"""
find_coordinated_phosphoresidues_all_rsa.py

Searches the PDB for HUMAN X-ray structures (≤ 2.0 Å resolution) containing
phosphoserine (SEP), phosphothreonine (TPO), or phosphotyrosine (PTR) where
the phosphoresidue is coordinated by one or more positively charged residues
(Lys or Arg) within 4.0 Å (P atom to closest charged atom), regardless of
burial status.

NMR structures are excluded: SASA computation on NMR ensembles is unreliable
because BioPython's ShrakeRupley processes all conformers simultaneously,
causing inter-conformer screening artefacts that render RSA values meaningless.

RSA is still computed and reported for downstream stratification.
buried_in_monomer is also computed (RSA in isolated chain vs full assembly).

Protein name and UniProt accession are assigned per-chain by fetching all
polymer entities for each structure, correctly handling multi-chain co-crystals.

Requirements:
    pip install biopython requests tqdm

Usage:
    python find_coordinated_phosphoresidues_all_rsa.py

Output:
    coordinated_phosphoresidues_all_rsa.tsv
    failed_structures_all_rsa.txt
"""

import logging
import os
import requests
import gzip
import time
import csv
from io import StringIO
from tqdm import tqdm
from Bio.PDB import MMCIFParser, SASA

# ── Parameters ────────────────────────────────────────────────────────────────

PHOS_RESIDUES   = {"SEP", "TPO", "PTR"}   # phosphoresidue 3-letter codes
                                          # all carry -2 charge at pH 7.4
COORD_RESIDUES  = {"LYS", "ARG"}          # positively charged residues (+1 each)
DISTANCE_CUTOFF = 4.0                     # Å, ionic contact cutoff (Supp Table 4)
                                          # measured P atom to closest charged atom
RSA_CUTOFF      = 0.25                    # relative SASA; <= this = buried

# Quality filters
MAX_RESOLUTION_XRAY        = 2.0          # Å; high-resolution cutoff for atom-level analysis

# Phosphorus atom name in each modified residue
PHOS_ATOM = "P"

# Charged nitrogen/carbon atoms on coordinating residues
COORD_ATOMS = {
    "LYS": ["NZ"],
    "ARG": ["NH1", "NH2", "CZ"],
}

# Reference ASA values (Å²) for RSA calculation (Tien et al. 2013, empirical)
# Modified residues mapped to closest standard residue:
#   SEP -> SER, TPO -> THR, PTR -> TYR
REF_ASA = {
    "ALA": 121.0, "ARG": 265.0, "ASN": 187.0, "ASP": 187.0,
    "CYS": 148.0, "GLN": 214.0, "GLU": 214.0, "GLY": 97.0,
    "HIS": 216.0, "ILE": 195.0, "LEU": 191.0, "LYS": 230.0,
    "MET": 203.0, "PHE": 228.0, "PRO": 154.0, "SER": 143.0,
    "THR": 163.0, "TRP": 264.0, "TYR": 255.0, "VAL": 165.0,
    "SEP": 143.0,   # -> SER
    "TPO": 163.0,   # -> THR
    "PTR": 255.0,   # -> TYR
}

log = logging.getLogger(__name__)

OUTPUT_FILE = "coordinated_phosphoresidues_all_rsa.tsv"
FAILED_FILE = "failed_structures_all_rsa.txt"
CACHE_DIR   = "pdb_cache"   # local cache for CIFs and API responses



# ── Cache helpers ─────────────────────────────────────────────────────────────

def _cache_path(*parts):
    """Return path inside CACHE_DIR, creating subdirs as needed."""
    path = os.path.join(CACHE_DIR, *parts)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return path

def _cache_get(path):
    """Return file contents as str if cached, else None."""
    try:
        with open(path) as f:
            return f.read()
    except FileNotFoundError:
        return None

def _cache_put(path, text):
    """Write text to cache file."""
    with open(path, "w") as f:
        f.write(text)


# ── Step 1: Query RCSB for human structures + fetch metadata ──────────────────

def _query_rcsb(query):
    """POST a query to the RCSB Search API and return list of PDB IDs."""
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    response = requests.post(url, json=query, timeout=30)
    response.raise_for_status()
    results = response.json()
    return [hit["identifier"] for hit in results.get("result_set", [])]


def count_total_human_structures(max_resolution_xray):
    """
    Counts total qualifying human X-ray structures in the PDB,
    plus unique human proteins (by UniProt accession) represented in those
    structures. Used as denominators for the phosphosite prevalence stats.
    """
    print("Counting total qualifying human structures and proteins (denominators)...")

    base_human = {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": "rcsb_entity_source_organism.scientific_name",
            "operator": "exact_match",
            "value": "Homo sapiens"
        }
    }

    def _count(extra_nodes):
        q = {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [base_human] + extra_nodes
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {"start": 0, "rows": 1},
                "results_content_type": ["experimental"],
            }
        }
        url = "https://search.rcsb.org/rcsbsearch/v2/query"
        r = requests.post(url, json=q, timeout=30)
        r.raise_for_status()
        return r.json().get("total_count", 0)

    def _count_proteins(extra_nodes):
        """Count unique UniProt accessions via polymer_entity return type."""
        q = {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [base_human] + extra_nodes
            },
            "return_type": "polymer_entity",
            "request_options": {
                "paginate": {"start": 0, "rows": 1},
                "results_content_type": ["experimental"],
            }
        }
        url = "https://search.rcsb.org/rcsbsearch/v2/query"
        r = requests.post(url, json=q, timeout=30)
        r.raise_for_status()
        return r.json().get("total_count", 0)

    xray_filter = [
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_entry_info.resolution_combined",
            "operator": "less_or_equal", "value": max_resolution_xray}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_entry_info.experimental_method",
            "operator": "exact_match", "value": "X-ray"}}
    ]
    n_xray_structs   = _count(xray_filter)
    n_xray_proteins  = _count_proteins(xray_filter)

    print(f"  X-ray ≤ {max_resolution_xray} Å: {n_xray_structs} structures, "
          f"~{n_xray_proteins} polymer entities")
    return n_xray_structs, n_xray_proteins


def query_pdb_xray_ids(phos_residues):
    """
    Queries RCSB for human X-ray (≤ MAX_RESOLUTION_XRAY Å) structures
    containing SEP, TPO, or PTR. Returns a list of PDB IDs sorted by
    resolution (best first). No per-structure metadata is fetched here —
    that is done lazily in the main processing loop alongside the CIF download.
    """
    residue_nodes = [
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers.chem_comp_monomers",
            "operator": "exact_match", "value": code}}
        for code in phos_residues
    ]
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {"type": "group", "logical_operator": "or", "nodes": residue_nodes},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_entity_source_organism.scientific_name",
                    "operator": "exact_match", "value": "Homo sapiens"}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_entry_info.resolution_combined",
                    "operator": "less_or_equal", "value": MAX_RESOLUTION_XRAY}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_entry_info.experimental_method",
                    "operator": "exact_match", "value": "X-ray"}},
            ]
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": 10000},
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined",
                      "direction": "asc"}]
        }
    }
    pdb_ids = _query_rcsb(query)
    print(f"  Found {len(pdb_ids)} X-ray structures with SEP/TPO/PTR.")
    return pdb_ids


def fetch_structure_metadata(pdb_id):
    """
    Fetches resolution and per-chain protein name/UniProt for one structure.
    Called from the main processing loop alongside the CIF download.
    Returns a dict: {resolution, method, protein_name, uniprot_id,
                     entity_map, chain_to_entity}.
    Results are cached to CACHE_DIR/metadata/{pdb_id}.json.
    """
    import json
    cache_file = _cache_path("metadata", f"{pdb_id}.json")
    cached = _cache_get(cache_file)
    if cached:
        result = json.loads(cached)
        # Check for fields added after initial cache was written.
        # If is_human is missing from entity_map entries, the cache is stale.
        entity_map = result.get("entity_map", {})
        first_entity = next(iter(entity_map.values()), {})
        if "is_human" in first_entity:
            return result
        # Stale cache — fall through to re-fetch and overwrite
        try:
            os.remove(cache_file)
        except OSError:
            pass

    # Resolution and method from entry endpoint
    data_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    r = requests.get(data_url, timeout=15)
    r.raise_for_status()
    d = r.json()
    method     = d.get("rcsb_entry_info", {}).get("experimental_method", "unknown")
    _res = d.get("rcsb_entry_info", {}).get("resolution_combined", None)
    resolution = _res[0] if isinstance(_res, list) else _res
    n_entities = d.get("rcsb_entry_info", {}).get("polymer_entity_count_protein", 1)

    # Per-entity protein name and UniProt, keyed by auth chain
    entity_map = {}
    for eid in range(1, max(int(n_entities) + 1, 2)):
        try:
            poly_url = (f"https://data.rcsb.org/rest/v1/core/"
                        f"polymer_entity/{pdb_id}/{eid}")
            pr = requests.get(poly_url, timeout=15)
            if not pr.ok:
                break
            pd_ = pr.json()
            ename = (pd_.get("rcsb_polymer_entity", {})
                        .get("pdbx_description", "unknown"))
            ref_ids = (pd_.get("rcsb_polymer_entity_container_identifiers", {})
                          .get("reference_sequence_identifiers", []))
            euniprot = next((ref.get("database_accession") for ref in ref_ids
                             if ref.get("database_name") == "UniProt"), None)
            auth_chains = (pd_.get("rcsb_polymer_entity_container_identifiers", {})
                              .get("auth_asym_ids", []))
            # Is this entity from Homo sapiens?
            # rcsb_entity_source_organism is a list; we check all entries.
            organisms = pd_.get("rcsb_entity_source_organism", [])
            is_human = any(
                org.get("scientific_name", "").lower() == "homo sapiens"
                for org in (organisms if isinstance(organisms, list) else [organisms])
            )
            entity_map[str(eid)] = {"protein_name": ename,
                                    "uniprot_id":   euniprot,
                                    "auth_chains":  auth_chains,
                                    "is_human":     is_human}
            time.sleep(0.02)
        except Exception:
            break

    if not entity_map:
        entity_map["1"] = {"protein_name": "unknown",
                           "uniprot_id": None, "auth_chains": [],
                           "is_human": False}

    chain_to_entity = {ch: eid
                       for eid, edata in entity_map.items()
                       for ch in edata.get("auth_chains", [])}
    default = entity_map.get("1", list(entity_map.values())[0])

    result = {
        "resolution":      resolution,
        "method":          method,
        "protein_name":    default["protein_name"],
        "uniprot_id":      default["uniprot_id"],
        "entity_map":      entity_map,
        "chain_to_entity": chain_to_entity,
    }
    _cache_put(cache_file, json.dumps(result))
    return result



# ── Step 2: Download and parse each structure ──────────────────────────────────

def download_cif(pdb_id):
    """Downloads the mmCIF file for a PDB entry. Returns file content as string.
    CIF is cached to CACHE_DIR/cif/{pdb_id}.cif to avoid re-downloading."""
    cache_file = _cache_path("cif", f"{pdb_id.upper()}.cif")
    cached = _cache_get(cache_file)
    if cached:
        return cached
    url = f"https://files.rcsb.org/download/{pdb_id.lower()}.cif.gz"
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    text = gzip.decompress(response.content).decode("utf-8")
    _cache_put(cache_file, text)
    return text


def fetch_uniprot_phosphosites(uniprot_id):
    """
    Fetches phosphorylation annotations from UniProt for the given accession.

    Returns a dict: {residue_position (int): evidence_level (str)}
    evidence_level is one of:
        'experimental'  — backed by ECO:0000269 (experimental evidence)
        'by_similarity' — inferred by similarity (ECO:0000250)
        'predicted'     — computational or other lower-confidence annotation

    Returns None if the request fails, so callers can distinguish
    "not annotated" from "fetch error".

    Results are cached to CACHE_DIR/uniprot/{uniprot_id}.json.
    """
    import json
    cache_file = _cache_path("uniprot", f"{uniprot_id}.json")
    cached = _cache_get(cache_file)
    if cached:
        try:
            data = json.loads(cached)
        except Exception:
            data = None
    else:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        try:
            r = requests.get(url, timeout=15)
            if not r.ok:
                return None
            data = r.json()
            _cache_put(cache_file, json.dumps(data))
        except Exception:
            return None

    phos_map = {}
    for feat in data.get("features", []):
        if feat.get("type") != "Modified residue":
            continue
        desc = feat.get("description", "").lower()
        if "phospho" not in desc:
            continue
        pos = feat.get("location", {}).get("start", {}).get("value")
        if pos is None:
            continue

        # Classify evidence quality from ECO codes
        ev_level = "predicted"
        for ev in feat.get("evidences", []):
            # UniProt API returns evidenceCode as either a dict {"code": "ECO:..."}
            # or directly as a string "ECO:..." depending on API version
            ec = ev.get("evidenceCode", "")
            code = ec.get("code", "") if isinstance(ec, dict) else str(ec)
            if code == "ECO:0000269":          # experimental evidence
                ev_level = "experimental"
                break
            elif code == "ECO:0000250":        # by similarity
                if ev_level != "experimental":
                    ev_level = "by_similarity"
            elif code in ("ECO:0007744",       # combinatorial evidence
                          "ECO:0000305"):       # curator inference
                if ev_level == "predicted":
                    ev_level = "by_similarity"
        phos_map[int(pos)] = ev_level

    return phos_map  # empty dict = valid response but no phosphosites annotated



def fetch_sifts_mappings(pdb_id, retries=10, backoff=2.0):
    """
    Fetches SIFTS segment-level residue mappings from PDBe for a structure.
    Results are cached to CACHE_DIR/sifts/{pdb_id}.json.

    Returns a dict:
        {chain_id: {uniprot_acc: [{"auth_start": int, "auth_end": int,
                                   "unp_start": int,  "unp_end": int}]}}

    Each entry is one contiguous mapped segment. A chain that maps to a single
    UniProt accession in a single segment is a clean, uninterrupted fragment.
    Multiple segments for the same accession indicate insertions, deletions, or
    numbering breaks within the deposited chain.

    Uses PDBe's SIFTS API: /pdbe/api/mappings/uniprot_segments/{pdb_id}
    Returns {} on failure.

    author_residue_number is preferred for matching PDB auth coords; falls back
    to residue_number (SEQRES label number) when author_residue_number is null.
    The same fallback is used in sifts_uniprot_pos so the match is consistent.
    """
    import json
    cache_file = _cache_path("sifts", f"{pdb_id.upper()}.json")
    cached = _cache_get(cache_file)
    if cached:
        data = json.loads(cached)
    else:
        url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments/{pdb_id.lower()}"
        data = {}
        for attempt in range(retries):
            try:
                r = requests.get(url, timeout=30)
                if r.ok:
                    data = r.json().get(pdb_id.lower(), {}).get("UniProt", {})
                    break
                elif r.status_code == 404:
                    _cache_put(cache_file, json.dumps({}))
                    return {}   # structure genuinely has no SIFTS mapping
                # Other HTTP errors: retry
            except Exception:
                pass
            if attempt < retries - 1:
                time.sleep(backoff * (attempt + 1))
        if not data:
            return {}
        _cache_put(cache_file, json.dumps(data))

    # Restructure: chain_id -> uniprot_acc -> [segments]
    chain_map = {}
    for uniprot_acc, acc_data in data.items():
        for seg in acc_data.get("mappings", []):
            chain_id = seg.get("chain_id") or seg.get("struct_asym_id")

            # Store auth_start and auth_end independently: null means unavailable
            # for that bound. Do NOT replace both when only one is null — that
            # would discard valid auth data and mix coordinate systems.
            auth_start = seg.get("start", {}).get("author_residue_number")
            auth_end   = seg.get("end",   {}).get("author_residue_number")

            unp_start  = seg.get("unp_start")
            unp_end    = seg.get("unp_end")
            if None in (chain_id, unp_start, unp_end):
                continue
            # Need at least one auth bound to compute a position
            if auth_start is None and auth_end is None:
                continue
            chain_map.setdefault(chain_id, {}).setdefault(uniprot_acc, []).append({
                "auth_start": int(auth_start) if auth_start is not None else None,
                "auth_end":   int(auth_end)   if auth_end   is not None else None,
                "unp_start":  int(unp_start),
                "unp_end":    int(unp_end),
            })
    return chain_map


def sifts_uniprot_pos(phos_resnum, chain_id, uniprot_acc, sifts):
    """
    Converts a PDB author residue number (auth_seq_id) to the corresponding
    UniProt canonical position using SIFTS segment data.

    Each segment has auth_start and/or auth_end (either may be None if the
    depositor did not provide author numbering for that bound). The UniProt
    position is computed using whichever auth bounds are available:

      - Both bounds present: standard linear interpolation within the range.
      - Only auth_start present: unp_pos = unp_start + (phos_resnum - auth_start),
        validated against unp_end.
      - Only auth_end present: unp_pos = unp_end - (auth_end - phos_resnum),
        validated against unp_start.

    This correctly handles structures like 6GLC chain A (auth_end=null) and
    3TMP chains C/E (auth_start present, auth_end=null).

    Returns (uniprot_pos, chain_is_contiguous):
        uniprot_pos         — int if computed and within [unp_start, unp_end],
                              None if phos_resnum falls outside all segments
        chain_is_contiguous — True if the chain maps to exactly one segment,
                              False if multiple segments, None if no SIFTS data
    """
    acc_segs = sifts.get(chain_id, {}).get(uniprot_acc)
    if acc_segs is None:
        return None, None

    is_contiguous = (len(acc_segs) == 1)

    for seg in acc_segs:
        a_start = seg["auth_start"]   # int or None
        a_end   = seg["auth_end"]     # int or None
        unp_start = seg["unp_start"]
        unp_end   = seg["unp_end"]

        if a_start is not None and a_end is not None:
            # Both bounds known: standard range check + linear offset
            if a_start <= phos_resnum <= a_end:
                return unp_start + (phos_resnum - a_start), is_contiguous

        elif a_start is not None:
            # Only start known: compute from start, validate against unp_end
            unp_pos = unp_start + (phos_resnum - a_start)
            if unp_start <= unp_pos <= unp_end:
                return unp_pos, is_contiguous

        elif a_end is not None:
            # Only end known: compute from end, validate against unp_start
            unp_pos = unp_end - (a_end - phos_resnum)
            if unp_start <= unp_pos <= unp_end:
                return unp_pos, is_contiguous

    return None, is_contiguous


def find_coordinated_sites(pdb_id, cif_content, phos_residues, coord_residues,
                           coord_atoms, distance_cutoff, rsa_cutoff):
    """
    Parses a structure and returns one dict per (phosphoresidue, coordinating
    residue) pair, deduplicated at the residue level: if a single Arg
    contributes contacts via both NH1 and NH2, only the closest atom distance
    is reported and it counts as one coordinating residue.

    This means a phosphoresidue coordinated by two Arg residues will appear
    on two rows, allowing downstream identification of fully coordinated
    (>= 2 positive charges) sites.
    """
    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure(pdb_id, StringIO(cif_content))
    except Exception:
        return []

    # RSA on full assembly (used for RSA column and burial check)
    try:
        sr = SASA.ShrakeRupley()
        sr.compute(structure, level="R")
    except Exception:
        return []

    # Build per-chain monomer RSA: isolate each chain and recompute SASA
    # so we can flag whether burial is intrinsic to the monomer rather than
    # arising from a protein–protein interface.
    from Bio.PDB import Model as PDBModel, Chain as PDBChain
    monomer_rsa = {}  # (chain_id, res_id) -> rsa_in_monomer
    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            # Build a minimal structure containing only this chain
            tmp_struct = structure.__class__("tmp")
            tmp_model  = PDBModel.Model(0)
            tmp_chain  = chain.copy()
            tmp_model.add(tmp_chain)
            tmp_struct.add(tmp_model)
            try:
                sr_mono = SASA.ShrakeRupley()
                sr_mono.compute(tmp_struct, level="R")
                for res in tmp_chain.get_residues():
                    abs_sasa_mono = getattr(res, "sasa", None)
                    resname = res.get_resname().strip()
                    ref = REF_ASA.get(resname)
                    if abs_sasa_mono is not None and ref:
                        monomer_rsa[(chain_id, res.get_id())] = \
                            min(round(abs_sasa_mono / ref, 3), 1.0)
            except Exception:
                pass  # monomer RSA unavailable for this chain; flag will be None

    results = []

    for model in structure:
        for chain in model:
            residues = list(chain.get_residues())
            chain_id = chain.get_id()

            # Collect ALL phosphoresidues (no RSA burial filter — include surface-exposed too)
            phos_sites = []
            for res in residues:
                resname = res.get_resname().strip()
                if resname not in phos_residues:
                    continue
                if PHOS_ATOM not in res:
                    continue
                ref = REF_ASA.get(resname)
                if not ref:
                    continue
                abs_sasa = getattr(res, "sasa", None)
                if abs_sasa is None:
                    continue
                rsa = min(round(abs_sasa / ref, 3), 1.0)  # cap at 1.0 (empirical refs can give >1)
                phos_sites.append((res, rsa))  # RSA retained for output, not used as filter

            if not phos_sites:
                continue

            # Collect coordinating residues (all atoms per residue)
            coord_res_list = []
            for res in residues:
                resname = res.get_resname().strip()
                if resname in coord_residues:
                    coord_res_list.append(res)

            for phos_res, rsa in phos_sites:
                p_coord = phos_res[PHOS_ATOM].get_vector()
                phos_resnum = phos_res.get_id()[1]
                phos_resname = phos_res.get_resname().strip()

                # For each coordinating residue, find closest atom distance
                # (residue-level deduplication)
                contacts_by_coord_res = {}  # coord_resnum -> (resname, atom, dist)
                for coord_res in coord_res_list:
                    if coord_res == phos_res:
                        continue
                    coord_resname = coord_res.get_resname().strip()
                    coord_resnum  = coord_res.get_id()[1]
                    best_atom = None
                    best_dist = float("inf")
                    for atom_name in coord_atoms.get(coord_resname, []):
                        if atom_name in coord_res:
                            dist = (p_coord - coord_res[atom_name].get_vector()).norm()
                            if dist < best_dist:
                                best_dist = dist
                                best_atom = atom_name
                    if best_dist <= distance_cutoff:
                        contacts_by_coord_res[coord_resnum] = (
                            coord_resname, best_atom, round(best_dist, 2)
                        )

                if not contacts_by_coord_res:
                    continue

                n_coord_residues = len(contacts_by_coord_res)

                # One row per coordinating residue
                for coord_resnum, (coord_resname, coord_atom, dist) in \
                        contacts_by_coord_res.items():
                    # Composite ranking score: lower = more buried + closer contact
                    score = round((rsa / rsa_cutoff) + (dist / distance_cutoff), 4)
                    # Is burial intrinsic to the monomer (not from an interface)?
                    mono_rsa = monomer_rsa.get((chain_id, phos_res.get_id()))
                    buried_in_monomer = (
                        bool(mono_rsa <= rsa_cutoff) if mono_rsa is not None else None
                    )
                    results.append({
                        "burial_contact_score_lower_is_better": score,
                        "buried_in_monomer":  buried_in_monomer,
                        "pdb_id":             pdb_id,
                        "chain":              chain_id,
                        "chain_length":       len(residues),
                        "phos_resname":       phos_resname,
                        "phos_resnum":        phos_resnum,
                        "rsa":                rsa,
                        "n_coord_residues":   n_coord_residues,
                        "coord_res":          coord_resname,
                        "coord_resnum":       coord_resnum,
                        "coord_atom":         coord_atom,
                        "distance_A":         dist,
                    })

    return results


# ── Step 3: Main ───────────────────────────────────────────────────────────────

def main():
    import sys

    # Get denominator counts first (quick queries, no structure parsing)
    n_xray_structs, n_xray_proteins = \
        count_total_human_structures(MAX_RESOLUTION_XRAY)

    print("Querying RCSB for X-ray structures with SEP/TPO/PTR...")
    pdb_ids = query_pdb_xray_ids(PHOS_RESIDUES)
    print(f"\n── Phosphosite prevalence in human PDB (X-ray only) ──")
    print(f"  X-ray ≤ {MAX_RESOLUTION_XRAY} Å with SEP/TPO/PTR: "
          f"{len(pdb_ids)} structures ({100*len(pdb_ids)/n_xray_structs:.1f}% "
          f"of {n_xray_structs})")
    print()

    all_hits = []
    failed   = []

    for pdb_id in tqdm(pdb_ids, desc="Processing structures"):
        try:
            meta = fetch_structure_metadata(pdb_id)
            sifts = fetch_sifts_mappings(pdb_id)
            cif_content = download_cif(pdb_id)
            hits = find_coordinated_sites(
                pdb_id, cif_content,
                PHOS_RESIDUES, COORD_RESIDUES, COORD_ATOMS,
                DISTANCE_CUTOFF, RSA_CUTOFF
            )
            # Attach metadata: use chain-specific entity for protein/UniProt
            chain_to_entity = meta.get("chain_to_entity", {})
            entity_map      = meta.get("entity_map", {})
            for h in hits:
                eid   = chain_to_entity.get(h.get("chain", ""), "1")
                edata = entity_map.get(eid, entity_map.get("1", {}))
                h["resolution"]   = meta.get("resolution")
                h["method"]       = meta.get("method")
                h["protein_name"] = edata.get("protein_name",
                                       meta.get("protein_name", "unknown"))
                uid = edata.get("uniprot_id", meta.get("uniprot_id"))
                h["uniprot_id"]   = uid
                h["is_human"]     = edata.get("is_human", False)
    
                # SIFTS: map PDB author resnum -> UniProt canonical position
                # and flag whether the chain is a single contiguous fragment
                if uid and sifts:
                    unp_pos, contiguous = sifts_uniprot_pos(
                        h["phos_resnum"], h["chain"], uid, sifts)
                    h["uniprot_phos_resnum"]    = unp_pos
                    h["chain_is_contiguous"]    = contiguous
                else:
                    h["uniprot_phos_resnum"]    = None
                    h["chain_is_contiguous"]    = None

            # ── Flag and filter substrate peptide chains and non-human entries ──
            # Short chains (< 30 residues) or "peptide" in name are substrate
            # peptides used in co-crystal binding studies; their RSA in the
            # deposited (bound) conformation is not meaningful.
            # Non-human chains (viral, designed, synthetic) are flagged via
            # is_human=False from the entity organism field.
            # All hits are kept in the output for transparency, but flagged
            # with is_canonical_human=False so they can be excluded from stats.
            n_before = len(hits)
            for h in hits:
                is_peptide = (
                    h.get("chain_length", 9999) < 30
                    or "peptide" in h.get("protein_name", "").lower()
                    or not h.get("uniprot_id")
                )
                h["is_canonical_human"] = (
                    h.get("is_human", False) and not is_peptide
                )
            n_flagged = sum(1 for h in hits if not h["is_canonical_human"])
            if n_flagged:
                print(f"    [{pdb_id}] flagged {n_flagged} hit(s) as non-canonical "
                      f"(non-human or peptide chain)")

            all_hits.extend(hits)
            time.sleep(0.05)
        except Exception as e:
            failed.append((pdb_id, str(e)))

    if failed:
        print(f"  ({len(failed)} structures failed — see {FAILED_FILE})")
        with open(FAILED_FILE, "w") as f:
            for pid, err in failed:
                f.write(f"{pid}\t{err}\n")



    # ── UniProt phosphosite validation ──────────────────────────────────────────
    # For each hit, check whether the phosphosite position is annotated in
    # UniProt as a phosphorylation site. Results are cached per UniProt accession
    # to avoid redundant API calls (many PDB structures share the same protein).
    #
    # Caveat: UniProt uses canonical sequence numbering; PDB author numbering
    # (auth_seq_id) is used here. These usually match but may differ for
    # truncated constructs or non-standard numbering. uniprot_phos_evidence=
    # 'not_annotated' may therefore reflect a numbering offset rather than a
    # truly unvalidated site.
    print("\nValidating phosphosites against UniProt annotations...")
    uniprot_cache = {}  # uniprot_id -> {pos: evidence_level} or None (fetch error)

    for h in all_hits:
        uid = h.get("uniprot_id")
        if not uid:
            h["uniprot_phos_evidence"] = "no_uniprot_id"
            continue
        if uid not in uniprot_cache:
            uniprot_cache[uid] = fetch_uniprot_phosphosites(uid)
            time.sleep(0.05)
        phos_map = uniprot_cache[uid]
        if phos_map is None:
            h["uniprot_phos_evidence"] = "fetch_error"
        else:
            # Use SIFTS-mapped UniProt position if available, else fall back
            # to PDB author numbering (with caveat about potential offset)
            lookup_pos = h.get("uniprot_phos_resnum") or h["phos_resnum"]
            ev = phos_map.get(int(lookup_pos))
            if ev is not None:
                h["uniprot_phos_evidence"] = ev
            elif h.get("uniprot_phos_resnum") is None:
                # SIFTS couldn't map this residue - numbering unknown
                h["uniprot_phos_evidence"] = "sifts_unmapped"
            else:
                h["uniprot_phos_evidence"] = "not_annotated"

    # Print validation summary
    from collections import Counter
    ev_counts = Counter(h["uniprot_phos_evidence"] for h in all_hits)
    # Unique (uniprot_id, phos_resnum) sites by evidence level
    site_ev = {}
    for h in all_hits:
        key = (h.get("uniprot_id"), h["phos_resnum"])
        # Upgrade evidence if multiple rows for the same site disagree
        prev = site_ev.get(key, "not_annotated")
        curr = h["uniprot_phos_evidence"]
        priority = {"experimental": 4, "by_similarity": 3, "predicted": 2,
                    "not_annotated": 1, "fetch_error": 0, "no_uniprot_id": -1}
        site_ev[key] = curr if priority.get(curr, 0) > priority.get(prev, 0) else prev

    ev_site_counts = Counter(site_ev.values())
    print(f"\nUniProt phosphosite validation (unique uniprot+resnum sites):")
    for ev in ["experimental", "by_similarity", "predicted",
               "not_annotated", "sifts_unmapped", "fetch_error", "no_uniprot_id"]:
        n = ev_site_counts.get(ev, 0)
        if n:
            print(f"  {ev:20s}: {n} sites")
    print(f"  NOTE: 'not_annotated' uses SIFTS-corrected UniProt positions — "
          f"these are genuinely absent from UniProt, not numbering artefacts.")
    print(f"  'sifts_unmapped' = residue falls outside SIFTS-mapped segments "
          f"(e.g. expression tags, engineered residues).")

    # Contiguity summary
    total_hits = len(all_hits)
    contiguous     = sum(1 for h in all_hits if h.get("chain_is_contiguous") is True)
    non_contiguous = sum(1 for h in all_hits if h.get("chain_is_contiguous") is False)
    no_sifts       = sum(1 for h in all_hits if h.get("chain_is_contiguous") is None)
    print(f"\nChain contiguity (single uninterrupted UniProt segment):")
    print(f"  Contiguous:     {contiguous} hits")
    print(f"  Non-contiguous: {non_contiguous} hits (insertions/breaks in chain)")
    print(f"  No SIFTS data:  {no_sifts} hits")

    n_structures = len({h["pdb_id"] for h in all_hits})
    site_n_coord = {}
    for h in all_hits:
        key = (h["pdb_id"], h["chain"], h["phos_resnum"])
        rsa     = h.get("rsa") or h.get("rsa_0to1_lower_is_more_buried")
        n_coord = h.get("n_coord_residues") or h.get("n_coord_residues_higher_is_more_coordinated")
        site_n_coord[key] = (n_coord, rsa)
    n_sites_total  = len(site_n_coord)
    fully_coord    = sum(1 for n, _   in site_n_coord.values() if n is not None and n >= 2)
    n_buried       = sum(1 for _, rsa in site_n_coord.values() if rsa is not None and rsa <= RSA_CUTOFF)
    n_surface      = sum(1 for _, rsa in site_n_coord.values() if rsa is not None and rsa >  RSA_CUTOFF)
    fully_coord_buried  = sum(1 for n, rsa in site_n_coord.values()
                              if n is not None and rsa is not None and n >= 2 and rsa <= RSA_CUTOFF)
    fully_coord_surface = sum(1 for n, rsa in site_n_coord.values()
                              if n is not None and rsa is not None and n >= 2 and rsa >  RSA_CUTOFF)
    print(f"\nFound {n_sites_total} coordinated phosphoresidue sites "
          f"across {n_structures} structures (all hits, including flagged).")
    print(f"  Buried (RSA ≤ {RSA_CUTOFF}):   {n_buried} sites "
          f"({fully_coord_buried} fully coordinated with ≥ 2 K/R)")
    print(f"  Surface (RSA > {RSA_CUTOFF}):  {n_surface} sites "
          f"({fully_coord_surface} fully coordinated with ≥ 2 K/R)")
    print(f"  Total fully coordinated (≥ 2 K/R contacts): {fully_coord}")

    # ── Canonical-human-only statistics (for rebuttal) ────────────────────────
    # Deduplication: unique sites are defined by (uniprot_id, canonical_resnum),
    # where canonical_resnum = uniprot_phos_resnum if available (SIFTS-mapped),
    # else phos_resnum (PDB deposited number). This matches the deduplication
    # used in classify_kinase_activation_loops.py and the rebuttal numbers.
    canon = [h for h in all_hits if h.get("is_canonical_human")]
    n_flagged_total = len(all_hits) - len(canon)

    # Build per-site aggregation
    site_data = {}  # (uniprot_id, canonical_resnum) -> dict of aggregated values
    for h in canon:
        uid       = h.get("uniprot_id")
        up_res    = h.get("uniprot_phos_resnum")
        pdb_res   = h.get("phos_resnum")
        canon_res = int(up_res) if up_res is not None else int(pdb_res)
        key       = (uid, canon_res)

        rsa     = h.get("rsa") or h.get("rsa_0to1_lower_is_more_buried")
        n_coord = h.get("n_coord_residues") or h.get("n_coord_residues_higher_is_more_coordinated", 0) or 0
        bim     = h.get("buried_in_monomer")  # True / False / None
        evid    = h.get("uniprot_phos_evidence", "")
        pdb_id  = h["pdb_id"]

        if key not in site_data:
            site_data[key] = {
                "max_ncoord":     n_coord,
                "buried_in_any":  bim is True,
                "surface_in_any": bim is False,
                "evidence":       evid,
                "pdb_ids_bur":    set(),
                "pdb_ids_surf":   set(),
            }
        else:
            site_data[key]["max_ncoord"]     = max(site_data[key]["max_ncoord"], n_coord)
            if bim is True:
                site_data[key]["buried_in_any"]  = True
            if bim is False:
                site_data[key]["surface_in_any"] = True

        # Track which structures each site appears in as buried / surface
        if bim is True:
            site_data[key]["pdb_ids_bur"].add(pdb_id)
        if bim is False:
            site_data[key]["pdb_ids_surf"].add(pdb_id)

        # Upgrade evidence level
        priority = {"experimental": 4, "by_similarity": 3, "predicted": 2,
                    "sifts_unmapped": 1, "not_annotated": 0, "fetch_error": -1}
        if priority.get(evid, -2) > priority.get(site_data[key]["evidence"], -2):
            site_data[key]["evidence"] = evid

    # Structure-level counts
    n_canon_structs   = len({h["pdb_id"] for h in canon})
    bur_structs_set   = {h["pdb_id"] for h in canon if h.get("buried_in_monomer") is True}
    surf_structs_set  = {h["pdb_id"] for h in canon if h.get("buried_in_monomer") is False}
    bur_fc_structs_set = {h["pdb_id"] for h in canon
                          if h.get("buried_in_monomer") is True and
                          (h.get("n_coord_residues") or
                           h.get("n_coord_residues_higher_is_more_coordinated") or 0) >= 2}

    # Site-level counts
    all_sites   = list(site_data.values())
    n_sites     = len(all_sites)
    n_proteins  = len({k[0] for k in site_data})

    bur_sites   = [s for s in all_sites if s["buried_in_any"]]
    surf_sites  = [s for s in all_sites if s["surface_in_any"]]
    surf_only   = [s for s in all_sites if s["surface_in_any"] and not s["buried_in_any"]]

    bur_fc      = [s for s in bur_sites  if s["max_ncoord"] >= 2]
    bur_sc      = [s for s in bur_sites  if s["max_ncoord"] == 1]
    surf_fc     = [s for s in surf_sites if s["max_ncoord"] >= 2]
    surf_sc     = [s for s in surf_sites if s["max_ncoord"] == 1]
    all_fc      = [s for s in all_sites  if s["max_ncoord"] >= 2]
    all_fc_surf_only = [s for s in surf_only if s["max_ncoord"] >= 2]

    bur_fc_proteins = len({k[0] for k, s in site_data.items()
                           if s["buried_in_any"] and s["max_ncoord"] >= 2})
    bur_proteins    = len({k[0] for k, s in site_data.items() if s["buried_in_any"]})
    surf_proteins   = len({k[0] for k, s in site_data.items() if s["surface_in_any"]})

    up_supported = sum(1 for s in all_sites
                       if s["evidence"] in ("experimental", "by_similarity"))

    print(f"\n{'='*65}")
    print(f"CANONICAL HUMAN STATISTICS — REBUTTAL NUMBERS")
    print(f"({'='*65})")
    print(f"(Excludes {n_flagged_total} non-canonical/peptide hits; "
          f"sites deduplicated by UniProt position)")
    print()
    print(f"OVERALL")
    print(f"  Structures:                         {n_canon_structs}")
    print(f"  Unique coordinated phosphosites:    {n_sites}")
    print(f"  Unique proteins (UniProt IDs):       {n_proteins}")
    print(f"  UniProt-supported sites             ")
    print(f"    (experimental or by_similarity):  {up_supported}")
    print()
    print(f"SURFACE (RSA > {RSA_CUTOFF}, present in ≥1 structure)")
    print(f"  Structures:                         {len(surf_structs_set)}")
    print(f"  Unique sites:                       {len(surf_sites)}")
    print(f"  Unique proteins:                    {surf_proteins}")
    print(f"  Fully coordinated (≥2 K/R):         {len(surf_fc)}")
    print(f"  Singly coordinated (=1 K/R):        {len(surf_sc)}")
    print(f"  Check: {len(surf_fc)}+{len(surf_sc)} = {len(surf_fc)+len(surf_sc)} "
          f"(should equal {len(surf_sites)})")
    print()
    print(f"BURIED (RSA ≤ {RSA_CUTOFF}, present in ≥1 structure)")
    print(f"  Structures:                         {len(bur_structs_set)}")
    print(f"  Unique sites:                       {len(bur_sites)}")
    print(f"  Unique proteins:                    {bur_proteins}")
    print(f"  Fully coordinated (≥2 K/R):")
    print(f"    Unique sites:                     {len(bur_fc)}")
    print(f"    Structures:                       {len(bur_fc_structs_set)}")
    print(f"    Proteins:                         {bur_fc_proteins}")
    print(f"  Singly coordinated (=1 K/R):        {len(bur_sc)}")
    print(f"  Check: {len(bur_fc)}+{len(bur_sc)} = {len(bur_fc)+len(bur_sc)} "
          f"(should equal {len(bur_sites)})")
    print()
    print(f"ALL FULLY COORDINATED (≥2 K/R, any burial status)")
    print(f"  Unique sites total:                 {len(all_fc)}")
    print(f"    of which buried (any struct):     {len(bur_fc)}")
    print(f"    of which surface-only:            {len(all_fc_surf_only)}")
    print(f"  Check: {len(bur_fc)}+{len(all_fc_surf_only)} = "
          f"{len(bur_fc)+len(all_fc_surf_only)} (should equal {len(all_fc)})")

    if failed:
        pass  # already written inside the else block above

    if all_hits:
        # Sort: monomer-buried first, then by composite score, then n_coord desc
        all_hits.sort(key=lambda h: (
            0 if h.get("buried_in_monomer") is True  else
            1 if h.get("buried_in_monomer") is None  else 2,
            h.get("burial_contact_score_lower_is_better", 0),
            -(h.get("n_coord_residues") or h.get("n_coord_residues_higher_is_more_coordinated") or 0),
        ))

        fieldnames = [
            "burial_contact_score_lower_is_better",
            "buried_in_monomer",
            "pdb_id", "protein_name", "uniprot_id", "resolution_A", "method",
            "chain", "chain_length", "phos_resname", "phos_resnum", "uniprot_phos_resnum",
            "rsa_0to1_lower_is_more_buried",
            "n_coord_residues_higher_is_more_coordinated",
            "coord_res", "coord_resnum", "coord_atom",
            "contact_distance_A_lower_is_closer",
            "uniprot_phos_evidence",
            "chain_is_contiguous",
            "is_human",
            "is_canonical_human",
        ]
        col_map = {
            "burial_contact_score_lower_is_better": "burial_contact_score_lower_is_better",
            "buried_in_monomer":       "buried_in_monomer",
            "pdb_id":                  "pdb_id",
            "protein_name":            "protein_name",
            "uniprot_id":              "uniprot_id",
            "resolution":              "resolution_A",
            "method":                  "method",
            "chain":                   "chain",
            "chain_length":            "chain_length",
            "phos_resname":            "phos_resname",
            "phos_resnum":             "phos_resnum",
            "uniprot_phos_resnum":     "uniprot_phos_resnum",
            "rsa":                     "rsa_0to1_lower_is_more_buried",
            "n_coord_residues":        "n_coord_residues_higher_is_more_coordinated",
            "coord_res":               "coord_res",
            "coord_resnum":            "coord_resnum",
            "coord_atom":              "coord_atom",
            "distance_A":              "contact_distance_A_lower_is_closer",
            "uniprot_phos_evidence":   "uniprot_phos_evidence",
            "chain_is_contiguous":     "chain_is_contiguous",
            "is_human":                "is_human",
            "is_canonical_human":      "is_canonical_human",
        }
        renamed_hits = [
            {col_map[k]: v for k, v in h.items()} for h in all_hits
        ]
        with open(OUTPUT_FILE, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(renamed_hits)
        print(f"\nResults written to: {OUTPUT_FILE}")
        print(f"\nTop 10 hits:")
        print("\t".join(fieldnames))
        for h in renamed_hits[:10]:
            print("\t".join(str(h[k]) for k in fieldnames))
    else:
        print("No hits found — try relaxing DISTANCE_CUTOFF or RSA_CUTOFF.")


if __name__ == "__main__":
    main()
