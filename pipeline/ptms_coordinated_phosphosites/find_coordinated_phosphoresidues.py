#!/usr/bin/env python3
"""
find_coordinated_phosphoresidues.py
-----------------------------------
Finds phosphoserine, phosphothreonine and phosphotyrosine residues that are
coordinated by lysine or arginine side chains in high-resolution human crystal
structures, and measures how buried each one is.

Written for Reviewer #2, major point 6, which asks whether the coordinated
buried phosphosites predicted in this work resemble experimentally determined
structures.

Definition of "coordinated"
    A Lys or Arg side chain of the same chain whose charged group comes within
    4.0 Å of the phosphate group, measured between the closest pair of atoms
    of those two groups. This is the ionic-contact cut-off used elsewhere in
    this work (Supplementary Table 4). The charged group is NZ for Lys and
    NE, CZ, NH1, NH2 for Arg; CZ is included as the centre of the delocalised
    guanidinium. "Fully coordinated" means two or more such side chains.

    The closest approach between all non-hydrogen atoms of the two residues is
    reported alongside, as dist_heavy, but is not what the cut-off is applied
    to.

Definition of "buried"
    Relative solvent accessibility of 0.25 or less (Levy 2010), consistent
    with the alphasa database definition used throughout this work.
    Accessibility is computed with Hekkel mkdssp and normalised against the
    Tien et al. 2013 Empirical reference values, as in the main analysis.

Single-chain restriction
    Burial and coordination are both assessed within one chain at a time, with
    DSSP run on that chain in isolation. This matches the monomeric AlphaFold
    analysis in the main paper: a residue buried in the isolated chain is
    buried whether or not the protein has a binding partner, whereas a residue
    buried only by a crystallographic partner is not comparable. Accessibility
    in the deposited assembly is reported alongside for reference but is never
    used to call burial.

Reference accessibility for phosphoresidues
    Tien et al. 2013 cover the twenty standard residues only. The values used
    for SEP, TPO and PTR were derived by repeating their Gly-X-Gly enumeration;
    see phospho_reference_asa.py. Absolute accessibility is reported alongside
    relative, so that no conclusion depends on that derivation alone.

Input
-----
None. Structures, metadata, residue mappings and modification annotations are
retrieved from the RCSB Search and Data APIs, PDBe SIFTS and UniProt, and
cached under pdb_cache/ so that re-runs are offline and fast.

The search ignores structures released after RELEASE_CUTOFF below, so that it
returns the same set of entries however long after publication it is repeated.
Raising that date brings the survey up to date.

Output
------
TSV with one row per phosphoresidue-to-Lys/Arg contact, and one row with
n_contacts 0 for each phosphoresidue that has no such contact, so that every
entry in which a site is resolved can be counted:
    pdb_id, resolution, method, chain, chain_length, protein_name, acc, site,
    site_evidence, mapping, mapping_ok, single_segment, phos_res, phos_resnum,
    phos_icode, asa_chain, relasa_chain, buried, relasa_assembly, n_contacts,
    coord_res, coord_resnum, phos_atom, coord_atom, dist_charged, dist_heavy,
    is_human, is_peptide, include
Structures that could not be processed are listed in failed_structures.txt.

Usage
-----
    python find_coordinated_phosphoresidues.py
    python find_coordinated_phosphoresidues.py --out coordinated_phosphoresidues.tsv

Dependencies
------------
    pip install biopython requests tqdm
    mkdssp on the PATH (https://github.com/PDB-REDO/dssp)
"""

import argparse
import csv
import gzip
import json
import os
import subprocess
import tempfile
import time
import warnings
from collections import defaultdict
from datetime import date
from io import StringIO

import numpy as np
import requests
from Bio.PDB import MMCIFParser
from tqdm import tqdm

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DIST_CUTOFF        = 4.0   # Å, ionic contact (Supplementary Table 4)
RELASA_CUTOFF      = 0.25  # buried at or below this (Levy 2010)
MAX_RESOLUTION     = 2.0   # Å
MIN_CHAIN_LENGTH   = 30    # shorter chains are treated as substrate peptides

# Structures released after this date are ignored, so that the survey can be
# repeated and give the same answer as the PDB grows. Raise it to bring the
# survey up to date. As of this cut-off the search returns 1236 entries.
RELEASE_CUTOFF     = "2026-09-16"

PHOS_RESIDUES = frozenset({"SEP", "TPO", "PTR"})
PARENT_AA     = {"SEP": "S", "TPO": "T", "PTR": "Y"}

# Atoms carrying the negative charge on the phosphate. Both the modern (OP1)
# and the legacy (O1P) naming are accepted.
PHOSPHATE_ATOMS = frozenset({"P", "O1P", "O2P", "O3P", "OP1", "OP2", "OP3"})

# Atoms carrying the positive charge on the coordinating side chain.
CHARGED_ATOMS = {"LYS": frozenset({"NZ"}),
                 "ARG": frozenset({"NE", "CZ", "NH1", "NH2"})}

# Maximum accessible surface area used to normalise ASA.
#
# Standard residues: Tien et al. 2013, PLoS One 8(11):e80635, PMID 24278298,
# Empirical column, the same values used elsewhere in this work.
#
# Phosphorylated residues: derived by phospho_reference_asa.py, because no
# published values exist.
MAX_ASA: dict[str, float] = {
    "ALA": 121.0, "ARG": 265.0, "ASN": 187.0, "ASP": 187.0, "CYS": 148.0,
    "GLN": 214.0, "GLU": 214.0, "GLY":  97.0, "HIS": 216.0, "ILE": 195.0,
    "LEU": 191.0, "LYS": 230.0, "MET": 203.0, "PHE": 228.0, "PRO": 154.0,
    "SER": 143.0, "THR": 163.0, "TRP": 264.0, "TYR": 255.0, "VAL": 165.0,
    "SEP": 226.0, "TPO": 243.0, "PTR": 334.0,
}

# UniProt evidence codes, collapsed to three levels.
EXPERIMENTAL_ECO = frozenset({"ECO:0000269"})
SIMILARITY_ECO   = frozenset({"ECO:0000250", "ECO:0007744", "ECO:0000305"})

RCSB_SEARCH   = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_ENTRY    = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
RCSB_ENTITY   = "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
RCSB_CIF      = "https://files.rcsb.org/download/{pdb_id}.cif.gz"
SIFTS_URL     = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments/{pdb_id}"
UNIPROT_URL   = "https://rest.uniprot.org/uniprotkb/{acc}.json"

CACHE_DIR   = "pdb_cache"
OUTFILE     = "coordinated_phosphoresidues.tsv"
FAILED_FILE = "failed_structures.txt"

COLUMNS = ["pdb_id", "resolution", "method", "chain", "chain_length", "protein_name",
           "acc", "site", "site_evidence", "mapping", "mapping_ok", "single_segment",
           "phos_res", "phos_resnum", "phos_icode",
           "asa_chain", "relasa_chain", "buried", "relasa_assembly",
           "n_contacts", "coord_res", "coord_resnum", "phos_atom", "coord_atom",
           "dist_charged", "dist_heavy", "is_human", "is_peptide", "include"]

HEADER = f"""\
# Phosphoserine, phosphothreonine and phosphotyrosine residues coordinated by
# Lys or Arg in human X-ray structures at {MAX_RESOLUTION} Å resolution or better.
# One row per phosphoresidue-to-Lys/Arg contact, and one row with n_contacts 0 for each
# phosphoresidue with no such contact. Written by find_coordinated_phosphoresidues.py.
#
# structures_searched  {{structures_searched}}
# release_cutoff       {{release_cutoff}}
# retrieved            {{retrieved}}
#
# chain_length      residues in the polymer entity, or observed residues where unavailable
# acc, site         UniProt accession and position in its canonical sequence
# site_evidence     UniProt evidence for phosphorylation at that position
# mapping           how site was obtained: sifts, author_numbering, sifts_conflict or unmapped
# mapping_ok        the canonical sequence carries the expected Ser, Thr or Tyr at site
# single_segment    the chain maps to one contiguous UniProt segment
# phos_resnum       residue number as deposited, with phos_icode the insertion code
# asa_chain         accessible surface area (Å²) of the phosphoresidue, DSSP on the isolated chain
# relasa_chain      asa_chain over the reference maximum for that residue type
# buried            relasa_chain <= {RELASA_CUTOFF}
# relasa_assembly   the same quantity for the deposited assembly, reported but not used
# n_contacts        distinct Lys or Arg side chains of this chain coordinating the phosphate;
#                   0 for a phosphoresidue with none, whose coordination fields are empty
# dist_charged      Å between the closest atoms of the two charged groups, cut-off {DIST_CUTOFF}
# dist_heavy        Å between the closest non-hydrogen atoms of the two residues
# is_peptide        chain shorter than {MIN_CHAIN_LENGTH} residues, i.e. a substrate phosphopeptide
# include           is_human and mapping_ok and not is_peptide
"""


# ---------------------------------------------------------------------------
# Cached retrieval
# ---------------------------------------------------------------------------

def cache_path(*parts):
    """Path under CACHE_DIR, with the parent directory created."""
    path = os.path.join(CACHE_DIR, *parts)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return path


def cached_json(path, fetch):
    """Parsed JSON from path, fetching and storing it if absent."""
    if os.path.exists(path):
        try:
            with open(path) as fh:
                return json.load(fh)
        except ValueError:
            pass
    data = fetch()
    if data is not None:
        with open(path, "w") as fh:
            json.dump(data, fh)
    return data


def query_structures():
    """PDB IDs of human X-ray structures containing SEP, TPO or PTR."""
    residue_nodes = [
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers.chem_comp_monomers",
            "operator": "exact_match", "value": code}}
        for code in sorted(PHOS_RESIDUES)
    ]
    query = {
        "query": {"type": "group", "logical_operator": "and", "nodes": [
            {"type": "group", "logical_operator": "or", "nodes": residue_nodes},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_entity_source_organism.scientific_name",
                "operator": "exact_match", "value": "Homo sapiens"}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_entry_info.resolution_combined",
                "operator": "less_or_equal", "value": MAX_RESOLUTION}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_entry_info.experimental_method",
                "operator": "exact_match", "value": "X-ray"}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_accession_info.initial_release_date",
                "operator": "less_or_equal", "value": RELEASE_CUTOFF}},
        ]},
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": 10000},
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
        },
    }
    r = requests.post(RCSB_SEARCH, json=query, timeout=60)
    r.raise_for_status()
    return [hit["identifier"] for hit in r.json().get("result_set", [])]


def fetch_metadata(pdb_id):
    """Resolution, method, and per-chain protein name, accession and organism.

    Entity identifiers are read from the entry rather than assumed to run from
    1 to the number of protein entities, because nucleic acid and other
    polymer entities share the same numbering.
    """
    def fetch():
        r = requests.get(RCSB_ENTRY.format(pdb_id=pdb_id), timeout=30)
        r.raise_for_status()
        entry = r.json()
        info = entry.get("rcsb_entry_info", {})
        resolution = info.get("resolution_combined")
        meta = {"resolution": resolution[0] if isinstance(resolution, list) else resolution,
                "method": info.get("experimental_method", "unknown"),
                "entities": {}}
        entity_ids = entry.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])
        for entity_id in entity_ids:
            er = requests.get(RCSB_ENTITY.format(pdb_id=pdb_id, entity_id=entity_id), timeout=30)
            if not er.ok:
                continue
            entity = er.json()
            identifiers = entity.get("rcsb_polymer_entity_container_identifiers", {})
            accessions = [ref.get("database_accession")
                          for ref in identifiers.get("reference_sequence_identifiers", [])
                          if ref.get("database_name") == "UniProt"]
            organisms = entity.get("rcsb_entity_source_organism") or []
            meta["entities"][str(entity_id)] = {
                "name": entity.get("rcsb_polymer_entity", {}).get("pdbx_description", "unknown"),
                "accessions": accessions,
                "chains": identifiers.get("auth_asym_ids", []),
                "polymer_length": entity.get("entity_poly", {})
                                 .get("rcsb_sample_sequence_length"),
                "is_human": any(o.get("scientific_name", "").lower() == "homo sapiens"
                                for o in organisms),
            }
            time.sleep(0.02)
        return meta
    return cached_json(cache_path("metadata", f"{pdb_id}.json"), fetch)


def fetch_cif(pdb_id):
    """Structure in mmCIF format, from the cache if it is there."""
    path = cache_path("cif", f"{pdb_id.upper()}.cif")
    if os.path.exists(path):
        with open(path) as fh:
            return fh.read()
    r = requests.get(RCSB_CIF.format(pdb_id=pdb_id.lower()), timeout=90)
    r.raise_for_status()
    text = gzip.decompress(r.content).decode("utf-8")
    with open(path, "w") as fh:
        fh.write(text)
    return text


def fetch_sifts(pdb_id, retries=5):
    """SIFTS residue-level mapping, as {chain: {accession: [segments]}}."""
    def fetch():
        url = SIFTS_URL.format(pdb_id=pdb_id.lower())
        for attempt in range(retries):
            try:
                r = requests.get(url, timeout=60)
                if r.ok:
                    return r.json().get(pdb_id.lower(), {}).get("UniProt", {})
                if r.status_code == 404:
                    return {}
            except requests.RequestException:
                pass
            time.sleep(2.0 * (attempt + 1))
        return None

    raw = cached_json(cache_path("sifts", f"{pdb_id.upper()}.json"), fetch)
    if not raw:
        return {}
    mapping = defaultdict(lambda: defaultdict(list))
    for acc, block in raw.items():
        for seg in block.get("mappings", []):
            chain = seg.get("chain_id") or seg.get("struct_asym_id")
            unp_start, unp_end = seg.get("unp_start"), seg.get("unp_end")
            auth_start = seg.get("start", {}).get("author_residue_number")
            auth_end = seg.get("end", {}).get("author_residue_number")
            if chain is None or unp_start is None or unp_end is None:
                continue
            if auth_start is None and auth_end is None:
                continue
            mapping[chain][acc].append({"auth_start": auth_start, "auth_end": auth_end,
                                        "unp_start": int(unp_start), "unp_end": int(unp_end)})
    return {chain: dict(accs) for chain, accs in mapping.items()}


def fetch_uniprot(acc):
    """Canonical sequence, organism and annotated phosphosites for one accession."""
    def fetch():
        try:
            r = requests.get(UNIPROT_URL.format(acc=acc), timeout=30)
            return r.json() if r.ok else None
        except requests.RequestException:
            return None

    data = cached_json(cache_path("uniprot", f"{acc}.json"), fetch)
    if data is None:
        return None
    sites = {}
    for feature in data.get("features", []):
        if feature.get("type") != "Modified residue":
            continue
        if "phospho" not in feature.get("description", "").lower():
            continue
        position = feature.get("location", {}).get("start", {}).get("value")
        if position is None:
            continue
        codes = set()
        for evidence in feature.get("evidences", []):
            code = evidence.get("evidenceCode", "")
            codes.add(code.get("code", "") if isinstance(code, dict) else str(code))
        if codes & EXPERIMENTAL_ECO:
            sites[int(position)] = "experimental"
        elif codes & SIMILARITY_ECO:
            sites[int(position)] = "by_similarity"
        else:
            sites[int(position)] = "predicted"
    return {"sequence": data.get("sequence", {}).get("value", ""),
            "organism": data.get("organism", {}).get("scientificName", ""),
            "sites": sites}


# ---------------------------------------------------------------------------
# Residue mapping
# ---------------------------------------------------------------------------

def sifts_position(resnum, chain, acc, sifts):
    """Author residue number mapped to a position in the canonical sequence.

    A segment may be missing one of its author bounds, so the offset is taken
    from whichever bound is present and checked against the UniProt range.
    Returns (position, whether the chain maps as a single segment).
    """
    segments = sifts.get(chain, {}).get(acc)
    if not segments:
        return None, None
    single = len(segments) == 1
    for seg in segments:
        auth_start, auth_end = seg["auth_start"], seg["auth_end"]
        unp_start, unp_end = seg["unp_start"], seg["unp_end"]
        if auth_start is not None and auth_end is not None:
            if auth_start <= resnum <= auth_end:
                return unp_start + (resnum - auth_start), single
        elif auth_start is not None:
            position = unp_start + (resnum - auth_start)
            if unp_start <= position <= unp_end:
                return position, single
        elif auth_end is not None:
            position = unp_end - (auth_end - resnum)
            if unp_start <= position <= unp_end:
                return position, single
    return None, single


def assign_accession(record, entity, sifts, uniprot_cache):
    """Accession for this chain, and the phosphoresidue's position on it.

    A chain may reference more than one accession, either because the construct
    is a fusion or because SIFTS splits it. The accession chosen is the one
    whose SIFTS segment covers the phosphoresidue and whose canonical sequence
    carries the right residue there, which is also what catches mis-mapped
    positions.

    Returns (accession, position, single_segment, residue_matches, mapping).
    """
    expected = PARENT_AA[record["phos_res"]]
    candidates = entity.get("accessions") or []

    def sequence_of(acc):
        if acc not in uniprot_cache:
            uniprot_cache[acc] = fetch_uniprot(acc)
            time.sleep(0.05)
        entry = uniprot_cache[acc]
        return entry["sequence"] if entry else ""

    best = None
    for acc in candidates:
        position, single = sifts_position(record["phos_resnum"], record["chain"], acc, sifts)
        if position is None:
            continue
        sequence = sequence_of(acc)
        matches = bool(sequence) and 1 <= position <= len(sequence) \
            and sequence[position - 1] == expected
        score = (matches, single is True)
        if best is None or score > best[0]:
            best = (score, acc, position, single, matches)
    if best is not None and best[0][0]:
        _, acc, position, single, matches = best
        return acc, position, single, matches, "sifts"

    # SIFTS does not cover every deposited residue. Where it fails, the author
    # numbering is tested directly against the canonical sequence and accepted
    # only if the residue there is the right one, which rules out the offset
    # constructs that would otherwise split one site into two.
    for acc in candidates:
        sequence = sequence_of(acc)
        position = record["phos_resnum"]
        if sequence and 1 <= position <= len(sequence) and sequence[position - 1] == expected:
            return acc, position, None, True, "author_numbering"

    if best is not None:
        _, acc, position, single, matches = best
        return acc, position, single, matches, "sifts_conflict"
    acc = candidates[0] if candidates else None
    if acc:
        sequence_of(acc)
    return acc, None, None, False, "unmapped"


# ---------------------------------------------------------------------------
# Accessibility
# ---------------------------------------------------------------------------

def atom_site_chain_column(cif_text):
    """Index of the auth chain field within the _atom_site loop."""
    names, in_loop = [], False
    for line in cif_text.split("\n"):
        if line.startswith("_atom_site."):
            in_loop = True
            names.append(line.strip().split(".", 1)[1])
        elif in_loop:
            break
    columns = {name: i for i, name in enumerate(names)}
    return columns.get("auth_asym_id", columns.get("label_asym_id"))


def split_by_chain(cif_text, wanted):
    """One copy of the structure per wanted auth chain, as mmCIF text.

    Every line other than a coordinate line is kept, in place, so the file
    stays a valid mmCIF and only the atoms change. Doing all the chains in one
    pass keeps the cost down on entries with many of them.
    """
    chain_column = atom_site_chain_column(cif_text)
    if chain_column is None:
        return {}
    out = {chain: [] for chain in wanted}
    for line in cif_text.split("\n"):
        if line.startswith(("ATOM ", "HETATM")):
            fields = line.split()
            chain = fields[chain_column] if len(fields) > chain_column else None
            if chain in out:
                out[chain].append(line)
        else:
            for lines in out.values():
                lines.append(line)
    return {chain: "\n".join(lines) for chain, lines in out.items()}


def run_dssp(cif_text, workdir):
    """Accessible surface area per residue, keyed by "chain|resnum|icode".

    DSSP is the program used for accessibility throughout this work. It works
    on the polymer alone, so ordered waters and bound ligands do not shield
    the surface.
    """
    path = os.path.join(workdir, "in.cif")
    out = os.path.join(workdir, "out.dssp")
    with open(path, "w") as fh:
        fh.write(cif_text)
    if os.path.exists(out):
        os.remove(out)
    subprocess.run(["mkdssp", "--output-format", "dssp", path, out],
                   capture_output=True, text=True)
    if not os.path.exists(out):
        return {}
    with open(out) as fh:
        lines = fh.read().split("\n")
    start = next((i for i, line in enumerate(lines) if line.startswith("  #  RESIDUE")), None)
    if start is None:
        return {}
    asa = {}
    for line in lines[start + 1:]:
        if len(line) < 38 or line[13] == "!":
            continue
        try:
            resnum = int(line[5:10])
        except ValueError:
            continue
        asa[f"{line[11]}|{resnum}|{line[10].strip()}"] = int(line[34:38])
    os.remove(out)
    return asa


def structure_asa(pdb_id, cif_text, chains):
    """DSSP accessibility for the whole entry and for each chain in isolation.

    Results are cached on disk under pdb_cache/dssp/, which is what makes a
    re-run of the whole survey take seconds rather than minutes. The key "*"
    holds the deposited assembly; every other key is one auth chain.
    """
    path = cache_path("dssp", f"{pdb_id.upper()}.json")
    cached = {}
    if os.path.exists(path):
        try:
            with open(path) as fh:
                cached = json.load(fh)
        except ValueError:
            cached = {}

    missing = [key for key in ["*"] + sorted(chains) if key not in cached]
    if missing:
        by_chain = split_by_chain(cif_text, [key for key in missing if key != "*"])
        with tempfile.TemporaryDirectory() as workdir:
            for key in missing:
                text = cif_text if key == "*" else by_chain.get(key)
                cached[key] = run_dssp(text, workdir) if text else {}
        with open(path, "w") as fh:
            json.dump(cached, fh)
    return cached


# ---------------------------------------------------------------------------
# Contacts
# ---------------------------------------------------------------------------

def closest_approach(coords_a, coords_b):
    """Indices and distance of the closest pair of atoms between two sets."""
    distances = np.linalg.norm(coords_a[:, None, :] - coords_b[None, :, :], axis=-1)
    flat = int(np.argmin(distances))
    i, j = divmod(flat, distances.shape[1])
    return i, j, float(distances[i, j])


def analyse_structure(pdb_id, cif_text):
    """One record per phosphoresidue-to-Lys/Arg contact within a single chain."""
    try:
        structure = MMCIFParser(QUIET=True).get_structure(pdb_id, StringIO(cif_text))
        model = next(iter(structure))
    except Exception as exc:
        raise RuntimeError(f"parse failed: {exc}")

    phos_by_chain = {}
    for chain in model:
        residues = [r for r in chain if r.get_resname().strip() in PHOS_RESIDUES
                    and any(a.get_name() in PHOSPHATE_ATOMS for a in r)]
        if residues:
            phos_by_chain[chain.id] = residues
    if not phos_by_chain:
        return []

    asa = structure_asa(pdb_id, cif_text, phos_by_chain)
    assembly_asa = asa.get("*", {})

    records = []
    for chain in model:
        if chain.id not in phos_by_chain:
            continue
        chain_asa = asa.get(chain.id, {})
        basic = [r for r in chain if r.get_resname().strip() in CHARGED_ATOMS]
        observed_length = sum(1 for r in chain if r.get_resname().strip() in MAX_ASA)

        for residue in phos_by_chain[chain.id]:
            name = residue.get_resname().strip()
            resnum, icode = residue.get_id()[1], residue.get_id()[2].strip()
            key = f"{chain.id}|{resnum}|{icode}"
            if key not in chain_asa:
                continue

            phosphate_atoms = [a for a in residue if a.get_name() in PHOSPHATE_ATOMS]
            phosphate = np.array([a.get_coord() for a in phosphate_atoms])
            heavy = np.array([a.get_coord() for a in residue if a.element != "H"])

            contacts = {}
            for other in basic:
                if other is residue:
                    continue
                other_name = other.get_resname().strip()
                charged_atoms = [a for a in other if a.get_name() in CHARGED_ATOMS[other_name]]
                if not charged_atoms:
                    continue
                charged = np.array([a.get_coord() for a in charged_atoms])
                i, j, dist_charged = closest_approach(phosphate, charged)
                if dist_charged > DIST_CUTOFF:
                    continue
                other_heavy = np.array([a.get_coord() for a in other if a.element != "H"])
                _, _, dist_heavy = closest_approach(heavy, other_heavy)
                contacts[other.get_id()] = (other_name,
                                            phosphate_atoms[i].get_name(),
                                            charged_atoms[j].get_name(),
                                            dist_charged, dist_heavy)
            asa_chain = chain_asa[key]
            asa_assembly = assembly_asa.get(key)
            relasa_chain = min(asa_chain / MAX_ASA[name], 1.0)
            residue_fields = {
                "pdb_id": pdb_id,
                "chain": chain.id,
                "observed_length": observed_length,
                "phos_res": name,
                "phos_resnum": resnum,
                "phos_icode": icode,
                "asa_chain": asa_chain,
                "relasa_chain": round(relasa_chain, 3),
                "relasa_assembly": (round(min(asa_assembly / MAX_ASA[name], 1.0), 3)
                                    if asa_assembly is not None else ""),
                "buried": relasa_chain <= RELASA_CUTOFF,
                "n_contacts": len(contacts),
            }

            # A phosphoresidue with no Lys or Arg contact still gets one row, so
            # that every entry in which a site is resolved is counted when its
            # reproducibility is assessed, not only the entries where it happens
            # to be coordinated.
            if not contacts:
                records.append({**residue_fields, "coord_res": "", "coord_resnum": "",
                                "phos_atom": "", "coord_atom": "",
                                "dist_charged": "", "dist_heavy": ""})
                continue

            for other_id, (other_name, phos_atom, coord_atom,
                           dist_charged, dist_heavy) in contacts.items():
                records.append({**residue_fields,
                                "coord_res": other_name,
                                "coord_resnum": other_id[1],
                                "phos_atom": phos_atom,
                                "coord_atom": coord_atom,
                                "dist_charged": round(dist_charged, 2),
                                "dist_heavy": round(dist_heavy, 2)})
    return records


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--out", default=OUTFILE, help=f"output TSV (default {OUTFILE})")
    parser.add_argument("--limit", type=int, default=None,
                        help="process only the first N structures, for testing")
    args = parser.parse_args()

    print("Querying the PDB:")
    pdb_ids = query_structures()
    if args.limit:
        pdb_ids = pdb_ids[:args.limit]
    print(f" >> {len(pdb_ids)} structures")

    uniprot_cache = {}
    rows, failed = [], []

    for pdb_id in tqdm(pdb_ids, desc="Structures"):
        try:
            meta = fetch_metadata(pdb_id)
            sifts = fetch_sifts(pdb_id)
            records = analyse_structure(pdb_id, fetch_cif(pdb_id))
        except Exception as exc:
            failed.append((pdb_id, str(exc)))
            continue

        entity_of_chain = {chain: entity_id
                           for entity_id, entity in meta["entities"].items()
                           for chain in entity.get("chains", [])}

        for record in records:
            entity = meta["entities"].get(entity_of_chain.get(record["chain"], ""), {})
            acc, site, single_segment, mapping_ok, mapping = assign_accession(
                record, entity, sifts, uniprot_cache)

            entry = uniprot_cache.get(acc) if acc else None
            if entry is None:
                evidence = "no_uniprot_data"
            elif site is None:
                evidence = "unmapped"
            elif not mapping_ok:
                evidence = "mapping_conflict"
            else:
                evidence = entry["sites"].get(site, "not_annotated")

            # Humanness is judged on the accession the phosphoresidue itself
            # belongs to. A fusion construct can be part human and part not,
            # and in those the phosphoresidue may belong to the other half.
            is_human = bool(entry and entry.get("organism") == "Homo sapiens")

            # Chains too short to fold on their own are co-crystallised
            # substrate peptides, whose burial in the bound state says nothing
            # about burial in the protein they came from.
            chain_length = entity.get("polymer_length") or record["observed_length"]
            is_peptide = chain_length < MIN_CHAIN_LENGTH

            record.update({
                "protein_name": entity.get("name", "unknown"),
                "acc": acc,
                "site": site,
                "site_evidence": evidence,
                "mapping": mapping,
                "mapping_ok": mapping_ok,
                "single_segment": single_segment,
                "chain_length": chain_length,
                "resolution": meta["resolution"],
                "method": meta["method"],
                "is_human": is_human,
                "is_peptide": is_peptide,
                "include": is_human and not is_peptide and acc is not None and mapping_ok,
            })
            rows.append(record)

    if failed:
        with open(FAILED_FILE, "w") as fh:
            for pdb_id, error in failed:
                fh.write(f"{pdb_id}\t{error}\n")
        print(f" >> {len(failed)} structures failed, see {FAILED_FILE}")

    # Ties are broken on identity so that the row order is reproducible.
    rows.sort(key=lambda r: (r["relasa_chain"], -r["n_contacts"],
                             r["dist_charged"] if r["n_contacts"] else 0.0,
                             r["pdb_id"], r["chain"], r["phos_resnum"],
                             r["coord_resnum"] if r["n_contacts"] else 0))
    with open(args.out, "w", newline="") as fh:
        fh.write(HEADER.format(structures_searched=len(pdb_ids),
                               release_cutoff=RELEASE_CUTOFF,
                               retrieved=date.today().isoformat()))
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    contacts = [r for r in rows if r["n_contacts"]]
    kept = [r for r in contacts if r["include"]]
    print(f" >> {len(contacts)} contacts, {len(kept)} after filtering, written to {args.out}")


if __name__ == "__main__":
    main()
