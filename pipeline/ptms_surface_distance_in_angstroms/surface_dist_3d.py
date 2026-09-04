#!/usr/bin/env python3
"""
surface_dist_3d.py
------------------
Computes the 3D Euclidean distance (Å) from each buried PTM residue's
sidechain heavy atoms to the nearest surface-exposed atom in the
corresponding AlphaFold structure.

This is the structural depth metric requested by Reviewers #1 and #2,
complementing the sequence-space `surfdist` metric in distance_to_surface.R.

Definition of "surface-exposed atom"
    Any heavy atom in a residue whose relative SASA (computed with
    Hekkel mkdssp, normalised against Tien et al. 2013 Empirical reference
    values, 1.4 Å probe) exceeds 0.25 (Levy 2010 threshold), consistent
    with the alphasa database definition.

Definition of "depth" for a buried PTM residue
    The minimum Euclidean distance (Å) from any sidechain heavy atom of
    the PTM residue (CB included; backbone N/CA/C/O/OXT excluded;
    fallback to CA for Gly) to any surface-exposed heavy atom in the
    structure, queried via a cKDTree for efficiency.

Input
-----
tmp-qm.rds from alphasa_relasa_vs_ptms.R.
All rows are PTM sites (the SQL query hardcodes 'Modified'; no ptmbin filter).
Required columns (exact names, no aliases):
    acc      – UniProt accession (canonical isoform)
    site     – residue position, 1-indexed, matching AlphaFold chain A
    ptm      – PTM type string  (e.g. "S-p")
    source   – data source      (e.g. "Ochoa", "PhosphoSitePlus")
    relasa   – relative SASA from alphasa  (NOT relasa10)
    plddt    – per-residue pLDDT
    min_pae  – minimum PAE over all non-neighbour (|Δsite| > 1) contacts

Filters applied (matching the main paper criterion; see Methods and
alphasa_relasa_vs_ptms.R line ~604)
    relasa <= 0.25         buried (Levy 2010)
    plddt  >= 70           high local confidence
    min_pae <= 2           at least one confident non-neighbour contact
                           (PAE ≤ 2 Å; non-neighbour defined as |Δsite| > 1)

Structure download
    Two sources are tried in priority order:
    1. Local AlphaSync v0 files ({local_cif_dir}/AF-{acc}-F1-model_v0.cif.gz)
       These exactly match the structures used to build the alphasa database.
    2. AFDB GCS v4 via gcloud storage cp
       gs://public-datasets-deepmind-alphafold-v4/AF-{acc}-F1-model_v4.cif
       Requires: brew install --cask google-cloud-sdk && gcloud auth login
    Raises RuntimeError if neither source succeeds.
    CIFs are cached as {acc}.cif and gemmi-converted PDBs as {acc}.pdb.

Output
------
CSV with columns:
    acc, site, ptm, source, relasa, plddt, min_pae, surf_dist_3d_ang

Usage
-----
    python surface_dist_3d.py --rds tmp-qm.rds --out output-surfdist-3d.csv

Dependencies
------------
    pip install pyreadr MDAnalysis scipy tqdm gemmi matplotlib
    mkdssp next to this script (https://github.com/PDB-REDO/dssp)
    brew install --cask google-cloud-sdk  &&  gcloud auth login
"""

import argparse
import gzip
import logging
import os
import shutil
import subprocess
import threading
import warnings
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import pyreadr
from scipy.spatial import cKDTree
from tqdm import tqdm

import gemmi
import MDAnalysis as mda

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

AFDB_GCS_BUCKET    = "gs://public-datasets-deepmind-alphafold-v4"
AFDB_GCS_CIF_PATH  = AFDB_GCS_BUCKET + "/AF-{acc}-F1-model_v4.cif"

# Default path to local AlphaSync CIF files (AF 2.3.2, "v0").
# These take priority over all remote sources.
# Files are named:  AF-{acc}-F1-model_v0.cif.gz
DEFAULT_LOCAL_CIF_DIR = (
    "alphasync_cif"
)

# RSA threshold matching the alphasa / Levy 2010 definition used throughout
# the paper:  surface = relASA > 0.25,  buried = relASA <= 0.25.
SURFACE_RSA_THRESHOLD = 0.25

# Backbone atom names excluded when selecting the sidechain of a PTM residue.
BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O", "OXT"})

# Tien et al. (2013) Empirical maximum ASA values per residue type (Å²).
# Source: PLoS One 8(11):e80635  PMID: 24278298
# These are the same reference values used to populate the alphasa database
# (combine_fragments_dssp.py, "Empirical" column of max_asa_for_residues.tsv).
MAX_ASA_TIEN_EMPIRICAL: dict[str, float] = {
    "A": 121.0, "R": 265.0, "N": 187.0, "D": 187.0, "C": 148.0,
    "E": 214.0, "Q": 214.0, "G":  97.0, "H": 216.0, "I": 195.0,
    "L": 191.0, "K": 230.0, "M": 203.0, "F": 228.0, "P": 154.0,
    "S": 143.0, "T": 163.0, "W": 264.0, "Y": 255.0, "V": 165.0,
}

# mkdssp stderr patterns that are spurious AlphaFold citation warnings.
# Identical to the grep filter used in job_dssp.py.
_MKDSSP_IGNORABLE_WARNINGS = (
    "Links for citation_author:citation:1 are incomplete",
    "There are 33 items in citation_author that don",
    "Warning, the input file is not valid",
)

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# MDAnalysis emits very verbose INFO messages about segids, types, and masses
# for every PDB file loaded.  Silence these — errors and warnings still show.
for _mda_logger in ("MDAnalysis", "MDAnalysis.topology", "MDAnalysis.coordinates"):
    logging.getLogger(_mda_logger).setLevel(logging.WARNING)


# ---------------------------------------------------------------------------
# Data loading and filtering
# ---------------------------------------------------------------------------

REQUIRED_COLUMNS = ["acc", "site", "aa", "ptm", "source",
                    "relasa", "plddt", "min_pae"]


def load_and_filter_rds(rds_path: str) -> pd.DataFrame:
    """
    Load tmp-qm.rds and apply the canonical buried-PTM filter:
        relasa <= 0.25        AND
        plddt  >= 70          AND
        min_pae <= 2

    All rows are assumed to be PTM sites (no ptmbin column present).
    Column names must match exactly; no aliases are accepted.
    Raises ValueError if any required column is absent.
    Returns a deduplicated (acc, site) DataFrame.
    """
    log.info("Loading %s ...", rds_path)
    result = pyreadr.read_r(rds_path)
    df = result[None] if None in result else result[list(result.keys())[0]]
    log.info("  Loaded %d rows; columns: %s", len(df), list(df.columns))

    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(
            f"Required columns missing from RDS: {missing}\n"
            f"Available columns: {list(df.columns)}\n"
            "Use tmp-qm.rds (alphasa_relasa_vs_ptms.R), not q.rds — "
            "q.rds does not contain 'relasa' or 'min_pae'."
        )

    n0 = len(df)

    # multi-fragment proteins by seeing if F2+ structures exist.
    log.info("  Proceeding with structure-based fragment detection ...")

    # 2. Row-level filtering for analysis quality
    df = df[df["relasa"].astype(float) <= SURFACE_RSA_THRESHOLD]
    log.info("  After relasa <= %.2f:           %d rows", SURFACE_RSA_THRESHOLD, len(df))

    df = df[df["plddt"].astype(float) >= 70.0]
    log.info("  After plddt >= 70:             %d rows", len(df))

    # min_pae filter handles missing PAE rows (mixed or otherwise).
    df = df[df["min_pae"].astype(float) <= 2.0]
    log.info("  After min_pae <= 2:            %d rows", len(df))

    return df


# ---------------------------------------------------------------------------
# Structure download
# ---------------------------------------------------------------------------

def prefetch_structures(
    accs: list[str],
    struct_dir: Path,
    local_cif_dir: Path,
    max_structure_ts: float | None = None,
    skipped_proteins: dict[str, str] | None = None,
) -> set[str]:
    """
    Ensure CIF files are present in *struct_dir* for all accessions in *accs*
    before the main processing loop begins, so the loop itself never blocks on
    I/O.  Returns the set of accessions whose CIF could not be obtained.

    Strategy:
      1. Decompress all available local v0 .cif.gz files (fast, Python gzip).
      2. For everything still missing, issue a single  gcloud storage cp  call
         with all required GCS URIs at once — one subprocess startup instead
         of one per protein.

    Already-present .cif files (from a previous run) are skipped in both steps.
    """
    if skipped_proteins is None:
        skipped_proteins = {}

    # ── Step 1: decompress local v0 files ────────────────────────────────────
    need_gcs: list[str] = []
    n_local = 0
    n_date_skipped = 0
    del_counts = {"cif": 0, "pdb": 0, "dssp": 0}
    # Track WHY we need GCS for each accession to give better diagnostics on failure
    gcs_reasons: dict[str, str] = {} # acc -> reason
    
    for acc in tqdm(accs, desc="Local CIFs", unit="acc"):
        local_gz = local_cif_dir / f"AF-{acc}-F1-model_v0.cif.gz" if local_cif_dir else None
        cif_dest = struct_dir / f"{acc}.cif"
        
        use_local = False
        if local_gz and local_gz.exists():
            mtime = local_gz.stat().st_mtime
            if max_structure_ts is None or mtime <= max_structure_ts:
                use_local = True
            else:
                log.debug("  -> Local v0 for %s is newer than cutoff, falling back to GCS.", acc)
                n_date_skipped += 1
                gcs_reasons[acc] = "local v0 too new (> cutoff)"
        else:
            gcs_reasons[acc] = "missing on local drive"

        if use_local:
            # Check for local fragments
            if (local_cif_dir / f"AF-{acc}-F2-model_v0.cif.gz").exists():
                log.info("  -> Skipping %s: Multi-fragment protein detected (found local -F2)", acc)
                skipped_proteins[acc] = "multi-fragment protein (found local -F2)"
                if acc in gcs_reasons: del gcs_reasons[acc]
                continue

            if not cif_dest.exists() or cif_dest.stat().st_size == 0:
                try:
                    with gzip.open(local_gz, "rb") as f_in:
                        cif_dest.write_bytes(f_in.read())
                    # Preserve timestamp from local file
                    os.utime(cif_dest, (mtime, mtime))
                    n_local += 1
                    if acc in gcs_reasons: del gcs_reasons[acc]
                except Exception as e:
                    log.error("Failed to decompress local CIF for %s: %s", acc, e)
                    cif_dest.unlink(missing_ok=True)
                    need_gcs.append(acc)
                    gcs_reasons[acc] = "local decompression failed"
            else:
                n_local += 1
                if acc in gcs_reasons: del gcs_reasons[acc]
            continue
        else:
            # Not using local (either absent or too new). Check cache or add to GCS list.
            if not cif_dest.exists() or cif_dest.stat().st_size == 0:
                need_gcs.append(acc)
                # reason already set above
            elif max_structure_ts is not None and cif_dest.stat().st_mtime > max_structure_ts:
                # Cached file is too new. We'll mark for re-fetch from GCS
                # and delete the current files to ensure a clean refresh.
                log.info("  -> Cached CIF for %s is newer than cutoff, re-fetching.", acc)
                cif_dest.unlink()
                del_counts["cif"] += 1
                for ext in ["pdb", "dssp"]:
                    fpath = struct_dir / f"{acc}.{ext}"
                    if fpath.exists():
                        fpath.unlink()
                        del_counts[ext] += 1
                need_gcs.append(acc)
                gcs_reasons[acc] = "cached file too new (> cutoff)"
            else:
                # Cached file is fine
                n_local += 1
                if acc in gcs_reasons: del gcs_reasons[acc]

    log.info(
        "Prefetch step 1/2 done: %d copied/cached, %d need GCS.",
        n_local, len(need_gcs),
    )
    if n_date_skipped > 0:
        log.info("  -> %d local v0 files skipped (newer than cutoff).", n_date_skipped)
    if sum(del_counts.values()) > 0:
        log.info("  -> Cleaned up %d outdated cached files (cif=%d, pdb=%d, dssp=%d).",
                 sum(del_counts.values()), del_counts["cif"], del_counts["pdb"], del_counts["dssp"])

    # ── Step 2: bulk GCS download for everything still missing ───────────────
    failed: dict[str, str] = {} # acc -> reason
    if need_gcs:
        log.info(
            "Prefetch step 2/2: downloading %d CIFs from GCS in one batch ...",
            len(need_gcs),
        )
        existing_v4_files = list(struct_dir.glob("AF-*-F1-model_v4.cif*"))
        acc_set = set(need_gcs)
        for f in existing_v4_files:
            fname = f.name
            if not fname.startswith("AF-"): continue
            parts = fname.split("-")
            if len(parts) > 1 and parts[1] in acc_set:
                try: f.unlink()
                except OSError: pass

        gcs_uris_str = "\n".join(AFDB_GCS_CIF_PATH.format(acc=acc).replace("-F1-", "-F*-") for acc in need_gcs)
        
        # Check local fragments again for anything that became need_gcs due to date skip
        final_skip = set()
        for acc in need_gcs:
             if (local_cif_dir / f"AF-{acc}-F2-model_v0.cif.gz").exists():
                 final_skip.add(acc)
        if final_skip:
            log.info("  -> Discarding %d additional proteins found to have -F2 fragments locally.", len(final_skip))
            for acc in final_skip:
                skipped_proteins[acc] = "multi-fragment protein (found local -F2)"
            need_gcs = [a for a in need_gcs if a not in final_skip]
            if not need_gcs: return set()
            gcs_uris_str = "\n".join(AFDB_GCS_CIF_PATH.format(acc=acc).replace("-F1-", "-F*-") for acc in need_gcs)

        gcloud_result: dict = {}
        def _run_gcloud() -> None:
            try:
                proc = subprocess.run(
                    ["gcloud", "-q", "storage", "cp", "--preserve-posix", "-I", str(struct_dir) + "/"],
                    input=gcs_uris_str, text=True, check=False, capture_output=True,
                )
                gcloud_result["returncode"] = proc.returncode
            except FileNotFoundError:
                gcloud_result["error"] = "gcloud not found."

        gcloud_thread = threading.Thread(target=_run_gcloud, daemon=True)
        gcloud_thread.start()

        def _is_valid(f):
            try: return f.stat().st_size > 0
            except OSError: return False

        with tqdm(total=len(need_gcs), desc="GCS download", unit="cif") as pbar:
            seen = 0
            while gcloud_thread.is_alive():
                current = sum(1 for f in struct_dir.glob("AF-*-F1-model_v4.cif") if _is_valid(f))
                if current > seen:
                    pbar.update(current - seen)
                    seen = current
                gcloud_thread.join(timeout=0.5)
            current = sum(1 for f in struct_dir.glob("AF-*-F1-model_v4.cif") if _is_valid(f))
            if current > seen: pbar.update(current - seen)

        gcloud_thread.join()
        if "error" in gcloud_result: raise RuntimeError(gcloud_result["error"])
        
        n_gcs_ok = 0
        for acc in need_gcs:
            cif_dest = struct_dir / f"{acc}.cif"
            f1_gcs = struct_dir / f"AF-{acc}-F1-model_v4.cif"
            higher_fragments = list(struct_dir.glob(f"AF-{acc}-F[2-9]-model_v4.cif"))
            
            if higher_fragments:
                log.info("  -> Skipping %s: Multi-fragment protein detected on GCS (F2+ found)", acc)
                skipped_proteins[acc] = "multi-fragment protein on GCS"
                for df in higher_fragments + [f1_gcs]:
                    try: df.unlink()
                    except: pass
                continue

            if f1_gcs.exists() and f1_gcs.stat().st_size > 0:
                f1_gcs.rename(cif_dest)
            
            if cif_dest.exists() and cif_dest.stat().st_size > 0:
                n_gcs_ok += 1
            else:
                reason = gcs_reasons.get(acc, "unknown")
                log.error("No CIF obtained for %s (%s, and not on GCS)", acc, reason)
                failed[acc] = reason

        log.info("Prefetch step 2/2 done: %d/%d downloaded from GCS, %d failed.",
                 n_gcs_ok, len(need_gcs), len(failed))

    n_total = len(accs)
    log.info("Prefetch complete: %d/%d structures ready (%d copies, %d GCS, %d unavailable).",
             n_total - len(failed), n_total, n_local, n_total - n_local - len(failed), len(failed))
    
    if failed:
        log.warning("Unavailable structures (first 20):")
        for acc, reason in list(failed.items())[:20]:
            log.warning("  - %s: %s", acc, reason)
    return set(failed.keys())


def download_structure(acc: str, struct_dir: Path) -> Path:
    """
    Convert the already-cached {acc}.cif to {acc}.pdb via gemmi and return
    the PDB path.  Assumes prefetch_structures() has already been called.
    Already-present PDB files are returned immediately.
    """
    pdb_dest = struct_dir / f"{acc}.pdb"
    if pdb_dest.exists() and pdb_dest.stat().st_size > 0:
        return pdb_dest

    cif_dest = struct_dir / f"{acc}.cif"
    if not cif_dest.exists() or cif_dest.stat().st_size == 0:
        raise RuntimeError(
            f"CIF not found for {acc} — prefetch_structures() may have failed"
        )

    try:
        structure = gemmi.read_structure(str(cif_dest))
        structure.write_pdb(str(pdb_dest))
    except Exception as e:
        raise RuntimeError(
            f"gemmi CIF→PDB conversion failed for {acc}: {e}"
        ) from e

    return pdb_dest


# ---------------------------------------------------------------------------
# mkdssp runner
# ---------------------------------------------------------------------------

def _run_mkdssp(input_path: Path, dssp_path: Path) -> None:
    """
    Invoke mkdssp on *input_path* (mmCIF), writing DSSP output to *dssp_path*.

    mkdssp v4 (Hekkel/PDB-REDO) requires mmCIF input — PDB is not accepted.
    AlphaFold-specific citation warnings in stderr are suppressed (matching
    the filter used in job_dssp.py).  All other stderr output is forwarded
    to the DEBUG log.  Raises RuntimeError on non-zero exit or empty output.
    """
    # Prefer the dssp422 conda environment; then an mkdssp executable sitting next to this script; fall back to PATH.
    conda_dssp = Path.home() / ".local" / "share" / "mamba" / "envs" / "dssp422" / "bin" / "mkdssp"
    local_bin  = Path(__file__).resolve().parent / "mkdssp"
    
    if conda_dssp.exists():
        mkdssp_bin = str(conda_dssp)
    elif local_bin.exists():
        mkdssp_bin = str(local_bin)
    else:
        mkdssp_bin = "mkdssp"

    # mkdssp (via Miniforge/conda) needs LIBCIFPP_DATA_DIR to locate its
    # dictionary files (mmcif_pdbx.dic, mmcif_ma.dic, etc.).  When called
    # outside the conda environment the variable is unset and mkdssp fails.
    # Derive the data directory from the resolved executable location:
    #   .../miniforge/bin/mkdssp  →  .../miniforge/share/libcifpp
    env = os.environ.copy()
    if "LIBCIFPP_DATA_DIR" not in env:
        script_dir = Path(__file__).resolve().parent
        # First: look for a libcifpp folder next to the script itself
        local_data = script_dir / "libcifpp"
        # Second: derive from the mkdssp binary location (e.g. Miniforge)
        #   .../miniforge/bin/mkdssp  →  .../miniforge/share/libcifpp
        conda_data = Path(mkdssp_bin).resolve().parent.parent / "share" / "libcifpp"
        for candidate in (local_data, conda_data):
            if candidate.is_dir():
                env["LIBCIFPP_DATA_DIR"] = str(candidate)
                log.debug("Set LIBCIFPP_DATA_DIR=%s", candidate)
                break

    try:
        proc = subprocess.run(
            [mkdssp_bin, str(input_path), str(dssp_path)],
            capture_output=True,
            text=True,
            check=False,
            env=env,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(
            "mkdssp executable not found on PATH. "
            "Install from https://github.com/PDB-REDO/dssp"
        ) from exc

    # Forward non-trivial stderr (excluding known AlphaFold citation noise)
    noisy_lines = [
        line for line in proc.stderr.splitlines()
        if line.strip() and not any(p in line for p in _MKDSSP_IGNORABLE_WARNINGS)
    ]
    if noisy_lines:
        log.debug("mkdssp stderr for %s:\n%s", input_path.name, "\n".join(noisy_lines))

    if proc.returncode != 0:
        raise RuntimeError(
            f"mkdssp exited with code {proc.returncode} for {input_path.name}. "
            f"stderr: {proc.stderr[:500]}"
        )

    if not dssp_path.exists() or dssp_path.stat().st_size == 0:
        raise RuntimeError(
            f"mkdssp produced no output for {input_path.name}"
        )


# ---------------------------------------------------------------------------
# DSSP parser
# ---------------------------------------------------------------------------

def _parse_dssp(dssp_path: Path) -> list[tuple[str, int, str, int]]:
    """
    Parse a DSSP output file (Kabsch & Sander 1983 fixed-width format).

    Returns a list of (chain, resnum, aa_1letter, asa_absolute_Å2) tuples,
    one per residue, in sequence order.  Chain-break pseudo-residues ('!')
    are skipped.

    Fixed-width column positions used (0-indexed, per Kabsch & Sander 1983):
        [5:10]  PDB residue sequence number
        [11]    chain ID
        [13]    amino acid (1-letter code)
        [34:38] solvent-accessible surface area (Å²)
    """
    records: list[tuple[str, int, str, int]] = []
    in_data = False

    with open(dssp_path) as fh:
        for line in fh:
            # The data section begins after the column-header line
            if line.startswith("  #  RESIDUE"):
                in_data = True
                continue
            if not in_data or len(line) < 38:
                continue
            # Chain-break marker inserted by DSSP between chain segments
            if line[13] == "!":
                continue
            resnum_str = line[5:10].strip()
            if not resnum_str:
                continue
            try:
                resnum = int(resnum_str)
                chain  = line[11]
                aa     = line[13]
                asa    = int(line[34:38])
            except ValueError:
                continue

            records.append((chain, resnum, aa, asa))

    return records


# ---------------------------------------------------------------------------
# Surface residue identification
# ---------------------------------------------------------------------------

def get_surface_resids(pdb_path: Path) -> set[int]:
    """
    Run mkdssp on the sibling {stem}.cif (mkdssp v4 requires mmCIF input),
    cache the output as {stem}.dssp, normalise the per-residue absolute ASA
    with Tien et al. (2013) Empirical maxima, and return the set of residue
    numbers whose relASA exceeds SURFACE_RSA_THRESHOLD (> 0.25).

    Non-standard amino acids (not in MAX_ASA_TIEN_EMPIRICAL) are skipped.
    Raises RuntimeError if mkdssp fails or no surface residues are found.
    """
    cif_path  = pdb_path.with_suffix(".cif")
    dssp_path = pdb_path.with_suffix(".dssp")
    if not dssp_path.exists() or dssp_path.stat().st_size == 0:
        _run_mkdssp(cif_path, dssp_path)

    records = _parse_dssp(dssp_path)
    if not records:
        raise RuntimeError(
            f"DSSP parse yielded no residues for {pdb_path.name}"
        )

    surface_resids: set[int] = set()
    for (_chain, resnum, aa, asa_abs) in records:
        max_asa = MAX_ASA_TIEN_EMPIRICAL.get(aa)
        if max_asa is None or max_asa <= 0:
            continue
        if (asa_abs / max_asa) > SURFACE_RSA_THRESHOLD:
            surface_resids.add(resnum)

    if not surface_resids:
        raise RuntimeError(
            f"No surface-exposed residues (relASA > {SURFACE_RSA_THRESHOLD}) "
            f"found in {pdb_path.name}"
        )

    return surface_resids


# ---------------------------------------------------------------------------
# Surface atom coordinate array
# ---------------------------------------------------------------------------

def get_surface_atom_coords(
    u: mda.Universe,
    surface_resids: set[int],
    pdb_name: str,
) -> np.ndarray:
    """
    Return an (N, 3) float32 array of heavy-atom coordinates for all residues
    in *surface_resids*, selected from the already-loaded Universe *u*.

    AlphaFold/AlphaSync PDB files are single-chain, so residue numbers are
    matched by resid alone without chain filtering.

    Raises RuntimeError if the coordinate array is empty.
    """
    surf_coords: list[np.ndarray] = []
    for atom in u.select_atoms("not name H*"):
        if atom.resid in surface_resids:
            surf_coords.append(atom.position)

    if not surf_coords:
        raise RuntimeError(
            f"Surface atom selection is empty for {pdb_name} "
            f"({len(surface_resids)} surface residues, "
            f"{len(u.select_atoms('not name H*'))} heavy atoms total)"
        )

    return np.array(surf_coords, dtype=np.float32)


# ---------------------------------------------------------------------------
# Per-residue sidechain coordinates
# ---------------------------------------------------------------------------

def sidechain_heavy_coords(
    u: mda.Universe, resid: int, acc: str
) -> np.ndarray | None:
    """
    Return (N, 3) float32 sidechain heavy-atom positions for *resid*
    (1-indexed) in Universe *u*.

    Sidechain = all heavy atoms except backbone N, CA, C, O, OXT.
    Falls back to CA for Glycine (no sidechain heavy atoms).
    Returns None and logs a warning if the residue is not found.
    """
    backbone_sel = " or name ".join(sorted(BACKBONE_ATOMS))
    sel = u.select_atoms(
        f"resid {resid} and not (name {backbone_sel}) and not name H*"
    )
    if len(sel) == 0:
        # Glycine fallback
        sel = u.select_atoms(f"resid {resid} and name CA")
    if len(sel) == 0:
        log.warning("Residue %d not found in structure for %s", resid, acc)
        return None
    return sel.positions.astype(np.float32)


# ---------------------------------------------------------------------------
# QC reporting
# ---------------------------------------------------------------------------

# Distance thresholds used to flag unusual values.
QC_ZERO_DIST    = 0.0   # exact zero — sidechain atom at same coords as surface atom
QC_NEAR_ZERO    = 1.0   # < 1 Å — right at the buried/exposed boundary
QC_LARGE_DIST   = 15.0  # > 15 Å — unusually deep burial


def report_qc(qc_events: list[dict]) -> None:
    """
    Log a structured QC summary from the *qc_events* list accumulated
    during the main processing loop.

    Categories tracked:
        zero_dist            surf_dist_3d_ang == 0.0 Å
                             Sidechain atom at identical coordinates to a
                             surface atom — physically impossible, likely a
                             boundary residue or gemmi PDB conversion artefact.

        near_zero_dist       0 < surf_dist_3d_ang < 1.0 Å
                             Residue is at the buried/exposed boundary.
                             alphasa called it buried (relASA ≤ 0.25) but the
                             sidechain barely clears the surface.

        large_dist           surf_dist_3d_ang > 15.0 Å
                             Unusually deeply buried residue — rare but possible
                             in large multi-domain proteins.

        dssp_alphasa_conflict  DSSP (mkdssp, Tien et al.) classifies this
                             residue as surface-exposed (relASA > 0.25), but
                             alphasa (the database source of the input filter)
                             called it buried.  Indicates a discrepancy between
                             the two SASA calculations, e.g. from different DSSP
                             versions, reference values, or multi-fragment
                             averaging in alphasa.  The distance is still
                             computed (treating the residue as buried per the
                             input data) but is likely to be zero or very small.
    """
    if not qc_events:
        log.info("QC report: no anomalies detected.")
        return

    by_cat: dict[str, list[dict]] = {}
    for e in qc_events:
        by_cat.setdefault(e["category"], []).append(e)

    total = len(qc_events)
    log.warning("QC report — %d anomalous (acc, site) observations:", total)

    descriptions = {
        "missing_residue":       "residue not found in structure (version mismatch?)",
        "aa_mismatch":           "amino acid mismatch between input data and structure",
        "zero_dist":             "exact zero distance (sidechain overlaps surface atom)",
        "near_zero_dist":        "near-zero distance < 1.0 Å (boundary residue)",
        "large_dist":            "large distance > 15.0 Å (unusually deep burial)",
        "dssp_alphasa_conflict": "DSSP/alphasa conflict (DSSP says surface, alphasa says buried)",
    }
    for cat, events in sorted(by_cat.items(), key=lambda x: -len(x[1])):
        desc = descriptions.get(cat, cat)
        examples = ", ".join(
            f"{e['acc']}:{e['site']}({e['ptm']})" for e in events[:5]
        )
        suffix = f" ... (+{len(events) - 5} more)" if len(events) > 5 else ""
        log.warning(
            "  %-28s  n=%d   examples: %s%s",
            desc, len(events), examples, suffix,
        )




# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main(args: argparse.Namespace) -> None:
    struct_dir    = Path(args.struct_dir)
    out_path      = Path(args.out)
    local_cif_dir = Path(args.local_cif_dir)

    # ── --reset: wipe struct_dir and output CSV, then exit ───────────────────
    if args.reset:
        if out_path.exists():
            out_path.unlink()
            log.info("Deleted %s", out_path)
        if struct_dir.exists():
            counts: dict[str, int] = {}
            for f in struct_dir.iterdir():
                ext = f.suffix.lower()
                counts[ext] = counts.get(ext, 0) + 1
            shutil.rmtree(struct_dir)
            summary = ", ".join(
                f"{n} {ext}" for ext, n in sorted(counts.items())
            ) or "0 files"
            log.info("Deleted %s/ (%s)", struct_dir, summary)
        log.info("Reset complete. Continuing with a fresh run ...")

    struct_dir.mkdir(parents=True, exist_ok=True)

    df = load_and_filter_rds(args.rds)

    if "ptmbin" not in df.columns:
        log.info("No ptmbin column found in input RDS. Defaulting all rows to 'Modified'.")
        df["ptmbin"] = "Modified"

    # ── Resumability: skip proteins already fully written to the output CSV ──
    done_accs: set[str] = set()
    if out_path.exists() and out_path.stat().st_size > 0:
        try:
            existing  = pd.read_csv(out_path, usecols=["acc"])
            done_accs = set(existing["acc"].unique())
            log.info(
                "Resuming — %d proteins already in %s, skipping them.",
                len(done_accs), out_path,
            )
        except Exception as e:
            log.warning("Could not read existing output file (%s); starting fresh.", e)

    unique_accs = [a for a in df["acc"].unique() if a not in done_accs]
    log.info(
        "Processing %d unique proteins (%d remaining after resume skip) ...",
        df["acc"].nunique(), len(unique_accs),
    )

    # ── Bulk prefetch all structures before processing ───────────────────────
    skipped_proteins:  dict[str, str] = {}   # acc → reason

    if args.max_structure_date:
        max_ts = time.mktime(time.strptime(args.max_structure_date, "%Y-%m-%d"))
        max_ts += 86399  # include the whole day
    else:
        max_ts = None

    prefetch_failed = prefetch_structures(
        unique_accs, struct_dir, local_cif_dir, max_ts, skipped_proteins
    )

    # ── Open output CSV for incremental appending ────────────────────────────
    OUTPUT_COLS  = ["acc", "site", "aa", "ptm", "source", "ptmbin", "dis",
                    "relasa", "plddt", "min_pae", "surf_dist_3d_ang"]
    write_header = not (out_path.exists() and out_path.stat().st_size > 0)

    failed_structures: list[str]   = list(prefetch_failed)
    failed_residues:   list[tuple] = []
    qc_events:         list[dict]  = []
    total_written = len(done_accs)

    # Build a 1-letter AA lookup from the input data for residue verification.
    # The RDS has an "aa" column with the expected amino acid per site.
    has_aa_col = "aa" in df.columns

    # Skip proteins whose CIF could not be obtained during prefetch, OR were skipped due to date
    accs_to_process = [a for a in unique_accs if a not in prefetch_failed and a not in skipped_proteins]

    with open(out_path, "a", newline="") as out_fh:
        for i, acc in enumerate(tqdm(accs_to_process, desc="Proteins", unit="protein")):
            protein_rows = df[df["acc"] == acc]

            # ── 1. Convert cached CIF → PDB via gemmi ─────────────────────────
            try:
                pdb_path = download_structure(acc, struct_dir)
            except RuntimeError as e:
                log.error("STRUCTURE UNAVAILABLE -- %s: %s", acc, e)
                failed_structures.append(acc)
                continue

            # ── 2. Load MDAnalysis Universe ────────────────────────────────────
            try:
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", message=".*CRYST1.*")
                    u = mda.Universe(str(pdb_path))
            except Exception as e:
                log.error("MDAnalysis load failed -- %s: %s", acc, e)
                failed_structures.append(acc)
                continue

            # ── 3. Run mkdssp, identify surface residues, build KDTree ─────────
            try:
                surface_resids = get_surface_resids(pdb_path)
                surf_coords    = get_surface_atom_coords(u, surface_resids, acc)
                surf_tree      = cKDTree(surf_coords)
            except RuntimeError as e:
                log.error("DSSP/SASA FAILED -- %s: %s", acc, e)
                failed_structures.append(acc)
                continue

            # ── 4. Pre-scan: check all sites for AA mismatch / missing residue ─
            #     If any site in a protein triggers a structural mismatch, the
            #     whole protein is skipped.  This avoids mixing valid and invalid
            #     rows for the same protein — the structure is likely a newer
            #     AlphaSync version than the one used for the alphasa database.
            protein_skip_reason: str | None = None
            for _, row in protein_rows.iterrows():
                site    = int(row["site"])
                ptm_str = str(row["ptm"])
                res_sel = u.select_atoms(f"resid {site}")
                if len(res_sel) == 0:
                    protein_skip_reason = (
                        f"residue {site} ({ptm_str}) not found in structure"
                    )
                    qc_events.append({"category": "missing_residue",
                                      "acc": acc, "site": site,
                                      "ptm": ptm_str, "value": float("nan")})
                    break
                if has_aa_col:
                    expected_aa = str(row["aa"]).upper()
                    observed_aa = gemmi.find_tabulated_residue(
                        res_sel[0].resname
                    ).one_letter_code.upper()
                    if observed_aa != expected_aa:
                        protein_skip_reason = (
                            f"AA mismatch at site {site} ({ptm_str}): "
                            f"expected {expected_aa}, structure has {observed_aa}"
                        )
                        qc_events.append({"category": "aa_mismatch",
                                          "acc": acc, "site": site,
                                          "ptm": ptm_str,
                                          "value": f"{expected_aa}→{observed_aa}"})
                        break

            if protein_skip_reason is not None:
                skipped_proteins[acc] = protein_skip_reason
                log.info(
                    "SKIP  %s — %s (structure likely newer than alphasa)",
                    acc, protein_skip_reason,
                )
                continue

            # ── 5. Compute depth for each PTM residue ─────────────────────────
            protein_results = []
            for _, row in protein_rows.iterrows():
                site    = int(row["site"])
                ptm_str = str(row["ptm"])

                sc_coords = sidechain_heavy_coords(u, site, acc)

                if sc_coords is None:
                    failed_residues.append((acc, site))
                    dist = float("nan")
                else:
                    atom_dists, _ = surf_tree.query(sc_coords, workers=-1)
                    dist = float(atom_dists.min())

                # ── QC checks ─────────────────────────────────────────────────
                if not np.isnan(dist):
                    if dist == QC_ZERO_DIST:
                        qc_events.append({"category": "zero_dist",
                                          "acc": acc, "site": site,
                                          "ptm": ptm_str, "value": dist})
                        # Expected: a sidechain atom of this "buried" residue
                        # (RSA ≤ 0.25) is itself part of the solvent-accessible
                        # surface — boundary residue, not an error.
                        log.debug(
                            "ZERO DIST  %s site %d (%s) — sidechain atom "
                            "is itself on the surface (boundary residue)",
                            acc, site, ptm_str,
                        )
                    elif dist < QC_NEAR_ZERO:
                        qc_events.append({"category": "near_zero_dist",
                                          "acc": acc, "site": site,
                                          "ptm": ptm_str, "value": dist})
                    elif dist > QC_LARGE_DIST:
                        qc_events.append({"category": "large_dist",
                                          "acc": acc, "site": site,
                                          "ptm": ptm_str, "value": dist})

                if site in surface_resids:
                    qc_events.append({"category": "dssp_alphasa_conflict",
                                      "acc": acc, "site": site,
                                      "ptm": ptm_str, "value": dist})
                    log.debug(
                        "DSSP/ALPHASA CONFLICT  %s site %d (%s) — "
                        "DSSP classifies this residue as surface-exposed "
                        "(relASA > 0.25) but alphasa called it buried",
                        acc, site, ptm_str,
                    )

                out_dict = {
                    "acc":              acc,
                    "site":             site,
                    "aa":               str(row.get("aa", "")),
                    "ptm":              ptm_str,
                    "source":           str(row.get("source", "")),
                    "ptmbin":           str(row.get("ptmbin", "Modified")),
                    "dis":              str(row.get("dis", "")),
                    "relasa":           float(row.get("relasa", np.nan)),
                    "plddt":            float(row.get("plddt", np.nan)),
                    "surf_dist_3d_ang": dist,
                }
                
                # min_pae might be missing/NaN
                mpa = row.get("min_pae", np.nan)
                out_dict["min_pae"] = np.nan if pd.isna(mpa) else float(mpa)
                
                protein_results.append(out_dict)

            # ── 6. Append this protein's rows to the CSV immediately ──────────
            chunk = pd.DataFrame(protein_results, columns=OUTPUT_COLS)
            chunk.to_csv(
                out_fh,
                index=False,
                header=(write_header and i == 0),
            )
            out_fh.flush()

            # ── 7. Per-protein success log ─────────────────────────────────────
            valid_dists = [r["surf_dist_3d_ang"] for r in protein_results
                           if not np.isnan(r["surf_dist_3d_ang"])]
            if valid_dists:
                log.info(
                    "OK  %s — %d site(s)  median=%.1f Å  min=%.1f Å  max=%.1f Å",
                    acc, len(valid_dists),
                    float(np.median(valid_dists)), min(valid_dists), max(valid_dists),
                )
            else:
                log.info("OK  %s — 0 sites with valid distances", acc)
            total_written += 1

    # ── Final summary ─────────────────────────────────────────────────────────
    log.info("Done. Wrote results for %d proteins to %s", total_written, out_path)
    out_df = pd.read_csv(out_path)
    valid  = out_df["surf_dist_3d_ang"].dropna()
    if len(valid) > 0:
        log.info(
            "surf_dist_3d_ang:  n=%d  median=%.1f Å  mean=%.1f Å  "
            "min=%.1f Å  max=%.1f Å",
            len(valid), valid.median(), valid.mean(), valid.min(), valid.max(),
        )

    if failed_structures:
        log.warning(
            "%d proteins had no usable structure: %s%s",
            len(failed_structures),
            ", ".join(failed_structures[:10]),
            " ..." if len(failed_structures) > 10 else "",
        )
    if failed_residues:
        log.warning(
            "%d (acc, site) pairs had no matching residue in the structure",
            len(failed_residues),
        )

    # ── Skipped proteins (structure/sequence version mismatch) ────────────────
    if skipped_proteins:
        log.warning(
            "%d proteins skipped due to AA mismatch or missing residues "
            "(structure likely newer than alphasa database):",
            len(skipped_proteins),
        )
        for skip_acc, reason in sorted(skipped_proteins.items()):
            log.warning("  %-12s  %s", skip_acc, reason)

    # ── QC report ─────────────────────────────────────────────────────────────
    report_qc(qc_events)

    zero_dist_count = sum(1 for e in qc_events if e["category"] == "zero_dist")
    if zero_dist_count > 0:
        log.info("")
        log.info("=" * 60)
        log.info("DIAGNOSTICS: Zero-Distance Boundary Atoms")
        log.info("Found %d residues with exactly 0.0 A distance.", zero_dist_count)
        log.info("These are surface-touching atoms (often the modifiable terminal atom)")
        log.info("on residues that are otherwise overall buried (RSA <= 0.25).")
        log.info("Their exact 'relasa' values are preserved in the output CSV.")
        log.info("=" * 60)

    # ── 7. Generate R Plots ───────────────────────────────────────────────────
    log.info("")
    log.info("Generating density plots via Rscript plot_surfdist_3d.R ...")
    try:
        subprocess.run(
            ["Rscript", "plot_surfdist_3d.R", str(out_path)],
            check=True
        )
    except Exception as e:
        log.error("Failed to execute Rscript: %s", e)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Compute 3D surface distance (Angstroms) for buried PTM residues "
            "from AlphaFold structures."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--rds",
        required=True,
        help=(
            "Path to tmp-qm.rds (from alphasa_relasa_vs_ptms.R). "
            "Required columns (exact names): "
            "acc, site, ptm, source, relasa, plddt, min_pae."
        ),
    )
    parser.add_argument(
        "--out",
        default="output-surfdist-3d.csv",
        help="Output CSV file.",
    )
    parser.add_argument(
        "--struct_dir",
        default="af_structures",
        help=(
            "Directory for caching CIF, PDB, and DSSP files. "
            "Created automatically if absent."
        ),
    )
    parser.add_argument(
        "--local_cif_dir",
        type=Path,
        default=DEFAULT_LOCAL_CIF_DIR,
        help=(
            "Directory containing local AlphaSync v0 CIF files "
            "(AF-{acc}-F1-model_v0.cif.gz). "
            "These take priority over all remote sources."
        ),
    )
    parser.add_argument(
        "--reset",
        action="store_true",
        help=(
            "Delete --struct_dir and --out, then exit. "
            "Use this to start completely fresh. "
            "Re-run without --reset afterwards to process everything from scratch."
        ),
    )
    parser.add_argument(
        "--max_structure_date",
        default=None,
        help=(
            "Maximum allowed structure file modification date (YYYY-MM-DD). "
            "If provided, local files newer than this will be ignored. "
            "Useful for ensuring structures match a specific database version."
        ),
    )
    args = parser.parse_args()
    main(args)

