#!/usr/bin/env python3
"""
classify_kinase_activation_loops.py

Classifies phosphosites from coordinated_phosphoresidues_all_rsa.tsv as:
  - kinase / non-kinase
  - activation loop / other kinase region / outside kinase domain / not in alignment

Uses the Kincore/Modi & Dunbrack 2019 structure-based multiple sequence alignment
of 497 human protein kinase domains to define activation loop boundaries.

Reference:
  Modi & Dunbrack, Sci Rep 2019, 9, 19790
  https://doi.org/10.1038/s41598-019-56499-4

Activation loop is defined as the union of aligned blocks ALN and ALC
(columns 1331-1351 and 1904-1920 in the alignment, 1-indexed),
plus biologically confirmed unaligned inserts between them (columns 1352-1903)
for dual-phosphorylation motifs (JAK1/2, SYK).

Usage:
  python classify_kinase_activation_loops.py \
      --tsv coordinated_phosphoresidues_all_rsa.tsv \
      --alignment Human-PK-alignment.fasta \
      [--canonical-only] \
      [--output results_with_kinase_classification.tsv]

Obtaining the alignment file:
  The Human-PK-alignment.fasta file is available for download from the Kincore
  Alignment page at https://dunbrack.fccc.edu/kincore/alignment
  (select "FASTA-formatted file"). The Clustal format (.aln) is also available
  from the same page.

  Reference:
    Modi & Dunbrack, Sci Rep 2019, 9, 19790
    https://doi.org/10.1038/s41598-019-56499-4

  Kincore web resource:
    Modi & Dunbrack, Nucleic Acids Res 2022, 50, D654-D664
    https://doi.org/10.1093/nar/gkab920
"""

import argparse
import re
import sys
import pandas as pd


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Activation loop column ranges in the Kincore alignment (1-indexed)
ALN_COLS = (1331, 1351)   # N-terminal activation loop (contains DFG motif)
ALC_COLS = (1904, 1920)   # C-terminal activation loop (contains APE region)

# Unaligned insert between ALN and ALC: columns 1352-1903
# This region is too variable to align structurally but biologically contains
# the middle of the activation loop. We include specific well-confirmed cases
# where the phosphosite falls here (dual-Tyr activation loop motifs).
UNALIGNED_AL_COLS = (1352, 1903)

# Sites in the unaligned insert that are confirmed activation loop by biology.
# These are the second Tyr of dual-Tyr activation loop motifs where the first
# Tyr aligns to ALN col 1351 and the second falls one position outside in the
# unaligned insert. Add (uniprot_id, resnum) tuples here if new cases arise.
CONFIRMED_UNALIGNED_AL_SITES = {
    ("Q9NWZ3", 345),   # IRAK4 pThr345 (activation loop, col 1356)
    ("O60674", 1008),  # JAK2 pTyr1008 (dual activation loop Tyr, col 1354)
    ("P23458", 1035),  # JAK1 pTyr1035 (dual activation loop Tyr, col 1354)
    ("P43405", 526),   # SYK pTyr526   (dual activation loop Tyr, col 1354)
}

# Known SIFTS mapping errors where the TSV maps a phospho-residue to a
# UniProt position with a chemically incompatible amino acid (e.g. pThr→Glu).
# These are excluded from activation loop classification.
KNOWN_SIFTS_ERRORS = {
    ("Q16512", 780),   # PKN1: pThr780 maps to UniProt Glu780 — impossible
}

# Regex patterns used to identify kinase entries by protein_name.
# We use a conservative set of patterns that match curated kinase names
# without over-matching (e.g. "kinase-like", "kinase inhibitor").
KINASE_PATTERNS = [
    r"\bkinase\b",
    r"\btyrosine-protein kinase\b",
    r"\bserine[/\-]threonine[- ]protein kinase\b",
]
KINASE_REGEX = re.compile("|".join(KINASE_PATTERNS), re.IGNORECASE)


# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------

def parse_fasta(path: str) -> dict:
    """Parse a FASTA file and return {header: sequence}."""
    seqs = {}
    name = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                name = line[1:].strip()
                seqs[name] = ""
            elif name is not None:
                seqs[name] += line
    return seqs


# ---------------------------------------------------------------------------
# Build UniProt → alignment lookup
# ---------------------------------------------------------------------------

def build_uniprot_lookup(aligned: dict) -> dict:
    """
    Parse aligned FASTA headers to extract UniProt accession and domain
    start/end positions.

    Header format: FAMILY_GENE/start-end PROTEIN_ID GENE UNIPROT
    e.g. CMGC_CDK2/4-286 CDK2_HUMAN CDK2 P24941

    Returns {uniprot_id: (aligned_sequence, domain_start, domain_end)}.
    """
    lookup = {}
    # UniProt accession pattern: letter, digit, 3 alphanums, digit
    acc_re = re.compile(r"^[A-Z][0-9][A-Z0-9]{3}[0-9]$")
    range_re = re.compile(r"/(\d+)-(\d+)")

    for header, seq in aligned.items():
        uniprot = None
        for token in header.split():
            if acc_re.match(token):
                uniprot = token
                break
        m = range_re.search(header)
        if uniprot and m:
            dom_start = int(m.group(1))
            dom_end = int(m.group(2))
            lookup[uniprot] = (seq, dom_start, dom_end)

    return lookup


# ---------------------------------------------------------------------------
# Column lookup
# ---------------------------------------------------------------------------

def get_alignment_col(aligned_seq: str, dom_start: int, dom_end: int,
                      target_resnum: int):
    """
    Return (column_1indexed, amino_acid_char) for the given protein residue
    number, or (None, None) if the residue is outside the domain.

    Gap characters ('-', '.') are skipped; residue counting uses both
    upper- and lower-case letters (upper = aligned block, lower = unaligned
    insert in Kincore format).
    """
    if target_resnum < dom_start or target_resnum > dom_end:
        return None, None

    res_count = dom_start - 1
    for col_idx, char in enumerate(aligned_seq):
        if char not in "-.":
            res_count += 1
            if res_count == target_resnum:
                return col_idx + 1, char  # 1-indexed column

    return None, None


# ---------------------------------------------------------------------------
# Activation loop classification
# ---------------------------------------------------------------------------

def classify_site(uniprot_id: str, resnum: int,
                  uniprot_lookup: dict) -> str:
    """
    Classify a single (uniprot_id, resnum) site.

    Returns one of:
      'activation_loop_ALN'     — falls in ALN aligned block
      'activation_loop_ALC'     — falls in ALC aligned block
      'activation_loop_insert'  — falls in unaligned insert, confirmed AL
      'not_activation_loop'     — within kinase domain but outside AL
      'outside_domain'          — outside the annotated kinase domain
      'not_in_alignment'        — UniProt ID not found in the alignment
      'sifts_error'             — known SIFTS mapping error, excluded
    """
    if (uniprot_id, resnum) in KNOWN_SIFTS_ERRORS:
        return "sifts_error"

    if uniprot_id not in uniprot_lookup:
        return "not_in_alignment"

    seq, dom_start, dom_end = uniprot_lookup[uniprot_id]
    col, char = get_alignment_col(seq, dom_start, dom_end, resnum)

    if col is None:
        return "outside_domain"

    if ALN_COLS[0] <= col <= ALN_COLS[1]:
        return "activation_loop_ALN"
    if ALC_COLS[0] <= col <= ALC_COLS[1]:
        return "activation_loop_ALC"
    if (UNALIGNED_AL_COLS[0] <= col <= UNALIGNED_AL_COLS[1]
            and (uniprot_id, resnum) in CONFIRMED_UNALIGNED_AL_SITES):
        return "activation_loop_insert"

    return "not_activation_loop"


def is_activation_loop(classification: str) -> bool:
    return classification.startswith("activation_loop")


# ---------------------------------------------------------------------------
# Kinase detection
# ---------------------------------------------------------------------------

def is_kinase(protein_name: str) -> bool:
    """Return True if protein_name matches kinase patterns."""
    if pd.isna(protein_name):
        return False
    return bool(KINASE_REGEX.search(protein_name))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Classify phosphosites as kinase/non-kinase and "
                    "activation loop / other, using the Kincore alignment."
    )
    parser.add_argument(
        "--tsv", default="coordinated_phosphoresidues_all_rsa.tsv",
        help="Path to coordinated_phosphoresidues_all_rsa.tsv "
             "(default: %(default)s)"
    )
    parser.add_argument(
        "--alignment", default="Human-PK-alignment.fasta",
        help="Path to Human-PK-alignment.fasta (Kincore/Modi & Dunbrack 2019) "
             "(default: %(default)s)"
    )
    parser.add_argument(
        "--canonical-only", action="store_true", default=True,
        help="Restrict analysis to is_canonical_human == True rows "
             "(default: True; use --no-canonical-only to disable)"
    )
    parser.add_argument(
        "--no-canonical-only", dest="canonical_only", action="store_false",
        help="Include all rows, not just canonical human"
    )
    parser.add_argument(
        "--output", default="results_with_kinase_classification.tsv",
        help="Output TSV path (default: %(default)s)"
    )
    args = parser.parse_args()

    # ── Load data ──────────────────────────────────────────────────────────
    print(f"Loading TSV: {args.tsv}", file=sys.stderr)
    df = pd.read_csv(args.tsv, sep="\t")

    if args.canonical_only:
        df = df[df["is_canonical_human"] == True].copy()
        print(f"  Restricted to canonical human: {len(df)} rows", file=sys.stderr)

    print(f"Loading alignment: {args.alignment}", file=sys.stderr)
    aligned = parse_fasta(args.alignment)
    uniprot_lookup = build_uniprot_lookup(aligned)
    print(f"  Kinases in alignment: {len(uniprot_lookup)}", file=sys.stderr)

    # ── Classify at unique-site level ──────────────────────────────────────
    # Deduplicate to (uniprot_id, phos_resnum) — one row per unique site.
    # Use uniprot_phos_resnum where available (canonical UniProt position),
    # falling back to phos_resnum (PDB deposited number).
    df["_site_resnum"] = df["uniprot_phos_resnum"].where(
        df["uniprot_phos_resnum"].notna(), df["phos_resnum"]
    ).astype(int)

    sites = (
        df.groupby(["uniprot_id", "_site_resnum"])
        .agg(
            protein_name=("protein_name", "first"),
            phos_resname=("phos_resname", "first"),
            phos_resnum_pdb=("phos_resnum", "first"),
            uniprot_phos_resnum=("uniprot_phos_resnum", "first"),
            min_rsa=("rsa_0to1_lower_is_more_buried", "min"),
            max_rsa=("rsa_0to1_lower_is_more_buried", "max"),
            max_ncoord=("n_coord_residues_higher_is_more_coordinated", "max"),
            buried_in_any=("buried_in_monomer", lambda x: any(x == True)),
            surface_in_any=("buried_in_monomer", lambda x: any(x == False)),
            n_structures=("pdb_id", "nunique"),
            evidence=("uniprot_phos_evidence", "first"),
        )
        .reset_index()
        .rename(columns={"_site_resnum": "site_resnum"})
    )

    # Add kinase flag
    sites["is_kinase"] = sites["protein_name"].apply(is_kinase)

    # Classify activation loop for kinase sites only
    al_classifications = []
    for _, row in sites.iterrows():
        if not row["is_kinase"]:
            al_classifications.append("not_a_kinase")
        else:
            cl = classify_site(
                row["uniprot_id"],
                int(row["site_resnum"]),
                uniprot_lookup
            )
            al_classifications.append(cl)

    sites["al_classification"] = al_classifications
    sites["is_activation_loop"] = sites["al_classification"].apply(
        lambda c: is_activation_loop(c)
    )

    # ── Summary statistics ─────────────────────────────────────────────────
    total        = len(sites)
    n_kinase     = int(sites["is_kinase"].sum())
    n_non_kinase = total - n_kinase
    n_al         = int(sites["is_activation_loop"].sum())
    n_kinase_not_al = n_kinase - n_al

    buried           = sites[sites["buried_in_any"]]
    n_buried         = len(buried)
    n_buried_kinase  = int(buried["is_kinase"].sum())
    n_buried_al      = int(buried["is_activation_loop"].sum())
    bur_fc           = buried[buried["max_ncoord"] >= 2]
    bur_sc           = buried[buried["max_ncoord"] == 1]

    surf_any         = sites[sites["surface_in_any"]]   # surface in ≥1 struct
    surf_only        = sites[sites["surface_in_any"] & ~sites["buried_in_any"]]
    n_surf_any       = len(surf_any)
    n_surf_kinase    = int(surf_any["is_kinase"].sum())
    n_surf_al        = int(surf_any["is_activation_loop"].sum())
    surf_fc          = surf_any[surf_any["max_ncoord"] >= 2]
    surf_sc          = surf_any[surf_any["max_ncoord"] == 1]

    all_fc           = sites[sites["max_ncoord"] >= 2]
    all_fc_surf_only = all_fc[all_fc["surface_in_any"] & ~all_fc["buried_in_any"]]

    print("\n=== CLASSIFICATION SUMMARY ===")
    print(f"Total unique canonical phosphosites:  {total}")
    print(f"  Kinase sites:                       {n_kinase} ({100*n_kinase/total:.0f}%)")
    print(f"    of which activation loop:          {n_al} ({100*n_al/n_kinase:.0f}% of kinase)")
    print(f"    of which other kinase region:      {n_kinase_not_al}")
    print(f"  Non-kinase sites:                   {n_non_kinase}")
    print(f"\nActivation loop as % of ALL sites:    {100*n_al/total:.0f}%")
    print(f"\nBuried sites (RSA ≤ 0.25 in ≥1 structure): {n_buried}")
    print(f"  Kinase activation loop:             {n_buried_al} ({100*n_buried_al/n_buried:.0f}%)")
    print(f"  Other kinase:                       {n_buried_kinase - n_buried_al}")
    print(f"  Non-kinase:                         {n_buried - n_buried_kinase}")
    print(f"  Fully coordinated (≥2 K/R):         {len(bur_fc)}")
    print(f"  Singly coordinated (=1 K/R):        {len(bur_sc)}")
    print(f"  Check: {len(bur_fc)}+{len(bur_sc)} = {len(bur_fc)+len(bur_sc)} "
          f"(should equal {n_buried})")
    print(f"\nSurface sites (RSA > 0.25 in ≥1 structure): {n_surf_any}")
    print(f"  Kinase sites:                       {n_surf_kinase}")
    print(f"    of which activation loop:          {n_surf_al} "
          f"({100*n_surf_al//n_surf_kinase if n_surf_kinase else 0}% of surface kinase)")
    print(f"  Fully coordinated (≥2 K/R):         {len(surf_fc)}")
    print(f"  Singly coordinated (=1 K/R):        {len(surf_sc)}")
    print(f"  Check: {len(surf_fc)}+{len(surf_sc)} = {len(surf_fc)+len(surf_sc)} "
          f"(should equal {n_surf_any})")
    print(f"\nAll fully coordinated (≥2 K/R, any burial): {len(all_fc)}")
    print(f"  Buried (any struct):                {len(bur_fc)}")
    print(f"  Surface-only (never buried):        {len(all_fc_surf_only)}")
    print(f"  Check: {len(bur_fc)}+{len(all_fc_surf_only)} = "
          f"{len(bur_fc)+len(all_fc_surf_only)} (should equal {len(all_fc)})")

    # AL classification breakdown for kinase sites
    print("\n=== ACTIVATION LOOP CLASSIFICATION DETAIL (kinase sites) ===")
    kinase_sites = sites[sites["is_kinase"]]
    for cl, grp in kinase_sites.groupby("al_classification"):
        print(f"  {cl}: {len(grp)}")

    # List non-kinase buried sites (most relevant for the paper)
    print("\n=== NON-KINASE BURIED SITES ===")
    nk_buried = sites[sites["buried_in_any"] & ~sites["is_kinase"]].sort_values("min_rsa")
    for _, r in nk_buried.iterrows():
        print(f"  {r['uniprot_id']}  pos {r['site_resnum']}  "
              f"RSA_min={r['min_rsa']:.3f}  {r['phos_resname']}  "
              f"{r['evidence']}  {r['protein_name'][:55]}")

    # ── Output ────────────────────────────────────────────────────────────
    if args.output:
        out_cols = [
            "uniprot_id", "site_resnum", "phos_resnum_pdb", "uniprot_phos_resnum",
            "protein_name", "phos_resname", "min_rsa", "max_rsa",
            "buried_in_any", "surface_in_any", "n_structures",
            "evidence", "is_kinase", "al_classification", "is_activation_loop",
        ]
        sites[out_cols].to_csv(args.output, sep="\t", index=False)
        print(f"\nFull site-level results written to: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
