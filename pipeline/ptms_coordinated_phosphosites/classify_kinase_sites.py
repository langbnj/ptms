#!/usr/bin/env python3
"""
classify_kinase_sites.py
------------------------
Labels each phosphosite from summarise_phosphosites.py as belonging to a
protein kinase or not and, for kinases, says where in the kinase domain it
lies. Supersedes classify_kinase_activation_loops.py.

This is what allows the response to Reviewer #2, major point 6 to say how many
of the buried, coordinated phosphosites are activation-loop sites in kinases
and how many are not.

Why an alignment rather than the protein name
    Kinase membership is taken from the Kincore structure-based alignment of
    human protein kinase domains, not from the PDB or UniProt protein name.
    Names are inconsistent, and matching on the word "kinase" also picks up
    things like "cyclin-dependent kinase inhibitor", which is not a kinase.

Where the activation loop is
    The activation segment runs from the DFG motif to the APE motif. In this
    alignment the DFG aspartate sits in column 1338 (Asp in 469 of the 497
    domains) and the APE glutamate in column 1916 (Glu in 461 of them), so
    those two columns define the segment for every kinase at once.

    Within it, two blocks are structurally alignable and are annotated ALN and
    ALC in the alignment's own annotation row: columns 1338 to 1351 from the
    DFG motif, and columns 1904 to 1916 up to the APE motif. Between them lies
    an insert of variable length that cannot be aligned across the family.
    Residues there are still between the two motifs in their own sequence, so
    they are part of the activation loop; this is where the second site of a
    dual-phosphorylation pair falls, as in JAK1 pTyr1035, JAK2 pTyr1008,
    JAK3 pTyr981, SYK pTyr526 and IRAK4 pThr345.

    Anchoring on the motifs rather than on the extent of the two aligned
    blocks matters at both ends. Columns 1331 to 1337 precede the DFG motif
    and belong to the preceding strand, which is where DYRK3 Ser350 lies, and
    columns 1917 to 1920 follow the APE motif and belong to the P+1 loop.

Input
-----
phosphosites.tsv from summarise_phosphosites.py, and the Kincore alignment
Human-PK-alignment.fasta (https://dunbrack.fccc.edu/kincore/alignment).
Modi V, Dunbrack RL Jr (2019) Sci Rep 9:19790, PMID 31875044.

Output
------
The input table with three columns added: is_kinase, region and
is_activation_loop.

Usage
-----
    python classify_kinase_sites.py
    python classify_kinase_sites.py --sites phosphosites.tsv --out phosphosites_classified.tsv

Dependencies
------------
    pip install pandas
"""

import argparse
import re
from collections import defaultdict

import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Alignment columns, counting from 1. See "Where the activation loop is" above.
DFG_COLUMN = 1338         # DFG aspartate, start of the activation segment
APE_COLUMN = 1916         # APE glutamate, end of the activation segment
ALN_END    = 1351         # last alignable column of the N-terminal block
ALC_START  = 1904         # first alignable column of the C-terminal block

SITES_FILE     = "phosphosites.tsv"
ALIGNMENT_FILE = "Human-PK-alignment.fasta"
OUTFILE        = "phosphosites_classified.tsv"

# UniProt accession, in both the six- and the ten-character form.
ACCESSION    = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]"
                          r"|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")
DOMAIN_RANGE = re.compile(r"/(\d+)-(\d+)")
GAP          = "-."

HEADER = """\
# Coordinated phosphosites with their kinase assignment, from phosphosites.tsv.
# Written by classify_kinase_sites.py.
#
# is_kinase           the protein has a domain in the Kincore alignment of human
#                     protein kinase domains
# region              activation_loop_N, activation_loop_insert or activation_loop_C
#                     between the DFG and APE motifs; kinase_domain_other elsewhere
#                     in the domain; outside_kinase_domain for another part of a
#                     kinase; not_a_kinase otherwise
# is_activation_loop  region is one of the three activation loop segments
"""


# ---------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------

def read_alignment(path):
    """{accession: [(aligned sequence, domain start, domain end)]} from the Kincore FASTA.

    A few kinases contribute two domains, such as the pseudokinase and kinase
    domains of the JAK family, so the domains are kept as a list rather than
    one per accession. The annotation row that labels the alignment blocks is
    skipped.
    """
    sequences, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                name = line[1:].strip()
                sequences[name] = ""
            elif name:
                sequences[name] += line

    domains = defaultdict(list)
    for header, sequence in sequences.items():
        if "ANNOTATION" in header:
            continue
        acc = next((token for token in header.split() if ACCESSION.match(token)), None)
        bounds = DOMAIN_RANGE.search(header)
        if acc and bounds:
            domains[acc].append((sequence, int(bounds.group(1)), int(bounds.group(2))))
    return dict(domains)


def alignment_column(sequence, start, end, resnum):
    """Column number, counting from 1, of a residue within an aligned domain."""
    if not start <= resnum <= end:
        return None
    position = start - 1
    for i, char in enumerate(sequence):
        if char not in GAP:
            position += 1
            if position == resnum:
                return i + 1
    return None


def classify(acc, site, domains):
    """Where in the kinase domain a site lies, by alignment column."""
    if acc not in domains:
        return "not_in_alignment"
    column = next((c for c in (alignment_column(sequence, start, end, site)
                               for sequence, start, end in domains[acc]) if c), None)
    if column is None:
        return "outside_kinase_domain"
    if not DFG_COLUMN <= column <= APE_COLUMN:
        return "kinase_domain_other"
    if column <= ALN_END:
        return "activation_loop_N"
    if column < ALC_START:
        return "activation_loop_insert"
    return "activation_loop_C"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--sites", default=SITES_FILE, help=f"site table (default {SITES_FILE})")
    parser.add_argument("--alignment", default=ALIGNMENT_FILE,
                        help=f"Kincore alignment (default {ALIGNMENT_FILE})")
    parser.add_argument("--out", default=OUTFILE, help=f"output TSV (default {OUTFILE})")
    args = parser.parse_args()

    sites = pd.read_csv(args.sites, sep="\t", comment="#")
    domains = read_alignment(args.alignment)
    print(f"Kincore alignment:\n >> {sum(len(d) for d in domains.values())} kinase domains "
          f"in {len(domains)} human proteins")

    sites["is_kinase"] = sites.acc.isin(domains)
    sites["region"] = [classify(row.acc, int(row.site), domains) if row.acc in domains
                       else "not_a_kinase" for row in sites.itertuples()]
    sites["is_activation_loop"] = sites.region.str.startswith("activation_loop")

    subsets = [
        ("All coordinated sites", sites),
        ("Buried and fully coordinated in at least one entry",
         sites[sites.buried_and_coordinated]),
        ("Buried and fully coordinated in every entry",
         sites[sites.buried_and_coordinated & (sites.fraction_of_entries == 1.0)]),
    ]
    for label, subset in subsets:
        total = len(subset)
        if not total:
            continue
        kinases = int(subset.is_kinase.sum())
        loops = int(subset.is_activation_loop.sum())
        print(f"\n{label}: {total}")
        print(f" >> {kinases} in protein kinases ({100 * kinases / total:.0f}%)")
        print(f" >> {loops} at activation loops ({100 * loops / total:.0f}% of all"
              + (f", {100 * loops / kinases:.0f}% of kinases)" if kinases else ")"))
        for region, count in subset.region.value_counts().items():
            print(f"      {region:<24} {count}")

    with open(args.out, "w") as fh:
        fh.write(HEADER)
    sites.to_csv(args.out, sep="\t", index=False, mode="a")
    print(f"\n >> Written to {args.out}")


if __name__ == "__main__":
    main()
