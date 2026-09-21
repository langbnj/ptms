#!/usr/bin/env python3
"""
summarise_phosphosites.py
-------------------------
Collapses the per-contact table from find_coordinated_phosphoresidues.py into
one row per unique phosphosite, and reports the counts quoted in the response
to Reviewer #2, major point 6.

Definition of a site
    One UniProt accession and one position in its canonical sequence. Only
    contacts whose PDB-to-UniProt mapping was confirmed against the canonical
    sequence are used, so a site is never split in two by inconsistent
    deposited numbering. PRKACA pThr198, for instance, is deposited as
    residue 197 in some entries and 252 in others.

Burial and coordination together
    A site qualifies when one chain of one entry is both buried and fully
    coordinated at the same time. Taking the lowest accessibility from one
    entry and the contact count from another would overstate the case.

Reproducibility
    A phosphosite usually appears in more than one entry, and burial can
    depend on the crystal form or on what is bound. Three levels are therefore
    reported: qualifying in at least one entry, in at least half of them, and
    in every entry. The denominator is every entry in which the phosphorylated
    residue is resolved, including those where it has no Lys or Arg contact at
    all.

Input
-----
coordinated_phosphoresidues.tsv from find_coordinated_phosphoresidues.py.

Output
------
TSV with one row per site: accession, position, protein, residue type, the
lowest accessibility seen in any chain, the largest number of coordinating
side chains, the entry counts, and the best UniProt evidence for the
modification.

Usage
-----
    python summarise_phosphosites.py
    python summarise_phosphosites.py --tsv coordinated_phosphoresidues.tsv --out phosphosites.tsv

Dependencies
------------
    pip install pandas
"""

import argparse

import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

RELASA_CUTOFF = 0.25  # buried at or below this (Levy 2010)
MIN_CONTACTS  = 2     # fully coordinated at or above this many Lys/Arg side chains

INFILE  = "coordinated_phosphoresidues.tsv"
OUTFILE = "phosphosites.tsv"

# UniProt evidence for the modification itself, best first.
EVIDENCE_RANK = {"experimental": 4, "by_similarity": 3, "predicted": 2, "not_annotated": 1,
                 "unmapped": 0, "mapping_conflict": -1, "no_uniprot_data": -2}
SUPPORTED = ("experimental", "by_similarity")

HEADER = f"""\
# Unique phosphosites coordinated by Lys or Arg, collapsed from
# coordinated_phosphoresidues.tsv. Written by summarise_phosphosites.py.
#
# acc, site           UniProt accession and position in its canonical sequence
# min_relasa          lowest relative accessibility of the site in any chain where it is
#                     coordinated by Lys or Arg
# min_asa             the same in Å², before normalisation
# max_contacts        largest number of coordinating Lys or Arg side chains seen
# buried_and_coordinated  some chain of some entry is both buried (<= {RELASA_CUTOFF}) and fully
#                     coordinated (>= {MIN_CONTACTS} side chains) at the same time
# entries_total       deposited entries in which the site is resolved, coordinated or not,
#                     and entries_buried_and_coordinated how many of those qualify
# fraction_of_entries the proportion of entries that qualify
# evidence            best UniProt evidence for phosphorylation at this position
"""


# ---------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------

def best_evidence(values):
    """The strongest UniProt evidence level among several entries for a site."""
    return max(values, key=lambda v: EVIDENCE_RANK.get(v, -3))


def most_common(values):
    """The commonest protein name for a site, entries being inconsistent."""
    counts = values.mode()
    return counts.iat[0] if len(counts) else values.iat[0]


def collapse(rows):
    """One row per coordinated site, with the reproducibility counts.

    rows holds every resolved phosphoresidue, coordinated or not. A site is
    reported if it is coordinated in at least one entry, and entries_total
    counts every entry in which it is resolved, coordinated or not.
    """
    rows = rows.copy()
    rows["qualifies"] = ((rows.relasa_chain <= RELASA_CUTOFF)
                         & (rows.n_contacts >= MIN_CONTACTS))
    coordinated = rows[rows.n_contacts >= 1]

    chains = coordinated.groupby(["acc", "site", "pdb_id", "chain"]).agg(
        relasa=("relasa_chain", "min"),
        asa=("asa_chain", "min"),
        contacts=("n_contacts", "max"),
        qualifies=("qualifies", "any"),
    ).reset_index()

    sites = chains.groupby(["acc", "site"]).agg(
        min_relasa=("relasa", "min"),
        min_asa=("asa", "min"),
        max_contacts=("contacts", "max"),
        buried_and_coordinated=("qualifies", "any"),
    ).reset_index()

    entries = rows.groupby(["acc", "site", "pdb_id"]).agg(
        qualifies=("qualifies", "any")).reset_index()
    counts = entries.groupby(["acc", "site"]).agg(
        entries_total=("pdb_id", "nunique"),
        entries_buried_and_coordinated=("qualifies", "sum")).reset_index()

    annotation = coordinated.groupby(["acc", "site"]).agg(
        protein_name=("protein_name", most_common),
        phos_res=("phos_res", "first"),
        evidence=("site_evidence", best_evidence),
    ).reset_index()

    sites = sites.merge(counts, on=["acc", "site"]).merge(annotation, on=["acc", "site"])
    sites["fraction_of_entries"] = sites.entries_buried_and_coordinated / sites.entries_total
    return sites[["acc", "site", "protein_name", "phos_res", "min_relasa", "min_asa",
                  "max_contacts", "buried_and_coordinated", "entries_total",
                  "entries_buried_and_coordinated", "fraction_of_entries", "evidence"]]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--tsv", default=INFILE, help=f"contact table (default {INFILE})")
    parser.add_argument("--out", default=OUTFILE, help=f"output TSV (default {OUTFILE})")
    args = parser.parse_args()

    rows = pd.read_csv(args.tsv, sep="\t", comment="#")
    kept = rows[rows["include"] == True].copy()
    contacts, kept_contacts = rows[rows.n_contacts >= 1], kept[kept.n_contacts >= 1]
    print(f"Contacts in {args.tsv}:")
    print(f" >> {len(contacts)} total, in {contacts.pdb_id.nunique()} structures")
    print(f" >> {len(kept_contacts)} after filtering, in {kept_contacts.pdb_id.nunique()} structures")
    print(f" >> plus {int((kept.n_contacts == 0).sum())} phosphoresidues with no Lys/Arg contact, "
          f"used only to count the entries in which each site is resolved")

    sites = collapse(kept)
    qualifying = sites[sites.buried_and_coordinated]
    most = qualifying[qualifying.fraction_of_entries >= 0.5]
    every = qualifying[qualifying.fraction_of_entries == 1.0]

    print("\nCoordinated phosphosites:")
    print(f" >> {len(sites)} unique sites in {sites.acc.nunique()} proteins")
    print(f" >> {int(sites.evidence.isin(SUPPORTED).sum())} annotated as phosphorylated in UniProt")
    print(f" >> {int((sites.max_contacts >= MIN_CONTACTS).sum())} fully coordinated "
          f"(>= {MIN_CONTACTS} side chains)")
    print(f" >> {int((sites.min_relasa <= RELASA_CUTOFF).sum())} buried "
          f"(relative accessibility <= {RELASA_CUTOFF}) in some chain")

    print("\nBuried and fully coordinated in the same chain:")
    print(f" >> {len(qualifying):3d} sites in {qualifying.acc.nunique()} proteins, "
          f"in at least one entry")
    print(f" >> {len(most):3d} sites in {most.acc.nunique()} proteins, "
          f"in at least half the entries")
    print(f" >> {len(every):3d} sites in {every.acc.nunique()} proteins, in every entry")

    with open(args.out, "w") as fh:
        fh.write(HEADER)
    sites.sort_values(["min_relasa", "acc", "site"]).to_csv(
        args.out, sep="\t", index=False, mode="a")
    print(f"\n >> Written to {args.out}")


if __name__ == "__main__":
    main()
