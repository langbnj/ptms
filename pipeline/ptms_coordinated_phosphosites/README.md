# Coordinated phosphosites in the PDB

Scripts written for Reviewer #2, major point 6, which asks whether the buried
phosphosites coordinated by positively charged side chains that are predicted in
this work are also seen in experimentally determined structures.

They survey every human X-ray structure in the PDB at 2.0 Å resolution or better
that contains phosphoserine (SEP), phosphothreonine (TPO) or phosphotyrosine
(PTR), and collect the phosphoresidues that are coordinated by lysine or
arginine, measuring how buried each one is.

Burial and coordination are both assessed within one protein chain at a time, so
that the survey is comparable with the monomeric AlphaFold analysis used
elsewhere in this work. A residue buried in the isolated chain is buried whether
or not the protein has a binding partner, whereas a residue buried only by a
crystallographic partner is not. Accessibility is computed with `mkdssp` and
normalised against the Tien et al. (2013) empirical reference values, as in the
main analysis.

## Running

```bash
python run.py
```

runs the four steps in order:

```bash
python find_coordinated_phosphoresidues.py # the survey itself
python summarise_phosphosites.py           # one row per site, with the counts
python classify_kinase_sites.py            # kinase and activation loop labels
python build_supplementary_table.py        # Supplementary Table 6b
```

`phospho_reference_asa.py` derives the maximum accessible surface area of SEP,
TPO and PTR, for which no published values exist. Its results are already
written into `find_coordinated_phosphoresidues.py` as constants, so it only needs
running to reproduce them. It takes under a minute.

Nothing else needs preparing. Structures, metadata, residue-level SIFTS
mappings and modification annotations are downloaded as needed and cached under
`pdb_cache/`, including the DSSP results. The first run takes from about 15
minutes to over an hour depending on the connection, and a second run about
three minutes. The cache is roughly 750 MB.

## Reproducibility

The search is limited to structures released on or before 16 September 2026
(`RELEASE_CUTOFF` in `find_coordinated_phosphoresidues.py`), which returns 1236
entries. Repeating it returns the same entries as the PDB grows, unless any are
later made obsolete. Raise the date to bring the survey up to date. UniProt
annotations and SIFTS mappings are retrieved live and can change over time; the
cache under `pdb_cache/` keeps them fixed for re-runs.

Accessibility has been computed on two machines with different DSSP builds,
mkdssp 4.2.2 on x86-64 and 4.0.4 on ARM64, with identical results. The values
in Supplementary Table 6b also agree with the DSSP server at PDB-REDO
(https://pdb-redo.eu/dssp, DSSP 4.6.1) run on the same chains.

Expected counts, printed by `summarise_phosphosites.py`:

| | |
| --- | --- |
| Unique coordinated phosphosites | 101, in 71 proteins |
| Fully coordinated (two or more Lys/Arg) | 69 |
| Buried within the chain | 29 |
| Buried and fully coordinated in the same chain | 23, in 19 proteins |
| of which in every entry in which the site is resolved | 9 |
| in at least half of those entries | 16 |

## Files

| Script | What it does |
| --- | --- |
| `run.py` | Runs the four steps below in order. |
| `phospho_reference_asa.py` | Derives reference accessibilities for SEP, TPO and PTR by repeating the Gly-X-Gly enumeration of Tien et al. |
| `find_coordinated_phosphoresidues.py` | Searches the PDB and writes one row per phosphoresidue-to-Lys/Arg contact. |
| `summarise_phosphosites.py` | Collapses those contacts to one row per unique phosphosite and reports the counts quoted in the response. |
| `classify_kinase_sites.py` | Labels each site as kinase or not, and locates it in the kinase domain, using `Human-PK-alignment.fasta`. |
| `build_supplementary_table.py` | Builds Supplementary Table 6b as a TSV and a formatted workbook, with the matching paragraph of the legend. Given the curated literature table with `--table6a`, also assembles the two-tab Supplementary Table 6. |

| Input | Contents |
| --- | --- |
| `Human-PK-alignment.fasta` | The Kincore structure-based alignment of 497 human protein kinase domains (Modi & Dunbrack, 2019), from https://dunbrack.fccc.edu/kincore/alignment. |

| Output | Contents |
| --- | --- |
| `coordinated_phosphoresidues.tsv` | One row per contact, with accessibility, distances and the UniProt mapping; phosphoresidues with no contact get one row, so that every entry in which a site is resolved is counted. |
| `failed_structures.txt` | Structures that could not be processed, with the reason. |
| `phosphosites.tsv` | One row per unique phosphosite. |
| `phosphosites_classified.tsv` | The same, with the kinase labels added. |
| `Supplementary_Table_S6b.tsv`, `.xlsx` | The supplementary table. |
| `Supplementary_Table_S6b_legend.txt` | Its legend and the notes that define its columns, with every number filled in from the data. The same notes are written beneath the table in the workbook and as a commented header in the TSV. |
| `phospho_reference_asa.tsv` | The derived reference values. |

Each output file carries a commented header explaining its columns.

## Criteria

| | |
| --- | --- |
| Structures | human, X-ray, 2.0 Å resolution or better, released by 16 September 2026 |
| Buried | relative solvent accessibility ≤ 0.25 (Levy 2010), DSSP on the isolated chain |
| Ionic contact | ≤ 4.0 Å between the closest atoms of the phosphate group and the Lys or Arg charged group (Supplementary Table 4) |
| Fully coordinated | two or more such side chains, in the same chain |
| Mapping | each residue mapped to UniProt individually through SIFTS, and accepted only where the canonical sequence has the matching Ser, Thr or Tyr |
| Ubiquitin | encoded by UBB, UBC, UBA52 and RPS27A with an identical sequence; its sites are numbered within the 76-residue ubiquitin unit and counted once, under UBC (P0CG48) |
| Excluded | chains under 30 residues, which are co-crystallised substrate peptides, chains whose phosphoresidue does not belong to a human UniProt entry, and phosphoresidues that SIFTS does not map to UniProt, such as those in expression tags |

## Dependencies

```bash
pip install biopython numpy openpyxl pandas requests tqdm
```

and `mkdssp` on the `PATH` (https://github.com/PDB-REDO/dssp). Versions 4.0.4
and 4.2.2 have been tested and give identical results; see Reproducibility
above.

## References

Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO (2013) Maximum allowed
solvent accessibilites of residues in proteins. *PLoS One* 8(11):e80635.
PMID 24278298.

Levy ED (2010) A simple definition of structural regions in proteins and its use
in analyzing interface evolution. *J Mol Biol* 403:660-670. PMID 20868694.

Modi V, Dunbrack RL Jr (2019) A structurally-validated multiple sequence
alignment of 497 human protein kinase domains. *Sci Rep* 9:19790. PMID 31875044.
