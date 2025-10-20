# ptms
Code for manuscript:

"Protein dynamics at different timescales unlock access to hidden post-translational modification sites"

Preprint: https://www.biorxiv.org/content/10.1101/2025.06.25.661537v1

## Contents

This project also includes code and uses SQL tables from AlphaSync (https://alphasync.stjude.org, code at https://github.com/langbnj/alphasync, preprint at https://www.biorxiv.org/content/10.1101/2025.03.12.642845v2).

**analyses**: R scripts used to perform statistical analyses and produce figures.

**include**: Utility functions.

**mysql_tables**: MySQL table structures (CREATE statements).

**pipeline**: Scripts used to extract and analyse data from data tables (SQL).

**scripts**: Additional utility scripts.

**update**: Scripts used to fill data tables (SQL).

For "pipeline" and "update", first run "download.py" (or "download.pl") to download datasets, followed by "run.py" (or "run.pl").

## Update sequence

Create MySQL tables using "mysql_tables/create_statements.sql". Then, fill MySQL tables in this sequence:

- update/uniprot (UniProt version 2022_04)
- update/ensembl (Ensembl version 108)
- update/tax (NCBI taxonomy tables)
- update/compara (Ensembl Compara tables)
- update/comparanopara (Ensembl Compara tables, paralog-filtered)
- update/evolutionary_rate_rate4site (Rate4Site evolutionary rates)
- update/evolutionary_rate_lichtarge (real-valued Evolutionary Trace score, rvET)
- update/evolutionary_rate_capra (Capra et al. conservation scores, Jensen-Shannon Divergence)
- update/evolutionary_rate_sporadic_protein_removal (Remove a handful of proteins not fully covered by all conservation scoring methods)
- update/evolutionary_rate_normalize (Normalize conservation scores)

...followed by the other tables.


## Note

Please don't hesitate to contact me/us for questions and assistance with these scripts (<benjamin.lang@cantab.net>). They have grown historically as these analyses grew far more complex than originally envisioned, though we have done our best to comment them well enough to follow along and use them. Thank you!
