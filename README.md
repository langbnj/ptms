# ptms
Code for manuscript:

**"Protein dynamics at different timescales unlock access to hidden post-translational modification sites"**

**Article preprint: https://www.biorxiv.org/content/10.1101/2025.06.25.661537v1**

## AlphaSync

This project also uses code and SQL tables from AlphaSync (https://alphasync.stjude.org). 

AlphaSync code: https://github.com/langbnj/alphasync.

AlphaSync article: https://www.nature.com/articles/s41594-025-01719-x.

## Contents

**analyses**: R scripts used to perform statistical analyses and produce figures.

**include**: Utility functions.

**mysql_tables**: MySQL table structures (CREATE statements).

**pipeline**: Scripts used to extract and analyse data from data tables (SQL).

**scripts**: Additional utility scripts.

**update**: Scripts used to fill data tables (SQL).

For "pipeline" and "update", first run "download.py" (or "download.pl") to download datasets, followed by "run.py" (or "run.pl").

## Install dependencies

### Python
`pip install biopython requests sqlalchemy ipdb natsort numpy pandas pymysql scipy tabulate tqdm`
### Perl
`cpan -i Bundle::CPAN Term::ReadKey Term::ReadLine::Perl CPAN::DistnameInfo XML::Parser Statistics::Descriptive Statistics::R Text::CSV List::Compare Sort::Key::Natural Carp::Always Math::Round Statistics::Robust::Scale List::Vectorize Statistics::Multtest Statistics::RankCorrelation DBD::mysql`
### SQL
A MySQL or MariaDB server.

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

Please don't hesitate to contact us for questions and assistance with these scripts (<benjamin.lang@cantab.net>, <madan.babu@stjude.org>). They have grown historically as these analyses grew far more complex than originally envisioned, though we have done our best to comment them well enough to follow along and use them. Thank you!
