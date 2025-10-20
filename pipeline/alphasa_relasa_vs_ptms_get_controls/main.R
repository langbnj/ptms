#!/usr/bin/env Rscript --vanilla

library("RMySQL")
library("tibble")

if (exists("superreservedcon"))
{
  try(dbDisconnect(superreservedcon), silent=T)
  rm(superreservedcon)
}

# Insert server name here:
superreservedcon = dbConnect(dbDriver("MySQL"), host="[server name]", dbname="blang")
Query <- function(query)
{
  a = as_tibble(dbGetQuery(superreservedcon, statement=query))
  return(a)
}

# Load packages
suppressPackageStartupMessages(library("tidyverse"))
# suppressPackageStartupMessages(library("scales"))
# suppressPackageStartupMessages(library("magrittr"))
# suppressPackageStartupMessages(library("ggrepel"))
# suppressPackageStartupMessages(library("measurements"))
# suppressPackageStartupMessages(library("ggbeeswarm"))
# suppressPackageStartupMessages(library("tictoc"))

# Convenience function like Python's f-strings
f <- stringr::str_glue

# myquery <- "SELECT 2+2"
# Incorrect version with alphaccs instead of accs:
# myquery <- "SELECT a.acc, a.site, a.aa, '0' AS ptm, '0' AS source, '0' AS subset, '0' AS scale, a.relasa, a.relasa10, a.dis, a.dis10, a.plddt, a.plddt10, a.membrane, IFNULL(COUNT(c.pae), 0) AS contacts, MIN(c.pae) AS min_pae, AVG(c.pae) AS avg_pae, MAX(c.pae) AS max_pae FROM unimod_control m, uniseq s, alphamap am, alphasa a LEFT OUTER JOIN alphacon c ON a.acc=c.acc AND a.afdb=c.afdb AND (a.site=c.site1 OR a.site=c.site2) WHERE m.species='human' AND m.acc=s.acc AND s.type='UniProt' AND LENGTH(s.seq)>=16 AND m.acc=am.value AND am.type='uniprot' AND am.version='2022_04' AND am.best=1 AND am.map=a.acc AND am.afdb=a.afdb AND m.site=a.site AND a.membrane IS NULL GROUP BY a.acc, a.site, a.aa"
# Corrected version:
myquery <- "SELECT m.acc, m.site, m.aa, '0' AS ptm, '0' AS source, '0' AS subset, '0' AS scale, a.relasa, a.relasa10, a.dis, a.dis10, a.plddt, a.plddt10, a.membrane, IFNULL(COUNT(DISTINCT c.site1, c.site2), 0) AS contacted_residues, MIN(c.pae) AS min_pae, AVG(c.pae) AS avg_pae, MAX(c.pae) AS max_pae FROM unimod_control m, uniseq s, alphamap am, alphasa a LEFT OUTER JOIN alphacon c ON a.acc=c.acc AND a.afdb=c.afdb AND (a.site=c.site1 OR a.site=c.site2) WHERE m.species='human' AND m.acc=s.acc AND s.type='UniProt' AND LENGTH(s.seq)>=16 AND m.acc=am.value AND am.type='uniprot' AND am.version='2022_04' AND am.best=1 AND am.map=a.acc AND am.afdb=a.afdb AND m.site=a.site AND a.membrane IS NULL GROUP BY m.acc, m.site, m.aa"
print(f("Running query:\n\n{myquery}\n\n"))
system.time(qc <- Query(myquery))
print(f("Done!"))
write_rds(qc, "tmp-qc.rds")

print(f("Wrote to 'tmp-qc.rds'"))
