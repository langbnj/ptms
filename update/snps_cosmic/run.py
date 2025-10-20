#!/usr/bin/env python3
"""
Run: Run entire pipeline
"""

# Initialize

from blang_mysql import *
from blang import *

table = "snps_cosmic"
rel = 99
vep_infiles = (
    f"input/Cosmic_GenomeScreensMutant_Normal_v{rel}_GRCh38.vcf.gz", 
    f"input/Cosmic_NonCodingVariants_Normal_v{rel}_GRCh38.vcf.gz",
    f"input/Cosmic_CompleteTargetedScreensMutant_Normal_v{rel}_GRCh38.vcf.gz", 
    f"input/CellLinesProject_GenomeScreensMutant_Normal_v{rel}_GRCh38.vcf.gz", 
    f"input/CellLinesProject_NonCodingVariants_Normal_v{rel}_GRCh38.vcf.gz", 
    )

# Download
# Run("Download", "download.py")

# Start
Clear(table)

# VEP (only run if output file doesn't exist yet)
for vep_infile in vep_infiles:
    vep_outfile = re.sub(r"\.vcf\.gz$", ".vep.vcf.gz", vep_infile)
    if not Exists(vep_outfile):
        # 32 cores, 0.25 GB of memory per core (8 GB in total, of which ~5.5 GB will be used)
        Run("Run VEP to get protein-level consequence information", f'~/scripts/qsubcm.sh 32 0.25 ~/update/ensembl-vep/vep -i {vep_infile} -o {vep_outfile} --cache --offline --force -v --fork 32 --buffer_size 10000 --format vcf --vcf --compress_output gzip --coding_only --fields="Allele,Gene,Feature_type,Feature,Consequence,HGVSp,Protein_position,Amino_acids,CHECK_REF"')
# >> ~3 hours
Waitforjobs()

# Only keep header lines and "missense_variants"
for vep_infile in vep_infiles:
    vep_outfile = re.sub(r"\.vcf\.gz$", ".vep.vcf.gz", vep_infile)
    vep_outfile_missense = re.sub(r"\.vcf\.gz$", ".vep.missense.vcf.gz", vep_infile)
    if not Exists(vep_outfile_missense):
        Run(f"Filter '{vep_outfile}' and keep only header and missense_variants", f"zcat {vep_outfile} | grep -P '^# |missense_variant' | gzip > {vep_outfile_missense}")

# Disable keys during inserts
Query(f"ALTER TABLE {table} DISABLE KEYS")

# Insert into SQL
Starttime()
for vep_infile in vep_infiles:
    vep_outfile_missense = re.sub(r"\.vcf\.gz$", ".vep.missense.vcf.gz", vep_infile)
    # Run("Main", f"main.py {vep_outfile_missense}")
    Run("Main", f"~/scripts/qsub.sh main.py {vep_outfile_missense}")
Waitforjobs()
Stoptime()

# Re-enable keys
print(f"Enabling keys for table '{table}'...")
Starttime()
Query(f"ALTER TABLE {table} ENABLE KEYS")
Stoptime()

# Sort table
print(f"Sorting table '{table}' by species, name, acc, site, original, variant, source...")
Starttime()
Query(f"ALTER TABLE {table} ORDER BY species, name, acc, site, original, variant, source")
Stoptime()

# Optimize
Starttime()
Optimize(table)
Stoptime()

# To do? Add in all the other information (cancer type etc.) from non-VCF files? I can match on COSV IDs.
# To do? Can also add sample information ("primary" etc.) by matching on COSS IDs.
# To do? However, see download.py for an analysis of how useful/complete these annotations are. They're not particularly necessary.

# Diagnostic queries

Run("Query: Counts & Averages", '~/scripts/query.pl "SELECT COUNT(*), COUNT(DISTINCT acc), COUNT(DISTINCT acc, site), COUNT(DISTINCT acc, site, original), COUNT(DISTINCT acc, site, original, variant), MIN(ac), AVG(ac), STD(ac), MAX(ac) FROM snps_cosmic" -h')
Run("Query: Top 100 mutations (probably BRAF V600E etc.)", '~/scripts/query.pl "SELECT DISTINCT source, cell_line, name, site, original, variant, ac, cgc FROM snps_cosmic ORDER BY ac DESC LIMIT 100" -h')
Run("Query: Allele counts (actually sample counts here, but the column is AC for simplicity)", '~/scripts/query.pl "SELECT ac, COUNT(*) FROM snps_cosmic GROUP BY ac ORDER BY ac LIMIT 10" -h')

print("Done!")
