#!/usr/bin/env python3
"""
Run: Run entire pipeline
"""

# Initialize

from blang_mysql import *
from blang import *

apitoken = "yzCtuDYU3AeC9ldf7uAJPVMeC3OlxNUpPE8Ic5KlP8gQQnRVfWPV1SWd1qbIkUS4"

rel = "2024.01.03"
vep_infile = f"input/mastermind_cited_variants_reference-{rel}-grch38.vcf.gz"
vep_outfile = f"input/mastermind_cited_variants_reference-{rel}-grch38.vep.vcf.gz"
vep_outfile_missense = f"input/mastermind_cited_variants_reference-{rel}-grch38.vep.missense.vcf.gz"

# Start

# Download
# Run("Download", "download.py")

# VEP (only run if output file doesn't exist yet)
if not Exists(vep_outfile):
    # 32 cores, 0.25 GB of memory per core (8 GB in total, of which ~5.5 GB will be used)
    Run("Run VEP to get protein-level consequence information", f'~/scripts/qsubcm.sh 32 0.25 ~/update/ensembl-vep/vep -i {vep_infile} -o {vep_outfile} --cache --offline --force -v --fork 32 --buffer_size 10000 --format vcf --vcf --compress_output gzip --coding_only --fields="Allele,Gene,Feature_type,Feature,Consequence,HGVSp,Protein_position,Amino_acids,CHECK_REF"')
    Waitforjobs()
    # >> 2.25 hours

# Only keep header lines and "missense_variants"
if not Exists(vep_outfile_missense):
    Run(f"Filter '{vep_outfile}' and keep only header and missense_variants", f"zcat {vep_outfile} | grep -P '^# |missense_variant' | gzip > {vep_outfile_missense}")

# Run
Run("Main", f"main.py {vep_outfile_missense}")



# Diagnostic queries

Run("Query: Counts & Averages", '~/scripts/query.pl "SELECT COUNT(*), COUNT(DISTINCT acc), COUNT(DISTINCT acc, site), COUNT(DISTINCT acc, site, original), COUNT(DISTINCT acc, site, original, variant), MIN(papers), AVG(papers), STD(papers), MAX(papers) FROM snps_mastermind" -h')
Run("Query: Top 10 mutations (BRAF V600E etc.)", '~/scripts/query.pl "SELECT * FROM snps_mastermind ORDER BY papers DESC LIMIT 10" -h')
Run("Query: Occurrences per acc|site|original|variant, descending (should be max 1 since main.py prevents multiple inserts)", '~/scripts/query.pl "SELECT acc, site, original, variant, COUNT(*) AS c FROM snps_mastermind GROUP BY acc, site, original, variant ORDER BY c DESC, acc, site, original, variant LIMIT 3" -h')

# There are no variants larger than 3 bases that only affect a single amino acid (which they technically could, due to the +3 wobble position in a preceding codon):
# zcat input/mastermind_cited_variants_reference-2024.01.03-grch38.vep.missense.vcf.gz | g -o "CSQ=.+" | perl -ne 's/^CSQ=//; @a = split(/\|/); ($orig, $var) = split(/\//, $a[7]); if ((length($a[0]) > 3) and (length($var) == 1)) { print; }'
# >> None
# SELECT * FROM snps_mastermind WHERE LENGTH(variantbase)>3 AND LENGTH(variant)=1;
# >> None
# >> Most of the plentiful 3-base variants must actually come from the protein level, and be putative genomic variants, as the documentation says. 2-base ones are probably truncated (simplified), but also based on protein-level variants.

print("Done!")
