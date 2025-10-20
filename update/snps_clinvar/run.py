#!/usr/bin/env python3
"""
Run: Run entire pipeline
"""

# Initialize

from blang_mysql import *
from blang import *

# Start
table = "snps_clinvar"

# For variant_summary (TSV) and XML
# rel = "2023-06"
rel = "2023-07"
# For VCF
rel2 = "20230722"   

# VEP (only run if output file doesn't exist yet)
vep_infile = f"input/clinvar_{rel2}_combined.vcf.gz"
vep_outfile = f"input/clinvar_{rel2}_combined.vep.vcf.gz"
if not Exists(vep_outfile):
    # 32 cores
    Run("Run VEP to get protein-level consequence information", f'~/scripts/qsubcm.sh 32 3 ~/update/ensembl-vep/vep --cache --offline --force -v --fork 32 --buffer_size 10000 --format vcf --vcf --compress_output gzip --coding_only --fields="Allele,Gene,Feature_type,Feature,Consequence,HGVSp,Protein_position,Amino_acids,CHECK_REF" -i {vep_infile} -o {vep_outfile}')
    # >> Actually only took 13 minutes for these 2,248,538 variants
    # Note: VEP returns slightly fewer variants. This is because it drops 984 variants that don't have an alternative allele (only ".") (cat first7_vcf_only.vcf | cut -f5 | suq)
    # Waitforjobs()

# Run("Parse flat file release using Perl script", f"main_tsv.pl human {rel}")
# # Run("Parse XML release using new Python script", f"main_xml.py human {rel}")

# VEP-based
Run("Parse VEP VCF file", f"main.py {rel} {rel2}")


# Validating the VEP-based version:
# e.g. KDM3B_HUMAN
# >> The overlap between the ClinVar website (102 variants):
# https://www.ncbi.nlm.nih.gov/clinvar/?term=kdm3b%5Bgene%5D (display options 200, then download)
# ...and snps_clinvar (98 variants):
# query "SELECT CONCAT(original, site, variant) FROM snps_clinvar WHERE name='KDM3B_HUMAN' AND acc NOT LIKE '%-%'"
# ...is good. All snps_clinvar variants are also on the ClinVar website.
# The ClinVar website additionally has 4 more variants for KDM3B: A709P, N351H, P1254R, R245I.
# I can't see an obvious reason why they weren't included on e.g. https://www.ncbi.nlm.nih.gov/clinvar/variation/1676521/ (perhaps it's because the "Genomic location" field is blank).
# https://github.com/SeqOne/clinvcf also mentioned that some variants are mysteriously absent in the official ClinVar VCF files.
# 
# >> Since they are the official ClinVar VCF files, however, I'm going to trust them. The MacArthur lab also uses the offical VCFs: https://github.com/macarthur-lab/clinvar.
#
# >> OK, use VEP-based version!

# Validating KDM3B_HUMAN vs. UniProt:
# https://www.uniprot.org/uniprotkb/Q7LBC6/variant-viewer
# Oddly, UniProt only has very few ClinVar variants (only 14). However, all 14 are in snps_clinvar.
#
# >> OK, use VEP-based version!



print("Done!")
