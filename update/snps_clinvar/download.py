#!/usr/bin/env python3
"""
Download: Download source data
"""

# Initialize

# from blang_mysql import *
from blang import Run

# For variant_summary (TSV) and XML
rel = "2023-07"

# For VCF
rel2 = "20230722"   

# Download

# https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/clinvar_variation/


# # Flat file release (TSV) (doesn't have NP_... protein IDs, but has binary ClinSigSimple, which can't simply be derived from ClinicalSignificance)
# Run("Download", f"wget -O input/variant_summary_{rel}.txt.gz 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz'");


# XML release (has NP_... protein IDs) (no parser implemented)
Run("Download", f"wget -O input/ClinVarVariationRelease_{rel}.xml.gz 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/clinvar_variation/ClinVarVariationRelease_{rel}.xml.gz'");


# VCF release (main chromosomes)
# >> Doesn't have any effect annotation (only genomic locations)
# >> Could run VEP on this, but:
#   >> Also doesn't have the simple 0/1 clinsig (which is needed to clear up some cases)
#     >> Can pretty easily get that and other annotation from the variant_summary though
Run("Download", f"wget -O input/clinvar_{rel2}.vcf.gz 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_{rel2}.vcf.gz'");
# VCF release (patches etc., but really just contains the pseudoautosomal part of Y) (can concatenate according to https://ftp.ncbi.nlm.nih.gov/pub/clinvar/README_VCF.txt)
Run("Download", f"wget -O input/clinvar_{rel2}_papu.vcf.gz 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_{rel2}_papu.vcf.gz'");
# Combine the two (concatenate)
# Note - clinvar_..._combined.vcf.gz looks smaller than .vcf.gz, but this is only because it is compressed better.
Run("Get header and body from main VCF", f"zcat input/clinvar_{rel2}.vcf.gz > input/clinvar_{rel2}_combined.vcf");
Run("Get only body from extra VCF", f"zcat input/clinvar_{rel2}_papu.vcf.gz | grep -v '^#' >> input/clinvar_{rel2}_combined.vcf");
Run("Compress", f"gzip -f input/clinvar_{rel2}_combined.vcf");


# Show directory
Run("Show directory", "ls -lah input")

print("Done!")
