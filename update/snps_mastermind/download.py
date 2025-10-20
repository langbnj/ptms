#!/usr/bin/env python3
"""
Download: Download source data
"""

# Initialize

# from blang_mysql import *
from blang import *

# Now available at: https://www.genomenon.com/cvr-2/

# Download (wget doesn't manage to get the file name, so perhaps just download by hand)
# VCF
# https://mastermind.genomenon.com/cvr/download?build=grch38
# CSV
# https://mastermind.genomenon.com/cvr/download?build=grch38&format=csv
# rel = "2023.07.02"
rel = "2024.01.03"

# License (PDF)
Run("Download", f"wget -O input/Genomenon-CVR-License-11-09-18.pdf 'https://www.genomenon.com/wp-content/uploads/2018/11/Genomenon-CVR-License-11-09-18.pdf'")
# VCF format
Run("Download", f"wget -O input/mastermind_cited_variants_reference-{rel}-grch38-vcf.zip 'https://mastermind.genomenon.com/cvr/download?build=grch38'")
# # CSV format
# Run("Download", f"wget -O input/mastermind_cited_variants_reference-{rel}-grch38-csv.zip 'https://mastermind.genomenon.com/cvr/download?build=grch38&format=csv'")

# Convert .zip to .gz (a bit problematic since there are multiple files in the .zips)
# VCF format
Run("Convert .zip to .gz", f"unzip -p input/mastermind_cited_variants_reference-{rel}-grch38-vcf.zip '*.vcf' | gzip > input/mastermind_cited_variants_reference-{rel}-grch38.vcf.gz");
# # CSV format
# Run("Convert .zip to .gz", f"unzip -p input/mastermind_cited_variants_reference-{rel}-grch38-csv.zip '*.csv' | gzip > input/mastermind_cited_variants_reference-{rel}-grch38.csv.gz");

# # Unzip
# Cd("input")
# # VCF format
# Run("Decompress (.zip)", f"unzip -o input/mastermind_cited_variants_reference-{rel}-grch38-vcf.zip");
# # CSV format
# Run("Decompress (.zip)", f"unzip -o input/mastermind_cited_variants_reference-{rel}-grch38-csv.zip");
# Cd("..")
# 
# # Recompress
# Run("Compress (.gz)", "gzip")

# Clean up
Run("Clean up", "rm -fv input/mastermind_cited_variants_reference-{rel}-grch38-vcf.zip")
# Run("Clean up", "rm -fv input/mastermind_cited_variants_reference-{rel}-grch38-csv.zip")

# Show directory
Run("Show directory", "ls -lah input")

print("Done!")
