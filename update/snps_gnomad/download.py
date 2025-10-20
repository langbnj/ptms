#!/usr/bin/env python3
"""
Download: Download source data
"""

# Initialize

# from blang_mysql import *
from blang import Run

# gnomAD release versions
rel2 = "2.1.1"
rel3 = "3.1.2"
rel4 = "4.0"

# Download from Google Cloud

# v2 exomes (WES)
# Seems to have 17,205,543 lines (including headers), which is pretty similar to snps_univar (11 million, which took VEP only 45 minutes)
# >> runtime should be okay
# Single file (~60 GB):
# Run(f"Download v{rel2} WES", f"gsutil -m cp gs://gcp-public-data--gnomad/release/{rel2}/liftover_grch38/vcf/exomes/gnomad.exomes.r{rel2}.sites.liftover_grch38.vcf.bgz* input")
# Individual chromosomes
Run(f"Download v{rel2} WES", f"gsutil -m cp gs://gcp-public-data--gnomad/release/{rel2}/liftover_grch38/vcf/exomes/gnomad.exomes.r{rel2}.sites.*.liftover_grch38.vcf.bgz* input")

# v2 genomes (WGS)
Run(f"Download v{rel2} WGS", f"gsutil -m cp gs://gcp-public-data--gnomad/release/{rel2}/liftover_grch38/vcf/genomes/gnomad.genomes.r{rel2}.sites.*.liftover_grch38.vcf.bgz* input")

# v3 genomes (WGS) (no WES available)
# Another 10 times bigger than v2 genomes, but it can be parallelised, so hopefully only a few hours
Run(f"Download v{rel3} WGS", f"gsutil -m cp gs://gcp-public-data--gnomad/release/{rel3}/vcf/genomes/gnomad.genomes.v{rel3}.sites.chr*.vcf.bgz* input")

# # v4 exomes (WES)
# # File sizes are okay.
# # 730,947 exomes, which is amazing. Gigantic. "Minimum viable product" release, which means the files are actually smaller and have less annotation. That's okay for me.
# Run(f"Download v{rel4} WES", f"gsutil -m cp gs://gcp-public-data--gnomad/release/{rel4}/vcf/exomes/gnomad.exomes.v{rel4}.sites.chr*.vcf.bgz* input")
# 
# # v4 genomes (WGS)
# # File sizes are okay.
# # 76,215 genomes, which are actually just the v3 genomes (https://gnomad.broadinstitute.org/news/2023-11-gnomad-v4-0/).
# # They say they changed these genomes only slightly: "Differences between gnomAD v3 and v4: The 76,156 genomes from gnomAD v3.1.2 have been included in gnomAD v4. We have not updated the quality control (QC) of these samples, but we have updated some of the frequencies within the HGDP and 1000 Genomes (HGDP/1KG) subset."
# # >> Limited value in getting these.
# Run(f"Download v{rel4} WGS", f"gsutil -m cp gs://gcp-public-data--gnomad/release/{rel4}/vcf/genomes/gnomad.genomes.v{rel4}.sites.chr*.vcf.bgz* input")



# Notes
# Q: Does gnomAD contain all of ExAC?
# A: Almost (~90% of samples/individuals):
# https://gnomad.broadinstitute.org/help
# "Why is a particular variant found in some versions of ExAC/gnomAD but not others? Likely because of differences between the projects in sample inclusion and variant quality control.
# Sample QC: Not all data from previous call sets is included in each gnomAD release. Most of the samples from prior releases are included during the initial joint calling; however, we make changes to our sample QC process in each release, which always results in the removal of some samples that previously passed. For instance, approximately 10% of the samples in ExAC are missing from gnomAD v2, so we expect some variants in ExAC, particularly those that were found at low frequencies, to be absent in gnomAD v2."




# Note: There are also "constraint" scores, but they are at the gene and transcript level and therefore useless for individual residues and protein regions (https://gnomad.broadinstitute.org/help/constraint, gsutil -m cat gs://gcp-public-data--gnomad/release/4.0/constraint/gnomad.v4.0.constraint_metrics.tsv | b).



# Show directory
Run("Show directory", "ls -lah input")

print("Done!")
