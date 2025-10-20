#!/usr/bin/env python3
"""
Main: Main script
"""

# Initialize
# import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
import gzip
import pickle
from blang_mysql import *
from blang import *

table = "snps_gnomad"

# Fields of interest (defines the order) (also selected via regexps below)
myfields = (
    "ac",
    "ac_xx",
    "ac_xy",
    "ac_grpmax",

    "an",
    "an_xx",
    "an_xy",
    "an_grpmax",

    "af",
    "af_xx",
    "af_xy",
    "af_grpmax",

    "grpmax",
    "n_alt_alleles",
    "allele_type",
    "variant_type",
    "was_mixed",

    "nhomalt",
    "nhomalt_xx",
    "nhomalt_xy",
    "nhomalt_grpmax",
    "inbreeding_coeff",

    "af_controls_and_biobanks",
    "af_non_cancer",
    "af_non_neuro",
    "af_non_topmed",
    "af_non_v2",

    "af_afr",
    "af_ami",
    "af_amr",
    "af_asj",
    "af_eas",
    "af_fin",
    "af_mid",
    "af_nfe",
    "af_sas",
    "af_remaining",

    "fafmax95",
    "fafmax_faf95_max",
    "faf95",
    "faf95_afr",
    "faf95_amr",
    "faf95_eas",
    "faf95_nfe",
    "faf95_sas",

    "as_fs",
    "as_sor",
    "as_pab_max",
    "fs",
    "sor",
    "lcr",
    "non_par",
    "segdup",
)

infile = Args(1, "[input file]", "../output/gnomad.exomes.r2.1.1.sites.liftover_grch38.vep.missense.vcf.gz")
outfile = re.sub(r"\.(b)?gz$", ".tmp_sql_insert.txt", infile)

# Clear(table)



# # Get uniens
# print("Getting enst-to-name/acc/canon/species/ensg mappings from table 'uniens':")
# uniens = {}
# for (enst,) in tq(Query(f"SELECT DISTINCT enst FROM uniens WHERE species='human'")):
#     uniens[enst] = FetchList(Query(f"SELECT DISTINCT name, acc, canon, species, ensg FROM uniens WHERE enst='{enst}'"))
# 
# # Get uniseq
# print("\nGetting acc-to-UniProt-sequence mappings from table 'uniseq':")
# uniseqs = {}
# for (acc,) in tq(Query(f"SELECT DISTINCT acc FROM uniseq WHERE species='human' AND type IN ('UniProt', 'UniIso')")):
#     uniseqs[acc] = FetchOne(Query(f"SELECT DISTINCT seq FROM uniseq WHERE acc='{acc}' AND type IN ('UniProt', 'UniIso')"))

# Get uniens (cached query)
print("Getting enst-to-name/acc/canon/species/ensg mappings from table 'uniens':")
tmpfile = "../tmp-uniens.pkl"
if not Exists(tmpfile):
    uniens = {}
    for (enst,) in tq(Query(f"SELECT DISTINCT enst FROM uniens WHERE species='human'")):
        uniens[enst] = FetchList(Query(f"SELECT DISTINCT name, acc, canon, species, ensg FROM uniens WHERE enst='{enst}'"))
    # Write to file
    with open(tmpfile, "wb") as f:
        pickle.dump(uniens, f)
    print(f"Wrote to '{tmpfile}'")
else:
    # Load from file
    with open(tmpfile, "rb") as f:
        uniens = pickle.load(f)
    print(f"Loaded from '{tmpfile}'")

# Get uniseq (cached query)
print("\nGetting acc-to-UniProt-sequence mappings from table 'uniseq':")
tmpfile = "../tmp-uniseq.pkl"
if not Exists(tmpfile):
    uniseqs = {}
    for (acc,) in tq(Query(f"SELECT DISTINCT acc FROM uniseq WHERE species='human' AND type IN ('UniProt', 'UniIso')")):
        uniseqs[acc] = FetchOne(Query(f"SELECT DISTINCT seq FROM uniseq WHERE acc='{acc}' AND type IN ('UniProt', 'UniIso')"))
    # Write to file
    with open(tmpfile, "wb") as f:
        pickle.dump(uniseqs, f)
    print(f"Wrote to '{tmpfile}'")
else:
    with open(tmpfile, "rb") as f:
        uniseqs = pickle.load(f)
    print(f"Loaded from '{tmpfile}'")



# Start

# Determine source subset (v2 exomes, v2 genomes, v3 genomes)
# gnomad.exomes.r2.1.1.sites.liftover_grch38.vep.vcf.bgz
# gnomad.genomes.r2.1.1.sites.liftover_grch38.vep.vcf.bgz
# gnomad.genomes.v3.1.2.sites.chr10.vep.vcf.bgz
# Also get maximum allele numbers: https://gnomad.broadinstitute.org:
# "The v2.1.1 data set (GRCh37/hg19) provided on this website spans 125,748 (125748) exome sequences and 15,708 (15708) whole-genome sequences from unrelated individuals sequenced as part of various disease-specific and population genetic studies."
# "The v3.1.2 data set (GRCh38) spans 76,156 (76156) genomes of diverse ancestries, selected as in v2."
# This also matches "SELECT source, MAX(an) FROM snps_gnomad GROUP BY source ORDER BY source" on the filled table.
# "The gnomAD v4.0.0 data set contains data from 730,947 (730947) exomes and 76,215 (76215) whole genomes, all mapped to the GRCh38 reference sequence." (https://gnomad.broadinstitute.org/downloads#v4)
# "Today, we are delighted to announce the release of gnomAD v4, which includes data from 807,162 (807162) total individuals." (https://gnomad.broadinstitute.org/news/2023-11-gnomad-v4-0/)



# Diagnostic queries:
# SELECT COUNT(*) FROM snps_gnomad;
# # >> 32,213,842
# SELECT COUNT(*) FROM snps_gnomad WHERE (source='v2_wes' AND an>=251496/2) OR (source='v2_wgs' AND an>=31416/2) OR (source='v3_wgs' AND an>=152312/2);
# # >> 31,970,409 passing >=50% max AN filter
# # 243,433 are filtered out (that's ok, small fraction)
# SELECT COUNT(*) FROM snps_gnomad WHERE (source='v2_wes' AND an<251496/2) OR (source='v2_wgs' AND an<31416/2) OR (source='v3_wgs' AND an<152312/2);
# RENAME TABLE snps_gnomad TO snps_gnomad_backup;
# CREATE TABLE snps_gnomad LIKE snps_gnomad_backup;
# INSERT INTO snps_gnomad SELECT * FROM snps_gnomad_backup WHERE (source='v2_wes' AND an>=251496/2) OR (source='v2_wgs' AND an>=31416/2) OR (source='v3_wgs' AND an>=152312/2);

# Number of males/XY for different sources (estimated from missense variants, as well as from complete chrY VCF files (see also: test_number_of_xy/commands.sh)):

# v2_wes: 125,748 individuals (i.e. expected max AN=251,496)
# v2_wes: zcat gnomad.exomes.r2.1.1.sites.*.liftover_grch38.vep.missense.vcf.bgz | g -o "AN_male=\d+" | g -o "\d+" | max            # >> AN_male=135,922 (autosomes in males) >> number of males should be / 2 = 67,961
# v2_wes: zcat gnomad.exomes.r2.1.1.sites.1.liftover_grch38.vep.missense.vcf.bgz | g -o ";AN_female=\d+" | g -o "\d+" | max         # >> AN_female=115,574 (autosomes in females) >> number of females should be / 2 = 57,787 >> OK, matches up!
# chrY complete VCF, v2_wes: zcat gnomad.exomes.r2.1.1.sites.Y.liftover_grch38.vcf.bgz | g -o "AN_male=\d+" | g -o "\d+" | max      # >> AN=AN_male=67,953 (chrY)

# v2_wgs: 15,708 individuals (i.e. expected max AN=31,416)
# v2_wgs: zcat gnomad.genomes.r2.1.1.sites.*.liftover_grch38.vep.missense.vcf.bgz | g -o "AN_male=\d+" | g -o "\d+" | max           # >> AN_male=17,482 (autosomes in males) >> number of males should be / 2 = 8,741
# v2_wgs: zcat gnomad.genomes.r2.1.1.sites.1.liftover_grch38.vep.missense.vcf.bgz | g -o ";AN_female=\d+" | g -o "\d+" | max        # >> AN_female=13,934 (autosomes in females) >> number of females should be / 2 = 6,967 >> OK, matches up!

# v3_wgs: 76,156 individuals (i.e. expected max AN=152,312)
# v3_wgs: zcat gnomad.genomes.v3.1.2.sites.chr*.vep.missense.vcf.bgz | g -o "AN_XY=\d+" | g -o "\d+" | max                          # >> AN_XY=74,418 (autosomes in males) >> number of males should be / 2 = 37,209
# v3_wgs: zcat gnomad.genomes.v3.1.2.sites.chr1.vep.missense.vcf.bgz | g -o ";AN_XX=\d+" | g -o "\d+" | max                         # >> AN_XX=77,894 (autosomes in females) >> number of females should be / 2 = 38,947 >> OK, matches up!
# chrY complete VCF, v3_wgs: zcat gnomad.genomes.v3.1.2.sites.chrY.vcf.bgz | g -o "AN_XY=\d+" | g -o "\d+" | max                    # >> AN=AN_XY=37,209 (chrY)

# Not joint
# v4_wes: 730,947 individuals (i.e. expected max AN=1,461,894)
# v4_wes: zcat gnomad.exomes.v4.0.sites.chr*.vep.missense.vcf.bgz | g -o "AN_XY=\d+" | g -o "\d+" | datamash max 1                  # >> AN_XY=727,248 (autosomes in males) >> number of males should be / 2 = 363,624
# v4_wes: zcat gnomad.exomes.v4.0.sites.chr1.vep.missense.vcf.bgz | g -o ";AN_XX=\d+" | g -o "\d+" | datamash max 1                 # >> AN_XX=734,646 (autosomes in females) >> number of females should be / 2 = 367,323 >> OK, matches up!
# chrY complete VCF, v4_wes: zcat gnomad.exomes.v4.0.sites.chrY.vcf.bgz | g -o "AN_XY=\d+" | g -o "\d+" | max                       # >> AN=AN_XY=363,593 (chrY)

# Not joint
# v4_wgs: 76,215 individuals (i.e. expected max AN=152,430)
# v4_wgs: zcat gnomad.genomes.v4.0.sites.chr*.vep.missense.vcf.bgz | g -o "AN_XY=\d+" | g -o "\d+" | datamash max 1                 # >> AN_XY=74,546 (autosomes in males) >> number of males should be / 2 = 37,273
# v4_wgs: zcat gnomad.genomes.v4.0.sites.chr1.vep.missense.vcf.bgz | g -o ";AN_XX=\d+" | g -o "\d+" | datamash max 1                # >> AN_XX=77,894 (autosomes in females) >> number of females should be / 2 = 38,947 >> adds up to 76,220 (5 too many, presumably due to XXY etc.)
# chrY complete VCF, v4_wgs: zcat gnomad.genomes.v4.0.sites.chrY.vcf.bgz | g -o "AN_XY=\d+" | g -o "\d+" | max                      # >> AN=AN_XY=37,269 (chrY)
# v4_wgs is the only source where the XY and XX numbers don't add up perfectly, but this makes sense:
# >> "The final list of possible sex karyotypes is: X, XX, XXX, XXXY, XXY, XXYY, XY, XYY, and ambiguous." (https://gnomad.broadinstitute.org/news/2023-11-gnomad-v4-0/)

# joint
# v4_wes: 807,162 individuals (i.e. expected max AN_joint=1,614,324) (actual max AN_joint is 1,614,320, i.e. 4 lower) (...and actual max AC_joint is 1,614,165, i.e. 159 lower - it of course should normally be lower since the reference allele should be present, too)
# v4_wes: zcat gnomad.exomes.v4.0.sites.chr*.vep.missense.vcf.bgz | g -o ";AN_joint_XY=\d+" | g -o "\d+" | datamash max 1           # >> AN_joint_XY=801,794 (autosomes in males) >> number of males should be / 2 = 400,897
# v4_wes: zcat gnomad.exomes.v4.0.sites.chr1.vep.missense.vcf.bgz | g -o ";AN_joint_XX=\d+" | g -o "\d+" | datamash max 1           # >> AN_joint_XX=812,540 (autosomes in females) >> number of females should be / 2 = 406,270 >> adds up to 807,167 (5 too many, presumably due to XXY etc.)
# chrY complete VCF, v4_wes: zcat gnomad.exomes.v4.0.sites.chrY.vcf.bgz | g -o ";AN_joint_XY=\d+" | g -o "\d+" | max                # >> AN_joint=AN_joint_XY=398,004 (chrY)
# For joint v4, the XY and XX numbers don't add up perfectly, but this makes sense:
# >> "The final list of possible sex karyotypes is: X, XX, XXX, XXXY, XXY, XXYY, XY, XYY, and ambiguous." (https://gnomad.broadinstitute.org/news/2023-11-gnomad-v4-0/)

# joint
# v4_wgs: 807,162 individuals (i.e. expected max AN_joint=1,614,324) (actual max AN_joint is 1,614,320, i.e. 4 lower) (...and actual max AC_joint is 1,614,165, i.e. 159 lower - it of course should normally be lower since the reference allele should be present, too)
# v4_wgs: zcat gnomad.genomes.v4.0.sites.chr*.vep.missense.vcf.bgz | g -o ";AN_joint_XY=\d+" | g -o "\d+" | datamash max 1          # >> AN_joint_XY=801,794 (autosomes in males) >> number of males should be / 2 = 400,897
# v4_wgs: zcat gnomad.genomes.v4.0.sites.chr1.vep.missense.vcf.bgz | g -o ";AN_joint_XX=\d+" | g -o "\d+" | datamash max 1          # >> AN_joint_XX=812,540 (autosomes in females) >> number of females should be / 2 = 406,270 >> adds up to 807,167 (5 too many, presumably due to XXY etc.)
# chrY complete VCF, v4_wgs: zcat gnomad.genomes.v4.0.sites.chrY.vcf.bgz | g -o ";AN_joint_XY=\d+" | g -o "\d+" | max               # >> AN_joint=AN_joint_XY=398,004 (chrY)
# For joint v4, the XY and XX numbers don't add up perfectly, but this makes sense:
# >> "The final list of possible sex karyotypes is: X, XX, XXX, XXXY, XXY, XXYY, XY, XYY, and ambiguous." (https://gnomad.broadinstitute.org/news/2023-11-gnomad-v4-0/)

# chrX AN_male/AN_XY:
# v2_wes:     135,922 >> OK!
# v2_wgs:     17,482  >> OK!
# v3_wgs:     74,418  >> OK!
# v4_wes:     727,248 >> OK!
# v4_wgs:     74,546  >> OK!
# v4_joint:   801,716 >> (AN_joint) Expected 78 more (801,794), but that's ok
# chrX AN:
# v2_wes:     251,496 >> OK!
# v2_wgs:      31,416 >> OK!
# v3_wgs:     152,312 >> OK!
# v4_wes:   1,461,890 >> Expected 4 more (1,461,894), but that's almost identical
# v4_wgs:     152,430 >> OK!
# v4_joint: 1,614,172 >> (AN_joint) Expected 152 more (1,614,324), but that's ok

# >> Can use the number of males/XY as max AN for chrY (which then gets divided by 2 below to arrive at a 50% threshold):
# v2_wes:    67,961 = 135,922 / 2   ( 67961 = 135922 / 2)
# v2_wgs:     8,741 =  17,482 / 2   (  8741 =  17482 / 2)
# v3_wgs:    37,209 =  74,418 / 2   ( 37209 =  74418 / 2)
# v4_wes:   363,624 = 727,248 / 2   (363624 = 727248 / 2)
# v4_wgs:    37,273 =  74,546 / 2   ( 37273 =  74546 / 2)
# v4_joint: 400,897 = 801,794 / 2   (400897 = 801794 / 2)

# >> Actually, the number of males also needs to be halved (males are haploid for Y):
# >> Thresholds:
# v2_wes:    33,980.5 =  67,961 / 2   ( 67961 / 2)
# v2_wgs:     4,370.5 =   8,741 / 2   (  8741 / 2)
# v3_wgs:    18,604.5 =  37,209 / 2   ( 37209 / 2)
# v4_wes:   181,812   = 363,624 / 2   (363624 / 2)
# v4_wgs:    18,636.5 =  37,273 / 2   ( 37273 / 2)
# v4_joint: 200,448.5 = 400,897 / 2   (400897 / 2)
    
# >> OK! I verified these values in queries, and the AC ≥ 50% max(AN) thresholds got applied correctly, including the special ones for chrY.

if "gnomad.exomes.r2." in infile:
    source = "v2_wes"
    anmax = 125748 * 2
    anmaxy = 67961
elif "gnomad.genomes.r2." in infile:
    source = "v2_wgs"
    anmax = 15708 * 2
    anmaxy = 8741
elif "gnomad.genomes.v3." in infile:
    source = "v3_wgs"
    anmax = 76156 * 2
    anmaxy = 37209
elif "gnomad.exomes.v4." in infile:
    source = "v4_wes"
    # anmax = 730947 * 2
    # anmaxy = 363624
    anmax = 807162 * 2   # v4 wes/wgs joint value
    anmaxy = 400897      # v4 wes/wgs joint value
elif "gnomad.genomes.v4." in infile:
    source = "v4_wgs"
    # anmax = 76215 * 2
    # anmaxy = 37273
    anmax = 807162 * 2   # v4 wes/wgs joint value
    anmaxy = 400897      # v4 wes/wgs joint value
else:
    Die(f"Error: Unhandled gnomAD release for '{infile}'")

# print(f"\nReading '{infile}' and inserting into table '{table}':")
print(f"\nReading '{infile}' and writing to temporary file '{outfile}':")

# # Disable keys during inserts
# Query(f"ALTER TABLE {table} DISABLE KEYS")

# e.g.
# chr1	1675544	rs1	C	G	179.81	.	Sources=ExAC,gnomAD,ClinGen

# Read line by line
Time(1)
with gzip.open(infile, 'rt') as f:
    with open(outfile, 'w') as out:

        # Print header (not really required, LOAD DATA LOCAL INFILE will ignore it below - just for clarity)
        # print("id\tname\tacc\tcanon\tensg\tenst\tspecies\tsite\toriginal\tvariant\tsource\t" + "\t".join(myfields), file=out)
        print("id\tchr\tpos\toriginalbase\tvariantbase\trsid\tname\tacc\tcanon\tensg\tenst\tspecies\tsite\toriginal\tvariant\tsource\t" + "\t".join(myfields), file=out)

        for line in tq(f):

            # Skip header
            if line.startswith("#"):
                continue

            # Interesting fields exclusively from v2 exomes:

            # Interesting fields exclusively from v2 genomes:
            # Note: XX and XY are called female and male here instead
            # AS_pab_max is called pab_max (definition is identical) and AS_SOR is called SOR
            ##INFO=<ID=pab_max,Number=A,Type=Float,Description="Maximum p-value over callset for binomial test of observed allele balance for a heterozygous genotype, given expectation of AB=0.5">
            # Could "decoy" be from liftover? Not sure:
            ##INFO=<ID=decoy,Number=0,Type=Flag,Description="Variant falls within a reference decoy region">
            ##INFO=<ID=has_star,Number=0,Type=Flag,Description="Variant locus coincides with a spanning deletion (represented by a star) observed elsewhere in the callset">


            # Interesting fields (from v3 genomes):

            # Subfields that exist for most interesting fields:
            # geographic populations: afr, ami, amr, asj, eas, fin, mid, nfe, oth, sas
            # XX, XY, (female), (male), controls_and_biobanks, non_cancer, non_neuro, non_topmed, non_v2

            # Subsets:
            # https://gnomad.broadinstitute.org/news/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/
            # """
            # By popular request, we’ve computed allele counts and allele frequencies for major subsets of gnomAD v3.1, as follows:
            # 
            # Non-v2: only genome samples that are new to the v3 or v3.1 release and not included in the v2 genomes
            # Non-TOPMed: only samples that are not present in the Trans-Omics for Precision Medicine (TOPMed)/BRAVO release. The allele counts in this subset can thus be added to those of BRAVO to enable federated use of both datasets
            # Non-cancer: only samples from individuals who were not ascertained for having cancer in a cancer study
            # Controls and biobanks: only samples collected specifically as controls for disease studies, or samples belonging to biobanks (e.g. BioMe, Genizon) or general population studies (e.g., 1000 Genomes, HGDP, PAGE)
            # Non-neuro: only samples that were not collected as part of a neurologic or psychiatric case/control study, or samples collected as part of a neurologic or psychiatric case/control study but designated as controls                
            # """

            ##INFO=<ID=AC,Number=A,Type=Integer,Description="Alternate allele count">
            ##INFO=<ID=AC_XX,Number=A,Type=Integer,Description="Alternate allele count for XX samples">
            ##INFO=<ID=AC_XY,Number=A,Type=Integer,Description="Alternate allele count for XY samples">
            ##INFO=<ID=AC_popmax,Number=A,Type=Integer,Description="Allele count in the population with the maximum allele frequency">

            ##INFO=<ID=AC_controls_and_biobanks,Number=A,Type=Integer,Description="Alternate allele count in controls_and_biobanks subset">
            ##INFO=<ID=AC_controls_and_biobanks_XX,Number=A,Type=Integer,Description="Alternate allele count for XX samples in controls_and_biobanks subset">
            ##INFO=<ID=AC_controls_and_biobanks_XY,Number=A,Type=Integer,Description="Alternate allele count for XY samples in controls_and_biobanks subset">
            ##INFO=<ID=AC_non_cancer,Number=A,Type=Integer,Description="Alternate allele count in non_cancer subset">
            ##INFO=<ID=AC_non_cancer,Number=A,Type=Integer,Description="Alternate allele count in non_cancer subset">
            ##INFO=<ID=AC_non_cancer_XX,Number=A,Type=Integer,Description="Alternate allele count for XX samples in non_cancer subset">
            ##INFO=<ID=AC_non_cancer_XY,Number=A,Type=Integer,Description="Alternate allele count for XY samples in non_cancer subset">
            ##INFO=<ID=AC_non_neuro,Number=A,Type=Integer,Description="Alternate allele count in non_neuro subset">
            ##INFO=<ID=AC_non_neuro_XX,Number=A,Type=Integer,Description="Alternate allele count for XX samples in non_neuro subset">
            ##INFO=<ID=AC_non_neuro_XY,Number=A,Type=Integer,Description="Alternate allele count for XY samples in non_neuro subset">
            ##INFO=<ID=AC_non_topmed,Number=A,Type=Integer,Description="Alternate allele count in non_topmed subset">
            ##INFO=<ID=AC_non_topmed_XX,Number=A,Type=Integer,Description="Alternate allele count for XX samples in non_topmed subset">
            ##INFO=<ID=AC_non_topmed_XY,Number=A,Type=Integer,Description="Alternate allele count for XY samples in non_topmed subset">
            ##INFO=<ID=AC_non_v2,Number=A,Type=Integer,Description="Alternate allele count in non_v2 subset">
            ##INFO=<ID=AC_non_v2_XX,Number=A,Type=Integer,Description="Alternate allele count for XX samples in non_v2 subset">
            ##INFO=<ID=AC_non_v2_XY,Number=A,Type=Integer,Description="Alternate allele count for XY samples in non_v2 subset">

            ##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency">
            ##INFO=<ID=AF_XX,Number=A,Type=Float,Description="Alternate allele frequency in XX samples">
            ##INFO=<ID=AF_XY,Number=A,Type=Float,Description="Alternate allele frequency in XY samples">
            ##INFO=<ID=AF_popmax,Number=A,Type=Float,Description="Maximum allele frequency across populations">

            ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles">
            ##INFO=<ID=AN_XX,Number=1,Type=Integer,Description="Total number of alleles in XX samples">
            ##INFO=<ID=AN_XY,Number=1,Type=Integer,Description="Total number of alleles in XY samples">
            ##INFO=<ID=AN_popmax,Number=A,Type=Integer,Description="Total number of alleles in the population with the maximum allele frequency">

            ##INFO=<ID=popmax,Number=A,Type=String,Description="Population with maximum allele frequency">

            ##INFO=<ID=AS_SOR,Number=A,Type=Float,Description="Allele-specific strand bias estimated by the symmetric odds ratio test">
            ##INFO=<ID=AS_pab_max,Number=A,Type=Float,Description="Maximum p-value over callset for binomial test of observed allele balance for a heterozygous genotype, given expectation of 0.5">
            ##INFO=<ID=AS_FS,Number=A,Type=Float,Description="Allele-specific phred-scaled p-value of Fisher's exact test for strand bias">
            ##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value of Fisher's exact test for strand bias">

            ##INFO=<ID=InbreedingCoeff,Number=A,Type=Float,Description="Inbreeding coefficient, the excess heterozygosity at a variant site, computed as 1 - (the number of heterozygous genotypes)/(the number of heterozygous genotypes expected under Hardy-Weinberg equilibrium)">

            ##INFO=<ID=lcr,Number=0,Type=Flag,Description="Variant falls within a low complexity region">

            ##INFO=<ID=monoallelic,Number=0,Type=Flag,Description="All samples are homozygous alternate for the variant">
            ##INFO=<ID=n_alt_alleles,Number=1,Type=Integer,Description="Total number of alternate alleles observed at variant locus">

            ##INFO=<ID=nhomalt,Number=A,Type=Integer,Description="Count of homozygous individuals">
            ##INFO=<ID=nhomalt_XX,Number=A,Type=Integer,Description="Count of homozygous individuals in XX samples">
            ##INFO=<ID=nhomalt_XY,Number=A,Type=Integer,Description="Count of homozygous individuals in XY samples">
            ##INFO=<ID=nhomalt_popmax,Number=A,Type=Integer,Description="Count of homozygous individuals in the population with the maximum allele frequency">

            ##INFO=<ID=nonpar,Number=0,Type=Flag,Description="Variant (on sex chromosome) falls outside a pseudoautosomal region">
            ##INFO=<ID=segdup,Number=0,Type=Flag,Description="Variant falls within a segmental duplication region">

            ##INFO=<ID=transmitted_singleton,Number=0,Type=Flag,Description="Variant was a callset-wide doubleton that was transmitted within a family from a parent to a child (i.e., a singleton amongst unrelated samples in cohort)">

            ##INFO=<ID=allele_type,Number=1,Type=String,Description="Allele type (snv, insertion, deletion, or mixed)">
            ##INFO=<ID=variant_type,Number=1,Type=String,Description="Variant type (snv, indel, multi-snv, multi-indel, or mixed)">
            ##INFO=<ID=was_mixed,Number=0,Type=Flag,Description="Variant type was mixed">

            ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Gene|Feature_type|Feature|Consequence|HGVSp|Protein_position|Amino_acids|CHECK_REF">

            # Filter column values:
            ##FILTER=<ID=AC0,Description="Allele count is zero after filtering out low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het calls)">
            ##FILTER=<ID=AS_VQSR,Description="Failed VQSR filtering thresholds of -2.7739 for SNPs and -1.0606 for indels">
            ##FILTER=<ID=InbreedingCoeff,Description="InbreedingCoeff < -0.3">
            ##FILTER=<ID=PASS,Description="Passed all variant filters">

            # Example line 1 (rs62635282, noncoding SNV):
            # v2 exomes (WES)
            # chr1	12198	rs62635282	G	C	9876.24	AC0	AC=0;AC_afr=0;AC_afr_female=0;AC_afr_male=0;AC_amr=0;AC_amr_female=0;AC_amr_male=0;AC_asj=0;AC_asj_female=0;AC_asj_male=0;AC_eas=0;AC_eas_female=0;AC_eas_jpn=0;AC_eas_kor=0;AC_eas_male=0;AC_eas_oea=0;AC_female=0;AC_fin=0;AC_fin_female=0;AC_fin_male=0;AC_male=0;AC_nfe=0;AC_nfe_bgr=0;AC_nfe_est=0;AC_nfe_female=0;AC_nfe_male=0;AC_nfe_nwe=0;AC_nfe_onf=0;AC_nfe_seu=0;AC_nfe_swe=0;AC_oth=0;AC_oth_female=0;AC_oth_male=0;AC_raw=227;AC_sas=0;AC_sas_female=0;AC_sas_male=0;AF_raw=4.57108e-02;AN=0;AN_afr=0;AN_afr_female=0;AN_afr_male=0;AN_amr=0;AN_amr_female=0;AN_amr_male=0;AN_asj=0;AN_asj_female=0;AN_asj_male=0;AN_eas=0;AN_eas_female=0;AN_eas_jpn=0;AN_eas_kor=0;AN_eas_male=0;AN_eas_oea=0;AN_female=0;AN_fin=0;AN_fin_female=0;AN_fin_male=0;AN_male=0;AN_nfe=0;AN_nfe_bgr=0;AN_nfe_est=0;AN_nfe_female=0;AN_nfe_male=0;AN_nfe_nwe=0;AN_nfe_onf=0;AN_nfe_seu=0;AN_nfe_swe=0;AN_oth=0;AN_oth_female=0;AN_oth_male=0;AN_raw=4966;AN_sas=0;AN_sas_female=0;AN_sas_male=0;BaseQRankSum=0.00000e+00;ClippingRankSum=3.58000e-01;DP=9204;FS=0.00000e+00;InbreedingCoeff=9.80000e-03;MQ=2.30400e+01;MQRankSum=7.36000e-01;OriginalContig=1;OriginalStart=12198;QD=1.39500e+01;ReadPosRankSum=7.36000e-01;SOR=3.02000e-01;VQSLOD=1.01000e+00;VQSR_culprit=MQ;ab_hist_alt_bin_freq=0|0|0|0|1|0|2|0|2|0|10|0|1|28|0|3|0|0|0|0;age_hist_het_bin_freq=0|0|0|0|0|0|0|0|0|0;age_hist_het_n_larger=0;age_hist_het_n_smaller=0;age_hist_hom_bin_freq=0|0|0|0|0|0|0|0|0|0;age_hist_hom_n_larger=0;age_hist_hom_n_smaller=0;allele_type=snv;controls_AC=0;controls_AC_afr=0;controls_AC_afr_female=0;controls_AC_afr_male=0;controls_AC_amr=0;controls_AC_amr_female=0;controls_AC_amr_male=0;controls_AC_asj=0;controls_AC_asj_female=0;controls_AC_asj_male=0;controls_AC_eas=0;controls_AC_eas_female=0;controls_AC_eas_jpn=0;controls_AC_eas_kor=0;controls_AC_eas_male=0;controls_AC_eas_oea=0;controls_AC_female=0;controls_AC_fin=0;controls_AC_fin_female=0;controls_AC_fin_male=0;controls_AC_male=0;controls_AC_nfe=0;controls_AC_nfe_bgr=0;controls_AC_nfe_est=0;controls_AC_nfe_female=0;controls_AC_nfe_male=0;controls_AC_nfe_nwe=0;controls_AC_nfe_onf=0;controls_AC_nfe_seu=0;controls_AC_nfe_swe=0;controls_AC_oth=0;controls_AC_oth_female=0;controls_AC_oth_male=0;controls_AC_raw=109;controls_AC_sas=0;controls_AC_sas_female=0;controls_AC_sas_male=0;controls_AF_raw=4.66610e-02;controls_AN=0;controls_AN_afr=0;controls_AN_afr_female=0;controls_AN_afr_male=0;controls_AN_amr=0;controls_AN_amr_female=0;controls_AN_amr_male=0;controls_AN_asj=0;controls_AN_asj_female=0;controls_AN_asj_male=0;controls_AN_eas=0;controls_AN_eas_female=0;controls_AN_eas_jpn=0;controls_AN_eas_kor=0;controls_AN_eas_male=0;controls_AN_eas_oea=0;controls_AN_female=0;controls_AN_fin=0;controls_AN_fin_female=0;controls_AN_fin_male=0;controls_AN_male=0;controls_AN_nfe=0;controls_AN_nfe_bgr=0;controls_AN_nfe_est=0;controls_AN_nfe_female=0;controls_AN_nfe_male=0;controls_AN_nfe_nwe=0;controls_AN_nfe_onf=0;controls_AN_nfe_seu=0;controls_AN_nfe_swe=0;controls_AN_oth=0;controls_AN_oth_female=0;controls_AN_oth_male=0;controls_AN_raw=2336;controls_AN_sas=0;controls_AN_sas_female=0;controls_AN_sas_male=0;controls_faf95=0.00000e+00;controls_faf95_afr=0.00000e+00;controls_faf95_amr=0.00000e+00;controls_faf95_eas=0.00000e+00;controls_faf95_nfe=0.00000e+00;controls_faf95_sas=0.00000e+00;controls_faf99=0.00000e+00;controls_faf99_afr=0.00000e+00;controls_faf99_amr=0.00000e+00;controls_faf99_eas=0.00000e+00;controls_faf99_nfe=0.00000e+00;controls_faf99_sas=0.00000e+00;controls_nhomalt=0;controls_nhomalt_afr=0;controls_nhomalt_afr_female=0;controls_nhomalt_afr_male=0;controls_nhomalt_amr=0;controls_nhomalt_amr_female=0;controls_nhomalt_amr_male=0;controls_nhomalt_asj=0;controls_nhomalt_asj_female=0;controls_nhomalt_asj_male=0;controls_nhomalt_eas=0;controls_nhomalt_eas_female=0;controls_nhomalt_eas_jpn=0;controls_nhomalt_eas_kor=0;controls_nhomalt_eas_male=0;controls_nhomalt_eas_oea=0;controls_nhomalt_female=0;controls_nhomalt_fin=0;controls_nhomalt_fin_female=0;controls_nhomalt_fin_male=0;controls_nhomalt_male=0;controls_nhomalt_nfe=0;controls_nhomalt_nfe_bgr=0;controls_nhomalt_nfe_est=0;controls_nhomalt_nfe_female=0;controls_nhomalt_nfe_male=0;controls_nhomalt_nfe_nwe=0;controls_nhomalt_nfe_onf=0;controls_nhomalt_nfe_seu=0;controls_nhomalt_nfe_swe=0;controls_nhomalt_oth=0;controls_nhomalt_oth_female=0;controls_nhomalt_oth_male=0;controls_nhomalt_raw=44;controls_nhomalt_sas=0;controls_nhomalt_sas_female=0;controls_nhomalt_sas_male=0;dp_hist_all_bin_freq=125724|24|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;dp_hist_all_n_larger=0;dp_hist_alt_bin_freq=130|7|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;dp_hist_alt_n_larger=0;faf95=0.00000e+00;faf95_afr=0.00000e+00;faf95_amr=0.00000e+00;faf95_eas=0.00000e+00;faf95_nfe=0.00000e+00;faf95_sas=0.00000e+00;faf99=0.00000e+00;faf99_afr=0.00000e+00;faf99_amr=0.00000e+00;faf99_eas=0.00000e+00;faf99_nfe=0.00000e+00;faf99_sas=0.00000e+00;gq_hist_all_bin_freq=1898|511|26|28|8|4|0|5|2|0|0|0|1|0|0|0|0|0|0|0;gq_hist_alt_bin_freq=14|78|1|25|7|4|0|5|2|0|0|0|1|0|0|0|0|0|0|0;n_alt_alleles=1;nhomalt=0;nhomalt_afr=0;nhomalt_afr_female=0;nhomalt_afr_male=0;nhomalt_amr=0;nhomalt_amr_female=0;nhomalt_amr_male=0;nhomalt_asj=0;nhomalt_asj_female=0;nhomalt_asj_male=0;nhomalt_eas=0;nhomalt_eas_female=0;nhomalt_eas_jpn=0;nhomalt_eas_kor=0;nhomalt_eas_male=0;nhomalt_eas_oea=0;nhomalt_female=0;nhomalt_fin=0;nhomalt_fin_female=0;nhomalt_fin_male=0;nhomalt_male=0;nhomalt_nfe=0;nhomalt_nfe_bgr=0;nhomalt_nfe_est=0;nhomalt_nfe_female=0;nhomalt_nfe_male=0;nhomalt_nfe_nwe=0;nhomalt_nfe_onf=0;nhomalt_nfe_seu=0;nhomalt_nfe_swe=0;nhomalt_oth=0;nhomalt_oth_female=0;nhomalt_oth_male=0;nhomalt_raw=90;nhomalt_sas=0;nhomalt_sas_female=0;nhomalt_sas_male=0;non_cancer_AC=0;non_cancer_AC_afr=0;non_cancer_AC_afr_female=0;non_cancer_AC_afr_male=0;non_cancer_AC_amr=0;non_cancer_AC_amr_female=0;non_cancer_AC_amr_male=0;non_cancer_AC_asj=0;non_cancer_AC_asj_female=0;non_cancer_AC_asj_male=0;non_cancer_AC_eas=0;non_cancer_AC_eas_female=0;non_cancer_AC_eas_jpn=0;non_cancer_AC_eas_kor=0;non_cancer_AC_eas_male=0;non_cancer_AC_eas_oea=0;non_cancer_AC_female=0;non_cancer_AC_fin=0;non_cancer_AC_fin_female=0;non_cancer_AC_fin_male=0;non_cancer_AC_male=0;non_cancer_AC_nfe=0;non_cancer_AC_nfe_bgr=0;non_cancer_AC_nfe_est=0;non_cancer_AC_nfe_female=0;non_cancer_AC_nfe_male=0;non_cancer_AC_nfe_nwe=0;non_cancer_AC_nfe_onf=0;non_cancer_AC_nfe_seu=0;non_cancer_AC_nfe_swe=0;non_cancer_AC_oth=0;non_cancer_AC_oth_female=0;non_cancer_AC_oth_male=0;non_cancer_AC_raw=227;non_cancer_AC_sas=0;non_cancer_AC_sas_female=0;non_cancer_AC_sas_male=0;non_cancer_AF_raw=4.57293e-02;non_cancer_AN=0;non_cancer_AN_afr=0;non_cancer_AN_afr_female=0;non_cancer_AN_afr_male=0;non_cancer_AN_amr=0;non_cancer_AN_amr_female=0;non_cancer_AN_amr_male=0;non_cancer_AN_asj=0;non_cancer_AN_asj_female=0;non_cancer_AN_asj_male=0;non_cancer_AN_eas=0;non_cancer_AN_eas_female=0;non_cancer_AN_eas_jpn=0;non_cancer_AN_eas_kor=0;non_cancer_AN_eas_male=0;non_cancer_AN_eas_oea=0;non_cancer_AN_female=0;non_cancer_AN_fin=0;non_cancer_AN_fin_female=0;non_cancer_AN_fin_male=0;non_cancer_AN_male=0;non_cancer_AN_nfe=0;non_cancer_AN_nfe_bgr=0;non_cancer_AN_nfe_est=0;non_cancer_AN_nfe_female=0;non_cancer_AN_nfe_male=0;non_cancer_AN_nfe_nwe=0;non_cancer_AN_nfe_onf=0;non_cancer_AN_nfe_seu=0;non_cancer_AN_nfe_swe=0;non_cancer_AN_oth=0;non_cancer_AN_oth_female=0;non_cancer_AN_oth_male=0;non_cancer_AN_raw=4964;non_cancer_AN_sas=0;non_cancer_AN_sas_female=0;non_cancer_AN_sas_male=0;non_cancer_faf95=0.00000e+00;non_cancer_faf95_afr=0.00000e+00;non_cancer_faf95_amr=0.00000e+00;non_cancer_faf95_eas=0.00000e+00;non_cancer_faf95_nfe=0.00000e+00;non_cancer_faf95_sas=0.00000e+00;non_cancer_faf99=0.00000e+00;non_cancer_faf99_afr=0.00000e+00;non_cancer_faf99_amr=0.00000e+00;non_cancer_faf99_eas=0.00000e+00;non_cancer_faf99_nfe=0.00000e+00;non_cancer_faf99_sas=0.00000e+00;non_cancer_nhomalt=0;non_cancer_nhomalt_afr=0;non_cancer_nhomalt_afr_female=0;non_cancer_nhomalt_afr_male=0;non_cancer_nhomalt_amr=0;non_cancer_nhomalt_amr_female=0;non_cancer_nhomalt_amr_male=0;non_cancer_nhomalt_asj=0;non_cancer_nhomalt_asj_female=0;non_cancer_nhomalt_asj_male=0;non_cancer_nhomalt_eas=0;non_cancer_nhomalt_eas_female=0;non_cancer_nhomalt_eas_jpn=0;non_cancer_nhomalt_eas_kor=0;non_cancer_nhomalt_eas_male=0;non_cancer_nhomalt_eas_oea=0;non_cancer_nhomalt_female=0;non_cancer_nhomalt_fin=0;non_cancer_nhomalt_fin_female=0;non_cancer_nhomalt_fin_male=0;non_cancer_nhomalt_male=0;non_cancer_nhomalt_nfe=0;non_cancer_nhomalt_nfe_bgr=0;non_cancer_nhomalt_nfe_est=0;non_cancer_nhomalt_nfe_female=0;non_cancer_nhomalt_nfe_male=0;non_cancer_nhomalt_nfe_nwe=0;non_cancer_nhomalt_nfe_onf=0;non_cancer_nhomalt_nfe_seu=0;non_cancer_nhomalt_nfe_swe=0;non_cancer_nhomalt_oth=0;non_cancer_nhomalt_oth_female=0;non_cancer_nhomalt_oth_male=0;non_cancer_nhomalt_raw=90;non_cancer_nhomalt_sas=0;non_cancer_nhomalt_sas_female=0;non_cancer_nhomalt_sas_male=0;non_neuro_AC=0;non_neuro_AC_afr=0;non_neuro_AC_afr_female=0;non_neuro_AC_afr_male=0;non_neuro_AC_amr=0;non_neuro_AC_amr_female=0;non_neuro_AC_amr_male=0;non_neuro_AC_asj=0;non_neuro_AC_asj_female=0;non_neuro_AC_asj_male=0;non_neuro_AC_eas=0;non_neuro_AC_eas_female=0;non_neuro_AC_eas_jpn=0;non_neuro_AC_eas_kor=0;non_neuro_AC_eas_male=0;non_neuro_AC_eas_oea=0;non_neuro_AC_female=0;non_neuro_AC_fin=0;non_neuro_AC_fin_female=0;non_neuro_AC_fin_male=0;non_neuro_AC_male=0;non_neuro_AC_nfe=0;non_neuro_AC_nfe_bgr=0;non_neuro_AC_nfe_est=0;non_neuro_AC_nfe_female=0;non_neuro_AC_nfe_male=0;non_neuro_AC_nfe_nwe=0;non_neuro_AC_nfe_onf=0;non_neuro_AC_nfe_seu=0;non_neuro_AC_nfe_swe=0;non_neuro_AC_oth=0;non_neuro_AC_oth_female=0;non_neuro_AC_oth_male=0;non_neuro_AC_raw=225;non_neuro_AC_sas=0;non_neuro_AC_sas_female=0;non_neuro_AC_sas_male=0;non_neuro_AF_raw=4.70908e-02;non_neuro_AN=0;non_neuro_AN_afr=0;non_neuro_AN_afr_female=0;non_neuro_AN_afr_male=0;non_neuro_AN_amr=0;non_neuro_AN_amr_female=0;non_neuro_AN_amr_male=0;non_neuro_AN_asj=0;non_neuro_AN_asj_female=0;non_neuro_AN_asj_male=0;non_neuro_AN_eas=0;non_neuro_AN_eas_female=0;non_neuro_AN_eas_jpn=0;non_neuro_AN_eas_kor=0;non_neuro_AN_eas_male=0;non_neuro_AN_eas_oea=0;non_neuro_AN_female=0;non_neuro_AN_fin=0;non_neuro_AN_fin_female=0;non_neuro_AN_fin_male=0;non_neuro_AN_male=0;non_neuro_AN_nfe=0;non_neuro_AN_nfe_bgr=0;non_neuro_AN_nfe_est=0;non_neuro_AN_nfe_female=0;non_neuro_AN_nfe_male=0;non_neuro_AN_nfe_nwe=0;non_neuro_AN_nfe_onf=0;non_neuro_AN_nfe_seu=0;non_neuro_AN_nfe_swe=0;non_neuro_AN_oth=0;non_neuro_AN_oth_female=0;non_neuro_AN_oth_male=0;non_neuro_AN_raw=4778;non_neuro_AN_sas=0;non_neuro_AN_sas_female=0;non_neuro_AN_sas_male=0;non_neuro_faf95=0.00000e+00;non_neuro_faf95_afr=0.00000e+00;non_neuro_faf95_amr=0.00000e+00;non_neuro_faf95_eas=0.00000e+00;non_neuro_faf95_nfe=0.00000e+00;non_neuro_faf95_sas=0.00000e+00;non_neuro_faf99=0.00000e+00;non_neuro_faf99_afr=0.00000e+00;non_neuro_faf99_amr=0.00000e+00;non_neuro_faf99_eas=0.00000e+00;non_neuro_faf99_nfe=0.00000e+00;non_neuro_faf99_sas=0.00000e+00;non_neuro_nhomalt=0;non_neuro_nhomalt_afr=0;non_neuro_nhomalt_afr_female=0;non_neuro_nhomalt_afr_male=0;non_neuro_nhomalt_amr=0;non_neuro_nhomalt_amr_female=0;non_neuro_nhomalt_amr_male=0;non_neuro_nhomalt_asj=0;non_neuro_nhomalt_asj_female=0;non_neuro_nhomalt_asj_male=0;non_neuro_nhomalt_eas=0;non_neuro_nhomalt_eas_female=0;non_neuro_nhomalt_eas_jpn=0;non_neuro_nhomalt_eas_kor=0;non_neuro_nhomalt_eas_male=0;non_neuro_nhomalt_eas_oea=0;non_neuro_nhomalt_female=0;non_neuro_nhomalt_fin=0;non_neuro_nhomalt_fin_female=0;non_neuro_nhomalt_fin_male=0;non_neuro_nhomalt_male=0;non_neuro_nhomalt_nfe=0;non_neuro_nhomalt_nfe_bgr=0;non_neuro_nhomalt_nfe_est=0;non_neuro_nhomalt_nfe_female=0;non_neuro_nhomalt_nfe_male=0;non_neuro_nhomalt_nfe_nwe=0;non_neuro_nhomalt_nfe_onf=0;non_neuro_nhomalt_nfe_seu=0;non_neuro_nhomalt_nfe_swe=0;non_neuro_nhomalt_oth=0;non_neuro_nhomalt_oth_female=0;non_neuro_nhomalt_oth_male=0;non_neuro_nhomalt_raw=89;non_neuro_nhomalt_sas=0;non_neuro_nhomalt_sas_female=0;non_neuro_nhomalt_sas_male=0;non_topmed_AC=0;non_topmed_AC_afr=0;non_topmed_AC_afr_female=0;non_topmed_AC_afr_male=0;non_topmed_AC_amr=0;non_topmed_AC_amr_female=0;non_topmed_AC_amr_male=0;non_topmed_AC_asj=0;non_topmed_AC_asj_female=0;non_topmed_AC_asj_male=0;non_topmed_AC_eas=0;non_topmed_AC_eas_female=0;non_topmed_AC_eas_jpn=0;non_topmed_AC_eas_kor=0;non_topmed_AC_eas_male=0;non_topmed_AC_eas_oea=0;non_topmed_AC_female=0;non_topmed_AC_fin=0;non_topmed_AC_fin_female=0;non_topmed_AC_fin_male=0;non_topmed_AC_male=0;non_topmed_AC_nfe=0;non_topmed_AC_nfe_bgr=0;non_topmed_AC_nfe_est=0;non_topmed_AC_nfe_female=0;non_topmed_AC_nfe_male=0;non_topmed_AC_nfe_nwe=0;non_topmed_AC_nfe_onf=0;non_topmed_AC_nfe_seu=0;non_topmed_AC_nfe_swe=0;non_topmed_AC_oth=0;non_topmed_AC_oth_female=0;non_topmed_AC_oth_male=0;non_topmed_AC_raw=218;non_topmed_AC_sas=0;non_topmed_AC_sas_female=0;non_topmed_AC_sas_male=0;non_topmed_AF_raw=4.59334e-02;non_topmed_AN=0;non_topmed_AN_afr=0;non_topmed_AN_afr_female=0;non_topmed_AN_afr_male=0;non_topmed_AN_amr=0;non_topmed_AN_amr_female=0;non_topmed_AN_amr_male=0;non_topmed_AN_asj=0;non_topmed_AN_asj_female=0;non_topmed_AN_asj_male=0;non_topmed_AN_eas=0;non_topmed_AN_eas_female=0;non_topmed_AN_eas_jpn=0;non_topmed_AN_eas_kor=0;non_topmed_AN_eas_male=0;non_topmed_AN_eas_oea=0;non_topmed_AN_female=0;non_topmed_AN_fin=0;non_topmed_AN_fin_female=0;non_topmed_AN_fin_male=0;non_topmed_AN_male=0;non_topmed_AN_nfe=0;non_topmed_AN_nfe_bgr=0;non_topmed_AN_nfe_est=0;non_topmed_AN_nfe_female=0;non_topmed_AN_nfe_male=0;non_topmed_AN_nfe_nwe=0;non_topmed_AN_nfe_onf=0;non_topmed_AN_nfe_seu=0;non_topmed_AN_nfe_swe=0;non_topmed_AN_oth=0;non_topmed_AN_oth_female=0;non_topmed_AN_oth_male=0;non_topmed_AN_raw=4746;non_topmed_AN_sas=0;non_topmed_AN_sas_female=0;non_topmed_AN_sas_male=0;non_topmed_faf95=0.00000e+00;non_topmed_faf95_afr=0.00000e+00;non_topmed_faf95_amr=0.00000e+00;non_topmed_faf95_eas=0.00000e+00;non_topmed_faf95_nfe=0.00000e+00;non_topmed_faf95_sas=0.00000e+00;non_topmed_faf99=0.00000e+00;non_topmed_faf99_afr=0.00000e+00;non_topmed_faf99_amr=0.00000e+00;non_topmed_faf99_eas=0.00000e+00;non_topmed_faf99_nfe=0.00000e+00;non_topmed_faf99_sas=0.00000e+00;non_topmed_nhomalt=0;non_topmed_nhomalt_afr=0;non_topmed_nhomalt_afr_female=0;non_topmed_nhomalt_afr_male=0;non_topmed_nhomalt_amr=0;non_topmed_nhomalt_amr_female=0;non_topmed_nhomalt_amr_male=0;non_topmed_nhomalt_asj=0;non_topmed_nhomalt_asj_female=0;non_topmed_nhomalt_asj_male=0;non_topmed_nhomalt_eas=0;non_topmed_nhomalt_eas_female=0;non_topmed_nhomalt_eas_jpn=0;non_topmed_nhomalt_eas_kor=0;non_topmed_nhomalt_eas_male=0;non_topmed_nhomalt_eas_oea=0;non_topmed_nhomalt_female=0;non_topmed_nhomalt_fin=0;non_topmed_nhomalt_fin_female=0;non_topmed_nhomalt_fin_male=0;non_topmed_nhomalt_male=0;non_topmed_nhomalt_nfe=0;non_topmed_nhomalt_nfe_bgr=0;non_topmed_nhomalt_nfe_est=0;non_topmed_nhomalt_nfe_female=0;non_topmed_nhomalt_nfe_male=0;non_topmed_nhomalt_nfe_nwe=0;non_topmed_nhomalt_nfe_onf=0;non_topmed_nhomalt_nfe_seu=0;non_topmed_nhomalt_nfe_swe=0;non_topmed_nhomalt_oth=0;non_topmed_nhomalt_oth_female=0;non_topmed_nhomalt_oth_male=0;non_topmed_nhomalt_raw=87;non_topmed_nhomalt_sas=0;non_topmed_nhomalt_sas_female=0;non_topmed_nhomalt_sas_male=0;pab_max=1.00000e+00;rf_label=FP;rf_negative_label;rf_tp_probability=8.36542e-01;rf_train;segdup;variant_type=snv;vep=C|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000423562|unprocessed_pseudogene||||||||||rs62635282|1|2165|-1||SNV|1|HGNC|38034||||||||||||||||||||||||||||||||||||||||||,C|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000438504|unprocessed_pseudogene||||||||||rs62635282|1|2165|-1||SNV|1|HGNC|38034|YES|||||||||||||||||||||||||||||||||||||||||,C|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene|2/6||ENST00000450305.2:n.68G>C||68|||||rs62635282|1||1||SNV|1|HGNC|37102||||||||||||||||||||||||||||||||||||||||||,C|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript|1/3||ENST00000456328.2:n.330G>C||330|||||rs62635282|1||1||SNV|1|HGNC|37102|YES|||||||||||||||||||||||||||||||||||||||||,C|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||||||||||rs62635282|1|2206|-1||SNV|1|HGNC|38034||||||||||||||||||||||||||||||||||||||||||,C|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000515242|transcribed_unprocessed_pseudogene|1/3||ENST00000515242.2:n.327G>C||327|||||rs62635282|1||1||SNV|1|HGNC|37102||||||||||||||||||||||||||||||||||||||||||,C|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000518655|transcribed_unprocessed_pseudogene|1/4||ENST00000518655.2:n.325G>C||325|||||rs62635282|1||1||SNV|1|HGNC|37102||||||||||||||||||||||||||||||||||||||||||,C|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000538476|unprocessed_pseudogene||||||||||rs62635282|1|2213|-1||SNV|1|HGNC|38034||||||||||||||||||||||||||||||||||||||||||,C|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000541675|unprocessed_pseudogene||||||||||rs62635282|1|2165|-1||SNV|1|HGNC|38034||||||||||||||||||||||||||||||||||||||||||,C|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00001576075|CTCF_binding_site||||||||||rs62635282|1||||SNV|1||||||||||||||||||||||||||||||||||||||||||||
            # v2 genomes (WGS)
            # chr1	12198	rs62635282	G	C	62.26	AC0;RF	AC=0;AC_afr=0;AC_afr_female=0;AC_afr_male=0;AC_amr=0;AC_amr_female=0;AC_amr_male=0;AC_asj=0;AC_asj_female=0;AC_asj_male=0;AC_eas=0;AC_eas_female=0;AC_eas_male=0;AC_female=0;AC_fin=0;AC_fin_female=0;AC_fin_male=0;AC_male=0;AC_nfe=0;AC_nfe_est=0;AC_nfe_female=0;AC_nfe_male=0;AC_nfe_nwe=0;AC_nfe_onf=0;AC_nfe_seu=0;AC_oth=0;AC_oth_female=0;AC_oth_male=0;AC_raw=4;AF_raw=2.66667e-02;AN=0;AN_afr=0;AN_afr_female=0;AN_afr_male=0;AN_amr=0;AN_amr_female=0;AN_amr_male=0;AN_asj=0;AN_asj_female=0;AN_asj_male=0;AN_eas=0;AN_eas_female=0;AN_eas_male=0;AN_female=0;AN_fin=0;AN_fin_female=0;AN_fin_male=0;AN_male=0;AN_nfe=0;AN_nfe_est=0;AN_nfe_female=0;AN_nfe_male=0;AN_nfe_nwe=0;AN_nfe_onf=0;AN_nfe_seu=0;AN_oth=0;AN_oth_female=0;AN_oth_male=0;AN_raw=150;DP=225;FS=0.00000e+00;InbreedingCoeff=-4.04000e-02;MQ=2.60200e+01;OriginalContig=1;OriginalStart=12198;QD=1.55600e+01;SOR=3.25800e+00;VQSLOD=-9.77700e+01;VQSR_culprit=MQ;ab_hist_alt_bin_freq=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;age_hist_het_bin_freq=0|0|0|0|0|0|0|0|0|0;age_hist_het_n_larger=0;age_hist_het_n_smaller=0;age_hist_hom_bin_freq=0|0|0|0|0|0|0|0|0|0;age_hist_hom_n_larger=0;age_hist_hom_n_smaller=0;allele_type=snv;controls_AC=0;controls_AC_afr=0;controls_AC_afr_female=0;controls_AC_afr_male=0;controls_AC_amr=0;controls_AC_amr_female=0;controls_AC_amr_male=0;controls_AC_asj=0;controls_AC_asj_female=0;controls_AC_asj_male=0;controls_AC_eas=0;controls_AC_eas_female=0;controls_AC_eas_male=0;controls_AC_female=0;controls_AC_fin=0;controls_AC_fin_female=0;controls_AC_fin_male=0;controls_AC_male=0;controls_AC_nfe=0;controls_AC_nfe_est=0;controls_AC_nfe_female=0;controls_AC_nfe_male=0;controls_AC_nfe_nwe=0;controls_AC_nfe_onf=0;controls_AC_nfe_seu=0;controls_AC_oth=0;controls_AC_oth_female=0;controls_AC_oth_male=0;controls_AC_raw=0;controls_AF_raw=0.00000e+00;controls_AN=0;controls_AN_afr=0;controls_AN_afr_female=0;controls_AN_afr_male=0;controls_AN_amr=0;controls_AN_amr_female=0;controls_AN_amr_male=0;controls_AN_asj=0;controls_AN_asj_female=0;controls_AN_asj_male=0;controls_AN_eas=0;controls_AN_eas_female=0;controls_AN_eas_male=0;controls_AN_female=0;controls_AN_fin=0;controls_AN_fin_female=0;controls_AN_fin_male=0;controls_AN_male=0;controls_AN_nfe=0;controls_AN_nfe_est=0;controls_AN_nfe_female=0;controls_AN_nfe_male=0;controls_AN_nfe_nwe=0;controls_AN_nfe_onf=0;controls_AN_nfe_seu=0;controls_AN_oth=0;controls_AN_oth_female=0;controls_AN_oth_male=0;controls_AN_raw=38;controls_faf95=0.00000e+00;controls_faf95_afr=0.00000e+00;controls_faf95_amr=0.00000e+00;controls_faf95_eas=0.00000e+00;controls_faf95_nfe=0.00000e+00;controls_faf99=0.00000e+00;controls_faf99_afr=0.00000e+00;controls_faf99_amr=0.00000e+00;controls_faf99_eas=0.00000e+00;controls_faf99_nfe=0.00000e+00;controls_nhomalt=0;controls_nhomalt_afr=0;controls_nhomalt_afr_female=0;controls_nhomalt_afr_male=0;controls_nhomalt_amr=0;controls_nhomalt_amr_female=0;controls_nhomalt_amr_male=0;controls_nhomalt_asj=0;controls_nhomalt_asj_female=0;controls_nhomalt_asj_male=0;controls_nhomalt_eas=0;controls_nhomalt_eas_female=0;controls_nhomalt_eas_male=0;controls_nhomalt_female=0;controls_nhomalt_fin=0;controls_nhomalt_fin_female=0;controls_nhomalt_fin_male=0;controls_nhomalt_male=0;controls_nhomalt_nfe=0;controls_nhomalt_nfe_est=0;controls_nhomalt_nfe_female=0;controls_nhomalt_nfe_male=0;controls_nhomalt_nfe_nwe=0;controls_nhomalt_nfe_onf=0;controls_nhomalt_nfe_seu=0;controls_nhomalt_oth=0;controls_nhomalt_oth_female=0;controls_nhomalt_oth_male=0;controls_nhomalt_raw=0;dp_hist_all_bin_freq=15708|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;dp_hist_all_n_larger=0;dp_hist_alt_bin_freq=2|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;dp_hist_alt_n_larger=0;faf95=0.00000e+00;faf95_afr=0.00000e+00;faf95_amr=0.00000e+00;faf95_eas=0.00000e+00;faf95_nfe=0.00000e+00;faf99=0.00000e+00;faf99_afr=0.00000e+00;faf99_amr=0.00000e+00;faf99_eas=0.00000e+00;faf99_nfe=0.00000e+00;gq_hist_all_bin_freq=71|4|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;gq_hist_alt_bin_freq=0|2|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;n_alt_alleles=1;nhomalt=0;nhomalt_afr=0;nhomalt_afr_female=0;nhomalt_afr_male=0;nhomalt_amr=0;nhomalt_amr_female=0;nhomalt_amr_male=0;nhomalt_asj=0;nhomalt_asj_female=0;nhomalt_asj_male=0;nhomalt_eas=0;nhomalt_eas_female=0;nhomalt_eas_male=0;nhomalt_female=0;nhomalt_fin=0;nhomalt_fin_female=0;nhomalt_fin_male=0;nhomalt_male=0;nhomalt_nfe=0;nhomalt_nfe_est=0;nhomalt_nfe_female=0;nhomalt_nfe_male=0;nhomalt_nfe_nwe=0;nhomalt_nfe_onf=0;nhomalt_nfe_seu=0;nhomalt_oth=0;nhomalt_oth_female=0;nhomalt_oth_male=0;nhomalt_raw=2;non_neuro_AC=0;non_neuro_AC_afr=0;non_neuro_AC_afr_female=0;non_neuro_AC_afr_male=0;non_neuro_AC_amr=0;non_neuro_AC_amr_female=0;non_neuro_AC_amr_male=0;non_neuro_AC_asj=0;non_neuro_AC_asj_female=0;non_neuro_AC_asj_male=0;non_neuro_AC_eas=0;non_neuro_AC_eas_female=0;non_neuro_AC_eas_male=0;non_neuro_AC_female=0;non_neuro_AC_fin=0;non_neuro_AC_fin_female=0;non_neuro_AC_fin_male=0;non_neuro_AC_male=0;non_neuro_AC_nfe=0;non_neuro_AC_nfe_est=0;non_neuro_AC_nfe_female=0;non_neuro_AC_nfe_male=0;non_neuro_AC_nfe_nwe=0;non_neuro_AC_nfe_onf=0;non_neuro_AC_nfe_seu=0;non_neuro_AC_oth=0;non_neuro_AC_oth_female=0;non_neuro_AC_oth_male=0;non_neuro_AC_raw=2;non_neuro_AF_raw=1.88679e-02;non_neuro_AN=0;non_neuro_AN_afr=0;non_neuro_AN_afr_female=0;non_neuro_AN_afr_male=0;non_neuro_AN_amr=0;non_neuro_AN_amr_female=0;non_neuro_AN_amr_male=0;non_neuro_AN_asj=0;non_neuro_AN_asj_female=0;non_neuro_AN_asj_male=0;non_neuro_AN_eas=0;non_neuro_AN_eas_female=0;non_neuro_AN_eas_male=0;non_neuro_AN_female=0;non_neuro_AN_fin=0;non_neuro_AN_fin_female=0;non_neuro_AN_fin_male=0;non_neuro_AN_male=0;non_neuro_AN_nfe=0;non_neuro_AN_nfe_est=0;non_neuro_AN_nfe_female=0;non_neuro_AN_nfe_male=0;non_neuro_AN_nfe_nwe=0;non_neuro_AN_nfe_onf=0;non_neuro_AN_nfe_seu=0;non_neuro_AN_oth=0;non_neuro_AN_oth_female=0;non_neuro_AN_oth_male=0;non_neuro_AN_raw=106;non_neuro_faf95=0.00000e+00;non_neuro_faf95_afr=0.00000e+00;non_neuro_faf95_amr=0.00000e+00;non_neuro_faf95_eas=0.00000e+00;non_neuro_faf95_nfe=0.00000e+00;non_neuro_faf99=0.00000e+00;non_neuro_faf99_afr=0.00000e+00;non_neuro_faf99_amr=0.00000e+00;non_neuro_faf99_eas=0.00000e+00;non_neuro_faf99_nfe=0.00000e+00;non_neuro_nhomalt=0;non_neuro_nhomalt_afr=0;non_neuro_nhomalt_afr_female=0;non_neuro_nhomalt_afr_male=0;non_neuro_nhomalt_amr=0;non_neuro_nhomalt_amr_female=0;non_neuro_nhomalt_amr_male=0;non_neuro_nhomalt_asj=0;non_neuro_nhomalt_asj_female=0;non_neuro_nhomalt_asj_male=0;non_neuro_nhomalt_eas=0;non_neuro_nhomalt_eas_female=0;non_neuro_nhomalt_eas_male=0;non_neuro_nhomalt_female=0;non_neuro_nhomalt_fin=0;non_neuro_nhomalt_fin_female=0;non_neuro_nhomalt_fin_male=0;non_neuro_nhomalt_male=0;non_neuro_nhomalt_nfe=0;non_neuro_nhomalt_nfe_est=0;non_neuro_nhomalt_nfe_female=0;non_neuro_nhomalt_nfe_male=0;non_neuro_nhomalt_nfe_nwe=0;non_neuro_nhomalt_nfe_onf=0;non_neuro_nhomalt_nfe_seu=0;non_neuro_nhomalt_oth=0;non_neuro_nhomalt_oth_female=0;non_neuro_nhomalt_oth_male=0;non_neuro_nhomalt_raw=1;non_topmed_AC=0;non_topmed_AC_afr=0;non_topmed_AC_afr_female=0;non_topmed_AC_afr_male=0;non_topmed_AC_amr=0;non_topmed_AC_amr_female=0;non_topmed_AC_amr_male=0;non_topmed_AC_asj=0;non_topmed_AC_asj_female=0;non_topmed_AC_asj_male=0;non_topmed_AC_eas=0;non_topmed_AC_eas_female=0;non_topmed_AC_eas_male=0;non_topmed_AC_female=0;non_topmed_AC_fin=0;non_topmed_AC_fin_female=0;non_topmed_AC_fin_male=0;non_topmed_AC_male=0;non_topmed_AC_nfe=0;non_topmed_AC_nfe_est=0;non_topmed_AC_nfe_female=0;non_topmed_AC_nfe_male=0;non_topmed_AC_nfe_nwe=0;non_topmed_AC_nfe_onf=0;non_topmed_AC_nfe_seu=0;non_topmed_AC_oth=0;non_topmed_AC_oth_female=0;non_topmed_AC_oth_male=0;non_topmed_AC_raw=4;non_topmed_AF_raw=2.81690e-02;non_topmed_AN=0;non_topmed_AN_afr=0;non_topmed_AN_afr_female=0;non_topmed_AN_afr_male=0;non_topmed_AN_amr=0;non_topmed_AN_amr_female=0;non_topmed_AN_amr_male=0;non_topmed_AN_asj=0;non_topmed_AN_asj_female=0;non_topmed_AN_asj_male=0;non_topmed_AN_eas=0;non_topmed_AN_eas_female=0;non_topmed_AN_eas_male=0;non_topmed_AN_female=0;non_topmed_AN_fin=0;non_topmed_AN_fin_female=0;non_topmed_AN_fin_male=0;non_topmed_AN_male=0;non_topmed_AN_nfe=0;non_topmed_AN_nfe_est=0;non_topmed_AN_nfe_female=0;non_topmed_AN_nfe_male=0;non_topmed_AN_nfe_nwe=0;non_topmed_AN_nfe_onf=0;non_topmed_AN_nfe_seu=0;non_topmed_AN_oth=0;non_topmed_AN_oth_female=0;non_topmed_AN_oth_male=0;non_topmed_AN_raw=142;non_topmed_faf95=0.00000e+00;non_topmed_faf95_afr=0.00000e+00;non_topmed_faf95_amr=0.00000e+00;non_topmed_faf95_eas=0.00000e+00;non_topmed_faf95_nfe=0.00000e+00;non_topmed_faf99=0.00000e+00;non_topmed_faf99_afr=0.00000e+00;non_topmed_faf99_amr=0.00000e+00;non_topmed_faf99_eas=0.00000e+00;non_topmed_faf99_nfe=0.00000e+00;non_topmed_nhomalt=0;non_topmed_nhomalt_afr=0;non_topmed_nhomalt_afr_female=0;non_topmed_nhomalt_afr_male=0;non_topmed_nhomalt_amr=0;non_topmed_nhomalt_amr_female=0;non_topmed_nhomalt_amr_male=0;non_topmed_nhomalt_asj=0;non_topmed_nhomalt_asj_female=0;non_topmed_nhomalt_asj_male=0;non_topmed_nhomalt_eas=0;non_topmed_nhomalt_eas_female=0;non_topmed_nhomalt_eas_male=0;non_topmed_nhomalt_female=0;non_topmed_nhomalt_fin=0;non_topmed_nhomalt_fin_female=0;non_topmed_nhomalt_fin_male=0;non_topmed_nhomalt_male=0;non_topmed_nhomalt_nfe=0;non_topmed_nhomalt_nfe_est=0;non_topmed_nhomalt_nfe_female=0;non_topmed_nhomalt_nfe_male=0;non_topmed_nhomalt_nfe_nwe=0;non_topmed_nhomalt_nfe_onf=0;non_topmed_nhomalt_nfe_seu=0;non_topmed_nhomalt_oth=0;non_topmed_nhomalt_oth_female=0;non_topmed_nhomalt_oth_male=0;non_topmed_nhomalt_raw=2;rf_label=FP;rf_negative_label;rf_tp_probability=6.55868e-02;rf_train;segdup;variant_type=snv;vep=C|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000423562|unprocessed_pseudogene||||||||||rs62635282|1|2165|-1||SNV|1|HGNC|38034||||||||||||||||||||||||||||||||||||||||||,C|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000438504|unprocessed_pseudogene||||||||||rs62635282|1|2165|-1||SNV|1|HGNC|38034|YES|||||||||||||||||||||||||||||||||||||||||,C|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene|2/6||ENST00000450305.2:n.68G>C||68|||||rs62635282|1||1||SNV|1|HGNC|37102||||||||||||||||||||||||||||||||||||||||||,C|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript|1/3||ENST00000456328.2:n.330G>C||330|||||rs62635282|1||1||SNV|1|HGNC|37102|YES|||||||||||||||||||||||||||||||||||||||||,C|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||||||||||rs62635282|1|2206|-1||SNV|1|HGNC|38034||||||||||||||||||||||||||||||||||||||||||,C|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000515242|transcribed_unprocessed_pseudogene|1/3||ENST00000515242.2:n.327G>C||327|||||rs62635282|1||1||SNV|1|HGNC|37102||||||||||||||||||||||||||||||||||||||||||,C|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000518655|transcribed_unprocessed_pseudogene|1/4||ENST00000518655.2:n.325G>C||325|||||rs62635282|1||1||SNV|1|HGNC|37102||||||||||||||||||||||||||||||||||||||||||,C|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000538476|unprocessed_pseudogene||||||||||rs62635282|1|2213|-1||SNV|1|HGNC|38034||||||||||||||||||||||||||||||||||||||||||,C|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000541675|unprocessed_pseudogene||||||||||rs62635282|1|2165|-1||SNV|1|HGNC|38034||||||||||||||||||||||||||||||||||||||||||,C|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00001576075|CTCF_binding_site||||||||||rs62635282|1||||SNV|1||||||||||||||||||||||||||||||||||||||||||||
            # v3 genomes (WGS)
            # chr1	12198	rs62635282	G	C	.	AC0;AS_VQSR	AC=0;AN=0;AC_non_v2_XX=0;AN_non_v2_XX=0;nhomalt_non_v2_XX=0;AC_non_cancer_fin_XX=0;AN_non_cancer_fin_XX=0;nhomalt_non_cancer_fin_XX=0;AC_non_neuro_nfe=0;AN_non_neuro_nfe=0;nhomalt_non_neuro_nfe=0;AC_non_neuro_afr_XY=0;AN_non_neuro_afr_XY=0;nhomalt_non_neuro_afr_XY=0;AC_non_neuro_nfe_XY=0;AN_non_neuro_nfe_XY=0;nhomalt_non_neuro_nfe_XY=0;AC_controls_and_biobanks_eas_XY=0;AN_controls_and_biobanks_eas_XY=0;nhomalt_controls_and_biobanks_eas_XY=0;AC_non_neuro_sas_XX=0;AN_non_neuro_sas_XX=0;nhomalt_non_neuro_sas_XX=0;AC_non_v2=0;AN_non_v2=0;nhomalt_non_v2=0;AC_non_topmed_nfe_XX=0;AN_non_topmed_nfe_XX=0;nhomalt_non_topmed_nfe_XX=0;AC_non_v2_mid=0;AN_non_v2_mid=0;nhomalt_non_v2_mid=0;AC_non_topmed_sas=0;AN_non_topmed_sas=0;nhomalt_non_topmed_sas=0;AC_non_cancer_eas_XX=0;AN_non_cancer_eas_XX=0;nhomalt_non_cancer_eas_XX=0;AC_amr_XY=0;AN_amr_XY=0;nhomalt_amr_XY=0;AC_non_v2_nfe_XX=0;AN_non_v2_nfe_XX=0;nhomalt_non_v2_nfe_XX=0;AC_controls_and_biobanks_XY=0;AN_controls_and_biobanks_XY=0;nhomalt_controls_and_biobanks_XY=0;AC_non_neuro_asj_XY=0;AN_non_neuro_asj_XY=0;nhomalt_non_neuro_asj_XY=0;AC_oth=0;AN_oth=0;nhomalt_oth=0;AC_non_topmed_mid_XY=0;AN_non_topmed_mid_XY=0;nhomalt_non_topmed_mid_XY=0;AC_non_cancer_asj_XX=0;AN_non_cancer_asj_XX=0;nhomalt_non_cancer_asj_XX=0;AC_sas_XY=0;AN_sas_XY=0;nhomalt_sas_XY=0;AC_non_neuro_fin=0;AN_non_neuro_fin=0;nhomalt_non_neuro_fin=0;AC_non_topmed_amr_XY=0;AN_non_topmed_amr_XY=0;nhomalt_non_topmed_amr_XY=0;AC_non_neuro_XX=0;AN_non_neuro_XX=0;nhomalt_non_neuro_XX=0;AC_fin_XX=0;AN_fin_XX=0;nhomalt_fin_XX=0;AC_controls_and_biobanks_asj_XX=0;AN_controls_and_biobanks_asj_XX=0;nhomalt_controls_and_biobanks_asj_XX=0;AC_non_v2_raw=4;AN_non_v2_raw=302;AF_non_v2_raw=0.0132450;nhomalt_non_v2_raw=2;AC_non_v2_asj=0;AN_non_v2_asj=0;nhomalt_non_v2_asj=0;AC_nfe_XX=0;AN_nfe_XX=0;nhomalt_nfe_XX=0;AC_controls_and_biobanks_raw=0;AN_controls_and_biobanks_raw=36;AF_controls_and_biobanks_raw=0.00000;nhomalt_controls_and_biobanks_raw=0;AC_controls_and_biobanks_ami=0;AN_controls_and_biobanks_ami=0;nhomalt_controls_and_biobanks_ami=0;AC_non_topmed_eas=0;AN_non_topmed_eas=0;nhomalt_non_topmed_eas=0;AC_non_v2_amr=0;AN_non_v2_amr=0;nhomalt_non_v2_amr=0;AC_non_neuro_sas=0;AN_non_neuro_sas=0;nhomalt_non_neuro_sas=0;AC_non_cancer_fin_XY=0;AN_non_cancer_fin_XY=0;nhomalt_non_cancer_fin_XY=0;AC_non_cancer_nfe_XY=0;AN_non_cancer_nfe_XY=0;nhomalt_non_cancer_nfe_XY=0;AC_non_v2_oth=0;AN_non_v2_oth=0;nhomalt_non_v2_oth=0;AC_ami=0;AN_ami=0;nhomalt_ami=0;AC_non_cancer_XY=0;AN_non_cancer_XY=0;nhomalt_non_cancer_XY=0;AC_non_v2_sas=0;AN_non_v2_sas=0;nhomalt_non_v2_sas=0;AC_non_topmed_afr_XX=0;AN_non_topmed_afr_XX=0;nhomalt_non_topmed_afr_XX=0;AC_sas=0;AN_sas=0;nhomalt_sas=0;AC_non_neuro_nfe_XX=0;AN_non_neuro_nfe_XX=0;nhomalt_non_neuro_nfe_XX=0;AC_non_topmed_ami_XX=0;AN_non_topmed_ami_XX=0;nhomalt_non_topmed_ami_XX=0;AC_ami_XY=0;AN_ami_XY=0;nhomalt_ami_XY=0;AC_oth_XX=0;AN_oth_XX=0;nhomalt_oth_XX=0;AC_non_cancer_eas=0;AN_non_cancer_eas=0;nhomalt_non_cancer_eas=0;AC_non_topmed_XY=0;AN_non_topmed_XY=0;nhomalt_non_topmed_XY=0;AC_non_v2_ami=0;AN_non_v2_ami=0;nhomalt_non_v2_ami=0;AC_non_neuro=0;AN_non_neuro=0;nhomalt_non_neuro=0;AC_amr_XX=0;AN_amr_XX=0;nhomalt_amr_XX=0;AC_controls_and_biobanks_nfe_XY=0;AN_controls_and_biobanks_nfe_XY=0;nhomalt_controls_and_biobanks_nfe_XY=0;AC_controls_and_biobanks_eas=0;AN_controls_and_biobanks_eas=0;nhomalt_controls_and_biobanks_eas=0;AC_XX=0;AN_XX=0;nhomalt_XX=0;AC_non_cancer_oth_XY=0;AN_non_cancer_oth_XY=0;nhomalt_non_cancer_oth_XY=0;AC_non_v2_XY=0;AN_non_v2_XY=0;nhomalt_non_v2_XY=0;AC_non_topmed_amr_XX=0;AN_non_topmed_amr_XX=0;nhomalt_non_topmed_amr_XX=0;AC_fin=0;AN_fin=0;nhomalt_fin=0;AC_controls_and_biobanks_nfe_XX=0;AN_controls_and_biobanks_nfe_XX=0;nhomalt_controls_and_biobanks_nfe_XX=0;AC_controls_and_biobanks_afr=0;AN_controls_and_biobanks_afr=0;nhomalt_controls_and_biobanks_afr=0;AC_asj_XX=0;AN_asj_XX=0;nhomalt_asj_XX=0;AC_non_topmed_mid=0;AN_non_topmed_mid=0;nhomalt_non_topmed_mid=0;AC_non_cancer_sas_XY=0;AN_non_cancer_sas_XY=0;nhomalt_non_cancer_sas_XY=0;AC_sas_XX=0;AN_sas_XX=0;nhomalt_sas_XX=0;AC_non_topmed=0;AN_non_topmed=0;nhomalt_non_topmed=0;AC_non_v2_oth_XX=0;AN_non_v2_oth_XX=0;nhomalt_non_v2_oth_XX=0;AC_non_neuro_ami_XY=0;AN_non_neuro_ami_XY=0;nhomalt_non_neuro_ami_XY=0;AC_controls_and_biobanks_afr_XY=0;AN_controls_and_biobanks_afr_XY=0;nhomalt_controls_and_biobanks_afr_XY=0;AC_controls_and_biobanks_amr_XX=0;AN_controls_and_biobanks_amr_XX=0;nhomalt_controls_and_biobanks_amr_XX=0;AC_non_topmed_amr=0;AN_non_topmed_amr=0;nhomalt_non_topmed_amr=0;AC_controls_and_biobanks_sas_XX=0;AN_controls_and_biobanks_sas_XX=0;nhomalt_controls_and_biobanks_sas_XX=0;AC_controls_and_biobanks_amr=0;AN_controls_and_biobanks_amr=0;nhomalt_controls_and_biobanks_amr=0;AC_non_neuro_fin_XX=0;AN_non_neuro_fin_XX=0;nhomalt_non_neuro_fin_XX=0;AC_non_cancer_raw=6;AN_non_cancer_raw=444;AF_non_cancer_raw=0.0135135;nhomalt_non_cancer_raw=3;AC_non_neuro_mid=0;AN_non_neuro_mid=0;nhomalt_non_neuro_mid=0;AC_non_v2_asj_XY=0;AN_non_v2_asj_XY=0;nhomalt_non_v2_asj_XY=0;AC_non_v2_afr=0;AN_non_v2_afr=0;nhomalt_non_v2_afr=0;AC_non_neuro_fin_XY=0;AN_non_neuro_fin_XY=0;nhomalt_non_neuro_fin_XY=0;AC_non_cancer_afr=0;AN_non_cancer_afr=0;nhomalt_non_cancer_afr=0;AC_non_topmed_sas_XY=0;AN_non_topmed_sas_XY=0;nhomalt_non_topmed_sas_XY=0;AC_mid_XY=0;AN_mid_XY=0;nhomalt_mid_XY=0;AC_non_v2_oth_XY=0;AN_non_v2_oth_XY=0;nhomalt_non_v2_oth_XY=0;AC_controls_and_biobanks_fin=0;AN_controls_and_biobanks_fin=0;nhomalt_controls_and_biobanks_fin=0;AC_non_neuro_eas_XY=0;AN_non_neuro_eas_XY=0;nhomalt_non_neuro_eas_XY=0;AC_non_topmed_eas_XX=0;AN_non_topmed_eas_XX=0;nhomalt_non_topmed_eas_XX=0;AC_non_v2_afr_XX=0;AN_non_v2_afr_XX=0;nhomalt_non_v2_afr_XX=0;AC_non_neuro_amr_XX=0;AN_non_neuro_amr_XX=0;nhomalt_non_neuro_amr_XX=0;AC_non_cancer_ami=0;AN_non_cancer_ami=0;nhomalt_non_cancer_ami=0;AC_XY=0;AN_XY=0;nhomalt_XY=0;AC_non_topmed_asj_XX=0;AN_non_topmed_asj_XX=0;nhomalt_non_topmed_asj_XX=0;AC_non_topmed_eas_XY=0;AN_non_topmed_eas_XY=0;nhomalt_non_topmed_eas_XY=0;AC_non_v2_eas_XY=0;AN_non_v2_eas_XY=0;nhomalt_non_v2_eas_XY=0;AC_eas=0;AN_eas=0;nhomalt_eas=0;AC_asj_XY=0;AN_asj_XY=0;nhomalt_asj_XY=0;AC_non_v2_eas_XX=0;AN_non_v2_eas_XX=0;nhomalt_non_v2_eas_XX=0;AC_controls_and_biobanks_mid_XY=0;AN_controls_and_biobanks_mid_XY=0;nhomalt_controls_and_biobanks_mid_XY=0;AC_fin_XY=0;AN_fin_XY=0;nhomalt_fin_XY=0;AC_non_topmed_nfe=0;AN_non_topmed_nfe=0;nhomalt_non_topmed_nfe=0;AC_amr=0;AN_amr=0;nhomalt_amr=0;AC_non_neuro_ami=0;AN_non_neuro_ami=0;nhomalt_non_neuro_ami=0;AC_non_cancer_nfe_XX=0;AN_non_cancer_nfe_XX=0;nhomalt_non_cancer_nfe_XX=0;AC_non_cancer_mid=0;AN_non_cancer_mid=0;nhomalt_non_cancer_mid=0;AC_non_v2_mid_XY=0;AN_non_v2_mid_XY=0;nhomalt_non_v2_mid_XY=0;AC_controls_and_biobanks_amr_XY=0;AN_controls_and_biobanks_amr_XY=0;nhomalt_controls_and_biobanks_amr_XY=0;AC_non_cancer_ami_XY=0;AN_non_cancer_ami_XY=0;nhomalt_non_cancer_ami_XY=0;AC_non_neuro_asj_XX=0;AN_non_neuro_asj_XX=0;nhomalt_non_neuro_asj_XX=0;AC_afr=0;AN_afr=0;nhomalt_afr=0;AC_non_v2_sas_XX=0;AN_non_v2_sas_XX=0;nhomalt_non_v2_sas_XX=0;AC_non_neuro_afr_XX=0;AN_non_neuro_afr_XX=0;nhomalt_non_neuro_afr_XX=0;AC_non_cancer_sas=0;AN_non_cancer_sas=0;nhomalt_non_cancer_sas=0;AC_non_topmed_fin=0;AN_non_topmed_fin=0;nhomalt_non_topmed_fin=0;AC_non_cancer_asj_XY=0;AN_non_cancer_asj_XY=0;nhomalt_non_cancer_asj_XY=0;AC_non_cancer_mid_XY=0;AN_non_cancer_mid_XY=0;nhomalt_non_cancer_mid_XY=0;AC_raw=6;AN_raw=444;AF_raw=0.0135135;nhomalt_raw=3;AC_non_topmed_XX=0;AN_non_topmed_XX=0;nhomalt_non_topmed_XX=0;AC_ami_XX=0;AN_ami_XX=0;nhomalt_ami_XX=0;AC_eas_XY=0;AN_eas_XY=0;nhomalt_eas_XY=0;AC_controls_and_biobanks_mid=0;AN_controls_and_biobanks_mid=0;nhomalt_controls_and_biobanks_mid=0;AC_non_v2_nfe_XY=0;AN_non_v2_nfe_XY=0;nhomalt_non_v2_nfe_XY=0;AC_controls_and_biobanks_sas=0;AN_controls_and_biobanks_sas=0;nhomalt_controls_and_biobanks_sas=0;AC_non_v2_eas=0;AN_non_v2_eas=0;nhomalt_non_v2_eas=0;AC_mid=0;AN_mid=0;nhomalt_mid=0;AC_oth_XY=0;AN_oth_XY=0;nhomalt_oth_XY=0;AC_non_cancer_nfe=0;AN_non_cancer_nfe=0;nhomalt_non_cancer_nfe=0;AC_non_neuro_eas_XX=0;AN_non_neuro_eas_XX=0;nhomalt_non_neuro_eas_XX=0;AC_non_neuro_sas_XY=0;AN_non_neuro_sas_XY=0;nhomalt_non_neuro_sas_XY=0;AC_non_cancer_ami_XX=0;AN_non_cancer_ami_XX=0;nhomalt_non_cancer_ami_XX=0;AC_mid_XX=0;AN_mid_XX=0;nhomalt_mid_XX=0;AC_non_topmed_asj=0;AN_non_topmed_asj=0;nhomalt_non_topmed_asj=0;AC_non_v2_asj_XX=0;AN_non_v2_asj_XX=0;nhomalt_non_v2_asj_XX=0;nhomalt=0;AC_non_v2_amr_XY=0;AN_non_v2_amr_XY=0;nhomalt_non_v2_amr_XY=0;AC_non_cancer_amr_XX=0;AN_non_cancer_amr_XX=0;nhomalt_non_cancer_amr_XX=0;AC_controls_and_biobanks_afr_XX=0;AN_controls_and_biobanks_afr_XX=0;nhomalt_controls_and_biobanks_afr_XX=0;AC_asj=0;AN_asj=0;nhomalt_asj=0;AC_non_topmed_asj_XY=0;AN_non_topmed_asj_XY=0;nhomalt_non_topmed_asj_XY=0;AC_non_v2_fin_XX=0;AN_non_v2_fin_XX=0;nhomalt_non_v2_fin_XX=0;AC_non_topmed_ami=0;AN_non_topmed_ami=0;nhomalt_non_topmed_ami=0;AC_controls_and_biobanks_eas_XX=0;AN_controls_and_biobanks_eas_XX=0;nhomalt_controls_and_biobanks_eas_XX=0;AC_controls_and_biobanks_fin_XX=0;AN_controls_and_biobanks_fin_XX=0;nhomalt_controls_and_biobanks_fin_XX=0;AC_non_topmed_raw=2;AN_non_topmed_raw=178;AF_non_topmed_raw=0.0112360;nhomalt_non_topmed_raw=1;AC_non_cancer_eas_XY=0;AN_non_cancer_eas_XY=0;nhomalt_non_cancer_eas_XY=0;AC_non_cancer=0;AN_non_cancer=0;nhomalt_non_cancer=0;AC_controls_and_biobanks_ami_XY=0;AN_controls_and_biobanks_ami_XY=0;nhomalt_controls_and_biobanks_ami_XY=0;AC_controls_and_biobanks_mid_XX=0;AN_controls_and_biobanks_mid_XX=0;nhomalt_controls_and_biobanks_mid_XX=0;AC_non_v2_afr_XY=0;AN_non_v2_afr_XY=0;nhomalt_non_v2_afr_XY=0;AC_non_v2_sas_XY=0;AN_non_v2_sas_XY=0;nhomalt_non_v2_sas_XY=0;AC_non_v2_fin=0;AN_non_v2_fin=0;nhomalt_non_v2_fin=0;AC_non_neuro_oth=0;AN_non_neuro_oth=0;nhomalt_non_neuro_oth=0;AC_non_cancer_sas_XX=0;AN_non_cancer_sas_XX=0;nhomalt_non_cancer_sas_XX=0;AC_non_neuro_asj=0;AN_non_neuro_asj=0;nhomalt_non_neuro_asj=0;AC_non_topmed_afr=0;AN_non_topmed_afr=0;nhomalt_non_topmed_afr=0;AC_non_topmed_afr_XY=0;AN_non_topmed_afr_XY=0;nhomalt_non_topmed_afr_XY=0;AC_non_neuro_eas=0;AN_non_neuro_eas=0;nhomalt_non_neuro_eas=0;AC_afr_XX=0;AN_afr_XX=0;nhomalt_afr_XX=0;AC_non_neuro_mid_XY=0;AN_non_neuro_mid_XY=0;nhomalt_non_neuro_mid_XY=0;AC_non_topmed_fin_XX=0;AN_non_topmed_fin_XX=0;nhomalt_non_topmed_fin_XX=0;AC_non_cancer_amr=0;AN_non_cancer_amr=0;nhomalt_non_cancer_amr=0;AC_non_v2_ami_XX=0;AN_non_v2_ami_XX=0;nhomalt_non_v2_ami_XX=0;AC_afr_XY=0;AN_afr_XY=0;nhomalt_afr_XY=0;AC_non_v2_mid_XX=0;AN_non_v2_mid_XX=0;nhomalt_non_v2_mid_XX=0;AC_non_topmed_fin_XY=0;AN_non_topmed_fin_XY=0;nhomalt_non_topmed_fin_XY=0;AC_non_neuro_amr_XY=0;AN_non_neuro_amr_XY=0;nhomalt_non_neuro_amr_XY=0;AC_non_topmed_mid_XX=0;AN_non_topmed_mid_XX=0;nhomalt_non_topmed_mid_XX=0;AC_controls_and_biobanks_asj_XY=0;AN_controls_and_biobanks_asj_XY=0;nhomalt_controls_and_biobanks_asj_XY=0;AC_non_v2_fin_XY=0;AN_non_v2_fin_XY=0;nhomalt_non_v2_fin_XY=0;AC_controls_and_biobanks_ami_XX=0;AN_controls_and_biobanks_ami_XX=0;nhomalt_controls_and_biobanks_ami_XX=0;AC_eas_XX=0;AN_eas_XX=0;nhomalt_eas_XX=0;AC_non_cancer_amr_XY=0;AN_non_cancer_amr_XY=0;nhomalt_non_cancer_amr_XY=0;AC_non_neuro_ami_XX=0;AN_non_neuro_ami_XX=0;nhomalt_non_neuro_ami_XX=0;AC_controls_and_biobanks=0;AN_controls_and_biobanks=0;nhomalt_controls_and_biobanks=0;AC_controls_and_biobanks_oth=0;AN_controls_and_biobanks_oth=0;nhomalt_controls_and_biobanks_oth=0;AC_nfe_XY=0;AN_nfe_XY=0;nhomalt_nfe_XY=0;AC_non_cancer_afr_XX=0;AN_non_cancer_afr_XX=0;nhomalt_non_cancer_afr_XX=0;AC_controls_and_biobanks_sas_XY=0;AN_controls_and_biobanks_sas_XY=0;nhomalt_controls_and_biobanks_sas_XY=0;AC_non_cancer_oth=0;AN_non_cancer_oth=0;nhomalt_non_cancer_oth=0;AC_non_topmed_oth=0;AN_non_topmed_oth=0;nhomalt_non_topmed_oth=0;AC_non_topmed_nfe_XY=0;AN_non_topmed_nfe_XY=0;nhomalt_non_topmed_nfe_XY=0;AC_non_topmed_sas_XX=0;AN_non_topmed_sas_XX=0;nhomalt_non_topmed_sas_XX=0;AC_non_v2_nfe=0;AN_non_v2_nfe=0;nhomalt_non_v2_nfe=0;AC_non_topmed_oth_XX=0;AN_non_topmed_oth_XX=0;nhomalt_non_topmed_oth_XX=0;AC_non_cancer_mid_XX=0;AN_non_cancer_mid_XX=0;nhomalt_non_cancer_mid_XX=0;AC_controls_and_biobanks_nfe=0;AN_controls_and_biobanks_nfe=0;nhomalt_controls_and_biobanks_nfe=0;AC_controls_and_biobanks_oth_XY=0;AN_controls_and_biobanks_oth_XY=0;nhomalt_controls_and_biobanks_oth_XY=0;AC_controls_and_biobanks_fin_XY=0;AN_controls_and_biobanks_fin_XY=0;nhomalt_controls_and_biobanks_fin_XY=0;AC_non_v2_amr_XX=0;AN_non_v2_amr_XX=0;nhomalt_non_v2_amr_XX=0;AC_non_cancer_asj=0;AN_non_cancer_asj=0;nhomalt_non_cancer_asj=0;AC_non_cancer_oth_XX=0;AN_non_cancer_oth_XX=0;nhomalt_non_cancer_oth_XX=0;AC_non_neuro_amr=0;AN_non_neuro_amr=0;nhomalt_non_neuro_amr=0;AC_non_cancer_XX=0;AN_non_cancer_XX=0;nhomalt_non_cancer_XX=0;AC_non_v2_ami_XY=0;AN_non_v2_ami_XY=0;nhomalt_non_v2_ami_XY=0;AC_non_neuro_raw=4;AN_non_neuro_raw=382;AF_non_neuro_raw=0.0104712;nhomalt_non_neuro_raw=2;AC_non_neuro_afr=0;AN_non_neuro_afr=0;nhomalt_non_neuro_afr=0;AC_non_topmed_ami_XY=0;AN_non_topmed_ami_XY=0;nhomalt_non_topmed_ami_XY=0;AC_non_neuro_oth_XY=0;AN_non_neuro_oth_XY=0;nhomalt_non_neuro_oth_XY=0;AC_non_neuro_oth_XX=0;AN_non_neuro_oth_XX=0;nhomalt_non_neuro_oth_XX=0;AC_controls_and_biobanks_XX=0;AN_controls_and_biobanks_XX=0;nhomalt_controls_and_biobanks_XX=0;AC_non_cancer_afr_XY=0;AN_non_cancer_afr_XY=0;nhomalt_non_cancer_afr_XY=0;AC_non_cancer_fin=0;AN_non_cancer_fin=0;nhomalt_non_cancer_fin=0;AC_controls_and_biobanks_asj=0;AN_controls_and_biobanks_asj=0;nhomalt_controls_and_biobanks_asj=0;AC_non_topmed_oth_XY=0;AN_non_topmed_oth_XY=0;nhomalt_non_topmed_oth_XY=0;AC_non_neuro_mid_XX=0;AN_non_neuro_mid_XX=0;nhomalt_non_neuro_mid_XX=0;AC_controls_and_biobanks_oth_XX=0;AN_controls_and_biobanks_oth_XX=0;nhomalt_controls_and_biobanks_oth_XX=0;AC_non_neuro_XY=0;AN_non_neuro_XY=0;nhomalt_non_neuro_XY=0;AC_nfe=0;AN_nfe=0;nhomalt_nfe=0;faf95_sas=0.00000;faf99_sas=0.00000;faf95_eas=0.00000;faf99_eas=0.00000;faf95_amr=0.00000;faf99_amr=0.00000;faf95_afr=0.00000;faf99_afr=0.00000;faf95=0.00000;faf99=0.00000;faf95_nfe=0.00000;faf99_nfe=0.00000;age_hist_het_bin_freq=0|0|0|0|0|0|0|0|0|0;age_hist_het_n_smaller=0;age_hist_het_n_larger=0;age_hist_hom_bin_freq=0|0|0|0|0|0|0|0|0|0;age_hist_hom_n_smaller=0;age_hist_hom_n_larger=0;FS=.;MQ=24.5080;MQRankSum=-0.736000;QUALapprox=368;QD=26.2857;ReadPosRankSum=-0.736000;VarDP=14;AS_FS=.;AS_MQ=24.5080;AS_MQRankSum=-0.736000;AS_pab_max=1.00000;AS_QUALapprox=368;AS_QD=26.2857;AS_ReadPosRankSum=-0.736000;AS_SB_TABLE=1,0|13,0;AS_SOR=3.91202;InbreedingCoeff=1.00000;AS_culprit=AS_SOR;AS_VQSLOD=-13.8706;allele_type=snv;n_alt_alleles=1;variant_type=snv;segdup;gq_hist_alt_bin_freq=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;gq_hist_all_bin_freq=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;dp_hist_alt_bin_freq=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;dp_hist_alt_n_larger=0;dp_hist_all_bin_freq=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;dp_hist_all_n_larger=0;ab_hist_alt_bin_freq=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;cadd_raw_score=0.853610;cadd_phred=9.93800;vep=C|non_coding_transcript_exon_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene|2/6||ENST00000450305.2:n.68G>C||68|||||1||1|SNV||HGNC|HGNC:37102|||||||||||||||||||||,C|non_coding_transcript_exon_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript|1/3||ENST00000456328.2:n.330G>C||330|||||1||1|SNV||HGNC|HGNC:37102|YES|1|||||||||||||||||||,C|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||||||||||1|2206|-1|SNV||HGNC|HGNC:38034|YES||||||||||||||||||||,C|downstream_gene_variant|MODIFIER|WASH7P|653635|Transcript|NR_024540.1|transcribed_pseudogene||||||||||1|2164|-1|SNV||EntrezGene|HGNC:38034|YES||||||||||||||||||||,C|non_coding_transcript_exon_variant|MODIFIER|DDX11L1|100287102|Transcript|NR_046018.2|transcribed_pseudogene|1/3||NR_046018.2:n.325G>C||325|||||1||1|SNV||EntrezGene|HGNC:37102|YES||||||||||||||||||||

            # Example line 2 (coding SNV):
            # v2 exomes (WES)
            # chr1	69428	rs140739101	T	G	7394499.09	PASS	AC=4068;AC_afr=42;AC_afr_female=30;AC_afr_male=12;AC_amr=185;AC_amr_female=102;AC_amr_male=83;AC_asj=72;AC_asj_female=39;AC_asj_male=33;AC_eas=0;AC_eas_female=0;AC_eas_jpn=0;AC_eas_kor=0;AC_eas_male=0;AC_eas_oea=0;AC_female=1758;AC_fin=566;AC_fin_female=231;AC_fin_male=335;AC_male=2310;AC_nfe=2975;AC_nfe_bgr=62;AC_nfe_est=10;AC_nfe_female=1272;AC_nfe_male=1703;AC_nfe_nwe=1150;AC_nfe_onf=739;AC_nfe_seu=168;AC_nfe_swe=846;AC_oth=114;AC_oth_female=51;AC_oth_male=63;AC_popmax=2975;AC_raw=4512;AC_sas=114;AC_sas_female=33;AC_sas_male=81;AF=2.45436e-02;AF_afr=3.34821e-03;AF_afr_female=3.86399e-03;AF_afr_male=2.51046e-03;AF_amr=9.05886e-03;AF_amr_female=8.51135e-03;AF_amr_male=9.83645e-03;AF_asj=1.39643e-02;AF_asj_female=1.55131e-02;AF_asj_male=1.24905e-02;AF_eas=0.00000e+00;AF_eas_female=0.00000e+00;AF_eas_jpn=0.00000e+00;AF_eas_kor=0.00000e+00;AF_eas_male=0.00000e+00;AF_eas_oea=0.00000e+00;AF_female=2.31566e-02;AF_fin=4.83761e-02;AF_fin_female=4.15020e-02;AF_fin_male=5.46136e-02;AF_male=2.57158e-02;AF_nfe=4.17497e-02;AF_nfe_bgr=4.18919e-02;AF_nfe_est=6.49351e-02;AF_nfe_female=4.04657e-02;AF_nfe_male=4.27632e-02;AF_nfe_nwe=4.32884e-02;AF_nfe_onf=3.84495e-02;AF_nfe_seu=2.41657e-02;AF_nfe_swe=5.01007e-02;AF_oth=2.93512e-02;AF_oth_female=2.75676e-02;AF_oth_male=3.09735e-02;AF_popmax=4.17497e-02;AF_raw=2.15376e-02;AF_sas=4.97165e-03;AF_sas_female=5.74913e-03;AF_sas_male=4.71204e-03;AN=165746;AN_afr=12544;AN_afr_female=7764;AN_afr_male=4780;AN_amr=20422;AN_amr_female=11984;AN_amr_male=8438;AN_asj=5156;AN_asj_female=2514;AN_asj_male=2642;AN_eas=17852;AN_eas_female=9066;AN_eas_jpn=136;AN_eas_kor=3576;AN_eas_male=8786;AN_eas_oea=14140;AN_female=75918;AN_fin=11700;AN_fin_female=5566;AN_fin_male=6134;AN_male=89828;AN_nfe=71258;AN_nfe_bgr=1480;AN_nfe_est=154;AN_nfe_female=31434;AN_nfe_male=39824;AN_nfe_nwe=26566;AN_nfe_onf=19220;AN_nfe_seu=6952;AN_nfe_swe=16886;AN_oth=3884;AN_oth_female=1850;AN_oth_male=2034;AN_popmax=71258;AN_raw=209494;AN_sas=22930;AN_sas_female=5740;AN_sas_male=17190;BaseQRankSum=1.75000e+00;ClippingRankSum=5.50000e-02;DP=5250818;FS=1.12090e+01;InbreedingCoeff=3.21200e-01;MQ=2.55500e+01;MQRankSum=-1.24700e+00;OriginalContig=1;OriginalStart=69428;QD=1.38900e+01;ReadPosRankSum=1.06000e+00;SOR=1.19600e+00;VQSLOD=-6.28500e-01;VQSR_culprit=MQ;ab_hist_alt_bin_freq=0|0|1|15|38|89|149|255|323|133|81|15|8|46|1|9|4|11|0|0;age_hist_het_bin_freq=24|35|60|76|90|78|73|59|54|31;age_hist_het_n_larger=23;age_hist_het_n_smaller=17;age_hist_hom_bin_freq=47|58|101|99|124|130|108|109|79|65;age_hist_hom_n_larger=28;age_hist_hom_n_smaller=31;allele_type=snv;controls_AC=1686;controls_AC_afr=19;controls_AC_afr_female=15;controls_AC_afr_male=4;controls_AC_amr=72;controls_AC_amr_female=36;controls_AC_amr_male=36;controls_AC_asj=16;controls_AC_asj_female=5;controls_AC_asj_male=11;controls_AC_eas=0;controls_AC_eas_female=0;controls_AC_eas_jpn=0;controls_AC_eas_kor=0;controls_AC_eas_male=0;controls_AC_eas_oea=0;controls_AC_female=743;controls_AC_fin=337;controls_AC_fin_female=144;controls_AC_fin_male=193;controls_AC_male=943;controls_AC_nfe=1139;controls_AC_nfe_bgr=20;controls_AC_nfe_est=5;controls_AC_nfe_female=501;controls_AC_nfe_male=638;controls_AC_nfe_nwe=395;controls_AC_nfe_onf=236;controls_AC_nfe_seu=67;controls_AC_nfe_swe=416;controls_AC_oth=42;controls_AC_oth_female=19;controls_AC_oth_male=23;controls_AC_popmax=1139;controls_AC_raw=1856;controls_AC_sas=61;controls_AC_sas_female=23;controls_AC_sas_male=38;controls_AF=2.31530e-02;controls_AF_afr=3.40624e-03;controls_AF_afr_female=4.53995e-03;controls_AF_afr_male=1.75901e-03;controls_AF_amr=7.19568e-03;controls_AF_amr_female=5.97213e-03;controls_AF_amr_male=9.04977e-03;controls_AF_asj=1.40351e-02;controls_AF_asj_female=8.25083e-03;controls_AF_asj_male=2.05993e-02;controls_AF_eas=0.00000e+00;controls_AF_eas_female=0.00000e+00;controls_AF_eas_jpn=0.00000e+00;controls_AF_eas_kor=0.00000e+00;controls_AF_eas_male=0.00000e+00;controls_AF_eas_oea=0.00000e+00;controls_AF_female=2.19589e-02;controls_AF_fin=4.65598e-02;controls_AF_fin_female=4.20807e-02;controls_AF_fin_male=5.05765e-02;controls_AF_male=2.41894e-02;controls_AF_nfe=4.20481e-02;controls_AF_nfe_bgr=5.10204e-02;controls_AF_nfe_est=1.31579e-01;controls_AF_nfe_female=4.19177e-02;controls_AF_nfe_male=4.21512e-02;controls_AF_nfe_nwe=4.29815e-02;controls_AF_nfe_onf=3.75199e-02;controls_AF_nfe_seu=2.33775e-02;controls_AF_nfe_swe=5.00481e-02;controls_AF_oth=3.39806e-02;controls_AF_oth_female=2.86145e-02;controls_AF_oth_male=4.02098e-02;controls_AF_popmax=4.20481e-02;controls_AF_raw=2.02757e-02;controls_AF_sas=5.18972e-03;controls_AF_sas_female=7.15619e-03;controls_AF_sas_male=4.44965e-03;controls_AN=72820;controls_AN_afr=5578;controls_AN_afr_female=3304;controls_AN_afr_male=2274;controls_AN_amr=10006;controls_AN_amr_female=6028;controls_AN_amr_male=3978;controls_AN_asj=1140;controls_AN_asj_female=606;controls_AN_asj_male=534;controls_AN_eas=8780;controls_AN_eas_female=4646;controls_AN_eas_jpn=102;controls_AN_eas_kor=1764;controls_AN_eas_male=4134;controls_AN_eas_oea=6914;controls_AN_female=33836;controls_AN_fin=7238;controls_AN_fin_female=3422;controls_AN_fin_male=3816;controls_AN_male=38984;controls_AN_nfe=27088;controls_AN_nfe_bgr=392;controls_AN_nfe_est=38;controls_AN_nfe_female=11952;controls_AN_nfe_male=15136;controls_AN_nfe_nwe=9190;controls_AN_nfe_onf=6290;controls_AN_nfe_seu=2866;controls_AN_nfe_swe=8312;controls_AN_oth=1236;controls_AN_oth_female=664;controls_AN_oth_male=572;controls_AN_popmax=27088;controls_AN_raw=91538;controls_AN_sas=11754;controls_AN_sas_female=3214;controls_AN_sas_male=8540;controls_faf95=2.22333e-02;controls_faf95_afr=2.23026e-03;controls_faf95_amr=5.85936e-03;controls_faf95_eas=0.00000e+00;controls_faf95_nfe=4.00191e-02;controls_faf95_sas=4.14649e-03;controls_faf99=2.22333e-02;controls_faf99_afr=2.23000e-03;controls_faf99_amr=5.85932e-03;controls_faf99_eas=0.00000e+00;controls_faf99_nfe=4.00198e-02;controls_faf99_sas=4.14654e-03;controls_nhomalt=650;controls_nhomalt_afr=7;controls_nhomalt_afr_female=6;controls_nhomalt_afr_male=1;controls_nhomalt_amr=28;controls_nhomalt_amr_female=14;controls_nhomalt_amr_male=14;controls_nhomalt_asj=7;controls_nhomalt_asj_female=2;controls_nhomalt_asj_male=5;controls_nhomalt_eas=0;controls_nhomalt_eas_female=0;controls_nhomalt_eas_jpn=0;controls_nhomalt_eas_kor=0;controls_nhomalt_eas_male=0;controls_nhomalt_eas_oea=0;controls_nhomalt_female=282;controls_nhomalt_fin=142;controls_nhomalt_fin_female=61;controls_nhomalt_fin_male=81;controls_nhomalt_male=368;controls_nhomalt_nfe=431;controls_nhomalt_nfe_bgr=6;controls_nhomalt_nfe_est=2;controls_nhomalt_nfe_female=185;controls_nhomalt_nfe_male=246;controls_nhomalt_nfe_nwe=150;controls_nhomalt_nfe_onf=89;controls_nhomalt_nfe_seu=26;controls_nhomalt_nfe_swe=158;controls_nhomalt_oth=15;controls_nhomalt_oth_female=7;controls_nhomalt_oth_male=8;controls_nhomalt_popmax=431;controls_nhomalt_raw=688;controls_nhomalt_sas=20;controls_nhomalt_sas_female=7;controls_nhomalt_sas_male=13;controls_popmax=nfe;dp_hist_all_bin_freq=40005|2733|572|276|434|371|25490|39940|9514|2882|826|357|230|156|120|89|98|92|75|101;dp_hist_all_n_larger=1387;dp_hist_alt_bin_freq=176|126|36|23|17|30|52|66|83|118|131|115|76|76|54|52|66|69|60|85;dp_hist_alt_n_larger=1334;faf95=2.39136e-02;faf95_afr=2.54576e-03;faf95_amr=7.99100e-03;faf95_eas=0.00000e+00;faf95_nfe=4.04978e-02;faf95_sas=4.23127e-03;faf99=2.39136e-02;faf99_afr=2.54575e-03;faf99_amr=7.99128e-03;faf99_eas=0.00000e+00;faf99_nfe=4.04986e-02;faf99_sas=4.23070e-03;gq_hist_all_bin_freq=11880|6528|1156|1299|692|233|284|188|79|96|113|79|258|45|179|63|274|60|270|80971;gq_hist_alt_bin_freq=17|83|24|52|20|27|50|33|12|8|21|16|12|11|19|13|26|20|23|2358;n_alt_alleles=1;nhomalt=1552;nhomalt_afr=14;nhomalt_afr_female=10;nhomalt_afr_male=4;nhomalt_amr=73;nhomalt_amr_female=42;nhomalt_amr_male=31;nhomalt_asj=30;nhomalt_asj_female=17;nhomalt_asj_male=13;nhomalt_eas=0;nhomalt_eas_female=0;nhomalt_eas_jpn=0;nhomalt_eas_kor=0;nhomalt_eas_male=0;nhomalt_eas_oea=0;nhomalt_female=666;nhomalt_fin=237;nhomalt_fin_female=97;nhomalt_fin_male=140;nhomalt_male=886;nhomalt_nfe=1120;nhomalt_nfe_bgr=24;nhomalt_nfe_est=4;nhomalt_nfe_female=471;nhomalt_nfe_male=649;nhomalt_nfe_nwe=438;nhomalt_nfe_onf=275;nhomalt_nfe_seu=63;nhomalt_nfe_swe=316;nhomalt_oth=43;nhomalt_oth_female=19;nhomalt_oth_male=24;nhomalt_popmax=1120;nhomalt_raw=1667;nhomalt_sas=35;nhomalt_sas_female=10;nhomalt_sas_male=25;non_cancer_AC=3771;non_cancer_AC_afr=39;non_cancer_AC_afr_female=29;non_cancer_AC_afr_male=10;non_cancer_AC_amr=177;non_cancer_AC_amr_female=102;non_cancer_AC_amr_male=75;non_cancer_AC_asj=69;non_cancer_AC_asj_female=36;non_cancer_AC_asj_male=33;non_cancer_AC_eas=0;non_cancer_AC_eas_female=0;non_cancer_AC_eas_jpn=0;non_cancer_AC_eas_kor=0;non_cancer_AC_eas_male=0;non_cancer_AC_eas_oea=0;non_cancer_AC_female=1639;non_cancer_AC_fin=564;non_cancer_AC_fin_female=229;non_cancer_AC_fin_male=335;non_cancer_AC_male=2132;non_cancer_AC_nfe=2703;non_cancer_AC_nfe_bgr=61;non_cancer_AC_nfe_est=10;non_cancer_AC_nfe_female=1164;non_cancer_AC_nfe_male=1539;non_cancer_AC_nfe_nwe=1081;non_cancer_AC_nfe_onf=572;non_cancer_AC_nfe_seu=153;non_cancer_AC_nfe_swe=826;non_cancer_AC_oth=107;non_cancer_AC_oth_female=48;non_cancer_AC_oth_male=59;non_cancer_AC_popmax=2703;non_cancer_AC_raw=4154;non_cancer_AC_sas=112;non_cancer_AC_sas_female=31;non_cancer_AC_sas_male=81;non_cancer_AF=2.41316e-02;non_cancer_AF_afr=3.39307e-03;non_cancer_AF_afr_female=4.08106e-03;non_cancer_AF_afr_male=2.27894e-03;non_cancer_AF_amr=8.75025e-03;non_cancer_AF_amr_female=8.58441e-03;non_cancer_AF_amr_male=8.98634e-03;non_cancer_AF_asj=1.41509e-02;non_cancer_AF_asj_female=1.54242e-02;non_cancer_AF_asj_male=1.29819e-02;non_cancer_AF_eas=0.00000e+00;non_cancer_AF_eas_female=0.00000e+00;non_cancer_AF_eas_jpn=0.00000e+00;non_cancer_AF_eas_kor=0.00000e+00;non_cancer_AF_eas_male=0.00000e+00;non_cancer_AF_eas_oea=0.00000e+00;non_cancer_AF_female=2.31608e-02;non_cancer_AF_fin=4.82381e-02;non_cancer_AF_fin_female=4.12019e-02;non_cancer_AF_fin_male=5.46136e-02;non_cancer_AF_male=2.49351e-02;non_cancer_AF_nfe=4.19786e-02;non_cancer_AF_nfe_bgr=4.35093e-02;non_cancer_AF_nfe_est=1.11111e-01;non_cancer_AF_nfe_female=4.18886e-02;non_cancer_AF_nfe_male=4.20469e-02;non_cancer_AF_nfe_nwe=4.34451e-02;non_cancer_AF_nfe_onf=3.74591e-02;non_cancer_AF_nfe_seu=2.38466e-02;non_cancer_AF_nfe_swe=5.05818e-02;non_cancer_AF_oth=3.01239e-02;non_cancer_AF_oth_female=2.83353e-02;non_cancer_AF_oth_male=3.17546e-02;non_cancer_AF_popmax=4.19786e-02;non_cancer_AF_raw=2.11079e-02;non_cancer_AF_sas=4.89853e-03;non_cancer_AF_sas_female=5.43288e-03;non_cancer_AF_sas_male=4.72083e-03;non_cancer_AN=156268;non_cancer_AN_afr=11494;non_cancer_AN_afr_female=7106;non_cancer_AN_afr_male=4388;non_cancer_AN_amr=20228;non_cancer_AN_amr_female=11882;non_cancer_AN_amr_male=8346;non_cancer_AN_asj=4876;non_cancer_AN_asj_female=2334;non_cancer_AN_asj_male=2542;non_cancer_AN_eas=17172;non_cancer_AN_eas_female=8698;non_cancer_AN_eas_jpn=112;non_cancer_AN_eas_kor=3536;non_cancer_AN_eas_male=8474;non_cancer_AN_eas_oea=13524;non_cancer_AN_female=70766;non_cancer_AN_fin=11692;non_cancer_AN_fin_female=5558;non_cancer_AN_fin_male=6134;non_cancer_AN_male=85502;non_cancer_AN_nfe=64390;non_cancer_AN_nfe_bgr=1402;non_cancer_AN_nfe_est=90;non_cancer_AN_nfe_female=27788;non_cancer_AN_nfe_male=36602;non_cancer_AN_nfe_nwe=24882;non_cancer_AN_nfe_onf=15270;non_cancer_AN_nfe_seu=6416;non_cancer_AN_nfe_swe=16330;non_cancer_AN_oth=3552;non_cancer_AN_oth_female=1694;non_cancer_AN_oth_male=1858;non_cancer_AN_popmax=64390;non_cancer_AN_raw=196798;non_cancer_AN_sas=22864;non_cancer_AN_sas_female=5706;non_cancer_AN_sas_male=17158;non_cancer_faf95=2.34883e-02;non_cancer_faf95_afr=2.55139e-03;non_cancer_faf95_amr=7.69649e-03;non_cancer_faf95_eas=0.00000e+00;non_cancer_faf95_nfe=4.06589e-02;non_cancer_faf95_sas=4.16222e-03;non_cancer_faf99=2.34888e-02;non_cancer_faf99_afr=2.55064e-03;non_cancer_faf99_amr=7.69683e-03;non_cancer_faf99_eas=0.00000e+00;non_cancer_faf99_nfe=4.06585e-02;non_cancer_faf99_sas=4.16247e-03;non_cancer_nhomalt=1441;non_cancer_nhomalt_afr=14;non_cancer_nhomalt_afr_female=10;non_cancer_nhomalt_afr_male=4;non_cancer_nhomalt_amr=69;non_cancer_nhomalt_amr_female=42;non_cancer_nhomalt_amr_male=27;non_cancer_nhomalt_asj=29;non_cancer_nhomalt_asj_female=16;non_cancer_nhomalt_asj_male=13;non_cancer_nhomalt_eas=0;non_cancer_nhomalt_eas_female=0;non_cancer_nhomalt_eas_jpn=0;non_cancer_nhomalt_eas_kor=0;non_cancer_nhomalt_eas_male=0;non_cancer_nhomalt_eas_oea=0;non_cancer_nhomalt_female=627;non_cancer_nhomalt_fin=236;non_cancer_nhomalt_fin_female=96;non_cancer_nhomalt_fin_male=140;non_cancer_nhomalt_male=814;non_cancer_nhomalt_nfe=1019;non_cancer_nhomalt_nfe_bgr=24;non_cancer_nhomalt_nfe_est=4;non_cancer_nhomalt_nfe_female=436;non_cancer_nhomalt_nfe_male=583;non_cancer_nhomalt_nfe_nwe=414;non_cancer_nhomalt_nfe_onf=212;non_cancer_nhomalt_nfe_seu=57;non_cancer_nhomalt_nfe_swe=308;non_cancer_nhomalt_oth=40;non_cancer_nhomalt_oth_female=18;non_cancer_nhomalt_oth_male=22;non_cancer_nhomalt_popmax=1019;non_cancer_nhomalt_raw=1544;non_cancer_nhomalt_sas=34;non_cancer_nhomalt_sas_female=9;non_cancer_nhomalt_sas_male=25;non_cancer_popmax=nfe;non_neuro_AC=3129;non_neuro_AC_afr=42;non_neuro_AC_afr_female=30;non_neuro_AC_afr_male=12;non_neuro_AC_amr=144;non_neuro_AC_amr_female=76;non_neuro_AC_amr_male=68;non_neuro_AC_asj=41;non_neuro_AC_asj_female=24;non_neuro_AC_asj_male=17;non_neuro_AC_eas=0;non_neuro_AC_eas_female=0;non_neuro_AC_eas_jpn=0;non_neuro_AC_eas_kor=0;non_neuro_AC_eas_male=0;non_neuro_AC_eas_oea=0;non_neuro_AC_female=1351;non_neuro_AC_fin=412;non_neuro_AC_fin_female=170;non_neuro_AC_fin_male=242;non_neuro_AC_male=1778;non_neuro_AC_nfe=2295;non_neuro_AC_nfe_bgr=8;non_neuro_AC_nfe_est=10;non_neuro_AC_nfe_female=980;non_neuro_AC_nfe_male=1315;non_neuro_AC_nfe_nwe=947;non_neuro_AC_nfe_onf=691;non_neuro_AC_nfe_seu=161;non_neuro_AC_nfe_swe=478;non_neuro_AC_oth=81;non_neuro_AC_oth_female=38;non_neuro_AC_oth_male=43;non_neuro_AC_popmax=2295;non_neuro_AC_raw=3456;non_neuro_AC_sas=114;non_neuro_AC_sas_female=33;non_neuro_AC_sas_male=81;non_neuro_AF=2.26559e-02;non_neuro_AF_afr=3.35678e-03;non_neuro_AF_afr_female=3.86997e-03;non_neuro_AF_afr_male=2.52101e-03;non_neuro_AF_amr=8.05099e-03;non_neuro_AF_amr_female=7.07767e-03;non_neuro_AF_amr_male=9.51315e-03;non_neuro_AF_asj=1.29501e-02;non_neuro_AF_asj_female=1.50754e-02;non_neuro_AF_asj_male=1.08005e-02;non_neuro_AF_eas=0.00000e+00;non_neuro_AF_eas_female=0.00000e+00;non_neuro_AF_eas_jpn=0.00000e+00;non_neuro_AF_eas_kor=0.00000e+00;non_neuro_AF_eas_male=0.00000e+00;non_neuro_AF_eas_oea=0.00000e+00;non_neuro_AF_female=2.13922e-02;non_neuro_AF_fin=4.53744e-02;non_neuro_AF_fin_female=4.15242e-02;non_neuro_AF_fin_male=4.85359e-02;non_neuro_AF_male=2.37206e-02;non_neuro_AF_nfe=4.06612e-02;non_neuro_AF_nfe_bgr=3.44828e-02;non_neuro_AF_nfe_est=6.94444e-02;non_neuro_AF_nfe_female=3.90781e-02;non_neuro_AF_nfe_male=4.19271e-02;non_neuro_AF_nfe_nwe=4.23032e-02;non_neuro_AF_nfe_onf=3.96853e-02;non_neuro_AF_nfe_seu=2.41960e-02;non_neuro_AF_nfe_swe=4.97192e-02;non_neuro_AF_oth=2.59449e-02;non_neuro_AF_oth_female=2.48042e-02;non_neuro_AF_oth_male=2.70440e-02;non_neuro_AF_popmax=4.06612e-02;non_neuro_AF_raw=1.99478e-02;non_neuro_AF_sas=4.97252e-03;non_neuro_AF_sas_female=5.75113e-03;non_neuro_AF_sas_male=4.71259e-03;non_neuro_AN=138110;non_neuro_AN_afr=12512;non_neuro_AN_afr_female=7752;non_neuro_AN_afr_male=4760;non_neuro_AN_amr=17886;non_neuro_AN_amr_female=10738;non_neuro_AN_amr_male=7148;non_neuro_AN_asj=3166;non_neuro_AN_asj_female=1592;non_neuro_AN_asj_male=1574;non_neuro_AN_eas=12976;non_neuro_AN_eas_female=6630;non_neuro_AN_eas_jpn=134;non_neuro_AN_eas_kor=3574;non_neuro_AN_eas_male=6346;non_neuro_AN_eas_oea=9268;non_neuro_AN_female=63154;non_neuro_AN_fin=9080;non_neuro_AN_fin_female=4094;non_neuro_AN_fin_male=4986;non_neuro_AN_male=74956;non_neuro_AN_nfe=56442;non_neuro_AN_nfe_bgr=232;non_neuro_AN_nfe_est=144;non_neuro_AN_nfe_female=25078;non_neuro_AN_nfe_male=31364;non_neuro_AN_nfe_nwe=22386;non_neuro_AN_nfe_onf=17412;non_neuro_AN_nfe_seu=6654;non_neuro_AN_nfe_swe=9614;non_neuro_AN_oth=3122;non_neuro_AN_oth_female=1532;non_neuro_AN_oth_male=1590;non_neuro_AN_popmax=56442;non_neuro_AN_raw=173252;non_neuro_AN_sas=22926;non_neuro_AN_sas_female=5738;non_neuro_AN_sas_male=17188;non_neuro_faf95=2.19937e-02;non_neuro_faf95_afr=2.55227e-03;non_neuro_faf95_amr=6.97927e-03;non_neuro_faf95_eas=0.00000e+00;non_neuro_faf95_nfe=3.92753e-02;non_neuro_faf95_sas=4.23109e-03;non_neuro_faf99=2.19936e-02;non_neuro_faf99_afr=2.55238e-03;non_neuro_faf99_amr=6.97941e-03;non_neuro_faf99_eas=0.00000e+00;non_neuro_faf99_nfe=3.92748e-02;non_neuro_faf99_sas=4.23138e-03;non_neuro_nhomalt=1192;non_neuro_nhomalt_afr=14;non_neuro_nhomalt_afr_female=10;non_neuro_nhomalt_afr_male=4;non_neuro_nhomalt_amr=59;non_neuro_nhomalt_amr_female=31;non_neuro_nhomalt_amr_male=28;non_neuro_nhomalt_asj=17;non_neuro_nhomalt_asj_female=10;non_neuro_nhomalt_asj_male=7;non_neuro_nhomalt_eas=0;non_neuro_nhomalt_eas_female=0;non_neuro_nhomalt_eas_jpn=0;non_neuro_nhomalt_eas_kor=0;non_neuro_nhomalt_eas_male=0;non_neuro_nhomalt_eas_oea=0;non_neuro_nhomalt_female=509;non_neuro_nhomalt_fin=171;non_neuro_nhomalt_fin_female=72;non_neuro_nhomalt_fin_male=99;non_neuro_nhomalt_male=683;non_neuro_nhomalt_nfe=866;non_neuro_nhomalt_nfe_bgr=3;non_neuro_nhomalt_nfe_est=4;non_neuro_nhomalt_nfe_female=362;non_neuro_nhomalt_nfe_male=504;non_neuro_nhomalt_nfe_nwe=358;non_neuro_nhomalt_nfe_onf=260;non_neuro_nhomalt_nfe_seu=60;non_neuro_nhomalt_nfe_swe=181;non_neuro_nhomalt_oth=30;non_neuro_nhomalt_oth_female=14;non_neuro_nhomalt_oth_male=16;non_neuro_nhomalt_popmax=866;non_neuro_nhomalt_raw=1269;non_neuro_nhomalt_sas=35;non_neuro_nhomalt_sas_female=10;non_neuro_nhomalt_sas_male=25;non_neuro_popmax=nfe;non_topmed_AC=3990;non_topmed_AC_afr=30;non_topmed_AC_afr_female=22;non_topmed_AC_afr_male=8;non_topmed_AC_amr=183;non_topmed_AC_amr_female=100;non_topmed_AC_amr_male=83;non_topmed_AC_asj=68;non_topmed_AC_asj_female=37;non_topmed_AC_asj_male=31;non_topmed_AC_eas=0;non_topmed_AC_eas_female=0;non_topmed_AC_eas_jpn=0;non_topmed_AC_eas_kor=0;non_topmed_AC_eas_male=0;non_topmed_AC_eas_oea=0;non_topmed_AC_female=1713;non_topmed_AC_fin=566;non_topmed_AC_fin_female=231;non_topmed_AC_fin_male=335;non_topmed_AC_male=2277;non_topmed_AC_nfe=2921;non_topmed_AC_nfe_bgr=62;non_topmed_AC_nfe_est=10;non_topmed_AC_nfe_female=1245;non_topmed_AC_nfe_male=1676;non_topmed_AC_nfe_nwe=1117;non_topmed_AC_nfe_onf=723;non_topmed_AC_nfe_seu=165;non_topmed_AC_nfe_swe=844;non_topmed_AC_oth=108;non_topmed_AC_oth_female=45;non_topmed_AC_oth_male=63;non_topmed_AC_popmax=2921;non_topmed_AC_raw=4427;non_topmed_AC_sas=114;non_topmed_AC_sas_female=33;non_topmed_AC_sas_male=81;non_topmed_AF=2.47882e-02;non_topmed_AF_afr=3.25027e-03;non_topmed_AF_afr_female=3.90071e-03;non_topmed_AF_afr_male=2.22841e-03;non_topmed_AF_amr=9.00502e-03;non_topmed_AF_amr_female=8.38364e-03;non_topmed_AF_amr_male=9.88802e-03;non_topmed_AF_asj=1.32968e-02;non_topmed_AF_asj_female=1.48594e-02;non_topmed_AF_asj_male=1.18140e-02;non_topmed_AF_eas=0.00000e+00;non_topmed_AF_eas_female=0.00000e+00;non_topmed_AF_eas_jpn=0.00000e+00;non_topmed_AF_eas_kor=0.00000e+00;non_topmed_AF_eas_male=0.00000e+00;non_topmed_AF_eas_oea=0.00000e+00;non_topmed_AF_female=2.35018e-02;non_topmed_AF_fin=4.83843e-02;non_topmed_AF_fin_female=4.15169e-02;non_topmed_AF_fin_male=5.46136e-02;non_topmed_AF_male=2.58527e-02;non_topmed_AF_nfe=4.17381e-02;non_topmed_AF_nfe_bgr=4.20054e-02;non_topmed_AF_nfe_est=6.66667e-02;non_topmed_AF_nfe_female=4.06305e-02;non_topmed_AF_nfe_male=4.26008e-02;non_topmed_AF_nfe_nwe=4.31808e-02;non_topmed_AF_nfe_onf=3.85600e-02;non_topmed_AF_nfe_seu=2.39339e-02;non_topmed_AF_nfe_swe=5.01009e-02;non_topmed_AF_oth=2.81397e-02;non_topmed_AF_oth_female=2.47525e-02;non_topmed_AF_oth_male=3.11881e-02;non_topmed_AF_popmax=4.17381e-02;non_topmed_AF_raw=2.17340e-02;non_topmed_AF_sas=4.97165e-03;non_topmed_AF_sas_female=5.74913e-03;non_topmed_AF_sas_male=4.71204e-03;non_topmed_AN=160964;non_topmed_AN_afr=9230;non_topmed_AN_afr_female=5640;non_topmed_AN_afr_male=3590;non_topmed_AN_amr=20322;non_topmed_AN_amr_female=11928;non_topmed_AN_amr_male=8394;non_topmed_AN_asj=5114;non_topmed_AN_asj_female=2490;non_topmed_AN_asj_male=2624;non_topmed_AN_eas=17848;non_topmed_AN_eas_female=9066;non_topmed_AN_eas_jpn=136;non_topmed_AN_eas_kor=3576;non_topmed_AN_eas_male=8782;non_topmed_AN_eas_oea=14136;non_topmed_AN_female=72888;non_topmed_AN_fin=11698;non_topmed_AN_fin_female=5564;non_topmed_AN_fin_male=6134;non_topmed_AN_male=88076;non_topmed_AN_nfe=69984;non_topmed_AN_nfe_bgr=1476;non_topmed_AN_nfe_est=150;non_topmed_AN_nfe_female=30642;non_topmed_AN_nfe_male=39342;non_topmed_AN_nfe_nwe=25868;non_topmed_AN_nfe_onf=18750;non_topmed_AN_nfe_seu=6894;non_topmed_AN_nfe_swe=16846;non_topmed_AN_oth=3838;non_topmed_AN_oth_female=1818;non_topmed_AN_oth_male=2020;non_topmed_AN_popmax=69984;non_topmed_AN_raw=203690;non_topmed_AN_sas=22930;non_topmed_AN_sas_female=5740;non_topmed_AN_sas_male=17190;non_topmed_faf95=2.41458e-02;non_topmed_faf95_afr=2.33920e-03;non_topmed_faf95_amr=7.93860e-03;non_topmed_faf95_eas=0.00000e+00;non_topmed_faf95_nfe=4.04759e-02;non_topmed_faf95_sas=4.23127e-03;non_topmed_faf99=2.41458e-02;non_topmed_faf99_afr=2.33918e-03;non_topmed_faf99_amr=7.93807e-03;non_topmed_faf99_eas=0.00000e+00;non_topmed_faf99_nfe=4.04757e-02;non_topmed_faf99_sas=4.23070e-03;non_topmed_nhomalt=1522;non_topmed_nhomalt_afr=10;non_topmed_nhomalt_afr_female=7;non_topmed_nhomalt_afr_male=3;non_topmed_nhomalt_amr=72;non_topmed_nhomalt_amr_female=41;non_topmed_nhomalt_amr_male=31;non_topmed_nhomalt_asj=28;non_topmed_nhomalt_asj_female=16;non_topmed_nhomalt_asj_male=12;non_topmed_nhomalt_eas=0;non_topmed_nhomalt_eas_female=0;non_topmed_nhomalt_eas_jpn=0;non_topmed_nhomalt_eas_kor=0;non_topmed_nhomalt_eas_male=0;non_topmed_nhomalt_eas_oea=0;non_topmed_nhomalt_female=650;non_topmed_nhomalt_fin=237;non_topmed_nhomalt_fin_female=97;non_topmed_nhomalt_fin_male=140;non_topmed_nhomalt_male=872;non_topmed_nhomalt_nfe=1099;non_topmed_nhomalt_nfe_bgr=24;non_topmed_nhomalt_nfe_est=4;non_topmed_nhomalt_nfe_female=462;non_topmed_nhomalt_nfe_male=637;non_topmed_nhomalt_nfe_nwe=425;non_topmed_nhomalt_nfe_onf=269;non_topmed_nhomalt_nfe_seu=62;non_topmed_nhomalt_nfe_swe=315;non_topmed_nhomalt_oth=41;non_topmed_nhomalt_oth_female=17;non_topmed_nhomalt_oth_male=24;non_topmed_nhomalt_popmax=1099;non_topmed_nhomalt_raw=1636;non_topmed_nhomalt_sas=35;non_topmed_nhomalt_sas_female=10;non_topmed_nhomalt_sas_male=25;non_topmed_popmax=nfe;pab_max=1.00000e+00;popmax=nfe;rf_label=FP;rf_negative_label;rf_tp_probability=8.47237e-01;rf_train;segdup;variant_type=snv;vep=G|missense_variant|MODERATE|OR4F5|ENSG00000186092|Transcript|ENST00000335137|protein_coding|1/1||ENST00000335137.3:c.338T>G|ENSP00000334393.3:p.Phe113Cys|338|338|113|F/C|tTt/tGt|rs140739101|1||1||SNV|1|HGNC|14825|YES|||CCDS30547.1|ENSP00000334393|Q8NH21||UPI0000041BC1||deleterious(0.01)|possibly_damaging(0.568)|Transmembrane_helices:TMhelix&Prints_domain:PR00237&Superfamily_domains:SSF81321&Gene3D:1.20.1070.10&PROSITE_patterns:PS00237&hmmpanther:PTHR26451&hmmpanther:PTHR26451:SF72&PROSITE_profiles:PS50262||G:0.0190|G:0.004888|G:0.0015|G:0.036|G:0.003|G:0.0497|G:0.0153|G:0.0037|G:0.0457|G:0.008045|G:0.022|G:0.002553|G:0.02462|G:0|G:0.04624|G:0.04058|G:0.02716||||||||||||,G|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00000278218|open_chromatin_region||||||||||rs140739101|1||||SNV|1||||||||||||||||G:0.0190|G:0.004888|G:0.0015|G:0.036|G:0.003|G:0.0497|G:0.0153|G:0.0037|G:0.0457|G:0.008045|G:0.022|G:0.002553|G:0.02462|G:0|G:0.04624|G:0.04058|G:0.02716||||||||||||;CSQ=G|ENSG00000186092|Transcript|ENST00000641515|missense_variant||134|F/C|
            # v2 genomes (WGS)
            # chr1	69428	rs140739101	T	G	40085.55	PASS	AC=103;AC_afr=9;AC_afr_female=3;AC_afr_male=6;AC_amr=0;AC_amr_female=0;AC_amr_male=0;AC_asj=1;AC_asj_female=0;AC_asj_male=1;AC_eas=0;AC_eas_female=0;AC_eas_male=0;AC_female=41;AC_fin=10;AC_fin_female=8;AC_fin_male=2;AC_male=62;AC_nfe=83;AC_nfe_est=30;AC_nfe_female=30;AC_nfe_male=53;AC_nfe_nwe=42;AC_nfe_onf=11;AC_nfe_seu=0;AC_oth=0;AC_oth_female=0;AC_oth_male=0;AC_popmax=83;AC_raw=457;AF=7.35504e-03;AF_afr=1.82556e-03;AF_afr_female=1.43266e-03;AF_afr_male=2.11566e-03;AF_amr=0.00000e+00;AF_amr_female=0.00000e+00;AF_amr_male=0.00000e+00;AF_asj=1.06383e-02;AF_asj_female=0.00000e+00;AF_asj_male=1.47059e-02;AF_eas=0.00000e+00;AF_eas_female=0.00000e+00;AF_eas_male=0.00000e+00;AF_female=6.92100e-03;AF_fin=8.46024e-03;AF_fin_female=1.25392e-02;AF_fin_male=3.67647e-03;AF_male=7.67327e-03;AF_nfe=1.43997e-02;AF_nfe_est=2.23881e-02;AF_nfe_female=1.25733e-02;AF_nfe_male=1.56898e-02;AF_nfe_nwe=1.19658e-02;AF_nfe_onf=1.25571e-02;AF_nfe_seu=0.00000e+00;AF_oth=0.00000e+00;AF_oth_female=0.00000e+00;AF_oth_male=0.00000e+00;AF_popmax=1.43997e-02;AF_raw=2.06694e-02;AN=14004;AN_afr=4930;AN_afr_female=2094;AN_afr_male=2836;AN_amr=328;AN_amr_female=142;AN_amr_male=186;AN_asj=94;AN_asj_female=26;AN_asj_male=68;AN_eas=1310;AN_eas_female=464;AN_eas_male=846;AN_female=5924;AN_fin=1182;AN_fin_female=638;AN_fin_male=544;AN_male=8080;AN_nfe=5764;AN_nfe_est=1340;AN_nfe_female=2386;AN_nfe_male=3378;AN_nfe_nwe=3510;AN_nfe_onf=876;AN_nfe_seu=38;AN_oth=396;AN_oth_female=174;AN_oth_male=222;AN_popmax=5764;AN_raw=22110;BaseQRankSum=-5.55000e-01;ClippingRankSum=7.40000e-02;DP=189138;FS=5.24000e-01;InbreedingCoeff=3.77700e-01;MQ=2.50000e+01;MQRankSum=-7.65000e-01;OriginalContig=1;OriginalStart=69428;QD=1.02500e+01;ReadPosRankSum=6.04000e-01;SOR=6.40000e-01;VQSLOD=-5.41900e+01;VQSR_culprit=MQ;ab_hist_alt_bin_freq=0|0|1|11|13|22|22|12|19|3|13|4|1|4|2|5|7|2|0|0;age_hist_het_bin_freq=9|4|5|6|11|3|10|4|2|0;age_hist_het_n_larger=0;age_hist_het_n_smaller=9;age_hist_hom_bin_freq=0|1|0|1|1|1|0|1|0|0;age_hist_hom_n_larger=0;age_hist_hom_n_smaller=2;allele_type=snv;controls_AC=36;controls_AC_afr=1;controls_AC_afr_female=0;controls_AC_afr_male=1;controls_AC_amr=0;controls_AC_amr_female=0;controls_AC_amr_male=0;controls_AC_asj=0;controls_AC_asj_female=0;controls_AC_asj_male=0;controls_AC_eas=0;controls_AC_eas_female=0;controls_AC_eas_male=0;controls_AC_female=18;controls_AC_fin=4;controls_AC_fin_female=4;controls_AC_fin_male=0;controls_AC_male=18;controls_AC_nfe=31;controls_AC_nfe_est=24;controls_AC_nfe_female=14;controls_AC_nfe_male=17;controls_AC_nfe_nwe=5;controls_AC_nfe_onf=2;controls_AC_nfe_seu=0;controls_AC_oth=0;controls_AC_oth_female=0;controls_AC_oth_male=0;controls_AC_popmax=31;controls_AC_raw=186;controls_AF=7.81929e-03;controls_AF_afr=6.98324e-04;controls_AF_afr_female=0.00000e+00;controls_AF_afr_male=1.25628e-03;controls_AF_amr=0.00000e+00;controls_AF_amr_female=0.00000e+00;controls_AF_amr_male=0.00000e+00;controls_AF_asj=0.00000e+00;controls_AF_asj_female=0.00000e+00;controls_AF_asj_male=0.00000e+00;controls_AF_eas=0.00000e+00;controls_AF_eas_female=0.00000e+00;controls_AF_eas_male=0.00000e+00;controls_AF_female=8.75486e-03;controls_AF_fin=1.00000e-02;controls_AF_fin_female=1.94175e-02;controls_AF_fin_male=0.00000e+00;controls_AF_male=7.06436e-03;controls_AF_nfe=1.79191e-02;controls_AF_nfe_est=2.03735e-02;controls_AF_nfe_female=1.76768e-02;controls_AF_nfe_male=1.81237e-02;controls_AF_nfe_nwe=1.45349e-02;controls_AF_nfe_onf=1.07527e-02;controls_AF_nfe_seu=0.00000e+00;controls_AF_oth=0.00000e+00;controls_AF_oth_female=0.00000e+00;controls_AF_oth_male=0.00000e+00;controls_AF_popmax=1.79191e-02;controls_AF_raw=2.49731e-02;controls_AN=4604;controls_AN_afr=1432;controls_AN_afr_female=636;controls_AN_afr_male=796;controls_AN_amr=108;controls_AN_amr_female=54;controls_AN_amr_male=54;controls_AN_asj=20;controls_AN_asj_female=12;controls_AN_asj_male=8;controls_AN_eas=782;controls_AN_eas_female=302;controls_AN_eas_male=480;controls_AN_female=2056;controls_AN_fin=400;controls_AN_fin_female=206;controls_AN_fin_male=194;controls_AN_male=2548;controls_AN_nfe=1730;controls_AN_nfe_est=1178;controls_AN_nfe_female=792;controls_AN_nfe_male=938;controls_AN_nfe_nwe=344;controls_AN_nfe_onf=186;controls_AN_nfe_seu=22;controls_AN_oth=132;controls_AN_oth_female=54;controls_AN_oth_male=78;controls_AN_popmax=1730;controls_AN_raw=7448;controls_faf95=5.80532e-03;controls_faf95_afr=3.50000e-05;controls_faf95_amr=0.00000e+00;controls_faf95_eas=0.00000e+00;controls_faf95_nfe=1.29730e-02;controls_faf99=5.80584e-03;controls_faf99_afr=3.50000e-05;controls_faf99_amr=0.00000e+00;controls_faf99_eas=0.00000e+00;controls_faf99_nfe=1.29730e-02;controls_nhomalt=1;controls_nhomalt_afr=0;controls_nhomalt_afr_female=0;controls_nhomalt_afr_male=0;controls_nhomalt_amr=0;controls_nhomalt_amr_female=0;controls_nhomalt_amr_male=0;controls_nhomalt_asj=0;controls_nhomalt_asj_female=0;controls_nhomalt_asj_male=0;controls_nhomalt_eas=0;controls_nhomalt_eas_female=0;controls_nhomalt_eas_male=0;controls_nhomalt_female=0;controls_nhomalt_fin=0;controls_nhomalt_fin_female=0;controls_nhomalt_fin_male=0;controls_nhomalt_male=1;controls_nhomalt_nfe=1;controls_nhomalt_nfe_est=1;controls_nhomalt_nfe_female=0;controls_nhomalt_nfe_male=1;controls_nhomalt_nfe_nwe=0;controls_nhomalt_nfe_onf=0;controls_nhomalt_nfe_seu=0;controls_nhomalt_oth=0;controls_nhomalt_oth_female=0;controls_nhomalt_oth_male=0;controls_nhomalt_popmax=1;controls_nhomalt_raw=62;controls_popmax=nfe;dp_hist_all_bin_freq=5866|2696|3019|1797|1203|530|351|200|44|2|0|0|0|0|0|0|0|0|0|0;dp_hist_all_n_larger=0;dp_hist_alt_bin_freq=105|88|42|36|19|7|1|0|1|0|0|0|0|0|0|0|0|0|0|0;dp_hist_alt_n_larger=0;faf95=6.20402e-03;faf95_afr=9.51830e-04;faf95_amr=0.00000e+00;faf95_eas=0.00000e+00;faf95_nfe=1.19015e-02;faf99=6.20436e-03;faf99_afr=9.51800e-04;faf99_amr=0.00000e+00;faf99_eas=0.00000e+00;faf99_nfe=1.19019e-02;gq_hist_all_bin_freq=733|365|334|884|1212|712|1341|1137|486|786|600|284|692|154|319|137|261|37|160|421;gq_hist_alt_bin_freq=1|71|30|35|27|9|11|12|10|3|6|2|5|9|6|6|5|3|3|45;n_alt_alleles=1;nhomalt=9;nhomalt_afr=0;nhomalt_afr_female=0;nhomalt_afr_male=0;nhomalt_amr=0;nhomalt_amr_female=0;nhomalt_amr_male=0;nhomalt_asj=0;nhomalt_asj_female=0;nhomalt_asj_male=0;nhomalt_eas=0;nhomalt_eas_female=0;nhomalt_eas_male=0;nhomalt_female=1;nhomalt_fin=0;nhomalt_fin_female=0;nhomalt_fin_male=0;nhomalt_male=8;nhomalt_nfe=9;nhomalt_nfe_est=1;nhomalt_nfe_female=1;nhomalt_nfe_male=8;nhomalt_nfe_nwe=7;nhomalt_nfe_onf=1;nhomalt_nfe_seu=0;nhomalt_oth=0;nhomalt_oth_female=0;nhomalt_oth_male=0;nhomalt_popmax=9;nhomalt_raw=158;non_neuro_AC=72;non_neuro_AC_afr=4;non_neuro_AC_afr_female=1;non_neuro_AC_afr_male=3;non_neuro_AC_amr=0;non_neuro_AC_amr_female=0;non_neuro_AC_amr_male=0;non_neuro_AC_asj=1;non_neuro_AC_asj_female=0;non_neuro_AC_asj_male=1;non_neuro_AC_eas=0;non_neuro_AC_eas_female=0;non_neuro_AC_eas_male=0;non_neuro_AC_female=27;non_neuro_AC_fin=4;non_neuro_AC_fin_female=4;non_neuro_AC_fin_male=0;non_neuro_AC_male=45;non_neuro_AC_nfe=63;non_neuro_AC_nfe_est=25;non_neuro_AC_nfe_female=22;non_neuro_AC_nfe_male=41;non_neuro_AC_nfe_nwe=31;non_neuro_AC_nfe_onf=7;non_neuro_AC_nfe_seu=0;non_neuro_AC_oth=0;non_neuro_AC_oth_female=0;non_neuro_AC_oth_male=0;non_neuro_AC_popmax=63;non_neuro_AC_raw=370;non_neuro_AF=7.81250e-03;non_neuro_AF_afr=2.08551e-03;non_neuro_AF_afr_female=1.00000e-03;non_neuro_AF_afr_male=3.26797e-03;non_neuro_AF_amr=0.00000e+00;non_neuro_AF_amr_female=0.00000e+00;non_neuro_AF_amr_male=0.00000e+00;non_neuro_AF_asj=1.42857e-02;non_neuro_AF_asj_female=0.00000e+00;non_neuro_AF_asj_male=1.92308e-02;non_neuro_AF_eas=0.00000e+00;non_neuro_AF_eas_female=0.00000e+00;non_neuro_AF_eas_male=0.00000e+00;non_neuro_AF_female=6.77031e-03;non_neuro_AF_fin=1.00000e-02;non_neuro_AF_fin_female=1.94175e-02;non_neuro_AF_fin_male=0.00000e+00;non_neuro_AF_male=8.60750e-03;non_neuro_AF_nfe=1.25199e-02;non_neuro_AF_nfe_est=2.10438e-02;non_neuro_AF_nfe_female=1.05871e-02;non_neuro_AF_nfe_male=1.38795e-02;non_neuro_AF_nfe_nwe=9.88520e-03;non_neuro_AF_nfe_onf=1.02041e-02;non_neuro_AF_nfe_seu=0.00000e+00;non_neuro_AF_oth=0.00000e+00;non_neuro_AF_oth_female=0.00000e+00;non_neuro_AF_oth_male=0.00000e+00;non_neuro_AF_popmax=1.25199e-02;non_neuro_AF_raw=2.48991e-02;non_neuro_AN=9216;non_neuro_AN_afr=1918;non_neuro_AN_afr_female=1000;non_neuro_AN_afr_male=918;non_neuro_AN_amr=214;non_neuro_AN_amr_female=114;non_neuro_AN_amr_male=100;non_neuro_AN_asj=70;non_neuro_AN_asj_female=18;non_neuro_AN_asj_male=52;non_neuro_AN_eas=1310;non_neuro_AN_eas_female=464;non_neuro_AN_eas_male=846;non_neuro_AN_female=3988;non_neuro_AN_fin=400;non_neuro_AN_fin_female=206;non_neuro_AN_fin_male=194;non_neuro_AN_male=5228;non_neuro_AN_nfe=5032;non_neuro_AN_nfe_est=1188;non_neuro_AN_nfe_female=2078;non_neuro_AN_nfe_male=2954;non_neuro_AN_nfe_nwe=3136;non_neuro_AN_nfe_onf=686;non_neuro_AN_nfe_seu=22;non_neuro_AN_oth=272;non_neuro_AN_oth_female=108;non_neuro_AN_oth_male=164;non_neuro_AN_popmax=5032;non_neuro_AN_raw=14860;non_neuro_faf95=6.36172e-03;non_neuro_faf95_afr=7.12130e-04;non_neuro_faf95_amr=0.00000e+00;non_neuro_faf95_eas=0.00000e+00;non_neuro_faf95_nfe=1.00427e-02;non_neuro_faf99=6.36188e-03;non_neuro_faf99_afr=7.12220e-04;non_neuro_faf99_amr=0.00000e+00;non_neuro_faf99_eas=0.00000e+00;non_neuro_faf99_nfe=1.00423e-02;non_neuro_nhomalt=6;non_neuro_nhomalt_afr=0;non_neuro_nhomalt_afr_female=0;non_neuro_nhomalt_afr_male=0;non_neuro_nhomalt_amr=0;non_neuro_nhomalt_amr_female=0;non_neuro_nhomalt_amr_male=0;non_neuro_nhomalt_asj=0;non_neuro_nhomalt_asj_female=0;non_neuro_nhomalt_asj_male=0;non_neuro_nhomalt_eas=0;non_neuro_nhomalt_eas_female=0;non_neuro_nhomalt_eas_male=0;non_neuro_nhomalt_female=0;non_neuro_nhomalt_fin=0;non_neuro_nhomalt_fin_female=0;non_neuro_nhomalt_fin_male=0;non_neuro_nhomalt_male=6;non_neuro_nhomalt_nfe=6;non_neuro_nhomalt_nfe_est=1;non_neuro_nhomalt_nfe_female=0;non_neuro_nhomalt_nfe_male=6;non_neuro_nhomalt_nfe_nwe=5;non_neuro_nhomalt_nfe_onf=0;non_neuro_nhomalt_nfe_seu=0;non_neuro_nhomalt_oth=0;non_neuro_nhomalt_oth_female=0;non_neuro_nhomalt_oth_male=0;non_neuro_nhomalt_popmax=6;non_neuro_nhomalt_raw=131;non_neuro_popmax=nfe;non_topmed_AC=83;non_topmed_AC_afr=8;non_topmed_AC_afr_female=3;non_topmed_AC_afr_male=5;non_topmed_AC_amr=0;non_topmed_AC_amr_female=0;non_topmed_AC_amr_male=0;non_topmed_AC_asj=0;non_topmed_AC_asj_female=0;non_topmed_AC_asj_male=0;non_topmed_AC_eas=0;non_topmed_AC_eas_female=0;non_topmed_AC_eas_male=0;non_topmed_AC_female=38;non_topmed_AC_fin=10;non_topmed_AC_fin_female=8;non_topmed_AC_fin_male=2;non_topmed_AC_male=45;non_topmed_AC_nfe=65;non_topmed_AC_nfe_est=29;non_topmed_AC_nfe_female=27;non_topmed_AC_nfe_male=38;non_topmed_AC_nfe_nwe=27;non_topmed_AC_nfe_onf=9;non_topmed_AC_nfe_seu=0;non_topmed_AC_oth=0;non_topmed_AC_oth_female=0;non_topmed_AC_oth_male=0;non_topmed_AC_popmax=65;non_topmed_AC_raw=359;non_topmed_AF=6.98183e-03;non_topmed_AF_afr=1.65769e-03;non_topmed_AF_afr_female=1.46056e-03;non_topmed_AF_afr_male=1.80375e-03;non_topmed_AF_amr=0.00000e+00;non_topmed_AF_amr_female=0.00000e+00;non_topmed_AF_amr_male=0.00000e+00;non_topmed_AF_asj=0.00000e+00;non_topmed_AF_asj_female=0.00000e+00;non_topmed_AF_asj_male=0.00000e+00;non_topmed_AF_eas=0.00000e+00;non_topmed_AF_eas_female=0.00000e+00;non_topmed_AF_eas_male=0.00000e+00;non_topmed_AF_female=7.10015e-03;non_topmed_AF_fin=8.46024e-03;non_topmed_AF_fin_female=1.25392e-02;non_topmed_AF_fin_male=3.67647e-03;non_topmed_AF_male=6.88494e-03;non_topmed_AF_nfe=1.68394e-02;non_topmed_AF_nfe_est=2.17391e-02;non_topmed_AF_nfe_female=1.43008e-02;non_topmed_AF_nfe_male=1.92698e-02;non_topmed_AF_nfe_nwe=1.44077e-02;non_topmed_AF_nfe_onf=1.46580e-02;non_topmed_AF_nfe_seu=0.00000e+00;non_topmed_AF_oth=0.00000e+00;non_topmed_AF_oth_female=0.00000e+00;non_topmed_AF_oth_male=0.00000e+00;non_topmed_AF_popmax=1.68394e-02;non_topmed_AF_raw=1.90917e-02;non_topmed_AN=11888;non_topmed_AN_afr=4826;non_topmed_AN_afr_female=2054;non_topmed_AN_afr_male=2772;non_topmed_AN_amr=312;non_topmed_AN_amr_female=132;non_topmed_AN_amr_male=180;non_topmed_AN_asj=54;non_topmed_AN_asj_female=22;non_topmed_AN_asj_male=32;non_topmed_AN_eas=1278;non_topmed_AN_eas_female=448;non_topmed_AN_eas_male=830;non_topmed_AN_female=5352;non_topmed_AN_fin=1182;non_topmed_AN_fin_female=638;non_topmed_AN_fin_male=544;non_topmed_AN_male=6536;non_topmed_AN_nfe=3860;non_topmed_AN_nfe_est=1334;non_topmed_AN_nfe_female=1888;non_topmed_AN_nfe_male=1972;non_topmed_AN_nfe_nwe=1874;non_topmed_AN_nfe_onf=614;non_topmed_AN_nfe_seu=38;non_topmed_AN_oth=376;non_topmed_AN_oth_female=170;non_topmed_AN_oth_male=206;non_topmed_AN_popmax=3860;non_topmed_AN_raw=18804;non_topmed_faf95=5.77088e-03;non_topmed_faf95_afr=8.24500e-04;non_topmed_faf95_amr=0.00000e+00;non_topmed_faf95_eas=0.00000e+00;non_topmed_faf95_nfe=1.35569e-02;non_topmed_faf99=5.77022e-03;non_topmed_faf99_afr=8.24260e-04;non_topmed_faf99_amr=0.00000e+00;non_topmed_faf99_eas=0.00000e+00;non_topmed_faf99_nfe=1.35563e-02;non_topmed_nhomalt=6;non_topmed_nhomalt_afr=0;non_topmed_nhomalt_afr_female=0;non_topmed_nhomalt_afr_male=0;non_topmed_nhomalt_amr=0;non_topmed_nhomalt_amr_female=0;non_topmed_nhomalt_amr_male=0;non_topmed_nhomalt_asj=0;non_topmed_nhomalt_asj_female=0;non_topmed_nhomalt_asj_male=0;non_topmed_nhomalt_eas=0;non_topmed_nhomalt_eas_female=0;non_topmed_nhomalt_eas_male=0;non_topmed_nhomalt_female=1;non_topmed_nhomalt_fin=0;non_topmed_nhomalt_fin_female=0;non_topmed_nhomalt_fin_male=0;non_topmed_nhomalt_male=5;non_topmed_nhomalt_nfe=6;non_topmed_nhomalt_nfe_est=1;non_topmed_nhomalt_nfe_female=1;non_topmed_nhomalt_nfe_male=5;non_topmed_nhomalt_nfe_nwe=4;non_topmed_nhomalt_nfe_onf=1;non_topmed_nhomalt_nfe_seu=0;non_topmed_nhomalt_oth=0;non_topmed_nhomalt_oth_female=0;non_topmed_nhomalt_oth_male=0;non_topmed_nhomalt_popmax=6;non_topmed_nhomalt_raw=123;non_topmed_popmax=nfe;pab_max=1.00000e+00;popmax=nfe;rf_label=FP;rf_negative_label;rf_tp_probability=8.59190e-01;rf_train;segdup;variant_type=snv;vep=G|missense_variant|MODERATE|OR4F5|ENSG00000186092|Transcript|ENST00000335137|protein_coding|1/1||ENST00000335137.3:c.338T>G|ENSP00000334393.3:p.Phe113Cys|338|338|113|F/C|tTt/tGt|rs140739101|1||1||SNV|1|HGNC|14825|YES|||CCDS30547.1|ENSP00000334393|Q8NH21||UPI0000041BC1||deleterious(0.01)|possibly_damaging(0.568)|Transmembrane_helices:TMhelix&Prints_domain:PR00237&Superfamily_domains:SSF81321&Gene3D:1.20.1070.10&PROSITE_patterns:PS00237&hmmpanther:PTHR26451&hmmpanther:PTHR26451:SF72&PROSITE_profiles:PS50262||G:0.0190|G:0.004888|G:0.0015|G:0.036|G:0.003|G:0.0497|G:0.0153|G:0.0037|G:0.0457|G:0.008045|G:0.022|G:0.002553|G:0.02462|G:0|G:0.04624|G:0.04058|G:0.02716||||||||||||,G|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00000278218|open_chromatin_region||||||||||rs140739101|1||||SNV|1||||||||||||||||G:0.0190|G:0.004888|G:0.0015|G:0.036|G:0.003|G:0.0497|G:0.0153|G:0.0037|G:0.0457|G:0.008045|G:0.022|G:0.002553|G:0.02462|G:0|G:0.04624|G:0.04058|G:0.02716||||||||||||;CSQ=G|ENSG00000186092|Transcript|ENST00000641515|missense_variant||134|F/C|
            # v3 genomes (WGS)
            # chr1	69428	rs140739101	T	G	.	PASS	AC=430;AN=75720;AF=0.00567882;popmax=nfe;faf95_popmax=0.0100955;AC_non_v2_XX=174;AN_non_v2_XX=30666;AF_non_v2_XX=0.00567404;nhomalt_non_v2_XX=17;AC_non_cancer_fin_XX=5;AN_non_cancer_fin_XX=922;AF_non_cancer_fin_XX=0.00542299;nhomalt_non_cancer_fin_XX=0;AC_non_neuro_nfe=327;AN_non_neuro_nfe=30318;AF_non_neuro_nfe=0.0107857;nhomalt_non_neuro_nfe=35;AC_non_neuro_afr_XY=18;AN_non_neuro_afr_XY=8420;AF_non_neuro_afr_XY=0.00213777;nhomalt_non_neuro_afr_XY=0;AC_non_neuro_nfe_XY=143;AN_non_neuro_nfe_XY=12342;AF_non_neuro_nfe_XY=0.0115865;nhomalt_non_neuro_nfe_XY=17;AC_controls_and_biobanks_eas_XY=0;AN_controls_and_biobanks_eas_XY=1188;AF_controls_and_biobanks_eas_XY=0.00000;nhomalt_controls_and_biobanks_eas_XY=0;AC_non_neuro_sas_XX=0;AN_non_neuro_sas_XX=572;AF_non_neuro_sas_XX=0.00000;nhomalt_non_neuro_sas_XX=0;AC_non_v2=326;AN_non_v2=56262;AF_non_v2=0.00579432;nhomalt_non_v2=28;AC_non_topmed_nfe_XX=43;AN_non_topmed_nfe_XX=3788;AF_non_topmed_nfe_XX=0.0113516;nhomalt_non_topmed_nfe_XX=4;AC_non_v2_mid=0;AN_non_v2_mid=136;AF_non_v2_mid=0.00000;nhomalt_non_v2_mid=0;AC_non_topmed_sas=3;AN_non_topmed_sas=2298;AF_non_topmed_sas=0.00130548;nhomalt_non_topmed_sas=0;AC_non_cancer_eas_XX=0;AN_non_cancer_eas_XX=1814;AF_non_cancer_eas_XX=0.00000;nhomalt_non_cancer_eas_XX=0;AC_amr_XY=5;AN_amr_XY=3462;AF_amr_XY=0.00144425;nhomalt_amr_XY=1;AC_non_v2_nfe_XX=160;AN_non_v2_nfe_XX=15466;AF_non_v2_nfe_XX=0.0103453;nhomalt_non_v2_nfe_XX=17;AC_controls_and_biobanks_XY=37;AN_controls_and_biobanks_XY=9110;AF_controls_and_biobanks_XY=0.00406147;nhomalt_controls_and_biobanks_XY=3;AC_non_neuro_asj_XY=4;AN_non_neuro_asj_XY=610;AF_non_neuro_asj_XY=0.00655738;nhomalt_non_neuro_asj_XY=0;AC_oth=5;AN_oth=930;AF_oth=0.00537634;nhomalt_oth=0;AC_non_topmed_mid_XY=0;AN_non_topmed_mid_XY=56;AF_non_topmed_mid_XY=0.00000;nhomalt_non_topmed_mid_XY=0;AC_non_cancer_asj_XX=1;AN_non_cancer_asj_XX=652;AF_non_cancer_asj_XX=0.00153374;nhomalt_non_cancer_asj_XX=0;AC_sas_XY=3;AN_sas_XY=1742;AF_sas_XY=0.00172216;nhomalt_sas_XY=0;AC_non_neuro_fin=14;AN_non_neuro_fin=1952;AF_non_neuro_fin=0.00717213;nhomalt_non_neuro_fin=0;AC_non_topmed_amr_XY=3;AN_non_topmed_amr_XY=3010;AF_non_topmed_amr_XY=0.000996678;nhomalt_non_topmed_amr_XY=0;AC_non_neuro_XX=202;AN_non_neuro_XX=36494;AF_non_neuro_XX=0.00553516;nhomalt_non_neuro_XX=18;AC_fin_XX=5;AN_fin_XX=922;AF_fin_XX=0.00542299;nhomalt_fin_XX=0;AC_controls_and_biobanks_asj_XX=0;AN_controls_and_biobanks_asj_XX=34;AF_controls_and_biobanks_asj_XX=0.00000;nhomalt_controls_and_biobanks_asj_XX=0;AC_non_v2_raw=1395;AN_non_v2_raw=77212;AF_non_v2_raw=0.0180671;nhomalt_non_v2_raw=499;AC_non_v2_asj=4;AN_non_v2_asj=1200;AF_non_v2_asj=0.00333333;nhomalt_non_v2_asj=0;AC_nfe_XX=195;AN_nfe_XX=18612;AF_nfe_XX=0.0104771;nhomalt_nfe_XX=20;AC_controls_and_biobanks_raw=288;AN_controls_and_biobanks_raw=22576;AF_controls_and_biobanks_raw=0.0127569;nhomalt_controls_and_biobanks_raw=103;AC_controls_and_biobanks_ami=0;AN_controls_and_biobanks_ami=32;AF_controls_and_biobanks_ami=0.00000;nhomalt_controls_and_biobanks_ami=0;AC_non_topmed_eas=0;AN_non_topmed_eas=2930;AF_non_topmed_eas=0.00000;nhomalt_non_topmed_eas=0;AC_non_v2_amr=8;AN_non_v2_amr=5642;AF_non_v2_amr=0.00141794;nhomalt_non_v2_amr=1;AC_non_neuro_sas=3;AN_non_neuro_sas=2314;AF_non_neuro_sas=0.00129646;nhomalt_non_neuro_sas=0;AC_non_cancer_fin_XY=16;AN_non_cancer_fin_XY=2296;AF_non_cancer_fin_XY=0.00696864;nhomalt_non_cancer_fin_XY=0;AC_non_cancer_nfe_XY=146;AN_non_cancer_nfe_XY=12248;AF_non_cancer_nfe_XY=0.0119203;nhomalt_non_cancer_nfe_XY=15;AC_non_v2_oth=5;AN_non_v2_oth=832;AF_non_v2_oth=0.00600962;nhomalt_non_v2_oth=0;AC_ami=0;AN_ami=448;AF_ami=0.00000;nhomalt_ami=0;AC_non_cancer_XY=200;AN_non_cancer_XY=34492;AF_non_cancer_XY=0.00579845;nhomalt_non_cancer_XY=16;AC_non_v2_sas=2;AN_non_v2_sas=1758;AF_non_v2_sas=0.00113766;nhomalt_non_v2_sas=0;AC_non_topmed_afr_XX=8;AN_non_topmed_afr_XX=6648;AF_non_topmed_afr_XX=0.00120337;nhomalt_non_topmed_afr_XX=0;AC_sas=3;AN_sas=2314;AF_sas=0.00129646;nhomalt_sas=0;AC_non_neuro_nfe_XX=184;AN_non_neuro_nfe_XX=17976;AF_non_neuro_nfe_XX=0.0102359;nhomalt_non_neuro_nfe_XX=18;AC_non_topmed_ami_XX=0;AN_non_topmed_ami_XX=32;AF_non_topmed_ami_XX=0.00000;nhomalt_non_topmed_ami_XX=0;AC_ami_XY=0;AN_ami_XY=236;AF_ami_XY=0.00000;nhomalt_ami_XY=0;AC_oth_XX=2;AN_oth_XX=474;AF_oth_XX=0.00421941;nhomalt_oth_XX=0;AC_non_cancer_eas=0;AN_non_cancer_eas=3992;AF_non_cancer_eas=0.00000;nhomalt_non_cancer_eas=0;AC_non_topmed_XY=100;AN_non_topmed_XY=21668;AF_non_topmed_XY=0.00461510;nhomalt_non_topmed_XY=6;AC_non_v2_ami=0;AN_non_v2_ami=448;AF_non_v2_ami=0.00000;nhomalt_non_v2_ami=0;AC_non_neuro=389;AN_non_neuro=67672;AF_non_neuro=0.00574832;nhomalt_non_neuro=36;AC_amr_XX=4;AN_amr_XX=2908;AF_amr_XX=0.00137552;nhomalt_amr_XX=0;AC_controls_and_biobanks_nfe_XY=20;AN_controls_and_biobanks_nfe_XY=1652;AF_controls_and_biobanks_nfe_XY=0.0121065;nhomalt_controls_and_biobanks_nfe_XY=3;AC_controls_and_biobanks_eas=0;AN_controls_and_biobanks_eas=2070;AF_controls_and_biobanks_eas=0.00000;nhomalt_controls_and_biobanks_eas=0;AC_XX=218;AN_XX=39842;AF_XX=0.00547161;nhomalt_XX=20;AC_non_cancer_oth_XY=3;AN_non_cancer_oth_XY=422;AF_non_cancer_oth_XY=0.00710900;nhomalt_non_cancer_oth_XY=0;AC_non_v2_XY=152;AN_non_v2_XY=25596;AF_non_v2_XY=0.00593843;nhomalt_non_v2_XY=11;AC_non_topmed_amr_XX=4;AN_non_topmed_amr_XX=2228;AF_non_topmed_amr_XX=0.00179533;nhomalt_non_topmed_amr_XX=0;AC_fin=21;AN_fin=3218;AF_fin=0.00652579;nhomalt_fin=0;AC_controls_and_biobanks_nfe_XX=14;AN_controls_and_biobanks_nfe_XX=1414;AF_controls_and_biobanks_nfe_XX=0.00990099;nhomalt_controls_and_biobanks_nfe_XX=1;AC_controls_and_biobanks_afr=4;AN_controls_and_biobanks_afr=5128;AF_controls_and_biobanks_afr=0.000780031;nhomalt_controls_and_biobanks_afr=0;AC_asj_XX=1;AN_asj_XX=720;AF_asj_XX=0.00138889;nhomalt_asj_XX=0;AC_non_topmed_mid=0;AN_non_topmed_mid=126;AF_non_topmed_mid=0.00000;nhomalt_non_topmed_mid=0;AC_non_cancer_sas_XY=3;AN_non_cancer_sas_XY=1732;AF_non_cancer_sas_XY=0.00173210;nhomalt_non_cancer_sas_XY=0;AC_sas_XX=0;AN_sas_XX=572;AF_sas_XX=0.00000;nhomalt_sas_XX=0;AC_non_topmed=159;AN_non_topmed=37468;AF_non_topmed=0.00424362;nhomalt_non_topmed=10;AC_non_v2_oth_XX=2;AN_non_v2_oth_XX=438;AF_non_v2_oth_XX=0.00456621;nhomalt_non_v2_oth_XX=0;AC_non_neuro_ami_XY=0;AN_non_neuro_ami_XY=232;AF_non_neuro_ami_XY=0.00000;nhomalt_non_neuro_ami_XY=0;AC_controls_and_biobanks_afr_XY=2;AN_controls_and_biobanks_afr_XY=2494;AF_controls_and_biobanks_afr_XY=0.000801925;nhomalt_controls_and_biobanks_afr_XY=0;AC_controls_and_biobanks_amr_XX=2;AN_controls_and_biobanks_amr_XX=1160;AF_controls_and_biobanks_amr_XX=0.00172414;nhomalt_controls_and_biobanks_amr_XX=0;AC_non_topmed_amr=7;AN_non_topmed_amr=5238;AF_non_topmed_amr=0.00133639;nhomalt_non_topmed_amr=0;AC_controls_and_biobanks_sas_XX=0;AN_controls_and_biobanks_sas_XX=444;AF_controls_and_biobanks_sas_XX=0.00000;nhomalt_controls_and_biobanks_sas_XX=0;AC_controls_and_biobanks_amr=3;AN_controls_and_biobanks_amr=2188;AF_controls_and_biobanks_amr=0.00137112;nhomalt_controls_and_biobanks_amr=0;AC_non_neuro_fin_XX=3;AN_non_neuro_fin_XX=246;AF_non_neuro_fin_XX=0.0121951;nhomalt_non_neuro_fin_XX=0;AC_non_cancer_raw=1804;AN_non_cancer_raw=101230;AF_non_cancer_raw=0.0178208;nhomalt_non_cancer_raw=644;AC_non_neuro_mid=0;AN_non_neuro_mid=138;AF_non_neuro_mid=0.00000;nhomalt_non_neuro_mid=0;AC_non_v2_asj_XY=3;AN_non_v2_asj_XY=552;AF_non_v2_asj_XY=0.00543478;nhomalt_non_v2_asj_XY=0;AC_non_v2_afr=25;AN_non_v2_afr=17078;AF_non_v2_afr=0.00146387;nhomalt_non_v2_afr=0;AC_non_neuro_fin_XY=11;AN_non_neuro_fin_XY=1706;AF_non_neuro_fin_XY=0.00644783;nhomalt_non_neuro_fin_XY=0;AC_non_cancer_afr=35;AN_non_cancer_afr=24660;AF_non_cancer_afr=0.00141930;nhomalt_non_cancer_afr=0;AC_non_topmed_sas_XY=3;AN_non_topmed_sas_XY=1726;AF_non_topmed_sas_XY=0.00173812;nhomalt_non_topmed_sas_XY=0;AC_mid_XY=0;AN_mid_XY=58;AF_mid_XY=0.00000;nhomalt_mid_XY=0;AC_non_v2_oth_XY=3;AN_non_v2_oth_XY=394;AF_non_v2_oth_XY=0.00761421;nhomalt_non_v2_oth_XY=0;AC_controls_and_biobanks_fin=12;AN_controls_and_biobanks_fin=1524;AF_controls_and_biobanks_fin=0.00787402;nhomalt_controls_and_biobanks_fin=0;AC_non_neuro_eas_XY=0;AN_non_neuro_eas_XY=2298;AF_non_neuro_eas_XY=0.00000;nhomalt_non_neuro_eas_XY=0;AC_non_topmed_eas_XX=0;AN_non_topmed_eas_XX=1146;AF_non_topmed_eas_XX=0.00000;nhomalt_non_topmed_eas_XX=0;AC_non_v2_afr_XX=7;AN_non_v2_afr_XX=9470;AF_non_v2_afr_XX=0.000739176;nhomalt_non_v2_afr_XX=0;AC_non_neuro_amr_XX=4;AN_non_neuro_amr_XX=2834;AF_non_neuro_amr_XX=0.00141143;nhomalt_non_neuro_amr_XX=0;AC_non_cancer_ami=0;AN_non_cancer_ami=448;AF_non_cancer_ami=0.00000;nhomalt_non_cancer_ami=0;AC_XY=212;AN_XY=35878;AF_XY=0.00590891;nhomalt_XY=18;AC_non_topmed_asj_XX=0;AN_non_topmed_asj_XX=116;AF_non_topmed_asj_XX=0.00000;nhomalt_non_topmed_asj_XX=0;AC_non_topmed_eas_XY=0;AN_non_topmed_eas_XY=1784;AF_non_topmed_eas_XY=0.00000;nhomalt_non_topmed_eas_XY=0;AC_non_v2_eas_XY=0;AN_non_v2_eas_XY=1090;AF_non_v2_eas_XY=0.00000;nhomalt_non_v2_eas_XY=0;AC_eas=0;AN_eas=4204;AF_eas=0.00000;nhomalt_eas=0;AC_asj_XY=4;AN_asj_XY=626;AF_asj_XY=0.00638978;nhomalt_asj_XY=0;AC_non_v2_eas_XX=0;AN_non_v2_eas_XX=1052;AF_non_v2_eas_XX=0.00000;nhomalt_non_v2_eas_XX=0;AC_controls_and_biobanks_mid_XY=0;AN_controls_and_biobanks_mid_XY=46;AF_controls_and_biobanks_mid_XY=0.00000;nhomalt_controls_and_biobanks_mid_XY=0;AC_fin_XY=16;AN_fin_XY=2296;AF_fin_XY=0.00696864;nhomalt_fin_XY=0;AC_non_topmed_nfe=105;AN_non_topmed_nfe=8760;AF_non_topmed_nfe=0.0119863;nhomalt_non_topmed_nfe=10;AC_amr=9;AN_amr=6370;AF_amr=0.00141287;nhomalt_amr=1;AC_non_neuro_ami=0;AN_non_neuro_ami=430;AF_non_neuro_ami=0.00000;nhomalt_non_neuro_ami=0;AC_non_cancer_nfe_XX=187;AN_non_cancer_nfe_XX=17938;AF_non_cancer_nfe_XX=0.0104248;nhomalt_non_cancer_nfe_XX=18;AC_non_cancer_mid=0;AN_non_cancer_mid=132;AF_non_cancer_mid=0.00000;nhomalt_non_cancer_mid=0;AC_non_v2_mid_XY=0;AN_non_v2_mid_XY=54;AF_non_v2_mid_XY=0.00000;nhomalt_non_v2_mid_XY=0;AC_controls_and_biobanks_amr_XY=1;AN_controls_and_biobanks_amr_XY=1028;AF_controls_and_biobanks_amr_XY=0.000972763;nhomalt_controls_and_biobanks_amr_XY=0;AC_non_cancer_ami_XY=0;AN_non_cancer_ami_XY=236;AF_non_cancer_ami_XY=0.00000;nhomalt_non_cancer_ami_XY=0;AC_non_neuro_asj_XX=1;AN_non_neuro_asj_XX=712;AF_non_neuro_asj_XX=0.00140449;nhomalt_non_neuro_asj_XX=0;AC_afr=35;AN_afr=24884;AF_afr=0.00140653;nhomalt_afr=0;AC_non_v2_sas_XX=0;AN_non_v2_sas_XX=356;AF_non_v2_sas_XX=0.00000;nhomalt_non_v2_sas_XX=0;AC_non_neuro_afr_XX=8;AN_non_neuro_afr_XX=11510;AF_non_neuro_afr_XX=0.000695048;nhomalt_non_neuro_afr_XX=0;AC_non_cancer_sas=3;AN_non_cancer_sas=2294;AF_non_cancer_sas=0.00130776;nhomalt_non_cancer_sas=0;AC_non_topmed_fin=20;AN_non_topmed_fin=3180;AF_non_topmed_fin=0.00628931;nhomalt_non_topmed_fin=0;AC_non_cancer_asj_XY=3;AN_non_cancer_asj_XY=608;AF_non_cancer_asj_XY=0.00493421;nhomalt_non_cancer_asj_XY=0;AC_non_cancer_mid_XY=0;AN_non_cancer_mid_XY=52;AF_non_cancer_mid_XY=0.00000;nhomalt_non_cancer_mid_XY=0;AC_raw=1871;AN_raw=103974;AF_raw=0.0179949;nhomalt_raw=668;AC_non_topmed_XX=59;AN_non_topmed_XX=15800;AF_non_topmed_XX=0.00373418;nhomalt_non_topmed_XX=4;AC_ami_XX=0;AN_ami_XX=212;AF_ami_XX=0.00000;nhomalt_ami_XX=0;AC_eas_XY=0;AN_eas_XY=2298;AF_eas_XY=0.00000;nhomalt_eas_XY=0;AC_controls_and_biobanks_mid=0;AN_controls_and_biobanks_mid=114;AF_controls_and_biobanks_mid=0.00000;nhomalt_controls_and_biobanks_mid=0;AC_non_v2_nfe_XY=106;AN_non_v2_nfe_XY=9550;AF_non_v2_nfe_XY=0.0110995;nhomalt_non_v2_nfe_XY=10;AC_controls_and_biobanks_sas=2;AN_controls_and_biobanks_sas=1614;AF_controls_and_biobanks_sas=0.00123916;nhomalt_controls_and_biobanks_sas=0;AC_non_v2_eas=0;AN_non_v2_eas=2142;AF_non_v2_eas=0.00000;nhomalt_non_v2_eas=0;AC_mid=0;AN_mid=140;AF_mid=0.00000;nhomalt_mid=0;AC_oth_XY=3;AN_oth_XY=456;AF_oth_XY=0.00657895;nhomalt_oth_XY=0;AC_non_cancer_nfe=333;AN_non_cancer_nfe=30186;AF_non_cancer_nfe=0.0110316;nhomalt_non_cancer_nfe=33;AC_non_neuro_eas_XX=0;AN_non_neuro_eas_XX=1906;AF_non_neuro_eas_XX=0.00000;nhomalt_non_neuro_eas_XX=0;AC_non_neuro_sas_XY=3;AN_non_neuro_sas_XY=1742;AF_non_neuro_sas_XY=0.00172216;nhomalt_non_neuro_sas_XY=0;AC_non_cancer_ami_XX=0;AN_non_cancer_ami_XX=212;AF_non_cancer_ami_XX=0.00000;nhomalt_non_cancer_ami_XX=0;AC_mid_XX=0;AN_mid_XX=82;AF_mid_XX=0.00000;nhomalt_mid_XX=0;AC_non_topmed_asj=2;AN_non_topmed_asj=396;AF_non_topmed_asj=0.00505051;nhomalt_non_topmed_asj=0;AC_non_v2_asj_XX=1;AN_non_v2_asj_XX=648;AF_non_v2_asj_XX=0.00154321;nhomalt_non_v2_asj_XX=0;nhomalt=38;AC_non_v2_amr_XY=5;AN_non_v2_amr_XY=3080;AF_non_v2_amr_XY=0.00162338;nhomalt_non_v2_amr_XY=1;AC_non_cancer_amr_XX=4;AN_non_cancer_amr_XX=2872;AF_non_cancer_amr_XX=0.00139276;nhomalt_non_cancer_amr_XX=0;AC_controls_and_biobanks_afr_XX=2;AN_controls_and_biobanks_afr_XX=2634;AF_controls_and_biobanks_afr_XX=0.000759301;nhomalt_controls_and_biobanks_afr_XX=0;AC_asj=5;AN_asj=1346;AF_asj=0.00371471;nhomalt_asj=0;AC_non_topmed_asj_XY=2;AN_non_topmed_asj_XY=280;AF_non_topmed_asj_XY=0.00714286;nhomalt_non_topmed_asj_XY=0;AC_non_v2_fin_XX=1;AN_non_v2_fin_XX=380;AF_non_v2_fin_XX=0.00263158;nhomalt_non_v2_fin_XX=0;AC_non_topmed_ami=0;AN_non_topmed_ami=50;AF_non_topmed_ami=0.00000;nhomalt_non_topmed_ami=0;AC_controls_and_biobanks_eas_XX=0;AN_controls_and_biobanks_eas_XX=882;AF_controls_and_biobanks_eas_XX=0.00000;nhomalt_controls_and_biobanks_eas_XX=0;AC_controls_and_biobanks_fin_XX=2;AN_controls_and_biobanks_fin_XX=206;AF_controls_and_biobanks_fin_XX=0.00970874;nhomalt_controls_and_biobanks_fin_XX=0;AC_non_topmed_raw=768;AN_non_topmed_raw=54594;AF_non_topmed_raw=0.0140675;nhomalt_non_topmed_raw=273;AC_non_cancer_eas_XY=0;AN_non_cancer_eas_XY=2178;AF_non_cancer_eas_XY=0.00000;nhomalt_non_cancer_eas_XY=0;AC_non_cancer=410;AN_non_cancer=73334;AF_non_cancer=0.00559086;nhomalt_non_cancer=34;AC_controls_and_biobanks_ami_XY=0;AN_controls_and_biobanks_ami_XY=14;AF_controls_and_biobanks_ami_XY=0.00000;nhomalt_controls_and_biobanks_ami_XY=0;AC_controls_and_biobanks_mid_XX=0;AN_controls_and_biobanks_mid_XX=68;AF_controls_and_biobanks_mid_XX=0.00000;nhomalt_controls_and_biobanks_mid_XX=0;AC_non_v2_afr_XY=18;AN_non_v2_afr_XY=7608;AF_non_v2_afr_XY=0.00236593;nhomalt_non_v2_afr_XY=0;AC_non_v2_sas_XY=2;AN_non_v2_sas_XY=1402;AF_non_v2_sas_XY=0.00142653;nhomalt_non_v2_sas_XY=0;AC_non_v2_fin=16;AN_non_v2_fin=2010;AF_non_v2_fin=0.00796020;nhomalt_non_v2_fin=0;AC_non_neuro_oth=5;AN_non_neuro_oth=892;AF_non_neuro_oth=0.00560538;nhomalt_non_neuro_oth=0;AC_non_cancer_sas_XX=0;AN_non_cancer_sas_XX=562;AF_non_cancer_sas_XX=0.00000;nhomalt_non_cancer_sas_XX=0;AC_non_neuro_asj=5;AN_non_neuro_asj=1322;AF_non_neuro_asj=0.00378215;nhomalt_non_neuro_asj=0;AC_non_topmed_afr=19;AN_non_topmed_afr=13836;AF_non_topmed_afr=0.00137323;nhomalt_non_topmed_afr=0;AC_non_topmed_afr_XY=11;AN_non_topmed_afr_XY=7188;AF_non_topmed_afr_XY=0.00153033;nhomalt_non_topmed_afr_XY=0;AC_non_neuro_eas=0;AN_non_neuro_eas=4204;AF_non_neuro_eas=0.00000;nhomalt_non_neuro_eas=0;AC_afr_XX=11;AN_afr_XX=13434;AF_afr_XX=0.000818818;nhomalt_afr_XX=0;AC_non_neuro_mid_XY=0;AN_non_neuro_mid_XY=56;AF_non_neuro_mid_XY=0.00000;nhomalt_non_neuro_mid_XY=0;AC_non_topmed_fin_XX=4;AN_non_topmed_fin_XX=900;AF_non_topmed_fin_XX=0.00444444;nhomalt_non_topmed_fin_XX=0;AC_non_cancer_amr=9;AN_non_cancer_amr=6258;AF_non_cancer_amr=0.00143816;nhomalt_non_cancer_amr=1;AC_non_v2_ami_XX=0;AN_non_v2_ami_XX=212;AF_non_v2_ami_XX=0.00000;nhomalt_non_v2_ami_XX=0;AC_afr_XY=24;AN_afr_XY=11450;AF_afr_XY=0.00209607;nhomalt_afr_XY=0;AC_non_v2_mid_XX=0;AN_non_v2_mid_XX=82;AF_non_v2_mid_XX=0.00000;nhomalt_non_v2_mid_XX=0;AC_non_topmed_fin_XY=16;AN_non_topmed_fin_XY=2280;AF_non_topmed_fin_XY=0.00701754;nhomalt_non_topmed_fin_XY=0;AC_non_neuro_amr_XY=5;AN_non_neuro_amr_XY=3338;AF_non_neuro_amr_XY=0.00149790;nhomalt_non_neuro_amr_XY=1;AC_non_topmed_mid_XX=0;AN_non_topmed_mid_XX=70;AF_non_topmed_mid_XX=0.00000;nhomalt_non_topmed_mid_XX=0;AC_controls_and_biobanks_asj_XY=0;AN_controls_and_biobanks_asj_XY=14;AF_controls_and_biobanks_asj_XY=0.00000;nhomalt_controls_and_biobanks_asj_XY=0;AC_non_v2_fin_XY=15;AN_non_v2_fin_XY=1630;AF_non_v2_fin_XY=0.00920245;nhomalt_non_v2_fin_XY=0;AC_controls_and_biobanks_ami_XX=0;AN_controls_and_biobanks_ami_XX=18;AF_controls_and_biobanks_ami_XX=0.00000;nhomalt_controls_and_biobanks_ami_XX=0;AC_eas_XX=0;AN_eas_XX=1906;AF_eas_XX=0.00000;nhomalt_eas_XX=0;AC_non_cancer_amr_XY=5;AN_non_cancer_amr_XY=3386;AF_non_cancer_amr_XY=0.00147667;nhomalt_non_cancer_amr_XY=1;AC_non_neuro_ami_XX=0;AN_non_neuro_ami_XX=198;AF_non_neuro_ami_XX=0.00000;nhomalt_non_neuro_ami_XX=0;AC_controls_and_biobanks=57;AN_controls_and_biobanks=16160;AF_controls_and_biobanks=0.00352723;nhomalt_controls_and_biobanks=4;AC_controls_and_biobanks_oth=2;AN_controls_and_biobanks_oth=376;AF_controls_and_biobanks_oth=0.00531915;nhomalt_controls_and_biobanks_oth=0;AC_nfe_XY=157;AN_nfe_XY=13254;AF_nfe_XY=0.0118455;nhomalt_nfe_XY=17;AC_non_cancer_afr_XX=11;AN_non_cancer_afr_XX=13326;AF_non_cancer_afr_XX=0.000825454;nhomalt_non_cancer_afr_XX=0;AC_controls_and_biobanks_sas_XY=2;AN_controls_and_biobanks_sas_XY=1170;AF_controls_and_biobanks_sas_XY=0.00170940;nhomalt_controls_and_biobanks_sas_XY=0;AC_non_cancer_oth=5;AN_non_cancer_oth=886;AF_non_cancer_oth=0.00564334;nhomalt_non_cancer_oth=0;AC_non_topmed_oth=3;AN_non_topmed_oth=654;AF_non_topmed_oth=0.00458716;nhomalt_non_topmed_oth=0;AC_non_topmed_nfe_XY=62;AN_non_topmed_nfe_XY=4972;AF_non_topmed_nfe_XY=0.0124698;nhomalt_non_topmed_nfe_XY=6;AC_non_topmed_sas_XX=0;AN_non_topmed_sas_XX=572;AF_non_topmed_sas_XX=0.00000;nhomalt_non_topmed_sas_XX=0;AC_non_v2_nfe=266;AN_non_v2_nfe=25016;AF_non_v2_nfe=0.0106332;nhomalt_non_v2_nfe=27;AC_non_topmed_oth_XX=0;AN_non_topmed_oth_XX=300;AF_non_topmed_oth_XX=0.00000;nhomalt_non_topmed_oth_XX=0;AC_non_cancer_mid_XX=0;AN_non_cancer_mid_XX=80;AF_non_cancer_mid_XX=0.00000;nhomalt_non_cancer_mid_XX=0;AC_controls_and_biobanks_nfe=34;AN_controls_and_biobanks_nfe=3066;AF_controls_and_biobanks_nfe=0.0110894;nhomalt_controls_and_biobanks_nfe=4;AC_controls_and_biobanks_oth_XY=2;AN_controls_and_biobanks_oth_XY=186;AF_controls_and_biobanks_oth_XY=0.0107527;nhomalt_controls_and_biobanks_oth_XY=0;AC_controls_and_biobanks_fin_XY=10;AN_controls_and_biobanks_fin_XY=1318;AF_controls_and_biobanks_fin_XY=0.00758725;nhomalt_controls_and_biobanks_fin_XY=0;AC_non_v2_amr_XX=3;AN_non_v2_amr_XX=2562;AF_non_v2_amr_XX=0.00117096;nhomalt_non_v2_amr_XX=0;AC_non_cancer_asj=4;AN_non_cancer_asj=1260;AF_non_cancer_asj=0.00317460;nhomalt_non_cancer_asj=0;AC_non_cancer_oth_XX=2;AN_non_cancer_oth_XX=464;AF_non_cancer_oth_XX=0.00431034;nhomalt_non_cancer_oth_XX=0;AC_non_neuro_amr=9;AN_non_neuro_amr=6172;AF_non_neuro_amr=0.00145820;nhomalt_non_neuro_amr=1;AC_non_cancer_XX=210;AN_non_cancer_XX=38842;AF_non_cancer_XX=0.00540652;nhomalt_non_cancer_XX=18;AC_non_v2_ami_XY=0;AN_non_v2_ami_XY=236;AF_non_v2_ami_XY=0.00000;nhomalt_non_v2_ami_XY=0;AC_non_neuro_raw=1701;AN_non_neuro_raw=91914;AF_non_neuro_raw=0.0185064;nhomalt_non_neuro_raw=606;AC_non_neuro_afr=26;AN_non_neuro_afr=19930;AF_non_neuro_afr=0.00130457;nhomalt_non_neuro_afr=0;AC_non_topmed_ami_XY=0;AN_non_topmed_ami_XY=18;AF_non_topmed_ami_XY=0.00000;nhomalt_non_topmed_ami_XY=0;AC_non_neuro_oth_XY=3;AN_non_neuro_oth_XY=434;AF_non_neuro_oth_XY=0.00691244;nhomalt_non_neuro_oth_XY=0;AC_non_neuro_oth_XX=2;AN_non_neuro_oth_XX=458;AF_non_neuro_oth_XX=0.00436681;nhomalt_non_neuro_oth_XX=0;AC_controls_and_biobanks_XX=20;AN_controls_and_biobanks_XX=7050;AF_controls_and_biobanks_XX=0.00283688;nhomalt_controls_and_biobanks_XX=1;AC_non_cancer_afr_XY=24;AN_non_cancer_afr_XY=11334;AF_non_cancer_afr_XY=0.00211752;nhomalt_non_cancer_afr_XY=0;AC_non_cancer_fin=21;AN_non_cancer_fin=3218;AF_non_cancer_fin=0.00652579;nhomalt_non_cancer_fin=0;AC_controls_and_biobanks_asj=0;AN_controls_and_biobanks_asj=48;AF_controls_and_biobanks_asj=0.00000;nhomalt_controls_and_biobanks_asj=0;AC_non_topmed_oth_XY=3;AN_non_topmed_oth_XY=354;AF_non_topmed_oth_XY=0.00847458;nhomalt_non_topmed_oth_XY=0;AC_non_neuro_mid_XX=0;AN_non_neuro_mid_XX=82;AF_non_neuro_mid_XX=0.00000;nhomalt_non_neuro_mid_XX=0;AC_controls_and_biobanks_oth_XX=0;AN_controls_and_biobanks_oth_XX=190;AF_controls_and_biobanks_oth_XX=0.00000;nhomalt_controls_and_biobanks_oth_XX=0;AC_non_neuro_XY=187;AN_non_neuro_XY=31178;AF_non_neuro_XY=0.00599782;nhomalt_non_neuro_XY=18;AC_nfe=352;AN_nfe=31866;AF_nfe=0.0110463;nhomalt_nfe=37;AC_popmax=352;AN_popmax=31866;AF_popmax=0.0110463;nhomalt_popmax=37;faf95_sas=0.000352970;faf99_sas=0.000187690;faf95_eas=0.00000;faf99_eas=0.00000;faf95_amr=0.000736940;faf99_amr=0.000550510;faf95_afr=0.00103881;faf99_afr=0.000912740;faf95=0.00523513;faf99=0.00506114;faf95_nfe=0.0100955;faf99_nfe=0.00972244;age_hist_het_bin_freq=11|5|9|12|20|16|21|15|9|3;age_hist_het_n_smaller=15;age_hist_het_n_larger=1;age_hist_hom_bin_freq=1|1|0|2|3|3|1|2|1|0;age_hist_hom_n_smaller=2;age_hist_hom_n_larger=1;FS=4.94928;MQ=31.7700;MQRankSum=-1.07600;QUALapprox=288675;QD=12.6131;ReadPosRankSum=0.715000;VarDP=22887;AS_FS=4.94928;AS_MQ=31.7625;AS_MQRankSum=-1.09200;AS_pab_max=1.00000;AS_QUALapprox=288641;AS_QD=12.6143;AS_ReadPosRankSum=0.711000;AS_SB_TABLE=6448,4706|5893,5838;AS_SOR=0.433605;InbreedingCoeff=0.708817;AS_culprit=AS_MQ;AS_VQSLOD=-2.49130;NEGATIVE_TRAIN_SITE;allele_type=snv;n_alt_alleles=2;variant_type=multi-snv;segdup;gq_hist_alt_bin_freq=0|0|0|0|2|9|16|13|10|22|15|20|12|14|18|22|19|8|9|183;gq_hist_all_bin_freq=0|0|0|0|22737|7005|4472|1754|685|460|220|100|69|45|39|33|27|11|17|186;dp_hist_alt_bin_freq=0|0|155|140|59|26|9|3|0|0|0|0|0|0|0|0|0|0|0|0;dp_hist_alt_n_larger=0;dp_hist_all_bin_freq=0|0|14696|14332|5329|2524|897|70|9|0|0|0|2|0|1|0|0|0|0|0;dp_hist_all_n_larger=0;ab_hist_alt_bin_freq=0|0|0|0|69|78|71|50|32|14|26|8|3|1|1|1|0|0|0|0;cadd_raw_score=3.44120;cadd_phred=24.5000;revel_score=0.0550000;splice_ai_max_ds=0.00000;splice_ai_consequence=no_consequence;primate_ai_score=0.371419;vep=G|missense_variant|MODERATE|OR4F5|ENSG00000186092|Transcript|ENST00000335137|protein_coding|1/1||ENST00000335137.4:c.338T>G|ENSP00000334393.3:p.Phe113Cys|374|338|113|F/C|tTt/tGt|1||1|SNV||HGNC|HGNC:14825|YES||P1|CCDS30547.1|ENSP00000334393|||||deleterious(0.010)|probably_damaging(0.984)|Prints:PR00237&Gene3D:1&Pfam:PF13853&PROSITE_patterns:PS00237&PROSITE_profiles:PS50262&Superfamily:SSF81321&Transmembrane_helices:TMhelix&CDD:cd15226|||||||||,G|missense_variant|MODERATE|OR4F5|ENSG00000186092|Transcript|ENST00000641515|protein_coding|3/3||ENST00000641515.2:c.401T>G|ENSP00000493376.2:p.Phe134Cys|461|401|134|F/C|tTt/tGt|1||1|SNV||HGNC|HGNC:14825|||||ENSP00000493376|||||deleterious(0.020)|probably_damaging(1.000)|Transmembrane_helices:TMhelix&CDD:cd15226&PANTHER:PTHR26451&PANTHER:PTHR26451&Pfam:PF13853&PROSITE_profiles:PS50262&Gene3D:1&PROSITE_patterns:PS00237&Superfamily:SSF81321|||||||||,G|missense_variant|MODERATE|OR4F5|79501|Transcript|NM_001005484.1|protein_coding|1/1||NM_001005484.1:c.338T>G|NP_001005484.1:p.Phe113Cys|338|338|113|F/C|tTt/tGt|1||1|SNV||EntrezGene|HGNC:14825|YES||||NP_001005484.1|||||deleterious(0.010)|probably_damaging(0.984)||||||||||,G|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00000918279|TF_binding_site||||||||||1|||SNV||||||||||||||||||||||||;CSQ=G|ENSG00000186092|Transcript|ENST00000641515|missense_variant||134|F/C|

            # Fields that exist in v2 exomes, v2 genomes and v3 genomes for both a coding and a noncoding variant (60), with example values from v3 genomes:
            # AC=430
            # AC_afr=35
            # AC_amr=9
            # AC_asj=5
            # AC_eas=0
            # AC_fin=21
            # AC_nfe=352
            # AC_oth=5
            # AC_raw=1871
            # AF_raw=0.0179949
            # AN=75720
            # AN_afr=24884
            # AN_amr=6370
            # AN_asj=1346
            # AN_eas=4204
            # AN_fin=3218
            # AN_nfe=31866
            # AN_oth=930
            # AN_raw=103974
            # FS=4.94928
            # AS_FS=4.94928
            # InbreedingCoeff=0.708817
            # MQ=31.7700
            # AS_MQ=31.7625
            # QD=12.6131
            # AS_QD=12.6143
            # ab_hist_alt_bin_freq=0|0|0|0|69|78|71|50|32|14|26|8|3|1|1|1|0|0|0|0
            # age_hist_het_bin_freq=11|5|9|12|20|16|21|15|9|3
            # age_hist_het_n_larger=1
            # age_hist_het_n_smaller=15
            # age_hist_hom_bin_freq=1|1|0|2|3|3|1|2|1|0
            # age_hist_hom_n_larger=1
            # age_hist_hom_n_smaller=2
            # allele_type=snv
            # dp_hist_all_bin_freq=0|0|14696|14332|5329|2524|897|70|9|0|0|0|2|0|1|0|0|0|0|0
            # dp_hist_all_n_larger=0
            # dp_hist_alt_bin_freq=0|0|155|140|59|26|9|3|0|0|0|0|0|0|0|0|0|0|0|0
            # dp_hist_alt_n_larger=0
            # faf95=0.00523513
            # faf95_afr=0.00103881
            # faf95_amr=0.000736940
            # faf95_eas=0.00000
            # faf95_nfe=0.0100955
            # faf99=0.00506114
            # faf99_afr=0.000912740
            # faf99_amr=0.000550510
            # faf99_eas=0.00000
            # faf99_nfe=0.00972244
            # gq_hist_all_bin_freq=0|0|0|0|22737|7005|4472|1754|685|460|220|100|69|45|39|33|27|11|17|186
            # gq_hist_alt_bin_freq=0|0|0|0|2|9|16|13|10|22|15|20|12|14|18|22|19|8|9|183
            # n_alt_alleles=2
            # nhomalt=38
            # nhomalt_afr=0
            # nhomalt_amr=1
            # nhomalt_asj=0
            # nhomalt_eas=0
            # nhomalt_fin=0
            # nhomalt_nfe=37
            # nhomalt_oth=0
            # nhomalt_raw=668
            # segdup
            # variant_type=multi-snv
            # vep=G|missense_variant|MODERATE|OR4F5|ENSG00000186092|Transcript|ENST00000335137|protein_coding|1/1||ENST00000335137.4:c.338T>G|ENSP00000334393.3:p.Phe113Cys|374|338|113|F/C|tTt/tGt|1||1|SNV||HGNC|HGNC:14825|YES||P1|CCDS30547.1|ENSP00000334393|||||deleterious(0.010)|probably_damaging(0.984)|Prints:PR00237&Gene3D:1&Pfam:PF13853&PROSITE_patterns:PS00237&PROSITE_profiles:PS50262&Superfamily:SSF81321&Transmembrane_helices:TMhelix&CDD:cd15226|||||||||,G|missense_variant|MODERATE|OR4F5|ENSG00000186092|Transcript|ENST00000641515|protein_coding|3/3||ENST00000641515.2:c.401T>G|ENSP00000493376.2:p.Phe134Cys|461|401|134|F/C|tTt/tGt|1||1|SNV||HGNC|HGNC:14825|||||ENSP00000493376|||||deleterious(0.020)|probably_damaging(1.000)|Transmembrane_helices:TMhelix&CDD:cd15226&PANTHER:PTHR26451&PANTHER:PTHR26451&Pfam:PF13853&PROSITE_profiles:PS50262&Gene3D:1&PROSITE_patterns:PS00237&Superfamily:SSF81321|||||||||,G|missense_variant|MODERATE|OR4F5|79501|Transcript|NM_001005484.1|protein_coding|1/1||NM_001005484.1:c.338T>G|NP_001005484.1:p.Phe113Cys|338|338|113|F/C|tTt/tGt|1||1|SNV||EntrezGene|HGNC:14825|YES||||NP_001005484.1|||||deleterious(0.010)|probably_damaging(0.984)||||||||||,G|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00000918279|TF_binding_site||||||||||1|||SNV||||||||||||||||||||||||

            # Fields that exist in v2 exomes, v2 genomes and v3 genomes only for a coding variant (16), with example values from v3 genomes:
            # AC_popmax=352
            # AF=0.00567882
            # AF_afr=0.00140653
            # AF_amr=0.00141287
            # AF_asj=0.00371471
            # AF_eas=0.00000
            # AF_fin=0.00652579
            # AF_nfe=0.0110463
            # AF_oth=0.00537634
            # AF_popmax=0.0110463
            # AN_popmax=31866
            # CSQ=G|ENSG00000186092|Transcript|ENST00000641515|missense_variant||134|F/C|
            # MQRankSum=-1.07600
            # AS_MQRankSum=-1.09200
            # ReadPosRankSum=0.715000
            # AS_ReadPosRankSum=0.711000
            # nhomalt_popmax=37
            # popmax=nfe
            # faf95_popmax=0.0100955
            # AC_popmax=352
            # AN_popmax=31866
            # AF_popmax=0.0110463
            # nhomalt_popmax=37

            # v3 has a lot more interesting fields (see "Interesting fields" above), so limiting ourselves to fields present in both v2 and v3 makes no sense. Should probably parse all fields.
            # On the other hand, there are ~900 fields if I were to parse all. Should probably limit it to fields of interest (but include all the subpopulations for these).
            


            # v4 exomes fields (for a missense variant):
            # zcat gnomad.exomes.v4.0.sites.chr1.vcf.bgz | head -n 10000 | g -v "^#" | g missense | head -n1
            # chr1	65571	.	A	C	.	AC0;AS_VQSR	AC=0;AN=6768;AF=0.00000;AC_XX=0;AF_XX=0.00000;AN_XX=3126;nhomalt_XX=0;AC_XY=0;AF_XY=0.00000;AN_XY=3642;nhomalt_XY=0;nhomalt=0;AC_afr_XX=0;AF_afr_XX=0.00000;AN_afr_XX=196;nhomalt_afr_XX=0;AC_afr_XY=0;AF_afr_XY=0.00000;AN_afr_XY=162;nhomalt_afr_XY=0;AC_afr=0;AF_afr=0.00000;AN_afr=358;nhomalt_afr=0;AC_amr_XX=0;AF_amr_XX=0.00000;AN_amr_XX=166;nhomalt_amr_XX=0;AC_amr_XY=0;AF_amr_XY=0.00000;AN_amr_XY=108;nhomalt_amr_XY=0;AC_amr=0;AF_amr=0.00000;AN_amr=274;nhomalt_amr=0;AC_asj_XX=0;AF_asj_XX=0.00000;AN_asj_XX=48;nhomalt_asj_XX=0;AC_asj_XY=0;AF_asj_XY=0.00000;AN_asj_XY=58;nhomalt_asj_XY=0;AC_asj=0;AF_asj=0.00000;AN_asj=106;nhomalt_asj=0;AC_eas_XX=0;AF_eas_XX=0.00000;AN_eas_XX=266;nhomalt_eas_XX=0;AC_eas_XY=0;AF_eas_XY=0.00000;AN_eas_XY=236;nhomalt_eas_XY=0;AC_eas=0;AF_eas=0.00000;AN_eas=502;nhomalt_eas=0;AC_fin_XX=0;AF_fin_XX=0.00000;AN_fin_XX=112;nhomalt_fin_XX=0;AC_fin_XY=0;AF_fin_XY=0.00000;AN_fin_XY=116;nhomalt_fin_XY=0;AC_fin=0;AF_fin=0.00000;AN_fin=228;nhomalt_fin=0;AC_mid_XX=0;AF_mid_XX=0.00000;AN_mid_XX=8;nhomalt_mid_XX=0;AC_mid_XY=0;AF_mid_XY=0.00000;AN_mid_XY=16;nhomalt_mid_XY=0;AC_mid=0;AF_mid=0.00000;AN_mid=24;nhomalt_mid=0;AC_nfe_XX=0;AF_nfe_XX=0.00000;AN_nfe_XX=1790;nhomalt_nfe_XX=0;AC_nfe_XY=0;AF_nfe_XY=0.00000;AN_nfe_XY=1470;nhomalt_nfe_XY=0;AC_nfe=0;AF_nfe=0.00000;AN_nfe=3260;nhomalt_nfe=0;AC_non_ukb_XX=0;AF_non_ukb_XX=0.00000;AN_non_ukb_XX=3126;nhomalt_non_ukb_XX=0;AC_non_ukb_XY=0;AF_non_ukb_XY=0.00000;AN_non_ukb_XY=3642;nhomalt_non_ukb_XY=0;AC_non_ukb=0;AF_non_ukb=0.00000;AN_non_ukb=6768;nhomalt_non_ukb=0;AC_non_ukb_afr_XX=0;AF_non_ukb_afr_XX=0.00000;AN_non_ukb_afr_XX=196;nhomalt_non_ukb_afr_XX=0;AC_non_ukb_afr_XY=0;AF_non_ukb_afr_XY=0.00000;AN_non_ukb_afr_XY=162;nhomalt_non_ukb_afr_XY=0;AC_non_ukb_afr=0;AF_non_ukb_afr=0.00000;AN_non_ukb_afr=358;nhomalt_non_ukb_afr=0;AC_non_ukb_amr_XX=0;AF_non_ukb_amr_XX=0.00000;AN_non_ukb_amr_XX=166;nhomalt_non_ukb_amr_XX=0;AC_non_ukb_amr_XY=0;AF_non_ukb_amr_XY=0.00000;AN_non_ukb_amr_XY=108;nhomalt_non_ukb_amr_XY=0;AC_non_ukb_amr=0;AF_non_ukb_amr=0.00000;AN_non_ukb_amr=274;nhomalt_non_ukb_amr=0;AC_non_ukb_asj_XX=0;AF_non_ukb_asj_XX=0.00000;AN_non_ukb_asj_XX=48;nhomalt_non_ukb_asj_XX=0;AC_non_ukb_asj_XY=0;AF_non_ukb_asj_XY=0.00000;AN_non_ukb_asj_XY=58;nhomalt_non_ukb_asj_XY=0;AC_non_ukb_asj=0;AF_non_ukb_asj=0.00000;AN_non_ukb_asj=106;nhomalt_non_ukb_asj=0;AC_non_ukb_eas_XX=0;AF_non_ukb_eas_XX=0.00000;AN_non_ukb_eas_XX=266;nhomalt_non_ukb_eas_XX=0;AC_non_ukb_eas_XY=0;AF_non_ukb_eas_XY=0.00000;AN_non_ukb_eas_XY=236;nhomalt_non_ukb_eas_XY=0;AC_non_ukb_eas=0;AF_non_ukb_eas=0.00000;AN_non_ukb_eas=502;nhomalt_non_ukb_eas=0;AC_non_ukb_fin_XX=0;AF_non_ukb_fin_XX=0.00000;AN_non_ukb_fin_XX=112;nhomalt_non_ukb_fin_XX=0;AC_non_ukb_fin_XY=0;AF_non_ukb_fin_XY=0.00000;AN_non_ukb_fin_XY=116;nhomalt_non_ukb_fin_XY=0;AC_non_ukb_fin=0;AF_non_ukb_fin=0.00000;AN_non_ukb_fin=228;nhomalt_non_ukb_fin=0;AC_non_ukb_mid_XX=0;AF_non_ukb_mid_XX=0.00000;AN_non_ukb_mid_XX=8;nhomalt_non_ukb_mid_XX=0;AC_non_ukb_mid_XY=0;AF_non_ukb_mid_XY=0.00000;AN_non_ukb_mid_XY=16;nhomalt_non_ukb_mid_XY=0;AC_non_ukb_mid=0;AF_non_ukb_mid=0.00000;AN_non_ukb_mid=24;nhomalt_non_ukb_mid=0;AC_non_ukb_nfe_XX=0;AF_non_ukb_nfe_XX=0.00000;AN_non_ukb_nfe_XX=1790;nhomalt_non_ukb_nfe_XX=0;AC_non_ukb_nfe_XY=0;AF_non_ukb_nfe_XY=0.00000;AN_non_ukb_nfe_XY=1470;nhomalt_non_ukb_nfe_XY=0;AC_non_ukb_nfe=0;AF_non_ukb_nfe=0.00000;AN_non_ukb_nfe=3260;nhomalt_non_ukb_nfe=0;AC_non_ukb_raw=6;AF_non_ukb_raw=9.54223e-06;AN_non_ukb_raw=628784;nhomalt_non_ukb_raw=1;AC_non_ukb_remaining_XX=0;AF_non_ukb_remaining_XX=0.00000;AN_non_ukb_remaining_XX=220;nhomalt_non_ukb_remaining_XX=0;AC_non_ukb_remaining_XY=0;AF_non_ukb_remaining_XY=0.00000;AN_non_ukb_remaining_XY=204;nhomalt_non_ukb_remaining_XY=0;AC_non_ukb_remaining=0;AF_non_ukb_remaining=0.00000;AN_non_ukb_remaining=424;nhomalt_non_ukb_remaining=0;AC_non_ukb_sas_XX=0;AF_non_ukb_sas_XX=0.00000;AN_non_ukb_sas_XX=320;nhomalt_non_ukb_sas_XX=0;AC_non_ukb_sas_XY=0;AF_non_ukb_sas_XY=0.00000;AN_non_ukb_sas_XY=1272;nhomalt_non_ukb_sas_XY=0;AC_non_ukb_sas=0;AF_non_ukb_sas=0.00000;AN_non_ukb_sas=1592;nhomalt_non_ukb_sas=0;AC_raw=6;AF_raw=9.54223e-06;AN_raw=628784;nhomalt_raw=1;AC_remaining_XX=0;AF_remaining_XX=0.00000;AN_remaining_XX=220;nhomalt_remaining_XX=0;AC_remaining_XY=0;AF_remaining_XY=0.00000;AN_remaining_XY=204;nhomalt_remaining_XY=0;AC_remaining=0;AF_remaining=0.00000;AN_remaining=424;nhomalt_remaining=0;AC_sas_XX=0;AF_sas_XX=0.00000;AN_sas_XX=320;nhomalt_sas_XX=0;AC_sas_XY=0;AF_sas_XY=0.00000;AN_sas_XY=1272;nhomalt_sas_XY=0;AC_sas=0;AF_sas=0.00000;AN_sas=1592;nhomalt_sas=0;AC_joint_XX=0;AF_joint_XX=0.00000;AN_joint_XX=3126;nhomalt_joint_XX=0;AC_joint_XY=0;AF_joint_XY=0.00000;AN_joint_XY=3642;nhomalt_joint_XY=0;AC_joint=0;AF_joint=0.00000;AN_joint=6768;nhomalt_joint=0;AC_joint_afr_XX=0;AF_joint_afr_XX=0.00000;AN_joint_afr_XX=196;nhomalt_joint_afr_XX=0;AC_joint_afr_XY=0;AF_joint_afr_XY=0.00000;AN_joint_afr_XY=162;nhomalt_joint_afr_XY=0;AC_joint_afr=0;AF_joint_afr=0.00000;AN_joint_afr=358;nhomalt_joint_afr=0;AC_joint_ami_XX=0;AN_joint_ami_XX=0;nhomalt_joint_ami_XX=0;AC_joint_ami_XY=0;AN_joint_ami_XY=0;nhomalt_joint_ami_XY=0;AC_joint_ami=0;AN_joint_ami=0;nhomalt_joint_ami=0;AC_joint_amr_XX=0;AF_joint_amr_XX=0.00000;AN_joint_amr_XX=166;nhomalt_joint_amr_XX=0;AC_joint_amr_XY=0;AF_joint_amr_XY=0.00000;AN_joint_amr_XY=108;nhomalt_joint_amr_XY=0;AC_joint_amr=0;AF_joint_amr=0.00000;AN_joint_amr=274;nhomalt_joint_amr=0;AC_joint_asj_XX=0;AF_joint_asj_XX=0.00000;AN_joint_asj_XX=48;nhomalt_joint_asj_XX=0;AC_joint_asj_XY=0;AF_joint_asj_XY=0.00000;AN_joint_asj_XY=58;nhomalt_joint_asj_XY=0;AC_joint_asj=0;AF_joint_asj=0.00000;AN_joint_asj=106;nhomalt_joint_asj=0;AC_joint_eas_XX=0;AF_joint_eas_XX=0.00000;AN_joint_eas_XX=266;nhomalt_joint_eas_XX=0;AC_joint_eas_XY=0;AF_joint_eas_XY=0.00000;AN_joint_eas_XY=236;nhomalt_joint_eas_XY=0;AC_joint_eas=0;AF_joint_eas=0.00000;AN_joint_eas=502;nhomalt_joint_eas=0;AC_joint_fin_XX=0;AF_joint_fin_XX=0.00000;AN_joint_fin_XX=112;nhomalt_joint_fin_XX=0;AC_joint_fin_XY=0;AF_joint_fin_XY=0.00000;AN_joint_fin_XY=116;nhomalt_joint_fin_XY=0;AC_joint_fin=0;AF_joint_fin=0.00000;AN_joint_fin=228;nhomalt_joint_fin=0;AC_joint_mid_XX=0;AF_joint_mid_XX=0.00000;AN_joint_mid_XX=8;nhomalt_joint_mid_XX=0;AC_joint_mid_XY=0;AF_joint_mid_XY=0.00000;AN_joint_mid_XY=16;nhomalt_joint_mid_XY=0;AC_joint_mid=0;AF_joint_mid=0.00000;AN_joint_mid=24;nhomalt_joint_mid=0;AC_joint_nfe_XX=0;AF_joint_nfe_XX=0.00000;AN_joint_nfe_XX=1790;nhomalt_joint_nfe_XX=0;AC_joint_nfe_XY=0;AF_joint_nfe_XY=0.00000;AN_joint_nfe_XY=1470;nhomalt_joint_nfe_XY=0;AC_joint_nfe=0;AF_joint_nfe=0.00000;AN_joint_nfe=3260;nhomalt_joint_nfe=0;AC_joint_raw=6;AF_joint_raw=9.54223e-06;AN_joint_raw=628784;nhomalt_joint_raw=1;AC_joint_remaining_XX=0;AF_joint_remaining_XX=0.00000;AN_joint_remaining_XX=220;nhomalt_joint_remaining_XX=0;AC_joint_remaining_XY=0;AF_joint_remaining_XY=0.00000;AN_joint_remaining_XY=204;nhomalt_joint_remaining_XY=0;AC_joint_remaining=0;AF_joint_remaining=0.00000;AN_joint_remaining=424;nhomalt_joint_remaining=0;AC_joint_sas_XX=0;AF_joint_sas_XX=0.00000;AN_joint_sas_XX=320;nhomalt_joint_sas_XX=0;AC_joint_sas_XY=0;AF_joint_sas_XY=0.00000;AN_joint_sas_XY=1272;nhomalt_joint_sas_XY=0;AC_joint_sas=0;AF_joint_sas=0.00000;AN_joint_sas=1592;nhomalt_joint_sas=0;nhomalt_grpmax=2147483647;nhomalt_grpmax_non_ukb=2147483647;faf95=0.00000;faf95_afr=0.00000;faf95_amr=0.00000;faf95_eas=0.00000;faf95_nfe=0.00000;faf95_non_ukb=0.00000;faf95_non_ukb_afr=0.00000;faf95_non_ukb_amr=0.00000;faf95_non_ukb_eas=0.00000;faf95_non_ukb_nfe=0.00000;faf95_non_ukb_sas=0.00000;faf95_sas=0.00000;faf99=0.00000;faf99_afr=0.00000;faf99_amr=0.00000;faf99_eas=0.00000;faf99_nfe=0.00000;faf99_non_ukb=0.00000;faf99_non_ukb_afr=0.00000;faf99_non_ukb_amr=0.00000;faf99_non_ukb_eas=0.00000;faf99_non_ukb_nfe=0.00000;faf99_non_ukb_sas=0.00000;faf99_sas=0.00000;faf95_joint=0.00000;faf95_joint_afr=0.00000;faf95_joint_amr=0.00000;faf95_joint_eas=0.00000;faf95_joint_nfe=0.00000;faf95_joint_sas=0.00000;faf99_joint=0.00000;faf99_joint_afr=0.00000;faf99_joint_amr=0.00000;faf99_joint_eas=0.00000;faf99_joint_nfe=0.00000;faf99_joint_sas=0.00000;age_hist_het_bin_freq=0|0|0|0|0|0|0|0|0|0;age_hist_het_n_smaller=0;age_hist_het_n_larger=0;age_hist_hom_bin_freq=0|0|0|0|0|0|0|0|0|0;age_hist_hom_n_smaller=0;age_hist_hom_n_larger=0;MQ=29.2923;MQRankSum=0.00000;QUALapprox=370;QD=13.7037;ReadPosRankSum=0.967000;SOR=0.764518;VarDP=27;AS_MQ=29.2923;AS_MQRankSum=0.00000;AS_pab_max=1.00000;AS_QUALapprox=370;AS_QD=13.7037;AS_ReadPosRankSum=0.967000;AS_SB_TABLE=13,0|14,0;AS_SOR=0.764518;AS_VarDP=27;inbreeding_coeff=0.333327;AS_culprit=AS_MQ;AS_VQSLOD=-1.75460;allele_type=snv;n_alt_alleles=1;variant_type=snv;segdup;outside_ukb_capture_region;gq_hist_alt_bin_freq=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;gq_hist_all_bin_freq=0|0|0|0|723|0|2291|0|370|0|0|0|0|0|0|0|0|0|0|0;dp_hist_alt_bin_freq=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;dp_hist_alt_n_larger=0;dp_hist_all_bin_freq=0|0|3014|304|51|13|2|0|0|0|0|0|0|0|0|0|0|0|0|0;dp_hist_all_n_larger=0;ab_hist_alt_bin_freq=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;cadd_raw_score=1.54599;cadd_phred=16.0300;pangolin_largest_ds=-0.0700000;phylop=0.433000;sift_max=0.930000;polyphen_max=0.00000;VRS_Allele_IDs=ga4gh:VA.RR0YR51puLix-IUintOOBYT0yaB92fmZ,ga4gh:VA.yRO4oXggC6MKMgq5ZSaj8s4yUxP5de_F;VRS_Starts=65570,65570;VRS_Ends=65571,65571;VRS_States=A,C;vep=C|downstream_gene_variant|MODIFIER|OR4G11P|ENSG00000240361|Transcript|ENST00000492842|transcribed_unprocessed_pseudogene||||||||||1|1684|1||SNV|HGNC|HGNC:31276|YES||||||||Ensembl|||||||||||||||,C|missense_variant&splice_region_variant|MODERATE|OR4F5|ENSG00000186092|Transcript|ENST00000641515|protein_coding|2/3||ENST00000641515.2:c.7A>C|ENSP00000493376.2:p.Lys3Gln|67|7|3|K/Q|Aag/Cag|1||1||SNV|HGNC|HGNC:14825|YES|NM_001005484.2|||P1||ENSP00000493376||Ensembl|||||||||||||||,C|downstream_gene_variant|MODIFIER|OR4G11P|ENSG00000240361|Transcript|ENST00000642116|processed_transcript||||||||||1|1455|1||SNV|HGNC|HGNC:31276|||||||||Ensembl|||||||||||||||,C|missense_variant&splice_region_variant|MODERATE|OR4F5|79501|Transcript|NM_001005484.2|protein_coding|2/3||NM_001005484.2:c.7A>C|NP_001005484.2:p.Lys3Gln|67|7|3|K/Q|Aag/Cag|1||1||SNV|EntrezGene|HGNC:14825|YES|ENST00000641515.2|||||NP_001005484.2||RefSeq|||||||||||||||



            # v4 genomes fields (for a missense variant):
            # zcat gnomad.genomes.v4.0.sites.chr1.vcf.bgz | head -n 100000 | g -v "^#" | g missense | head -n1
            # chr1	69037	rs1639970567	G	A	.	AC0;AS_VQSR	AC=0;AN=11022;AF=0.00000;AC_XX=0;AF_XX=0.00000;AN_XX=6170;nhomalt_XX=0;AC_XY=0;AF_XY=0.00000;AN_XY=4852;nhomalt_XY=0;nhomalt=0;AC_afr_XX=0;AF_afr_XX=0.00000;AN_afr_XX=3560;nhomalt_afr_XX=0;AC_afr_XY=0;AF_afr_XY=0.00000;AN_afr_XY=2714;nhomalt_afr_XY=0;AC_afr=0;AF_afr=0.00000;AN_afr=6274;nhomalt_afr=0;AC_ami_XX=0;AF_ami_XX=0.00000;AN_ami_XX=12;nhomalt_ami_XX=0;AC_ami_XY=0;AF_ami_XY=0.00000;AN_ami_XY=8;nhomalt_ami_XY=0;AC_ami=0;AF_ami=0.00000;AN_ami=20;nhomalt_ami=0;AC_amr_XX=0;AF_amr_XX=0.00000;AN_amr_XX=360;nhomalt_amr_XX=0;AC_amr_XY=0;AF_amr_XY=0.00000;AN_amr_XY=390;nhomalt_amr_XY=0;AC_amr=0;AF_amr=0.00000;AN_amr=750;nhomalt_amr=0;AC_asj_XX=0;AF_asj_XX=0.00000;AN_asj_XX=42;nhomalt_asj_XX=0;AC_asj_XY=0;AF_asj_XY=0.00000;AN_asj_XY=28;nhomalt_asj_XY=0;AC_asj=0;AF_asj=0.00000;AN_asj=70;nhomalt_asj=0;AC_eas_XX=0;AF_eas_XX=0.00000;AN_eas_XX=440;nhomalt_eas_XX=0;AC_eas_XY=0;AF_eas_XY=0.00000;AN_eas_XY=534;nhomalt_eas_XY=0;AC_eas=0;AF_eas=0.00000;AN_eas=974;nhomalt_eas=0;AC_fin_XX=0;AF_fin_XX=0.00000;AN_fin_XX=26;nhomalt_fin_XX=0;AC_fin_XY=0;AF_fin_XY=0.00000;AN_fin_XY=34;nhomalt_fin_XY=0;AC_fin=0;AF_fin=0.00000;AN_fin=60;nhomalt_fin=0;AC_mid_XX=0;AF_mid_XX=0.00000;AN_mid_XX=20;nhomalt_mid_XX=0;AC_mid_XY=0;AF_mid_XY=0.00000;AN_mid_XY=24;nhomalt_mid_XY=0;AC_mid=0;AF_mid=0.00000;AN_mid=44;nhomalt_mid=0;AC_nfe_XX=0;AF_nfe_XX=0.00000;AN_nfe_XX=1456;nhomalt_nfe_XX=0;AC_nfe_XY=0;AF_nfe_XY=0.00000;AN_nfe_XY=768;nhomalt_nfe_XY=0;AC_nfe=0;AF_nfe=0.00000;AN_nfe=2224;nhomalt_nfe=0;AC_raw=1;AF_raw=1.08826e-05;AN_raw=91890;nhomalt_raw=0;AC_remaining_XX=0;AF_remaining_XX=0.00000;AN_remaining_XX=76;nhomalt_remaining_XX=0;AC_remaining_XY=0;AF_remaining_XY=0.00000;AN_remaining_XY=78;nhomalt_remaining_XY=0;AC_remaining=0;AF_remaining=0.00000;AN_remaining=154;nhomalt_remaining=0;AC_sas_XX=0;AF_sas_XX=0.00000;AN_sas_XX=178;nhomalt_sas_XX=0;AC_sas_XY=0;AF_sas_XY=0.00000;AN_sas_XY=274;nhomalt_sas_XY=0;AC_sas=0;AF_sas=0.00000;AN_sas=452;nhomalt_sas=0;AC_joint_XX=0;AF_joint_XX=0.00000;AN_joint_XX=6170;nhomalt_joint_XX=0;AC_joint_XY=0;AF_joint_XY=0.00000;AN_joint_XY=4852;nhomalt_joint_XY=0;AC_joint=0;AF_joint=0.00000;AN_joint=11022;nhomalt_joint=0;AC_joint_afr_XX=0;AF_joint_afr_XX=0.00000;AN_joint_afr_XX=3560;nhomalt_joint_afr_XX=0;AC_joint_afr_XY=0;AF_joint_afr_XY=0.00000;AN_joint_afr_XY=2714;nhomalt_joint_afr_XY=0;AC_joint_afr=0;AF_joint_afr=0.00000;AN_joint_afr=6274;nhomalt_joint_afr=0;AC_joint_ami_XX=0;AF_joint_ami_XX=0.00000;AN_joint_ami_XX=12;nhomalt_joint_ami_XX=0;AC_joint_ami_XY=0;AF_joint_ami_XY=0.00000;AN_joint_ami_XY=8;nhomalt_joint_ami_XY=0;AC_joint_ami=0;AF_joint_ami=0.00000;AN_joint_ami=20;nhomalt_joint_ami=0;AC_joint_amr_XX=0;AF_joint_amr_XX=0.00000;AN_joint_amr_XX=360;nhomalt_joint_amr_XX=0;AC_joint_amr_XY=0;AF_joint_amr_XY=0.00000;AN_joint_amr_XY=390;nhomalt_joint_amr_XY=0;AC_joint_amr=0;AF_joint_amr=0.00000;AN_joint_amr=750;nhomalt_joint_amr=0;AC_joint_asj_XX=0;AF_joint_asj_XX=0.00000;AN_joint_asj_XX=42;nhomalt_joint_asj_XX=0;AC_joint_asj_XY=0;AF_joint_asj_XY=0.00000;AN_joint_asj_XY=28;nhomalt_joint_asj_XY=0;AC_joint_asj=0;AF_joint_asj=0.00000;AN_joint_asj=70;nhomalt_joint_asj=0;AC_joint_eas_XX=0;AF_joint_eas_XX=0.00000;AN_joint_eas_XX=440;nhomalt_joint_eas_XX=0;AC_joint_eas_XY=0;AF_joint_eas_XY=0.00000;AN_joint_eas_XY=534;nhomalt_joint_eas_XY=0;AC_joint_eas=0;AF_joint_eas=0.00000;AN_joint_eas=974;nhomalt_joint_eas=0;AC_joint_fin_XX=0;AF_joint_fin_XX=0.00000;AN_joint_fin_XX=26;nhomalt_joint_fin_XX=0;AC_joint_fin_XY=0;AF_joint_fin_XY=0.00000;AN_joint_fin_XY=34;nhomalt_joint_fin_XY=0;AC_joint_fin=0;AF_joint_fin=0.00000;AN_joint_fin=60;nhomalt_joint_fin=0;AC_joint_mid_XX=0;AF_joint_mid_XX=0.00000;AN_joint_mid_XX=20;nhomalt_joint_mid_XX=0;AC_joint_mid_XY=0;AF_joint_mid_XY=0.00000;AN_joint_mid_XY=24;nhomalt_joint_mid_XY=0;AC_joint_mid=0;AF_joint_mid=0.00000;AN_joint_mid=44;nhomalt_joint_mid=0;AC_joint_nfe_XX=0;AF_joint_nfe_XX=0.00000;AN_joint_nfe_XX=1456;nhomalt_joint_nfe_XX=0;AC_joint_nfe_XY=0;AF_joint_nfe_XY=0.00000;AN_joint_nfe_XY=768;nhomalt_joint_nfe_XY=0;AC_joint_nfe=0;AF_joint_nfe=0.00000;AN_joint_nfe=2224;nhomalt_joint_nfe=0;AC_joint_raw=1;AF_joint_raw=1.08826e-05;AN_joint_raw=91890;nhomalt_joint_raw=0;AC_joint_remaining_XX=0;AF_joint_remaining_XX=0.00000;AN_joint_remaining_XX=76;nhomalt_joint_remaining_XX=0;AC_joint_remaining_XY=0;AF_joint_remaining_XY=0.00000;AN_joint_remaining_XY=78;nhomalt_joint_remaining_XY=0;AC_joint_remaining=0;AF_joint_remaining=0.00000;AN_joint_remaining=154;nhomalt_joint_remaining=0;AC_joint_sas_XX=0;AF_joint_sas_XX=0.00000;AN_joint_sas_XX=178;nhomalt_joint_sas_XX=0;AC_joint_sas_XY=0;AF_joint_sas_XY=0.00000;AN_joint_sas_XY=274;nhomalt_joint_sas_XY=0;AC_joint_sas=0;AF_joint_sas=0.00000;AN_joint_sas=452;nhomalt_joint_sas=0;faf95=0.00000;faf95_afr=0.00000;faf95_amr=0.00000;faf95_eas=0.00000;faf95_nfe=0.00000;faf95_sas=0.00000;faf99=0.00000;faf99_afr=0.00000;faf99_amr=0.00000;faf99_eas=0.00000;faf99_nfe=0.00000;faf99_sas=0.00000;faf95_joint=0.00000;faf95_joint_afr=0.00000;faf95_joint_amr=0.00000;faf95_joint_eas=0.00000;faf95_joint_nfe=0.00000;faf95_joint_sas=0.00000;faf99_joint=0.00000;faf99_joint_afr=0.00000;faf99_joint_amr=0.00000;faf99_joint_eas=0.00000;faf99_joint_nfe=0.00000;faf99_joint_sas=0.00000;age_hist_het_bin_freq=0|0|0|0|0|0|0|0|0|0;age_hist_het_n_smaller=0;age_hist_het_n_larger=0;age_hist_hom_bin_freq=0|0|0|0|0|0|0|0|0|0;age_hist_hom_n_smaller=0;age_hist_hom_n_larger=0;FS=0.00000;MQ=25.0000;MQRankSum=0.988000;QUALapprox=71;QD=11.8333;ReadPosRankSum=0.406000;SOR=0.132000;VarDP=6;AS_FS=0.00000;AS_MQ=25.0000;AS_MQRankSum=0.988000;AS_pab_max=1.00000;AS_QUALapprox=71;AS_QD=11.8333;AS_ReadPosRankSum=0.406000;AS_SB_TABLE=0,3|1,2;AS_SOR=0.131576;AS_VarDP=6;inbreeding_coeff=-1.08827e-05;AS_culprit=AS_MQ;AS_VQSLOD=-5.09250;allele_type=snv;n_alt_alleles=1;variant_type=snv;segdup;gq_hist_alt_bin_freq=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;gq_hist_all_bin_freq=0|0|0|0|4405|672|346|71|11|4|2|0|0|0|0|0|0|0|0|0;dp_hist_alt_bin_freq=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;dp_hist_alt_n_larger=0;dp_hist_all_bin_freq=0|0|3472|1575|301|131|30|2|0|0|0|0|0|0|0|0|0|0|0|0;dp_hist_all_n_larger=0;ab_hist_alt_bin_freq=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;cadd_raw_score=1.08152;cadd_phred=12.5800;pangolin_largest_ds=-0.180000;phylop=4.17300;sift_max=0.280000;polyphen_max=0.0580000;VRS_Allele_IDs=ga4gh:VA.1_1P_5X3Rc-sgXDfp-Irohwly7w35AcI,ga4gh:VA.HHpBCMkQw6mGBbl4P5ox-BrBt_I_SDKs;VRS_Starts=69036,69036;VRS_Ends=69037,69037;VRS_States=G,A;vep=A|missense_variant&splice_region_variant|MODERATE|OR4F5|ENSG00000186092|Transcript|ENST00000641515|protein_coding|3/3||ENST00000641515.2:c.10G>A|ENSP00000493376.2:p.Val4Ile|70|10|4|V/I|Gta/Ata|1||1||SNV|HGNC|HGNC:14825|YES|NM_001005484.2|||P1||ENSP00000493376||Ensembl|||||||||||||||,A|downstream_gene_variant|MODIFIER|OR4G11P|ENSG00000240361|Transcript|ENST00000642116|processed_transcript||||||||||1|4921|1||SNV|HGNC|HGNC:31276|||||||||Ensembl|||||||||||||||,A|missense_variant&splice_region_variant|MODERATE|OR4F5|79501|Transcript|NM_001005484.2|protein_coding|3/3||NM_001005484.2:c.10G>A|NP_001005484.2:p.Val4Ile|70|10|4|V/I|Gta/Ata|1||1||SNV|EntrezGene|HGNC:14825|YES|ENST00000641515.2|||||NP_001005484.2||RefSeq|||||||||||||||



            # Split line
            e = line.rstrip().split("\t")

            # Skip any variants that haven't passed all filters
            if e[6] != "PASS":
                continue
            
            # Parse tags
            csq = None

            tags = e[7].split(";")
            fields = {}
            for tag in tags:
                # CSQ (Ensembl VEP consequences)
                if tag.startswith("CSQ="):
                    csq = tag[4:]
                # Flags (no value)
                elif '=' not in tag:
                    fields[tag] = 1
                # Fields
                else:
                    m = rx(r'^([^=]+)=(.+?)$', tag)
                    if not m:
                        Die(f"Error: Couldn't match tag:\n\n{tag}\n\n")
                    (field, value) = m
                    fields[field] = value
            
            # VEP's CSQ string (only exists if there is any coding consequence for the variant)
            if csq is None:
                # Log("no consequence for chr|pos|originalbase|variantbase (skipped)", f"{chr}|{pos}|{originalbase}|{variantbase}")
                continue
            elif csq.find("missense_variant") == -1:
                # Log("no missense_variant in csq for chr|pos|originalbase|variantbase (skipped)", f"{chr}|{pos}|{originalbase}|{variantbase}")
                continue

            # Select fields of interest
            newfields = {}
            for field in fields:
                if (
                    # Allele counts
                    rx(r'^AC_?(female|male|XX|XY|popmax|grpmax)?$', field) or 
                    rx(r'^AC_?(grpmax)?(_joint)?(_XX|_XY)?$', field) or # for v4 joint
                    # Allele frequencies
                    rx(r'^AF_?(female|male|XX|XY|popmax|grpmax)?$', field) or 
                    rx(r'^AF_?(controls_and_biobanks|non_cancer|non_neuro|non_topmed|non_v2)?$', field) or 
                    rx(r'^(controls|non_cancer|non_neuro|non_topmed)_AF$', field) or # v2
                    rx(r'^AF_(afr|ami|amr|asj|eas|fin|mid|nfe|sas|oth|remaining)$', field) or 
                    rx(r'^AF_?(grpmax)?(_joint)?(_(XX|XY|afr|ami|amr|asj|eas|fin|mid|nfe|sas|oth|remaining))?$', field) or # for v4 joint
                    # Filtering allele frequency (FAF), 95% confidence, in all non-bottlenecked groups (to determine the one with the highest FAF, and also compare to the overall faf95)
                    rx(r'^faf95_?(afr|amr|eas|nfe|sas)?$', field) or # for v2 and v3
                    rx(r'^faf95_joint_?(afr|amr|eas|nfe|sas)?$', field) or # for v4 joint
                    rx(r'^fafmax_faf95_max_joint$', field) or # for v4 joint (not available if faf95_joint_(afr|amr|eas|nfe|sas) are all 0)
                    # Group with maximum allele frequency
                    rx(r'^(popmax|grpmax)$', field) or 
                    rx(r'^(popmax|grpmax)(_joint)?$', field) or # for v4 joint
                    # Allele numbers (total)
                    rx(r'^AN_?(female|male|XX|XY|popmax|grpmax)?$', field) or 
                    rx(r'^AN_?(grpmax)?(_joint)?(_XX|_XY)?$', field) or # for v4 joint
                    # Strand bias
                    rx(r'^(AS_)?(FS|pab_max|SOR)$', field) or 
                    # Homo-/heterozygosity, inbreeding coefficient etc.
                    rx(r'^(monoallelic|n_alt_alleles|transmitted_singleton|InbreedingCoeff|inbreeding_coeff)$', field) or 
                    # Number of homozygous individuals
                    rx(r'^nhomalt_?(female|male|XX|XY|popmax|grpmax)?$', field) or 
                    rx(r'^nhomalt_?(grpmax)?(_joint)?(_XX|_XY)?$', field) or # for v4 joint
                    # Region characteristics: low-complexity, non-pseudoautosomal region of a sex chromosome, in a segmental duplication etc.
                    rx(r'^(lcr|nonpar|non_par|segdup)$', field) or 
                    # Variant type
                    rx(r'^(allele_type|variant_type|was_mixed)$', field)
                    ):

                    # Rename some fields to make v2 match v3
                    newfield = field
                    newfield = re.sub(r'_female', r'_XX', newfield)
                    newfield = re.sub(r'_male', r'_XY', newfield)
                    newfield = re.sub(r'^controls_AF$', r'AF_controls_and_biobanks', newfield)
                    newfield = re.sub(r'^(\w+)_AF$', r'AF_\1', newfield)    # e.g. non_neuro_AF (v2) -> AF_non_neuro (v3)
                    # newfield = re.sub(r'^FS$', r'AS_FS', newfield)        # in v3 and v4, these are actually distinct (keeping both)
                    # newfield = re.sub(r'^SOR$', r'AS_SOR', newfield)      # in v4, these are actually distinct (keeping both)
                    newfield = re.sub(r'^pab_max$', r'AS_pab_max', newfield)

                    # Rename some fields to make v2 and v3 match v4
                    # ac_popmax
                    # an_popmax
                    # af_popmax
                    # popmax
                    # nhomalt_popmax
                    newfield = re.sub(r'popmax', r'grpmax', newfield)
                    # InbreedingCoeff
                    newfield = re.sub(r'InbreedingCoeff', r'inbreeding_coeff', newfield)
                    # af_controls_and_biobanks
                    # af_non_cancer
                    # af_non_neuro
                    # af_non_topmed
                    # af_non_v2
                    # >> missing from v4
                    # af_oth
                    newfield = re.sub(r'AF_oth', r'AF_remaining', newfield)
                    # nonpar
                    newfield = re.sub(r'^nonpar$', r'^non_par$', newfield)
                    
                    # v4 "joint" fields:
                    # apparently "called jointly across the gnomAD v4 exome and genome datasets" (https://gnomad.broadinstitute.org/news/2023-11-ga4gh-gks/)
                    # The AN numbers match up for this as well:
                    # Estimated on chr1 (VEP, which explains the difference to 1,614,320 above (yes, 1,614,320 is the result with chr*)):
                    # zcat input/gnomad.exomes.v4.0.sites.chr1.vep.vcf.bgz | g -o ";AN_joint=\d+" | g -o "\d+" | max
                    # >> 1,614,310 / 2 = 807,155 diploid individuals
                    # zcat input/gnomad.exomes.v4.0.sites.chr1.vcf.bgz | g -o ";AN_joint=\d+" | g -o "\d+" | max
                    # >> 1,614,310 / 2 = 807,155 diploid individuals
                    # Should be 730,947 exomes + 76,215 genomes = 807,162 total (https://gnomad.broadinstitute.org/downloads#v4)
                    # >> Almost matches - presumably looking across all chromosomes might yield a higher AN_joint value that matches perfectly? See above for this (chr*).
                    # zcat input/gnomad.exomes.v4.0.sites.chr*.vcf.bgz | g -o ";AN_joint=\d+" | g -o "\d+" | max
                    # >> 1,614,320
                    # >> Yes, almost perfect match. Should be 1,614,324 (4 more), but that's close enough.


                    # Also matches here:
                    # zcat gnomad.genomes.v4.0.sites.chrY.vep.vcf.bgz | g -o ";AN_joint=\d+" | g -o "\d+" | max # >> 398,004
                    # zcat gnomad.genomes.v4.0.sites.chrY.vep.vcf.bgz | g -o ";AN=\d+" | g -o "\d+" | max       # >> 37,269
                    # zcat gnomad.exomes.v4.0.sites.chrY.vep.vcf.bgz | g -o ";AN_joint=\d+" | g -o "\d+" | max  # >> 398,004
                    # zcat gnomad.exomes.v4.0.sites.chrY.vep.vcf.bgz | g -o ";AN=\d+" | g -o "\d+" | max        # >> 363,593
                    # >> genomes and exomes chrY AN_joint maximum is identical

                    # >> Since v4 is by far the biggest datasets and includes the previous versions (I think all), I should simply use the joint fields.
                    # >> These will give me the most accurate allele frequencies etc.

                    # Rename v4 "joint" fields (26) to be used instead of the regular fields
                    # AC_grpmax_joint
                    # AC_joint
                    # AC_joint_XX
                    # AC_joint_XY
                    # AF_grpmax_joint
                    # AF_joint
                    # AF_joint_afr
                    # AF_joint_ami
                    # AF_joint_amr
                    # AF_joint_asj
                    # AF_joint_eas
                    # AF_joint_fin
                    # AF_joint_mid
                    # AF_joint_nfe
                    # AF_joint_sas
                    # AF_joint_XX
                    # AF_joint_XY
                    # AN_grpmax_joint
                    # AN_joint
                    # AN_joint_XX
                    # AN_joint_XY
                    # grpmax_joint
                    # nhomalt_grpmax_joint
                    # nhomalt_joint
                    # nhomalt_joint_XX
                    # nhomalt_joint_XY
                    if source.startswith('v4_'):
                        # Skip non-joint version of these fields
                        if (
                            rx(r'^AC_?(female|male|XX|XY|popmax|grpmax)?$', field) or 
                            # Allele frequencies
                            rx(r'^AF_?(female|male|XX|XY|popmax|grpmax)?$', field) or 
                            rx(r'^AF_?(controls_and_biobanks|non_cancer|non_neuro|non_topmed|non_v2)?$', field) or 
                            rx(r'^AF_(afr|ami|amr|asj|eas|fin|mid|nfe|sas|oth|remaining)$', field) or 
                            # Filtering allele frequency (FAF), 95% confidence, in all non-bottlenecked groups (to determine the one with the highest FAF, and also compare to the overall faf95)
                            rx(r'^faf95_?(afr|amr|eas|nfe|sas)?$', field) or # for v2 and v3
                            # Group with maximum allele frequency
                            rx(r'^(popmax|grpmax)$', field) or 
                            # Allele numbers (total)
                            rx(r'^AN_?(female|male|XX|XY|popmax|grpmax)?$', field) or
                            # Number of homozygous individuals
                            rx(r'^nhomalt_?(female|male|XX|XY|popmax|grpmax)?$', field)
                        ):
                            continue
                        # Use "joint" field as newfield
                        newfield = re.sub(r'_joint', r'', newfield)
                        
                    # If renamed: replace old field with new
                    if newfield != field:
                        Log("field renamed for oldfield|newfield", f"{field}|{newfield}")
                        Log("field renamed for oldfield", field)
                        Log("field renamed for newfield", newfield)

                    # Also make all fields lower case (but do not log these renames)
                    newfield = newfield.lower()

                    # Add to newfields dict
                    if newfield in newfields:
                        Die(f"Error: newfield '{newfield}' was already defined as '{field}' in newfields '{newfields}'\n\nline: {line}\n\n")
                    newfields[newfield] = fields[field]
                    Log("field kept for field", field)

                else:

                    Log("field ignored for field", field)
            
            # Verify field content
            # popmax can occasionally be missing, e.g. for CSQ=C|ENSG00000184895|Transcript|ENST00000383070|missense_variant||119|M/R| (which has AC=2 and AC_ami=2, but every other subpopulation AC_...=0)
            # if "popmax" not in newfields:
            #     Die(f"Error: No 'popmax' field in tags:\n\n{tags}\n\n")
            # if newfields["popmax"] not in ("afr", "ami", "amr", "asj", "eas", "fin", "mid", "nfe", "oth", "sas"):
            #     Die(f"Error: popmax is '{newfields['popmax']}'")

            # fafmax95: Process faf95 fields (not included) to populate a single field (included)
            # v2_wes: Use max(faf95, faf95_(afr|amr|eas|nfe|sas)), but set faf95_… to 0 first for group singletons (AC_…=1) to match what's in the single file (https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/README.txt). Don't just use faf95 since that's not what gets used on the website. e.g. faf95=1.77680e-01;faf95_afr=0.00000e+00;faf95_amr=1.77680e-01;faf95_eas=0.00000e+00;faf95_nfe=0.00000e+00;faf95_sas=0.00000e+00 (matches https://gnomad.broadinstitute.org/variant/8-7439430-A-G?dataset=gnomad_r2_1) (faf95 seems to be the max of the other faf95_…s, but it isn't - see below) (also, faf95_max_gen_anc doesn't exist in v2)
            # v2_wgs: Use max(faf95, faf95_(afr|amr|eas|nfe)), but set faf95_… to 0 first for group singletons (AC_…=1) to match what's in the single file (https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/README.txt). Don't just use faf95 since that's not what gets used on the website. e.g. faf95=4.08837e-02;faf95_afr=3.20500e-03;faf95_amr=0.00000e+00;faf95_eas=0.00000e+00;faf95_nfe=8.88401e-02; (matches https://gnomad.broadinstitute.org/variant/8-7439430-A-G?dataset=gnomad_r2_1) (faf95 seems to be the max of the other faf95_…s, but it isn't - see below) (also, faf95_max_gen_anc doesn't exist in v2)
            # v3_wgs: Use max(faf95, faf95_(afr|amr|eas|nfe|sas)). Singleton MAF is already always 0, so I don't need to set it to 0 myself (OK).
            # v4_wes: Use fafmax_faf95_max_joint. Peculiarity: fafmax_faf95_max_joint doesn't exist if AC_joint=1 (singletons), should consider it 0 in those cases and use faf95_joint instead. e.g. fafmax_faf95_max=1.00200e-05;fafmax_faf95_max_gen_anc=afr (matches https://gnomad.broadinstitute.org/variant/22-15528195-A-G?dataset=gnomad_r4)
            # v4_wgs: Use fafmax_faf95_max_joint. Peculiarity: fafmax_faf95_max_joint doesn't exist if AC_joint=1 (singletons), should consider it 0 in those cases and use faf95_joint instead. e.g. fafmax_faf95_max=4.74800e-05;fafmax_faf95_max_gen_anc=afr (matches https://gnomad.broadinstitute.org/variant/22-15528195-A-G?dataset=gnomad_r4)
            # v2
            # v3
            # v4: fafmax_faf95_max_joint
            
            # Set groups to use for FAF95 (non-bottlenecked groups)
            if source == 'v2_wgs':
                # no sas
                groups = ('afr', 'amr', 'eas', 'nfe')
            else:
                groups = ('afr', 'amr', 'eas', 'nfe', 'sas')

            # set faf95_… to 0 first for group singletons (AC_…=1) to match what's in the single file (which, unlike the individual chromosome files for v2, has been updated)
            if source.startswith('v2_'):
                for group in groups:
                    if fields[f'AC_{group}'] == 1:
                        fields[f'faf95_{group}'] = 0
            
            # Calculate fafmax95 (maximum across non-bottlenecked FAF95s, and the global FAF95)
            if source.startswith('v4_'):
                if field == 'fafmax_faf95_max_joint' in fields:
                    newfields['fafmax95'] = max(float(fields['fafmax_faf95_max_joint']), float(fields['faf95_joint']))
                else:
                    # fafmax_faf95_max_joint is only ever missing if all group FAF95s are 0, in which case I'll need to use faf95_joint
                    newfields['fafmax95'] = fields['faf95_joint']
                    newfields['fafmax_faf95_max'] = max(float(fields[f'faf95_{group}']) for group in groups)    # Will be 0 (otherwise fafmax_faf95_max wouldn't be missing, i.e. NULL - this way it becomes 0)
            else:
                newfields['fafmax95'] = max(float(fields['faf95']), max(float(fields[f'faf95_{group}']) for group in groups))
            Log("field added for field", "fafmax95")

            # Calculate fafmax_faf95_max for v2 and v3
            if source.startswith('v2_') or source.startswith('v3_'):
                newfields['fafmax_faf95_max'] = max(float(fields[f'faf95_{group}']) for group in groups)
                Log("field added for field", "fafmax_faf95_max")


            # Construct string using newfields
            # List of fields of interest
            fieldstring = ""
            for field in myfields:
                if field in newfields:
                    # Append field
                    fieldstring += "\t" + str(newfields[field])
                else:
                    # Not defined: print \N (for NULL)
                    fieldstring += "\t\\N"


                    

            # Check if variants are larger than a single nucleotide (just for logging)
            
            # See snps_gnomad.sql for more details
            # zcat input/gnomad.exomes.v4.0.sites.chr1.vcf.bgz | perl -ne '@a = split(/\t/); if (((length($a[3]) == 1) and (length($a[4]) == 1)) and ($_ =~ /missense_variant/)) { print }' | wcl
            # >> The vast majority are true SNVs (only affecting a single base)
            # zcat input/gnomad.exomes.v4.0.sites.chr1.vcf.bgz | perl -ne '@a = split(/\t/); if (((length($a[3]) > 1) or (length($a[4]) > 1)) and ($_ =~ /missense_variant/)) { print }' | wcl
            # >> Very few variants are larger. In principle, since I don't care about the mutational mechanism of the variant, but only about the amino acid change, I'll keep these larger variants (but only if they only affect a single amino acid).
            # zcat gnomad.exomes.v4.0.sites.chr1.vcf.bgz | g -v ^# | cut -f4,5 | perl -ne 'chomp; @a = split(/\t/); if ((length($a[0]) == length($a[1])) and (length($a[0]) > 1)) { print "EQUAL LENGTH AND LENGTH > 1\t$_\n" }'
            # >> However, in gnomAD, whenever the original and variant are larger than a single base, their lengths are unequal (so they'll lead to frameshifts rather than missense). So there won't be any larger variants that only affect a single amino acid. See snps_gnomad.sql for confirmation.
            
            # Also, there's always only a single alternative allele in each line, so I don't need to handle these:
            # zcat input/gnomad.exomes.v4.0.sites.chr1.vcf.bgz | g -v "^#" | cut -f4,5 | g ','
            # >> None
            chr = e[0]
            pos = e[1]
            rsids = e[2]
            originalbase = e[3]
            variantbase = e[4]
            if len(originalbase) == 1 and len(variantbase) == 1:
                Log("original and variant are a single nucleotide (i.e. a true SNV) for chr|pos|originalbase|variantbase (kept)", f"{chr}|{pos}|{originalbase}|{variantbase}")
                Log("original and variant are a single nucleotide (i.e. a true SNV) for originalbase|variantbase (kept)", f"{originalbase}|{variantbase}")
            if len(originalbase) != len(variantbase):
                Log("original and variant lengths differ (must be a frameshift) for chr|pos|originalbase|variantbase (kept)", f"{chr}|{pos}|{originalbase}|{variantbase}")
                Log("original and variant lengths differ (must be a frameshift) for originalbase|variantbase (kept)", f"{originalbase}|{variantbase}")
            if len(originalbase) > 1:
                Log("original larger than a single nucleotide for chr|pos|originalbase|variantbase (kept)", f"{chr}|{pos}|{originalbase}|{variantbase}")
                Log("original larger than a single nucleotide for originalbase|variantbase (kept)", f"{originalbase}|{variantbase}")
                # continue
            if len(variantbase) > 1:
                Log("variant larger than a single nucleotide for chr|pos|originalbase|variantbase (kept)", f"{chr}|{pos}|{originalbase}|{variantbase}")
                Log("variant larger than a single nucleotide for originalbase|variantbase (kept)", f"{originalbase}|{variantbase}")
                # continue



            # Filter out variants where AN is below 50% of the possible maximum
            # "This variant is covered in fewer than 50% of individuals in [...]. This may indicate a low-quality site." (e.g. https://gnomad.broadinstitute.org/variant/8-7581908-A-G?dataset=gnomad_r3)

            # Skip any variant with AN (i.e. called genotypes) lower than 50% of samples, as flagged on gnomAD
            if chr == 'chrY':
                # chrY: use 2 * number of males 
                # Note: on chrX, maximum AN_male will also be anmaxy, but filtering needs to be done using anmax since maximum AN is twice as high.
                if newfields["an"] < anmaxy * 0.5:
                    Log(f"AN is lower than 50% of maximum ({Comma(anmaxy)}) for chrY|pos|originalbase|variantbase (skipped)", f"{chr}|{pos}|{originalbase}|{variantbase}")
                    continue
            else:
                # Other chromosomes
                if newfields["an"] < anmax * 0.5:
                    Log(f"AN is lower than 50% of maximum ({Comma(anmax)}) for chr|pos|originalbase|variantbase (skipped)", f"{chr}|{pos}|{originalbase}|{variantbase}")
                    continue



            # Process rsIDs (these are an int field in MySQL, without the "rs")
            # Usually there's only one rsID, but occasionally it's e.g. "rs1583111101;rs57725689;rs959273279"
            for rsid in nsort(rsids.split(";")):
                if rsid == '.':
                    rsid = "\\N"
                else:
                    m = rx(r"^rs(\d+)$", rsid)
                    if m:
                        rsid = m[0]
                    else:
                        Die(f"Error: Couldn't match rsid '{rsid}' in line:\n\n{line}\n\n")


                # Parse CSQs (Ensembl VEP output)
                # Split by comma
                for con in csq.split(","):
                    # Split by pipe
                    con = con.split("|")
                    # --fields="Allele,Gene,Feature_type,Feature,Consequence,HGVSp,Protein_position,Amino_acids,CHECK_REF"
                    if len(con) != 9:
                        Die("Error: Expected 9 fields in CSQ, but got:\n" + str(con))
                    
                    # Parse
                    allele = con[0]
                    ensg = con[1]
                    feature_type = con[2]
                    enst = con[3]
                    consequence = con[4]
                    hgvsp = con[5]
                    site = con[6]
                    aas = con[7]
                    check_ref = con[8]

                    # Verify
                    if feature_type != 'Transcript':
                        Die("Error: Expected feature_type to be 'Transcript' for missense_variants, but got:\n" + feature_type)

                    # Consequence
                    # e.g. missense_variant&NMD_transcript_variant
                    # if consequence != 'missense_variant':
                    #     Die("Error: Expected consequence to be 'missense_variant', but got:\n" + consequence)
                    #     # if 'missense_variant' in consequence:
                    #     #     Log("consequence is not 'missense_variant' for consequence (kept)", consequence)
                    #     # else:
                    #     #     Die("Error: Expected consequence to be 'missense_variant', or at least to contain it, but got:\n" + consequence)
                    consequences = consequence.split("&")
                    if 'missense_variant' not in consequences:
                        # Die("Error: Expected consequence to contain 'missense_variant', but got:\n" + consequence)
                        Log("consequences do not contain 'missense_variant' for this allele for consequences (skipped)", consequence)
                        Log("consequences do not contain 'missense_variant' for this allele for chr|pos|originalbase|variantbase (skipped)", f"{chr}|{pos}|{originalbase}|{variantbase}")
                        continue

                    # Verify feature type
                    if feature_type != 'Transcript':
                        Die("Error: Expected feature_type to be 'Transcript' for missense_variants, but got:\n" + feature_type)

                    # Check if site is a simple residue number or a range
                    # e.g. "169" or "169-170"
                    if rx(r'^\d+$', site):
                        site = int(site)
                    else:
                        Log("site isn't a single residue for site (skipped)", site)
                        Log("site isn't a single residue for originalbaselength|variantbaselength (skipped)", f"{len(originalbase)}|{len(variantbase)}")
                        Log("site isn't a single residue for originalbase|variantbase (skipped)", f"{originalbase}|{variantbase}")
                        Log("site isn't a single residue for csq (skipped)", csq)
                        continue

                    if hgvsp != '':
                        Die("Error: Expected hgvsp to be '', but got:\n" + hgvsp)
                    if check_ref != '':
                        Die("Error: Expected check_ref to be '', but got:\n" + check_ref)

                    # Parse original and variant from aas
                    m = rx(r'^([A-Z])/([A-Z])$', aas)
                    if m:
                        original = m[0]
                        variant = m[1]
                    else:
                        Die(f"Error: Couldn't parse aa string '{aas}' in line:\n" + line)

                    # Get UniProt accession etc.
                    # via SQL
                    # query = Query(f"SELECT DISTINCT name, acc, canon, species, ensg FROM uniens WHERE enst='{enst}'")
                    # if Numrows(query) == 0:
                    #     Log("no acc in uniens for enst (skipped)", enst)
                    #     continue
                    # elif Numrows(query) > 1:
                    #     Log("multiple accs in uniens for enst (kept)", enst)
                    #     # continue
                    # Log("found at least one acc in uniens for enst (kept)", enst)
                    # via dict
                    Log("total ensts that were looked for in uniens", enst)
                    if enst not in uniens:
                        Log("no acc in uniens for enst (skipped)", enst)
                        Log("no acc in uniens for chr|pos|originalbase|variantbase (skipped)", f"{chr}|{pos}|{originalbase}|{variantbase}")
                        continue
                    elif len(uniens[enst]) > 1:
                        Log("multiple accs in uniens for enst (kept)", enst)
                        Log("multiple accs in uniens for chr|pos|originalbase|variantbase (kept)", f"{chr}|{pos}|{originalbase}|{variantbase}")
                        # continue
                    Log("found at least one acc in uniens for enst (kept)", enst)
                    Log("found at least one acc in uniens for chr|pos|originalbase|variantbase (kept)", f"{chr}|{pos}|{originalbase}|{variantbase}")

                    # for (name, acc, canon, species, uniens_ensg) in query:
                    for (name, acc, canon, species, uniens_ensg) in uniens[enst]:

                        if uniens_ensg != ensg:
                            Die(f"Error: Expected uniens ensg for enst '{enst}' to be '{uniens_ensg}', but got '{ensg}'")

                        # Get UniProt sequence
                        # via SQL
                        # query = Query(f"SELECT DISTINCT seq FROM uniseq WHERE acc='{acc}' AND type IN ('UniProt', 'UniIso')")
                        # (uniseq) = FetchOne(query)
                        # via dict
                        # This would crash with a key error if the acc didn't have a sequence attached
                        uniseq = uniseqs[acc]

                        # Check if position is within UniProt sequence
                        if site > len(uniseq):
                            # Log("site is outside UniProt sequence (skipped) for acc", acc)
                            # Log("site is outside UniProt sequence (skipped) for acc|site|original|variant", f"{acc}|{site}|{original}|{variant}")
                            # continue
                            Die(f"Error: Site is outside UniProt sequence (skipped) for acc '{acc}' site '{site}' original '{original}' variant '{variant}' (shouldn't happen since this is based on Ensembl VEP 108, which should match up perfectly with uniens)")

                        # Check if original matches UniProt sequence
                        if original != uniseq[site - 1]:
                            # Log("original doesn't match UniProt sequence (skipped)", original, uniseq[site - 1])
                            # continue
                            Die(f"Error: Original doesn't match UniProt sequence: Should be '{uniseq[site - 1]}', but found '{original}' for acc '{acc}' site '{site}' original '{original}' variant '{variant}' (shouldn't happen since this is based on Ensembl VEP 108, which should match up perfectly with uniens)")

                        # Insert into MySQL
                        # Query(f"INSERT INTO {table} SET name='{name}', acc='{acc}', canon='{canon}', ensg='{ensg}', enst='{enst}', species='{species}', site='{site}', original='{original}', variant='{variant}', uniprot='{uniprot}', cosmic_curated='{cosmic_curated}', dbsnp='{dbsnp}', tcga='{tcga}', thousand_genomes='{thousand_genomes}', esp='{esp}', clinvar='{clinvar}', topmed='{topmed}', exac='{exac}', gnomad='{gnomad}', clingen='{clingen}', sources='{sourcecount}'")
                        # Write to temporary file (for LOAD DATA INFILE) (with NULL for the primary key, id)
                        # print(f"NULL\t{name}\t{acc}\t{canon}\t{ensg}\t{enst}\t{species}\t{site}\t{original}\t{variant}\t{source}{fieldstring}", file=out)
                        print(f"0\t{chr}\t{pos}\t{originalbase}\t{variantbase}\t{rsid}\t{name}\t{acc}\t{canon}\t{ensg}\t{enst}\t{species}\t{site}\t{original}\t{variant}\t{source}{fieldstring}", file=out)
                        

                        Log("successfully inserted variant for chr|pos|originalbase|variantbase", f"{chr}|{pos}|{originalbase}|{variantbase}")
                        Log("successfully inserted variant for acc|site|original|variant", f"{acc}|{site}|{original}|{variant}")
                        Log("successfully inserted variant for acc|site", f"{acc}|{site}")
                        Log("successfully inserted variant for acc", acc)
                        Log("successfully inserted variant for canon|site|original|variant", f"{canon}|{site}|{original}|{variant}")
                        Log("successfully inserted variant for canon|site", f"{canon}|{site}")
                        Log("successfully inserted variant for canon", canon)
                        Log("successfully inserted variant for ensg", ensg)
                        Log("successfully inserted variant for enst", enst)
                        Log("successfully inserted variant for species", species)
                        Log("successfully inserted variant for originalbaselength|variantbaselength", f"{len(originalbase)}|{len(variantbase)}")
                        Log("successfully inserted variant for originalbase|variantbase", f"{originalbase}|{variantbase}")
                    
Time(1)
Show(lim=50, sort=True)

print("\nDone!")
