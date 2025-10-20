#!/usr/bin/env python3
"""
Run: Run entire pipeline
"""

# Initialize

import pickle
from blang_mysql import *
from blang import *

table = "snps_gnomad"

# Start

# Download
Run("Download", "download.py")

# Note: While gnomAD provides .bgz files, there is no benefit to .bgz either for running VEP on it or for any later step here. It'll only produce larger files.
# Previously this pipeline was using .bgz everywhere, but I've switched to .gz for the final output to save space.

# VEP (only run if output file doesn't exist yet)
# Go through all .vcf.bgz files in input directory
Cd("tmp")
# for vep_infile in glob.glob("../input/gnomad.exomes.r2.1.1.sites.*.liftover_grch38.vcf.bgz"):         # v2 exomes individual chromos only
# for vep_infile in glob.glob("../input/gnomad.genomes.r2.1.1.sites.*.liftover_grch38.vcf.bgz"):        # v2 genomes individual chromos only
# for vep_infile in glob.glob("../input/gnomad.genomes.v3.1.2.sites.chr*.vcf.bgz"):                     # v3 genomes individual chromos only
# for vep_infile in glob.glob("../input/gnomad.exomes.v4.0.sites.chrY.vcf.bgz"):                          # v4 exomes chrY only
for vep_infile in glob.glob("../input/*.vcf.bgz"):
    vep_outfile = vep_infile.replace("../input/", "../output/").replace(".vcf.bgz", ".vep.vcf.gz")
    if not Exists(vep_outfile):
        # # 32 cores (fastest, but difficult to find nodes right now, apparently)
        # Run(f"Run VEP on '{vep_infile}' to get protein-level consequence information", f'~/scripts/qsubcm.sh 32 3 ~/update/ensembl-vep/vep -i {vep_infile} -o {vep_outfile} --cache --offline --force -v --fork 32 --buffer_size 10000 --format vcf --vcf --compress_output gzip --coding_only --fields="Allele,Gene,Feature_type,Feature,Consequence,HGVSp,Protein_position,Amino_acids,CHECK_REF"')
        # # 8 cores (still kind of difficult)
        # Run(f"Run VEP on '{vep_infile}' to get protein-level consequence information", f'~/scripts/qsubcm.sh 8 3 ~/update/ensembl-vep/vep -i {vep_infile} -o {vep_outfile} --cache --offline --force -v --fork 8 --buffer_size 10000 --format vcf --vcf --compress_output gzip --coding_only --fields="Allele,Gene,Feature_type,Feature,Consequence,HGVSp,Protein_position,Amino_acids,CHECK_REF"')
        # 4 cores
        Run(f"Run VEP on '{vep_infile}' to get protein-level consequence information", f'~/scripts/qsubcm.sh 4 3 ~/update/ensembl-vep/vep -i {vep_infile} -o {vep_outfile} --cache --offline --force -v --fork 4 --buffer_size 10000 --format vcf --vcf --compress_output gzip --coding_only --fields="Allele,Gene,Feature_type,Feature,Consequence,HGVSp,Protein_position,Amino_acids,CHECK_REF"')
Cd("..")
Waitforjobs()


# Only keep header lines and "missense_variants"
Cd("tmp")
# for vep_infile in glob.glob("../output/gnomad.exomes.r2.1.1.sites.*.liftover_grch38.vep.vcf.gz"):    # v2 exomes individual chromos only
# for vep_infile in glob.glob("../output/gnomad.genomes.r2.1.1.sites.*.vep.vcf.gz"):                   # v2 genomes only
# for vep_infile in glob.glob("../output/gnomad.genomes.v3.1.2.sites.chr*.vep.vcf.gz"):                # v3 genomes individual chromos only
# for vep_infile in glob.glob("../output/gnomad.*.v4.0.sites.chr*.vep.vcf.gz"):                        # v4 genomes and exomes, individual chromos only
# for vep_infile in glob.glob("../output/gnomad.exomes.v4.0.sites.chrY.vep.vcf.gz"):                   # v4 exomes chrY only
for vep_infile in glob.glob("../output/*.vep.vcf.gz"):
    vep_outfile = vep_infile.replace(".vep.vcf.gz", ".vep.missense.vcf.gz")
    if not Exists(vep_outfile):
        Run(f"Filter '{vep_infile}' and keep only header and missense_variants", f"~/scripts/qsub.sh zcat {vep_infile} \| grep -P '^#\|missense_variant' \| gzip \> {vep_outfile}")
Cd("..")
Waitforjobs()

# Run cached queries here (to avoid any clashes in the individual jobs)
# Write uniens (cached query) if the pickle doesn't exist yet
tmpfile = "tmp-uniens.pkl"
if not Exists(tmpfile):
    print("Getting enst-to-name/acc/canon/species/ensg mappings from table 'uniens':")
    uniens = {}
    for (enst,) in tq(Query(f"SELECT DISTINCT enst FROM uniens WHERE species='human'")):
        uniens[enst] = FetchList(Query(f"SELECT DISTINCT name, acc, canon, species, ensg FROM uniens WHERE enst='{enst}'"))
    # Write to file
    with open(tmpfile, "wb") as f:
        pickle.dump(uniens, f)
    print(f"Wrote to '{tmpfile}'")
# Write uniseq (cached query) if the pickle doesn't exist yet
tmpfile = "tmp-uniseq.pkl"
if not Exists(tmpfile):
    print("\nGetting acc-to-UniProt-sequence mappings from table 'uniseq':")
    uniseqs = {}
    for (acc,) in tq(Query(f"SELECT DISTINCT acc FROM uniseq WHERE species='human' AND type IN ('UniProt', 'UniIso')")):
        uniseqs[acc] = FetchOne(Query(f"SELECT DISTINCT seq FROM uniseq WHERE acc='{acc}' AND type IN ('UniProt', 'UniIso')"))
    # Write to file
    with open(tmpfile, "wb") as f:
        pickle.dump(uniseqs, f)
    print(f"Wrote to '{tmpfile}'")

# Run main.py, which writes SQL insert queries to temporary files (queries will be run below)
Time(1)
Cd("tmp")
# for infile in nsort(glob.glob("../output/gnomad.exomes.r2.1.1.sites.*.liftover_grch38.vep.missense.vcf.gz")):    # v2 exomes individual chromos only (the rest are already run)
# for infile in nsort(glob.glob("../output/gnomad.*.v4.0.sites.chr*.vep.missense.vcf.gz")):                        # v4 genomes and exomes, individual chromos only
# for infile in nsort(glob.glob("../output/gnomad.exomes.v4.0.sites.chrY.vep.missense.vcf.gz")):                   # v4 exomes chrY only
for infile in nsort(glob.glob("../output/*.vep.missense.vcf.gz")):
    outfile = infile.replace(".vep.missense.vcf.gz", ".vep.missense.vcf.tmp_sql_insert.txt")
    if not Exists(outfile):
        # Run(f"Parse output from '{infile}' into temporary file", f"~/scripts/qsub.sh ../main.py {infile}")
        Run(f"Parse output from '{infile}' into temporary file", f"~/scripts/qsub_standard.sh ../main.py {infile}")
Cd("..")
Waitforjobs()
Time(1)

# Remove complete .vep.vcf.gz files
Run("Remove complete .vep.vcf.gz files (huge, only need to keep missense)", "rm -fv output/*.vep.vcf.gz")

# Clear table
Time(2)
Clear("snps_gnomad")
Time(2)

# Disable keys for inserts
Time(3)
Query(f"ALTER TABLE {table} DISABLE KEYS")
Time(3)

# Do SQL inserts via LOAD DATA LOCAL INFILE
Time(4)
Cd("tmp")
# for outfile in nsort(glob.glob("../output/gnomad.*.v4.0.sites.chr*.vep.missense.vcf.tmp_sql_insert.txt")):                    # v4 genomes and exomes, individual chromos only
for outfile in nsort(glob.glob("../output/*.vep.missense.vcf.tmp_sql_insert.txt")):

    # Load data from temporary files
    print(f"Importing data from temporary file '{outfile}' into table '{table}'...")

    # Load data from temporary file
    Query(f"LOAD DATA LOCAL INFILE '{outfile}' INTO TABLE {table} FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES")
    Time(4)

Cd("..")

# Re-enable keys
print(f"Enabling keys for table '{table}'...")
Time(5)
Query(f"ALTER TABLE {table} ENABLE KEYS")
Time(5)

# Reorder table
print(f"Reordering table '{table}'...")
Time(6)
Query(f"ALTER TABLE {table} ORDER BY species, name, acc, site, original, variant, source")
# Query(f"ALTER TABLE {table} ORDER BY species, chr, pos, originalbase, variantbase, name, acc, site, original, variant, source")   # No benefit (I won't be querying by chr|pos)
Time(6)

# Reset id column
print(f"Resetting id column in table '{table}'...")
Time(7)
Query(f"ALTER TABLE {table} DROP COLUMN id, DROP PRIMARY KEY, ADD COLUMN id int unsigned NOT NULL AUTO_INCREMENT FIRST, ADD PRIMARY KEY (id)")
Time(7)

Time(8)
Optimize(table)
Time(8)



print("Done!")
