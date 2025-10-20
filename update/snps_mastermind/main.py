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

table = "snps_mastermind"
chrsynfile = "../../.vep/homo_sapiens/108_GRCh38/chr_synonyms.txt"

infile = Args(1, "[input file]", "input/mastermind_cited_variants_reference-2024.01.03-grch38.vep.missense.vcf.gz")
# outfile_raw = f"tmp-sql-insert-raw.txt"
outfile = f"tmp-sql-insert.txt"

Clear(table)

# Get chromosome synonyms from Ensembl VEP cache
# cat ../../.vep/homo_sapiens/108_GRCh38/chr_synonyms.txt | \grep --color=auto -P "chr(\d+|M|X|Y)\t"
# >> None, i.e. the mapping is always from cryptic chromosome IDs such as NC_000001.11 to e.g. chr1. chr1 is never mapped to anything else.
# cat ../../.vep/homo_sapiens/108_GRCh38/chr_synonyms.txt | g "^NC_0000\d+" | nsort
# >> All of these only have one mapping, apart from X and Y (to chrX and X etc.)
# Load ~/.vep/homo_sapiens/108_GRCh38chrsynonyms.txt:
# NC_000024.10	chrY
# NC_000004.12	chr4
# NC_000001.11	chr1
print("Getting chromosome synonyms from Ensembl VEP cache file '{chrsynfile}':")
chrsyn = {}
with open(chrsynfile) as f:
    for line in tq(f):
        (value, map) = line.rstrip().split("\t")
        # Only retain mappings to chr1 etc.
        if rx(r"^chr(\d{1,2}|M|X|Y)$", map):
            chrsyn[value] = map

# Get uniens
print("\nGetting enst-to-name/acc/canon/species/ensg mappings from table 'uniens':")
uniens = {}
for (enst,) in tq(Query(f"SELECT DISTINCT enst FROM uniens WHERE species='human'")):
    uniens[enst] = FetchList(Query(f"SELECT DISTINCT name, acc, canon, species, ensg FROM uniens WHERE enst='{enst}'"))

# Get uniseq
print("\nGetting acc-to-UniProt-sequence mappings from table 'uniseq':")
uniseqs = {}
for (acc,) in tq(Query(f"SELECT DISTINCT acc FROM uniseq WHERE species='human' AND type IN ('UniProt', 'UniIso')")):
    uniseqs[acc] = FetchOne(Query(f"SELECT DISTINCT seq FROM uniseq WHERE acc='{acc}' AND type IN ('UniProt', 'UniIso')"))



# Start
# print(f"\nReading '{infile}' and inserting into table '{table}':")
print(f"\nReading '{infile}' and writing to temporary file '{outfile}':")

# # Disable keys during inserts
# Query(f"ALTER TABLE {table} DISABLE KEYS")

# Read line by line
Time(1)
paper1_stats = []
paper2_stats = []
paper_stats = []
with gzip.open(infile, 'rt') as f:
    with open(outfile, 'w') as out:

        # Print header (not really required, LOAD DATA LOCAL INFILE will ignore it below)
        # For VEP tables, this is: (chr pos originalbase variantbase rsid) name acc canon ensg enst species site original variant
        print("id\tchr\tpos\toriginalbase\tvariantbase\tname\tacc\tcanon\tensg\tenst\tspecies\tsite\toriginal\tvariant\tpapers1\tpapers2\tpapers", file=out)

        for line in tq(f):

            # Skip header
            if line.startswith("#"):
                continue

            e = line.rstrip().split("\t")

            rsid = e[2]

            # Parse tags
            papers1 = None
            papers2 = None
            papers = None
            csq = None

            tags = e[7].split(";")
            for tag in tags:
                if tag.startswith("MMCNT1="):
                    papers1 = tag[7:]
                elif tag.startswith("MMCNT2="):
                    papers2 = tag[7:]
                elif tag.startswith("MMCNT3="):
                    papers = tag[7:]
                elif tag.startswith("CSQ="):
                    csq = tag[4:]
            
            # Official documentation (input/README-Nov-2018.pdf):
            # MMCNT1 (most specific) – cDNA-level exact matches. This is the number of articles that mention the variant at the nucleotide level in either the title/abstract or the full-text.
            # >> cdna
            # MMCNT2 – cDNA-level possible matches. This is the number of articles with nucleotide-level matches (from 1) plus articles with protein-level matches in which the publication did not specify the cDNA-level change, meaning they could be referring to this nucleotide-level variant but there is insufficient data in these articles to determine conclusively.
            # >> cdna_or_protein
            # MMCNT3 (most sensitive) – This is the number of articles citing any variant resulting in the same biological effect as this variant. This includes the articles from MMCNT1 and MMCNT2 plus articles with alternative cDNA-level variants that result in the same protein effect.
            # >> cdna_or_protein_or_same_effect
            # >> This is the most relevant count for me since these citations are identical at the protein level.

            # In other words (from https://raw.githubusercontent.com/ensembl-variation/VEP_plugins/master/Mastermind.pm, an Ensembl VEP plugin for Mastermind):
            # 'MMCNT1' is the count of Mastermind articles with cDNA matches for a specific variant;
            # 'MMCNT2' is the count of Mastermind articles with variants either explicitly matching at the cDNA level or given only at protein level;
            # 'MMCNT3' is the count of Mastermind articles including other DNA-level variants resulting in the same amino acid change;
            # 'MMID3' is the Mastermind variant identifier(s), as gene:key. Link to the Genomenon Mastermind Genomic Search Engine;   
            # >> Since we only care about the protein-level change, MMCNT3 is the most relevant count to use.     

            # cdna = e['MMCNT1']
            # cdna_or_protein = e['MMCNT2']
            # cdna_or_protein_or_same_effect = e['MMCNT3']
            # papers = e['MMCNT3']

            # Papers (MMCNT1)
            if papers1 is None:
                Die(f"Error: Missing MMCNT1 (papers1) tag in line:\n" + line)
            # Papers (MMCNT2)
            if papers2 is None:
                Die(f"Error: Missing MMCNT2 (papers2) tag in line:\n" + line)
            # Papers (MMCNT3)
            if papers is None:
                Die(f"Error: Missing MMCNT3 (papers) tag in line:\n" + line)



            # Check if variants are larger than a single nucleotide
            
            # zcat input/mastermind_cited_variants_reference-2024.01.03-grch38.vep.missense.vcf.gz | perl -ne '@a = split(/\t/); if (((length($a[3]) == 1) and (length($a[4]) == 1)) and ($_ =~ /missense_variant/)) { print }' | wcl
            # >> 4,576,252 single-base missense variants
            # >> Only aroun a third are true SNVs (only affecting a single base)
            # zcat input/mastermind_cited_variants_reference-2024.01.03-grch38.vep.missense.vcf.gz | perl -ne '@a = split(/\t/); if (((length($a[3]) > 1) or (length($a[4]) > 1)) and ($_ =~ /missense_variant/)) { print }' | wcl
            # >> 12,663,956 missense or indel variants that affect more than one base
            # >> Two thirds of variants are larger. In principle, since I don't care about the mutational mechanism of the variant, but only about the amino acid change, I'll keep these larger variants (but only if they only affect a single amino acid).
            # zcat input/mastermind_cited_variants_reference-2024.01.03-grch38.vep.missense.vcf.gz | g -v ^# | cut -f4,5 | perl -ne 'chomp; @a = split(/\t/); if ((length($a[0]) == length($a[1])) and (length($a[0]) > 1)) { print "EQUAL LENGTH AND LENGTH > 1\t$_\n" }' | wc -l
            # >> 12,663,956
            # >> There are many missense variants where more than one base is affected but the length of originalbase and variantbase is >1 and equal. In fact, for all of them.
            # zcat input/mastermind_cited_variants_reference-2024.01.03-grch38.vep.missense.vcf.gz | perl -ne '@a = split(/\t/); if (((length($a[3]) > 3) or (length($a[4]) > 3)) and ($_ =~ /missense_variant/)) { print }' | wcl
            # >> Only 4,088 variants are larger than 3 bases (multi-AA).
            # zcat input/mastermind_cited_variants_reference-2024.01.03-grch38.vep.missense.vcf.gz | perl -ne '@a = split(/\t/); if (((length($a[3]) > 3) or (length($a[4]) > 3)) and ($_ =~ /missense_variant/)) { print }' | head 
            # >> Looking at these, they tend to be extremely large, actually. Skipping these.
            # Also, there's always only a single alternative allele in each line, so I don't need to handle these:
            # zcat input/mastermind_cited_variants_reference-2024.01.03-grch38.vep.missense.vcf.gz | g -v "^#" | cut -f4,5 | g ','
            # >> None! >> OK

            chr = e[0]
            pos = e[1]
            rsids = e[2]
            originalbase = e[3]
            variantbase = e[4]
            # Use readable chromosome synonym instead of the provided chromosome accessions (e.g. use chr1 instead of NC_000001.11)
            chr = chrsyn[chr]
            # if len(originalbase) == 1 and len(variantbase) == 1:
            #     Log("original and variant are a single nucleotide (i.e. a true SNV) for chr|pos|originalbase|variantbase (kept)", f"{chr}|{pos}|{originalbase}|{variantbase}")
            #     Log("original and variant are a single nucleotide (i.e. a true SNV) for originalbase|variantbase (kept)", f"{originalbase}|{variantbase}")
            if len(originalbase) != len(variantbase):
                Log("original and variant lengths differ (must be a frameshift) for chr|pos|originalbase|variantbase (skipped)", f"{chr}|{pos}|{originalbase}|{variantbase}")
                Log("original and variant lengths differ (must be a frameshift) for originalbase|variantbase (skipped)", f"{originalbase}|{variantbase}")
                continue
            # if len(originalbase) > 5:
            #     Log("original larger than 5 nucleotides (must affect multiple aas) for chr|pos|originalbase|variantbase (skipped)", f"{chr}|{pos}|{originalbase}|{variantbase}")
            #     Log("original larger than 5 nucleotides (must affect multiple aas) for originalbase|variantbase (skipped)", f"{originalbase}|{variantbase}")
            #     continue
            # if len(variantbase) > 5:
            #     Log("variant larger than 5 nucleotides (must affect multiple aas) for chr|pos|originalbase|variantbase (skipped)", f"{chr}|{pos}|{originalbase}|{variantbase}")
            #     Log("variant larger than 5 nucleotides (must affect multiple aas) for originalbase|variantbase (skipped)", f"{originalbase}|{variantbase}")
            #     continue
            
            
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
                    # Not every consequence is "missense_variant", but at least one will be, e.g.:

                    # HGVSG=NC_000008.11:g.39071318G>A;
                    # GENE=ADAM9;
                    # MMCNT1=0;
                    # MMCNT2=1;
                    # MMCNT3=1;
                    # MMID3=ADAM9:D538N;
                    # MMURI3=https://mastermind.genomenon.com/detail?mutation=NC_000008.11%3Ag.39071318G%3EA&ref=cvr;
                    # CSQ=
                    # A|ENSG00000168615|Transcript|ENST00000379917|missense_variant&NMD_transcript_variant||538|D/N|,
                    # A|ENSG00000168615|Transcript|ENST00000468065|missense_variant&NMD_transcript_variant||538|D/N|,
                    # A|ENSG00000168615|Transcript|ENST00000481873|missense_variant&NMD_transcript_variant||538|D/N|,
                    # A|ENSG00000168615|Transcript|ENST00000487273|missense_variant||538|D/N|,
                    # A|ENSG00000168615|Transcript|ENST00000676617|missense_variant&NMD_transcript_variant||538|D/N|,
                    # A|ENSG00000168615|Transcript|ENST00000676643|missense_variant||538|D/N|,
                    # A|ENSG00000168615|Transcript|ENST00000676669|missense_variant&NMD_transcript_variant||538|D/N|,
                    # A|ENSG00000168615|Transcript|ENST00000676765|missense_variant||538|D/N|,
                    # A|ENSG00000168615|Transcript|ENST00000677004|missense_variant||569|D/N|,
                    # A|ENSG00000168615|Transcript|ENST00000677137|missense_variant&NMD_transcript_variant||538|D/N|,
                    # A|ENSG00000168615|Transcript|ENST00000677165|missense_variant&NMD_transcript_variant||538|D/N|,
                    # A|ENSG00000168615|Transcript|ENST00000677359|missense_variant&NMD_transcript_variant||538|D/N|,
                    # A|ENSG00000168615|Transcript|ENST00000677582|missense_variant||538|D/N|,
                    # A|ENSG00000168615|Transcript|ENST00000677908|synonymous_variant&NMD_transcript_variant||480|K|

                    # Die("Error: Expected consequence to contain 'missense_variant', but got:\n" + consequence)
                    # Log("consequences do not contain 'missense_variant' for consequences (skipped)", consequence)
                    # Log("consequences do not contain 'missense_variant' for chr|pos|originalbase|variantbase (skipped)", f"{chr}|{pos}|{originalbase}|{variantbase}")
                    continue

                # Verify feature type
                if feature_type != 'Transcript':
                    Die("Error: Expected feature_type to be 'Transcript' for missense_variants, but got:\n" + feature_type)

                # Check if site is a simple residue number or a range
                # e.g. "169" or "169-170"
                if rx(r'^\d+$', site):
                    site = int(site)
                else:
                    # Note: This check should work, but Mastermind sometimes seems to give a single residue position (the starting position) rather than a range. Added aa sequence length check below.
                    Log("site isn't a single residue for chr|pos|originalbase|variantbase (skipped)", f"{chr}|{pos}|{originalbase}|{variantbase}")
                    Log("site isn't a single residue for originalbase|variantbase (skipped)", f"{originalbase}|{variantbase}")
                    Log("site isn't a single residue for originalbaselength|variantbaselength (skipped)", f"{len(originalbase)}|{len(variantbase)}")
                    Log("site isn't a single residue for original|variant (skipped)", f"{original}|{variant}")
                    Log("site isn't a single residue for originallength|variantlength (skipped)", f"{len(original)}|{len(variant)}")
                    continue

                if hgvsp != '':
                    Die("Error: Expected hgvsp to be '', but got:\n" + hgvsp)
                if check_ref != '':
                    Die("Error: Expected check_ref to be '', but got:\n" + check_ref)

                # Parse original and variant from aas
                # m = rx(r'^([A-Z])/([A-Z])$', aas)
                # More tolerant parsing (nonsense mutations (*) will be skipped below)
                m = rx(r'^([A-Z\*]+)/([A-Z\*]+)$', aas)
                if m:
                    original = m[0]
                    variant = m[1]
                    if (len(original) != 1) or (len(variant) != 1):
                        if len(original) == len(variant):
                            Log("variant affects multiple amino acids (same number for original and variant) for chr|pos|originalbase|variantbase (skipped)", f"{chr}|{pos}|{originalbase}|{variantbase}")
                            Log("variant affects multiple amino acids (same number for original and variant) for originalbase|variantbase (skipped)", f"{originalbase}|{variantbase}")
                            Log("variant affects multiple amino acids (same number for original and variant) for original|variant (skipped)", f"{original}|{variant}")
                            continue
                        else:
                            Log("variant affects multiple amino acids (different number for original and variant) for chr|pos|originalbase|variantbase (skipped)", f"{chr}|{pos}|{originalbase}|{variantbase}")
                            Log("variant affects multiple amino acids (different number for original and variant) for originalbase|variantbase (skipped)", f"{originalbase}|{variantbase}")
                            Log("variant affects multiple amino acids (different number for original and variant) for original|variant (skipped)", f"{original}|{variant}")
                            continue
                    if not Aa(original):
                        Log("variant has a non-aa original residue for chr|pos|originalbase|variantbase (skipped)", f"{chr}|{pos}|{originalbase}|{variantbase}")
                        Log("variant has a non-aa original residue for originalbase|variantbase (skipped)", f"{originalbase}|{variantbase}")
                        Log("variant has a non-aa original residue for original|variant (skipped)", f"{original}|{variant}")
                        Log("variant has a non-aa original residue for original (skipped)", original)
                        continue
                    if not Aa(variant):
                        Log("variant has a non-aa variant residue for chr|pos|originalbase|variantbase (skipped)", f"{chr}|{pos}|{originalbase}|{variantbase}")
                        Log("variant has a non-aa variant residue for originalbase|variantbase (skipped)", f"{originalbase}|{variantbase}")
                        Log("variant has a non-aa variant residue for original|variant (skipped)", f"{original}|{variant}")
                        Log("variant has a non-aa variant residue for variant (skipped)", variant)
                        continue
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
                if enst not in uniens:
                    Log("no acc in uniens for enst (skipped)", enst)
                    continue
                elif len(uniens[enst]) > 1:
                    Log("multiple accs in uniens for enst (kept)", enst)
                    # continue
                Log("found at least one acc in uniens for enst (kept)", enst)

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

                    # Insert into table 'snps_mastermind'
                    # Query(f"INSERT INTO {table} SET name='{name}', acc='{acc}', canon='{canon}', ensg='{ensg}', enst='{enst}', species='{species}', site='{site}', original='{original}', variant='{variant}', papers='{papers}'")
                    # Write to temporary file (for LOAD DATA INFILE) (with NULL for the primary key, id)
                    # print(f"0\t{name}\t{acc}\t{canon}\t{ensg}\t{enst}\t{species}\t{site}\t{original}\t{variant}\t{papers}", file=out)
                    print(f"0\t{chr}\t{pos}\t{originalbase}\t{variantbase}\t{name}\t{acc}\t{canon}\t{ensg}\t{enst}\t{species}\t{site}\t{original}\t{variant}\t{papers1}\t{papers2}\t{papers}", file=out)
                    

                    Log("successfully inserted variant for acc|site|original|variant", f"{acc}|{site}|{original}|{variant}")
                    Log("successfully inserted variant for acc|site", f"{acc}|{site}")
                    Log("successfully inserted variant for acc", acc)
                    Log("successfully inserted variant for name|site|original|variant", f"{name}|{site}|{original}|{variant}")
                    Log("successfully inserted variant for name|site", f"{name}|{site}")
                    Log("successfully inserted variant for name", name)
                    Log("successfully inserted variant for canon|site|original|variant", f"{canon}|{site}|{original}|{variant}")
                    Log("successfully inserted variant for canon|site", f"{canon}|{site}")
                    Log("successfully inserted variant for canon", canon)
                    Log("successfully inserted variant for ensg", ensg)
                    Log("successfully inserted variant for enst", enst)
                    Log("successfully inserted variant for species", species)
                    Log("successfully inserted variant for originalbaselength|variantbaselength", f"{len(originalbase)}|{len(variantbase)}")
                    Log("successfully inserted variant for originalbase|variantbase", f"{originalbase}|{variantbase}")
                    Log("successfully inserted variant for originallength|variantlength", f"{len(original)}|{len(variant)}")
                    Log("successfully inserted variant for original|variant", f"{original}|{variant}")
                    
                    paper1_stats.append(int(papers1))
                    paper2_stats.append(int(papers2))
                    paper_stats.append(int(papers))

Time(1)
Show(lim=50, sort=True)

# Show paper statistics
print("\nPapers per variant (MMCNT1, papers1):")
Characterise(paper1_stats)
print()
print("\nPapers per variant (MMCNT2, papers2):")
Characterise(paper2_stats)
print()
print("\nPapers per variant (MMCNT3, papers):")
Characterise(paper_stats)
print()

# # Unique outfile_raw
# Starttime()
# print(f"Make temporary file '{outfile_raw}' unique, writing to '{outfile}'...")
# Run("Unique", f"cat {outfile_raw} | sort | uniq > {outfile}")
# Stoptime()

# Load data from temporary file
print(f"Importing data from temporary file '{outfile}' into table '{table}'...")
Starttime()
# # Disable keys during inserts
# Query(f"ALTER TABLE {table} DISABLE KEYS")
# Load data from temporary file
Query(f"LOAD DATA LOCAL INFILE '{outfile}' INTO TABLE {table} FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES")
Stoptime()

# # Re-enable keys
# print(f"Enabling keys for table '{table}'...")
# Starttime()
# Query(f"ALTER TABLE {table} ENABLE KEYS")
# Stoptime()

# Sort table
print(f"Sorting table '{table}' by species, name, acc, site, original, variant...")
Starttime()
Query(f"ALTER TABLE {table} ORDER BY species, name, acc, site, original, variant")
Stoptime()

Optimize(table)

# Run("Remove temporary output file", f"rm -fv {outfile_raw}")
Run("Remove temporary output file", f"rm -fv {outfile}")

print("\nDone!")
