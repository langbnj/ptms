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
from blang_mysql import *
from blang import *

table = "snps_clinvar"

# Args(0, "", "")
(rel, rel2) = Args(2, "[release version (rel)] [release date (rel2)]", "2023-07 20230722")
summaryfile = f"input/variant_summary_{rel}.txt.gz"
infile = f"input/clinvar_{rel2}_combined.vep.vcf.gz"

Clear(table)

# Start



# Read variant_summary file into memory
variant_summary = {}
print(f"Loading '{summaryfile}' into memory:")
with gzip.open(summaryfile, 'rt') as f:
    for line in tq(f):
        #AlleleID	Type	Name	GeneID	GeneSymbol	HGNC_ID	ClinicalSignificance	ClinSigSimple	LastEvaluated	RS# (dbSNP)	nsv/esv (dbVar)	RCVaccession	PhenotypeIDS	PhenotypeList	Origin	OriginSimple	Assembly	ChromosomeAccession	Chromosome	Start	Stop	ReferenceAllele	AlternateAllele	Cytogenetic	ReviewStatus	NumberSubmitters	Guidelines	TestedInGTR	OtherIDs	SubmitterCategories	VariationID	PositionVCF	ReferenceAlleleVCF	AlternateAlleleVCF
        # 15045	single nucleotide variant	NM_017547.4(FOXRED1):c.1289A>G (p.Asn430Ser)	55572	FOXRED1	HGNC:26927	Pathogenic	1	Oct 01, 2010	267606830	-	RCV000000016	MONDO:MONDO:0032624,MedGen:C4748791,OMIM:618241	Mitochondrial complex 1 deficiency, nuclear type 19	germline	germline	GRCh37	NC_000011.9	11	126147412	126147412	na	na	11q24.2	no assertion criteria provided	1	-	N	ClinGen:CA113794,UniProtKB:Q96CU9#VAR_064571,OMIM:613622.0002	1	6	126147412	A	G

        # Skip header
        if line.startswith("#"):
            continue

        e = line.rstrip().split("\t")
        alleleid = e[0]

        variant_summary[alleleid] = e


# Read VEP VCF
print(f"\nReading '{infile}' and inserting into table '{table}':")

# e.g.
# # 1	69134	2205837	A	G	.	.	ALLELEID=2193183;CLNDISDB=MeSH:D030342,MedGen:C0950123;CLNDN=Inborn_genetic_diseases;CLNHGVS=NC_000001.11:g.69134A>G;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Likely_benign;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;GENEINFO=OR4F5:79501;MC=SO:0001583|missense_variant;ORIGIN=1;CSQ=G|ENSG00000186092|Transcript|ENST00000641515|missense_variant||36|E/G|
# ALLELEID=2193183
# CLNDISDB=MeSH:D030342,MedGen:C0950123
# CLNDN=Inborn_genetic_diseases
# CLNHGVS=NC_000001.11:g.69134A>G
# CLNREVSTAT=criteria_provided,_single_submitter
# CLNSIG=Likely_benign
# CLNVC=single_nucleotide_variant
# CLNVCSO=SO:0001483
# GENEINFO=OR4F5:79501
# MC=SO:0001583|missense_variant
# ORIGIN=1
# CSQ=G|ENSG00000186092|Transcript|ENST00000641515|missense_variant||36|E/G|

# Read line by line
mafs = []
macs = []
with gzip.open(infile, 'rt') as f:
    for line in tq(f):
        # Skip header
        if line.startswith("#"):
            continue

        e = line.rstrip().split("\t")

        # Parse tags
        alleleid = None
        csq = None

        tags = e[7].split(";")
        for tag in tags:
            # Allele ID, in order to retrieve its information from the variant_summary file
            if tag.startswith("ALLELEID="):
                alleleid = tag[9:]
            elif tag.startswith("CSQ="):
                csq = tag[4:]
        
        # ALLELEID (for mapping to variant_summary file)
        if alleleid is None:
            Die(f"Error: Missing ALLELEID '{alleleid}' in line:\n" + line)
        Log("total alleleids", alleleid)

        # Check if variants are larger than a single nucleotide
        
        # zcat input/clinvar_20230722_combined.vep.vcf.gz | perl -ne '@a = split(/\t/); if (((length($a[3]) == 1) and (length($a[4]) == 1)) and ($_ =~ /missense_variant/)) { print }' | wcl
        # >> 1,153,075 variants (the vast majority) are true SNVs (only affecting a single base)
        # zcat input/clinvar_20230722_combined.vep.vcf.gz | perl -ne '@a = split(/\t/); if (((length($a[3]) > 1) or (length($a[4]) > 1)) and ($_ =~ /missense_variant/)) { print }' | wcl
        # >> 5,608 variants are larger (mostly 2 nt, e.g. CA>TG)
        # >> Since I don't care about the mutational mechanism of the variant, but only about the amino acid change, I'll keep these larger variants (but only if they only affect a single amino acid).
        
        # Also, ther's always only a single alternative allele in each line, so I don't need to handle these:
        # zcat clinvar_20230722_combined.vep.vcf.gz | g -v "^#" | cut -f4,5 | g ','
        # >> None
        originalbase = e[3]
        variantbase = e[4]
        if len(originalbase) == 1 and len(variantbase) == 1:
            Log("original and variant are a single nucleotide (i.e. a true SNV) for alleleid (kept)", alleleid)
        if len(originalbase) != len(variantbase):
            Log("original and variant lengths differ (must be a frameshift) for alleleid (kept)", alleleid)
        if len(originalbase) > 1:
            Log("original larger than a single nucleotide for alleleid (kept)", alleleid)
            # continue
        if len(variantbase) > 1:
            Log("variant larger than a single nucleotide for alleleid (kept)", alleleid)
            # continue

        # VEP's CSQ string (only exists if there is any coding consequence for the variant)
        if csq is None:
            Log("no consequence for alleleid (skipped)", alleleid)
            continue
        elif csq.find("missense_variant") == -1:
            Log("no missense_variant in csq for alleleid (skipped)", alleleid)
            continue

        # Get annotation from variant_summary
        ann = variant_summary[alleleid]

        # Parse annotation

        # Clinical significance plain text
        # zcat variant_summary_2023-06.txt.gz | cut -f7 | suq
        #     305 Affects
        #     312 conflicting data from submitters
        #     714 association
        #     949 risk factor
        #    1412 no interpretation for the single variant
        #    3233 other
        #    3792 drug response
        #   21972 not provided
        #   38320 Pathogenic/Likely pathogenic
        #   73812 Benign/Likely benign
        #  135682 Likely pathogenic
        #  201273 Conflicting interpretations of pathogenicity
        #  295057 Pathogenic
        #  412783 Benign
        # 1203289 Likely benign
        # 2070399 Uncertain significance	
        # >> The ones I should capture are: 
        # Pathogenic/Likely pathogenic
        # Likely pathogenic
        # Pathogenic
        # ...but actually: I can simply use ClinSigSimple (= 1)
        # This is actually superior since some "Uncertain significance" entries have a ClinSigSimple of 1, while others are 0. Would lose the 1s otherwise.
        clin = ann[6]

    	# Clinical significance yes/no
        clinsig = ann[7]

        # Date last evaluated
        evaluated = ann[8]
        if evaluated == '-':
            evaluated = 'NULL'
        else:
            # evaluated = f"evaluated=STR_TO_DATE('{evaluated}','%b %d, %Y'), ";
            evaluated = f"STR_TO_DATE('{evaluated}','%b %d, %Y')";

    	# rsid
        rsid = ann[9]
        if rsid == '-1':
            # Blank
            rsid = ''
        else:
            # Add "rs" to the beginning
            rsid = f"rs{rsid}";

        # Phenotype IDs
        # e.g. MedGen:CN517202|MeSH:D030342,MedGen:C0950123 (comma groups are equivalent terms, apparently, while pipes divide term groups)
        phenoids_str = ann[12]
        phenoids = []
        for phenoid in nsort(phenoids_str.split('|')):
            # I've checked that MedGen:CN517202 ("not provided") doesn't co-occur with any other equivalent (comma-separated) terms:
            # zcat variant_summary_2023-07.txt.gz | cut -f13 | suq
            # zcat variant_summary_2023-07.txt.gz | cut -f13 | suq | g ",MedGen:CN517202"
            # zcat variant_summary_2023-07.txt.gz | cut -f13 | suq | g "MedGen:CN517202,"
            if phenoid != 'MedGen:CN517202':
                phenoids.append(phenoid)
        phenoids = '|'.join(phenoids)
        
        # Phenotype plain text
        # e.g. not provided|Cardiovascular phenotype
        phenotypes_str = ann[13]
        phenotypes = []
        for phenotype in nsort(phenotypes_str.split('|')):
            # I've checked that "not provided" doesn't co-occur with any other equivalent (comma-separated) terms:
            # zcat variant_summary_2023-07.txt.gz | cut -f14 | suq
            # zcat variant_summary_2023-07.txt.gz | cut -f14 | suq | g ",(not provided|not specified)"
            # zcat variant_summary_2023-07.txt.gz | cut -f14 | suq | g "(not provided|not specified),"
            if phenotype != 'not provided' and phenotype != 'not specified':
                phenotypes.append(phenotype)
        phenotypes = '|'.join(phenotypes)

        # Origin
        origin = ann[14]
        
        # OriginSimple
        # zcat input/variant_summary_2023-06.txt.gz | cut -f16 | suq
        #       1 OriginSimple
        #      65 tested-inconclusive
        #    4647 germline/somatic
        #    6289 not applicable
        #   12889 somatic
        #   34717 not provided
        #   97946 unknown
        # 4308279 germline
        originsimple = ann[15]

        # Review status
        status = ann[24]
        #      1 ReviewStatus
        #     46 practice guideline
        #   1284 no interpretation for the single variant
        #  18484 reviewed by expert panel
        #  21712 no assertion provided
        #  36841 criteria provided, conflicting interpretations
        # 106019 criteria provided, multiple submitters, no conflicts
        # 122536 no assertion criteria provided
        # 503533 criteria provided, single submitter
        # >> Not really interpretable
        # $conflicts = 0;
        # $conflicts = 1 if ($status =~ /conflicting interpretations/);

        submitters = ann[25]
        
        # Included in the Genetic Testing Registry (GTR) panel yes/no
        genetictest = ann[27]
        if genetictest == 'Y':
            genetictest = 1
        elif genetictest == 'N':
            genetictest = 0



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
                # Check if CSQ contains missense_variant
                Log("non-missense_variant for consequence (skipped)", consequence)
                Log("non-missense_variant for con (one of the others will be a missense_variant for this alleleid) (skipped)", '|'.join(con))
                continue

            # Verify feature type
            if feature_type != 'Transcript':
                Die("Error: Expected feature_type to be 'Transcript' for missense_variants, but got:\n" + feature_type)

            # Check if site is a simple residue number or a range
            # e.g. "169" or "169-170"
            if rx(r'^\d+$', site):
                site = int(site)
            else:
                Log("site isn't a single residue for alleleid (skipped)", alleleid)
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
            query = Query(f"SELECT DISTINCT name, acc, canon, species, ensg FROM uniens WHERE enst='{enst}'")
            if Numrows(query) == 0:
                Log("no acc in uniens for enst (skipped)", enst)
                continue
            elif Numrows(query) > 1:
                Log("multiple accs in uniens for enst (kept)", enst)
                # continue
            Log("found at least one acc in uniens for enst (kept)", enst)

            # (name, acc, canon, species, ensembl_ensg) = FetchOne(query)
            for (name, acc, canon, species, uniens_ensg) in query:

                if uniens_ensg != ensg:
                    Die(f"Error: Expected uniens ensg for enst '{enst}' to be '{uniens_ensg}', but got '{ensg}'")

                # Get UniProt sequence
                query = Query(f"SELECT DISTINCT seq FROM uniseq WHERE acc='{acc}' AND type IN ('UniProt', 'UniIso')")
                (uniseq) = FetchOne(query)

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

                # Insert into table
                q = f"INSERT INTO {table} SET name='{name}', acc='{acc}', canon='{canon}', ensg='{ensg}', enst='{enst}', species='{species}', site='{site}', original='{original}', variant='{variant}', clinsig='{clinsig}', clin='{clin}', genetictest='{genetictest}', submitters='{submitters}', evaluated={evaluated}, rsid='{rsid}', phenoids='{phenoids}', phenotype='{Esc(phenotypes)}', origin='{origin}', originsimple='{originsimple}', status='{status}'"
                q = q.replace("=''", "=NULL")
                Query(q)

                Log("successfully inserted variant for acc|site|original|variant", f"{acc}|{site}|{original}|{variant}")
                Log("successfully inserted variant for acc|site", f"{acc}|{site}")
                Log("successfully inserted variant for acc", acc)
                Log("successfully inserted variant for canon|site|original|variant", f"{canon}|{site}|{original}|{variant}")
                Log("successfully inserted variant for canon|site", f"{canon}|{site}")
                Log("successfully inserted variant for canon", canon)
                Log("successfully inserted variant for ensg", ensg)
                Log("successfully inserted variant for enst", enst)
                Log("successfully inserted variant for species", species)
                Log("successfully inserted variant for alleleid", alleleid)
                Log("successfully inserted variant for originalbaselength|variantbaselength", f"{len(originalbase)}|{len(variantbase)}")
                Log("successfully inserted variant for originalbase|variantbase", f"{originalbase}|{variantbase}")



Show(lim=20, sort=True)
# Show("non-missense_variant for consequences (skipped)")

# Reorder table
print(f"\nReordering table '{table}'...")
Query(f"ALTER TABLE {table} ORDER BY species, name, acc, site, original, variant")

Optimize(table)

print("\nDone!")
