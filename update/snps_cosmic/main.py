#!/usr/bin/env python3
"""
Main: Main script
"""

# Initialize
import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
import gzip
import pickle
from blang_mysql import *
from blang import *

rel = 99
ensembl_rel = 108

table = "snps_cosmic"
cgcfile = f"input/Cosmic_CancerGeneCensus_v{rel}_GRCh38.tsv"
chrsynfile = f"../../.vep/homo_sapiens/{ensembl_rel}_GRCh38/chr_synonyms.txt"

infile = Args(1, "[input file]", f"input/Cosmic_GenomeScreensMutant_Normal_v{rel}_GRCh38.vep.missense.vcf.gz")
source = infile
# Remove input/
source = re.sub(r"^.+/", "", source)
# Remove any file endings
source = re.sub(r"\..+$", "", source)
# Remove _v99_GRCh38 at end
source = re.sub(rf"_v{rel}_GRCh38$", "", source)
# Set cell_line to 1 if this input file is from the cell lines project
cell_line = 0
if rx(r"^cell", source, re.IGNORECASE):
    cell_line = 1
outfile = f"tmp-sql-insert-{source}.txt"

# def norm (startsite, stopsite, original, variant):
#     print(startsite, stopsite, original, variant)
# 
#     # Normalise start (removing identical amino acids at the beginning, e.g. in AP>AS)
#     while original[0] == variant[0]:
#         # Remove this position from original and variant
#         original = original[1:]
#         variant = variant[1:]
#         # Increase startsite
#         startsite += 1
# 
#     # Normalise stop (removing identical amino acids at the end, e.g. in PA>SA)
#     while original[-1] == variant[-1]:
#         # Remove this position from original and variant
#         original = original[:-1]
#         variant = variant[:-1]
#         # Decrease stopsite
#         stopsite -= 1
# 
#     return startsite, stopsite, original, variant
# print(norm(687, 688, "AP", "AS"))
# print(norm(687, 689, "APA", "ASA"))
# d()
# sys.exit()

# Clear(table)

# Get Cancer Gene Census (CGC) gene list. These are confirmed cancer driver genes.
# Read cgcfile using pandas
t = read_tsv(cgcfile)

# Symbols
# Contains 743 symbols
# >> 728 map to 'uniprot' 'symbol'
# via hgnc_aliases:
# FetchSet(Query("SELECT DISTINCT p.symbol FROM hgnc_aliases a, uniprot p WHERE p.species='human' AND a.symbol=p.symbol AND a.alias IN ('" + "', '".join(nsort(cgc["Gene Symbol"])) + "')"))
# >> 735 symbols in uniprot
# A1CF, ABI1, ABL1, ABL2, ACKR3, ACSL3, ACSL6, ACVR1, ACVR1B, ACVR2A, AFDN, AFF1, AFF3, AFF4, AKAP9, AKT1, AKT2, AKT3, ALDH2, ALK, AMER1, ANK1, APC, APOBEC3B, AR, ARAF, ARHGAP5, ARHGAP26, ARHGAP35, ARHGEF10, ARHGEF10L, ARHGEF12, ARID1A, ARID1B, ARID2, ARNT, ASPM, ASPSCR1, ASXL1, ASXL2, ATF1, ATIC, ATM, ATP1A1, ATP2B3, ATR, ATRX, AXIN1, AXIN2, B2M, BAP1, BARD1, BAX, BAZ1A, BCL2, BCL2L12, BCL3, BCL6, BCL7A, BCL9, BCL9L, BCL10, BCL11A, BCL11B, BCLAF1, BCOR, BCORL1, BCR, BIRC3, BIRC6, BLM, BMP5, BMPR1A, BRAF, BRCA1, BRCA2, BRD3, BRD4, BRIP1, BTG1, BTG2, BTK, BUB1B, CACNA1D, CALR, CAMTA1, CANT1, CARD11, CARS1, CASP3, CASP8, CASP9, CBFA2T3, CBFB, CBL, CBLB, CBLC, CCDC6, CCNB1IP1, CCNC, CCND1, CCND2, CCND3, CCNE1, CCR4, CCR7, CD28, CD74, CD79A, CD79B, CD209, CD274, CDC73, CDH1, CDH10, CDH11, CDH17, CDK4, CDK6, CDK12, CDKN1A, CDKN1B, CDKN2A, CDKN2C, CDX2, CEBPA, CEP43, CEP89, CHCHD7, CHD2, CHD4, CHEK2, CHIC2, CHST11, CIC, CIITA, CLIP1, CLP1, CLTC, CLTCL1, CNBD1, CNBP, CNOT3, CNTNAP2, CNTRL, COL1A1, COL2A1, COL3A1, COX6C, CPEB3, CREB1, CREB3L1, CREB3L2, CREBBP, CRLF2, CRNKL1, CRTC1, CRTC3, CSF1R, CSF3R, CSMD3, CTCF, CTNNA1, CTNNA2, CTNNB1, CTNND1, CTNND2, CUL3, CUX1, CXCR4, CYLD, CYP2C8, CYSLTR2, DAXX, DCAF12L2, DCC, DCTN1, DDB2, DDIT3, DDR2, DDX3X, DDX5, DDX6, DDX10, DEK, DGCR8, DICER1, DNAJB1, DNM2, DNMT3A, DROSHA, EBF1, ECT2L, EED, EGFR, EIF1AX, EIF3E, EIF4A2, ELF3, ELF4, ELK4, ELL, ELN, EML4, EP300, EPAS1, EPHA3, EPHA7, EPS15, ERBB2, ERBB3, ERBB4, ERC1, ERCC2, ERCC3, ERCC4, ERCC5, ERG, ESR1, ETNK1, ETV1, ETV4, ETV5, ETV6, EWSR1, EXT1, EXT2, EZH2, EZR, FADD, FAM47C, FAM131B, FAM135B, FANCA, FANCC, FANCD2, FANCE, FANCF, FANCG, FAS, FAT1, FAT3, FAT4, FBLN2, FBXO11, FBXW7, FCGR2B, FCRL4, FEN1, FES, FEV, FGFR1, FGFR2, FGFR3, FGFR4, FH, FHIT, FIP1L1, FKBP9, FLCN, FLI1, FLNA, FLT3, FLT4, FNBP1, FOXA1, FOXL2, FOXO1, FOXO3, FOXO4, FOXP1, FOXR1, FSTL3, FUBP1, FUS, GAS7, GATA1, GATA2, GATA3, GLI1, GMPS, GNA11, GNAQ, GNAS, GOLGA5, GOLPH3, GOPC, GPC3, GPC5, GPHN, GRIN2A, GRM3, GSK3B, H3-3A, HERPUD1, HEY1, HGF, HIF1A, HIP1, HLA-A, HLF, HMGA1, HMGA2, HMGN2P46, HNF1A, HNRNPA2B1, HOOK3, HOXA9, HOXA11, HOXA13, HOXC11, HOXC13, HOXD11, HOXD13, HRAS, HSP90AA1, HSP90AB1, ID3, IDH1, IDH2, IGF2BP2, IKBKB, IKZF1, IKZF3, IL2, IL6ST, IL7R, IL21R, IRF4, IRS4, ISX, ITGAV, ITK, JAK1, JAK2, JAK3, JAZF1, JUN, KAT6A, KAT6B, KAT7, KCNJ5, KDM5A, KDM5C, KDM6A, KDR, KDSR, KEAP1, KIAA1549, KIF5B, KIT, KLF4, KLF6, KLK2, KMT2A, KMT2C, KMT2D, KNL1, KNSTRN, KRAS, KTN1, LARP4B, LASP1, LATS1, LATS2, LCK, LCP1, LEF1, LEPROTL1, LHFPL6, LIFR, LMNA, LMO1, LMO2, LPP, LRIG3, LRP1B, LSM14A, LYL1, LYN, LZTR1, MACC1, MAF, MAFB, MALT1, MAML2, MAP2K1, MAP2K2, MAP2K4, MAP3K1, MAP3K13, MAPK1, MAX, MB21D2, MDM2, MDM4, MDS2, MECOM, MED12, MEN1, MET, MGMT, MITF, MLF1, MLH1, MLLT1, MLLT3, MLLT6, MLLT10, MLLT11, MN1, MNX1, MPL, MRTFA, MSH2, MSH6, MSI2, MSN, MTCP1, MTOR, MUC1, MUC4, MUC6, MUC16, MUTYH, MYB, MYC, MYCL, MYCN, MYD88, MYH9, MYH11, MYO5A, MYOD1, N4BP2, NAB2, NACA, NBEA, NBN, NCKIPSD, NCOA1, NCOA2, NCOA4, NCOR1, NCOR2, NDRG1, NF1, NF2, NFATC2, NFE2L2, NFIB, NFKB2, NFKBIE, NIN, NKX2-1, NONO, NOTCH1, NOTCH2, NPM1, NR4A3, NRAS, NRG1, NSD1, NSD2, NSD3, NT5C2, NTHL1, NTRK1, NTRK2, NTRK3, NUMA1, NUP98, NUP214, NUTM1, NUTM2B, NUTM2D, OLIG2, OMD, P2RY8, PABPC1, PAFAH1B2, PALB2, PATZ1, PAX3, PAX5, PAX7, PAX8, PBRM1, PBX1, PCBP1, PCM1, PDCD1LG2, PDE4DIP, PDGFB, PDGFRA, PDGFRB, PER1, PHF6, PHOX2B, PICALM, PIERCE2, PIK3CA, PIK3CB, PIK3R1, PIM1, PLAG1, PLCG1, PML, PMS1, PMS2, POLD1, POLE, POLG, POLQ, POLR2A, POT1, POU2AF1, POU5F1, PPARG, PPFIBP1, PPM1D, PPP2R1A, PPP6C, PRCC, PRDM1, PRDM2, PRDM16, PREX2, PRF1, PRKACA, PRKAR1A, PRKCB, PRKD1, PRPF40B, PRRX1, PSIP1, PTCH1, PTEN, PTK6, PTPN6, PTPN11, PTPN13, PTPRB, PTPRC, PTPRD, PTPRK, PTPRT, PWWP2A, QKI, RABEP1, RAC1, RAD17, RAD21, RAD50, RAD51B, RAF1, RALGDS, RANBP2, RAP1B, RAP1GDS1, RARA, RB1, RBM10, RBM15, RECQL4, REL, RET, RFWD3, RGPD3, RGS7, RHOA, RHOH, RMI2, RNF43, RNF213, ROBO2, ROS1, RPL5, RPL10, RPL22, RPN1, RRAS2, RSPO2, RSPO3, RUNX1, RUNX1T1, S100A7, SALL4, SBDS, SDC4, SDHA, SDHAF2, SDHB, SDHC, SDHD, SEPTIN5, SEPTIN6, SEPTIN9, SET, SETBP1, SETD1B, SETD2, SETDB1, SF3B1, SFPQ, SFRP4, SGK1, SH2B3, SH3GL1, SHTN1, SIRPA, SIX1, SIX2, SKI, SLC34A2, SLC45A3, SMAD2, SMAD3, SMAD4, SMARCA4, SMARCB1, SMARCD1, SMARCE1, SMC1A, SMO, SND1, SNX29, SOCS1, SOX2, SOX21, SPECC1, SPEN, SPOP, SRC, SRGAP3, SRSF2, SRSF3, SS18, SS18L1, SSX1, SSX2, SSX4, STAG1, STAG2, STAT3, STAT5B, STAT6, STIL, STK11, STRN, SUB1, SUFU, SUZ12, SYK, TAF15, TAL1, TAL2, TBL1XR1, TBX3, TCEA1, TCF3, TCF7L2, TCF12, TCL1A, TEC, TENT5C, TERT, TET1, TET2, TFE3, TFEB, TFG, TFPT, TFRC, TGFBR2, THRAP3, TLX1, TLX3, TMEM127, TMPRSS2, TMSB4X, TNC, TNFAIP3, TNFRSF14, TNFRSF17, TOP1, TP53, TP63, TPM3, TPM4, TPR, TRA, TRAF7, TRB, TRIM24, TRIM27, TRIM33, TRIP11, TRRAP, TSC1, TSC2, TSHR, U2AF1, UBR5, USP6, USP8, USP9X, USP44, VAV1, VHL, VTI1A, WAS, WDCP, WIF1, WNK2, WRN, WT1, WWTR1, XPA, XPC, XPO1, YWHAE, ZBTB16, ZCCHC8, ZEB1, ZFHX3, ZMYM2, ZMYM3, ZNF331, ZNF384, ZNF429, ZNF479, ZNF521, ZNRF3, ZRSR2

# ENSGs
# cat cancer_gene_census.csv | g -v "ENSG\d+"
# Contains 734 ENSGs

# >> 9 rows without ENSG:
# IGH (immunoglobulin, skip)
# IGK (immunoglobulin, skip)
# IGL (immunoglobulin, skip)
# TRA (T-cell receptor, skip)
# TRB (T-cell receptor, skip)
# TRD (T-cell receptor, skip)
# MDS2 (found via uniprot symbol)
# MALAT1 (lncRNA, no protein)
# HMGN2P46 (found via uniprot symbol)
# >> OK! Getting MDS2 and HMGN2P46 via uniprot symbol.
# According to UniProt, MDS2 and HMGN2P46 are proteins:
# SELECT * FROM uniprot WHERE symbol IN ('MDS2', 'HMGN2P46');
# ...but according to Ensembl, they are a lncRNA and a pseudogene, with no protein entry in table 'ensembl':
# SELECT * FROM ensembl_gff3_gene g, ensembl e WHERE g.symbol IN ('MDS2', 'HMGN2P46') AND g.species='human' AND g.ensg=e.ensg;
# >> VEP 108 will in all probability not call protein-level consequences for these two.

# Out of the 734 listed ENSGs, 730 are in uniens (see ~/Documents/MySQL/snps_cosmic.sql).
# >> Going ahead with the 734 and ignoring the other 9. I can't map these specifically enough to use them.

# Get ENSGs from cgc data frame
print(f"\nGetting Cancer Gene Census (CGC) gene ENSG IDs from {cgcfile}:")
cgcs = {}
for i, e in tq(t.iterrows(), total=len(t)):
    syns = e["SYNONYMS"]
    symbol = e["GENE_SYMBOL"]
    tier = e["TIER"]
    Log("cgc: total cgc symbols", symbol)
    if pd.isnull(syns):
        Log("cgc: no ensg given for cgc symbol for symbol (skipped)", symbol)
        continue
    # Get all occurrences of ENSG\d+ (ignoring the ENSG version number) in synonyms
    for m in re.finditer(r"(ENSG\d+)", syns):
        ensg = m.group(0)
        cgcs[ensg] = tier
        Log(f"cgc: ensg given for cgc symbol for symbol (kept)", symbol)
        Log(f"cgc: ensg given for cgc symbol for ensg (kept)", ensg)
        Log(f"cgc: ensg given for tier {tier} cgc symbol for symbol (kept)", symbol)
        Log(f"cgc: ensg given for tier {tier} cgc symbol for ensg (kept)", ensg)
# Show()

# Get chromosome synonyms from Ensembl VEP cache
# cat ../../.vep/homo_sapiens/108_GRCh38/chr_synonyms.txt | \grep --color=auto -P "chr(\d+|M|X|Y)\t"
# >> None, i.e. the mapping is always from cryptic chromosome IDs such as NC_000001.11 to e.g. chr1. chr1 is never mapped to anything else.
# cat ../../.vep/homo_sapiens/108_GRCh38/chr_synonyms.txt | g "^NC_0000\d+" | nsort
# >> All of these only have one mapping, apart from X and Y (to chrX and X etc.)
# Load ~/.vep/homo_sapiens/108_GRCh38chrsynonyms.txt:
# NC_000024.10	chrY
# NC_000004.12	chr4
# NC_000001.11	chr1
print(f"\nGetting chromosome synonyms from Ensembl VEP cache file '{chrsynfile}':")
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
with gzip.open(infile, 'rt') as f:
    with open(outfile, 'w') as out:

        # Print header (not really required, LOAD DATA LOCAL INFILE will ignore it below)
        # For VEP tables, this is: (chr pos originalbase variantbase rsid) name acc canon ensg enst species site original variant
        print("id\tchr\tpos\toriginalbase\tvariantbase\tcosv\tsource\tcell_line\tname\tacc\tcanon\tensg\tenst\tspecies\tsite\toriginal\tvariant\tac\tan\taf\tcgc\tloc\tsubloc\thist\tsubhist\tpmid\tzyg\tloh\tsom\ttum", file=out)

        for line in tq(f):

            # Skip header
            if line.startswith("#"):
                continue

            e = line.rstrip().split("\t")

            cosv = e[2]

            # Parse tags
            ac = None
            csq = None

            tags = e[7].split(";")
            for tag in tags:
                if tag.startswith("SAMPLE_COUNT="):
                    ac = tag[13:]
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

            # SAMPLE_COUNT (allele count)
            if ac is None:
                Die(f"Error: Missing SAMPLE_COUNT (allele count) tag in line:\n" + line)



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
            cosv = e[2]
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
                    
                # Is this ENSG a CGC gene?
                cgc = 0
                # cgc = "\\N"
                if ensg in cgcs:
                    cgc = cgcs[ensg]

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

                # # Check if site is a simple residue number or a range
                # This check doesn't quite work - the "affected residues" here include silently affected ones, e.g.:
                # zcat input/Cosmic_CompleteTargetedScreensMutant_Normal_v99_GRCh38.vep.missense.vcf.gz | g "^16\t3778062\t[^\t]+\tGG\tAA\t"
                # 16	3778062	COSV105100373	GG	AA	.	.	GENE=CREBBP;TRANSCRIPT=ENST00000262367.9;STRAND=-;LEGACY_ID=COSM9656080;CDS=c.2061_2062delinsTT;AA=p.P688S;HGVSC=ENST00000262367.9:c.2061_2062delinsTT;HGVSP=ENSP00000262367.5:p.Pro688Ser;HGVSG=16:g.3778062_3778063delinsAA;SAMPLE_COUNT=1;IS_CANONICAL=y;TIER=1;SO_TERM=substitution;OLD_VARIANT=16:3778061:GGG/GAA;CSQ=AA|ENSG00000005339|Transcript|ENST00000262367|missense_variant||687-688|AP/AS|,AA|ENSG00000005339|Transcript|ENST00000382070|missense_variant||649-650|AP/AS|,AA|ENSG00000005339|Transcript|ENST00000570939|missense_variant||222-223|AP/AS|,AA|ENSG00000005339|Transcript|ENST00000571826|missense_variant||37-38|AP/AS|,AA|ENSG00000005339|Transcript|ENST00000572134|missense_variant||125-126|AP/AS|
                # 16	3778062	COSV105100373	GG	AA	.	.	GENE=CREBBP;TRANSCRIPT=ENST00000382070.7;STRAND=-;LEGACY_ID=COSM9656080;CDS=c.1947_1948delinsTT;AA=p.P650S;HGVSC=ENST00000382070.7:c.1947_1948delinsTT;HGVSP=ENSP00000371502.3:p.Pro650Ser;HGVSG=16:g.3778062_3778063delinsAA;SAMPLE_COUNT=1;IS_CANONICAL=n;TIER=1;SO_TERM=substitution;OLD_VARIANT=16:3778061:GGG/GAA;CSQ=AA|ENSG00000005339|Transcript|ENST00000262367|missense_variant||687-688|AP/AS|,AA|ENSG00000005339|Transcript|ENST00000382070|missense_variant||649-650|AP/AS|,AA|ENSG00000005339|Transcript|ENST00000570939|missense_variant||222-223|AP/AS|,AA|ENSG00000005339|Transcript|ENST00000571826|missense_variant||37-38|AP/AS|,AA|ENSG00000005339|Transcript|ENST00000572134|missense_variant||125-126|AP/AS|
                # "missense_variant||687-688|AP/AS"
                # >> A remains the same (though it has that G>A mutation), only P to S is an actual missense mutation.
                # >> It's not a normalisation issue since GG>AA can't be normalised further. It does mean, however, that I shouldn't use the AA position information here.
                # e.g. "169" or "169-170"
                if rx(r'^\d+$', site):
                    # Single aa variant (OK)
                    startsite = int(site)
                    stopsite = int(site)
                else:
                    # Multi-aa variant (possibly silent in all but one position, will attempt to normalise it below)
                    m = rx(r'^(\d+)-(\d+)$', site)
                    if m:
                        startsite = m[0]
                        stopsite = m[1]
                    else:
                        Die(f"Error: Couldn't parse site '{site}' in line:\n" + line)
                    # # Note: This check should work, but Mastermind sometimes seems to give a single residue position (the starting position) rather than a range. Added aa sequence length check below.
                    # Log("site isn't a single residue for chr|pos|originalbase|variantbase (skipped)", f"{chr}|{pos}|{originalbase}|{variantbase}")
                    # Log("site isn't a single residue for originalbase|variantbase (skipped)", f"{originalbase}|{variantbase}")
                    # Log("site isn't a single residue for originalbaselength|variantbaselength (skipped)", f"{len(originalbase)}|{len(variantbase)}")
                    # Log("site isn't a single residue for original|variant (skipped)", f"{original}|{variant}")
                    # Log("site isn't a single residue for originallength|variantlength (skipped)", f"{len(original)}|{len(variant)}")
                    # continue

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
                    
                    # originalbaselength and variantbaselength are identical, so originallength and variantlength should definitely be, too:
                    if len(original) != len(variant):
                        Die(f"Error: Expected original and variant lengths to be equal, but got: '{aas}'")
                    
                    # Check if original and variant length match the length given in site
                    if (stopsite - startsite + 1) != len(original):
                        Die(f"Error: Expected original and variant length to match site range given, but got: '{aas}'")
                    
                    # Attempt to normalise multi-aa variants if they are silent in all but one position
                    # e.g. "missense_variant||687-688|AP/AS"
                    # >> A remains the same (though it has that G>A mutation), only P to S is an actual missense mutation.
                    # >> It's not a normalisation issue since GG>AA can't be normalised further.
                    # Note: Only the amino acid level is being normalised. The genomic information will remain as before, e.g.:
                    # 1	12847848	COSV58611869	GACGTTGA	TGCGCTGG	.	.	GENE=HNRNPCL1;TRANSCRIPT=ENST00000317869.6;STRAND=-;LEGACY_ID=COSM97001;CDS=c.435_442delinsCCAGCGCA;AA=p.L148I;HGVSC=ENST00000317869.6:c.435_442delinsCCAGCGCA;HGVSP=ENSP00000365370.3:p.Leu148Ile;HGVSG=1:g.12847848_12847855delinsTGCGCTGG;SAMPLE_COUNT=1;IS_CANONICAL=y;SO_TERM=substitution;OLD_VARIANT=1:12847847:AGACGTTGA/ATGCGCTGG;CSQ=TGCGCTGG|ENSG00000179172|Transcript|ENST00000317869|missense_variant||145-148|RQRL/RQRI|
                    # >> GACGTTGA to CSQ=TGCGCTGG|ENSG00000179172|Transcript|ENST00000317869|missense_variant||145-148|RQRL/RQRI|
                    # >> 145-148|RQRL/RQRI
                    # >> ...turns into 148|L/I, but chr|pos|originalbase|variantbase will remain the same (GACGTTGA/TGCGCTGG).
                    
                    # Normalise start (removing identical amino acids at the beginning, e.g. in AP>AS)
                    while original[0] == variant[0]:
                        # Remove this position from original and variant
                        original = original[1:]
                        variant = variant[1:]
                        # Increase startsite
                        startsite += 1

                    # Normalise stop (removing identical amino acids at the end, e.g. in PA>SA)
                    while original[-1] == variant[-1]:
                        # Remove this position from original and variant
                        original = original[:-1]
                        variant = variant[:-1]
                        # Decrease stopsite
                        stopsite -= 1
                        
                    # Check if this is now a single-aa variant
                    if stopsite != startsite:
                        # Die(f"Error: Expected original and variant length to be 1 after amino acid-level normalisation, but got: {original}>{variant} ({startsite}-{stopsite})")
                        pass
                    else:
                        site = startsite
                    
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

                    # Insert into table 'snps_univar'
                    # Write to temporary file (for LOAD DATA INFILE) (with NULL for the primary key, id)
                    # id\tchr\tpos\toriginalbase\tvariantbase\tcosv\tsource\tcell_line\tname\tacc\tcanon\tensg\tenst\tspecies\tsite\toriginal\tvariant\tac\tan\taf\tcgc\tloc\tsubloc\thist\tsubhist\tpmid\tzyg\tloh\tsom\ttum
                    print(f"\\N\t{chr}\t{pos}\t{originalbase}\t{variantbase}\t{cosv}\t{source}\t{cell_line}\t{name}\t{acc}\t{canon}\t{ensg}\t{enst}\t{species}\t{site}\t{original}\t{variant}\t{ac}\t\\N\t\\N\t{cgc}\t\\N\t\\N\t\\N\t\\N\t\\N\t\\N\t\\N\t\\N\t\\N", file=out)
                    
                    Log(f"successfully inserted variant for chr|pos|originalbase|variantbase", f"{chr}|{pos}|{originalbase}|{variantbase}")
                    Log(f"successfully inserted variant for acc|site|original|variant", f"{acc}|{site}|{original}|{variant}")
                    Log(f"successfully inserted variant for CGC {cgc} acc|site|original|variant", f"{acc}|{site}|{original}|{variant}")
                    Log(f"successfully inserted variant for CGC {cgc} acc|site|original|variant", f"{acc}|{site}|{original}|{variant}")
                    Log(f"successfully inserted variant for acc|site", f"{acc}|{site}")
                    Log(f"successfully inserted variant for acc", acc)
                    Log(f"successfully inserted variant for name|site|original|variant", f"{name}|{site}|{original}|{variant}")
                    Log(f"successfully inserted variant for name|site", f"{name}|{site}")
                    Log(f"successfully inserted variant for name", name)
                    Log(f"successfully inserted variant for canon|site|original|variant", f"{canon}|{site}|{original}|{variant}")
                    Log(f"successfully inserted variant for canon|site", f"{canon}|{site}")
                    Log(f"successfully inserted variant for canon", canon)
                    Log(f"successfully inserted variant for ensg", ensg)
                    Log(f"successfully inserted variant for enst", enst)
                    Log(f"successfully inserted variant for cosv", cosv)
                    Log(f"successfully inserted variant for source", source)
                    Log(f"successfully inserted variant for species", species)
                    Log(f"successfully inserted variant for originalbaselength|variantbaselength", f"{len(originalbase)}|{len(variantbase)}")
                    Log(f"successfully inserted variant for originalbase|variantbase", f"{originalbase}|{variantbase}")
                    Log(f"successfully inserted variant for originallength|variantlength", f"{len(original)}|{len(variant)}")
                    Log(f"successfully inserted variant for original|variant", f"{original}|{variant}")
                    
Time(1)
Show(lim=50, sort=True)

# Load data from temporary file
print(f"\nImporting data from temporary file '{outfile}' into table '{table}'...")
Starttime()
Query(f"LOAD DATA LOCAL INFILE '{outfile}' INTO TABLE {table} FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES")
Stoptime()

Run("Remove temporary output file", f"rm -fv {outfile}")

print("\nDone!")
