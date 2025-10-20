#!/usr/bin/env python3 -u
"""
Download: Download source data
"""

# Initialize

import json
# from blang_mysql import *
from blang import *

# Release number
# rel = 98        # 2023-05-23 (https://cancer.sanger.ac.uk/cosmic, https://cancer.sanger.ac.uk/cosmic/release_notes)
rel = 99        # 2023-11-28 (https://cancer.sanger.ac.uk/cosmic, https://cancer.sanger.ac.uk/cosmic/release_notes)



# Define download function
def DownloadFile(filename, url=None):

    Cd("input")

    # Build URL if it isn't set yet:
    if url is None:
        # # Old style
        # url = f"https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v{rel}/{filename}"
        # New style
        url = f"https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch38/cosmic/v{rel}/{filename}_v{rel}_GRCh38.tar&bucket=downloads"
    
    # Check if filename has a type:
    if not re.search(r"\.\w+$", filename):
        filename += ".tar"
    
    # Get JSON
    s = Return(f"""curl -H "Authorization: Basic {key}" "{url}" 2> /dev/null""")

    # Get final URL
    url = json.loads(s)["url"]
    
    # # Make basename directory from {filename}
    # Run("Make directory", f"mkdir -p {os.path.dirname(filename)}")
    
    # Download from URL
    # print(f"wget --no-check-certificate -O {filename} '{url}'")
    # Run("Download", f"wget --no-check-certificate -O {filename} '{url}'")

    # # Old style
    # Run("Download", f"wget --no-check-certificate -O {os.path.basename(filename)} '{url}'")
    # New style
    Run("Download", f"wget --no-check-certificate -O {os.path.basename(filename)} '{url}'")
    
    # Unpack if necessary
    if re.search("\.tar$", filename):
        Run("Unpack", f"tar -xavpf {os.path.basename(filename)}")
        Run("Clean up", f"rm -fv {os.path.basename(filename)}")
    # if re.search("\.gz$", filename):
    #     Run("Unpack", f"gunzip -f {filename}")
    
    Cd("..")


# Get access key
# Insert email and password here
email = ""
password = ""
key = Return(f"""echo "{email}:{password}" | base64""")




# Download files

# # Old style
# 
# # # Cancer Gene Census: Apparently the CGC is the foundation of COSMIC (COSMIC Classic). What I’ve previously used from COSMIC was this curated Cancer Gene Census.
# # # This file is just a list of genes.
# # DownloadFile("cancer_gene_census.csv")
# # >> 738 unique genes
# 
# # Cancer Gene Census: Apparently the CGC is the foundation of COSMIC (COSMIC Classic). What I’ve previously used from COSMIC was this curated Cancer Gene Census.
# # This file is just a list of genes.
# # DownloadFile(f"Cosmic_CancerGeneCensus_Tsv_v{rel}_GRCh38.tar")
# # DownloadFile(f"Cosmic_MutantCensus_Tsv_v{rel}_GRCh38.tar")
# DownloadFile("cancer_gene_census.csv")
# # >> 738 unique genes
# 
# # The CosmicMutantExport covers all genes, not just the CGC ones. This is what I am using now.
# # zcat CosmicMutantExport.tsv.gz | cut -f1 | perl -npe 's/_ENST\d+$//' | suq | wcl
# # >> ~20000 unique genes
# # "COSMIC Mutation Data: A tab separated table of all COSMIC coding point mutations from targeted and genome wide screens from the current release. "
# # >> This was its name when v98 came out, but it has since been renamed.
# # https://cancer.sanger.ac.uk/cosmic/archive-download
# # DownloadFile("CosmicMutantExport.tsv.gz")
# 
# # VCF coding
# DownloadFile("VCF/CosmicCodingMuts.normal.vcf.gz")
# # VCF noncoding (but might in fact be coding according to VEP)
# DownloadFile("VCF/CosmicNonCodingVariants.normal.vcf.gz")
# 
# # VCF coding "non-normalised"
# DownloadFile("VCF/CosmicCodingMuts.vcf.gz")
# # VCF noncoding (but might in fact be coding according to VEP)
# DownloadFile("VCF/CosmicNonCodingVariants.vcf.gz")

# Should probably use the "normalised" VCFs since they're newer?
# https://cosmic-blog.sanger.ac.uk/cosmic-release-v91/
# This calls the normalised versions "5'-shifted", while the HGVSG strings were kept unchanged. It sounds a little as if they're avoiding to admit a previous error.

# Run both normalised and non-normalised through VEP and see if the outcome is the same:
# >> More or less, yes! "Normalising" seems to simplify some indels (a few hundred) to SNVs. There's only a very small difference. The normalised versions are most likely better.
# >> Use normalised versions

# New style

DownloadFile("Cosmic_CancerGeneCensus_Tsv")
DownloadFile("VCF/Cosmic_NonCodingVariants_VcfNormal")
DownloadFile("VCF/Cosmic_CompleteTargetedScreensMutant_VcfNormal")
DownloadFile("VCF/Cosmic_GenomeScreensMutant_VcfNormal")



# What's the difference between the old and new style downloads?

# Old ("Archive Downloads"): https://cancer.sanger.ac.uk/cosmic/archive-download
# New: https://cancer.sanger.ac.uk/cosmic/download/cosmic/v99/cancergenecensus

# The file sizes are quite different...how about the content, though?
# For the Cancer Gene Census CSV/TSV: The "New style" TSV contains ENSGs (as "Synonyms") for 3 additional genes. It also seems to have newer gene names (descriptions) and more systematic column names. It clearly seems newer.
# The file dates show that the "Old style" files are about a week newer than the "New style" files. Perhaps that simply means they weren't a priority, though, and that they aren't the primary files.

# How about the variant content?
# The "New style" VCF looks like it has been tidied up ("CNT" is now called "SAMPLE_COUNT", and the gene symbol and ENST were split into two fields)
# Weirdly, the "fileDate" of the "Old style" file is from November while the "New style" one has September, but again, that might just be because the "Old style" ones aren't the primary files.
# >> Get the new ones!

# VCF tags present in the new VCF files (same as the old):
# ~/update/snps_cosmic/input >> zcat Cosmic_GenomeScreensMutant_Normal_v99_GRCh38.vcf.gz Cosmic_NonCodingVariants_Normal_v99_GRCh38.vcf.gz Cosmic_CompleteTargetedScreensMutant_Normal_v99_GRCh38.vcf.gz | g -v ^# | cut -f8 | perl -ne 'while (/(^|;)([^=]+)=/g) { print "$2\n" }' | suq
# 74,556,715 SO_TERM        
# 74,556,715 SAMPLE_COUNT   
# 74,556,715 LEGACY_ID      
# 74,556,715 HGVSG          
# 66,415,151 TRANSCRIPT     
# 66,415,151 STRAND         
# 66,415,151 IS_CANONICAL   
# 66,415,151 GENE           
# 40,458,308 HGVSC          
# 40,458,308 CDS            
# 40,458,308 AA             
# 13,992,233 HGVSP          
#  5,474,531 TIER           
#  3,953,919 OLD_VARIANT    
# (Note: CellLinesProject_GenomeScreensMutant_Normal_v99_GRCh38.vep.missense.vcf.gz and CellLinesProject_NonCodingVariants_Normal_v99_GRCh38.vep.missense.vcf.gz have the same tags.)
# All present in e.g.:
# GENE=SKI;TRANSCRIPT=ENST00000378536.4;STRAND=+;LEGACY_ID=COSM6190280;CDS=c.249_251del;AA=p.P84del;HGVSC=ENST00000378536.4:c.249_251del;HGVSP=ENSP00000367797.4:p.Pro84del;HGVSG=1:g.2229015_2229017del;SAMPLE_COUNT=1;IS_CANONICAL=y;TIER=2;SO_TERM=deletion;OLD_VARIANT=1:2229014:CGCC/C
# GENE=SKI;
# TRANSCRIPT=ENST00000378536.4;
# STRAND=+;
# LEGACY_ID=COSM6190280;
# CDS=c.249_251del;
# AA=p.P84del;
# HGVSC=ENST00000378536.4:c.249_251del;
# HGVSP=ENSP00000367797.4:p.Pro84del;
# HGVSG=1:g.2229015_2229017del;
# SAMPLE_COUNT=1;
# IS_CANONICAL=y;
# TIER=2;
# SO_TERM=deletion;
# OLD_VARIANT=1:2229014:CGCC/C

# >> The VCFs are missing all the cancer type annotation. The only thing they have is "SAMPLE_COUNT" and "TIER". Everything else is uninteresting.

# TIER:
# https://cosmic-blog.sanger.ac.uk/cancer-gene-census-hallmarks-and-new-tier-system/
# "We now have a gold standard, tier 1, where a gene must possess documented evidence of activity that may drive cancer, as well as evidence that changes to protein activity promote oncogenic transformation. We also take into account the existence of somatic mutation patterns in cancer samples typical for tumour suppressor genes (usually characterised by a broad range of inactivating mutations) or oncogenes (usually characterised by well-defined hotspots of missense mutations)."
# "Tier 2 contains genes that have strong indications of playing a role in cancer but with less expansive evidence than tier 1, it currently contains 41 genes however these are not reflected in the website or downloads due to data constraints with this release. The current Cancer Gene Census in COSMIC just contains tier 1 genes, however, as tier 2 is being expanded, with an initial planned release of about 200 genes in November 2017. "
# >> Tier 1 is "gold standard", Tier 2 is lower confidence

# TSV versions of the VCF files above (for additional cancer type etc. annotation)
# DownloadFile(f"https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch38/cosmic/v{rel}/Cosmic_GenomeScreensMutant_Tsv_v{rel}_GRCh38.tar&bucket=downloads")
DownloadFile(f"Cosmic_GenomeScreensMutant_Tsv")
DownloadFile(f"Cosmic_NonCodingVariants_Tsv")
DownloadFile(f"Cosmic_CompleteTargetedScreensMutant_Tsv")
# These have the information I previously parsed, plus:
# - LOH: loss of heterozygosity
# - PUBMED_PMID:
# zcat Cosmic_GenomeScreensMutant_v99_GRCh38.tsv.gz | cut -f19 | g -v "^$" | g -v "^\d+$"
# >> Always a single PMID.



# Sample information (need this for primary etc.)
DownloadFile(f"Cosmic_Sample_Tsv")
# >> Can map this via the sample ID (COSS...) in the TSV files.
# Contains:
# - TUMOUR_SOURCE (= tum):
# zcat Cosmic_Sample_v99_GRCh38.tsv.gz | cut -f12 | suq
#       1 TUMOUR_SOURCE
#       5 adenoma adjacent to primary tumour
#      36 hyperplasia adjacent to primary tumour
#    6665 secondary
#   17479 recurrent
#   32642 metastasis
#  407054 primary
# 1078989 NS
# >> Usually not specified (NS), and otherwise ~90% primary. Not useful.
# >> It was previously more commonly annotated, perhaps because it didn't include the targeted screens:
# SELECT tum, COUNT(*) AS c FROM snps_cosmic_backup_before_vep GROUP BY tum ORDER BY c DESC;
# 	5,070,985
# primary	3,930,670
# metastasis	775,196
# recurrent	108,232
# secondary	9,407
# hyperplasia adjacent to primary tumour	1,966
# adenoma adjacent to primary tumour	72
# >> Nothing else in the sample annotation file is interesting enough to justify the effort, and most things are super sporadically filled.
# >> There is "ethnicity", but it's very sporadically filled and at a fine-grained free-text level that wouldn't map easily onto the gnomAD groups.
# >> Ignore




# Cancer Mutation Census:
# single dataset
DownloadFile(f"CancerMutationCensus_AllData_Tsv_v{rel}_GRCh38.tar", url=f"https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=GRCh38/cmc/v{rel}/CancerMutationCensus_AllData_Tsv_v{rel}_GRCh38.tar&bucket=downloads")
# >> Covers all genes (not just CGC)
# >> This doesn't actually have any cancer type information. What it does have is:
# - ONC_TSG
#   24503 oncogene, TSG, fusion
#   30364 oncogene, TSG
#   31425 TSG, fusion
#   47554 fusion
#   60597 oncogene, fusion
#   71963 oncogene
#  154857 TSG
# 4804551 
# >> TSG = tumour suppressor gene
# >> Mostly blank >> not useful
# >> Ignore
# - "SHARED_AA - Number of mutations seen in this amino acid position":
# >> This doesn't match COSMIC_SAMPLE_MUTATED (I've seen it be 3 for only 2 samples), but that could be due to zygosity.
# >> Not useful, though. AC/AF are the interesting numbers.
# - The number of COSMIC samples a mutation was tested in (COSMIC_SAMPLE_TESTED) vs. the number it was found in (COSMIC_SAMPLE_MUTATED), i.e. essentially AC / AN = AF.
# >> This is nice, but I already have AC in the VCFs. AN and AF don't add too much more information. I can already judge recurrence based on the VCFs.
# >> Ignore for now (added to do note in folder)
# zcat CancerMutationCensus_AllData_v99_GRCh38.tsv.gz | cut -f22 | suq
# zcat CancerMutationCensus_AllData_v99_GRCh38.tsv.gz | cut -f22 | suq | g -v "^4988150 " | g -o "^\s*\d+" | g -o "\d+" | sum
# >> 4,988,150 were tested in 44,567 samples
# >> only 237,665 (~5%) were tested in a different number of samples (either in more or fewer)
# - DISEASE - Diseases with > 1% samples mutated (or frequency > 0.01), where disease = Primary site(tissue) / Primary histology / Sub-histology = Samples mutated / Samples tested = Frequency
# >> e.g. soft_tissue/gastrointestinal_stromal_tumour/NS=9/83=10.84%;upper_aerodigestive_tract/carcinoma/squamous_cell_carcinoma=23/1022=2.25%
# - WGS_DISEASE - Same as DISEASE, but for whole-genome screen data only
# >> e.g. soft_tissue/gastrointestinal_stromal_tumour/NS=6/62=9.68%;upper_aerodigestive_tract/carcinoma/squamous_cell_carcinoma=23/1022=2.25%
# >> This is more limited than simply DISEASE. Since the sample numbers were already so small, it's not as if the frequencies get more accurate by only using WGS screens.
# >> I suppose it avoids bias from targeted testing, though.
# >> Keep both.
# - CLINVAR_TRAIT: e.g.
# Epileptic encephalopathy
# Malignant tumor of prostate << cancer
# Primary ciliary dyskinesia
# >> I think I have this parsed better in snps_clinvar.
# >> Ignore
# - GERP++_RS   (conservation score of "missing substitutions" - left blank for silent mutations, not sure why)
# - MIN_SIFT_SCORE
# - MIN_SIFT_PRED
# 1342810 T
# 1863010 
# 2019994 D
# - DNDS_DISEASE_QVAL_SIG: e.g. colon=5.29e-07;endometrioid=1.93e-07;gastric=2.42e-09
# >> "DNDS_DISEASE_QVAL - dn/ds diseases with significant q-values (q-value < 0.05), analysed from TCGA whole-exome data in COSMIC. Diseases are classified into AML, HNSCC, NSCLC, bladder, breast, cervical, colon, endometrioid, gastric, glioma, kidney, liver, melanoma, ovary, pancreatic, prostate, sarcoma, thyroid"
# >> Purely based on TCGA, not actually on COSMIC itself.
# - MUTATION_SIGNIFICANCE_TIER
#     399 2
#    1972 1
#   49603 3
# 5173840 Other
# >> "MUTATION_SIGNIFICANCE_TIER - Mutation significance. 1 - high significance, 2 - medium significance, 3 - low significance, Other - No predicted significance (other mutations)"
# >> Not sure what this is based on, the README doesn't say.
# The rest is boring ExAC/gnomAD columns.
# I was hoping for a good breakdown of cancer type involvement, and I suppose DNDS_DISEASE_QVAL_SIG is good for this, but it doesn't say where a mutation was observed.
# >> I'll still need the other TSVs as well.

# CGC tier 1:
# ('ENSG00000002834', 'ENSG00000005073', 'ENSG00000005339', 'ENSG00000006468', 'ENSG00000007237', 'ENSG00000007312', 'ENSG00000009709', 'ENSG00000010671', 'ENSG00000012048', 'ENSG00000015285', 'ENSG00000018408', 'ENSG00000019582', 'ENSG00000019991', 'ENSG00000023445', 'ENSG00000026103', 'ENSG00000029725', 'ENSG00000037280', 'ENSG00000039068', 'ENSG00000046889', 'ENSG00000047410', 'ENSG00000047932', 'ENSG00000047936', 'ENSG00000048462', 'ENSG00000049618', 'ENSG00000051108', 'ENSG00000051341', 'ENSG00000051382', 'ENSG00000055609', 'ENSG00000057657', 'ENSG00000062822', 'ENSG00000064012', 'ENSG00000065361', 'ENSG00000065526', 'ENSG00000065559', 'ENSG00000066455', 'ENSG00000066468', 'ENSG00000067082', 'ENSG00000067560', 'ENSG00000067842', 'ENSG00000067955', 'ENSG00000068078', 'ENSG00000068323', 'ENSG00000069399', 'ENSG00000070371', 'ENSG00000070404', 'ENSG00000071564', 'ENSG00000072062', 'ENSG00000072364', 'ENSG00000072694', 'ENSG00000073282', 'ENSG00000073578', 'ENSG00000073584', 'ENSG00000073614', 'ENSG00000073803', 'ENSG00000073921', 'ENSG00000076242', 'ENSG00000076685', 'ENSG00000077150', 'ENSG00000077782', 'ENSG00000078399', 'ENSG00000078403', 'ENSG00000078674', 'ENSG00000079102', 'ENSG00000079432', 'ENSG00000079805', 'ENSG00000079999', 'ENSG00000080824', 'ENSG00000081237', 'ENSG00000082805', 'ENSG00000082898', 'ENSG00000083093', 'ENSG00000083168', 'ENSG00000083799', 'ENSG00000083857', 'ENSG00000084676', 'ENSG00000085185', 'ENSG00000085224', 'ENSG00000085276', 'ENSG00000085832', 'ENSG00000087088', 'ENSG00000087460', 'ENSG00000088038', 'ENSG00000088256', 'ENSG00000089280', 'ENSG00000091483', 'ENSG00000091831', 'ENSG00000092820', 'ENSG00000095002', 'ENSG00000095015', 'ENSG00000096384', 'ENSG00000096968', 'ENSG00000097007', 'ENSG00000099949', 'ENSG00000099956', 'ENSG00000100030', 'ENSG00000100105', 'ENSG00000100311', 'ENSG00000100345', 'ENSG00000100393', 'ENSG00000100503', 'ENSG00000100644', 'ENSG00000100697', 'ENSG00000100721', 'ENSG00000100814', 'ENSG00000100815', 'ENSG00000101096', 'ENSG00000101115', 'ENSG00000101213', 'ENSG00000101972', 'ENSG00000102034', 'ENSG00000102145', 'ENSG00000102974', 'ENSG00000103126', 'ENSG00000103197', 'ENSG00000103522', 'ENSG00000104320', 'ENSG00000104365', 'ENSG00000104408', 'ENSG00000104419', 'ENSG00000104517', 'ENSG00000104884', 'ENSG00000104903', 'ENSG00000105173', 'ENSG00000105221', 'ENSG00000105369', 'ENSG00000105568', 'ENSG00000105639', 'ENSG00000105656', 'ENSG00000105662', 'ENSG00000105810', 'ENSG00000105976', 'ENSG00000106031', 'ENSG00000106462', 'ENSG00000106483', 'ENSG00000107485', 'ENSG00000107779', 'ENSG00000107807', 'ENSG00000107882', 'ENSG00000108091', 'ENSG00000108375', 'ENSG00000108654', 'ENSG00000108821', 'ENSG00000108924', 'ENSG00000108946', 'ENSG00000108953', 'ENSG00000109132', 'ENSG00000109471', 'ENSG00000109670', 'ENSG00000109685', 'ENSG00000109906', 'ENSG00000110092', 'ENSG00000110367', 'ENSG00000110395', 'ENSG00000110619', 'ENSG00000110713', 'ENSG00000110777', 'ENSG00000110841', 'ENSG00000110987', 'ENSG00000111252', 'ENSG00000111276', 'ENSG00000111642', 'ENSG00000112039', 'ENSG00000112081', 'ENSG00000112531', 'ENSG00000112561', 'ENSG00000112576', 'ENSG00000113263', 'ENSG00000113360', 'ENSG00000113522', 'ENSG00000113594', 'ENSG00000113721', 'ENSG00000113916', 'ENSG00000114354', 'ENSG00000114423', 'ENSG00000114861', 'ENSG00000115170', 'ENSG00000115524', 'ENSG00000115808', 'ENSG00000116016', 'ENSG00000116044', 'ENSG00000116062', 'ENSG00000116128', 'ENSG00000116132', 'ENSG00000116251', 'ENSG00000116560', 'ENSG00000116990', 'ENSG00000117118', 'ENSG00000117400', 'ENSG00000117713', 'ENSG00000118046', 'ENSG00000118058', 'ENSG00000118260', 'ENSG00000118503', 'ENSG00000118513', 'ENSG00000118689', 'ENSG00000118971', 'ENSG00000119335', 'ENSG00000119397', 'ENSG00000119414', 'ENSG00000119508', 'ENSG00000119535', 'ENSG00000119537', 'ENSG00000119772', 'ENSG00000119866', 'ENSG00000120217', 'ENSG00000120457', 'ENSG00000121067', 'ENSG00000121741', 'ENSG00000121879', 'ENSG00000121966', 'ENSG00000121989', 'ENSG00000122025', 'ENSG00000122406', 'ENSG00000122512', 'ENSG00000122566', 'ENSG00000122779', 'ENSG00000123080', 'ENSG00000123268', 'ENSG00000123364', 'ENSG00000123388', 'ENSG00000123473', 'ENSG00000123983', 'ENSG00000124145', 'ENSG00000124181', 'ENSG00000124795', 'ENSG00000125618', 'ENSG00000125952', 'ENSG00000126012', 'ENSG00000126524', 'ENSG00000126746', 'ENSG00000126752', 'ENSG00000126777', 'ENSG00000126778', 'ENSG00000126883', 'ENSG00000126934', 'ENSG00000127152', 'ENSG00000127329', 'ENSG00000127616', 'ENSG00000127946', 'ENSG00000128052', 'ENSG00000128513', 'ENSG00000128602', 'ENSG00000128713', 'ENSG00000128714', 'ENSG00000129152', 'ENSG00000129204', 'ENSG00000129514', 'ENSG00000129993', 'ENSG00000130382', 'ENSG00000130396', 'ENSG00000130779', 'ENSG00000130844', 'ENSG00000131023', 'ENSG00000131653', 'ENSG00000131759', 'ENSG00000132002', 'ENSG00000132155', 'ENSG00000132170', 'ENSG00000132475', 'ENSG00000132781', 'ENSG00000133124', 'ENSG00000133392', 'ENSG00000133639', 'ENSG00000133703', 'ENSG00000133818', 'ENSG00000133895', 'ENSG00000134086', 'ENSG00000134250', 'ENSG00000134323', 'ENSG00000134352', 'ENSG00000134371', 'ENSG00000134574', 'ENSG00000134853', 'ENSG00000134899', 'ENSG00000134982', 'ENSG00000135100', 'ENSG00000135111', 'ENSG00000135363', 'ENSG00000135446', 'ENSG00000135503', 'ENSG00000135679', 'ENSG00000135903', 'ENSG00000135956', 'ENSG00000136238', 'ENSG00000136352', 'ENSG00000136492', 'ENSG00000136754', 'ENSG00000136826', 'ENSG00000136936', 'ENSG00000136997', 'ENSG00000137193', 'ENSG00000137265', 'ENSG00000137309', 'ENSG00000137497', 'ENSG00000137812', 'ENSG00000138081', 'ENSG00000138336', 'ENSG00000138363', 'ENSG00000138376', 'ENSG00000138413', 'ENSG00000138592', 'ENSG00000138698', 'ENSG00000138795', 'ENSG00000139083', 'ENSG00000139163', 'ENSG00000139219', 'ENSG00000139263', 'ENSG00000139618', 'ENSG00000139687', 'ENSG00000140262', 'ENSG00000140396', 'ENSG00000140464', 'ENSG00000140538', 'ENSG00000140577', 'ENSG00000140836', 'ENSG00000140937', 'ENSG00000141027', 'ENSG00000141367', 'ENSG00000141380', 'ENSG00000141510', 'ENSG00000141646', 'ENSG00000141736', 'ENSG00000141867', 'ENSG00000141985', 'ENSG00000142208', 'ENSG00000142273', 'ENSG00000142611', 'ENSG00000142867', 'ENSG00000143252', 'ENSG00000143294', 'ENSG00000143322', 'ENSG00000143437', 'ENSG00000143549', 'ENSG00000143924', 'ENSG00000144218', 'ENSG00000144476', 'ENSG00000144554', 'ENSG00000145012', 'ENSG00000145216', 'ENSG00000145675', 'ENSG00000145819', 'ENSG00000146232', 'ENSG00000146374', 'ENSG00000146648', 'ENSG00000147050', 'ENSG00000147065', 'ENSG00000147140', 'ENSG00000147257', 'ENSG00000147403', 'ENSG00000147548', 'ENSG00000147655', 'ENSG00000147862', 'ENSG00000147889', 'ENSG00000148053', 'ENSG00000148400', 'ENSG00000148737', 'ENSG00000149311', 'ENSG00000149948', 'ENSG00000150457', 'ENSG00000150907', 'ENSG00000151348', 'ENSG00000151702', 'ENSG00000152217', 'ENSG00000152894', 'ENSG00000153201', 'ENSG00000153944', 'ENSG00000154767', 'ENSG00000154803', 'ENSG00000156052', 'ENSG00000156076', 'ENSG00000156531', 'ENSG00000156650', 'ENSG00000156970', 'ENSG00000156976', 'ENSG00000157168', 'ENSG00000157388', 'ENSG00000157404', 'ENSG00000157554', 'ENSG00000157613', 'ENSG00000157764', 'ENSG00000157765', 'ENSG00000157873', 'ENSG00000158169', 'ENSG00000158711', 'ENSG00000158715', 'ENSG00000159216', 'ENSG00000160007', 'ENSG00000160201', 'ENSG00000160789', 'ENSG00000160867', 'ENSG00000160957', 'ENSG00000161405', 'ENSG00000161547', 'ENSG00000162367', 'ENSG00000162434', 'ENSG00000162613', 'ENSG00000162733', 'ENSG00000162775', 'ENSG00000162924', 'ENSG00000163026', 'ENSG00000163041', 'ENSG00000163161', 'ENSG00000163399', 'ENSG00000163497', 'ENSG00000163513', 'ENSG00000163518', 'ENSG00000163629', 'ENSG00000163902', 'ENSG00000163930', 'ENSG00000163939', 'ENSG00000164330', 'ENSG00000164362', 'ENSG00000164438', 'ENSG00000164683', 'ENSG00000164754', 'ENSG00000164985', 'ENSG00000165025', 'ENSG00000165392', 'ENSG00000165556', 'ENSG00000165671', 'ENSG00000165699', 'ENSG00000165731', 'ENSG00000166407', 'ENSG00000166710', 'ENSG00000166886', 'ENSG00000166888', 'ENSG00000166949', 'ENSG00000167258', 'ENSG00000167460', 'ENSG00000167548', 'ENSG00000167751', 'ENSG00000167985', 'ENSG00000168036', 'ENSG00000168092', 'ENSG00000168172', 'ENSG00000168421', 'ENSG00000168610', 'ENSG00000168646', 'ENSG00000168685', 'ENSG00000168702', 'ENSG00000168769', 'ENSG00000169032', 'ENSG00000169083', 'ENSG00000169184', 'ENSG00000169249', 'ENSG00000169696', 'ENSG00000169714', 'ENSG00000169925', 'ENSG00000170759', 'ENSG00000170791', 'ENSG00000170836', 'ENSG00000171094', 'ENSG00000171302', 'ENSG00000171456', 'ENSG00000171723', 'ENSG00000171735', 'ENSG00000171791', 'ENSG00000171843', 'ENSG00000171862', 'ENSG00000172175', 'ENSG00000172493', 'ENSG00000172936', 'ENSG00000173757', 'ENSG00000173821', 'ENSG00000174775', 'ENSG00000175054', 'ENSG00000175197', 'ENSG00000175387', 'ENSG00000175595', 'ENSG00000175643', 'ENSG00000175832', 'ENSG00000177084', 'ENSG00000177565', 'ENSG00000177606', 'ENSG00000178053', 'ENSG00000178104', 'ENSG00000178105', 'ENSG00000178568', 'ENSG00000178573', 'ENSG00000178691', 'ENSG00000179094', 'ENSG00000179218', 'ENSG00000179295', 'ENSG00000179348', 'ENSG00000179583', 'ENSG00000179750', 'ENSG00000180644', 'ENSG00000181163', 'ENSG00000181449', 'ENSG00000181555', 'ENSG00000181690', 'ENSG00000182054', 'ENSG00000182158', 'ENSG00000182162', 'ENSG00000182185', 'ENSG00000182197', 'ENSG00000182511', 'ENSG00000182866', 'ENSG00000182872', 'ENSG00000182944', 'ENSG00000183161', 'ENSG00000183337', 'ENSG00000183454', 'ENSG00000183508', 'ENSG00000183765', 'ENSG00000183770', 'ENSG00000184012', 'ENSG00000184384', 'ENSG00000184402', 'ENSG00000184481', 'ENSG00000184507', 'ENSG00000184634', 'ENSG00000184675', 'ENSG00000184937', 'ENSG00000185338', 'ENSG00000185499', 'ENSG00000185630', 'ENSG00000185811', 'ENSG00000185920', 'ENSG00000186051', 'ENSG00000186174', 'ENSG00000186575', 'ENSG00000186716', 'ENSG00000187098', 'ENSG00000187735', 'ENSG00000187741', 'ENSG00000188199', 'ENSG00000189079', 'ENSG00000189283', 'ENSG00000196090', 'ENSG00000196092', 'ENSG00000196159', 'ENSG00000196367', 'ENSG00000196498', 'ENSG00000196588', 'ENSG00000196712', 'ENSG00000196914', 'ENSG00000197122', 'ENSG00000197157', 'ENSG00000197299', 'ENSG00000197323', 'ENSG00000197535', 'ENSG00000197646', 'ENSG00000198286', 'ENSG00000198400', 'ENSG00000198625', 'ENSG00000198793', 'ENSG00000198795', 'ENSG00000198900', 'ENSG00000204103', 'ENSG00000204209', 'ENSG00000204370', 'ENSG00000204531', 'ENSG00000204713', 'ENSG00000204843', 'ENSG00000205755', 'ENSG00000205927', 'ENSG00000206503', 'ENSG00000213066', 'ENSG00000213190', 'ENSG00000213281', 'ENSG00000214562', 'ENSG00000214827', 'ENSG00000215301', 'ENSG00000221829', 'ENSG00000241476', 'ENSG00000244405', 'ENSG00000245848', 'ENSG00000257923', 'ENSG00000266412', 'ENSG00000268009', 'ENSG00000270647', 'ENSG00000274267', 'ENSG00000275023', 'ENSG00000276180')
# CGC tier 2:
# ('ENSG00000029363', 'ENSG00000029534', 'ENSG00000033030', 'ENSG00000036257', 'ENSG00000040731', 'ENSG00000041982', 'ENSG00000044115', 'ENSG00000044524', 'ENSG00000048471', 'ENSG00000049540', 'ENSG00000054118', 'ENSG00000064933', 'ENSG00000065057', 'ENSG00000066032', 'ENSG00000066117', 'ENSG00000066279', 'ENSG00000070756', 'ENSG00000072274', 'ENSG00000072501', 'ENSG00000073792', 'ENSG00000074266', 'ENSG00000074964', 'ENSG00000078061', 'ENSG00000078177', 'ENSG00000079112', 'ENSG00000082701', 'ENSG00000090659', 'ENSG00000100852', 'ENSG00000101343', 'ENSG00000104660', 'ENSG00000104728', 'ENSG00000105619', 'ENSG00000107864', 'ENSG00000107929', 'ENSG00000109220', 'ENSG00000110844', 'ENSG00000111087', 'ENSG00000111275', 'ENSG00000111679', 'ENSG00000112175', 'ENSG00000112237', 'ENSG00000113384', 'ENSG00000113387', 'ENSG00000115760', 'ENSG00000116731', 'ENSG00000117020', 'ENSG00000117318', 'ENSG00000118007', 'ENSG00000118515', 'ENSG00000121289', 'ENSG00000122642', 'ENSG00000122778', 'ENSG00000124486', 'ENSG00000124762', 'ENSG00000125285', 'ENSG00000125354', 'ENSG00000126353', 'ENSG00000126453', 'ENSG00000127083', 'ENSG00000127314', 'ENSG00000127914', 'ENSG00000128191', 'ENSG00000128487', 'ENSG00000128944', 'ENSG00000130675', 'ENSG00000132906', 'ENSG00000135333', 'ENSG00000135605', 'ENSG00000136014', 'ENSG00000136167', 'ENSG00000136504', 'ENSG00000138115', 'ENSG00000138448', 'ENSG00000139718', 'ENSG00000140521', 'ENSG00000141968', 'ENSG00000143379', 'ENSG00000143556', 'ENSG00000143970', 'ENSG00000145113', 'ENSG00000147130', 'ENSG00000147724', 'ENSG00000148516', 'ENSG00000148584', 'ENSG00000151532', 'ENSG00000152207', 'ENSG00000152942', 'ENSG00000153165', 'ENSG00000153707', 'ENSG00000153814', 'ENSG00000157933', 'ENSG00000159388', 'ENSG00000159784', 'ENSG00000160271', 'ENSG00000163435', 'ENSG00000163520', 'ENSG00000163655', 'ENSG00000164305', 'ENSG00000164398', 'ENSG00000164796', 'ENSG00000164919', 'ENSG00000165238', 'ENSG00000165323', 'ENSG00000165409', 'ENSG00000166501', 'ENSG00000168040', 'ENSG00000168411', 'ENSG00000168496', 'ENSG00000168542', 'ENSG00000169564', 'ENSG00000169862', 'ENSG00000170234', 'ENSG00000170430', 'ENSG00000170577', 'ENSG00000171310', 'ENSG00000172409', 'ENSG00000172915', 'ENSG00000173575', 'ENSG00000173674', 'ENSG00000174469', 'ENSG00000175329', 'ENSG00000176302', 'ENSG00000176571', 'ENSG00000178562', 'ENSG00000179362', 'ENSG00000179399', 'ENSG00000180611', 'ENSG00000181143', 'ENSG00000181222', 'ENSG00000182578', 'ENSG00000182901', 'ENSG00000183579', 'ENSG00000183722', 'ENSG00000183742', 'ENSG00000183813', 'ENSG00000184304', 'ENSG00000184640', 'ENSG00000184702', 'ENSG00000184956', 'ENSG00000185008', 'ENSG00000185177', 'ENSG00000187164', 'ENSG00000187239', 'ENSG00000187323', 'ENSG00000196220', 'ENSG00000196531', 'ENSG00000196924', 'ENSG00000197013', 'ENSG00000197880', 'ENSG00000198053', 'ENSG00000198173', 'ENSG00000198354', 'ENSG00000198561', 'ENSG00000198604', 'ENSG00000198822', 'ENSG00000203734', 'ENSG00000205542', 'ENSG00000213672', 'ENSG00000251562', 'ENSG00000254087', 'ENSG00000257103', 'ENSG00000261652')
# Note: 5 ENSGs from each set are missing in the VEP output and hence in snps_cosmic (I suppose they're not in Ensembl 108).

# Cell Lines Project (probably not useful since these will have accumulated so many additional mutations)
# coding
DownloadFile(f"CellLinesProject_GenomeScreensMutant_VcfNormal_v{rel}_GRCh38.tar", url=f"https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch38/cell_lines/v{rel}/VCF/CellLinesProject_GenomeScreensMutant_VcfNormal_v{rel}_GRCh38.tar&bucket=downloads")
# noncoding
DownloadFile(f"CellLinesProject_NonCodingVariants_VcfNormal_v{rel}_GRCh38.tar", url=f"https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch38/cell_lines/v{rel}/VCF/CellLinesProject_NonCodingVariants_VcfNormal_v{rel}_GRCh38.tar&bucket=downloads")

# TSV versions
DownloadFile(f"CellLinesProject_GenomeScreensMutant_Tsv_v{rel}_GRCh38.tar", url=f"https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch38/cell_lines/v{rel}/CellLinesProject_GenomeScreensMutant_Tsv_v{rel}_GRCh38.tar&bucket=downloads")
DownloadFile(f"CellLinesProject_NonCodingVariants_Tsv", url=f"https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch38/cell_lines/v{rel}/CellLinesProject_NonCodingVariants_Tsv_v{rel}_GRCh38.tar&bucket=downloads")


# Show directory
Run("Show directory", "ls -lah input")

print("Done!")
