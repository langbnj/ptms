#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

# our $superloudmysql = 1;

$comparahomology = "comparahomology";

our $usage = "$0 [-hc]\n\n -hc: Only keep high-confidence orthologs\n\nExample: $0 -hc";
args(0);



# start

Clear($comparahomology);

starttime();
@infiles = `ls input/hsapiens_gene_ensembl__homolog_*__dm.txt -1`;
%mart = ();
startme("Getting one2one orthologs from Ensembl Mart files for human proteins in other species and updating table '$comparahomology'", 0, scalar(@infiles));
foreach $infile (@infiles)
{
	chomp($infile);
	die if (!-e $infile);
	
	open(IN, $infile) or die;
	# state("Reading '$infile'");
	while (<IN>)
	{
		# next if !/one2one/;
		next if /^\\N\t/;

		chomp;
		
		# New (release 106) (15 columns)
		# 100.00	MT	564	4951	3914	ortholog_one2one	Euarchontoglires	57.3913	75	ENSMUSG00000064345	ENSMUSP00000080992	ENSP00000355046	mt-Nd2	57.0605	1
		# 99.35	5	13093	150493794	150446095	ortholog_one2one	Euarchontoglires	58.5161	100	ENSMUSG00000041147	ENSMUSP00000038576	ENSP00000369497	Brca2	56.9924	1
		# # C. elegans:
		# # V	29981	8606128	8605088	ortholog_one2one	Bilateria	60.8295	WBGene00018290	F41E6.9.1	ENSP00000223500	vps-60	60.274	1
		# Different number of columns! Need to account for that (13 instead of 15)
		# # It's missing the goc_score. It still seems to have 0/1 at the end, though, which I thought was high-confidence...
		# # ~/update/comparanopara >> head -n1 input/hsapiens_gene_ensembl__homolog_*__dm.txt | g -v '^==>' | g -v '^$' | perl -ne '$i=1; $i++ while /\t/g; print "$i\n"' | suq
		# >> Some species have 13, 14, or most commonly, 15 columns
		# head -n1 input/*
		# Perhaps theres a better version of this data somewhere? A better format?
		# Actually its a manageable format, but I need to account for the different column counts below!
		
		
		
		@a = split(/\t/);

		# Input fields (from Mart files)
		$ensp1 = '';
		$ensp2 = '';
		$ensg2 = '';
		$symbol2 = '';
		$homology = '';
		$hc = '';
		$clade = '';
		$target_percent_id = '';
		$query_percent_id = '';
		$gene_order_score = '';
		$whole_genome_align_score = '';
		$target_chr = '';
		$target_chr_start = '';
		$target_chr_stop = '';
		$target_chr_strand = '';

		# Additional fields (from MySQL)
		$ensg1 = '';
		$enst1 = '';
		$species1 = '';
		$fullspecies1 = '';
		$symbol1 = '';
		$enst2 = '';
		$species2 = '';
		$fullspecies2 = '';
		$acc1 = '';
		$name1 = '';
		$acc2 = '';
		$name2 = '';

		if (scalar(@a) == 13)
		{
			# state(" >> 13 columns");
			# head -vn1 input/hsapiens_gene_ensembl__homolog_*__dm.txt | perl -ne '@a = split(/\t/); print if ((scalar(@a) == 1) or (scalar(@a) == 13))' | g -v "^$" | g -B1 '^[^=]' | g -v '\-\-' | g -v '^[^=]' | perl -npe 's/^==> input\///; s/ <==$//'
			# Invertebrates
			# hsapiens_gene_ensembl__homolog_celegans__dm.txt
			# hsapiens_gene_ensembl__homolog_dmelanogaster__dm.txt
			# hsapiens_gene_ensembl__homolog_scerevisiae__dm.txt

			# g '^(([^\t]+)\t){12}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_*__dm.txt
			# g '^(([^\t]+)\t){12}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_*__dm.txt | perl -npe 's/^input\/hsapiens_gene_ensembl__homolog_\w+\.txt://;'
			# g '^(([^\t]+)\t){12}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_*__dm.txt | cut -f 13 | suq # >> $a[12] is always 0 or 1 and appears to be HC

			# C. elegans (13 columns)
			# I	68600	2489787	2485861	ortholog_one2one	Bilateria	39.4495	WBGene00001232	Y74C10AR.1a.1	ENSP00000503706	eif-3.I	29.4521	1
			# chromo	id	start	stop	homology	clade	target_id	homogene	homoensp	ensp	homosymbol	query_id	hc
			
			# S. cerevisiae (13 columns)
			# XIII	68600	558524	557481	ortholog_one2one	Opisthokonta	40.634	YMR146C	YMR146C	ENSP00000503706	TIF34	32.1918	1
			# chromo	id	start	stop	homology	clade	target_id	homogene	homoensp	ensp	homosymbol	query_id	hc
			
			$gene_order_score = '';
			$whole_genome_align_score = '';

			$target_chr = $a[0];
			# $?id? = $a[1];	# Some kind of running number (even exists in otherwise blank rows)
			$target_chr_start = $a[2];	# Get flipped for - strand, as opposed to elsewhere in Ensembl. Strand is not explicitly given.
			$target_chr_stop = $a[3];  # Get flipped for - strand, as opposed to elsewhere in Ensembl. Strand is not explicitly given.
			$homology = $a[4];
			$clade = $a[5];
			$target_percent_id = $a[6];
			$ensg2 = $a[7];
			$ensp2 = $a[8];
			$ensp1 = $a[9];
			$symbol2 = $a[10];
			$query_percent_id = $a[11];
			$hc = $a[12];			# High confidence (https://useast.ensembl.org/info/genome/compara/Ortholog_qc_manual.html)
			# "Orthologue pairs may be tagged as high confidence if they meet certain thresholds for these two measures [Gene Order Conservation score & Whole Genome Alignment score] as well as identity. The thresholds we use depend on the most-recent common ancestor of the species pair, according to the table below. The orthologue pair must satisfy all the criteria to be tagged as high confidence."

			# Clades	Min. GOC score	Min. WGA score	Min. %identity
			# Apes, Murinae	75	75	80
			# Mammalia, Aves, Percomorpha	75	75	50
			# Euteleostomi	50	50	25
			# Others	No threshold used	No threshold used	25	<<<
		}
		elsif ((scalar(@a) == 14) and ($a[0] =~ /^\d+\.\d+$/))
		# These have 14 columns, and they always have numbers in the first column, as well as in column $a[7].
		# They have WGA scores, but no GOC scores.
		{
			# state(" >> 14 columns, first numeric");
			# head -vn1 input/hsapiens_gene_ensembl__homolog_*__dm.txt | perl -ne '@a = split(/\t/); print if ((scalar(@a) == 1) or (scalar(@a) == 14))' | g -v "^$" | g -B1 '^[^=]' | g -v '\-\-' | g -v '^[^=]' | perl -npe 's/^==> input\///; s/ <==$//'
			# More distant relatives (non-mammals) (includes marsupial: Tasmanian devil)

			# These have 14 columns, and they always have numbers in the first column, as well as in column $a[7]:
			# cat input/hsapiens_gene_ensembl__homolog_cintestinalis__dm.txt input/hsapiens_gene_ensembl__homolog_cintestinalis__dm.txt input/hsapiens_gene_ensembl__homolog_csavignyi__dm.txt input/hsapiens_gene_ensembl__homolog_eburgeri__dm.txt input/hsapiens_gene_ensembl__homolog_pmarinus__dm.txt | cut -f 1 | perl -npe 's/^input\/hsapiens_gene_ensembl__homolog_\w+\.txt://;' | g -v '\\N' | g -v '^\d+\.\d+$' | suq
			# g '^(([^\t]+)\t){13}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_cintestinalis__dm.txt | cut -f 1 | g -v '\\N' | g -v '^\d+\.\d+$' | suq 	# a7 is WGA score
			# g '^(([^\t]+)\t){13}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_csavignyi__dm.txt | cut -f 1 | g -v '\\N' | g -v '^\d+\.\d+$' | suq 		# a7 is WGA score
			# g '^(([^\t]+)\t){13}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_eburgeri__dm.txt | cut -f 1 | g -v '\\N' | g -v '^\d+\.\d+$' | suq 		# a7 is WGA score
			# g '^(([^\t]+)\t){13}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_pmarinus__dm.txt | cut -f 1 | g -v '\\N' | g -v '^\d+\.\d+$' | suq 		# a7 is WGA score
			
 			# g '^(([^\t]+)\t){13}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_*__dm.txt | cut -f 8 | suq
 			# hc:
			# g '^(([^\t]+)\t){13}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_*__dm.txt | cut -f 14 | suq
			

			# Ciona intestinalis
			# cat input/hsapiens_gene_ensembl__homolog_cintestinalis__dm.txt | g -v '\\N' | head
			# 49.30	3	70635	2990193	2985217	ortholog_one2many	Chordata	14.5414	ENSCING00000005631	ENSCINP00000011665	ENSP00000478910	tektinA1	50.3876	0
			# 0.00	HT000823.1	73420	62896	51828	ortholog_many2many	Bilateria	5.8296	ENSCING00000009075	ENSCINP00000018418	ENSP00000318502	pbx	14.0541	0
			# 0.00	7	73588	2413385	2395065	ortholog_many2many	Bilateria	9.14286	ENSCING00000007901	ENSCINP00000016177	ENSP00000372547	soxC	23.5294	0
			# >>wga_score<<	chromo	id	start	stop	homology	clade	target_id	homogene	homoensp	ensp	homosymbol	query_id	hc

			$gene_order_score = '';

			$whole_genome_align_score = $a[0];
			$target_chr = $a[1];
			# $?id? = $a[2];	# Some kind of running number (even exists in otherwise blank rows)
			$target_chr_start = $a[3];	# Get flipped for - strand, as opposed to elsewhere in Ensembl. Strand is not explicitly given.
			$target_chr_stop = $a[4];  # Get flipped for - strand, as opposed to elsewhere in Ensembl. Strand is not explicitly given.
			$homology = $a[5];
			$clade = $a[6];
			$target_percent_id = $a[7];
			$ensg2 = $a[8];
			$ensp2 = $a[9];
			$ensp1 = $a[10];
			$symbol2 = $a[11];
			$query_percent_id = $a[12];
			$hc = $a[13];			# High confidence (https://useast.ensembl.org/info/genome/compara/Ortholog_qc_manual.html)
			# "Orthologue pairs may be tagged as high confidence if they meet certain thresholds for these two measures [Gene Order Conservation score & Whole Genome Alignment score] as well as identity. The thresholds we use depend on the most-recent common ancestor of the species pair, according to the table below. The orthologue pair must satisfy all the criteria to be tagged as high confidence."

			# Clades	Min. GOC score	Min. WGA score	Min. %identity
			# Apes, Murinae	75	75	80
			# Mammalia, Aves, Percomorpha	75	75	50
			# Euteleostomi	50	50	25		<<< ? This should really have WGA scores, then, though
			# Others	No threshold used	No threshold used	25
		}
		elsif ((scalar(@a) == 14) and ($a[0] !~ /^\d+\.\d+$/))
		# These have 14 columns, and they never have numbers (which would be WGA scores) in the first column (they have chromosome names instead), and they have 0/25/50/75/100 (GOC scores) in column $a[7]:
		# They have GOC scores, but no WGA scores.
		{
			# state(" >> 14 columns, first is a chromosome name");
			# head -vn1 input/hsapiens_gene_ensembl__homolog_*__dm.txt | perl -ne '@a = split(/\t/); print if ((scalar(@a) == 1) or (scalar(@a) == 14))' | g -v "^$" | g -B1 '^[^=]' | g -v '\-\-' | g -v '^[^=]' | perl -npe 's/^==> input\///; s/ <==$//'
			# More distant relatives (non-mammals) (includes marsupial: Tasmanian devil)
			
			# These have 14 columns, but they never have numbers (WGA scores) in the first column (they have chromosome names instead), and they have 0/25/50/75/100 (GOC scores) in column $a[7]:
			# cat input/hsapiens_gene_ensembl__homolog_acarolinensis__dm.txt input/hsapiens_gene_ensembl__homolog_amelanoleuca__dm.txt input/hsapiens_gene_ensembl__homolog_atestudineus__dm.txt input/hsapiens_gene_ensembl__homolog_cccarpio__dm.txt input/hsapiens_gene_ensembl__homolog_dlabrax__dm.txt input/hsapiens_gene_ensembl__homolog_falbicollis__dm.txt input/hsapiens_gene_ensembl__homolog_ggallus__dm.txt input/hsapiens_gene_ensembl__homolog_gmorhua__dm.txt input/hsapiens_gene_ensembl__homolog_marmatus__dm.txt input/hsapiens_gene_ensembl__homolog_mgallopavo__dm.txt input/hsapiens_gene_ensembl__homolog_omykiss__dm.txt input/hsapiens_gene_ensembl__homolog_sharrisii__dm.txt input/hsapiens_gene_ensembl__homolog_smaximus__dm.txt input/hsapiens_gene_ensembl__homolog_ssalar__dm.txt input/hsapiens_gene_ensembl__homolog_xtropicalis__dm.txt | cut -f 1 | perl -npe 's/^input\/hsapiens_gene_ensembl__homolog_\w+\.txt://;' | g -v '\\N' | g '^\d+\.\d+$' | suq
			# g '^(([^\t]+)\t){13}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_acarolinensis__dm.txt | cut -f 1 | g -v '\\N' | g -v '^\d+\.\d+$' | suq 	# a7 is GOC score
			# g '^(([^\t]+)\t){13}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_amelanoleuca__dm.txt | cut -f 1 | g -v '\\N' | g -v '^\d+\.\d+$' | suq 		# a7 is GOC score
			# g '^(([^\t]+)\t){13}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_atestudineus__dm.txt | cut -f 1 | g -v '\\N' | g -v '^\d+\.\d+$' | suq 		# a7 is GOC score
			# g '^(([^\t]+)\t){13}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_cccarpio__dm.txt | cut -f 1 | g -v '\\N' | g -v '^\d+\.\d+$' | suq 			# a7 is GOC score
			# g '^(([^\t]+)\t){13}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_dlabrax__dm.txt | cut -f 1 | g -v '\\N' | g -v '^\d+\.\d+$' | suq 			# a7 is GOC score
			# g '^(([^\t]+)\t){13}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_falbicollis__dm.txt | cut -f 1 | g -v '\\N' | g -v '^\d+\.\d+$' | suq 		# a7 is GOC score
			# g '^(([^\t]+)\t){13}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_ggallus__dm.txt | cut -f 1 | g -v '\\N' | g -v '^\d+\.\d+$' | suq 			# a7 is GOC score
			# g '^(([^\t]+)\t){13}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_gmorhua__dm.txt | cut -f 1 | g -v '\\N' | g -v '^\d+\.\d+$' | suq 			# a7 is GOC score
			# g '^(([^\t]+)\t){13}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_marmatus__dm.txt | cut -f 1 | g -v '\\N' | g -v '^\d+\.\d+$' | suq 			# a7 is GOC score
			# g '^(([^\t]+)\t){13}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_mgallopavo__dm.txt | cut -f 1 | g -v '\\N' | g -v '^\d+\.\d+$' | suq 		# a7 is GOC score
			# g '^(([^\t]+)\t){13}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_omykiss__dm.txt | cut -f 1 | g -v '\\N' | g -v '^\d+\.\d+$' | suq 			# a7 is GOC score
			# g '^(([^\t]+)\t){13}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_sharrisii__dm.txt | cut -f 1 | g -v '\\N' | g -v '^\d+\.\d+$' | suq 		# a7 is GOC score
			# g '^(([^\t]+)\t){13}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_smaximus__dm.txt | cut -f 1 | g -v '\\N' | g -v '^\d+\.\d+$' | suq 			# a7 is GOC score
			# g '^(([^\t]+)\t){13}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_ssalar__dm.txt | cut -f 1 | g -v '\\N' | g -v '^\d+\.\d+$' | suq 			# a7 is GOC score
			# g '^(([^\t]+)\t){13}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_xtropicalis__dm.txt | cut -f 1 | g -v '\\N' | g -v '^\d+\.\d+$' | suq 		# a7 is GOC score

			# MT	559	3633	2673	ortholog_one2one	Amniota	67.1875	50	ENSACAG00000028215	ENSACAP00000021683	ENSP00000354687	MT-ND1	67.6101	1
			# MT	564	4876	3841	ortholog_one2one	Vertebrata	46.6667	75	ENSACAG00000028219	ENSACAP00000021684	ENSP00000355046	MT-ND2	46.3977	1
			# MT	572	6802	5258	ortholog_one2one	Euteleostomi	86.965	100	ENSACAG00000028225	ENSACAP00000021685	ENSP00000354499	MT-CO1	87.1345	1



 			# g '^(([^\t]+)\t){13}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_*__dm.txt | cut -f 8 | suq
 			# hc:
			# g '^(([^\t]+)\t){13}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_*__dm.txt | cut -f 14 | suq
			

			# Anolis carolinensis
			# cat input/hsapiens_gene_ensembl__homolog_acarolinensis__dm.txt | g -v '\\N' | head
			# MT	559	3633	2673	ortholog_one2one	Amniota	67.1875	50	ENSACAG00000028215	ENSACAP00000021683	ENSP00000354687	MT-ND1	67.6101	1
			# MT	564	4876	3841	ortholog_one2one	Vertebrata	46.6667	75	ENSACAG00000028219	ENSACAP00000021684	ENSP00000355046	MT-ND2	46.3977	1
			# MT	572	6802	5258	ortholog_one2one	Euteleostomi	86.965	100	ENSACAG00000028225	ENSACAP00000021685	ENSP00000354499	MT-CO1	87.1345	1
			# chromo	id	start	stop	homology	clade	target_id	>>goc_score<<	homogene	homoensp	ensp	homosymbol	query_id	hc
			
			$whole_genome_align_score = '';

			$target_chr = $a[0];
			# $?id? = $a[1];	# Some kind of running number (even exists in otherwise blank rows)
			$target_chr_start = $a[2];	# Get flipped for - strand, as opposed to elsewhere in Ensembl. Strand is not explicitly given.
			$target_chr_stop = $a[3];  # Get flipped for - strand, as opposed to elsewhere in Ensembl. Strand is not explicitly given.
			$homology = $a[4];
			$clade = $a[5];
			$target_percent_id = $a[6];
			$gene_order_score = $a[7];	# Gene order conservation score (https://useast.ensembl.org/info/genome/compara/Ortholog_qc_manual.html)
			$ensg2 = $a[8];
			$ensp2 = $a[9];
			$ensp1 = $a[10];
			$symbol2 = $a[11];
			$query_percent_id = $a[12];
			$hc = $a[13];			# High confidence (https://useast.ensembl.org/info/genome/compara/Ortholog_qc_manual.html)
			# "Orthologue pairs may be tagged as high confidence if they meet certain thresholds for these two measures [Gene Order Conservation score & Whole Genome Alignment score] as well as identity. The thresholds we use depend on the most-recent common ancestor of the species pair, according to the table below. The orthologue pair must satisfy all the criteria to be tagged as high confidence."

			# Clades	Min. GOC score	Min. WGA score	Min. %identity
			# Apes, Murinae	75	75	80
			# Mammalia, Aves, Percomorpha	75	75	50
			# Euteleostomi	50	50	25		<<< ? This should really have WGA scores, then, though
			# Others	No threshold used	No threshold used	25
		}
		elsif (scalar(@a) == 15)
		# These have both GOC scores and WGA scores.
		{
			# state(" >> 15 columns");
			# head -vn1 input/hsapiens_gene_ensembl__homolog_*__dm.txt | perl -ne '@a = split(/\t/); print if ((scalar(@a) == 1) or (scalar(@a) == 15))' | g -v "^$" | g -B1 '^[^=]' | g -v '\-\-' | g -v '^[^=]' | perl -npe 's/^==> input\///; s/ <==$//'
			# Mammals etc.
			# e.g.
			# hsapiens_gene_ensembl__homolog_abrachyrhynchus__dm.txt
			# hsapiens_gene_ensembl__homolog_acalliptera__dm.txt
			# hsapiens_gene_ensembl__homolog_acchrysaetos__dm.txt
			# hsapiens_gene_ensembl__homolog_acitrinellus__dm.txt

			# GOC score:
			# g '^(([^\t]+)\t){14}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_*__dm.txt | cut -f 9 | suq
			# 0, 25, 50, 75, 100
			
			# WGA score:
			# g '^(([^\t]+)\t){14}([^\t]+)$' input/hsapiens_gene_ensembl__homolog_*__dm.txt | perl -npe 's/^input\/hsapiens_gene_ensembl__homolog_\w+\.txt://;' | cut -f 1 | g -v "^\\\\" | summary
			# nan	69.428429709263	91.59	100

			# Mus musculus (15 columns)
			# 100.00	4	68597	149392150	149260776	ortholog_one2one	Eutheria	93.5741	100	ENSMUSG00000063077	ENSMUSP00000158941	ENSP00000502065	Kif1b	97.0264	1
			# >>wga_score<<	chromo	id	start	stop	homology	clade	target_id	>>goc_score<<	homogene	homoensp	ensp	homosymbol	query_id	hc

			# A citrinellus
			# cat input/hsapiens_gene_ensembl__homolog_acitrinellus__dm.txt | g -v '\\N' | head
			# 0.00	CCOE01001026.1	70299	527693	516296	ortholog_many2many	Euteleostomi	67.4089	0	ENSACIG00000022724	ENSACIP00000029350	ENSP00000481127	prodhb	55.5	0
			# 88.79	CCOE01000490.1	70299	204272	194868	ortholog_many2many	Euteleostomi	67.2968	0	ENSACIG00000015775	ENSACIP00000020316	ENSP00000481127	prodha	59.3333	1
			# 86.63	CCOE01000977.1	70328	2045083	2042383	ortholog_one2many	Euteleostomi	59.0698	0	ENSACIG00000020714	ENSACIP00000026755	ENSP00000482514	dgcr6	57.7273	1
			# >>wga_score<<	chromo	id	start	stop	homology	clade	target_id	>>goc_score<<	homogene	homoensp	ensp	homosymbol	query_id	hc

			$whole_genome_align_score = $a[0];
			$target_chr = $a[1];
			# $?id? = $a[2];	# Some kind of running number (even exists in otherwise blank rows)
			$target_chr_start = $a[3];	# Get flipped for - strand, as opposed to elsewhere in Ensembl. Strand is not explicitly given.
			$target_chr_stop = $a[4];  # Get flipped for - strand, as opposed to elsewhere in Ensembl. Strand is not explicitly given.
			$homology = $a[5];
			$clade = $a[6];
			$target_percent_id = $a[7];
			$gene_order_score = $a[8];	# Gene order conservation score (https://useast.ensembl.org/info/genome/compara/Ortholog_qc_manual.html)
			$ensg2 = $a[9];
			$ensp2 = $a[10];
			$ensp1 = $a[11];
			$symbol2 = $a[12];
			$query_percent_id = $a[13];
			$hc = $a[14];			# High confidence (https://useast.ensembl.org/info/genome/compara/Ortholog_qc_manual.html)
			# "Orthologue pairs may be tagged as high confidence if they meet certain thresholds for these two measures [Gene Order Conservation score & Whole Genome Alignment score] as well as identity. The thresholds we use depend on the most-recent common ancestor of the species pair, according to the table below. The orthologue pair must satisfy all the criteria to be tagged as high confidence."
		}
		else
		{
			die("Error: Expected to have either 13, 14, or 15 columns, but got ".scalar(@a)." in line:\n\n$_\n\n");
		}
		
		# Get strand
		# Coordinates get flipped for - strand, as opposed to elsewhere in Ensembl. Strand is not explicitly given.
		die if ($target_chr_start == $target_chr_stop);
		if ($target_chr_start < $target_chr_stop)
		{
			$target_chr_strand = '+';
		}
		else
		{
			$target_chr_strand = '-';
		}


		# Get extra fields

		# ensg1
		# enst1
		# species1
		# fullspecies1
		# symbol1
		$query = Query("SELECT ensg, enst, species, fullspecies, symbol FROM ensembl WHERE ensp='$ensp1'");
		if (Numrows($query) == 0)
		{
			# # Skipping this for speed (these ensembl_gff3 queries were quite slow, and I don't really need non-coding transcripts)
			# Try ENSP as ENST (these are non-coding transcripts that I'm only keeping in the comparahomology table, not in any of the comparanopara output)
			# via ensembl (won't work since these are non-coding, and ensembl only contains protein_coding)
			# $query = Query("SELECT ensg, ensp, species, fullspecies, symbol FROM ensembl WHERE enst='$ensp1'");
			# via ensembl_gff3 tables
			$query = Query("SELECT g.ensg, NULL AS ensp, g.species, g.fullspecies, g.symbol FROM ensembl_gff3_gene g, ensembl_gff3_transcript t WHERE t.enst='$ensp1' AND t.ensg=g.ensg");
			if (Numrows($query) == 0)
			{
				addme("ensp not found in table 'ensembl' nor as an enst in table 'ensembl_gff3_transcript' (skipped)", $ensp1);
				next;
			}
			else
			{
				$enst1 = $ensp1;
				($ensg1, $ensp1, $species1, $fullspecies1, $symbol1) = FetchOne($query);
			}
		}
		elsif (Numrows($query) > 0)
		{
			($ensg1, $enst1, $species1, $fullspecies1, $symbol1) = FetchOne($query);
		}

		# enst2
		# species2
		# fullspecies2
		$query = Query("SELECT ensg, enst, species, fullspecies FROM ensembl WHERE ensp='$ensp2'");
		if (Numrows($query) == 0)
		{
			# # Skipping this for speed (these ensembl_gff3 queries were quite slow, and I don't really need non-coding transcripts)
			# Try ENSP as ENST (these are non-coding transcripts that I'm only keeping in the comparahomology table, not in any of the comparanopara output)
			# via ensembl (won't work since these are non-coding, and ensembl only contains protein_coding)
			# $query = Query("SELECT ensg, ensp, species, fullspecies, symbol FROM ensembl WHERE enst='$ensp2'");
			# via ensembl_gff3 tables
			$query = Query("SELECT g.ensg, NULL AS ensp, g.species, g.fullspecies, g.symbol FROM ensembl_gff3_gene g, ensembl_gff3_transcript t WHERE t.enst='$ensp2' AND t.ensg=g.ensg");
			if (Numrows($query) == 0)
			{
				addme("ensp not found in table 'ensembl' nor as an enst in table 'ensembl_gff3_transcript' (skipped)", $ensp2);
				next;
			}
			else
			{
				$enst2 = $ensp2;
				($ensg2, $ensp2, $species2, $fullspecies2, $symbol2) = FetchOne($query);
			}
		}
		elsif (Numrows($query) > 0)
		{
			($ensg2, $enst2, $species2, $fullspecies2) = FetchOne($query);
		}


		# acc1
		# name1
		$accquery1 = Query("SELECT NULL, NULL");
		if (defined($ensp1))
		{
			$accquery1 = Query("SELECT DISTINCT name, acc FROM uniens WHERE ensp='$ensp1'");
			if (Numrows($accquery1) == 0)
			{
				$accquery1 = Query("SELECT NULL, NULL");
			}
		}
		while (($name1, $acc1) = Fetch($accquery1))
		{
			# acc2
			# name2
			$accquery2 = Query("SELECT NULL, NULL");
			if (defined($ensp2))
			{
				$accquery2 = Query("SELECT DISTINCT name, acc FROM uniens WHERE ensp='$ensp2'");
				if (Numrows($accquery2) == 0)
				{
					$accquery2 = Query("SELECT NULL, NULL");
				}
			}
			while (($name2, $acc2) = Fetch($accquery2))
			{
				# Input fields (from Mart files)
				# if (!defined($ensp1)) { $ensp1 = ''; }
				# if (!defined($ensp2)) { $ensp2 = ''; }
				# if (!defined($ensg2)) { $ensg2 = ''; }
				# if (!defined($symbol2)) { $symbol2 = ''; }
				# if (!defined($homology)) { $homology = ''; }
				# if (!defined($hc)) { $hc = ''; }
				# if (!defined($clade)) { $clade = ''; }
				# if (!defined($target_percent_id)) { $target_percent_id = ''; }
				# if (!defined($query_percent_id)) { $query_percent_id = ''; }
				# if (!defined($gene_order_score)) { $gene_order_score = ''; }
				# if (!defined($whole_genome_align_score)) { $whole_genome_align_score = ''; }
				# if (!defined($target_chr)) { $target_chr = ''; }
				# if (!defined($target_chr_start)) { $target_chr_start = ''; }
				# if (!defined($target_chr_stop)) { $target_chr_stop = ''; }
				# if (!defined($target_chr_strand)) { $target_chr_strand = ''; }

				# Additional fields (from MySQL) (tolerate if these are blank)
				# if (!defined($ensg1)) { $ensg1 = ''; }
				# if (!defined($ensg2)) { $ensg2 = ''; }
				# if (!defined($enst1)) { $enst1 = ''; }
				# if (!defined($enst2)) { $enst2 = ''; }
				if (!defined($ensp1)) { $ensp1 = ''; }
				if (!defined($ensp2)) { $ensp2 = ''; }
				if (!defined($symbol1)) { $symbol1 = ''; }
				if (!defined($symbol2)) { $symbol2 = ''; }
				# if (!defined($fullspecies1)) { $fullspecies1 = ''; }
				# if (!defined($fullspecies2)) { $fullspecies2 = ''; }
				if (!defined($species1)) { $species1 = ''; }
				if (!defined($species2)) { $species2 = ''; }
				if (!defined($acc1)) { $acc1 = ''; }
				if (!defined($acc2)) { $acc2 = ''; }
				if (!defined($name1)) { $name1 = ''; }
				if (!defined($name2)) { $name2 = ''; }

				if ($symbol1 eq '\N') { $symbol1 = ''; }
				if ($symbol2 eq '\N') { $symbol2 = ''; }

				# Convert full species names to lower case (they're otherwise a mix of upper and lower case, e.g. "Homo sapiens" and "homo sapiens", apparently based on the file format above)
				$fullspecies1 = lc($fullspecies1);
				$fullspecies2 = lc($fullspecies2);

				# Insert into comparahomology table (for future use, not used yet)
				# state($_);
				$q = "INSERT INTO $comparahomology SET ensg1='$ensg1', enst1='$enst1', ensp1='$ensp1', acc1='$acc1', name1='$name1', species1='$species1', fullspecies1='$fullspecies1', symbol1='".esc($symbol1)."', ensg2='$ensg2', enst2='$enst2', ensp2='$ensp2', acc2='$acc2', name2='$name2', symbol2='".esc($symbol2)."', species2='$species2', fullspecies2='$fullspecies2', homology='$homology', hc='$hc', clade='$clade', target_percent_id='$target_percent_id', query_percent_id='$query_percent_id', gene_order_score='$gene_order_score', whole_genome_align_score='$whole_genome_align_score', target_chr='$target_chr', target_chr_start='$target_chr_start', target_chr_stop='$target_chr_stop', target_chr_strand='$target_chr_strand'";
				$q =~ s/=''/=NULL/g;
				$q =~ s/='\\N'/=NULL/g;
				Query($q);
			}
		}



		# Skip anything but one2one orthologs
		next if ($homology ne 'ortholog_one2one');

		# Skip anything that isn't high-confidence
		next if ($hc ne '1' and switch('hc'));

		# Skip ENSTs (the ones I looked at were non-coding RNAs)
		next if ($ensp1 =~ /^ENST\d+$/);
		
		# Handling for non-coding transcripts (which only end up in comparahomology - they aren't carried forward below)
		if (($ensp1 ne '') and ($ensp2 ne ''))
		{
			# Add to mart array (carried forward below)
			if (exists($mart{$ensp1}))
			{
				$mart{$ensp1} .= "$ensp2|";
			}
			else
			{
				$mart{$ensp1} = "$ensp2|";
			}
		}

		# Don't count non-coding transcripts (otherwise we'd get a blank ENSP entry)
		addme("total mart ensp1s", $ensp1) if ($ensp1 ne '');
		addme("total mart ensp2s", $ensp2) if ($ensp2 ne '');
	}
	
	stepme(1);
}
stopme();
stoptime();

state("Total ENSPs for Mart: ".scalar(keys(%mart)));

startme("Getting sequences and writing protein and CDS FASTA files");
foreach $ensp (sort(keys(%mart)))
{
	# Get human (reference) protein sequence from comparafasta
	# $query = Query("SELECT REPLACE(seq, '-', '') FROM comparafasta WHERE ensp='$ensp' AND species='homo_sapiens'");
	# Get it from ensembl instead
	$query = Query("SELECT seq, cds FROM ensembl WHERE ensp='$ensp' AND fullspecies='homo_sapiens'");
	if (Numrows($query) == 0)
	{
		# addme("no sequence in comparafasta for ensp", $ensp);
		addme("no sequence in ensembl for ensp", $ensp);
		next;
	}
	# $seq = FetchOne($query);
	($seq, $cds) = FetchOne($query);

	# # Get human (reference) CDS sequence from comparafasta
	# # $query = Query("SELECT REPLACE(seq, '-', '') FROM comparacds WHERE ensp='$ensp' AND species='homo_sapiens'");
	# # Get it from ensembl instead (which now contains CDS, and I've verified that it translates perfectly to the AA sequence)
	# $query = Query("SELECT cds FROM ensembl WHERE ensp='$ensp' AND species='homo_sapiens'");
	# if (Numrows($query) == 0)
	# {
	# 	# addme("no sequence in comparacds (but in comparafasta!) for ensp", $ensp);
	# 	die("Error: Found sequence, but not CDS in table ensembl for ensp '$ensp' (shouldn't happen, same row, so if it's missing from table 'ensembl' it should already have been skipped here)");
	# 	# next;
	# }
	# $cds = FetchOne($query);

	# #Can skip these verification steps for speed (ensembl cds and seq are already verified)
	# # Verify protein sequence (only accept the standard 20 AAs, else skip)
	# if (!aa($seq))
	# {
	# 	# addme("sequence in comparafasta contains non-AA characters for ensp (skipped)", $ensp);
	# 	addme("AA sequence in ensembl contains non-AA characters for ensp (skipped)", $ensp);
	# 	next;
	# }
	# # Verify CDS (only accept the standard 4 bases, else skip)
	# if (!dna($cds))
	# {
	# 	# addme("sequence in comparacds contains non-ACGT characters for ensp (skipped)", $ensp);
	# 	addme("CDS in ensembl contains non-ACGT characters for ensp (skipped)", $ensp);
	# 	next;
	# }
	# # Verify that CDS translates to protein sequence
	# if (translate($cds) ne $seq)
	# {
	# 	# addme("sequence in comparacds does not translate to comparafasta protein sequence for ensp (skipped)", $ensp);
	# 	addme("CDS in ensembl does not translate to comparafasta protein sequence for ensp (skipped)", $ensp);
	# 	next;
	# }

	# Print human (reference) protein sequence
	$outfile = "input/$ensp.txt";
	open(OUT, ">$outfile") or die("Error: Couldn't open '$outfile'");
	fastabreak();
	print OUT ">$ensp\n".split60($seq)."\n";

	# Print human (reference) CDS sequence
	$cdsfile = "input/$ensp.cds.txt";
	open(CDS, ">$cdsfile") or die("Error: Couldn't open '$cdsfile'");
	fastabreak();
	print CDS ">$ensp\n".split60($cds)."\n";
	
	# Write homolog sequences (rest of the homology cluster)
	foreach $homoensp (unique(split(/\|/, $mart{$ensp})))
	{
		# Protein

		# From comparafasta
		# $query = Query("SELECT REPLACE(seq, '-', '') FROM comparafasta WHERE ensp='$homoensp'");
		# From ensembl
		$query = Query("SELECT seq, cds FROM ensembl WHERE ensp='$homoensp'");
		if (Numrows($query) == 0)
		{
			# addme("no sequence in comparafasta for homoensp", $homoensp);
			addme("no sequence in ensembl for homoensp", $homoensp);
			next;
		}
		# $seq = FetchOne($query);
		($seq, $cds) = FetchOne($query);

		print OUT ">$homoensp\n".split60($seq)."\n";

		# # CDS
		# 
		# $query = Query("SELECT REPLACE(seq, '-', '') FROM comparacds WHERE ensp='$homoensp'");
		# if (Numrows($query) == 0)
		# {
		# 	addme("no sequence in comparacds (but in comparafasta!) for homoensp", $homoensp);
		# 	next;
		# }
		# $cds = FetchOne($query);

		print CDS ">$homoensp\n".split60($cds)."\n";
	}

	close(OUT);
	close(CDS);

	stepme(100);
}
stopme();
stoptime();

# showmesome();
showmeall();

done();
