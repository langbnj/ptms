#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
$dbptm_evid = 'ECO:0000312';		# "Another database"
$dbptm_pmid = '34788852';			# dbPTM paper

# our $superloudmysql = 1;

our $usage = "$0";
# ($var) = args(1);
args(0);

$infile = "input/unimod_dbptm.txt";
# $outfile = "output.txt";

open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");


# start

nl();
state("Deleting existing source='dbPTM' sites from table 'unimod'...", 1);
if (!switch('debug'))
{
	$query = Query("DELETE FROM unimod WHERE source='dbPTM'");
}
else
{
	$query = Query("SELECT id FROM unimod WHERE source='dbPTM'");
}
state("Rows affected: ".commify(Numrows($query)), 1);

startme("Reading '$infile' and inserting PTM sites into table 'unimod'", 0, chompme(`cat '$infile' | wc -l`));
starttime();

@numeric_sources = ();
while (<IN>)
{
	chomp;

	print "$_\n" if (switch('debug'));
	
	@a = split(/\t/);

	# ACTO_ACACA	P18281	1	Acetylation	2376577	----------MNPELQSAIGQ
	# RS3A_HUMAN	P61247	155	ADP-ribosylation	29954836	KRNNQIRKTSYAQHQQVRQIR
	# AMP1_PARCM	A0A2I8B346	61	Amidation	29338238	CCASCRRGPIYG---------
	# HEPT_APHF2	A0A0B0QJR1	109	AMPylation	33290744	SGLRNRLVHEYDDIDPNQVFM
	# BLAP_BACSU	C0H419	35	Biotinylation	16699181	GQEVAILESMKMEIPIVADRS
	# GLB2A_ANAIN	P14821	2	Blocked amino end	2599099	---------MVADAVAKVCGS
	# H4_HUMAN	P62805	6	Butyrylation	27105113	-----MSGRGKGGKGLGKGGA
	# ANS1A_HUMAN	Q92625	633	Carbamidation	17322306	LSKSDSDLLTCSPTEDATMGS
	# CBR1_HUMAN	P16152	239	Carboxyethylation	8421682	GWVRTDMAGPKATKSPEEGAE
	# RBL_ARATH	O03042	201	Carboxylation	29372894	ECLRGGLDFTKDDENVNSQPF
	# SHH_MOUSE	Q62226	198	Cholesterol ester	8824192	AENSVAAKSGGCFPGSATVHL
	# THOC4_MOUSE	O08583	140	Citrullination	24463520	TLKKAAVHYDRSGRSLGTADV
	# CO7_HUMAN	P10643	36	C-linked Glycosylation	10551839	VNCQWDFYAPWSECNGCTKTQ
	# H2A1H_HUMAN	Q96KK5	37	Crotonylation	21925322	PVGRVHRLLRKGNYAERVGAG
	# UBACT_YANXG	A0A0G0MF47	56	Deamidation	28087277	EQAKKKQPSRQ----------
	# CO1A1_HUMAN	P02452	170	Deamination	5529814	APQLSYGYDEKSTGGISVPGP
	# GHRL_HUMAN	Q9UBU3	26	Decanoylation	12414809;10604470	LWLDLAMAGSSFLSPEHQRVQ
	# THCL_BACCR	Q812G9	52	Decarboxylation	19196969	TCVCTCSCCTT----------
	# ABL1_HUMAN	P00519	393	Dephosphorylation	11163214	GLSRLMTGDTYTAHAGAKFPI
	# CUTI1_FUSVN	P00590	32	D-glucuronoylation	7398618	PAQELEARQLGRTTRDDLING
	# DEFA4_MOUSE	P28311	64	Disulfide bond	15595831	LHEKSLRGLLCYCRKGHCKRG
	# ASG2_ARATH	Q94BQ3	754	Farnesylation	26147561	FAEGNFHPFECTQS-------
	# GNAO_MOUSE	P18872	205	Formation of an isopeptide bond	23022564	LHFRLFDVGGQRSERKKWIHC
	# BS222_STAPS	A0A0P0C3P7	1	Formylation	26411997	----------MAGLLRFLLSK
	# INS3A_CONGE	A0A0B5A8P4	96	Gamma-carboxyglutamic acid	25605914	FLQRDGRGIVEVCCDNPCTVA
	# YPT1_YEAST	P01123	206	Geranylgeranylation	10071213	QSLTNTGGGCC----------
	# H2B1K_HUMAN	O60814	121	Glutarylation	31542297	AVSEGTKAVTKYTSAK-----
	# CLIC1_HUMAN	O00299	24	Glutathionylation	16286078;14613939;24253303;22833525	AGSDGAKIGNCPFSQRLFMVL
	# XPP2_PIG	Q95333	649	GPI-anchor	7744038	WLQRHTEPLSARAAPTTSLGS
	# INVO_HUMAN	P07476	79	Hydroxyceramide ester	9651377	GLPEQECEQQQKEPQEQELQQ
	# AGP10_ARATH	Q9M0S4	24	Hydroxylation	11006345	LIASSAIAQAPGPAPTRSPLP
	# THYG_HUMAN	P01266	2573	Iodination	2760035	YYSLEHSTDDYASFSRALENA
	# H2B1K_HUMAN	O60814	121	Lactoylation	31645732	AVSEGTKAVTKYTSAK-----
	# A0A384J5P0_BOTFB	A0A384J5P0	61	Lactylation	33193272	
	# GCVHL_STAAM	A0A0H3JT43	56	Lipoylation	26166706	DDEIVSIEASKTVIDVQTPLS
	# H2A2_YEAST	P04912	120	Malonylation	22389435	LPNIHQNLLPKKSAKTAKASQ
	# ACTO_ACACA	P18281	35	Methylation	2376577	APQIENVTVKKVDRSSFLEEV
	# EUDIS_EUDSF	P86455	1	Myristoylation	19053188	----------CPPLC------
	# HIP1_ECTMO	P83341	1	N-carbamoylation	12203680	----------AEKLEESSAEA
	# V9HW98_HUMAN	V9HW98	106	Neddylation	32015554	ICCDILDVLDKHLIPAANTGE
	# VIAAT_MOUSE	O35633	186	Nitration	16800626	DGEVVRVRDSYVAIANACCAP
	# HKT1_ARATH	Q84TI7	429	N-linked Glycosylation	11344270	QRDPINFNVLNITLEVISAYG
	# SPAN2_LAMBD	Q37935	20	N-palmitoylation	8626053	VMMLPLVVVGCTSKQSVSQCV
	# GHRL_RAT	Q9QYH7	26	Octanoylation	10801861;10604470	LWMDMAMAGSSFLSPEHQKAQ
	# RK15_ARATH	P25873	68	O-linked Glycosylation	18431481	FARPLVVVSQTAATSSAVVAP
	# WNT3A_HUMAN	P56704	209	O-palmitoleoylation	24292069;21244856;20826466;25731175	MHLKCKCHGLSGSCEVKTCWW
	# TX21A_PLETR	P0DL38	73	O-palmitoylation	23691198	SNCRCTGTKPSCGKR------
	# DHML_PARDE	P22619	114	Oxidation	1409575	PPGTKLATASWVASCYNPTDG
	# LGG1_CAEEL	Q09490	116	Phosphatidylethanolamine amidation	22767594;20550938	YIAYSDESVYGGEVEKKE---
	# MYSC_ACACA	P10569	311	Phosphorylation	2530230	TTGEQGRGRSSVYSCPQDPLG
	# H4_HUMAN	P62805	9	Propionylation	17267393	--MSGRGKGGKGLGKGGAKRH
	# TX85A_ETHRU	A0A023W0B6	95	Pyrrolidone carboxylic acid	24602922	KRLWRNWERRQVANEDDGEKA
	# HEM3_ARATH	Q43316	316	Pyrrolylation	23519422	RAFLETLDGSCRTPIAGYASK
	# PSD_PLAKH	B3L2V1	308	Pyruvate	25724650	GDEIGEFRMGSSIVVIFENKK
	# HCY_NATPH	P39442	25	S-archaeol	8195126	ATVAAATLAGCNGNGNGNGNG
	# HYPE_ECOLI	P24193	336	S-carbamoylation	15291820;12586941	LPHAEPLPRIC----------
	# HYPE_ECOLI	P24193	336	S-Cyanation	15291820;12586941	LPHAEPLPRIC----------
	# BGH3_HUMAN	Q15582	65	S-cysteinylation	27609313	IGTNRKYFTNCKQWYQRKICG
	# LLP_BPT5	Q38162	16	S-diacylglycerol	8057856	LAMAVVLLSACSTFGPKDIKC
	# RAB4A_MOUSE	P56371	72	Serotonylation	14697203	VKLQIWDTAGQERFRSVTRSY
	# ITIH1_HUMAN	P19827	60	S-linked Glycosylation	9425062;19782370	GVFIRSLKVNCKVTSRFAHYV
	# RIPL1_RAT	D3ZUQ0	47	S-nitrosylation	19607794	EFERVIDQHGCEAIARLMPKV
	# RGS10_HUMAN	O43665	66	S-palmitoylation	10608901	REFLKKEFSEENVLFWLACED
	# POLS_SFV	P03315	1248	Stearoylation	3143715	GAILVLVVVTCIGLRR-----
	# SODC_HUMAN	P00441	123	Succinylation	19608861;24140062	IIGRTLVVHEKADDLGKGGNE
	# PSK1_ORYSI	A2YFB4	80	Sulfation	9371850	SAKWEEFHTDYIYTQDVKNP-
	# CSDE_ECOLI	P0AGF2	61	Sulfhydration	15901727	LKAQAKEIAGCENRVWLGYTV
	# MAST4_HUMAN	O15021	1917	Sulfoxidation	21406390	ARSPGTVMESNPQQREGSSPK
	# ACTG_HUMAN	P63261	328	Sumoylation	19608861	TALAPSTMKIKIIAPPERKYS
	# SAMP1_HALVD	D4GUF6	87	Thiocarboxylation	21368171	ELALFPPVSGG----------
	# H2B1_ARATH	Q9LQQ4	144	Ubiquitination	17554311	AVSEGTKAVTKFTSS------
	# POLG_POL1M	P03300	1546	UMPylation	12502805;16840321	YKLFAGHQGAYTGLPNKKPNV

	# $name = $a[0];
	$thisacc = $a[1];
	$site = $a[2];
	$desc = $a[3];
	$sources = $a[4];
	$seq = $a[5];

	if (scalar(@a) != 6)
	{
		# die("\n\nError: Column count isn't 6, but ".scalar(@a).":\n\nLINE $_\n\n");
		addme("column count is ".scalar(@a)." for description|thisacc|site (skipped)", "$desc|$thisacc|$site");
		next;
	}

	if ($thisacc =~ / $/)
	{
		$thisacc =~ s/ $//;
		addme("cleaned up trailing space for thisacc", $thisacc);
	}

	# Validate listed accession (can be isoform) and update it to primary acc if it is outdated (including outdated isoform base accs)
	# Also get protein name
	if ($thisacc =~ /(-\d+)$/)
	{
		# Isoform
		$iso = $1;
		
		# First, check if this is a current isoform accession:
		$namequery = Query("SELECT name, canon, acc, species FROM uniiso WHERE acc='$thisacc'");
		# If not: try updating the base acc
		if (Numrows($namequery) == 0)
		{
			# Get base acc
			$thisbase = $thisacc;
			$thisbase =~ s/-\d+$//;
			
			# Get updated base acc
			$query = Query("SELECT canon FROM uniacc WHERE acc='$thisbase'");
			if (Numrows($query) == 0)
			{
				addme("listed isoform acc skipped (couldn't update base) for listed acc (skipped)", $thisacc);
				next;
			}
			else
			{
				if (Numrows($query) > 1)
				{
					# If there are multiple current accs for this canon: it's been demerged: skip
					addme("listed isoform acc skipped (tried updating base, but got multiple bases (must have been demerged)) for listed acc (skipped)", $thisacc);
					next;
				}
				else
				{
					($base) = FetchOne($query);

					if ($base eq $thisbase)
					{
						addme("listed isoform acc skipped (tried updating base, but base was up to date) for listed acc (skipped)", $thisacc);
						next;
					}
					$acc = $base.$iso;

					$namequery = Query("SELECT name, canon, acc, species FROM uniiso WHERE acc='$acc'");
					if (Numrows($namequery) == 0)
					{
						addme("listed isoform acc skipped (tried updating base, but still no match in uniiso) for listed acc (skipped)", $thisacc);
						next;
					}
				}
			}
		}
	}
	else
	{
		# Non-isoform
		$namequery = Query("SELECT DISTINCT name, canon, canon AS acc, species FROM uniacc WHERE canon='$thisacc'");
		if (Numrows($namequery) == 0)
		{
			# If this acc isn't a canonical, primary acc:
			# Get updated acc
			$namequery = Query("SELECT DISTINCT name, canon, canon AS acc, species FROM uniacc WHERE acc='$thisacc'");
			if (Numrows($namequery) == 0)
			{
				# If this still didn't work:
				addme("listed non-isoform acc skipped (no match in uniacc) for listed acc (skipped)", $thisacc);
				next;
			}
		}
	}

	# Get sequence (via canonical acc, since some are isoforms)
	while (($name, $canon, $acc, $species) = Fetch($namequery))
	{
		# Verify site in uniseq (is it within the protein?)
		die("Error: Couldn't match site '$site'") if ($site !~ /^\d+$/);
		$query = Query("SELECT DISTINCT seq FROM uniseq WHERE acc='$acc' AND type IN ('UniProt', 'UniIso')");
		($uniseq) = FetchOne($query);
		if ($site >= length($uniseq))
		{
			addme("site outside of uniseq for name (skipped)", $name);
			addme("site outside of uniseq for description|name|site (skipped)", "$desc|$name|$site");
			next;
		}

		# Check UniProt name
		if ($name =~ /^(\w{1,10})_(\w{3,5})$/)
		{
			$species = $2;
			addme("total species", $species);
		}
		else
		{
			die("Error: Couldn't match name '$name'\n\nLINE $_\n\n") ;
		}

		die("Error: Couldn't match site '$site'\n\nLINE $_\n\n") if ($site !~ /^\d+$/);

		# Verify sequence
		if (!defined($seq) or ($seq eq ''))
		{
			addme("sequence is blank for description|name|site (skipped)", "$desc|$name|$site");
			next;
		}
		if ($seq !~ /^[ACDEFGHIKLMNPQRSTVWY\-]+$/)
		{
			if ($seq =~ /^[ACDEFGHIKLMNPQRSTVWYBJZUX\-]+$/)
			{
				addme("sequence contains B/J/Z/U/X for description|name|site (kept)", "$desc|$name|$site");
				if ($seq =~ /U/)
				{
					addme("sequence contains selenocysteine (U) for description|name|site (kept)", "$desc|$name|$site");
				}
				if ($seq =~ /B/)
				{
					addme("sequence contains B for description|name|site (kept)", "$desc|$name|$site");
				}
				if ($seq =~ /J/)
				{
					addme("sequence contains J for description|name|site (kept)", "$desc|$name|$site");
				}
				if ($seq =~ /Z/)
				{
					addme("sequence contains Z for description|name|site (kept)", "$desc|$name|$site");
				}
				if ($seq =~ /X/)
				{
					addme("sequence contains X for description|name|site (kept)", "$desc|$name|$site");
				}
			}
			else
			{
				die("\n\nError: Sequence contains non-AA characters for sequence '$seq'\n\nLINE $_\n\n") ;
			}
		}
		$aa = substr($seq, 10, 1);
		addme("aa expected for $desc", $aa);
		# if ($seq =~ /^-+/)
		# {
		# }
		$expected = $seq;
		$expected =~ s/-//g;
		$start = max(0, $site-11);
		$stop = min($site+10, 21);
		$unisub = substr($uniseq, $start, $stop);
		if ($unisub ne $expected)
		{
			if (switch('debug'))
			{
				warn("Warning: Sequence mismatch: Was supposed to find $aa $desc '$expected', but got '".substr($unisub, 10, 1)."' '$unisub' ($name $aa$site >> ".($start + 1)."-".($start + $stop).") (skipped)\n\n$name:\n$uniseq\n\n");
			}
			addme("sequence mismatch for name (skipped)", $name);
			addme("sequence mismatch for description|name|site (skipped)", "$desc|$name|$site");
			next;
		}

		# Verify sources

		# To see problematic sources:
		# cat ~/update/uniprot/input/unimod_dbptm.txt | cut -f5 | g -v "^[\d;]+$" | g -v "UniProtKB CARBOHYD" | suq
		# g -l "UniProtKB CARBOHYD" ~/update/uniprot/input/dbptm/*
		# >> UniProtKB CARBOHYD only occurs in N-linked Glycosylation and O-linked Glycosylation.
		# >> Lots of other little typos, though ("N.N.", "-", DOIs, "PubMed", etc. for the references)
		# die("Error: Couldn't match source list '$sources'") if ($sources !~ /^[\d;]+$/);
		@sources = ();
		# Default evidence code (unless changed below): "other database"
		$evid = $dbptm_evid;
		# Filter sources
		foreach $source (split(/;/, $sources))
		{
			if ($source eq '0000269')
			{
				# ECO:0000269 (experimental evidence used in manual assertion). That's not a source.
				$source = '';
			}
			elsif ($source eq ' ECO')
			{
				$source = '';
			}
			elsif ($source eq ' E')
			{
				$source = '';
			}
			elsif ($source =~ /^\d+$/)
			{
				# Source is numeric
				if ($source < 10000000)
				{
					addme("DEBUG: total numeric sources below 10,000,000 (unlikely to be PubMed) (source removed)", $source);
					$source = '';
				}
				else
				{
					# Source is numeric (likely PMID)
					addme("total numeric (likely PMID) sources", $source);
					push(@numeric_sources, $source);
				}
			}
			elsif ($source =~ /^\s*(\d+)\s*$/)
			{
				$source = $1;
				# Source is numeric
				if ($source < 10000000)
				{
					# This source is very close to 10,000,000 and it's definitely not a PMID:
					# 10006251
					# https://pubmed.ncbi.nlm.nih.gov/10006251/
					# Phys Rev B Condens Matter
					# . 1993 Jan 15;47(4):2122-2129.  doi: 10.1103/physrevb.47.2122.
					# Self-consistent electronic structure of parabolic semiconductor quantum wells: Inhomogeneous-effective-mass and magnetic-field effects

					# 19,811 numeric sources get filtered out here. Some trash probably remains, but there's nothing I can do to fix this.

					addme("DEBUG: total numeric sources below 10,000,000 (unlikely to be PubMed) (source removed)", $source);
					$source = '';
				}
				else
				{
					# Source is numeric (likely PMID)
					addme("total numeric (likely PMID) sources", $source);
					push(@numeric_sources, $source);
				}
			}
			elsif ($source =~ /^\d+\?$/)
			{
				# Retracted article:
				# https://pubmed.ncbi.nlm.nih.gov/18321849/
				# Lots of fake gels, completely dubious
				# 15843167 has not been retracted (https://pubmed.ncbi.nlm.nih.gov/15843167/), but there's probably a reason for the question mark (skip this source)
				addme("probable retracted article (PMID with question mark) for description|name|site (source removed)", "$desc|$name|$site");
				addme("probable retracted article (PMID with question mark) for source (source removed)", $source);
				$source = '';
			}
			elsif ($source eq 'UniProtKB CARBOHYD')
			{
				# Original UniProt paper describing carbohydrate annotation
				# https://pubmed.ncbi.nlm.nih.gov/11680872/
				$source = '11680872';
				$evid = 'ECO:0000255';
				addme("UniProt-based glycosylation for description|name|site (kept)", "$desc|$name|$site");
			}
			elsif ($source eq '')
			{
				$source = '';
			}
			elsif ($source eq 'doi')
			{
				$source = '';
			}
			elsif ($source eq 'N.N.')
			{
				$source = '';
			}
			elsif ($source eq '-')
			{
				$source = '';
			}
			elsif ($source eq 'PubMed')
			{
				$source = '';
			}
			elsif ($source =~ /^(\d+)PubMed$/)
			{
				$source = $1;
			}
			elsif ($source =~ /^PubMed(\d+)$/)
			{
				# Doesn't actually occur
				$source = $1;
			}
			elsif (($source eq 'doi:10.1007/s13562-013-0225-7') or ($source eq '10.1007/s13562-013-0225-7'))
			{
				# There is no PubMed ID for this article, retaining DOI:
				# https://doi.org/10.1007/s13562-013-0225-7
				$source = "https://doi.org/10.1007/s13562-013-0225-7";
			}
			else
			{
				die("Error: Couldn't match source '$source' in sources '$sources'");
			}

			if ($source ne '')
			{
				push(@sources, $source);
			}
		}

		if (scalar(@sources) == 0)
		{
			addme("no non-dbPTM source given for description|name|site (skipped)", "$desc|$name|$site");
			next;
		}

		# Add dbPTM paper (2022) PMID to source list
		push(@sources, $dbptm_pmid);

		# $sources = join('|', unique(@sources));
		$tmpsources = join('|', @sources);
		# $sources .= "|34788852"


		# die("Error: Couldn't match source list '$sources'") if ($sources !~ /^[\d;]+$/);
		# die("Error: Couldn't match source list '$sources'") if ($sources !~ /^[httpsdoiorg\d;:\/\.\-\|]+$/);	# DOI-tolerant
		die("Error: Couldn't match source list '$tmpsources'") if ($tmpsources !~ /^([\d\|]|https:\/\/doi\.org\/10\.1007\/s13562-013-0225-7|\|)*$/);	# DOI-tolerant

		addme("total names inserted", $name);
		addme("total name|sites inserted", "$name|$site");
		addme("total accs inserted", $acc);
		addme("total acc|sites inserted", "$acc|$site");
		addme("total descriptions inserted", $desc);



		# PTM descriptions
		$type = 'modified residue';
		if ($desc eq "Phosphorylation") 	 { $desc = "$aa-p"; } 		# 1615054 sites
		if ($desc eq "Ubiquitination") 	 { $desc = "$aa-ub"; } 		# 348307 sites
		if ($desc eq "Acetylation") 	 { $desc = "$aa-ac"; } 		# 138169 sites
		if ($desc eq "N-linked Glycosylation") 	 { $desc = "$aa-gly"; } 		# 27361 sites
		if ($desc eq "Succinylation") 	 { $desc = "$aa-suc"; } 		# 17973 sites
		if ($desc eq "O-linked Glycosylation") 	 { $desc = "$aa-gly"; } 		# 16692 sites
		if ($desc eq "Methylation") 	 { $desc = "$aa-me"; } 		# 16114 sites
		if ($desc eq "Malonylation") 	 { $desc = "$aa-mal"; } 		# 12847 sites
		if ($desc eq "Sulfoxidation") 	 { $desc = "$aa-ox"; } 		# 7581 sites
		if ($desc eq "S-palmitoylation") 	 { $desc = "$aa-pal"; $type = 'lipid moiety-binding region'; } 		# 6501 sites
		if ($desc eq "Sumoylation") 	 { $desc = "$aa-sum"; } 		# 5889 sites
		if ($desc eq "S-nitrosylation") 	 { $desc = "$aa-nit"; } 		# 4172 sites
		if ($desc eq "Glutathionylation") 	 { $desc = "$aa-glt"; } 		# 4129 sites
		if ($desc eq "Amidation") 	 { $desc = "$aa-ami"; } 		# 3316 sites
		if ($desc eq "Hydroxylation") 	 { $desc = "$aa-hyd"; } 		# 2404 sites
		if ($desc eq "Neddylation") 	 { $desc = "$aa-ned"; } 		# 1676 sites
		if ($desc eq "Pyrrolidone carboxylic acid") 	 { $desc = "$aa-pyr"; } 		# 965 sites
		if ($desc eq "Glutarylation") 	 { $desc = "$aa-glu"; } 		# 952 sites
		if ($desc eq "Gamma-carboxyglutamic acid") 	 { $desc = "$aa-car"; } 		# 508 sites
		if ($desc eq "Oxidation") 	 { $desc = "$aa-ox"; } 		# 427 sites
		if ($desc eq "Crotonylation") 	 { $desc = "$aa-cro"; } 		# 373 sites
		if ($desc eq "Lactylation") 	 { $desc = "$aa-lac"; } 		# 336 sites
		# if ($desc eq "Myristoylation") 	 { $desc = "Myristoylation"; } 		# 312 sites
		if ($desc eq "Lactoylation") 	 { $desc = "$aa-lac"; } 		# 306 sites
		if ($desc eq "Formylation") 	 { $desc = "$aa-for"; } 		# 256 sites
		# if ($desc eq "Sulfation") 	 { $desc = "Sulfation"; } 		# 252 sites
		# if ($desc eq "C-linked Glycosylation") 	 { $desc = "C-linked Glycosylation"; } 		# 211 sites
		# if ($desc eq "Dephosphorylation") 	 { $desc = "Dephosphorylation"; } 		# 133 sites
		if ($desc eq "Citrullination") 	 { $desc = "$aa-cit"; } 		# 122 sites
		if ($desc eq "ADP-ribosylation") 	 { $desc = "$aa-adp"; } 		# 107 sites
		# if ($desc eq "GPI-anchor") 	 { $desc = "GPI-anchor"; } 		# 103 sites
		# if ($desc eq "Deamidation") 	 { $desc = "Deamidation"; } 		# 92 sites
		# if ($desc eq "Farnesylation") 	 { $desc = "Farnesylation"; } 		# 87 sites
		if ($desc eq "Butyrylation") 	 { $desc = "$aa-but"; } 		# 82 sites
		# if ($desc eq "Nitration") 	 { $desc = "Nitration"; } 		# 81 sites
		# if ($desc eq "Geranylgeranylation") 	 { $desc = "Geranylgeranylation"; } 		# 81 sites
		# if ($desc eq "N-palmitoylation") 	 { $desc = "N-palmitoylation"; } 		# 74 sites
		# if ($desc eq "S-diacylglycerol") 	 { $desc = "S-diacylglycerol"; } 		# 58 sites
		# if ($desc eq "Carboxylation") 	 { $desc = "Carboxylation"; } 		# 44 sites
		# if ($desc eq "Deamination") 	 { $desc = "Deamination"; } 		# 42 sites
		# if ($desc eq "Lipoylation") 	 { $desc = "Lipoylation"; } 		# 35 sites
		# if ($desc eq "Formation of an isopeptide bond") 	 { $desc = "Formation of an isopeptide bond"; } 		# 29 sites
		# if ($desc eq "AMPylation") 	 { $desc = "AMPylation"; } 		# 27 sites
		# if ($desc eq "Blocked amino end") 	 { $desc = "Blocked amino end"; } 		# 26 sites
		# if ($desc eq "Pyruvate") 	 { $desc = "Pyruvate"; } 		# 25 sites
		# if ($desc eq "Carbamidation") 	 { $desc = "Carbamidation"; } 		# 22 sites
		# if ($desc eq "Iodination") 	 { $desc = "Iodination"; } 		# 19 sites
		# if ($desc eq "Thiocarboxylation") 	 { $desc = "Thiocarboxylation"; } 		# 13 sites
		# if ($desc eq "Propionylation") 	 { $desc = "Propionylation"; } 		# 13 sites
		# if ($desc eq "Biotinylation") 	 { $desc = "Biotinylation"; } 		# 12 sites
		# if ($desc eq "UMPylation") 	 { $desc = "UMPylation"; } 		# 10 sites
		# if ($desc eq "Serotonylation") 	 { $desc = "Serotonylation"; } 		# 9 sites
		# if ($desc eq "Phosphatidylethanolamine amidation") 	 { $desc = "Phosphatidylethanolamine amidation"; } 		# 9 sites
		# if ($desc eq "Sulfhydration") 	 { $desc = "Sulfhydration"; } 		# 8 sites
		# if ($desc eq "Octanoylation") 	 { $desc = "Octanoylation"; } 		# 7 sites
		# if ($desc eq "O-palmitoleoylation") 	 { $desc = "O-palmitoleoylation"; } 		# 6 sites
		# if ($desc eq "Disulfide bond") 	 { $desc = "Disulfide bond"; } 		# 6 sites
		# if ($desc eq "Decanoylation") 	 { $desc = "Decanoylation"; } 		# 6 sites
		# if ($desc eq "S-linked Glycosylation") 	 { $desc = "S-linked Glycosylation"; } 		# 5 sites
		# if ($desc eq "Stearoylation") 	 { $desc = "Stearoylation"; } 		# 3 sites
		# if ($desc eq "O-palmitoylation") 	 { $desc = "O-palmitoylation"; } 		# 3 sites
		# if ($desc eq "Hydroxyceramide ester") 	 { $desc = "Hydroxyceramide ester"; } 		# 3 sites
		# if ($desc eq "S-cysteinylation") 	 { $desc = "S-cysteinylation"; } 		# 2 sites
		# if ($desc eq "Decarboxylation") 	 { $desc = "Decarboxylation"; } 		# 2 sites
		# if ($desc eq "Cholesterol ester") 	 { $desc = "Cholesterol ester"; } 		# 2 sites
		# if ($desc eq "S-Cyanation") 	 { $desc = "S-Cyanation"; } 		# 1 sites
		# if ($desc eq "S-carbamoylation") 	 { $desc = "S-carbamoylation"; } 		# 1 sites
		# if ($desc eq "S-archaeol") 	 { $desc = "S-archaeol"; } 		# 1 sites
		# if ($desc eq "Pyrrolylation") 	 { $desc = "Pyrrolylation"; } 		# 1 sites
		# if ($desc eq "N-carbamoylation") 	 { $desc = "N-carbamoylation"; } 		# 1 sites
		# if ($desc eq "D-glucuronoylation") 	 { $desc = "D-glucuronoylation"; } 		# 1 sites
		# if ($desc eq "Carboxyethylation") 	 { $desc = "Carboxyethylation"; } 		# 1 sites


		# Insert into unimod
		$q = "INSERT INTO unimod SET name='$name', acc='$acc', canon='$canon', species='$species', site='$site', ptm='', type='$type', description='$desc', source='dbPTM', subset='', scale='', evid='$evid', pmids='$tmpsources'";
		$q =~ s/=''/=NULL/g;
		if (switch('debug'))
		{
			print "$q\n" unless switch('quiet');
		}
		else
		{
			Query($q);
		}

		addme("added PTM site for acc|site", "$acc|$site");

	}

	stepme(1000);
}
stopme();
stoptime();

state("Numeric sources (likely PMIDs):");
characterise(@numeric_sources);

showme("total numeric sources below 10,000,000 (unlikely to be PubMed) (should probably be skipped)");
showme("column count is 5 for description|name|site (skipped)");
showme("column count is 4 for description|name|site (skipped)");
showme("column count is 3 for description|name|site (skipped)");
showme("column count is 2 for description|name|site (skipped)");
showme("column count is 1 for description|name|site (skipped)");
showme("column count is 0 for description|name|site (skipped)");

# showmeall(1);
showmesome(10);

Optimize('unimod');

done();
