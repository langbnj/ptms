#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
$evid = 'ECO:0000312';		# "Another database"
$pmid = '22135298';			# PhosphoSitePlus paper

our $usage = '';
args(0);

# ==> input/unimod_phosphositeplus_K-ac.txt <==
# 121417
# PhosphoSitePlus(R) (PSP) was created by Cell Signaling Technology Inc. It is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License. When using PSP data or analyses in printed publications or in online resources, the following acknowledgements must be included: (a) the words "PhosphoSitePlus(R), www.phosphosite.org" must be included at appropriate places in the text or webpage, and (b) the following citation must be included in the bibliography: "Hornbeck PV, Zhang B, Murray B, Kornhauser JM, Latham V, Skrzypek E PhosphoSitePlus, 2014: mutations, PTMs and recalibrations. Nucleic Acids Res. 2015 43:D512-20. PMID: 25514926."
# 
# GENE	PROTEIN	ACC_ID	HU_CHR_LOC	MOD_RSD	SITE_GRP_ID	ORGANISM	MW_kD	DOMAIN	SITE_+/-7_AA	LT_LIT	MS_LIT	MS_CST	CST_CAT#
# YWHAB	14-3-3 beta	P31946	20q13.12	K5-ac	33347661	human	28.08		___MtMDksELVQkA		2		
# Ywhab	14-3-3 beta	Q9CQV8	2|2 H3	K5-ac	33347661	mouse	28.09		___MtMDksELVQkA		1		
# Ywhab	14-3-3 beta	P35213	3q42	K5-ac	33347661	rat	28.05		___MTMDkSELVQkA		1		
# Ywhab	14-3-3 beta	P35213	3q42	K11-ac	36313495	rat	28.05	14-3-3	DkSELVQkAkLAEQA		1		
# YWHAB	14-3-3 beta	P31946	20q13.12	K13-ac	36297548	human	28.08	14-3-3	sELVQkAkLAEQAER		1		
# Ywhab	14-3-3 beta	P35213	3q42	K13-ac	36297548	rat	28.05	14-3-3	SELVQkAkLAEQAER		1		


nl();
state("Deleting existing source='PhosphoSitePlus' sites from table 'unimod'...", 1);
$query = Query("DELETE FROM unimod WHERE source='PhosphoSitePlus'");
state("Rows affected: ".commify(Numrows($query)), 1);


starttime();
# open(DIR, "ls input/unimod_phosphositeplus_K-ac.txt -1|");
open(DIR, "ls input/unimod_phosphositeplus_*.txt -1|");
while ($infile = <DIR>)
{
	chomp($infile);
	open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
	
	nl(); nl();
	startme("Reading sites from '$infile' (PhosphoSitePlus)", 0, chompme(`cat '$infile' | wc -l`));

	# Process headers
	$_ = <IN>; chomp;
	# $version = $_;
	<IN>; <IN>;
	$_ = <IN>; chomp;
	# die("Error: Header doesn't look as expected:\n\n$_\n\n") if ($_ ne 'NAME	ACC#	Gene Symbol	MOD_TYPE	RSD	SITE_GRP_ID	SPECIES	MW (kD)	IN_DOMAIN	MODSITE_SEQ	PUBMED_LTP	PUBMED_MS2	CST_MS2	CST_CAT#');
	# die("Error: Header doesn't look as expected:\n\n$_\n\n") if ($_ ne 'PROTEIN	ACC_ID	GENE	HU_CHR_LOC	MOD_RSD	SITE_GRP_ID	ORGANISM	MW_kD	DOMAIN	SITE_+/-7_AA	LT_LIT	MS_LIT	MS_CST	CST_CAT#');
	die("Error: Header doesn't look as expected:\n\n$_\n\n") if ($_ ne 'GENE	PROTEIN	ACC_ID	HU_CHR_LOC	MOD_RSD	SITE_GRP_ID	ORGANISM	MW_kD	DOMAIN	SITE_+/-7_AA	LT_LIT	MS_LIT	MS_CST	CST_CAT#');
	
	while (<IN>)
	{
		stepme(1000);

		chomp;
	
		@a = split(/\t/);

		# 2022_04:
		# GENE	PROTEIN	ACC_ID	HU_CHR_LOC	MOD_RSD	SITE_GRP_ID	ORGANISM	MW_kD	DOMAIN	SITE_+/-7_AA	LT_LIT	MS_LIT	MS_CST	CST_CAT#
		# YWHAB	14-3-3 beta	P31946	20q13.12	K5-ac	33347661	human	28.08		___MtMDksELVQkA		2		
		# Ywhab	14-3-3 beta	Q9CQV8	2|2 H3	K5-ac	33347661	mouse	28.09		___MtMDksELVQkA		1		
		# Ywhab	14-3-3 beta	P35213	3q42	K5-ac	33347661	rat	28.05		___MTMDkSELVQkA		1		
		# Ywhab	14-3-3 beta	P35213	3q42	K11-ac	36313495	rat	28.05	14-3-3	DkSELVQkAkLAEQA		1		
		# YWHAB	14-3-3 beta	P31946	20q13.12	K13-ac	36297548	human	28.08	14-3-3	sELVQkAkLAEQAER		1		
		# Ywhab	14-3-3 beta	P35213	3q42	K13-ac	36297548	rat	28.05	14-3-3	SELVQkAkLAEQAER		1		

        $thisacc = $a[2];
        # S149-p
        $ptmsite = $a[4];
        $ptmsite =~ /^([ACDEFGHIKLMNPQRSTVWY])(\d+)\-(\w+)$/ or die("Error: Couldn't parse MOD_RSD field '$ptmsite'");
        $aa = $1;
        $site = $2;
        $ptm = $3;
        $thisspecies = $a[6];
        $seq = $a[9];
        $pubmed_smallscale = $a[10];
        $pubmed_largescale = $a[11];
        $cst_largescale = $a[12];
		
		# Make sequence upper case
		# Lower-case letters are modified residues according to https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3245126/
		# However, on their own websites they sometimes leave them as capitals (e.g. P4 here: https://www.phosphosite.org/proteinAction.action?id=465&showAllSites=true))
		# List of all modified residues: ~/update/uniprot/input >> tail -qn+5 unimod_phosphositeplus_* | cut -f10 | grep -o '[a-z]' | suq
		# Not useful information for me (I could use it to verify that all lower-case letters are present with PTMs, but if not, there's nothing I can do to fix it since I don't have the type of PTM etc.)
		$seq = uc($seq);
		# Remove underscores at beginning and end
		$seq =~ s/^_+//g;
		$seq =~ s/_+$//g;
		
        # $site = substr($ptmsite, 1);      # Residue number
        # $aa = substr($ptmsite, 0, 1);     # AA type
		
		# These are just numbers of studies, blank/undefined if none:
		$pubmed_smallscale = '' if (!defined($pubmed_smallscale));
		$pubmed_largescale = '' if (!defined($pubmed_largescale));
		$cst_largescale = '' if (!defined($cst_largescale));
		
		# print "THISACC  = $thisacc\n";
		# print "PTM     = $ptm\n";
		# print "SITE    = $site\n";
		# print "SPECIES = $ptm\n";
		
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
		
		
		# I can't trust the species column. P03355, for example, is listed as MOUSE, but it's actually viral. I'll simply use the species UniProt lists for a given acc instead.
		# # Handle species
		# #      1 starfish
		# #      1 water buffalo
		# #      2 duck
		# #      2 marmoset
		# #      2 papillomavirus
		# #      2 torpedo
		# #      3 cat
		# #      4 turkey
		# #      5 guinea pig
		# #      6
		# #      6 quail
		# #      7 ORGANISM
		# #      9 goat
		# #     10 fruit fly
		# #     14 green monkey
		# #     14 horse
		# #     14 sheep
		# #     27 SARSCoV1
		# #     36 frog
		# #     41 hamster
		# #     49 dog
		# #     68 SARSCoV2
		# #    162 pig
		# #    201 rabbit
		# #    430 chicken
		# #    940 cow
		# #  48054 rat
		# # 147919 mouse
		# # 387825 human
		# if    ($species eq 'human') 		{ $species = 'HUMAN'; }
		# elsif ($species eq 'mouse') 		{ $species = 'MOUSE'; }
		# elsif ($species eq 'rat') 			{ $species = 'RAT'; }
		# elsif ($species eq 'cow') 			{ $species = 'BOVIN'; }
		# elsif ($species eq 'chicken') 		{ $species = 'CHICK'; }
		# elsif ($species eq 'rabbit') 		{ $species = 'RABIT'; }
		# elsif ($species eq 'pig') 			{ $species = 'PIG'; }
		# elsif ($species eq 'green monkey')	{ $species = 'CHLAE'; }
		# elsif ($species eq 'dog') 			{ $species = 'CANLF'; }
		# elsif ($species eq 'horse') 		{ $species = 'HORSE'; }
		# elsif ($species eq 'sheep') 		{ $species = 'SHEEP'; }
		# else
		# {
		# 	addme("rare, unhandled species skipped for species", $species);
		# 	addme("rare, unhandled species skipped for listed acc", $thisacc);
		# 	next;
		# }
		

		# Get sequence (via canonical acc, since some are isoforms)
		while (($name, $canon, $acc, $species) = Fetch($namequery))
		{
			# addme("pubmed_smallscale", $pubmed_smallscale);
			# addme("pubmed_largescale", $pubmed_largescale);
			# addme("cst_largescale", $cst_largescale);
	

		
			addme("total acc|aa|ptm|sites", "$acc|$aa|$ptm|$site");
			addme("total accs", $acc);
	        addme("total ptmsites", $ptmsite);
			addme("total ptms", $ptm);
			addme("total aa|ptms", "$aa|$ptm");
			addme("total $aa|$ptm", "$acc|$site");
			addme("total sites", $site);
			addme("total aas", $aa);
			addme("total listed species", $thisspecies);
			addme("total uniprot species", $species);
			addme("total listed|uniprot species", "$thisspecies|$species");
			addme("total sequences around PTM sites", $seq);



			# Check if site window (±7) matches
			# e.g. S6-p (and S9-p):
			# __MEEPQsDPsVEPP

			# SET @myseq = 'abcdef';
			# SET @myseq = '1234567890';
			# SET @mysite = 10;
			# SET @mywindow = 3;
			# SELECT @myseq, SUBSTRING(@myseq, GREATEST(@mysite - @mywindow, 1), 2 * @mywindow + LEAST(@mysite - @mywindow, 1));

			# # Verify that N-terminal sites are getting parsed
			# SELECT * FROM unimod WHERE source='PhosphoSitePlus' ORDER BY site;
			# # >> Only a handful have site=1 (cat unimod_phosphositeplus*| cut -f5 | g "[A-Z]1-"), and only for S-p and Y-p.
			# # >> These are problematic:
			# #  cat unimod_phosphositeplus*| g "[STY]1-p"
			# # SELECT * FROM uniseq WHERE acc IN ('P0C0S8', 'C0HKE1', 'Q96QV6', 'P04908', 'Q99878', 'Q9BTM1', 'Q93077', 'Q6FI13', 'P62805', 'NP_001139343', 'P20671', 'Q96KK5', 'Q8IUE6', 'Q16777', 'Q7L7L0', 'O70189', 'Q60908');
			# # >> Their methionine gets cleaved off. PhosphoSitePlus calls the subsequent serine "residue 1", but that's incorrect.
			# # >> Ignore these
			# # Verify that C-terminal sites are getting parsed
			# SELECT m.*, LENGTH(s.seq), LENGTH(s.seq) - m.site AS dist FROM unimod m, uniseq s WHERE source='PhosphoSitePlus' AND m.acc=s.acc AND s.type IN ('UniProt', 'UniIso') ORDER BY dist;
			# # >> OK
			
			$q = "SELECT SUBSTRING(seq, GREATEST($site - 7, 1), 2 * 7 + LEAST($site - 7, 1)) FROM uniseq WHERE acc='$acc' AND species='$species' AND type IN ('UniProt', 'UniIso') GROUP BY seq";

			$query = Query($q);
			if (Numrows($query) == 0)
			{
				die("Error: Couldn't get substring window at site '$site' for species '$species' acc '$acc' from table 'uniseq':\n$q\n");
			}
			if (Numrows($query) > 1)
			{
				die("Error: Multiple hits for substring window at site '$site' for species '$species' acc '$acc' from table 'uniseq':\n$q\n");
			}
			($uniseq) = FetchOne($query);
			if (!defined($uniseq))
			{
				addme("sequence around PTM site outside of string for acc (skipped)", $acc);
				addme("sequence around PTM site outside of string for acc|site (skipped)", "$acc|$site");
				next;
			}
			# if (length($uniseq) != 15)
			# {
			# 	addme("sequence around PTM site is wrong length (probably outside of string) for acc (skipped)", $acc);
			# 	addme("sequence around PTM site is wrong length (probably outside of string) for acc|site (skipped)", "$acc|$site");
			# 	addme("sequence around PTM site is wrong length (probably outside of string) for length (skipped)", length($uniseq));
			# 	next;
			# }
			# if ((length($uniseq) < 8) or (length($uniseq) > 15))
			if (length($uniseq) > 15)
			{
				die("Error: Unexpected length for sequence around PTM site uniseq '$uniseq' seq '$seq' (".length($uniseq).") for acc '$acc' and site '$site'");
			}
			if (length($uniseq) != length($seq))
			{
				addme("sequence around PTM site doesn't match for acc (skipped)", $acc);
				addme("sequence around PTM site doesn't match for acc|site (skipped)", "$acc|$site");
				addme("sequence around PTM site doesn't match & unexpected length for acc (skipped)", $acc);
				addme("sequence around PTM site doesn't match & unexpected length for acc|site (skipped)", "$acc|$site");
				next;
			}
			if ($seq ne $uniseq)
			{
				addme("sequence around PTM site doesn't match for acc (skipped)", $acc);
				addme("sequence around PTM site doesn't match for acc|site (skipped)", "$acc|$site");
				next;
			}
			else
			{
				addme("sequence around PTM site matches for acc", $acc);
				addme("sequence around PTM site matches for acc|site", "$acc|$site");
			}



			# Check if site matches
			$q = "SELECT SUBSTRING(seq, $site, 1) FROM uniseq WHERE acc='$acc' AND species='$species' AND type IN ('UniProt', 'UniIso') GROUP BY seq";
			$query = Query($q);
			if (Numrows($query) == 0)
			{
				die("Error: Couldn't get substring sequence at site '$site' for species '$species' acc '$acc' from table 'uniseq':\n$q\n");
			}
			if (Numrows($query) > 1)
			{
				die("Error: Multiple hits for substring sequence at site '$site' for species '$species' acc '$acc' from table 'uniseq':\n$q\n");
			}
			($uniaa) = FetchOne($query);
			if (!defined($uniaa))
			{
				addme("aa at PTM site wrong len string for acc (skipped)", $acc);
				addme("aa at PTM site outside of string for acc|site (skipped)", "$acc|$site");
				next;
			}
			if ($aa ne $uniaa)
			{
				addme("aa at PTM site doesn't match for acc (skipped)", $acc);
				addme("aa at PTM site doesn't match for acc|site (skipped)", "$acc|$site");
				next;
			}
			else
			{
				addme("aa at PTM site matches for acc", $acc);
				addme("aa at PTM site matches for acc|site", "$acc|$site");
			}
		
		

		
		
	        # cat input/unimod_phosphositeplus_*.txt | cut -f 5 | perl -ne 'if (/^([ACDEFGHIKLMNPQRSTVWY])(\d+)\-(\w+)$/) {print $1."-".$3."\n"}' | suq
		
			# Set type and description fields for unimod table
			$type = '';
			$description = '';
			# K-ac
	        # if (($aa eq 'K') and ($ptm eq 'ACETYLATION'))
			if (($aa eq 'K') and ($ptm eq 'ac'))
			{
				$type = 'modified residue';
				$description = 'N6-acetyllysine';
			}
			# elsif (($aa eq 'S') and ($ptm eq 'ACETYLATION'))
			# {
			# 	$type = 'modified residue';
			# 	$description = 'N-acetylserine';
			# }
			# elsif (($aa eq 'T') and ($ptm eq 'ACETYLATION'))
			# {
			# 	$type = 'modified residue';
			# 	$description = 'N-acetylthreonine';		# Could also be O-acetylthreonine? Not important though
			# }
			# elsif (($aa eq 'Y') and ($ptm eq 'ACETYLATION'))
			# {
			# 	$type = 'modified residue';
			# 	$description = 'N-acetyltyrosine';
			# }
			# K-me
	        # elsif (($aa eq 'K') and ($ptm eq 'METHYLATION'))
			elsif (($aa eq 'K') and ($ptm eq 'me'))
			{
				$type = 'modified residue';
				$description = 'N6-methylated lysine';				# Could also be N6-methyllysine, but I'm pretty sure this is what UniProt uses when it doesn't know whether it's mono, di or tri
			}
	        # elsif (($aa eq 'K') and ($ptm eq 'MONO-METHYLATION'))
			elsif (($aa eq 'K') and ($ptm eq 'm1'))
			{
				$type = 'modified residue';
				$description = 'N6-methyllysine';
			}
	        # elsif (($aa eq 'K') and ($ptm eq 'DI-METHYLATION'))
			elsif (($aa eq 'K') and ($ptm eq 'm2'))
			{
				$type = 'modified residue';
				$description = 'N6,N6-dimethyllysine';
			}
	        # elsif (($aa eq 'K') and ($ptm eq 'TRI-METHYLATION'))
			elsif (($aa eq 'K') and ($ptm eq 'm3'))
			{
				$type = 'modified residue';
				$description = 'N6,N6,N6-trimethyllysine';
			}
			# R-me
	        # elsif (($aa eq 'R') and ($ptm eq 'METHYLATION'))
			elsif (($aa eq 'R') and ($ptm eq 'me'))
			{
				$type = 'modified residue';
				$description = 'Omega-N-methylated arginine';		# Could also be Omega-N-methylarginine, but I'm pretty sure this is what UniProt uses when it doesn't know whether it's mono or di
			}
	        # elsif (($aa eq 'R') and ($ptm eq 'MONO-METHYLATION'))
			elsif (($aa eq 'R') and ($ptm eq 'm1'))
			{
				$type = 'modified residue';
				$description = 'Omega-N-methylarginine';
			}
	        # elsif (($aa eq 'R') and ($ptm eq 'DI-METHYLATION'))
			elsif (($aa eq 'R') and ($ptm eq 'm2'))
			{
				$type = 'modified residue';
				$description = 'Dimethylated arginine';			# I'm pretty sure this is what UniProt uses when it doesn't know whether it's symmetric or asymmetric
			}
			# K-sumo
	        # elsif (($aa eq 'K') and ($ptm eq 'SUMOYLATION'))
			elsif (($aa eq 'K') and ($ptm eq 'sm'))
			{
				$type = 'cross-link';
				$description = 'Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in SUMO)';
			}
			# K-ub
	        # elsif (($aa eq 'K') and ($ptm eq 'UBIQUITINATION'))
			elsif (($aa eq 'K') and ($ptm eq 'ub'))
			{
				$type = 'cross-link';
				$description = 'Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in ubiquitin)';
			}
			# O-linked GlcNAc (there is no N-linked glycosylation in PhosphoSitePlus)
			# SELECT type, description, COUNT(*) AS c FROM unimod WHERE type='glycosylation site' AND description LIKE '%O-linked%' AND (description LIKE '%glcnac%') GROUP BY type, description ORDER BY c DESC;
	        # elsif (($aa eq 'S') and ($ptm eq 'ACETYLATION'))
			elsif (($aa eq 'S') and ($ptm eq 'gl'))
			{
				$type = 'glycosylation site';
				$description = 'O-linked (GlcNAc) serine';
			}
	        # elsif (($aa eq 'T') and ($ptm eq 'ACETYLATION'))
			elsif (($aa eq 'T') and ($ptm eq 'gl'))
			{
				$type = 'glycosylation site';
				$description = 'O-linked (GlcNAc) threonine';
			}
			# O-GalNAc (there is no N-linked glycosylation in PhosphoSitePlus)
			# SELECT type, description, COUNT(*) AS c FROM unimod WHERE type='glycosylation site' AND description LIKE '%O-linked%' AND (description LIKE '%galnac%') GROUP BY type, description ORDER BY c DESC;
	        # elsif (($aa eq 'S') and ($ptm eq 'ACETYLATION'))
			elsif (($aa eq 'S') and ($ptm eq 'ga'))
			{
				$type = 'glycosylation site';
				$description = 'O-linked (GalNAc...) serine';
			}
	        # elsif (($aa eq 'T') and ($ptm eq 'ACETYLATION'))
			elsif (($aa eq 'T') and ($ptm eq 'ga'))
			{
				$type = 'glycosylation site';
				$description = 'O-linked (GalNAc...) threonine';
			}
			# S-p
	        # elsif (($aa eq 'S') and ($ptm eq 'PHOSPHORYLATION'))
			elsif (($aa eq 'S') and ($ptm eq 'p'))
			{
				$type = 'modified residue';
				$description = 'Phosphoserine';
			}
			# T-p
	        # elsif (($aa eq 'T') and ($ptm eq 'PHOSPHORYLATION'))
			elsif (($aa eq 'T') and ($ptm eq 'p'))
			{
				$type = 'modified residue';
				$description = 'Phosphothreonine';
			}
			# Y-p
	        # elsif (($aa eq 'Y') and ($ptm eq 'PHOSPHORYLATION'))
			elsif (($aa eq 'Y') and ($ptm eq 'p'))
			{
				$type = 'modified residue';
				$description = 'Phosphotyrosine';
			}
			else
			{
	            # addme("aa|ptm|acc|site with weird PTM types (presumably rare) (skipped)", "$aa|$ptm|$acc|$site");
	            # next;
				addme("aa|ptm|acc|site with weird PTM types (presumably rare) (kept)", "$aa|$ptm|$acc|$site");
				# $type = 'PhosphoSitePlus';
				$type = 'modified residue';
				$description = $aa.'-'.$ptm;
			}
		
			# Status
		
		# 	# $thisevidence = 'ECO:0000006';		# Experimental
		# 	# This has changed, ECO:0000006 is not in use anymore:
		# 	# SELECT * FROM unimod WHERE evidence='ECO:0000006';
		# 	# >> 0 rows
		# 	# Information on what UniProt now uses: https://www.uniprot.org/help/evidences#ECO:0000269
		# 	# Common evidence codes in current unimod:
		# 	# SELECT scale, evidence, COUNT(*) AS c FROM unimod GROUP BY scale, evidence ORDER BY c DESC;
		# 	# 		ECO:0000250	142670	# >> By similarity
		# 	# 					115877
		# 	#	 	ECO:0007744	68708	# Combinatorial evidence, manual assertion. "information inferred from a combination of experimental and computational evidence. It is currently used in UniProtKB for published large-scale proteomics data"
		# 	# 		ECO:0000255	45726	# Sequence model evidence, manual assertion
		# 	# small	ECO:0000269	26714	# manually curated information for which there is published experimental evidence, small-scale study
		# 	# large	ECO:0000269	7734	# manually curated information for which there is published experimental evidence, large-scale study
		# 	# small	ECO:0000305	1945
		# 	# 		ECO:0000269	1330
		# 	# small	ECO:0000269|ECO:0007744	1138
		# 	# 		ECO:0000250|ECO:0000255	972
		# 
		# 	# >> Should use ECO:0007829 for these mass spec data: Combinatorial evidence, automatic assertion.
		# 	$thisevidence = 'ECO:0007829';		# Experimental
		
		


		
		
			# Source studies
			
			# PhosphoSitePlus groups PTM site sources into three classes:
			# - 794,556   sites from small-scale low-throughput (LTP) studies
			# - 1,874,869 sites from large-scale high-throughput (HTP/MS) studies
			# - 2,656,908 sites from large-scale high-throughput (HTP/MS) studies done by Cell Signaling Technologies (CST) (the majority!)
			# They do give plenty of caveats about not trusting high-throughput mass spec studies: https://www.phosphosite.org/staticAboutPhosphosite.action#eleven

			# Set study scale

			$thisscale = '';
			if (($pubmed_smallscale eq '') and ($pubmed_largescale eq '') and ($cst_largescale eq ''))
			{
				# addme("No source information for uniacc|site", "$acc|$site");
				# next;
				# die("No source information for uniacc|site $acc|$site\n\n$_\n\n");
			
				# Actually, these sites are also okay. All of PhosphoSitePlus is experimental.
				# I looked up two of these sites on the website, and they were actually from high-throughput experiments.
				# Still, I can't be sure about scale, so I'm leaving it blank. This is okay since I've got lots of UniProt entries with source information, but without a scale statement.
				addme("scale unknown for uniacc|site", "$acc|$site");
				$thisscale = '';
			}
			elsif (($pubmed_smallscale ne '') and (($pubmed_largescale eq '') and ($cst_largescale eq '')))
			{
				$thisscale = 'small';
			}
			elsif (($pubmed_smallscale ne '') and (($pubmed_largescale ne '') or ($cst_largescale ne '')))
			{
				$thisscale = 'both';
			}
			elsif (($pubmed_smallscale eq '') and (($pubmed_largescale ne '') or ($cst_largescale ne '')))
			{
				$thisscale = 'large';
			}

		
			# Remove round bracket at end for LIKE query (to also match the "O-linked (GlcNAc...)" entries)
			$desclike = $description;
			$desclike =~ s/\)$//;
		
		# 	# Check if already in unimod
		# 	$query = Query("SELECT scale, sources, evidence FROM unimod WHERE acc='$acc' AND species='$species' AND site='$site' AND type='$type' AND description LIKE '$desclike%'");
		# 	if (Numrows($query) != 0)
		# 	{
		# 		while (($scale, $sources, $evidence) = Fetch($query))
		# 		{
		# 			if (!defined($scale))
		# 			{
		# 				$scale = $thisscale;
		# 			}
		# 			elsif (($scale eq 'small') and (($thisscale eq 'large') or ($thisscale eq 'both')))
		# 			{
		# 				$scale = 'both';
		# 			}
		# 			elsif (($scale eq 'large') and (($thisscale eq 'small') or ($thisscale eq 'both')))
		# 			{
		# 				$scale = 'both';
		# 			}
		# 			# ...else leave $scale as is.
		# 
		# 			# $sources = join("|", unique(split(/\|/, $thissource."|".$sources)));
		# 			if (!defined($sources))
		# 			{
		# 				$sources = $thissource;
		# 			}
		# 			else
		# 			{
		# 				$sources .= "|$thissource";
		# 			}
		# 
		# 			if (!defined($evidence))
		# 			{
		# 				$evidence = $thisevidence;	# experimental
		# 			}
		# 			else
		# 			{
		# 				$evidence .= "|$thisevidence";
		# 			}
		# 
		# 			addme("already known PTM site for acc|site", "$acc|$site");
		# 		
		# 			# Update unimod
		# 			$q = "UPDATE unimod SET scale='$scale', sources='$sources', evidence='$evidence' WHERE acc='$acc' AND species='$species' AND site='$site' AND type='$type' AND description LIKE '$desclike%'";
		# 			$q =~ s/=''/=NULL/g;
		# 			if (switch('debug'))
		# 			{
		# 				print "$q\n" unless switch('quiet');
		# 			}
		# 			else
		# 			{
		# 				Query($q);
		# 			}
		# 		}
		# 	}
		# 	else
		# 	{
				# Insert into unimod
				$q = "INSERT INTO unimod SET name='$name', acc='$acc', canon='$canon', species='$species', site='$site', ptm='', type='$type', description='$description', source='PhosphoSitePlus', subset='', scale='$thisscale', evid='$evid', pmids='$pmid'";
				$q =~ s/=''/=NULL/g;
				if (switch('debug'))
				{
					print "$q\n" unless switch('quiet');
				}
				else
				{
					Query($q);
				}

				addme("successfully added PTM site for acc|aa|ptm|site", "$acc|$aa|$ptm|$site");
			}
		
			# addme("total PTM sites for acc|sites", "$acc|$site");
		# }
	}
	stopme();
	
	# showmesome(20);
	# showmeall(1);
	showmeallsorted(1);
	clearall();							# clear all lists
	
	done();
}
stoptime();

Optimize('unimod');

done();
