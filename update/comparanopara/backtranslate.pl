#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [species] [MAFFT mode] [-tree] [-1para]\n\n Species: Reference species.\n MAFFT mode: linsi, einsi or ginsi (for L-INS-i, E-INS-i or G-INS-i)\n -tree: Use MAFFT results with species tree.\n -1para: Keeps the best-matching outparalog (by various metrics) for one2many and many2many cases.\n\nExample: $0 human linsi -tree";
($unispec, $mafftmode) = args(2);

if (switch('tree'))
{
	$mafftmode .= ".tree";
}

# Handling for 1para
if (switch('1parahc'))
{
	$label = '1parahc';
	$mafftmode .= '.'.$label;
}
elsif (switch('1para'))
{
	$label = '1para';
	$mafftmode .= '.'.$label;
}



# start

$query = Query("SELECT species FROM comparaspecies WHERE unispec='$unispec'");
($species) = FetchOne($query);

if (!switch('1para') and !switch('1parahc'))
{
	# normal
	# $mainquery = Query("SELECT c.ensp, c.aln FROM comparafasta c, uniens e, unimod m WHERE c.species='$species' AND e.ensp=c.ensp AND m.name=e.name AND m.ptm!='' GROUP BY c.ensp, c.aln ORDER BY c.ensp, c.aln");
	# $mainquery = Query("SELECT c.ensp, c.aln FROM comparafasta c, uniens e WHERE c.species='$species' AND e.ensp=c.ensp GROUP BY c.ensp, c.aln ORDER BY c.ensp, c.aln");
	$mainquery = Query("SELECT DISTINCT ensp1 FROM comparahomology WHERE species1='$unispec' AND homology='ortholog_one2one' AND hc=1 AND ensp1 IS NOT NULL AND ensp2 IS NOT NULL ORDER BY ensp1");
}
elsif (switch('1parahc'))
{
	# 1parahc
	# 1para
	if (!switch('debug'))
	{
		$mainquery = Query("SELECT DISTINCT ensp1 FROM comparahomology WHERE species1='$unispec' AND hc=1 AND ensp1 IS NOT NULL AND ensp2 IS NOT NULL ORDER BY ensp1");
	}
	else
	{
		$mainquery = Query("SELECT DISTINCT ensp1 FROM comparahomology WHERE species1='$unispec' AND hc=1 AND ensp1 IS NOT NULL AND ensp2 IS NOT NULL ORDER BY ensp1 LIMIT 1");
	}
}
else
{
	# 1para
	if (!switch('debug'))
	{
		$mainquery = Query("SELECT DISTINCT ensp1 FROM comparahomology WHERE species1='$unispec' AND ensp1 IS NOT NULL AND ensp2 IS NOT NULL ORDER BY ensp1");
	}
	else
	{
		$mainquery = Query("SELECT DISTINCT ensp1 FROM comparahomology WHERE species1='$unispec' AND ensp1 IS NOT NULL AND ensp2 IS NOT NULL ORDER BY ensp1 LIMIT 1");
	}
}

startme("Back-translating MAFFT alignments for '$unispec' ($species), mode $mafftmode to cDNA alignments", 0, Numrows($mainquery));
starttime();

while (($ensp) = Fetch($mainquery))
{
	stepme(10);
	# next if (getme() < 3860);
	
	$outfile = "output/cds.$ensp.$mafftmode.txt";
	# $outfile_pal2naltest = "tmp/cds.pal2naltest.$ensp.$mafftmode.txt";
	print(" >> $outfile\n") if (switch('debug'));
    
    if (-s $outfile)
    {
	    # Skip this ensp if output file already exists and is non-zero
        addme("output file already exists for ensp (skipped)", $ensp);
        next;
		# # Crash instead
		# die("Error: Alignment file '$outfile' already exists for ensp '$ensp' (should run clean.pl first)");
    }

	# Open input MAFFT file
	$mafftfile = "output/mafft.$ensp.$mafftmode.txt";
	if (!-s $mafftfile)
	{
		addme("couldn't find FASTA file for ensp (probably because it doesn't have any one-to-one orthologs in the orthologue mart files) (skipped)", $ensp);
		next;
	}
	

	# Open input file for this ENSP
	if (switch('1para') or switch('1parahc'))
	{
		$infile = "input/$label/$ensp.$label.txt";
	}
	else
	{
		$infile = "input/$ensp.txt";
	}
	print("   >> $infile\n") if (switch('debug'));
	if (!-s $infile)
	{
		addme("couldn't find FASTA file for ensp (skipped)", $ensp);
		next;
	}
	open(IN, $infile) or die("Error: Couldn't open '$infile'");
	fastabreak();
	
	# Read input file
	%seq = ();
	while (<IN>)
	{
		($title, $seq) = getfasta();
	
		$seq{$title} = $seq;
	}
	normalbreak();
	close(IN);

	# Open input CDS file
	if (switch('1para') or switch('1parahc'))
	{
		$cdsfile = "input/$label/$ensp.$label.cds.txt";
	}
	else
	{
		$cdsfile = "input/$ensp.cds.txt";
	}
	print("   >> $cdsfile\n") if (switch('debug'));
	if (!-s $cdsfile)
	{
		addme("couldn't find CDS FASTA file for ensp (skipped)", $ensp);
		next;
	}
	open(CDS, $cdsfile) or die("Error: Couldn't open '$cdsfile'");
	fastabreak();
	
	# Read CDS file
	%cds = ();
	while (<CDS>)
	{
		($title, $seq) = getfasta();
	
		$cds{$title} = $seq;
	}
	normalbreak();
	close(CDS);

	# Convert MAFFT protein alignment to the corresponding CDS alignment
	open(OUT, ">$outfile") or die("Error: Couldn't open '$outfile'");
	open(MAFFT, $mafftfile) or die("Error: Couldn't open '$mafftfile'");
	# Read MAFFT file and write to output file
	fastabreak();
	# $selenocysteine = 0;
	while (<MAFFT>)
	{
		($title, $seq) = getfasta();
		$pos = 0;
		$outcds = '';
		foreach $aa (split(//, $seq))
		{
			if ($aa eq '-')
			{
				# Gap
				$outcds .= '---';
			}
			else
			{
				# Coding
				$tmpcds3 = substr($cds{$title}, $pos * 3, 3);
				# $trancds3 = translate_loosely($tmpcds3);
				$trancds3 = translate($tmpcds3);
				if ($trancds3 ne $aa)
				{
					# # Tolerate selenocysteines
					# if (($tmpcds3 eq 'TGA') and ($aa eq 'U'))
					# {
					# 	addme("translation of TGA to U (selenocysteine) tolerated for ensp (kept)", $ensp);
					# 	$selenocysteine = 1;
					# }
					# else
					# {
						die("Error: Translation of codon '$tmpcds3' is '".translate_loosely($tmpcds3)."', but it should be aa '$aa'");
						# warn("Warning: Translation of codon '$tmpcds3' is '".translate_loosely($tmpcds3)."', but it should be aa '$aa'");
					# }
				}
				$outcds .= $tmpcds3;
				$pos++;
			}
		}

		# Check translation
		# $tranoutcds = translate_loosely($outcds);
		$tranoutcds = translate($outcds);
		$tmpseq = $seq;
		# if ($selenocysteine == 1) { $tmpseq =~ s/U/\*/g;	}
		# die("\n\nError: Translation of CDS\n\nOUTCDS $outcds\n\nis\n\nTRANSEQ $tranoutcds\n\n...but it should be sequence\n\nSEQ     $seq\n\n") if ($tranoutcds ne $tmpseq);
		warn("\n\nWarning:\n\nENSP $ensp\n\nMAFFTFILE $mafftfile\n\nTranslation of $ensp CDS\n\nOUTCDS $outcds\n\nis\n\nTRANSEQ $tranoutcds\n\n...but it should be sequence\n\nSEQ     $seq\n\n") if ($tranoutcds ne $tmpseq);

		print OUT ">$title\n".split60($outcds)."\n";
	}
	normalbreak();
	close(MAFFT);


# No longer necessary - used this to verify that my own CDS generation works perfectly.
# 	# Use PAL2NAL by Peer Bork's group for backtranslation (just for verification, now - building my own gapped CDS above)
# 	addme("ran pal2nal for ensp", $ensp);
# 	# addme("ran pal2nal for aln", $aln);
# 	$errorfile = "tmp/log-errors-pal2nal-cds.$ensp.$mafftmode.errors.txt";
#     # system("bin/pal2nal.pl $mafftfile $cdsfile -output fasta > $outfile 2> $errorfile");
# 	# system("bin/pal2nal.pl $mafftfile $cdsfile -output fasta -nomismatch > $outfile 2> $errorfile");
# 	# -nomismatch "removes mismatched codons from the output", but that would also lead to issues (the CDS would be too short).
# 	# state("bin/pal2nal.v14/pal2nal.pl $mafftfile $cdsfile -output fasta -nomismatch > $outfile_pal2naltest 2> $errorfile");
# 	# system("bin/pal2nal.v14/pal2nal.pl $mafftfile $cdsfile -output fasta -nomismatch > $outfile_pal2naltest 2> $errorfile");
# 	state("bin/pal2nal.v14/pal2nal.pl $mafftfile $cdsfile -output fasta > $outfile_pal2naltest 2> $errorfile");
# 	system("bin/pal2nal.v14/pal2nal.pl $mafftfile $cdsfile -output fasta > $outfile_pal2naltest 2> $errorfile");
# 
# 	# Check if PAL2NAL produced errors
# 	if (-z $errorfile)
# 	{
# 			addme("backtranslation via pal2nal completed successfully for ensp", $ensp);
# 
# 			# # Compare to my own CDS generation
# 			# # Only if there is no selenocysteine in the sequence (UGA -> U instead of *)
# 			# if ($selenocysteine == 0)
# 			# {
# 				# $diff = chompme(`cmp $outfile $outfile_pal2naltest`);
# 				$diff = chompme(`diff $outfile $outfile_pal2naltest`);
# 				if ($diff ne '')
# 				{
# 					die("\n\nError: My own CDS generation doesn't match PAL2NAL's:\n\nMINE $outfile\n\nPAL2NAL $outfile_pal2naltest\n\n");
# 				}
# 			# }
# 	}
# 	else
# 	{
# 		addme("backtranslation via pal2nal produced errors or warnings for ensp (skipped)", $ensp);
# 		addme("backtranslation via pal2nal produced errors or warnings for error file (skipped)", $errorfile);
# 		run("Remove pal2nal output file", "rm -f $outfile_pal2naltest", 1);
# 	}
}
stopme();
stoptime();

# showmeall();
showmesome(50);

print("\nNote: ENSP00000467141 is titin and always crashes MAFFT.\n\n");

done();
