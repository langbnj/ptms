#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');
use List::Util 'shuffle';

# initialize
# $nseqmin = 10;
$nseqmin = 1;
$factordefault = 30000;

our $usage = "$0 [mode: 'fast', 'slow'] [maximum GB of RAM to use (estimated)]";
($mode, $maxgb) = args(2);
die("Error: mode must be 'fast' or 'slow'") if (($mode ne 'fast') and ($mode ne 'slow'));

# start
$factor = $factordefault;
print " >> Attempting to get previous output...\n";
$outfile = "output-memtest-$mode.txt";
@f = ();
if (-e $outfile)
{
	open(IN, "$outfile") or die("\nError: Couldn't open '$outfile'\n\n");
	while (<IN>)
	{
		chomp;
		
		next if /^Alignment\t/;		# skip header
		
		@a = split(/\t/, $_);
		print "   >> $a[8]\n";
		push(@f, $a[8]);
		if ($factor < $a[8])
		{
			# set factor to max
			$factor = $a[8];
		}
	}
	close(IN);
}

if (scalar(@f) > 0)
{
	print " >> Using maximum factor of $factor bytes\n";
}
else
{
	print " >> Not found, using default factor of $factor bytes\n";
}
print "\n";


open(OUT, ">>$outfile") or die("\nError: Couldn't open '$outfile'\n\n");
if (tell(OUT) == 0)
{
	print OUT "Alignment\tNumber of sequences\tLength of alignment\tProduct of length and sequences\tReference ENSP ID\tExpected memory allocation\tExpected factor\tReal allocation\tReal factor\n";
}

@aln = `ls -1 tmp/*.fasta.txt`;
@aln = shuffle(@aln);
foreach $fastafile (@aln)
{
	chomp($fastafile);
	
	@a = `grep -E "^>ENSP[0-9]" $fastafile` or next;
	shuffle(@a);
	
	chomp($a[0]);
	$a[0] =~ s/^>//;
	$ensp = $a[0];
	
	$fastafile =~ /tmp\/(\d+)\.fasta\.txt/;
	$aln = $1;
	
	print " >> Alignment #$aln\n";
	print "   >> Reference sequence $ensp\n";
	
	@a = `grep $aln.fasta.txt output-check.txt` or die;
	die if (scalar(@a) != 1);
	chomp($a[0]);
	@a = split(/\t/, $a[0]);
	$nseq = $a[1];
	$lseq = $a[3];
	$score = $nseq * $lseq;
	$expect = $score * $factor;
	print "     >> $nseq sequences * $lseq alignment length = $score\n";
	
	if ($nseq < $nseqmin)
	{
		print "       >> too few sequences (<$nseqmin), skipping!\n";
		next;
	}
	
	
	print "     >> Expecting memory usage of ".sprintf("%.3f", ($expect) / 1024 / 1024 / 1024)." GB\n";
	
	$treefile = $fastafile;
	$treefile =~ s/\.fasta\./.tree./;
	
	$treeoutfile = "memtest/output-$mode-$aln-treeout.txt";
	$normalizedoutfile = "memtest/output-$mode-$aln-normalizedout.txt";
	$rate4siteoutfile = "memtest/output-$mode-$aln-mainout.txt";
	$logfile = "memtest/log-$mode-$aln.txt";
	
	if ($expect > $maxgb * 1024 * 1024 * 1024)
	{
		print "       >> Expected usage is >".sprintf("%.3f", $maxgb)." GB, skipping!\n";
		next;
	}
	
	print "       >> Running rate4site...\n";

	starttime();
	# print "valgrind rate4site.3.2.mysource_onlycleaned_$mode/sourceMar09/rate4site -bn -a $ensp -t $treefile -s $fastafile -x $treeoutfile -o $normalizedoutfile -y $rate4siteoutfile >& $logfile\n";
	# system("valgrind rate4site.3.2.mysource_onlycleaned_$mode/sourceMar09/rate4site -bn -a $ensp -t $treefile -s $fastafile -x $treeoutfile -o $normalizedoutfile -y $rate4siteoutfile >& $logfile");
	if ($mode eq 'fast')
	{
		system("valgrind rate4site_sources_3.3_fast/programs/rate4site/rate4site -bn -a $ensp -t $treefile -s $fastafile -x $treeoutfile -o $normalizedoutfile -y $rate4siteoutfile >& $logfile");
	}
	elsif ($mode eq 'slow')
	{
		system("valgrind rate4site_sources_3.3_slow/programs/rate4site/rate4site.doubleRep -bn -a $ensp -t $treefile -s $fastafile -x $treeoutfile -o $normalizedoutfile -y $rate4siteoutfile >& $logfile");
	}

	print "       >> ";
	stoptime(1);

	if (`grep -c "likelihoodComputation::getLofPos: likelihood of pos was zero!" $logfile` > 0)
	{
		print "         >> rate4site crashed! ('likelihood of pos was zero')\n";
		next;
	}

	if (`grep -c "Amino acid was not one of the following" $logfile` > 0)
	{
		print "         >> rate4site crashed! ('Amino acid was not one of the following')\n";
		next;
	}



	@a = `cat $logfile | grep allocated` or die("Error: Couldn't find number of bytes allocated in '$logfile'");
	die if (scalar(@a) != 1);
	chomp($a[0]);
	#==24195==   total heap usage: 247,987 allocs, 246,673 frees, 14,630,894 bytes allocated
	#==13906== malloc/free: 406,384 allocs, 405,070 frees, 25,452,887 bytes allocated.
	$a[0] =~ /: [\d,]+ allocs, [\d,]+ frees, ([\d,]+) bytes allocated/ or die("Error: Couldn't find number of bytes allocated in '$logfile', string '$a[0]'");
	$alloc = $1;
	$alloc =~ s/,//g;
	
	print "         >> Real memory usage was ".sprintf("%.3f", ($alloc) / 1024 / 1024 / 1024)." GB\n";

	print "         >> This implies a factor of ~".sprintf("%.0f", $alloc / $score)."\n";
	print "           >> Logging to '$outfile'\n";
	
	print OUT "$aln\t$nseq\t$lseq\t$score\t$ensp\t".sprintf("%.0f", $expect)."\t".sprintf("%.0f", $factor)."\t$alloc\t".sprintf("%.0f", ($alloc/$score))."\n";

	if ($alloc/$score > $factor)
	{
		$factor = $alloc/$score;
		print "           >> Adopting new HIGHER factor of $factor\n";
	}
	else
	{
		print "           >> Keeping current higher factor of $factor\n";
	}
}

done();