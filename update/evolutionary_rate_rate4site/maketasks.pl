#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');
use POSIX qw(ceil floor);

# initialize

# $nseqmin = 10;                # run only alignments with >= 10 sequences
$nseqmin = 1;					# run all alignments
$gbper = 16;					# RAM per node CPU: 16 GB
$maxgb = 1000;					# max RAM on nodes: 1024 GB (above that, I'll either skip them or try and run them locally)
$factordefault = 50000;			# default memory use factor

# our $usage = "$0 [species] [mode: 'fast', 'slow'] [-noread] [-ptmonly]\n\n -noread: Don't get memory usage factor from output-memtest*.txt, but default to $factordefault\n -ptmonly: Only run proteins that have PTM sites in the species\n\nExample: $0 human fast -noread";
our $usage = "$0 [species] [type] [-noread] [-ptmonly]\n\n -noread: Don't get memory usage factor from output-memtest*.txt, but default to $factordefault\n -ptmonly: Only run proteins that have PTM sites in the species\n\nExample: $0 human linsi_tree -noread";
($species, $mode) = args(2);
# die("Error: mode must be 'fast' or 'slow'") if (($mode ne 'fast') and ($mode ne 'slow'));

# start
run("Delete old task files", "rm -f tmp-task-*") if (!switch('debug'));

# Get Ensembl species 
$query = Query("SELECT species FROM comparaspecies WHERE unispec='$species'");
($comparaspec) = FetchOne($query);
state("UniProt species '$species' >> Ensembl species '$comparaspec'");

$factor = $factordefault;
print " >> Attempting to get previous output...\n";
$outfile = "output-memtest-$mode.txt";
@f = ();
if (-e $outfile and !switch('noread'))
{
	open(IN, "$outfile") or die("\nError: Couldn't open '$outfile'\n\n");
	while (<IN>)
	{
		chomp;
		
		next if /^Alignment\t/;		# skip header
		
		@a = split(/\t/, $_);
		# print "   >> $a[8]\n";
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

startme("Reading $species Compara sequences from 'tmp/*.fasta.txt'");
starttime();
%compara = ();
while (-e "tmp/".sprintf("%06d", getme() + 1).".fasta.txt")
{
	$infile = "tmp/".sprintf("%06d", getme() + 1).".fasta.txt";
	open(IN, "$infile") or die("\nError: Couldn't open '$infile'\n\n");
	fastabreak();
	while (<IN>)
	{
		($title, $seq) = getfasta();
		
		# Remove dashes from sequence (it's an alignment right now)
		$seq =~ s/-//g;
		
		# if $species sequence: store in hash

		# if ($title =~ /^ENSP\d+$/)
		# {
		# 	die("Error: '$title' sequence already exists") if (exists($compara{$title}));
		# 	
		# 	$compara{$title} = $seq;
		# }

		$query = Query("SELECT species FROM comparaenspspecies WHERE ensp='$title'");
		($tmpspec) = FetchOne($query);
		# If this is the desired species..
		if ($tmpspec eq $comparaspec)
		{
			die("Error: '$title' sequence already exists") if (exists($compara{$title}));
			
			$compara{$title} = $seq;
		}
	}
	normalbreak();
	close(IN);
	
	stepme(100);
}
stopme();
stoptime();

state("Read ".scalar(keys(%compara))." $species sequences");

if (switch('debug'))
{
	state("Reading alignments and making task files based on memory requirements");
}
else
{
	startme("Reading alignments and making task files based on memory requirements");
}
startme2();
startme3();

@aln = `ls -1 tmp/*.fasta.txt`;
@aln = unique(@aln);
foreach $fastafile (@aln)
{
	chomp($fastafile);
	
	$fastafile =~ /tmp\/(\d+)\.fasta\.txt/;
	$aln = $1;

	stepme2();

	# skip alignment if it doesn't contain any $species ENSPs
	# @ensps = `grep -E "^>ENSP[0-9]" $fastafile` or next;
	# foreach (@ensps)
	# {
	# 	chomp;
	# 	s/^>//;
	# }
	# @ensps = unique(@ensps);
	@tmp = `grep -E "^>\\S+" $fastafile` or next;
	@ensps = ();
	foreach $tmptitle (@tmp)
	{
		chomp($tmptitle);
		$tmptitle =~ s/^>//;
		
		$query = Query("SELECT species FROM comparaenspspecies WHERE ensp='$tmptitle'");
		($tmpspec) = FetchOne($query);
		# If this is the desired species..
		if ($tmpspec eq $comparaspec)
		{
			push(@ensps, $tmptitle);
		}
	}
	@ensps = unique(@ensps);
	
	print " >> Alignment #$aln\n" if (switch('debug'));

	@a = `grep $aln.fasta.txt output-check.txt` or next;
	die if (scalar(@a) != 1);
	chomp($a[0]);
	@a = split(/\t/, $a[0]);
	$nseq = $a[1];
	$lseq = $a[3];
	$score = $nseq * $lseq;
	$expect = $score * $factor;
	print "   >> $nseq sequences * $lseq alignment length = $score\n" if (switch('debug'));

	if ($nseq < $nseqmin)
	{
		print "     >> too few sequences (<$nseqmin), skipping!\n" if (switch('debug'));
		next;
	}

	$nodes = ceil($expect / ($gbper * 1024 * 1024 * 1024));
	
	print "     >> Expecting memory usage of ".sprintf("%.3f", ($expect) / 1024 / 1024 / 1024)." GB ($nodes nodes)\n" if (switch('debug'));
	
	if ($expect > $maxgb * 1024 * 1024 * 1024)
	{
		print "       >> Expected usage is >".sprintf("%.3f", $maxgb)." GB, assigning to 'extra' task!\n" if (switch('debug'));
		$nodes = 'extra';
	}

	foreach $ensp (@ensps)
	{
		stepme(100) if (!switch('debug'));

		print "   >> Reference sequence $ensp\n" if (switch('debug'));
		
		$mod = 0;
		$query = Query("SELECT name FROM uniens WHERE ensp='$ensp' AND species='$species'");
		if (Numrows($query) == 0)
		{
			addme("no names in uniens for ensp (kept)", $ensp);
            # next;
		}
		else
		{
		    $query = Query("SELECT m.id FROM unimod m, uniens e WHERE e.ensp='$ensp' AND e.name=m.name AND e.species='$species' AND m.ptm!=''");
    		if (Numrows($query) > 0)
    		{
    			$mod = 1;
    		}
		}

		if ($mod == 0)
		{
			print "     >> Doesn't have any PTMs according to unimod 'ptm' field (use switch -ptmonly to skip these)\n" if (switch('debug'));
			next if (switch('ptmonly'));
		}
		else
		{
			print "     >> Has PTMs\n" if (switch('debug'));
		}

		# Get Compara sequence
		if (!exists($compara{$ensp}))
		{
			die("Error: No sequence in Compara for ENSP '$ensp'");
			# addme("no sequence in compara for ensp", $ensp);
			# next;
		}
		
		# Try to get matching UniProt sequence
		$query = Query("SELECT s.seq FROM uniseq s, uniens e WHERE e.ensp='$ensp' AND e.name=s.name AND e.species='$species' AND s.type='UniProt' AND s.seq='$compara{$ensp}'");
		if (Numrows($query) == 0)
		{
			addme("couldn't find a uniprot name with matching sequence for ensp (kept)", $ensp);
            # next;
		}

		$treefile = $fastafile;
		$treefile =~ s/\.fasta\./.tree./;
		
		$outfile = "tmp-task-$nodes.txt";
		print "         >> Adding $aln|$ensp to task file '$outfile'\n" if (switch('debug'));
		addme("task $nodes", "$aln|$ensp");
		if (!switch('debug'))
		{
			open(OUT, ">>$outfile") or die("Error: Couldn't append to '$outfile'");
			print OUT "$aln|$ensp\n";
			close(OUT);
		}

		stepme3();

		addme("added ensps to tasks", $ensp);
		

		# $treeoutfile = "memtest/output-$mode-$aln-treeout.txt";
		# $normalizedoutfile = "memtest/output-$mode-$aln-normalizedout.txt";
		# $rate4siteoutfile = "memtest/output-$mode-$aln-mainout.txt";
		# $logfile = "memtest/log-$mode-$aln.txt";
		# 
		# print "       >> Running rate4site...\n";
		# 
		# starttime();
		# # print "valgrind rate4site.3.2.mysource_onlycleaned_$mode/sourceMar09/rate4site -bn -a $ensp -t $treefile -s $fastafile -x $treeoutfile -o $normalizedoutfile -y $rate4siteoutfile >& $logfile\n";
		# # system("valgrind rate4site.3.2.mysource_onlycleaned_$mode/sourceMar09/rate4site -bn -a $ensp -t $treefile -s $fastafile -x $treeoutfile -o $normalizedoutfile -y $rate4siteoutfile >& $logfile");
		# if ($mode eq 'fast')
		# {
		# 	system("valgrind rate4site_sources_3.3_fast/programs/rate4site/rate4site -bn -a $ensp -t $treefile -s $fastafile -x $treeoutfile -o $normalizedoutfile -y $rate4siteoutfile >& $logfile");
		# }
		# elsif ($mode eq 'slow')
		# {
		# 	system("valgrind rate4site_sources_3.3_slow/programs/rate4site/rate4site.doubleRep -bn -a $ensp -t $treefile -s $fastafile -x $treeoutfile -o $normalizedoutfile -y $rate4siteoutfile >& $logfile");
		# }
		# 
		# print "       >> ";
		# stoptime(1);
		# 
		# if (`grep -c "likelihoodComputation::getLofPos: likelihood of pos was zero!" $logfile` > 0)
		# {
		# 	print "         >> rate4site crashed! ('likelihood of pos was zero')\n";
		# 	next;
		# }
		# 
		# if (`grep -c "Amino acid was not one of the following" $logfile` > 0)
		# {
		# 	print "         >> rate4site crashed! ('Amino acid was not one of the following')\n";
		# 	next;
		# }
	}
}
stopme() if (!switch('debug'));

showme("added ensps to tasks");
print "\n";
showme("no names in uniens for ensp", 1);
showme("couldn't find a uniprot name with matching sequence for ensp", 1);
print "\n";
showme("task 1", 1);
showme("task 2", 1);
showme("task 3", 1);
showme("task 4", 1);
showme("task 5", 1);
showme("task 6", 1);
showme("task 7", 1);
showme("task 8", 1);
showme("task extra", 1);

done();