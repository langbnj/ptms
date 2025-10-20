#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
our $usage = "$0 [mode: e.g. linsi_tree, ginsi] [window size]\n\nExample: $0 linsi_tree 0";
($mode, $window) = args(2);

$table = 'evorate_capra'.$window.'_'.$mode;

$outfile = "output-check-$mode-$window.txt";
open(OUT, ">$outfile") or die("Error: Couldn't open '$outfile'");

chdir("tmp");

# start

startme("Checking mode '$mode' window '$window' output and writing to '$outfile'");
starttime();
@allensps = ();
@numseq = ();
@seqlengths = ();
while (-e $mode.".".sprintf("%06d", getme() + 1).".fasta.txt")
{
	@ensps = @theseseqlengths = @enspfails = ();
	%ensplen = ();
	
	# Get human ENSP IDs from alignment
	$infile = $mode.".".sprintf("%06d", getme() + 1).".fasta.txt";
	open(IN, "$infile") or die("\nError: Couldn't open '$infile'\n\n");
	fastabreak();
	while (<IN>)
	{
		($title, $seq) = getfasta();

		if ($title =~ /^(ENSP\d+)$/)
		{
			$ensp = $1;
			push(@ensps, $ensp);

			# Get length of this ENSP's sequence (without gaps)
			$tmpseq = $seq;
			$tmpseq =~ s/-//g;
			
			$ensplen{$ensp} = length($tmpseq);
		}

		push(@theseseqlengths, length($seq));
	}
	normalbreak();
	close(IN);

	# Skip this alignment if it doesn't contain any human sequences
	if (scalar(@ensps) == 0)
	{
		addme("no human sequences", $infile);
		stepme(100);
		next;
	}

	# print " >> $infile >> ".scalar(@ensps)." ENSPs\n";
	
	die("Error: ENSP found multiple times in alignment file '$infile'") if (scalar(@ensps) != scalar(unique(@ensps)));
	
	push(@numseq, scalar(@theseseqlengths));
	print OUT "$infile\t".scalar(@theseseqlengths)."\t".scalar(@ensps)."\t";
	
	@theseseqlengths = unique(@theseseqlengths);
	die("Error: Sequence lengths not equal in '$infile'\n\n".show(@theseseqlengths)) if (scalar(@theseseqlengths) != 1);
	push(@seqlengths, @theseseqlengths);
	print OUT "$theseseqlengths[0]\t";
	
	
	foreach $ensp (@ensps)
	{
		if (contains($ensp, @allensps))
		{
			addme("present in multiple alignments", $ensp);
		}
		push(@allensps, $ensp);
		
		# print " >> $ensp\n";
		$tmpfile = "../output/$mode.".sprintf("%06d", getme() + 1).".$ensp.capra$window.txt";
		# $tmpfile = "../output/$ensp.txt";
		
		if (!-e $tmpfile)
		{
			addme("output file not found", $tmpfile);
			addme("output files incomplete for alignment", $infile);
		}
		else
		{
			addme("enspsuccessfile", $ensp);
			addme("alignsuccessfile", $infile);
		}
		
		$query = Query("SELECT COUNT(rate) FROM $table WHERE ensp='$ensp'");
		if (Numrows($query) == 0)
		{
			push(@enspfails, $ensp);
		}
		while (($count) = Fetch($query))
		{
			if ($count > 0)
			{
				# addme("enspsuccessdb", $ensp);
				# addme("alignsuccessdb", $infile);

				# More sophisticated check (does the evorate table ($table) contain values for each residue, i.e. the complete protein?)
				if ($count == $ensplen{$ensp})
				{
					addme("enspsuccessdb", $ensp);
					addme("alignsuccessdb", $infile);
				}
				else
				{
					die("Error: $table contains $count residues for $ensp, but '$infile' contains ".$ensplen{$ensp});
				}
			}
			else
			{
				addme("no results in table $table for ensp", $ensp);
			}
		}
	}
	print OUT scalar(@enspfails)."\n";
	
	stepme(100);
}
stopme();
stoptime();

state("Wrote to '$outfile'");

# showme("no human sequences", !switch('debug'));
# showme("present in multiple alignments", !switch('debug'));
# showme("output file not found", !switch('debug'));
# showme("output files incomplete for alignment", !switch('debug'));

print "\nSequences per alignment:\n";
characterise(@numseq);

print "\nLength of sequences:\n";
characterise(@seqlengths);

nl();

# showme("enspsuccessfile", !switch('debug'));
# showme("alignsuccessfile", !switch('debug'));
# 
# showme("enspsuccessdb", !switch('debug'));
# showme("alignsuccessdb", !switch('debug'));

showmeall(1);

die("ERROR ERROR ERROR: 'present in multiple alignments' isn't zero") if (countme("present in multiple alignments") > 0);

done();
