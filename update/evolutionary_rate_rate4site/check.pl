#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
our $usage = "$0 [species] [type]\n\nExample: $0 human linsi_tree";
($species, $type) = args(2);

$outfile = "output-check.txt";
open(OUT, ">$outfile") or die("Error: Couldn't open '$outfile'");

chdir("tmp");

# start

# Get Ensembl species 
$query = Query("SELECT species FROM comparaspecies WHERE unispec='$species'");
($comparaspec) = FetchOne($query);
state("UniProt species '$species' >> Ensembl species '$comparaspec'");

# Output format:
# Alignment filename <TAB> Number of sequences <TAB> Number of $species sequences <TAB> Alignment length <TAB> Number of $species sequences without evolutionary rates <NEWLINE> 

startme("Checking output in table 'evorate_rate4site_$type'");
starttime();
@allensps = ();
@numseq = ();
@seqlengths = ();
while (-e sprintf("%06d", getme() + 1).".fasta.txt")
{
	@ensps = @theseseqlengths = @enspfails = ();
	
	# Get $species ENSP IDs from alignment
	$infile = sprintf("%06d", getme() + 1).".fasta.txt";
	open(IN, "$infile") or die("\nError: Couldn't open '$infile'\n\n");
	fastabreak();
	while (<IN>)
	{
		($title, $seq) = getfasta();
		
		# if ($title =~ /^(ENSP\d+)$/)
		# {
		# 	push(@ensps, $1);
		# }

		$query = Query("SELECT species FROM comparaenspspecies WHERE ensp='$title'");
		($tmpspec) = FetchOne($query);
		# If this is the desired species..
		if ($tmpspec eq $comparaspec)
		{
			push(@ensps, $title);
		}
		push(@theseseqlengths, length($seq));
	}
	normalbreak();
	close(IN);

	# Skip this alignment if it doesn't contain any $species sequences
	if (scalar(@ensps) == 0)
	{
		addme("no $species sequences", $infile);
		stepme(100);
		next;
	}

	# print " >> $infile >> ".scalar(@ensps)." ENSPs\n";
	
	die("Error: ENSP found multiple times in alignment file '$infile'") if (scalar(@ensps) != scalar(unique(@ensps)));
	
	push(@numseq, scalar(@theseseqlengths));
	print OUT "$infile\t".scalar(@theseseqlengths)."\t".scalar(@ensps)."\t";
	
	@theseseqlengths = unique(@theseseqlengths);
	die("Error: Sequence lengths not equal in '$infile'") if (scalar(@theseseqlengths) != 1);
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
		# $tmpfile = "../output/".sprintf("%06d", getme() + 1).".$ensp.Capra$window.txt";
		$tmpfile = "../output/$ensp.txt";
		
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
		
		# $query = Query("SELECT COUNT(rate) FROM evorate WHERE ensp='$ensp' AND type='$type' GROUP BY ensp");
		$query = Query("SELECT COUNT(rate) FROM `evorate_rate4site_$type` WHERE ensp='$ensp'");
		if (Numrows($query) == 0)
		{
			push(@enspfails, $ensp);
		}
		while (($count) = Fetch($query))
		{
			if ($count > 0)
			{
				addme("enspsuccessdb", $ensp);
				addme("alignsuccessdb", $infile);
			}
		}
	}
	print OUT scalar(@enspfails)."\n";
	
	stepme(100);
}
stopme();
stoptime();

state("Wrote to '$outfile'");

showme("no $species sequences", !switch('debug'));
showme("present in multiple alignments", !switch('debug'));
showme("output file not found", !switch('debug'));
showme("output files incomplete for alignment", !switch('debug'));

print "\nSequences per alignment:\n";
characterise(@numseq);

print "\nLength of sequences:\n";
characterise(@seqlengths);

showme("enspsuccessfile", !switch('debug'));
showme("alignsuccessfile", !switch('debug'));

showme("enspsuccessdb", !switch('debug'));
showme("alignsuccessdb", !switch('debug'));

die("ERROR: 'present in multiple alignments' isn't zero") if (countme("present in multiple alignments") > 0);

done();
