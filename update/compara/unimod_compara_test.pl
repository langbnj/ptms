#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

# start

startme("Reading human Compara sequences from 'tmp/*.fasta.txt'");
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
		
		# if human sequence: store in hash
		if ($title =~ /^ENSP\d+$/)
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

state("Read ".scalar(keys(%compara))." human sequences");

startme("Checking how many 'human' modified proteins from 'unimod' have the proper amino acids in their Compara sequences");
starttime();

$accquery = Query("SELECT acc FROM unimod WHERE ptm!='' AND species='human' GROUP BY acc ORDER BY acc");
while (($acc) = Fetch($accquery))
{
	stepme(100);
	
	# Get ENSPs
	$enspquery = Query("SELECT ensp FROM uniens WHERE acc='$acc' AND species='human'");
	if (Numrows($enspquery) == 0)
	{
		addme("no ensp in uniens for acc", $acc);
		next;
	}
	
	while (($ensp) = Fetch($enspquery))
	{
		@matches = ();
		
		# Get PTMs
		$ptmquery = Query("SELECT site, ptm FROM unimod WHERE acc='$acc' AND species='human' AND ptm!='' ORDER BY site, ptm");
		if (Numrows($ptmquery) == 0)
		{
			addme("no PTMs in unimod for acc", $acc);
			next;
		}
		while (($site, $ptm) = Fetch($ptmquery))
		{
			# Get expected amino acid
			$ptm =~ /^(\w)-.+/;
			$expected = $1;
			
			# Get Compara sequence
			if (!exists($compara{$ensp}))
			{
				# die("Error: No sequence in Compara for ENSP '$ensp'");
				addme("no sequence in compara for acc", $acc);
				addme("no sequence in compara for ensp", $ensp);
				next;
			}
			
			# Get site amino acid
			if ($site > length($compara{$ensp}))
			{
				# die("Error: No sequence in Compara for ENSP '$ensp'");
				addme("site outside sequence for acc", $acc);
				addme("site outside sequence for ptm site", "$acc|$site|$ptm");
				next;
			}
			$found = substr($compara{$ensp}, $site-1, 1);
			
			# Compare UniProt vs. Compara
			if ($found eq $expected)
			{
				addme("aa match for acc", $acc);
				addme("aa match for ptm site", "$acc|$site|$ptm");
				push(@matches, 1);
			}
			else
			{
				addme("aa mismatch for acc", $acc);
				addme("aa mismatch for ptm site", "$acc|$site|$ptm");
				push(@matches, 0);
			}
		}
		if ((scalar(@matches) > 0) and (min(@matches) == 1))
		{
			addme("at least one ensp matches perfectly for acc", $acc);
		}
	}
}
stopme();
stoptime();

showme("no PTMs in unimod for acc", 1);
print "\n";
showme("no ensp in uniens for acc", 1);
showme("no sequence in compara for acc", 1);
showme("no sequence in compara for ensp", 1);
print "\n";
showme("site outside sequence for acc", 1);
showme("site outside sequence for ptm site", 1);
print "\n";
showme("aa match for acc", 1);
showme("aa mismatch for acc", 1);
print "\n";
showme("aa match for ptm site", 1);
showme("aa mismatch for ptm site", 1);
print "\n";
showme("at least one ensp matches perfectly for acc", 1);

done();
