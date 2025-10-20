#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "";
args(0);


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

startme("Checking how many 'human' modified proteins from 'unimod' have identical 'uniseq' and Compara sequences");
starttime();

$accquery = Query("SELECT acc FROM unimod WHERE ptm!='' AND species='human' GROUP BY acc ORDER BY acc");
while (($acc) = Fetch($accquery))
{
	stepme(100);
	
	# Get UniProt sequence
	$query = Query("SELECT DISTINCT seq FROM uniseq WHERE acc='$acc' AND species='human' AND type IN ('UniProt', 'UniIso')");
	if (Numrows($query) == 0)
	{
		# die("Error: No sequence in 'uniseq' for acc '$acc'");
		addme("no sequence in uniseq for acc", $acc);
		next;
	}
	($seq) = FetchOne($query);
	
	# Get ENSPs
	$enspquery = Query("SELECT ensp FROM uniens WHERE acc='$acc' AND species='human'");
	if (Numrows($enspquery) == 0)
	{
		addme("no ensp in uniens for acc", $acc);
		next;
	}
	
	$onematch = 0;
	while (($ensp) = Fetch($enspquery))
	{
		# Get Compara sequence
		if (!exists($compara{$ensp}))
		{
			# die("Error: No sequence in Compara for ENSP '$ensp'");
			addme("no sequence in compara for ensp", $ensp);
			next;
		}
		
		# Compare UniProt vs. Compara
		if ($seq eq $compara{$ensp})
		{
			addme("sequence match for acc", $acc);
			addme("sequence match for ensp", $ensp);
			$onematch = 1;
		}
		else
		{
			addme("sequence mismatch for acc", $acc);
			addme("sequence mismatch for ensp", $ensp);
		}
	}
	
	if ($onematch == 1)
	{
		addme("at least one sequence match for acc", $acc);
	}
}
stopme();
stoptime();

showme("no ensp in uniens for acc", !switch('debug'));
print "\n";
showme("no sequence in uniseq for acc", !switch('debug'));
showme("no sequence in compara for ensp", !switch('debug'));
print "\n";
showme("sequence match for acc", !switch('debug'));
showme("sequence mismatch for acc", !switch('debug'));
print "\n";
showme("sequence match for ensp", !switch('debug'));
showme("sequence mismatch for ensp", !switch('debug'));
print "\n";
showme("at least one sequence match for acc", !switch('debug'));

done();
