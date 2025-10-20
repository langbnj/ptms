#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

# our $superloudmysql = 1;

$accquery = Query("SELECT acc FROM unimod WHERE ptm!='' AND species='human' GROUP BY acc ORDER BY acc");
startme("Checking how many 'human' modified proteins from 'unimod' have identical 'uniseq' and 'ensembl' sequences", 0, Numrows($accquery));
starttime();

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
		# Get ENSEMBL sequence
		$query = Query("SELECT seq FROM ensembl WHERE ensp='$ensp' AND species='human'");
		if (Numrows($query) == 0)
		{
			# die("Error: No sequence in 'ensembl' for ENSP '$ensp'");
			addme("no sequence in ensembl for ensp", $ensp);
			next;
		}
		($enseq) = FetchOne($query);
		
		# Compare UniProt vs. ENSEMBL
		if ($seq eq $enseq)
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

showme("no ensp in uniens for acc", 1);
print "\n";
showme("no sequence in uniseq for acc", 1);
showme("no sequence in ensembl for ensp", 1);
print "\n";
showme("sequence match for acc", 1);
showme("sequence mismatch for acc", 1);
print "\n";
showme("sequence match for ensp", 1);
showme("sequence mismatch for ensp", 1);
print "\n";
showme("at least one sequence match for acc", 1);

done();
