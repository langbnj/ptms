#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

# start

startme("Checking how many 'human' modified proteins from 'unimod' have the proper amino acids in their 'uniseq' sequences");
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
	
	# Get PTMs
	$ptmquery = Query("SELECT site, ptm FROM unimod WHERE acc='$acc' AND species='human' AND ptm!='' ORDER BY site, ptm");
	if (Numrows($ptmquery) == 0)
	{
		addme("no PTMs in unimod for acc", $acc);
		next;
	}
	
	$onematch = 0;
	while (($site, $ptm) = Fetch($ptmquery))
	{
		# Get expected amino acid
		$ptm =~ /^(\w)-.+/;
		$expected = $1;
		
		# Get site amino acid
		if ($site > length($seq))
		{
			# die("Error: No sequence in Compara for ENSP '$ensp'");
			addme("site outside sequence for acc", $acc);
			addme("site outside sequence for ptm site", "$acc|$site|$ptm");
			next;
		}
		$found = substr($seq, $site-1, 1);
		
		# Compare UniProt vs. Compara
		if ($found eq $expected)
		{
			addme("aa match for acc", $acc);
			addme("aa match for ptm site", "$acc|$site|$ptm");
			$onematch = 1;
		}
		else
		{
			addme("aa mismatch for acc", $acc);
			addme("aa mismatch for ptm site", "$acc|$site|$ptm");
		}
	}
	
	if ($onematch == 1)
	{
		addme("at least one aa match for acc", $acc);
	}
}
stopme();
stoptime();

showme("no PTMs in unimod for acc", 1);
print "\n";
showme("no sequence in uniseq for acc", 1);
showme("site outside sequence for acc", 1);
showme("site outside sequence for ptm site", 1);
print "\n";
showme("aa match for acc", 1);
showme("aa mismatch for acc", 1);
print "\n";
showme("aa match for ptm site", 1);
showme("aa mismatch for ptm site", 1);
print "\n";
showme("at least one aa match for acc", 1);

done();
