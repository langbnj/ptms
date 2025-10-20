#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$table = 'unimod_control';

our $usage = "$0 [species] [PTMs to look at] [PTMs to avoid] [window]\n\n".
" PTMs to look at:	'all' / 'hc' / 'K-ac' / ... (all: all PTMs including potential/probable/by similarity, hc: only PTMs where 'ptm' field is set)\n".
" PTMs to avoid:	'all' / 'hc' / 'K-ac' / ... (all: all PTMs including potential/probable/by similarity, hc: only PTMs where 'ptm' field is set)\n".
" Window: 			Stay more than this far away from any PTMs being avoided when picking control sites.\n\n".
"Example: $0 human hc all 0";
($species, $ptm, $avoid, $window) = args(4);

$species = uc($species);

# Clear($table);
Query("DELETE FROM `$table` WHERE species='$species'");
state("Cleared '$species' entries from $table");

# $infile = "input.txt";
#$outfile = "output.txt";

# open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
#open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");


state("Getting relevant amino acids for all '$species' PTM sites of type '$ptm':");

if ($ptm eq 'all')
{
	# All PTMs (including probable/potential/by similarity) ('all')
	$mainquery = Query("SELECT DISTINCT(SUBSTRING(s.seq, m.site, 1)) FROM unimod m, uniseq s
						WHERE m.species='$species' AND s.acc=m.acc AND s.type='UniProt'");
}
elsif ($ptm eq 'hc')
{
	# Experimental PTMs only ('hc')
	$mainquery = Query("SELECT DISTINCT(SUBSTRING(s.seq, m.site, 1)) FROM unimod m, uniseq s
						WHERE m.species='$species' AND m.ptm IS NOT NULL AND s.acc=m.acc AND s.type='UniProt'");
}
else
{
	# Specific experimental PTM only (e.g. 'K-ac')
	$mainquery = Query("DISTINCT(SUBSTRING(s.seq, m.site, 1))  FROM unimod m, uniseq s
						WHERE m.species='$species' AND m.ptm='$ptm' AND s.acc=m.acc AND s.type='UniProt'");
}
@aas = ();
while (($aa) = Fetch($mainquery))
{
	push(@aas, $aa);
}
@aas = unique(@aas);
print join('', @aas)."\n\n";


# Make temporary table for control sites
# startme("Filling table '$table' with control sites for all '$species' proteins with PTM sites of type '$ptm', while avoiding 'all' PTMs by at least '$window'");
startme("Filling table '$table' with control sites for all '$species' proteins with PTM sites of type '$ptm', while avoiding '$avoid' PTMs by at least '$window'");
startme2();
starttime();
if ($ptm eq 'all')
{
	# All PTMs (including probable/potential/by similarity) ('all')
	$ptmstr = "";
}
elsif ($ptm eq 'hc')
{
	# Experimental PTMs only ('hc')
	$ptmstr = "AND m.ptm IS NOT NULL ";
}
else
{
	# Specific experimental PTM only (e.g. 'K-ac')
	$ptmstr = "AND m.ptm='$ptm' ";
}

if ($avoid eq 'all')
{
	# All PTMs (including probable/potential/by similarity) ('all')
	$avoidstr = "";
}
elsif ($avoid eq 'hc')
{
	# Experimental PTMs only ('hc')
	$avoidstr = "AND n.ptm IS NOT NULL ";
}
else
{
	# Specific experimental PTM only (e.g. 'K-ac')
	$avoidstr = "AND n.ptm='$avoid' ";
}

# Limiting factor:
# only modified proteins (from unimod)
$mainquery = Query("SELECT m.name, m.acc, s.seq, GROUP_CONCAT(DISTINCT n.site ORDER BY n.site SEPARATOR '|') FROM unimod m
JOIN uniseq s ON s.acc=m.acc
RIGHT JOIN unimod n ON n.acc=m.acc
WHERE m.species='$species' ".$ptmstr.$avoidstr."AND s.acc=m.acc AND s.type='UniProt' GROUP BY m.acc");
while (($name, $acc, $seq, $ptmsites) = Fetch($mainquery))
{
	@ptmsites = split(/\|/, $ptmsites);
	foreach $site (positions('['.join('', @aas).']', $seq))
	{
		@avoidsites = ();
		$i = -$window;
		print " >> $site\n";
		while ($i <= $window)
		{
			if ((($site + $i) >= 1) and (($site + $i) <= length($seq)))
			{
				push(@avoidsites, $site + $i);
				# print "   >> ".($site + $i)."\n";
			}
			$i++;
		}

		# "AVOID" are the sites currently being looked at, not necessarily to be avoided yet (only if they overlap with PTM sites, including their +/- $window window)
		# print "AVOID=[".join(",", @avoidsites)."]\n";
		# print "PTM  =[".join(",", @ptmsites)."]\n";
		# print "INTER=[".join(",", intersection(\@avoidsites, \@ptmsites))."]\n";
		# # exit;
		
		if (scalar(intersection(\@avoidsites, \@ptmsites)) == 0)
		{
			Query("INSERT INTO `$table` SET name='$name', acc='$acc', species='$species', site='$site', aa='".substr($seq, $site-1, 1)."'");
			stepme2();
		}
	}
	stepme(100);
}
stopme();
stoptime();

state("Inserted ".getme2()." control sites for ".getme()." proteins");

done();
