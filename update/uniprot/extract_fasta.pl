#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [species]";
($species) = args(1);

# $infile = "input.txt";
$outfile = "output/$species.txt";

# open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");


# start

startme("Extracting all 'UniProt' type sequences for species '$species' from table 'uniseq' and writing output to '$outfile'");
starttime();
# $query = Query("SELECT name, seq FROM uniseq WHERE species='$species' AND type='UniProt' ORDER BY name");
$query = Query("SELECT name, seq FROM uniseq WHERE species='$species' AND type='UniProt'");
die("Error: Numrows ".Numrows($query)." for species '$species'") if (Numrows($query) == 0);
while (($title, $seq) = Fetch($query))
{
	print OUT ">$title\n".split60($seq)."\n";
	stepme(1000);
}
stopme();
stoptime();

done();
