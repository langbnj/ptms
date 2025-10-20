#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [species: e.g. human or 'all']\n\nExample: $0 all\nExample: $0 human";
($species) = args(1);
# args(0);

$seqtable = 'comparafasta';
$cdstable = 'comparacds';

# $infile = "input.txt";
# $outfile = "output.txt";

# open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");


# start

if ($species eq 'all')
{
	# All species
	$mainquery = Query("SELECT c.seq AS cds, f.seq AS seq FROM comparafasta f, comparacds c WHERE f.ensp=c.ensp");
}
else
{
	# Translate unispec to fullspecies
	$query = Query("SELECT species FROM comparaspecies WHERE unispec='$species'");
	($species) = FetchOne($query);

	# Specific species
	$mainquery = Query("SELECT REPLACE(c.seq, '-', '') AS cds, REPLACE(f.seq, '-', '') AS seq FROM comparafasta f, comparacds c WHERE f.ensp=c.ensp AND f.species='$species'");
}
startme("Getting CDS and AA sequences from tables '$seqtable' and '$cdstable' for species '$species' and checking whether the CDS translates to the AA sequence", 0, Numrows($mainquery));
starttime();
$match = 0;
$mismatch = 0;
while (($cds, $seq) = Fetch($mainquery))
{
	stepme(1000);

	# # Add stop codon to expected sequence (the CDS includes this)
	# $seq .= '*';

	# Replace selenocysteine (U) with its encoding stop codon (*) in expected sequence
	$seq =~ s/U/*/g;

	# Verify CDS
	# Check if CDS length is a multiple of 3
	if (length($cds) % 3 != 0)
	{
		addme("cds length isn't a multiple of 3 for expected seq (skipped)", $seq);
		$mismatch++;
		next;
	}
	# Check if CDS length divided by 3 matches expected seq length
	if (length($seq) != length($cds) / 3)
	{
		addme("cds length / 3 doesn't match expected seq length for expected seq (skipped)", $seq);
		$mismatch++;
		next;
	}

	# Translate CDS
	$transeq = translate_loosely($cds);

	# Check if CDS translation matches expected seq
	if ($transeq ne $seq)
	{
		addme("translation mismatch for expected seq", $seq);
		$mismatch++;
	}
	else
	{
		$match++;
	}
}
stopme();
stoptime();

showmeall(1);

state(commify($match)." matches, ".commify($mismatch)." mismatches");

done();
