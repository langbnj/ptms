#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [species: e.g. human or 'all']\n\nExample: $0 all\nExample: $0 human";
($species) = args(1);
# args(0);

# $infile = "input.txt";
# $outfile = "output.txt";

# open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");


# start

if ($species eq 'all')
{
	# All species
	$mainquery = Query("SELECT cds, seq FROM ensembl");
}
else
{
	# Specific species
	$mainquery = Query("SELECT cds, seq FROM ensembl WHERE species='$species'");
}
startme("Getting CDS and AA sequences from table 'ensembl' and checking whether the CDS translates to the AA sequence", 0, Numrows($mainquery));
starttime();
$match = 0;
$mismatch = 0;
while (($cds, $seq) = Fetch($mainquery))
{
	addme("total expected seqs", $seq);
	stepme(1000);

	# Verify that seq is fine (nothing except standard 20 AAs, i.e. no B/Z/U/X/*/etc.)
	if (!aa($seq))
	{
		addme("skipped sequence because it contained characters other than the standard 20 AA for seq", $seq);
		next;
	}

	# Verify that cds is fine (no thing except ACGT)
	if (!dna($cds))
	{
		addme("skipped sequence because it contained characters other than ACGT for seq", $seq);
		next;
	}

	# # Add stop codon to expected sequence (the CDS usually includes this)
	# $seq .= '*';
	# # Remove terminal stop codon from CDS (the CDS doesn't always include this, and the expected sequence never does)
	# $cds =~ s/(TAA|TGA|TAG)$//;

	# # Replace selenocysteine (U) with its encoding stop codon (*) in expected sequence
	# $seq =~ s/U/*/g;

	# Verify CDS
# 	# Check if CDS length is a multiple of 3
# 	if (length($cds) % 3 != 0)
# 	{
# 		# addme("cds length isn't a multiple of 3 for expected seq (skipped)", $seq);
# 		# $mismatch++;
# 		# next;
# 		$overhang = length($cds) - (length($seq) * 3);
# 		addme("cds length isn't a multiple of 3 for expected seq (trimmed)", $seq);
# 		addme("cds length isn't a multiple of 3 for expected seq (trimmed by $overhang)", $seq);
# 
# 		$cds = substr($cds, 0, length($seq) * 3);
# 	}
	# Check if CDS length divided by 3 matches expected seq length
	if (length($seq) != length($cds) / 3)
	{
		# addme("cds length / 3 doesn't match expected seq length for expected seq (skipped)", $seq);
		# $mismatch++;
		# next;

		$overhang = length($cds) - (length($seq) * 3);
		addme("cds length / 3 doesn't match expected seq length for expected seq (trimmed)", $seq);
		addme("cds length / 3 doesn't match expected seq length for expected seq (trimmed by $overhang)", $seq);

		$cds = substr($cds, 0, length($seq) * 3);
	}

	# Translate CDS
	# $transeq = translate_loosely($cds);
	$transeq = translate($cds);

	# # Remove terminal stop codon from translation (the CDS doesn't always include this, and the expected sequence never does)
	# $transeq =~ s/\*$//;

	# Check if CDS translation matches expected seq
	if ($transeq ne $seq)
	{
		addme("translation mismatch for expected seq", $seq);
		$mismatch++;
	}
	else
	{
		addme("translation match for expected seq", $seq);
		$match++;
	}
}
stopme();
stoptime();

showmeall(1);

state(commify($match)." matches, ".commify($mismatch)." mismatches");

done();
