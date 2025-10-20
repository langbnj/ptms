#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [mode: e.g. linsi, linsi_tree] [species: e.g. human or 'all']\n\nExample: $0 linsi_tree all\nExample: $0 linsi_tree human";
($mode, $species) = args(2);
# args(0);

$seqtable = 'comparafasta_'.$mode;
# $cdstable = 'comparacds'.$mode;		# There are no comparacds_linsi_tree etc. tables (only comparacds itself), therefore I'll get the CDS from the PAL2NAL output files

# $infile = "input.txt";
# $outfile = "output.txt";

# open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");



# Get CDSs from PAL2NAL output files
# e.g. cds.ENSP00000412673.linsi.tree.txt
$tmpmode = $mode;
$tmpmode =~ s/_/./g;	# Translate linsi_tree to linsi.tree
startme("Reading PAL2NAL CDS output files", 0, chompme(`find output -name 'cds.*.$tmpmode.txt' | wc -l`));
open(FIND, "find output -name 'cds.*.$tmpmode.txt'|") or die("Error: Couldn't list PAL2NAL output files (cds.*.$mode.txt)");
%cds = ();
while (<FIND>)
{
	chomp;
	$infile = $_;
	open(IN, $infile) or die("Error: Couldn't open '$infile'");
	fastabreak();
	while (<IN>)
	{
		($ensp, $cds) = getfasta();

		if ($species eq 'human')
		{
			next if ($ensp !~ /^ENSP\d+$/);
		}

		$cds{$ensp} = $cds;

		if ($species eq 'human')
		{
			last if ($ensp !~ /^ENSP\d+$/);
		}

		# stepme(1000);
	}
	normalbreak();
	close(IN);

	stepme(100);
}
stopme();
close(FIND);


# start

if ($species eq 'all')
{
	# All species
	$mainquery = Query("SELECT ensp, REPLACE(seq, '-', '') FROM comparafasta");
}
else
{
	# Translate unispec to fullspecies
	$query = Query("SELECT species FROM comparaspecies WHERE unispec='$species'");
	($species) = FetchOne($query);

	# Specific species
	$mainquery = Query("SELECT ensp, REPLACE(seq, '-', '') FROM comparafasta WHERE species='$species'");
}
startme("Getting AA sequences from table '$seqtable' and CDS from PAL2NAL output files for species '$species' and checking whether the CDS translates to the AA sequence", 0, Numrows($mainquery));
starttime();
$match = 0;
$mismatch = 0;
while (($ensp, $seq) = Fetch($mainquery))
{
	stepme(1000);

	# Get CDS from PAL2NAL output file
	if (!exists($cds{$ensp}))
	{
		# die("Error: No CDS for ENSP '$ensp'");
		addme("no cds in PAL2NAL output files for ensp (skipped)", $ensp);
		next;
	}
	$cds = $cds{$ensp};

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
