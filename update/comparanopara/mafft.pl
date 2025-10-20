#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [species] [MAFFT mode]\n\n Species: Reference species.\n MAFFT mode: linsi, einsi or ginsi (for L-INS-i, E-INS-i or G-INS-i)\n\nExample: $0 human linsi";
($unispec, $mafftmode) = args(2);



# start

# Translate UniProt species mnemonic to full Ensembl species name
$query = Query("SELECT species FROM comparaspecies WHERE unispec='$unispec'");
($species) = FetchOne($query);

# $mainquery = Query("SELECT c.ensp, c.aln FROM comparafasta c, uniens e, unimod m WHERE c.species='$species' AND e.ensp=c.ensp AND m.name=e.name AND m.ptm!='' GROUP BY c.ensp, c.aln ORDER BY c.ensp, c.aln");
# $mainquery = Query("SELECT c.ensp, c.aln FROM comparafasta c, uniens e WHERE c.species='$species' AND e.ensp=c.ensp GROUP BY c.ensp, c.aln ORDER BY c.ensp, c.aln");
# $mainquery = Query("SELECT c.ensp, c.aln FROM comparafasta c, uniens e WHERE c.species='$species' AND e.ensp=c.ensp AND e.ensp='ENSP00000199448' GROUP BY c.ensp, c.aln ORDER BY c.ensp, c.aln");
$mainquery = Query("SELECT DISTINCT ensp1 FROM comparahomology WHERE species1='$unispec' AND homology='ortholog_one2one' AND hc=1 AND ensp1 IS NOT NULL AND ensp2 IS NOT NULL ORDER BY ensp1");
startme("Aligning one2one orthologs for '$unispec' ($species) using MAFFT", 0, Numrows($mainquery));
starttime();

$free = 0;
while (($ensp) = Fetch($mainquery))
{
	stepme(1);

	# Open input file for this ENSP
	$infile = "input/$ensp.txt";
	
	# if input FASTA file doesn't exist or is zero
	if (!-s $infile)
	{
		addme("couldn't find FASTA file for ensp (probably because it doesn't have any one-to-one orthologs in the orthologue mart files) (skipped)", $ensp);
		next;
	}
	
	# Set MAFFT output file names
	$outfile = "output/mafft.$ensp.$mafftmode.txt";

	if (-s $outfile)
	{
		# Skip this ensp if output file already exists and is non-zero
		addme("alignment file already exists for ensp (skipped)", $ensp);
		next;
		# # Crash instead
		# die("Error: Alignment file '$outfile' already exists for ensp '$ensp' (should run clean.pl first)");
	}

	# Sequences all verified, run MAFFT!
	chdir("tmp"); # make sure log files end up in ./tmp/
	
	# Wait for grid nodes to free up
	while ($free <= 0)
	{
		$free = freenodes();
	}

	# Run MAFFT
	# run("MAFFT L-INS-i", "~/scripts/qsub.sh mafft --localpair --maxiterate 16 --reorder ../$infile \\> ../$outfile", 1);
	# run("MAFFT $mafftmode", "~/scripts/qsub.sh mafft --localpair --maxiterate 1000 --anysymbol ../$infile \\> ../$outfile", 1);

	state("Running MAFFT:", 1) if (switch('debug'));
	addme("ran MAFFT for ensp", $ensp);
    run("MAFFT $mafftmode", "~/scripts/qsub.sh mafft-$mafftmode --quiet --anysymbol ../$infile \\> ../$outfile", 1);
	$free--;
	
	# Done!
	chdir("..");
}
stopme();
stoptime();

# showmeall();
showmesome(50);

done();
