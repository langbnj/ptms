#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
our $usage = "$0 [MAFFT mode (linsi/ginsi/einsi)] [-tree]\n\nTakes its input from e.g. ../comparanopara/output_linsi_tree/Compara.protein_trees.nh.emf\n\nExample: $0 linsi -tree";
($mode) = args(1);
if (switch('tree')) { $mode .= "_tree"; }
if (switch('1para')) { $mode .= "_1para"; }

# Directory
$indir = "../comparanopara/output_".$mode;
if (!-d $indir)
{
	die("Error: Couldn't find directory '$indir'");
}


# Open input files
$fastafile = $indir."/Compara.protein_trees.aa.fasta.gz";
$treefile  = $indir."/Compara.protein_trees.nh.emf.gz";

open(FASTA, "zcat $fastafile|") or die("\nError: Couldn't zcat '$fastafile'\n\n");
open(TREES, "zcat $treefile|") or die("\nError: Couldn't zcat '$treefile'\n\n");


# start

startme("Splitting '$fastafile' into one-alignment chunks");
starttime();
$/ = "\n//\n\n";
%f = ();
while (<FASTA>)
{
	chomp;
	
	stepme(100);

	# replace spaces at the end of FASTA titles
	s/ +$//mg;

	# list FASTA titles
	$f{getme()} = '';
	while (/^>(.+)/mg)
	{
		$f{getme()} .= "$1|";
	}

	$outfile = "tmp/$mode.".sprintf("%06d", getme()).".fasta.txt";
	open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");
	print OUT $_;
	close(OUT);
}
close(FASTA);
stopme();
stoptime();

startme("Splitting '$treefile' into one-alignment chunks");
starttime();
$/ = "//\n\n";
%t = ();
while (<TREES>)
{
	chomp;
	
	stepme(100);

	# list FASTA titles
	$t{getme()} = '';
	while (/^SEQ \S+ (\S+) /mg)
	{
		$t{getme()} .= "$1|";
	}
	
	$f{getme()} = join("|", unique(split(/\|/, $f{getme()})));
	$t{getme()} = join("|", unique(split(/\|/, $t{getme()})));
	
	if ($f{getme()} ne $t{getme()})
	{
		die("Error: Proteins not the same in chunk ".getme().":\n\nFASTA\n".$f{getme()}."\nvs.\nTree\n".$t{getme()}."\n\n");
	}

	s/^SEQ.+?\n//msg;
	s/^DATA\n//msg;
	
	$outfile = "tmp/$mode.".sprintf("%06d", getme()).".tree.txt";
	open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");
	print OUT $_;
	close(OUT);
}
close(TREES);
stopme();
stoptime();

done();
