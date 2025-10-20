#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
$maxbranch = '100.0000';

our $usage = "$0 [source]\n\nExample: $0 linsi_tree";
($source) = args(1);

if ($source eq 'para')
{
	$fastafile = "../compara/input/Compara.protein_trees.aa.fasta.gz";
	$treefile = "../compara/input/Compara.protein_trees.nh.emf.gz";
}
else
{
	$fastafile = "../comparanopara/output_".$source."/Compara.protein_trees.aa.fasta.gz";
	$treefile = "../comparanopara/output_".$source."/Compara.protein_trees.nh.emf.gz";
}

open(FASTA, "zcat $fastafile|") or die("\nError: Couldn't zcat '$fastafile'\n\n");
open(TREES, "zcat $treefile|") or die("\nError: Couldn't zcat '$treefile'\n\n");


# start

startme("Splitting '$fastafile' into one-alignment chunks");
starttime();
$/ = "\n//\n\n";
%f = ();
@fastafiles = ();
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

	# $outfile = "tmp/".sprintf("%06d", getme()).".fasta.raw.txt";
	$outfile = "tmp/".sprintf("%06d", getme()).".fasta.txt";
	push(@fastafiles, $outfile);
	open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");
	print OUT $_;
	close(OUT);
}
$/ = "\n";
close(FASTA);
stopme();
stoptime();


# Not necessary anymore since all sequences now only contain the standard 20 amino acids
# # Reopen FASTA files to replace U (selenocysteine) with X
# foreach $infile (@fastafiles)
# {
#     $outfile = $infile;
#     $outfile =~ s/\.raw//;
#     
#     fastabreak();
# 	open(IN, "$infile") or die("\nError: Couldn't open '$infile'\n\n");
# 	open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");
#     while (<IN>)
#     {
#         ($title, $seq) = getfasta();
#     
#         $seq =~ s/U/X/g;
#     
#         print OUT ">$title\n".split60($seq)."\n";
#     }
# 	close(IN);
# 	close(OUT);
#     normalbreak();
# }



startme("Splitting '$treefile' into one-alignment chunks and setting branch lengths to a maximum of $maxbranch (will print a warning if this happens)");
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
	
	# Replace branch lengths greater than $maxbranch with $maxbranch
	while (/:(\d+\.\d+)/g)
	{
		if ($1 > $maxbranch)
		{
			warn("Warning: Branch $1 > $maxbranch for aln ".getme()." in '$treefile' (replaced with $maxbranch)");
			s/:$1/:$maxbranch/;
		}
	}
	
	$outfile = "tmp/".sprintf("%06d", getme()).".tree.txt";
	open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");
	print OUT $_;
	close(OUT);
}
$/ = "\n";
close(TREES);
stopme();
stoptime();

done();
