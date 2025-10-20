#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$nseqmin = 10;

our $usage = "$0 [MAFFT mode] [-tree] [-1para]\n\n -tree: Use MAFFT results with species tree.\n -1para: Keeps the best-matching outparalog (by various metrics) for one2many and many2many cases.\n\nExample: $0 linsi -tree (will test 'output_linsi_tree')";
($mafftmode) = args(1);
if (switch('tree'))
{
	$mafftmode .= "_tree";
}

# Handling for 1para
if (switch('1parahc'))
{
	$mafftmode .= "_1parahc";
}
elsif (switch('1para'))
{
	$mafftmode .= "_1para";
}

$fastafile = "output_$mafftmode/Compara.protein_trees.aa.fasta.gz";
$treefile = "output_$mafftmode/Compara.protein_trees.nh.emf.gz";

open(FASTA, "zcat $fastafile|") or die("\nError: Couldn't zcat '$fastafile'\n\n");
open(TREES, "zcat $treefile|") or die("\nError: Couldn't zcat '$treefile'\n\n");


# start

startme("Testing '$fastafile'");
starttime();
$/ = "\n//\n\n";
%f = ();
# $first = 1;
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

	# $outfile = "tmp/".sprintf("%06d", getme()).".fasta.txt";
	# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");
	# print OUT $_;
	# close(OUT);

	# if ($first == 1)
	# {
	# 	print $_;
	# 	$first = 0;
	# }
}
# print $_;
close(FASTA);
stopme();
stoptime();

startme("Testing '$treefile'");
starttime();
$/ = "//\n\n";
%t = ();
# $first = 1;
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
	
	$query = Query("SELECT species FROM comparaenspspecies WHERE ensp IN ('".join("', '", split(/\|/, $t{getme()}))."')");
	$human = 0;
	@species = ();
	while (($species) = Fetch($query))
	{
		if ($species eq 'homo_sapiens')
		{
			$human = 1;
		}
		push(@species, $species);
	}
	
	die("Error: Species duplicates for aln '".getme()."'") if (scalar(@species) != scalar(unique(@species)));
	die("Error: Sequence duplicates for aln '".getme()."'") if (scalar(split(/\|/, $t{getme()})) != scalar(unique(split(/\|/, $t{getme()}))));
	
	if (scalar(unique(split(/\|/, $t{getme()}))) < $nseqmin)
	{
		addme("aln contains less than $nseqmin sequences", getme());
	}
	else
	{
		addme("aln contains $nseqmin or more sequences", getme());
	}

	if (scalar(unique(@species)) < $nseqmin)
	{
		addme("aln contains less than $nseqmin species", getme());
	}
	else
	{
		addme("aln contains $nseqmin or more species", getme());
	}
	
	if ($human == 1)
	{
		addme("aln contains human", getme());
	}
	else
	{
		addme("aln doesn't contain human", getme());
	}

	s/^SEQ.+?\n//msg;
	s/^DATA\n//msg;
	
	# $outfile = "tmp/".sprintf("%06d", getme()).".tree.txt";
	# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");
	# print OUT $_;
	# close(OUT);
	
	# if ($first == 1)
	# {
	# 	print $_;
	# 	$first = 0;
	# }
}
close(TREES);
stopme();
stoptime();

showmeall(1);

done();
