#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

# $speciestable = 'comparaspecies';
$enspspeciestable = 'comparaenspspecies';
$fastatable = 'comparafasta';
$cdstable = 'comparacds';
$nhtable = 'comparanh';
$nhxtable = 'comparanhx';

our $usage = "$0 [Compara release number]\nExample: $0 107";
($rel) = args(1);

# start
$treefile = "input/Compara.$rel.protein_default.nh.emf.gz";
open(TREES, "zcat $treefile|") or die("\nError: Couldn't zcat '$treefile'\n\n");

Clear($enspspeciestable);
Clear($fastatable);
Clear($cdstable);
Clear($nhtable);
Clear($nhxtable);

# Get Compara species list:
# zcat input/Compara.protein_trees.nh.emf.gz | perl -ne '/^SEQ (\S+) /; print "$1\n";' | sort | uniq
# Longest species name:
# zcat input/Compara.protein_trees.nh.emf.gz | perl -ne '/^SEQ (\S+) /; print "$1\n";' | sort | uniq | wc -L


startme("Inserting Compara ENSP -> Species mappings from '$treefile' into table '$enspspeciestable'");
starttime();
while (<TREES>)
{
	next if (!/^SEQ/);
	
	/^SEQ (\S+) (\S+)/ or die("Error: Couldn't parse Species and ENSP from line '$_'");
	
	$species = $1;
	$ensp = $2;
	
	Query("INSERT INTO `$enspspeciestable` SET species='$species', ensp='$ensp'");
	
	addme("unique species", $species);
	addme("unique ensps", $ensp);

	stepme(10000);
}
close(TREES);
stopme();
stoptime();
state("Inserted ".getme()." ENSP -> Species mappings into '$enspspeciestable'");
showme("unique species", 1);
showme("unique ensps", 1);
clearme("unique species");
clearme("unique ensps");


startme("Inserting Compara FASTA alignments from 'tmp/*.fasta.txt' into table '$fastatable'");
startme2();
starttime();
while (-e "tmp/".sprintf("%06d", getme() + 1).".fasta.txt")
{
	$infile = "tmp/".sprintf("%06d", getme() + 1).".fasta.txt";
	open(IN, "$infile") or die("\nError: Couldn't open '$infile'\n\n");
	fastabreak();
	while (<IN>)
	{
		($ensp, $seq) = getfasta();
		
		$query = Query("SELECT species FROM `$enspspeciestable` WHERE ensp='$ensp'");
		($species) = FetchOne($query);
		
		Query("INSERT INTO `$fastatable` SET ensp='$ensp', species='$species', aln=".(getme() + 1).", seq='$seq'");
		stepme2();
		
		addme("unique species", $species);
		addme("unique ensps", $ensp);
		addme("unique alignments", getme()+1);
	}
	normalbreak();
	close(IN);
	
	stepme(100);
}
stopme();
stoptime();
state("Inserted alignment sequences for ".getme2()." ENSPs from ".getme()." alignments");
showme("unique species", 1);
showme("unique ensps", 1);
showme("unique alignments", 1);
clearme("unique species");
clearme("unique ensps");
clearme("unique alignments");

startme("Inserting Compara CDS FASTA alignments from 'tmp/*.cds.txt' into table '$cdstable'");
startme2();
starttime();
while (-e "tmp/".sprintf("%06d", getme() + 1).".cds.txt")
{
	$infile = "tmp/".sprintf("%06d", getme() + 1).".cds.txt";
	open(IN, "$infile") or die("\nError: Couldn't open '$infile'\n\n");
	fastabreak();
	while (<IN>)
	{
		($ensp, $seq) = getfasta();
		
		$query = Query("SELECT species FROM `$enspspeciestable` WHERE ensp='$ensp'");
		($species) = FetchOne($query);
		
		Query("INSERT INTO `$cdstable` SET ensp='$ensp', species='$species', aln=".(getme() + 1).", seq='$seq'");
		stepme2();
		
		addme("unique species", $species);
		addme("unique ensps", $ensp);
		addme("unique alignments", getme()+1);
	}
	normalbreak();
	close(IN);
	
	stepme(100);
}
stopme();
stoptime();
state("Inserted alignment sequences for ".getme2()." ENSPs from ".getme()." alignments");
showme("unique species", 1);
showme("unique ensps", 1);
showme("unique alignments", 1);
clearme("unique species");
clearme("unique ensps");
clearme("unique alignments");

startme("Inserting Compara NH trees from 'tmp/*.tree.txt' into table '$nhtable'");
starttime();
while (-e "tmp/".sprintf("%06d", getme() + 1).".tree.txt")
{
	$infile = "tmp/".sprintf("%06d", getme() + 1).".tree.txt";
	open(IN, "$infile") or die("\nError: Couldn't open '$infile'\n\n");
	$tree = "";
	while (<IN>)
	{
		$tree .= $_;
	}
	close(IN);
	
	Query("INSERT INTO `$nhtable` SET aln=".(getme() + 1).", tree='$tree'");

	addme("unique alignments", getme()+1);
	
	stepme(100);
}
stopme();
stoptime();
showme("unique alignments", 1);
clearme("unique alignments", 1);

startme("Inserting Compara NHX trees from 'tmp/*.treex.txt' into table '$nhxtable'");
starttime();
while (-e "tmp/".sprintf("%06d", getme() + 1).".treex.txt")
{
	$infile = "tmp/".sprintf("%06d", getme() + 1).".treex.txt";
	open(IN, "$infile") or die("\nError: Couldn't open '$infile'\n\n");
	$treex = "";
	while (<IN>)
	{
		$treex .= $_;
	}
	close(IN);
	
	Query("INSERT INTO `$nhxtable` SET aln=".(getme() + 1).", tree='$treex'");

	addme("unique alignments", getme()+1);
	
	stepme(100);
}
stopme();
stoptime();
showme("unique alignments", 1);
clearme("unique alignments");

done();
