#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$enspspeciestable = 'comparaenspspecies';

our $usage = "$0 [Compara release number]\nExample: $0 108";
($rel) = args(1);

# start
$treefile = "input/Compara.$rel.protein_default.nh.emf.gz";
open(TREES, "zcat $treefile|") or die("\nError: Couldn't zcat '$treefile'\n\n");

Clear($enspspeciestable);

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
state("Inserted ".commify(getme())." ENSP -> Species mappings into '$enspspeciestable'");
showmeall(1);

if (getme() != countme("unique ensps"))
{
	die("Error: ENSPs aren't unique (some are probably in multiple species)");
}
else
{
	state("All ENSPs are unique, and therefore uniquely assigned to one species.");
	run("Query", "~/scripts/query.pl 'SELECT COUNT(*), COUNT(DISTINCT ensp), COUNT(DISTINCT species, ensp) FROM comparaenspspecies' -h");
}

nl();
Optimize($enspspeciestable);

done();
