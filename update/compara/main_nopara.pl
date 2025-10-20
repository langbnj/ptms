#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [MAFFT mode (linsi/ginsi/einsi)] [-tree]\n\nExample: $0 linsi -tree";
($mode) = args(1);
if (switch('tree')) { $mode .= "_tree"; }
if (switch('1para')) { $mode .= "_1para"; }

$enspspeciestable = 'comparaenspspecies';		# Generic table (okay to use since paralog-free is just a reduced set of the generic Compara)
$fastatable = "comparafasta_$mode";
$nhtable = "comparanh_$mode";


# start

# Create tables if necessary
Query("CREATE TABLE IF NOT EXISTS `$fastatable` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `ensp` varchar(250) NOT NULL,
  `species` varchar(250) DEFAULT NULL,
  `aln` int unsigned DEFAULT NULL,
  `seq` mediumtext,
  PRIMARY KEY (`id`),
  KEY `Ensp` (`ensp`),
  KEY `Species` (`species`),
  KEY `Aln` (`aln`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Paralog-free Compara Alignments (FASTA format)'");

Query("CREATE TABLE IF NOT EXISTS `$nhtable` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `aln` int unsigned DEFAULT NULL,
  `tree` text NOT NULL,
  PRIMARY KEY (`id`),
  KEY `Aln` (`aln`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Paralog-free Compara Trees (Newick format, NH)'");

Clear($fastatable);
Clear($nhtable);


startme("Inserting Compara FASTA alignments from 'tmp/$mode.*.fasta.txt' into table '$fastatable'");
startme2();
starttime();
while (-e "tmp/$mode.".sprintf("%06d", getme() + 1).".fasta.txt")
{
	$infile = "tmp/$mode.".sprintf("%06d", getme() + 1).".fasta.txt";
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


startme("Inserting Compara NH trees from 'tmp/$mode.*.tree.txt' into table '$nhtable'");
starttime();
while (-e "tmp/$mode.".sprintf("%06d", getme() + 1).".tree.txt")
{
	$infile = "tmp/$mode.".sprintf("%06d", getme() + 1).".tree.txt";
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

done();
