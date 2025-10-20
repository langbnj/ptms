#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "$0 [type]";
($type) = args(1);

$newtype = "norm_$type";		# From evorate_lichtarge_linsi_tree to evorate_norm_lichtarge_linsi_tree

$table = "evorate_$type";
$newtable = "evorate_$newtype";

# start

Query("CREATE TABLE IF NOT EXISTS `$newtable` (
  `id` int(10) unsigned NOT NULL auto_increment,
  `ensp` char(15) default NULL,
  `site` mediumint(8) default NULL,
  `rate` float default NULL,
  PRIMARY KEY  (`id`),
  KEY `EnspSite` (`ensp`,`site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Evolutionary rates $newtype'
");

Clear($newtable);


$mainquery = Query("SELECT ensp, AVG(rate), STD(rate) FROM `$table` GROUP BY ensp ORDER BY ensp");
startme("Normalising evorate type '$type' to '$newtype' (per protein)", 0, Numrows($mainquery));
startme2();
starttime();
while (($ensp, $avg, $std) = Fetch($mainquery))
{
	# Workaround for histone H4 (100% conserved, resulting in ~30 evorates of exactly 1, with STDEV 0)
	if (($avg == 1) and ($std == 0))
	{
		$avg = 0;
		$std = 1;
		# >> Results in no change
	}
	if ($std == 0)
	{
		die("Error: Standard deviation is supposedly '$std' == 0 for '$ensp' (type '$type') (average '$avg')");
	}
	
	$query = Query("SELECT site, rate FROM `$table` WHERE ensp='$ensp'");
	while (($site, $rate) = Fetch($query))
	{
		# Normalize rate
		$newrate = ($rate - $avg) / $std;
		
		Query("INSERT INTO `$newtable` SET ensp='$ensp', site=$site, rate=$newrate");
		stepme2();
	}

	stepme(100);
}
stopme();
stoptime();

Optimize($newtable);

state("Inserted ".commify(getme2())." normalized rates for ".commify(getme())." proteins into '$newtable'");

done();
