#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
our $usage = "$0 [mode: e.g. linsi_tree, ginsi] [window size]\n\nExample: $0 linsi_tree 0";
($mode, $window) = args(2);

# $table = 'evorate_capra'.$window.'_'.$mode;
# Query("TRUNCATE TABLE $table");
# state("Cleared table '$table'");

# start

startme("Running Capra (mode '$mode' window '$window') on alignments from Compara");
starttime();
chdir("tmp");
while (-e $mode.".".sprintf("%06d", getme() + 1).".fasta.txt")
{
	# Cluster
	$free = freenodes();
	$mine = mynodes();

	if (($free <= 0) or ($mine >= 10000))
	{
		sleep(10);
	}

	while (($free > 0) and ($mine < 10000) and (-e $mode.".".sprintf("%06d", getme() + 1).".fasta.txt"))
	{
		stepme(100);
		run("Main", "~/scripts/qsub.sh ../main.pl ".getme()." $mode $window", 1);

		$free--;
		$mine++;
	}

	# # Just run locally
	# stepme(100);
	# run("Main", "../main.pl ".getme()." $mode $window", 1);
}
chdir("..");
stopme();
stoptime();

done();
