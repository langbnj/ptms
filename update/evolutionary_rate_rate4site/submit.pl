#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$maxnodes = 8;        # maximum number of nodes to use for 'extra' task
# $maxnodes = 1;        # maximum number of nodes to use for 'extra' task



# our $usage = "$0 [species] [mode: 'fast', 'slow'] [-clear] [-extra]\n\n -clear: remove existing evorates from MySQL table\n-extra: only run the tmp-task-extra.txt file (which is not run at all by default)\n\nExample: $0 human fast";
our $usage = "$0 [species] [type] [-clear] [-extra]\n\n -clear: remove existing evorates from MySQL table\n-extra: only run the tmp-task-extra.txt file\n\nExample: $0 human linsi_tree";
($species, $type) = args(2);
# die("Error: mode must be 'fast' or 'slow'") if (($mode ne 'fast') and ($mode ne 'slow'));

if (switch('clear'))
{
    Clear("evorate_rate4site_$type");
}

# start

chdir("tmp");

if (!switch('extra'))
{
	@tasks = `ls -1 ../tmp-task-*.txt`;
	@tasks = unique(@tasks);
}
else
{
	@tasks = ();
	push(@tasks, "../tmp-task-extra.txt");
}

foreach $task (@tasks)
{
	$task =~ /\.\.\/tmp-task-(\w+)\.txt/ or die("Error: Couldn't parse task file name '$task'");
	
	if (switch('extra'))
	{
		$nodes = $maxnodes;
	}
	else
	{
		$nodes = $1;
	}
	next if ($nodes eq 'extra');

	# Convert CPUs into memory (GB)
	$mem = $nodes * 16;
	# $mem .= "G";	# Not necessary, qsubm.sh will append this


	open(TASK, "$task") or die("Error: Couldn't open task file '$task'");
	
	startme("Running rate4site on task file '$task'");
	starttime();
	$free = freenodes();
	while (<TASK>)
	{
		chomp;
		
		($aln, $ensp) = split(/\|/);
		
		while ($free < $nodes)
		{
			sleep(5);
			$free = freenodes();
		}

        # run("Main", "~/scripts/qsub$nodes.sh ../job.pl $species $type $aln $ensp", 1);
        run("Main", "~/scripts/qsubm.sh $mem ../job.pl $species $type $aln $ensp", 1);



        # $free -= $nodes;
		$free--;

		stepme(1);
	}
	stopme();
	stoptime();
}

chdir("..");

starttime();
waitforjobs();
stoptime();

done();
