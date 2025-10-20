#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

our $usage = "[evorate type: capra/lichtarge/rate4site] [-para]\n\n -para: Run on 'evorate_..._para' tables (by default, run on non-para tables)\n\nExample: $0 capra";
($type) = args(1);
# args(0);

# $infile = "input.txt";
# $outfile = "output.txt";

# open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
# open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");


# start

state("Getting protein intersection that is present for all '$type' evorate types (i.e. none crashed):");
starttime();

# Get core list
$mainquery = Query("SHOW TABLES");
@corelist = ();
while (($table) = Fetch($mainquery))
{
	next if ($table !~ /^evorate_($type)/);              # evorate tables only (e.g. evorate_capra_...)
	if (switch('para'))
	{
		# -para: only _para
		next if ($table !~ /_para$/); # skip "para" tables (based on comparafasta, not the paralog-free e.g. comparafasta_linsi_tree)
	}
	else
	{
		# default: skip _para
		next if ($table =~ /_para$/); # skip "para" tables (based on comparafasta, not the paralog-free e.g. comparafasta_linsi_tree)
	}
	
	print " >> $table";
	$query = Query("SELECT ensp FROM `$table` GROUP BY ensp ORDER BY ensp");
	@thislist = ();
	while (($ensp) = Fetch($query))
	{
		push(@thislist, $ensp);
		# stepme(1000, 1);
	}
	# stopme(1);
	@thislist = unique(@thislist);
	print " >> ".scalar(@thislist)." proteins";
	
	if (scalar(@corelist) == 0)
	{
		@corelist = @thislist;
	}
	
	@corelist = intersection(\@corelist, \@thislist);

	print " >> core list ".scalar(@corelist)."\n";
}


state("Final core list: ".commify(scalar(@corelist))." proteins");



# Get sporadic list
$mainquery = Query("SHOW TABLES");
@deletelist = ();
while (($table) = Fetch($mainquery))
{
	next if ($table !~ /^evorate_($type)/);              # evorate tables only (e.g. evorate_capra_...)
	if (switch('para'))
	{
		# -para: only _para
		next if ($table !~ /_para$/); # skip "para" tables (based on comparafasta, not the paralog-free e.g. comparafasta_linsi_tree)
	}
	else
	{
		# default: skip _para
		next if ($table =~ /_para$/); # skip "para" tables (based on comparafasta, not the paralog-free e.g. comparafasta_linsi_tree)
	}
	
	print " >> $table";
	$query = Query("SELECT ensp FROM `$table` GROUP BY ensp ORDER BY ensp");
	@thislist = ();
	while (($ensp) = Fetch($query))
	{
		push(@thislist, $ensp);
		# stepme(1000, 1);
	}
	# stopme(1);
	@thislist = unique(@thislist);
	print " >> ".scalar(@thislist)." proteins";
	
	# foreach $missing (symdiff(\@corelist, \@thislist))
	foreach $missing (ronly(\@corelist, \@thislist))
	{
		if (!contains($missing, @deletelist))
		{
			print "\n   >> newly missing: '$missing'";

			push(@deletelist, $missing);

			addme("proteins to delete", $missing);
		}
	}
	
	print " >> sporadic list ".scalar(@deletelist)."\n";
}
stoptime();



nl();
state("Final core list: ".commify(scalar(@corelist))." proteins", 1);
state("Final sporadic list (to be deleted): ".commify(scalar(@deletelist))." proteins", 1);
state(join("\n", @deletelist), 1);
nl();



# Now cycle through all tables and delete sporadic proteins (@deletelist)
startme("Cycling through tables and deleting sporadic proteins");
starttime();

$mainquery = Query("SHOW TABLES");
$affected = 0;
@deletelist = unique(@deletelist);
while (($table) = Fetch($mainquery))
{
	next if ($table !~ /^evorate_($type)/);              # evorate tables only (e.g. evorate_capra_...)
	if (switch('para'))
	{
		# -para: only _para
		next if ($table !~ /_para$/); # skip "para" tables (based on comparafasta, not the paralog-free e.g. comparafasta_linsi_tree)
	}
	else
	{
		# default: skip _para
		next if ($table =~ /_para$/); # skip "para" tables (based on comparafasta, not the paralog-free e.g. comparafasta_linsi_tree)
	}

	foreach $missing (@deletelist)
	{
		if (!switch('debug'))
		{
			$query = Query("DELETE FROM `$table` WHERE ensp='$missing'");
		}
		else
		{
			$query = Query("SELECT id FROM `$table` WHERE ensp='$missing'");
		}
		
		if (Numrows($query) > 0)
		{
			addme("successfully deleted rows for table|protein", "$table|$missing");
			addme("successfully deleted rows for protein", $missing);
			addme("successfully deleted rows from table", $table);
		}
		
		$affected += Numrows($query);
	}
	
	Optimize($table, 1);
	# addme("optimized table", $table);

	stepme(1);
}
stopme();
stoptime();

# showmeall(1);
showmesome(50);

state("$affected rows affected");

state("Optimized all tables");

done();
