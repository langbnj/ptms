#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$table = 'tax';

# our $usage = "";
# ($var) = args(1);

$infile = "input/names.dmp";
#$outfile = "output.txt";

open(IN, $infile) or die("\nError: Couldn't open '$infile'\n\n");
#open(OUT, ">$outfile") or die("\nError: Couldn't open '$outfile'\n\n");

Clear($table);


# start

startme("Reading '$infile' and inserting into table '$table'");
starttime();

while (<IN>)
{
	chomp;
	
	@a = split(/\t/);
	
	$tax = $a[0];
	$name = $a[2];
	$type = $a[6];

	# next if ($type ne 'scientific name');
	
	Query("INSERT INTO `$table` SET tax='$tax', name='".esc($name)."', type='$type'");

	stepme(10000);
}
stopme();
stoptime();

Optimize($table);

done();
