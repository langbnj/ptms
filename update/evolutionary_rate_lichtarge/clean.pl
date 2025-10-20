#!/usr/bin/env perl -w

require('functions.inc.pl');
#require('mysql.inc.pl');


# Safety prompt

$text = "\nAlso remove output files? (y/n)\n";
print $text;

$input = <STDIN>;
chomp($input);
$clean = 'some';
if (lc($input) eq 'y') { $clean = 'all'; }



# Start

run("Clean up temporary files", "rm -f tmp-*");

run("Clean up temporary directory", "mv ./tmp ./deleting_tmp");
run("Clean up temporary directory", "rm -rf ./deleting_tmp &");
run("Recreate temporary directory", "mkdir ./tmp");


if ($clean eq 'all')
{
	run("Clean up output files", "rm -f output-*");

	run("Clean up output directory", "mv ./output ./deleting_output");
	run("Clean up output directory", "rm -rf ./deleting_output &");
	run("Recreate output directory", "mkdir ./output");
}



# $text = "\nAlso remove input files? (y/n)\n";
# print $text;
# 
# $input = <STDIN>;
# chomp($input);
# $clean = 'some';
# if (lc($input) eq 'y') { $clean = 'all'; }
# 
# 
# 
# if ($clean eq 'all')
# {
# 	# run("Clean up input files", "rm -f input-*");
# 
# 	run("Clean up input directory", "mv ./input ./deleting_input");
# 	run("Clean up input directory", "rm -rf ./deleting_input &");
# 	run("Recreate input directory", "mkdir ./input");
# }


done();
