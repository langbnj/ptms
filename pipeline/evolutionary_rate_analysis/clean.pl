#!/home/blang1/bin/perl -w

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

if ($clean eq 'all')
{
	run("Clean up temporary files", "rm -f tmp-*");

	run("Clean up temporary directory", "rm -rf ./tmp");
	run("Recreate temporary directory", "mkdir ./tmp");


	run("Clean up output files", "rm -f output-*");

	run("Clean up output directory", "rm -rf ./output");
	run("Recreate output directory", "mkdir ./output");
}






done();