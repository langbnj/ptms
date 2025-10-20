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

	run("Clean up output_linsi directory", "mv ./output_linsi ./deleting_output_linsi");
	run("Clean up output_linsi directory", "rm -rf ./deleting_output_linsi &");
	run("Recreate output_linsi directory", "mkdir ./output_linsi");

	run("Clean up output_ginsi directory", "mv ./output_ginsi ./deleting_output_ginsi");
	run("Clean up output_ginsi directory", "rm -rf ./deleting_output_ginsi &");
	run("Recreate output_ginsi directory", "mkdir ./output_ginsi");

	run("Clean up output_einsi directory", "mv ./output_einsi ./deleting_output_einsi");
	run("Clean up output_einsi directory", "rm -rf ./deleting_output_einsi &");
	run("Recreate output_einsi directory", "mkdir ./output_einsi");

	run("Clean up output_linsi_tree directory", "mv ./output_linsi_tree ./deleting_output_linsi_tree");
	run("Clean up output_linsi_tree directory", "rm -rf ./deleting_output_linsi_tree &");
	run("Recreate output_linsi_tree directory", "mkdir ./output_linsi_tree");

	run("Clean up output_ginsi_tree directory", "mv ./output_ginsi_tree ./deleting_output_ginsi_tree");
	run("Clean up output_ginsi_tree directory", "rm -rf ./deleting_output_ginsi_tree &");
	run("Recreate output_ginsi_tree directory", "mkdir ./output_ginsi_tree");

	run("Clean up output_einsi_tree directory", "mv ./output_einsi_tree ./deleting_output_einsi_tree");
	run("Clean up output_einsi_tree directory", "rm -rf ./deleting_output_einsi_tree &");
	run("Recreate output_einsi_tree directory", "mkdir ./output_einsi_tree");
}



$text = "\nAlso remove temporary files in input/ (while keeping input/hsapiens_gene_ensembl__homolog_*__dm.txt)? (y/n)\n";
print $text;

$input = <STDIN>;
chomp($input);
$clean = 'some';
if (lc($input) eq 'y') { $clean = 'all'; }



if ($clean eq 'all')
{

	run("Clean up input directory", "mv ./input ./deleting_input");
	run("Recreate input directory", "mkdir ./input");
	run("Keep hsapiens_gene_ensembl__homolog_*__dm.txt input files", "mv ./deleting_input/hsapiens_gene_ensembl__homolog_*__dm.txt ./input/");
	run("Clean up input directory", "rm -rf ./deleting_input &");
}


done();
