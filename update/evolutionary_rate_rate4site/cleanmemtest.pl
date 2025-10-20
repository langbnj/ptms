#!/usr/bin/env perl -w

require('functions.inc.pl');
#require('mysql.inc.pl');


# Start

run("Clean up temporary directory", "rm -rf ./memtest");
run("Recreate temporary directory", "mkdir ./memtest");

run("Clean up output files", "rm -f output-memtest-*.txt");






done();


