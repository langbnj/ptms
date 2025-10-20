#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
our $usage = "$0 [species] [type]\n\nExample: $0 human linsi_tree";
($species, $type) = args(2);

# run
starttime();
run("Split input files into single alignments", "split.pl $type");
run("Check input alignments for sequence length & number of sequences", "check.pl $species $type");
run("Submit jobs", "submit.pl $species $type -all");
stoptime();

starttime();
state("Waiting for jobs...");
waitforjobs();
stoptime();

run("Parse output", "parse.pl output $species $type -clear");

run("Check input alignments for sequence length & number of sequences", "check.pl $species $type");

done();
