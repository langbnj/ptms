#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
our $usage = "$0 [species] [type]\n\nExample: $0 human linsi_tree";
($species, $type) = args(2);

# run
starttime();
# run("Split input files into single alignments", "split.pl $type");
# run("Check & store alignment length and number of sequences", "check.pl $species $type");
# # run("Make task files containing alignments that need to be worked on (because they have $species PTMs), split according to their memory use estimates", "maketasks.pl $species $type -noread -ptmonly");
# run("Make task files containing alignments that need to be worked on (because they have $species proteins), split according to their memory use estimates", "maketasks.pl $species $type -noread");

run("Submit tasks/jobs (jobs update MySQL)", "submit.pl $species $type -clear");
# run("Submit tasks/jobs (jobs update MySQL)", "submit.pl $species $type");
stoptime();
starttime();
waitforjobs();
stoptime();
if (-e "tmp-task-extra.txt")
{
    run("Submit extra tasks/jobs (jobs update MySQL)", "submit.pl $species $type -extra");
    starttime();
    waitforjobs();
}

# print " >> Renaming rate4site_fast >> rate4site\n"; $query = Query("UPDATE evorate SET type='rate4site_$type' WHERE type='rate4site_fast'"); print "   >> ".Numrows($query)." rows affected\n";

stoptime();

starttime();
run("Check & store alignment length and number of sequences", "check.pl $species $type");
stoptime();

done();
