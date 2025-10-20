#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# run

# @modes = ('linsi', 'ginsi', 'einsi', 'linsi_tree', 'ginsi_tree', 'einsi_tree', 'para');
# @modes = ('linsi', 'ginsi', 'einsi', 'linsi_tree', 'ginsi_tree', 'einsi_tree');
@modes = ('einsi_tree_1para');
@windows = (0, 1, 3, 5);

starttime();

nl();
foreach $mode (@modes)
{
    foreach $window (@windows)
    {
        $type = 'capra'.$window.'_'.$mode;
        $table = 'evorate_'.$type;
        Query("TRUNCATE TABLE $table");
        state("Cleared table '$table'", 1);
    }
}
nl();

foreach $mode (@modes)
{
    run("Split input files into single alignments", "~/scripts/qsub.sh split.pl $mode");
}
waitforjobs();

foreach $mode (@modes)
{
    foreach $window (@windows)
    {
        run("Check (window size $window)", "~/scripts/qsub.sh check.pl $mode $window");
    }
}
waitforjobs();

foreach $mode (@modes)
{
    foreach $window (@windows)
    {
        run("Run jobs (window size $window)", "~/scripts/qsub.sh submit.pl $mode $window");
    }
}
waitforjobs();

foreach $mode (@modes)
{
    foreach $window (@windows)
    {
        run("Check (window size $window)", "~/scripts/qsub.sh check.pl $mode $window");
        # run("Check (window size $window)", "check.pl $mode $window");
    }
}
waitforjobs();



# No longer the case! Inverted below.
# warn("NOTE: Capra is a conservation score, unlike an evolutionary rate (variation score). High Capra means good conservation.
# Test case to see this:
# SELECT * FROM comparafasta_linsi_tree WHERE species='homo_sapiens' ORDER BY LENGTH(seq);
# SELECT * FROM comparafasta_linsi_tree WHERE aln=18120;
# SELECT * FROM evorate_capra0_linsi_tree WHERE ensp='ENSP00000502602';
# ");

# Invert Capra scores
state("Invert Capra scores (from a conservation score, high = good conservation, to a variation score, high = fast evolution, like the others)");
foreach $mode (@modes)
{
    foreach $window (@windows)
    {
        print " >> evorate_capra$window\_$mode";

        $query = Query("SELECT id FROM evorate_capra$window\_$mode WHERE rate<0");
        if (Numrows($query) > 0)
        {
            die("Error: Found rates < 0 in table 'evorate_capra$window\_$mode'");
        }
        $query = Query("SELECT id FROM evorate_capra$window\_$mode WHERE rate>1");
        if (Numrows($query) > 0)
        {
            die("Error: Found rates > 1 in table 'evorate_capra$window\_$mode'");
        }

        # Invert score
        $query = Query("UPDATE evorate_capra$window\_$mode SET rate=(1-rate)");
        print " >> ".commify(Numrows($query))." rows affected\n";
    }
}

# Check counts
state("Checking ENSP counts (in input and output):");
foreach $mode (@modes)
{
    $comparafasta = "comparafasta_$mode";
    $comparafasta = "comparafasta" if ($mode eq 'para');
    
    $query = Query("SELECT COUNT(DISTINCT ensp), SUM(LENGTH(REPLACE(seq, '-', ''))) AS residues FROM $comparafasta WHERE species='homo_sapiens'");
    ($ensps_comparafasta, $residues_comparafasta) = FetchOne($query);
    state(" >> $comparafasta >> ".commify($ensps_comparafasta)." ENSPs >> ".commify($residues_comparafasta)." residues", 1);
    foreach $window (@windows)
    {
        $query = Query("SELECT COUNT(DISTINCT ensp), COUNT(DISTINCT ensp, site), COUNT(*), MIN(rate), AVG(rate), MAX(rate) FROM evorate_capra$window\_$mode");
        ($ensps_evorate, $residues_evorate, $rows_evorate, $minrate, $avgrate, $maxrate) = FetchOne($query);
        state("   >> evorate_capra$window\_$mode >> ".commify($ensps_evorate)." ENSPs >> ".commify($residues_evorate)." residues >> MIN $minrate AVG $avgrate MAX $maxrate", 1);

        warn("Warning: 'evorate_capra$window\_$mode' should have ".commify($residues_comparafasta)." residues according to $comparafasta, but it has ".commify($residues_evorate)." residues") if ($residues_evorate != $residues_comparafasta);
        warn("Warning: 'evorate_capra$window\_$mode' has ".commify($residues_evorate)." residues, but ".commify($rows_evorate)." total rows") if ($residues_evorate != $rows_evorate);
        
        if ($mode eq 'para')
        {
            $query = Query("SELECT COUNT(DISTINCT ensp) FROM $comparafasta WHERE seq REGEXP 'X' AND species='homo_sapiens';");
            ($x_containing_ensps) = FetchOne($query);

            warn("Warning: 'evorate_capra$window\_$mode' is incomplete (".commify($ensps_evorate)." instead of the expected ".commify($ensps_comparafasta)." from $comparafasta) (".($ensps_comparafasta - $ensps_evorate)." missing, which is probably due to X-containing human ENSP sequences (there are $x_containing_ensps of them) for comparafasta (para))") if ($ensps_comparafasta ne $ensps_evorate);
        }
        else
        {
            warn("Warning: 'evorate_capra$window\_$mode' is incomplete (".commify($ensps_evorate)." instead of the expected ".commify($ensps_comparafasta)." from $comparafasta) (".($ensps_comparafasta - $ensps_evorate)." missing)") if ($ensps_comparafasta ne $ensps_evorate);
        }
    }
}

state("Optimizing tables:");
foreach $mode (@modes)
{
    foreach $window (@windows)
    {
        Optimize("evorate_capra$window\_$mode");
    }
}

stoptime();
done();
