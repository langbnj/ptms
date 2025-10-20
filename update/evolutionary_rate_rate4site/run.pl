#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize


# run

# @modes = ('linsi', 'ginsi', 'einsi', 'linsi_tree', 'ginsi_tree', 'einsi_tree', 'para');
# @modes = ('linsi', 'ginsi', 'einsi', 'linsi_tree', 'ginsi_tree', 'einsi_tree');
@modes = ('einsi_tree_1para');

foreach $mode (@modes)
{
	run("Main", "main.pl human $mode");
	waitforjobs();
	run("Backup", "~/scripts/backup.pl $mode");
}

# # Parse from "backups"
# run("Swap in linsi_tree", "mv backup_2023_03_17_linsi_tree/output backup_2023_03_17_linsi_tree/tmp .");
# run("Parse linsi_tree", "main.pl human linsi_tree");
# 
# run("Swap out linsi_tree", "mv output tmp backup_2023_03_17_linsi_tree");
# run("Swap in linsi", "mv backup_2023_03_17_linsi/output backup_2023_03_17_linsi/tmp .");
# run("Parse linsi", "main.pl human linsi");
# 
# run("Swap out linsi", "mv output tmp backup_2023_03_17_linsi");
# run("Swap in ginsi_tree", "mv backup_2023_03_17_ginsi_tree/output backup_2023_03_17_ginsi_tree/tmp .");
# run("Parse ginsi_tree", "main.pl human ginsi_tree");
# 
# run("Swap out ginsi_tree", "mv output tmp backup_2023_03_17_ginsi_tree");
# run("Swap in ginsi", "mv backup_2023_03_18_ginsi/output backup_2023_03_18_ginsi/tmp .");
# run("Parse ginsi", "main.pl human ginsi");
# 
# run("Swap out ginsi", "mv output tmp backup_2023_03_18_ginsi");
# run("Swap in einsi_tree", "mv backup_2023_03_18_einsi_tree/output backup_2023_03_18_einsi_tree/tmp .");
# run("Parse einsi_tree", "main.pl human einsi_tree");
# 
# run("Swap out einsi_tree", "mv output tmp backup_2023_03_18_einsi_tree");
# run("Swap in einsi", "mv backup_2023_03_18_einsi/output backup_2023_03_18_einsi/tmp .");
# run("Parse einsi", "main.pl human einsi");
# 
# run("Swap out einsi", "mv output tmp backup_2023_03_18_einsi");

# $mode = 'linsi_tree'; run("Main", "main.pl human $mode");   waitforjobs(); run("Backup", "~/scripts/backup.pl $mode");
# $mode = 'linsi'; run("Main", "main.pl human $mode");        waitforjobs(); run("Backup", "~/scripts/backup.pl $mode");
# $mode = 'ginsi_tree'; run("Main", "main.pl human $mode");   waitforjobs(); run("Backup", "~/scripts/backup.pl $mode");
# $mode = 'ginsi'; run("Main", "main.pl human $mode");        waitforjobs(); run("Backup", "~/scripts/backup.pl $mode");
# $mode = 'einsi_tree'; run("Main", "main.pl human $mode");   waitforjobs(); run("Backup", "~/scripts/backup.pl $mode");
# $mode = 'einsi'; run("Main", "main.pl human $mode");        waitforjobs(); run("Backup", "~/scripts/backup.pl $mode");
# Note: Need to clean between runs (hence the 'backup')
# run("Main", "main.pl human para");
# run("Main", "main.pl yeast linsi_tree");
# run("Main", "main.pl yeast para");

# Check counts
state("Checking ENSP counts (in input and output):");
foreach $mode (@modes)
{
    $comparafasta = "comparafasta_$mode";
    $comparafasta = "comparafasta" if ($mode eq 'para');
    

    $query = Query("SELECT COUNT(DISTINCT ensp), SUM(LENGTH(REPLACE(seq, '-', ''))) AS residues FROM $comparafasta WHERE species='homo_sapiens'");
    ($ensps_comparafasta, $residues_comparafasta) = FetchOne($query);
    state(" >> $comparafasta >> ".commify($ensps_comparafasta)." ENSPs >> ".commify($residues_comparafasta)." residues", 1);

    $query = Query("SELECT COUNT(DISTINCT ensp), COUNT(DISTINCT ensp, site), COUNT(*), MIN(rate), AVG(rate), MAX(rate) FROM evorate_rate4site_$mode");
    ($ensps_evorate, $residues_evorate, $rows_evorate, $minrate, $avgrate, $maxrate) = FetchOne($query);
    state("   >> evorate_rate4site_$mode >> ".commify($ensps_evorate)." ENSPs >> ".commify($residues_evorate)." residues >> MIN $minrate AVG $avgrate MAX $maxrate", 1);

    warn("Warning: 'evorate_rate4site_$mode' should have ".commify($residues_comparafasta)." residues according to $comparafasta, but it has ".commify($residues_evorate)." residues") if ($residues_evorate != $residues_comparafasta);
    warn("Warning: 'evorate_rate4site_$mode' has ".commify($residues_evorate)." residues, but ".commify($rows_evorate)." total rows") if ($residues_evorate != $rows_evorate);
    
    if ($mode eq 'para')
    {
        $query = Query("SELECT COUNT(DISTINCT ensp) FROM $comparafasta WHERE seq REGEXP 'X' AND species='homo_sapiens';");
        ($x_containing_ensps) = FetchOne($query);

        warn("Warning: 'evorate_rate4site_$mode' is incomplete (".commify($ensps_evorate)." instead of the expected ".commify($ensps_comparafasta)." from $comparafasta) (".($ensps_comparafasta - $ensps_evorate)." missing, which is probably due to X-containing human ENSP sequences (there are $x_containing_ensps of them) for comparafasta (para))") if ($ensps_comparafasta ne $ensps_evorate);
    }
    else
    {
        warn("Warning: 'evorate_rate4site_$mode' is incomplete (".commify($ensps_evorate)." instead of the expected ".commify($ensps_comparafasta)." from $comparafasta) (".($ensps_comparafasta - $ensps_evorate)." missing)") if ($ensps_comparafasta ne $ensps_evorate);
    }
}

state("Optimizing tables:");
foreach $mode (@modes)
{
    Optimize("evorate_rate4site_$mode");
}

# Longest jobs:
# in tmp: cat log-out* | g done | g elaps | sort -h
# >> ~ 2 hours 40 minutes, now (used to be days due to low memory)
# >> ~3500 jobs max simultaneously (all single-core)

done();
