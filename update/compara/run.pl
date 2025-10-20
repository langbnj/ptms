#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
$rel = 108;
our $usage = '$0 [-tests]\n\n -tests: Run tests only\n\nExample: $0 -tests';
args(0);

# run

# Download
# run("Download", "download.pl");

state("Steps:
1. Run ~/update/compara/run.pl
2. Run ~/update/compara/run.pl -tests         # All OK (now that uniens only contains mappings between UniProt and Ensembl where the sequences match perfectly)
3. Run ~/update/comparanopara/run.pl
4. Run ~/update/compara/run.pl -nopara");

if (!switch('nopara'))
{
    if (!switch('tests'))
    {
        # Update Compara
        run("Fill 'comparaenspspecies' table", "comparaenspspecies.pl $rel");
        run("Split", "split.pl $rel");
        run("Parse into compara tables", "main.pl $rel");
        run("Fill 'comparaspecies' table", "comparaspecies.pl");
    }
    else
    {
        # Testing
        run("NH tree vs. NHX tree agreement test", "tree_treex_test.pl");
        run("'uniens' test for UniProt/ENSEMBL ENSP ID agreement", "uniens_test.pl $rel");
        run("'uniseq' test for UniProt/ENSEMBL sequence agreement", "uniseq_ensembl_test.pl");
        run("'uniseq' test for UniProt/Compara sequence agreement", "uniseq_compara_test.pl");
        run("'unimod' test for UniProt PTM / UniProt site amino acid agreement", "unimod_uniseq_test.pl");
        run("'unimod' test for UniProt PTM / ENSEMBL site amino acid agreement", "unimod_ensembl_test.pl");
        run("'unimod' test for UniProt PTM / Compara site amino acid agreement", "unimod_compara_test.pl");
    }
}
# Update paralog-free Comparas (requires the other, at least for the comparaenspspecies table)
# Steps:
# 1. Run ~/update/compara/run.pl
# 2. Run ~/update/compara/run.pl -tests         # All OK (now that uniens only contains mappings between UniProt and Ensembl where the sequences match perfectly)
# 3. Run ~/update/comparanopara/run.pl
# 4. Run ~/update/compara/run.pl -nopara
elsif (switch('nopara'))
{
    # Running these locally since they're fast
    # run("Split paralog-free version", "split_nopara.pl einsi");
    # run("Split paralog-free version", "split_nopara.pl ginsi");
    # run("Split paralog-free version", "split_nopara.pl linsi");
    # run("Split paralog-free version", "split_nopara.pl einsi -tree");
    # run("Split paralog-free version", "split_nopara.pl ginsi -tree");
    # run("Split paralog-free version", "split_nopara.pl linsi -tree");
    run("Split paralog-free version", "split_nopara.pl einsi -tree -1para");

    # Running these on the cluster
    # run("Parse into paralog-free Compara tables", "~/scripts/qsub.sh main_nopara.pl einsi");
    # run("Parse into paralog-free Compara tables", "~/scripts/qsub.sh main_nopara.pl ginsi");
    # run("Parse into paralog-free Compara tables", "~/scripts/qsub.sh main_nopara.pl linsi");
    # run("Parse into paralog-free Compara tables", "~/scripts/qsub.sh main_nopara.pl einsi -tree");
    # run("Parse into paralog-free Compara tables", "~/scripts/qsub.sh main_nopara.pl ginsi -tree");
    # run("Parse into paralog-free Compara tables", "~/scripts/qsub.sh main_nopara.pl linsi -tree");
    run("Parse into paralog-free Compara tables", "~/scripts/qsub.sh main_nopara.pl einsi -tree -1para");
    waitforjobs();
}

done();
