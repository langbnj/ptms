#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize
our $usage = '';
args(0);


# Download

# run("Download sequences for alignments", "download.pl");

# Split

# # # run("Split into one-to-one ortholog clusters", "split.pl");
# # run("Split into one-to-one ortholog clusters (high-confidence only)", "split.pl -hc");
# # run("Split into ortholog clusters (high-confidence only), while keeping the best-matching outparalog by query_percent_id for one2many and many2many cases", "split_1para.pl -hc");
# run("Split into ortholog clusters, while keeping the best-matching outparalog for one2many and many2many cases (criteria: ORDER BY hc DESC, (gene_order_score + whole_genome_align_score) / 2 DESC, gene_order_score DESC, whole_genome_align_score DESC, query_percent_id DESC, target_percent_id DESC LIMIT 1)", "split_1para.pl");


# Run

if (!switch('analysisonly'))
{

    # run("Run MAFFT alignments", "~/scripts/qsub.sh mafft.pl human linsi");
    # run("Run MAFFT alignments", "~/scripts/qsub.sh mafft.pl human ginsi");
    # run("Run MAFFT alignments", "~/scripts/qsub.sh mafft.pl human einsi");
    # 
    # run("Run MAFFT alignments, with species tree", "~/scripts/qsub.sh mafft_tree.pl human linsi");
    # run("Run MAFFT alignments, with species tree", "~/scripts/qsub.sh mafft_tree.pl human ginsi");
    # run("Run MAFFT alignments, with species tree", "~/scripts/qsub.sh mafft_tree.pl human einsi");

    run("Run MAFFT alignments, with species tree, while retaining best paralog in one2many and many2many cases", "~/scripts/qsub.sh mafft_tree.pl human einsi -1para");

    waitforjobs();
    
    # run("Back-translate to cDNA alignments from MAFFT alignments", "~/scripts/qsub.sh backtranslate.pl human linsi");
    # run("Back-translate to cDNA alignments from MAFFT alignments", "~/scripts/qsub.sh backtranslate.pl human ginsi");
    # run("Back-translate to cDNA alignments from MAFFT alignments", "~/scripts/qsub.sh backtranslate.pl human einsi");
    # 
    # run("Back-translate to cDNA alignments from MAFFT alignments, with species tree", "~/scripts/qsub.sh backtranslate.pl human linsi -tree");
    # run("Back-translate to cDNA alignments from MAFFT alignments, with species tree", "~/scripts/qsub.sh backtranslate.pl human ginsi -tree");
    # run("Back-translate to cDNA alignments from MAFFT alignments, with species tree", "~/scripts/qsub.sh backtranslate.pl human einsi -tree");

    run("Back-translate to cDNA alignments from MAFFT alignments, with species tree, while retaining best paralog in one2many and many2many cases", "~/scripts/qsub.sh backtranslate.pl human einsi -tree -1para");

    waitforjobs();
    
    # run("Add NCBI taxon IDs to treebest backtrans CDS", "~/scripts/qsub.sh cds_addtaxonids.pl human linsi");
    # run("Add NCBI taxon IDs to treebest backtrans CDS", "~/scripts/qsub.sh cds_addtaxonids.pl human ginsi");
    # run("Add NCBI taxon IDs to treebest backtrans CDS", "~/scripts/qsub.sh cds_addtaxonids.pl human einsi");
    # 
    # run("Add NCBI taxon IDs to treebest backtrans CDS, for alignments generated with species tree", "~/scripts/qsub.sh cds_addtaxonids.pl human linsi -tree");
    # run("Add NCBI taxon IDs to treebest backtrans CDS, for alignments generated with species tree", "~/scripts/qsub.sh cds_addtaxonids.pl human ginsi -tree");
    # run("Add NCBI taxon IDs to treebest backtrans CDS, for alignments generated with species tree", "~/scripts/qsub.sh cds_addtaxonids.pl human einsi -tree");

    run("Add NCBI taxon IDs to treebest backtrans CDS, for alignments generated with species tree, while retaining best paralog in one2many and many2many cases", "~/scripts/qsub.sh cds_addtaxonids.pl human einsi -tree -1para");

    waitforjobs();
    
    # run("Build tree using TreeBeST", "~/scripts/qsub.sh treebest.pl human linsi");
    # run("Build tree using TreeBeST", "~/scripts/qsub.sh treebest.pl human ginsi");
    # run("Build tree using TreeBeST", "~/scripts/qsub.sh treebest.pl human einsi");
    # 
    # run("Build tree using TreeBeST, for alignments generated with species tree", "~/scripts/qsub.sh treebest.pl human linsi -tree");
    # run("Build tree using TreeBeST, for alignments generated with species tree", "~/scripts/qsub.sh treebest.pl human ginsi -tree");
    # run("Build tree using TreeBeST, for alignments generated with species tree", "~/scripts/qsub.sh treebest.pl human einsi -tree");

    run("Build tree using TreeBeST, for alignments generated with species tree, while retaining best paralog in one2many and many2many cases", "~/scripts/qsub.sh treebest.pl human einsi -tree -1para");

    waitforjobs();
    
    # run("Parse output", "~/scripts/qsub.sh parse.pl human linsi");
    # run("Parse output", "~/scripts/qsub.sh parse.pl human ginsi");
    # run("Parse output", "~/scripts/qsub.sh parse.pl human einsi");
    # 
    # run("Parse output, for alignments generated with species tree", "~/scripts/qsub.sh parse.pl human linsi -tree");
    # run("Parse output, for alignments generated with species tree", "~/scripts/qsub.sh parse.pl human ginsi -tree");
    # run("Parse output, for alignments generated with species tree", "~/scripts/qsub.sh parse.pl human einsi -tree");

    run("Parse output, for alignments generated with species tree, while retaining best paralog in one2many and many2many cases", "~/scripts/qsub.sh parse.pl human einsi -tree -1para");

    waitforjobs();
    


    # Tests

#     run("Test split output", "~/scripts/qsub.sh test_split.pl linsi");
#     run("Test split output", "~/scripts/qsub.sh test_split.pl ginsi");
#     run("Test split output", "~/scripts/qsub.sh test_split.pl einsi");
# 
#     run("Test split output, for alignments generated with species tree", "~/scripts/qsub.sh test_split.pl linsi -tree");
#     run("Test split output, for alignments generated with species tree", "~/scripts/qsub.sh test_split.pl ginsi -tree");
#     run("Test split output, for alignments generated with species tree", "~/scripts/qsub.sh test_split.pl einsi -tree");

    run("Test split output, for alignments generated with species tree, while retaining best paralog in one2many and many2many cases", "~/scripts/qsub.sh test_split.pl einsi -tree -1para");

    # run("Test CDS output", "~/scripts/qsub.sh test_cds.pl human");
    # No real need to wait for these, they're just for verification.

    waitforjobs();
}

# This would immediately run evorate calculations (and their analyses):
if (switch('analysis'))
{

    runanalysis("lichtarge", "linsi");
    runanalysis("lichtarge", "ginsi");
    runanalysis("lichtarge", "einsi");
    runanalysis("lichtarge", "linsi_tree");
    runanalysis("lichtarge", "ginsi_tree");
    runanalysis("lichtarge", "einsi_tree");

    runanalysis("capra", "linsi");
    runanalysis("capra", "ginsi");
    runanalysis("capra", "einsi");
    runanalysis("capra", "linsi_tree");
    runanalysis("capra", "ginsi_tree");
    runanalysis("capra", "einsi_tree");

    runanalysis("rate4site", "linsi");
    runanalysis("rate4site", "ginsi");
    runanalysis("rate4site", "einsi");
    runanalysis("rate4site", "linsi_tree");
    runanalysis("rate4site", "ginsi_tree");
    runanalysis("rate4site", "einsi_tree");
	
	runanalysis("lichtarge", "para");
	runanalysis("capra", "para");
	runanalysis("rate4site", "para");

}


done();


sub runanalysis
{
	my ($evorate, $mode) = @_;
	
	# Calculate evorates for this MAFFT mode / tree combination

	cd("~/update/evolutionary_rate_$evorate/");
	run("ls", "ls");
	run("clean", "~/clean.pl");
	run("clean", "echo y \| clean.pl");
	run("$evorate", "~/scripts/qsub.sh run.pl $mode");
	
	waitforalljobs();

	# # Run evorate analysis for this MAFFT mode / tree combination
    cd("~/pipeline/evolutionary_rate_analysis/");
    run("ls", "ls");
    run("clean", "~/clean.pl");
    run("clean", "echo y \| clean.pl");
    run("evorate analysis run", "~/scripts/qsub.sh run.pl");
    
    # Run evorate quantiles for this MAFFT mode / tree combination
    cd("~/pipeline/evolutionary_rate_quantile_logos/");
    run("ls", "ls");
    run("clean", "~/clean.pl");
    run("clean", "echo y \| clean.pl");
    run("evorate quantiles run", "~/scripts/qsub.sh run.pl");
    
    waitforalljobs();

	# Store results from evorate for this MAFFT mode / tree combination
	cd("~/update/evolutionary_rate_$evorate/");
	run("Make directory for this output", "mkdir output_$mode");
	run("Move", "mv log-* output_$mode");
	run("Move", "mv tmp-* output_$mode");
	run("Move", "mv output-* output_$mode");
	run("Move", "mv tmp output_$mode");
	run("Move", "mv output output_$mode");
	run("mkdir", "mkdir tmp");
	run("mkdir", "mkdir output");

    # Store results from evorate analysis for this MAFFT mode / tree combination
    cd("~/pipeline/evolutionary_rate_analysis/");
    run("Make directory for this output", "mkdir output_$mode");
    run("Move", "mv log-* output_$mode");
    run("Move", "mv tmp-* output_$mode");
    run("Move", "mv output-* output_$mode");
    
    # Store results from evorate quantiles for this MAFFT mode / tree combination
    cd("~/pipeline/evolutionary_rate_quantile_logos/");
    run("Make directory for this output", "mkdir output_$mode");
    run("Move", "mv log-* output_$mode");
    run("Move", "mv tmp-* output_$mode");
    run("Move", "mv output-* output_$mode");
    run("Move", "mv tmp output_$mode");
    run("Move", "mv output output_$mode");
    run("mkdir", "mkdir tmp");
    run("mkdir", "mkdir output");
	
	# Rename type in evorate table
	Query("UPDATE evorate SET type='".$evorate."_".$mode."' WHERE type='$evorate'");
}
