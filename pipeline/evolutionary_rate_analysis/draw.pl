#!/usr/bin/env perl -w

require('functions.inc.pl');
require('mysql.inc.pl');

# initialize

$colours = '"P" = "#4f81bd", "C" = "#bec1c0"';					# BLUE vs. GREY

$densityadjust = '1';											# kernel density estimate scaling factor (default: 1)

$significance = 'pmw_nonresampled_sig';						# Type of p-value to use to show significance in figures
# $significance = 'pmw_nonresampled_averaged_median_sig';		# Type of p-value to use to show significance in figures
# $significance = 'pmwsig';										# Type of p-value to use to show significance in figures
# $significance = 'pmediansig';									# Type of p-value to use to show significance in figures
# $significance = 'pmeansig';										# Type of p-value to use to show significance in figures

$directionality = 'median';									# use mean or median to determine directionality for *** vs. (***)? (only takes effect during recalc)
# $directionality = 'mean';										# use mean or median to determine directionality for *** vs. (***)? (only takes effect during recalc)

$lineplace = 'median';										# place vertical line at mean or median?
# $lineplace = 'mean';											# place vertical line at mean or median?

# $linetype = "solid";											# line type for median lines in density plots
# $linetype = "longdash";										# line type for median lines in density plots
$linetype = "dashed";											# line type for median lines in density plots
# $linetype = "dotted";											# line type for median lines in density plots

# $linesize = ', size=1, lineend="round"';						# size (thickness) for median lines in density plots
# $linesize = ', size=0.5';										# size (thickness) for median lines in density plots
$linesize = '';													# size (thickness) for median lines in density plots

our $usage = "$0 [species] [source: UniProt/Ochoa/PhosphoSitePlus/dbPTM/all] [evorate type] [disorder predictor] [bootstraps (e.g. 100)/all/averaged] [-allsites] [-bonferroni] [-nolog] [-noreplace] [-norecalc]\n\n".
" -allsites:     Show 'All sites' as a category (in addition to the breakdown).\n".
" -bonferroni:   Correct p-values for multiple testing (number of PTMs, not plots, since plots are just breakdowns).\n".
" -dismerge:     Merge disordered & structured (only distinguish by surface/core)\n".
# " -fastdebug:    Only do 100 repeats on ks.boot test\n".
# " -fastestdebug: Only do 10 repeats on ks.boot test\n".
" -filtered:     Don't expect surface/core information in input dataframe (only str/dis), but read from main tmp-dataframe (not -nosurf)\n".
" -nolog:        Don't log-transform the evolutionary rates (default for normalised values).\n".
" -noreplace     Sample controls without replacement, and skip cases if there aren't enough controls.\n".
" -nosurf:       Don't expect surface/core information in input dataframe (only str/dis).\n".
" -recalc:		 Recalculate all p-values, rather than just using an existing p-value R object if possible.\n".
" -surfmerge:    Merge surface & core (only distinguish by disorder)\n".
"\nExample: $0 human UniProt rate4site_einsi_tree AlphaFold";

# ($species, $source, $evorate, $predictor) = args(4);
($species, $source, $evorate, $predictor, $boots) = args(5);



$controlsamplingfunction = "sample_n_randomly_with_replacement";	# with replacement
$noreplacestr = '';
if (switch('noreplace'))
{
	$controlsamplingfunction = "sample_n_uniques_randomly";			# without replacement and skip cases if not enough controls
	$noreplacestr = '-noreplace';
}

$logstr = '';
if (switch('nolog'))
{
	$logstr = '-nolog';
}

$nosurfstr = '';
if (switch('nosurf'))
{
	$nosurfstr = '-nosurf';
}

$coresurfstr = '';
if (switch('coresurf'))
{
	$coresurfstr = '-coresurf';
}

$mergestr = '';
if (switch('dismerge') and switch('surfmerge'))
{
	die("Error: Switches -dismerge and -surfmerge can't both be active");
}
elsif (switch('dismerge'))
{
	$mergestr = '-dismerge';
	$coresurfstr = '-coresurf';
}
elsif (switch('surfmerge'))
{
	$mergestr = '-surfmerge';
	$coresurfstr = '-coresurf';
}
elsif (switch('coresurf'))
{
	$mergestr = '-coresurf';
	$coresurfstr = '-coresurf';
}

# $function = 'mean';
# if (switch('median'))
# {
# 	$function = 'median';
# }

# $infile = "tmp/tmp-dataframe-$boots-$species-$source-$evorate-$predictor$nosurfstr$noreplacestr.txt";
# $infile = "tmp/tmp-dataframe-$boots-$species-$source-$evorate-$predictor$coresurfstr$nosurfstr$noreplacestr.txt";
$infile = "tmp/tmp-dataframe-all-$species-$source-$evorate-$predictor$coresurfstr$nosurfstr$noreplacestr.txt";
$tmppvals = "output/tmp-R-pvalues-$boots-$species-$source-$evorate-$predictor$mergestr$nosurfstr$noreplacestr.txt";
$tmpevo = "output/tmp-R-processed-evoratedata-$boots-$species-$source-$evorate-$predictor$mergestr$nosurfstr$noreplacestr.txt";
$outdens = "output/output-densityplots-$boots-$species-$source-$evorate-$predictor$mergestr$logstr$nosurfstr$noreplacestr.pdf";
$outtable = "output/output-pvalues-$boots-$species-$source-$evorate-$predictor$mergestr$nosurfstr$noreplacestr.tsv";
$outevo = "output/output-processed-evoratedata-$boots-$species-$source-$evorate-$predictor$mergestr$nosurfstr$noreplacestr.tsv";
$outtest = "output/output-teststatistic-$boots-$species-$source-$evorate-$predictor$mergestr$nosurfstr$noreplacestr.tsv";
# $outbox = "output/output-boxplots-$boots-$species-$source-$evorate-$predictor$mergestr$logstr$nosurfstr$noreplacestr.pdf";

if (!-s $infile)
{
	die("Error: Couldn't open '$infile'");
}


# start

state("Drawing plots based on '$infile'");

$xlabel = "Variation score";
if ($evorate =~ /^lichtarge/)         { $xlabel = "Variation score (rvET)"; }
elsif ($evorate =~ /^capra/)          { $xlabel = "Variation score (inverse Jensen-Shannon divergence)"; }
elsif ($evorate =~ /^rate4site/)      { $xlabel = "Evolutionary rate (Rate4Site)"; }
elsif ($evorate =~ /^norm_lichtarge/) { $xlabel = "Normalised variation score (rvET)"; }
elsif ($evorate =~ /^norm_capra/)     { $xlabel = "Normalised variation score (inverse Jensen-Shannon divergence)"; }
elsif ($evorate =~ /^norm_rate4site/) { $xlabel = "Normalised evolutionary rate (Rate4Site)"; }

state("Initialising (loading libraries and data frame)...");
starttime();
startr();
state(runr('

options(warn = 1)		# always print R warnings
options(nwarnings = 10000)

# Quiet tidyverse startup message
options(tidyverse.quiet = T)
# Quiet dplyr summarise .groups message
options(dplyr.summarise.inform = F)

library("tidyverse")		# for plotting plots
library("scales")			# for plotting scales
library("glue")
library("magrittr")			# for pipes
# Not any faster than tibble filter lookups for accs
# library("data.table") 
# # Create an in-memory SQLite database
# con <- dbConnect(RSQLite::SQLite(), ":memory:")

', 1), 1);



# Set "dis" (disorder/residue type) plotting order
$dis_order = 'dis';

if (switch('dismerge'))
{
	# $dis_order = 'fct_relevel(dis, "All sites", "Core sites", "Surface sites")';
	# $dis_order = 'fct_relevel(dis, "Surface sites", "Core sites")';
	$dis_order = 'fct_relevel(dis, "Core sites", "Surface sites")';
}
elsif (switch('surfmerge'))
{
	# $dis_order = 'fct_relevel(dis, "All sites", "Disordered sites", "Structured sites")';
	# $dis_order = 'fct_relevel(dis, "Disordered sites", "Structured sites")';
	$dis_order = 'fct_relevel(dis, "Structured sites", "Disordered sites")';
}
elsif (switch('nosurf') or switch('filtered'))
{
	# $dis_order = 'fct_relevel(dis, "All sites", "Disordered sites", "Structured sites")';
	# $dis_order = 'fct_relevel(dis, "Disordered sites", "Structured sites")';
	$dis_order = 'fct_relevel(dis, "Structured sites", "Disordered sites")';
}
elsif (switch('coresurf'))
{
	# $dis_order = 'fct_relevel(dis, "All sites", "Buried", "Surface Structured", "Surface Disordered")';
	# $dis_order = 'fct_relevel(dis, "Buried", "Surface Structured", "Surface Disordered")';
	# $dis_order = 'fct_relevel(dis, "Surface Disordered", "Surface Structured", "Buried")';
	$dis_order = 'fct_relevel(dis, "Buried", "Surface Structured", "Surface Disordered")';
}
else
{
	# $dis_order = 'fct_relevel(dis, "All sites", "Disordered sites", "Structured sites")';
	# $dis_order = 'fct_relevel(dis, "Disordered sites", "Structured sites")';
	$dis_order = 'fct_relevel(dis, "Structured sites", "Disordered sites")';
}







if (switch('recalc') or (!-s $tmppvals))
{
	
	
	state(runr('
	# install.packages("doBy")
	library(doBy) 

	# evoratedata <- read.delim("'.$infile.'", header=T, quote="", stringsAsFactors = FALSE)
	evoratedata <- read.delim("'.$infile.'", header=T, quote="")
	evoratedata <- as_tibble(evoratedata)

	# evoratedata %<>% subset(ptm == "C-nit" & dis == "dissurf")

	# Convert all character columns to factors
	evoratedata %<>% mutate_if(is.character, as.factor)
	evoratedata

# 	# evoratedata$name <- NULL
# 	evoratedata$site <- NULL
# 
# 	# evoratedata <- evoratedata[evoratedata$dis != "A", ]
# 	levels(evoratedata$type) <- list("Control"="C", "PTM"="P")
# 	# evoratedata$rate <- log2(evoratedata$rate)

	', 1), 1);
	stoptime();
	
	if (switch('dismerge'))
	{
		state(runr('
		str(evoratedata$dis)
		evoratedata$dis <- as.character(evoratedata$dis)
		try(evoratedata[evoratedata$dis == "discore",]$dis <- "core", silent=T)
		try(evoratedata[evoratedata$dis == "dissurf",]$dis <- "surf", silent=T)
		try(evoratedata[evoratedata$dis == "strcore",]$dis <- "core", silent=T)
		try(evoratedata[evoratedata$dis == "strsurf",]$dis <- "surf", silent=T)
		evoratedata$dis <- as.factor(evoratedata$dis)
		str(evoratedata$dis)

		# levels(evoratedata$dis) <- list("Core"="core", "Surface"="surf", "Combined"="A")
		# levels(evoratedata$dis) <- list("All sites"="A", "Core"="core", "Surface"="surf")
		levels(evoratedata$dis) <- list("All sites"="A", "Core sites"="core", "Surface sites"="surf")
		str(evoratedata$dis)
		', 1), 1);
	}
	elsif (switch('surfmerge'))
	{
		state(runr('
		evoratedata$dis <- as.character(evoratedata$dis)
		try(evoratedata[evoratedata$dis == "discore",]$dis <- "dis", silent=T)
		try(evoratedata[evoratedata$dis == "dissurf",]$dis <- "dis", silent=T)
		try(evoratedata[evoratedata$dis == "strcore",]$dis <- "str", silent=T)
		try(evoratedata[evoratedata$dis == "strsurf",]$dis <- "str", silent=T)
		evoratedata$dis <- as.factor(evoratedata$dis)

		# levels(evoratedata$dis) <- list("Structured"="str", "Disordered"="dis", "Combined"="A")
		# levels(evoratedata$dis) <- list("All sites"="A", "Structured"="str", "Disordered"="dis")
		levels(evoratedata$dis) <- list("All sites"="A", "Disordered sites"="dis", "Structured sites"="str")
		str(evoratedata$dis)
		', 1), 1);
	}
	elsif (switch('nosurf') or switch('filtered'))
	{
		state(runr('
		str(evoratedata$dis)
		evoratedata$dis <- as.factor(evoratedata$dis)
		# levels(evoratedata$dis) <- list("Structured"="str", "Disordered"="dis", "Combined"="A")
		# levels(evoratedata$dis) <- list("All sites"="A", "Structured"="str", "Disordered"="dis")
		levels(evoratedata$dis) <- list("All sites"="A", "Disordered sites"="dis", "Structured sites"="str")
		str(evoratedata$dis)
		', 1), 1);
	}
	elsif (switch('coresurf'))
	{
		state(runr('
		str(evoratedata$dis)
		# levels(evoratedata$dis) <- list("Core / Structured"="strcore", "Core / Disordered"="discore", "Surface / Structured"="strsurf", "Surface / Disordered"="dissurf", "Combined"="A")
		# levels(evoratedata$dis) <- list("All sites"="A", "Core / Structured"="strcore", "Core / Disordered"="discore", "Surface / Structured"="strsurf", "Surface / Disordered"="dissurf")
		# levels(evoratedata$dis) <- list("All sites"="A", "Core / Disordered"="discore", "Core / Structured"="strcore", "Surface / Disordered"="dissurf", "Surface / Structured"="strsurf")
		# levels(evoratedata$dis) <- list("All sites"="A", "Disordered sites"="dis", "Structured sites"="str")

		# Consider Core / Disordered (very rare) to be Core / Structured (since Disorder uses a +/- 10 aa window, whereas Core/Surface does not - so Core / Disordered are clearly Core residues where the "Disorder" call comes from adjacent residues)
		# evoratedata <- evoratedata[evoratedata$dis != "discore", ]
		evoratedata$dis <- as.character(evoratedata$dis)
		try(evoratedata[evoratedata$dis == "discore",]$dis <- "strcore", silent=T)
		evoratedata$dis <- as.factor(evoratedata$dis)
		str(evoratedata$dis)

		# levels(evoratedata$dis) <- list("Core / Structured"="strcore", "Surface / Structured"="strsurf", "Surface / Disordered"="dissurf", "Combined"="A")
		# levels(evoratedata$dis) <- list("Core / Structured"="strcore", "Surface / Structured"="strsurf", "Surface / Disordered"="dissurf", "Combined"="A")
		# levels(evoratedata$dis) <- list("All sites"="A", "Core / Structured"="strcore", "Surface / Structured"="strsurf", "Surface / Disordered"="dissurf")
		print(levels(evoratedata$dis))
		levels(evoratedata$dis) <- list("All sites"="A", "Buried"="strcore", "Surface Structured"="strsurf", "Surface Disordered"="dissurf")
		str(evoratedata$dis)
		', 1), 1);
	}
	else
	{
		state(runr('
		str(evoratedata$dis)
		evoratedata$dis <- as.factor(evoratedata$dis)
		# levels(evoratedata$dis) <- list("Core / Structured"="strcore", "Core / Disordered"="discore", "Surface / Structured"="strsurf", "Surface / Disordered"="dissurf", "Combined"="A")
		# levels(evoratedata$dis) <- list("All sites"="A", "Core / Structured"="strcore", "Core / Disordered"="discore", "Surface / Structured"="strsurf", "Surface / Disordered"="dissurf")
		# levels(evoratedata$dis) <- list("All sites"="A", "Core / Disordered"="discore", "Core / Structured"="strcore", "Surface / Disordered"="dissurf", "Surface / Structured"="strsurf")
		levels(evoratedata$dis) <- list("All sites"="A", "Disordered sites"="dis", "Structured sites"="str")
		str(evoratedata$dis)

		# Remove Core / Disordered since it doesnt occur often enough?
		# evoratedata <- evoratedata[evoratedata$dis != "discore", ]
		# levels(evoratedata$dis) <- list("Core / Structured"="strcore", "Surface / Structured"="strsurf", "Surface / Disordered"="dissurf", "Combined"="A")
		# levels(evoratedata$dis) <- list("All sites"="A", "Core / Structured"="strcore", "Surface / Structured"="strsurf", "Surface / Disordered"="dissurf")
		', 1), 1);
	}

	# Calculate p-values
	state("Calculating p-values and writing them to '$outtable'...");
	starttime();
	state(runr('
	# library(Matching)		# for ks.boot test (Kolmogorov-Smirnov with 1000x bootstrapping, allows for ties in contrast to regular KS test)

	pvalues <- tibble()
	teststatistic <- tibble()
	
	evoratedata
	
	# print("Unique C-nit sites:")
	# # evoratedata %>% subset(ptm == "C-nit" & dis == "Surface Disordered") %>% select(acc, site) %>% unique %>% print
	# evoratedata %>% subset(ptm == "C-nit") %>% select(acc, site) %>% unique %>% print

	str(evoratedata$ptm)
	head(evoratedata$ptm)
	
	str(evoratedata$dis)
	head(evoratedata$dis)
	
	str(evoratedata$type)
	head(evoratedata$type)

	boots <- '.$boots.'

		', 1), 1);

	state(runr('

	# Custom function for sampling with replacement
	sample_with_replacement <- function(data_subset, sample_size) {
	  if(!is.null(data_subset) && nrow(data_subset) > 0) {
	    sampled_rows <- sample.int(n = nrow(data_subset), size = sample_size, replace = TRUE)
	    return(data_subset[sampled_rows, ])
	  } else {
	    return(data.frame())
	  }
	}

	'), 1);

# 	state(runr('
# 
# 	sample_one_randomly_and_return_it_n_times <- function (a, n) {  
# 	  b <- sample(a, 1)
# 	  return(rep(b, n))
# 	}
# 
# 	sample_n_uniques_randomly <- function (a, n) {
# 	  if (length(a) >= n) {
# 	    return(sample(a, n, replace=FALSE))
# 	  }
# 	  else {
# 	    return(rep(NA, n))
# 	  }
# 	}
# 
# 	sample_n_randomly_with_replacement <- function (a, n) {
# 	    return(sample(a, n, replace=TRUE))
# 	}
# 
# 	'), 1);
	
	state(runr(q(

	for (x in unique(evoratedata$ptm)) {
		cat(glue("x is now: {x}"))
		for (y in sort(unique(evoratedata$dis))) {
			cat(glue("y is now: {y}"))

			# ptmr <- subset(evoratedata, ptm == x & dis == y & evoratedata$type == "P")$rate
			# controlr <- subset(evoratedata, ptm == x & dis == y & evoratedata$type == "C")$rate
			
			ptmr <- subset(evoratedata, ptm == x & dis == y & evoratedata$type == "P")
			controlr <- subset(evoratedata, ptm == x & dis == y & evoratedata$type == "C")
			
			
			# Limit to proteins that have both PTM and control residues
			# i.e. skip any proteins that do not have control residues
			# length(unique(ptmr$acc))
			# length(unique(controlr$acc))
			# length(intersect(ptmr$acc, controlr$acc))
			# cat(glue("\nBEFORE filtering out proteins without control residues:\n\n"))
			# cat(ptmr %>% filter(ptm == 'C-nit' & dis == 'Surface Disordered' & type == "P") %>% select(acc, site) %>% unique %>% nrow)
			# cat(glue("\n\n"))
			ptmr <- subset(ptmr, acc %in% intersect(ptmr$acc, controlr$acc))
			controlr <- subset(controlr, acc %in% intersect(ptmr$acc, controlr$acc))
			# cat(glue("\nAFTER filtering out proteins without control residues:\n\n"))
			# cat(ptmr %>% filter(ptm == 'C-nit' & dis == 'Surface Disordered' & type == "P") %>% select(acc, site) %>% unique %>% nrow)
			# cat(glue("\n\n"))
			
			# Split controlr by 'acc'
			# controlr_split <- split(controlr, controlr$acc) # Using a factor to split doesn't work since unused levels get used as well
			controlr_split <- split(controlr, as.character(controlr$acc))
			# names(controlr_split) %>% unique %>% length
			# ptmr$acc %>% unique %>% length
			# controlr$acc %>% unique %>% length

			# Basic non-resampled p-values
			wilcox_nonresampled <- wilcox.test(ptmr$rate, controlr$rate, alternative="two.sided", paired=FALSE)
			ttest_nonresampled <- t.test(ptmr$rate, controlr$rate, alternative="two.sided", paired=FALSE)

			# Calculate means and medians for directionality
			cat(glue("\n\nPTM median and mean for {x} {y}:\n\n"))
			print(ptmr %>% summarise(median_rate = median(rate), mean_rate = mean(rate)))
			cat(glue("\n\nControl median and mean for {x} {y}:\n\n"))
			print(controlr %>% summarise(median_rate = median(rate), mean_rate = mean(rate)))
			
			# Average by accession
			ptmr %>% group_by(ptm, acc, type, dis) %>% summarise(median_rate = median(rate), mean_rate = mean(rate)) %>% ungroup %>% arrange(acc) -> ptmr_averaged_by_acc
			controlr %>% group_by(ptm, acc, type, dis) %>% summarise(median_rate = median(rate), mean_rate = mean(rate)) %>% ungroup %>% arrange(acc) -> controlr_averaged_by_acc

			# Verify that order is correct
			identical(ptmr_averaged_by_acc %>% select(ptm, acc, dis), controlr_averaged_by_acc %>% select(ptm, acc, dis))
			if (!identical(ptmr_averaged_by_acc %>% select(ptm, acc, dis), controlr_averaged_by_acc %>% select(ptm, acc, dis))) {
			stop("Error: ptm|acc|dis orders are not identical between ptmr_averaged_by_acc and controlr_averaged_by_acc")
			}
			
			# Paired
			# wilcox.test(ptmr_averaged_by_acc$median_rate, controlr_averaged_by_acc$median_rate, alternative="two.sided", paired=TRUE)$p.value
			# wilcox.test(ptmr_averaged_by_acc$mean_rate, controlr_averaged_by_acc$mean_rate, alternative="two.sided", paired=TRUE)$p.value
			wilcox_nonresampled_averaged_median <- wilcox.test(ptmr_averaged_by_acc$median_rate, controlr_averaged_by_acc$median_rate, alternative="two.sided", paired=TRUE)
			wilcox_nonresampled_averaged_mean <- wilcox.test(ptmr_averaged_by_acc$mean_rate, controlr_averaged_by_acc$mean_rate, alternative="two.sided", paired=TRUE)
			
			# Unpaired
			# wilcox.test(ptmr_averaged_by_acc$median_rate, controlr_averaged_by_acc$median_rate, alternative="two.sided", paired=FALSE)$p.value
			# wilcox.test(ptmr_averaged_by_acc$mean_rate, controlr_averaged_by_acc$mean_rate, alternative="two.sided", paired=FALSE)$p.value
			wilcox_nonresampled_averaged_median_unpaired <- wilcox.test(ptmr_averaged_by_acc$median_rate, controlr_averaged_by_acc$median_rate, alternative="two.sided", paired=FALSE)
			wilcox_nonresampled_averaged_mean_unpaired <- wilcox.test(ptmr_averaged_by_acc$mean_rate, controlr_averaged_by_acc$mean_rate, alternative="two.sided", paired=FALSE)



			# tibble(ptmr)
			# tibble(ptmr) %>% nrow
			# tibble(ptmr) %>% pull(acc) %>% unique %>% length
			# tibble(ptmr) %>% group_by(acc)
			# tibble(ptmr) %>% group_by(acc, site)
			# 
			# tibble(controlr)
			
			# # Convert controlr to data.table for faster acc lookup (data.table has indexes)
			# Not actually any faster than tibble acc lookup using filter
			# controlr <- as.data.table(controlr)
			# setkey(controltest, "acc")
			# 
			# # sqlite
			# Not actually faster either
			# # Load tibble into sqlite table
			# dbWriteTable(con, "controlr", controlr)
			# 
			# Benchmarking:
			# tibble:
			# controltest <- tibble(controlr); system.time(replicate(100, { controltest %>% filter(acc == "P51790") }))
			# >> ~200 ms
			# data.table:
			# controltest <- as.data.table(controlr)
			# setkey(controlr, "acc")
			# system.time(replicate(100, { controltest %>% filter(acc == "P51790") }))
			# >> ~200 ms (same as tibble - index makes no difference)
			# dbplyr (sqlite):
			# system.time(replicate(100, { tbl(con, "controlr") %>% filter(acc == "P51790") %>% collect() }))
			# >> ~2900 ms (super slow)
			# >> Just use tibble

			# Bootstrap

			mean_successes <- 0
			median_successes <- 0
			wilcox_successes <- 0
			ttest_successes <- 0

			for (i in 1:boots) {
				cat(glue(" >> {x} {y} sampling cycle {i}"))
				
				# Resample PTM sites
				ptm_sample <- sample_n(tibble(ptmr), size = nrow(ptmr), replace = T)
				# ptm_sample %>% arrange(acc, site)
				# tibble(ptmr) %>% arrange(acc, site)
				
				# Resample control sites proportionally from the same proteins as in the PTM sample

				# tibble(ptm_sample$acc)
				# tibble(controlr)
				# control_sample <- map(ptm_sample$acc, ~ controlr %>% filter(acc == .x) %>% sample_n(size = 1))
				# control_sample <- lapply(ptm_sample$acc, ~ controlr %>% filter(acc == .x) %>% sample_n(size = 1))
				# tibble(control_sample)
				
				# ptm_sample
				# ptm_sample %>% group_by(acc) %>% tally
				# ptm_sample %>% group_by(acc) %>% tally
				
				# for (my_acc, my_n in ptm_sample %>% group_by(acc) %>% tally) {
				#   print(my_acc, my_n)
				# }
				
				# Working method, but slightly slow (~25 secs):
				
				# Tally the number of rows for each acc in ptm_sample
				acc_counts <- ptm_sample %>% group_by(acc) %>% tally
				
				# # For each acc, sample rows from controlr
				# print(system.time(control_sample_list <- map2(acc_counts$acc, acc_counts$n,
				# 							~controlr %>%
				# 								filter(acc == .x) %>%
				# 								sample_n(size = .y, replace = TRUE))))

			    # Optimised approach:
			
			    # Create a named vector of sample sizes
			    sample_sizes <- setNames(acc_counts$n, acc_counts$acc)
			
			    # Ensure that every 'acc' in sample_sizes has a corresponding entry in controlr_split
			    # Fill missing 'acc' entries in controlr_split with empty data frames
			    my_controlr_split <- lapply(names(sample_sizes), function(acc) {
			    	if (acc %in% names(controlr_split)) {
						controlr_split[[acc]]
					} else {
						data.frame()
					}
				})

				# Main resampling step:
			    # Revised control sampling step
			    print(system.time(control_sample_list <- map2(my_controlr_split, acc_counts$n, sample_with_replacement)))
			
			    # Combine the list of sampled data frames into one data frame
			    control_sample <- bind_rows(control_sample_list)

				# sqlite experiment:
				# # For each acc, sample rows from controlr
				# print(system.time(control_sample_list <- map2(acc_counts$acc, acc_counts$n, 
				#                             ~tbl(con, "controlr") %>%
				#                               filter(acc == .x) %>% 
				#                               collect() %>% 
				#                               slice_sample(n = .y, replace = TRUE))))
				
				# Convert the list of tibbles into one tibble
				control_sample <- tibble(bind_rows(control_sample_list))
				
				# # Slightly faster method (15 secs) that doesn't work properly:
				# 
				# # Assuming acc_counts has been created as you suggested:
				# acc_counts <- ptm_sample %>% group_by(acc) %>% tally() %>% rename(num = n) %>% arrange(desc(num))
				# acc_counts
				# 
				# # Join with controlr, compute the sampling size, and then sample
				# system.time(control_sample <- map2_dfr(acc_counts$acc, acc_counts$num, ~ {
				#   acc_value <- .x
				#   sample_size <- .y
				#   available_rows <- nrow(filter(controlr, acc == acc_value))
				#   sample_rows <- min(sample_size, available_rows)
				#   
				#   controlr %>%
				#     filter(acc == acc_value) %>%
				#     slice_sample(n = sample_rows)
				# }))
				
				# Verify that the number of sites per acc is identical between ptm_sample and control_sample
				# acc_counts_control <- control_sample %>% group_by(acc) %>% tally() %>% rename(num = n) %>% arrange(desc(num))
				acc_counts_control <- control_sample %>% group_by(acc) %>% tally()
				# acc_counts_control
				# acc_counts
				# identical(acc_counts, acc_counts_control)
				if (!identical(acc_counts, acc_counts_control)) {
					stop("Error: acc_counts and acc_counts_control are not identical")
				}
				
				# Sort tibbles
				print(head(ptm_sample))
				ptm_sample %<>% arrange(acc, site)
				print(head(control_sample))
				control_sample %<>% arrange(acc, site)
				
				# Combine tibbles
				# ptm_sample
				# control_sample
				# Verify that order is correct
				identical(ptm_sample %>% select(ptm, acc, dis), control_sample %>% select(ptm, acc, dis))
				if (!identical(ptm_sample %>% select(ptm, acc, dis), control_sample %>% select(ptm, acc, dis))) {
					stop("Error: ptm|acc|dis orders are not identical between ptm_sample and control_sample")
				}
				# Combine in order to do paired wilcox.test below
				combined_sample <- ptm_sample
				combined_sample %<>% transmute(ptm, acc, dis, rate_ptm = rate)
				combined_sample$rate_control <- control_sample$rate
				# combined_sample

				# Count "successful" bootstraps
				# Wilcoxon signed-rank test (i.e. paired wilcox.test), one-tailed
				# ?wilcox.test
				# summary(combined_sample)
				sample_wilcox <- wilcox.test(combined_sample$rate_ptm, combined_sample$rate_control, alternative = "less", paired = T)
				if (sample_wilcox$p.value < 0.05) {
					wilcox_successes <- wilcox_successes + 1
				}
				# wilcox_successes
				# Student's t-test, paired, one-tailed
				?t.test
				sample_ttest <- t.test(combined_sample$rate_ptm, combined_sample$rate_control, alternative = "less", paired = T)
				if (sample_ttest$p.value < 0.05) {
					ttest_successes <- ttest_successes + 1
				}
				# ttest_successes
				
				
				# Calculate means and medians
				# Combine into one tidy tibble
				combined_sample_mean_and_median <- bind_rows(ptm_sample, control_sample)
				combined_sample_mean_and_median
				ptm_sample %>% group_by(type, dis) %>% tally
				control_sample %>% group_by(type, dis) %>% tally
				combined_sample_mean_and_median %>% group_by(type, dis) %>% tally
				combined_sample_mean_and_median %<>% group_by(ptm, type, dis) %>% summarise(mean_rate = mean(rate), median_rate = median(rate))
				# pivot_wider
				combined_sample_mean_and_median %<>% pivot_wider(id_cols = c(ptm, dis), names_from = type, values_from = c(mean_rate, median_rate))
				combined_sample_mean_and_median
				# Count "successful" bootstraps
				# First, verify that we only have a single tibble row (i.e. only one PTM-dis combination)
				if (nrow(combined_sample_mean_and_median) != 1) {
					stop(glue("Error: combined_sample_mean_and_median has nrow {nrow(combined_sample_mean_and_median)}"))
				}
				# Mean
				if (combined_sample_mean_and_median$mean_rate_P < combined_sample_mean_and_median$mean_rate_C) {
					mean_successes <- mean_successes + 1
				}
				# Median
				if (combined_sample_mean_and_median$median_rate_P < combined_sample_mean_and_median$median_rate_C) {
					median_successes <- median_successes + 1
				}
				
				# Record values from this bootstrap cycle
				teststatistic <- bind_rows(teststatistic, tibble_row(ptm=x, dis=y, i=i, nptm=nrow(ptm_sample), ncontrol=nrow(control_sample), medianptm=median(ptm_sample$rate), mediancontrol=median(control_sample$rate), meanptm=mean(ptm_sample$rate), meancontrol=mean(control_sample$rate), pmw=sample_wilcox$p.value, ptt=sample_ttest$p.value, mediandif=median(ptm_sample$rate)-median(control_sample$rate), meandif=mean(ptm_sample$rate)-mean(control_sample$rate)))
				# write.table(teststatistic, file="'.$outtest.'", sep="\t", quote=FALSE, row.names=FALSE)
				print(glue(" >> Test statistics from resample {i}:"))
				print(teststatistic)
			}

			# Calculate final resampling p-values (across resamples)
			teststatistic %>% filter(ptm == x & dis == y) -> mytest
			pvalues <- bind_rows(pvalues, tibble_row(ptm=x, dis=y, nptm=unique(mytest$nptm), ncontrol=unique(mytest$ncontrol), medianptm=median(mytest$medianptm), mediancontrol=median(mytest$mediancontrol), meanptm=mean(mytest$meanptm), meancontrol=mean(mytest$meancontrol), pmw_nonresampled=wilcox_nonresampled$p.value, pmw_nonresampled_averaged_median=wilcox_nonresampled_averaged_median$p.value, pmw_nonresampled_averaged_mean=wilcox_nonresampled_averaged_mean$p.value, pmw_nonresampled_averaged_median_unpaired=wilcox_nonresampled_averaged_median_unpaired$p.value, pmw_nonresampled_averaged_mean_unpaired=wilcox_nonresampled_averaged_mean_unpaired$p.value, ptt_nonresampled=ttest_nonresampled$p.value, pmw=1-(wilcox_successes/boots), ptt=1-(ttest_successes/boots), pmedian=1-(median_successes/boots), pmean=1-(mean_successes/boots), mediandif=round(median(mytest$mediandif), 3), meandif=round(mean(mytest$meandif), 3), resamples=boots, wilcox_successes=wilcox_successes, ttest_successes=ttest_successes, median_successes=median_successes, mean_successes=mean_successes))
		}
	}

# 	teststatistic
# 	pvalues
# 	pvalues %>% filter(medianptm < mediancontrol)
# 	pvalues %>% filter(meanptm < meancontrol)
# 	pvalues %>% filter(ptt < 0.05)
# 	pvalues %>% filter(pmw < 0.05)
# 	pvalues
	print(" >> Final p-values:")
	print(pvalues)

	)), 1);
	
	state(runr('
	# Calculate Bonferroni correction and add as extra columns
	pvalues$pmw_nonresampled_corr <- min(pvalues$pmw_nonresampled * nrow(pvalues), 1)
	pvalues$pmw_nonresampled_averaged_median_corr <- min(pvalues$pmw_nonresampled_averaged_median * nrow(pvalues), 1)
	pvalues$pmw_nonresampled_averaged_mean_corr <- min(pvalues$pmw_nonresampled_averaged_mean * nrow(pvalues), 1)
	pvalues$pmw_nonresampled_averaged_median_unpaired_corr <- min(pvalues$pmw_nonresampled_averaged_median_unpaired * nrow(pvalues), 1)
	pvalues$pmw_nonresampled_averaged_mean_unpaired_corr <- min(pvalues$pmw_nonresampled_averaged_mean_unpaired * nrow(pvalues), 1)
	pvalues$ptt_nonresampled_corr <- min(pvalues$ptt_nonresampled * nrow(pvalues), 1)
	pvalues$pmwcorr <- min(pvalues$pmw * nrow(pvalues), 1)
	pvalues$pttcorr <- min(pvalues$ptt * nrow(pvalues), 1)
	pvalues$pmediancorr <- min(pvalues$pmedian * nrow(pvalues), 1)
	pvalues$pmeancorr <- min(pvalues$pmean * nrow(pvalues), 1)
	', 1), 1);

	if (switch('bonferroni'))
	{
		state(runr('

		# Save uncorrected p-values
		pvalues$pmw_nonresampled_nocorr <- pvalues$pmw_nonresampled
		pvalues$pmw_nonresampled_averaged_median_nocorr <- pvalues$pmw_nonresampled_averaged_median
		pvalues$pmw_nonresampled_averaged_mean_nocorr <- pvalues$pmw_nonresampled_averaged_mean
		pvalues$pmw_nonresampled_averaged_median_unpaired_nocorr <- pvalues$pmw_nonresampled_averaged_median_unpaired
		pvalues$pmw_nonresampled_averaged_mean_unpaired_nocorr <- pvalues$pmw_nonresampled_averaged_mean_unpaired
		pvalues$ptt_nonresampled_nocorr <- pvalues$ptt_nonresampled
		pvalues$pmwnocorr <- pvalues$pmw
		pvalues$pttnocorr <- pvalues$ptt
		pvalues$pmediannocorr <- pvalues$pmedian
		pvalues$pmeannocorr <- pvalues$pmean

		# Apply Bonferroni correction to p-values
		pvalues$pmw_nonresampled <- pvalues$pmw_nonresampled_corr
		pvalues$pmw_nonresampled_averaged_median <- pvalues$pmw_nonresampled_averaged_median_corr
		pvalues$pmw_nonresampled_averaged_mean <- pvalues$pmw_nonresampled_averaged_mean_corr
		pvalues$pmw_nonresampled_averaged_median_unpaired <- pvalues$pmw_nonresampled_averaged_median_unpaired_corr
		pvalues$pmw_nonresampled_averaged_mean_unpaired <- pvalues$pmw_nonresampled_averaged_mean_unpaired_corr
		pvalues$ptt_nonresampled <- pvalues$ptt_nonresampled_corr
		pvalues$pmw <- pvalues$pmwcorr
		pvalues$ptt <- pvalues$pttcorr
		pvalues$pmedian <- pvalues$pmediancorr
		pvalues$pmean <- pvalues$pmeancorr
		', 1), 1);
	}

	state(runr('
	pvalues[,"pmw_nonresampled_sig"] <- ""
	try(pvalues[pvalues$pmw_nonresampled < 0.05,]$pmw_nonresampled_sig <- "*", silent=T)
	try(pvalues[pvalues$pmw_nonresampled < 0.01,]$pmw_nonresampled_sig <- "**", silent=T)
	try(pvalues[pvalues$pmw_nonresampled < 0.001,]$pmw_nonresampled_sig <- "***", silent=T)
	# pvalues$pmw_nonresampled_sig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$pmw_nonresampled_sig != "")] <- paste("+ (", pvalues$pmw_nonresampled_sig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$pmw_nonresampled_sig != "")], ")", sep = "")
	# pvalues$pmw_nonresampled_sig[(pvalues$medianptm < pvalues$mediancontrol) & (pvalues$pmw_nonresampled_sig != "")] <- paste("- (", pvalues$pmw_nonresampled_sig[(pvalues$medianptm < pvalues$mediancontrol) & (pvalues$pmw_nonresampled_sig != "")], ")", sep = "")
	# pvalues$pmw_nonresampled_sig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$pmw_nonresampled_sig != "")] <- paste("(", pvalues$pmw_nonresampled_sig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$pmw_nonresampled_sig != "")], ")", sep = "")
	pvalues$pmw_nonresampled_sig[(pvalues$'.$directionality.'ptm > pvalues$'.$directionality.'control) & (pvalues$pmw_nonresampled_sig != "")] <- paste("(", pvalues$pmw_nonresampled_sig[(pvalues$'.$directionality.'ptm > pvalues$'.$directionality.'control) & (pvalues$pmw_nonresampled_sig != "")], ")", sep = "")

	pvalues[,"pmw_nonresampled_averaged_median_sig"] <- ""
	try(pvalues[pvalues$pmw_nonresampled_averaged_median < 0.05,]$pmw_nonresampled_averaged_median_sig <- "*", silent=T)
	try(pvalues[pvalues$pmw_nonresampled_averaged_median < 0.01,]$pmw_nonresampled_averaged_median_sig <- "**", silent=T)
	try(pvalues[pvalues$pmw_nonresampled_averaged_median < 0.001,]$pmw_nonresampled_averaged_median_sig <- "***", silent=T)
	# pvalues$pmw_nonresampled_averaged_median_sig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$pmw_nonresampled_averaged_median_sig != "")] <- paste("+ (", pvalues$pmw_nonresampled_averaged_median_sig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$pmw_nonresampled_averaged_median_sig != "")], ")", sep = "")
	# pvalues$pmw_nonresampled_averaged_median_sig[(pvalues$medianptm < pvalues$mediancontrol) & (pvalues$pmw_nonresampled_averaged_median_sig != "")] <- paste("- (", pvalues$pmw_nonresampled_averaged_median_sig[(pvalues$medianptm < pvalues$mediancontrol) & (pvalues$pmw_nonresampled_averaged_median_sig != "")], ")", sep = "")
	# pvalues$pmw_nonresampled_averaged_median_sig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$pmw_nonresampled_averaged_median_sig != "")] <- paste("(", pvalues$pmw_nonresampled_averaged_median_sig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$pmw_nonresampled_averaged_median_sig != "")], ")", sep = "")
	pvalues$pmw_nonresampled_averaged_median_sig[(pvalues$'.$directionality.'ptm > pvalues$'.$directionality.'control) & (pvalues$pmw_nonresampled_averaged_median_sig != "")] <- paste("(", pvalues$pmw_nonresampled_averaged_median_sig[(pvalues$'.$directionality.'ptm > pvalues$'.$directionality.'control) & (pvalues$pmw_nonresampled_averaged_median_sig != "")], ")", sep = "")

	pvalues[,"pmw_nonresampled_averaged_mean_sig"] <- ""
	try(pvalues[pvalues$pmw_nonresampled_averaged_mean < 0.05,]$pmw_nonresampled_averaged_mean_sig <- "*", silent=T)
	try(pvalues[pvalues$pmw_nonresampled_averaged_mean < 0.01,]$pmw_nonresampled_averaged_mean_sig <- "**", silent=T)
	try(pvalues[pvalues$pmw_nonresampled_averaged_mean < 0.001,]$pmw_nonresampled_averaged_mean_sig <- "***", silent=T)
	# pvalues$pmw_nonresampled_averaged_mean_sig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$pmw_nonresampled_averaged_mean_sig != "")] <- paste("+ (", pvalues$pmw_nonresampled_averaged_mean_sig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$pmw_nonresampled_averaged_mean_sig != "")], ")", sep = "")
	# pvalues$pmw_nonresampled_averaged_mean_sig[(pvalues$medianptm < pvalues$mediancontrol) & (pvalues$pmw_nonresampled_averaged_mean_sig != "")] <- paste("- (", pvalues$pmw_nonresampled_averaged_mean_sig[(pvalues$medianptm < pvalues$mediancontrol) & (pvalues$pmw_nonresampled_averaged_mean_sig != "")], ")", sep = "")
	# pvalues$pmw_nonresampled_averaged_mean_sig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$pmw_nonresampled_averaged_mean_sig != "")] <- paste("(", pvalues$pmw_nonresampled_averaged_mean_sig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$pmw_nonresampled_averaged_mean_sig != "")], ")", sep = "")
	pvalues$pmw_nonresampled_averaged_mean_sig[(pvalues$'.$directionality.'ptm > pvalues$'.$directionality.'control) & (pvalues$pmw_nonresampled_averaged_mean_sig != "")] <- paste("(", pvalues$pmw_nonresampled_averaged_mean_sig[(pvalues$'.$directionality.'ptm > pvalues$'.$directionality.'control) & (pvalues$pmw_nonresampled_averaged_mean_sig != "")], ")", sep = "")

	pvalues[,"pmw_nonresampled_averaged_median_unpaired_sig"] <- ""
	try(pvalues[pvalues$pmw_nonresampled_averaged_median_unpaired < 0.05,]$pmw_nonresampled_averaged_median_unpaired_sig <- "*", silent=T)
	try(pvalues[pvalues$pmw_nonresampled_averaged_median_unpaired < 0.01,]$pmw_nonresampled_averaged_median_unpaired_sig <- "**", silent=T)
	try(pvalues[pvalues$pmw_nonresampled_averaged_median_unpaired < 0.001,]$pmw_nonresampled_averaged_median_unpaired_sig <- "***", silent=T)
	# pvalues$pmw_nonresampled_averaged_median_unpaired_sig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$pmw_nonresampled_averaged_median_unpaired_sig != "")] <- paste("+ (", pvalues$pmw_nonresampled_averaged_median_unpaired_sig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$pmw_nonresampled_averaged_median_unpaired_sig != "")], ")", sep = "")
	# pvalues$pmw_nonresampled_averaged_median_unpaired_sig[(pvalues$medianptm < pvalues$mediancontrol) & (pvalues$pmw_nonresampled_averaged_median_unpaired_sig != "")] <- paste("- (", pvalues$pmw_nonresampled_averaged_median_unpaired_sig[(pvalues$medianptm < pvalues$mediancontrol) & (pvalues$pmw_nonresampled_averaged_median_unpaired_sig != "")], ")", sep = "")
	# pvalues$pmw_nonresampled_averaged_median_unpaired_sig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$pmw_nonresampled_averaged_median_unpaired_sig != "")] <- paste("(", pvalues$pmw_nonresampled_averaged_median_unpaired_sig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$pmw_nonresampled_averaged_median_unpaired_sig != "")], ")", sep = "")
	pvalues$pmw_nonresampled_averaged_median_unpaired_sig[(pvalues$'.$directionality.'ptm > pvalues$'.$directionality.'control) & (pvalues$pmw_nonresampled_averaged_median_unpaired_sig != "")] <- paste("(", pvalues$pmw_nonresampled_averaged_median_unpaired_sig[(pvalues$'.$directionality.'ptm > pvalues$'.$directionality.'control) & (pvalues$pmw_nonresampled_averaged_median_unpaired_sig != "")], ")", sep = "")

	pvalues[,"pmw_nonresampled_averaged_mean_unpaired_sig"] <- ""
	try(pvalues[pvalues$pmw_nonresampled_averaged_mean_unpaired < 0.05,]$pmw_nonresampled_averaged_mean_unpaired_sig <- "*", silent=T)
	try(pvalues[pvalues$pmw_nonresampled_averaged_mean_unpaired < 0.01,]$pmw_nonresampled_averaged_mean_unpaired_sig <- "**", silent=T)
	try(pvalues[pvalues$pmw_nonresampled_averaged_mean_unpaired < 0.001,]$pmw_nonresampled_averaged_mean_unpaired_sig <- "***", silent=T)
	# pvalues$pmw_nonresampled_averaged_mean_unpaired_sig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$pmw_nonresampled_averaged_mean_unpaired_sig != "")] <- paste("+ (", pvalues$pmw_nonresampled_averaged_mean_unpaired_sig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$pmw_nonresampled_averaged_mean_unpaired_sig != "")], ")", sep = "")
	# pvalues$pmw_nonresampled_averaged_mean_unpaired_sig[(pvalues$medianptm < pvalues$mediancontrol) & (pvalues$pmw_nonresampled_averaged_mean_unpaired_sig != "")] <- paste("- (", pvalues$pmw_nonresampled_averaged_mean_unpaired_sig[(pvalues$medianptm < pvalues$mediancontrol) & (pvalues$pmw_nonresampled_averaged_mean_unpaired_sig != "")], ")", sep = "")
	# pvalues$pmw_nonresampled_averaged_mean_unpaired_sig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$pmw_nonresampled_averaged_mean_unpaired_sig != "")] <- paste("(", pvalues$pmw_nonresampled_averaged_mean_unpaired_sig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$pmw_nonresampled_averaged_mean_unpaired_sig != "")], ")", sep = "")
	pvalues$pmw_nonresampled_averaged_mean_unpaired_sig[(pvalues$'.$directionality.'ptm > pvalues$'.$directionality.'control) & (pvalues$pmw_nonresampled_averaged_mean_unpaired_sig != "")] <- paste("(", pvalues$pmw_nonresampled_averaged_mean_unpaired_sig[(pvalues$'.$directionality.'ptm > pvalues$'.$directionality.'control) & (pvalues$pmw_nonresampled_averaged_mean_unpaired_sig != "")], ")", sep = "")

	pvalues[,"ptt_nonresampled_sig"] <- ""
	try(pvalues[pvalues$ptt_nonresampled < 0.05,]$ptt_nonresampled_sig <- "*", silent=T)
	try(pvalues[pvalues$ptt_nonresampled < 0.01,]$ptt_nonresampled_sig <- "**", silent=T)
	try(pvalues[pvalues$ptt_nonresampled < 0.001,]$ptt_nonresampled_sig <- "***", silent=T)
	# pvalues$ptt_nonresampled_sig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$ptt_nonresampled_sig != "")] <- paste("+ (", pvalues$ptt_nonresampled_sig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$ptt_nonresampled_sig != "")], ")", sep = "")
	# pvalues$ptt_nonresampled_sig[(pvalues$medianptm < pvalues$mediancontrol) & (pvalues$ptt_nonresampled_sig != "")] <- paste("- (", pvalues$ptt_nonresampled_sig[(pvalues$medianptm < pvalues$mediancontrol) & (pvalues$ptt_nonresampled_sig != "")], ")", sep = "")
	# pvalues$ptt_nonresampled_sig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$ptt_nonresampled_sig != "")] <- paste("(", pvalues$ptt_nonresampled_sig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$ptt_nonresampled_sig != "")], ")", sep = "")
	pvalues$ptt_nonresampled_sig[(pvalues$'.$directionality.'ptm > pvalues$'.$directionality.'control) & (pvalues$ptt_nonresampled_sig != "")] <- paste("(", pvalues$ptt_nonresampled_sig[(pvalues$'.$directionality.'ptm > pvalues$'.$directionality.'control) & (pvalues$ptt_nonresampled_sig != "")], ")", sep = "")

	pvalues[,"pmwsig"] <- ""
	try(pvalues[pvalues$pmw < 0.05,]$pmwsig <- "*", silent=T)
	try(pvalues[pvalues$pmw < 0.01,]$pmwsig <- "**", silent=T)
	try(pvalues[pvalues$pmw < 0.001,]$pmwsig <- "***", silent=T)
	# pvalues$pmwsig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$pmwsig != "")] <- paste("+ (", pvalues$pmwsig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$pmwsig != "")], ")", sep = "")
	# pvalues$pmwsig[(pvalues$medianptm < pvalues$mediancontrol) & (pvalues$pmwsig != "")] <- paste("- (", pvalues$pmwsig[(pvalues$medianptm < pvalues$mediancontrol) & (pvalues$pmwsig != "")], ")", sep = "")
	# pvalues$pmwsig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$pmwsig != "")] <- paste("(", pvalues$pmwsig[(pvalues$medianptm > pvalues$mediancontrol) & (pvalues$pmwsig != "")], ")", sep = "")
	pvalues$pmwsig[(pvalues$'.$directionality.'ptm > pvalues$'.$directionality.'control) & (pvalues$pmwsig != "")] <- paste("(", pvalues$pmwsig[(pvalues$'.$directionality.'ptm > pvalues$'.$directionality.'control) & (pvalues$pmwsig != "")], ")", sep = "")

	pvalues[,"pttsig"] <- ""
	try(pvalues[pvalues$ptt < 0.05,]$pttsig <- "*", silent=T)
	try(pvalues[pvalues$ptt < 0.01,]$pttsig <- "**", silent=T)
	try(pvalues[pvalues$ptt < 0.001,]$pttsig <- "***", silent=T)
	# pvalues$pttsig[(pvalues$meanptm > pvalues$meancontrol) & (pvalues$pttsig != "")] <- paste("+ (", pvalues$pttsig[(pvalues$meanptm > pvalues$meancontrol) & (pvalues$pttsig != "")], ")", sep = "")
	# pvalues$pttsig[(pvalues$meanptm < pvalues$meancontrol) & (pvalues$pttsig != "")] <- paste("- (", pvalues$pttsig[(pvalues$meanptm < pvalues$meancontrol) & (pvalues$pttsig != "")], ")", sep = "")
	pvalues$pttsig[(pvalues$'.$directionality.'ptm > pvalues$'.$directionality.'control) & (pvalues$pttsig != "")] <- paste("(", pvalues$pttsig[(pvalues$'.$directionality.'ptm > pvalues$'.$directionality.'control) & (pvalues$pttsig != "")], ")", sep = "")

	pvalues[,"pmediansig"] <- ""
	try(pvalues[pvalues$pmedian < 0.05,]$pmediansig <- "*", silent=T)
	try(pvalues[pvalues$pmedian < 0.01,]$pmediansig <- "**", silent=T)
	try(pvalues[pvalues$pmedian < 0.001,]$pmediansig <- "***", silent=T)
	pvalues$pmediansig[(pvalues$'.$directionality.'ptm > pvalues$'.$directionality.'control) & (pvalues$pmediansig != "")] <- paste("(", pvalues$pmediansig[(pvalues$'.$directionality.'ptm > pvalues$'.$directionality.'control) & (pvalues$pmediansig != "")], ")", sep = "")

	pvalues[,"pmeansig"] <- ""
	try(pvalues[pvalues$pmean < 0.05,]$pmeansig <- "*", silent=T)
	try(pvalues[pvalues$pmean < 0.01,]$pmeansig <- "**", silent=T)
	try(pvalues[pvalues$pmean < 0.001,]$pmeansig <- "***", silent=T)
	pvalues$pmeansig[(pvalues$'.$directionality.'ptm > pvalues$'.$directionality.'control) & (pvalues$pmeansig != "")] <- paste("(", pvalues$pmeansig[(pvalues$'.$directionality.'ptm > pvalues$'.$directionality.'control) & (pvalues$pmeansig != "")], ")", sep = "")

	pvalues$sig <- pvalues$'.$significance.'

	write.table(pvalues, file="'.$outtable.'", sep="\t", quote=FALSE, row.names=FALSE)
	dput(pvalues, file = "'.$tmppvals.'")

	write.table(evoratedata, file="'.$outevo.'", sep="\t", quote=FALSE, row.names=FALSE)
	# write.table(subset(evoratedata, dis == "Combined"), file="'.$outevo.'", sep="\t", quote=FALSE, row.names=FALSE)
	dput(evoratedata, file = "'.$tmpevo.'")

	# dput(teststatistic, file = "'.$outtest.'")
	write.table(teststatistic, file="'.$outtest.'", sep="\t", quote=FALSE, row.names=FALSE)

	'), 1);
	stoptime();
}
else
{
	# Try to read existing p-value table
	state("Reading existing p-values from '$tmppvals'...");
	die("Error: p-value table '$tmppvals' doesn't exist") if (!-s $tmppvals);
	starttime();
	state(runr('
	# pvalues <- read.table(file="'.$outtable.'", sep="\t", quote="", header=TRUE)
	pvalues <- dget("'.$tmppvals.'")

	pvalues$sig <- pvalues$'.$significance.'

	'), 1);
	stoptime();

	# Try to read existing processed evorate data
	state("Reading processed evorate data from '$tmpevo'...");
	die("Error: Processed evorate table '$tmpevo' doesn't exist") if (!-s $tmpevo);
	starttime();
	state(runr('
	evoratedata <- dget("'.$tmpevo.'")
	'), 1);
	stoptime();
}



# Calculate p-values
state("Calculating disorder percentages for PTM types...");
starttime();
# state(runr('
# # Test
# subset(pvalues, ptm=="K-ac" & dis=="Disordered sites")$nptm
# ', 1), 1);
# state(runr('
# # Test
# subset(pvalues, ptm=="K-ac" & dis=="Combined")$nptm
# ', 1), 1);
# state(runr('str(evoratedata)', 1), 1);
# state(runr('head(evoratedata)', 1), 1);
# state(runr('str(pvalues)', 1), 1);
# state(runr('head(pvalues)', 1), 1);


state(runr('
percentages <- data.frame()
for (x in unique(evoratedata$ptm)) {
	for (y in sort(unique(evoratedata$dis))) {
		# print(paste("x (ptm)", x))
		# print(paste("y (dis)", y))
		pctthis <- subset(pvalues, ptm==x & dis==y)$nptm
		# print(paste("pctthis", pctthis))
		# pctall <- sum(subset(pvalues, ptm==x & dis!="Combined")$nptm)
		# print(paste("pctall", pctall))
		# print(paste(round((pctthis/pctall) * 100, 1)), "%", sep="")
		# percentages <- rbind(percentages, data.frame(ptm=x, dis=y, pct=paste(round((pctthis/pctall) * 100, 1), "%", sep="")))
		# percentages <- rbind(percentages, data.frame(ptm=x, dis=y, pct=paste(round((pctthis/pctall) * 100), "%", sep="")))
		# percentages <- rbind(percentages, data.frame(ptm=x, dis=y, pct=paste("n = ", pctthis, sep="")))

		percentages <- rbind(percentages, data.frame(ptm=x, dis=y, pct=paste("n=", pctthis, sep="")))

		# if (pctthis >= mincases) {
		# 	percentages <- rbind(percentages, data.frame(ptm=x, dis=y, pct=paste("n=", pctthis, sep="")))
		# }
		# else {
		# 	percentages <- rbind(percentages, data.frame(ptm=x, dis=y, pct=paste0("n<", mincases)))
		# }
	}
}
'), 1);

# state(runr('str(percentages)', 1), 1);
# state(runr('head(percentages)', 1), 1);
# state(runr('
# # Bigger test
# percentages
# ', 1), 1);
stoptime();

# state("Removing disorder/core/surface types with too few cases...");
# starttime();
# state(runr('
# for (y in evoratedata$dis) {
# 	if (nrows(evoratedata[evoratedata$dis==y],) == 0) {
# 		
# 	}
# }
# 
# '), 1);


$scalebox = '';
$scaledens = '';
$pvalpos = '-Inf';
$ymax = 4;
if (switch('nolog') or ($evorate =~ /^norm_/))
{
	# Linear
	# Normalised (i.e. z-scores)
	
	# $scalebox = '+ scale_y_continuous(trans = "identity", breaks = trans_breaks("identity", function(x) x'.$desiredbreaks.'), labels = trans_format("identity", math_format(.x)))';
	# $scaledens = '+ scale_x_continuous(trans = "identity", breaks = trans_breaks("identity", function(x) x'.$desiredbreaks.'), labels = trans_format("identity", math_format(.x)))';
	# $scalebox = '+ scale_y_continuous(breaks = pretty_breaks(20))';
	# $scaledens = '+ scale_x_continuous(breaks = pretty_breaks(20))';

	if ($evorate =~ /^norm_capra/)
	{
		# $scalebox = '+ scale_y_continuous(breaks = pretty_breaks(10))';
		# $scaledens = '+ scale_x_continuous(breaks = pretty_breaks(10))';
		# $ymax = 1.3;
	}
	if ($evorate =~ /^norm_lichtarge/)
	{
		# $scalebox = '+ scale_y_continuous(breaks = pretty_breaks(20))';
		# $scaledens = '+ scale_x_continuous(breaks = pretty_breaks(20))';
	}
	if ($evorate =~ /^norm_rate4site/)
	{
		# $s = '_continuous(breaks = pretty_breaks(20), limits=c(min(-2, min(evoratedata$rate)), max(evoratedata$rate)))';
		# $scalebox = '+ scale_y'.$s;
		# $scaledens = '+ scale_x'.$s;
		# $ymax = 2;
	}

	$scalebox = '';
	$scaledens = '';

	$pvalpos = '-Inf';
}
else
{
	# Log
	# Raw (i.e. raw scores)

	# $scalebox = '+ scale_y_log10()';
	# $scalebox = '+ scale_y_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x'.$scalebox.'), labels = trans_format("log10", math_format(10^.x)))';

	# $scaledens = '+ scale_x_log10()';
	# $scaledens = '+ scale_x_continuous(trans = "log2", breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x)))';
	# $scaledens = '+ scale_x_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x'.$scalebox.'), labels = trans_format("log10", math_format(10^.x)))';

	if ($evorate =~ /^capra/)
	{
		# $scalebox =  '+ scale_y_log10()';
		# $scaledens = '+ scale_x_log10()';

		# $s = '_continuous(trans = "log2", breaks = trans_breaks("log2", function(x) 2^x))';
		# $s = '_log10(breaks = pretty_breaks(5))';
		$s = '_log10(expand = c(0, 0), breaks = pretty_breaks(5))';
		
		$scalebox = '+ scale_y'.$s;
		$scaledens = '+ scale_x'.$s;

		# # $ymax = 6;
		# $ymax = 4;
	}
	if ($evorate =~ /^lichtarge/)
	{
		# $scalebox =  '+ scale_y_log10()';
		# $scaledens = '+ scale_x_log10()';

		# $s = '_log10(breaks = pretty_breaks(5))';
		# $s = '_continuous(trans = "log2", breaks = trans_breaks("log2", function(x) 2^x))';
		$s = '_continuous(expand = c(0, 0), trans = "log2", breaks = trans_breaks("log2", function(x) 2^x))';
		
		$scalebox = '+ scale_y'.$s;
		$scaledens = '+ scale_x'.$s;

		# $ymax = 1.2;
	}
	if ($evorate =~ /^rate4site/)
	{
		# $s = '_log10(breaks = pretty_breaks(4))';
		# $s = '_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x, n=10))';
		# $s = '_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x, n=4), labels = trans_format("log10", math_format(10^.x)))';
		# $s = '_continuous(expand = c(0, 0), trans = "log10", breaks = trans_breaks("log10", function(x) 10^x, n=4), labels = trans_format("log10", math_format(10^.x)))';
		# Fixed limits:
		$s = '_continuous(limits = c(2e-3, 9e1), expand = c(0, 0), trans = "log10", breaks = trans_breaks("log10", function(x) 10^x, n=4), labels = trans_format("log10", math_format(10^.x)))';
		$scalebox = '+ scale_y'.$s;
		$scaledens = '+ scale_x'.$s;

		# $scaledens += '+ scale_y_continuous(expand = expansion(mult = c(0, 0.1)))';

		# $ymax = 1.5;
	}

	$pvalpos = '0';
}


# # Boxplots
# if (switch('nolog') or ($evorate =~ /^norm_/))
# {
# 	# Linear
# 	state("Drawing linear (non-logged) boxplots to '$outbox'...");
# 	starttime();
# 	# Linear scale
# 	state(runr('
# 	# Draw boxplots
# 	# p <- ggplot(evoratedata, aes(type, rate)) + theme_bw() '.$scalebox.' + scale_x_discrete(expand=c(0.4, 0)) + geom_boxplot(aes(color = type, fill = type, alpha=0.7), notch = TRUE, notchwidth = 0.5, outlier.size = 1, outlier.shape = NA) + scale_colour_manual(values = c('.$colours.')) +  scale_fill_manual(values = c('.$colours.')) + facet_grid(ptm ~ dis) + xlab("") + ylab("'.$xlabel.'") + labs(title = "'."$evorate $predictor all".'") + theme(legend.position = "none", strip.background = element_rect(fill = NA, colour = NA), axis.title.y = element_text(angle=90), strip.text.y = element_text(angle=0)) + geom_text(data=pvalues, aes(x=Inf, y=-Inf, label=paste("\n", " ", pmwsig, sep=""), hjust=0, vjust=0.5), size = 3) + geom_text(data=percentages, aes(x=Inf, y=Inf, label=paste("\n", pct, " ", sep=""), hjust=1, vjust=0.5), size = 3)
# 	p <- ggplot(evoratedata, aes(type, rate)) + theme_bw() '.$scalebox.' + scale_x_discrete(expand=c(0.4, 0)) + geom_boxplot(aes(color = type, fill = type, alpha=0.7), notch = TRUE, notchwidth = 0.5, outlier.size = 1, outlier.shape = NA) + scale_colour_manual(values = c('.$colours.')) +  scale_fill_manual(values = c('.$colours.')) + facet_grid(factor(ptm, levels = pvalues %>% filter(dis == "All sites") %>% arrange(desc(nptm)) %>% pull(ptm)) ~ dis) + xlab("") + ylab("'.$xlabel.'") + labs(title = "'."$evorate $predictor all".'") + theme(legend.position = "none", strip.background = element_rect(fill = NA, colour = NA), axis.title.y = element_text(angle=90), strip.text.y = element_text(angle=0)) + geom_text(data=pvalues, aes(x=Inf, y=-Inf, label=paste("\n", " ", pmwsig, sep=""), hjust=0, vjust=0.5), size = 3) + geom_text(data=percentages, aes(x=Inf, y=Inf, label=paste("\n", pct, " ", sep=""), hjust=1, vjust=0.5), size = 3)
# 	# p <- p + geom_hline(data=pvalues, aes(yintercept='.$lineplace.'ptm), color = "#4f81bd", linetype = "'.$linetype.'"'.$linesize.')
# 	# p <- p + geom_hline(data=pvalues, aes(yintercept='.$lineplace.'control), color = "#bec1c0", linetype = "'.$linetype.'"'.$linesize.')
# 
# 	# compute lower and upper whiskers
# 	ylim1 = boxplot.stats(evoratedata$rate)$stats[c(1, 5)]
# 	ylim1
# 
# 	# Make boundaries symmetric
# 	ylim1 <- c(-max(abs(ylim1)), max(abs(ylim1)))
# 
# 	# scale y limits based on ylim1
# 	p <- p + coord_flip(ylim = ylim1*1.05)
# 
# 	# ggsave(p, file = "'.$outbox.'", width=234.48, height=154.28, units="mm")
# 	ggsave(p, file = "'.$outbox.'", width=183, height=247, units="mm")
# 	', 1), 1);
# 	# geom_text(data=pvalues, aes(x=-Inf, y=Inf, label=paste("\n", " ", pmwsig, sep=""), hjust=0, vjust=0.5), size = 3)
# }
# else
# {
# 	# Logged
# 	state("Drawing logged boxplots to '$outbox'...");
# 	starttime();
# 	# Log2
# 	# $scale = '+ scale_y_continuous(trans = "log2", breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x)))';
# 	# Log10
# 	# $scale = '+ scale_y_log10()';
# 	state(runr('
# 	# Draw boxplots
# 	# p <- ggplot(evoratedata, aes(type, rate)) + scale_x_discrete(expand=c(0.4, 0)) + theme_bw() '.$scalebox.' + geom_boxplot(aes(color = type, fill = type, alpha=0.7), notch = TRUE, notchwidth = 0.5, outlier.size = 1, outlier.shape = NA) + scale_colour_manual(values = c('.$colours.')) +  scale_fill_manual(values = c('.$colours.')) + facet_grid(ptm ~ dis) + xlab("") + ylab("'.$xlabel.'") + labs(title = "'."$evorate $predictor all".'") + theme(legend.position = "none", strip.background = element_rect(fill = NA, colour = NA), axis.title.y = element_text(angle=90), strip.text.y = element_text(angle=0)) + geom_text(data=pvalues, aes(x=Inf, y=0, label=paste("\n", " ", pmwsig, sep=""), hjust=0, vjust=0.5), size = 3) + geom_text(data=percentages, aes(x=Inf, y=Inf, label=paste("\n", pct, " ", sep=""), hjust=1, vjust=0.5), size = 3) + coord_flip()
# 	p <- ggplot(evoratedata, aes(type, rate)) + scale_x_discrete(expand=c(0.4, 0)) + theme_bw() '.$scalebox.' + geom_boxplot(aes(color = type, fill = type, alpha=0.7), notch = TRUE, notchwidth = 0.5, outlier.size = 1, outlier.shape = NA) + scale_colour_manual(values = c('.$colours.')) +  scale_fill_manual(values = c('.$colours.')) + facet_grid(factor(ptm, levels = pvalues %>% filter(dis == "All sites") %>% arrange(desc(nptm)) %>% pull(ptm)) ~ dis) + xlab("") + ylab("'.$xlabel.'") + labs(title = "'."$evorate $predictor all".'") + theme(legend.position = "none", strip.background = element_rect(fill = NA, colour = NA), axis.title.y = element_text(angle=90), strip.text.y = element_text(angle=0)) + geom_text(data=pvalues, aes(x=Inf, y=0, label=paste("\n", " ", pmwsig, sep=""), hjust=0, vjust=0.5), size = 3) + geom_text(data=percentages, aes(x=Inf, y=Inf, label=paste("\n", pct, " ", sep=""), hjust=1, vjust=0.5), size = 3) + coord_flip()
# 	# p <- p + geom_hline(data=pvalues, aes(yintercept='.$lineplace.'ptm), color = "#4f81bd", linetype = "'.$linetype.'"'.$linesize.')
# 	# p <- p + geom_hline(data=pvalues, aes(yintercept='.$lineplace.'control), color = "#bec1c0", linetype = "'.$linetype.'"'.$linesize.')
# 	# ggsave(p, file = "'.$outbox.'", width=234.48, height=154.28, units="mm")
# 	ggsave(p, file = "'.$outbox.'", width=183, height=247, units="mm")
# 	', 1), 1);
# 	# geom_text(data=pvalues, aes(x=0, y=Inf, label=paste("\n", " ", pttsig, sep=""), hjust=0, vjust=0.5), size = 3)
# }
# stoptime();


# Density plots
if (switch('nolog') or ($evorate =~ /^norm_/))
{
	# Linear
	state("Drawing linear (non-logged) density plots to '$outdens'...");

# 	state(runr('
# 
# 	# compute lower and upper whiskers
# 	xlim1 <- boxplot.stats(evoratedata$rate)$stats[c(1, 5)]
# 	xlim1
# 	# Make boundaries symmetric
# 	xlim1 <- c(-max(abs(xlim1)), max(abs(xlim1)))
# 
# 	', 1), 1);
}
else
{
	# Logged, including "All sites"
	state("Drawing logged density plots to '$outdens'...");

# 	state(runr('
# 
# 	xlim1 <- NULL
# 
# 	', 1), 1);
}

# Store PTM order, sorted by total number of sites ("All sites")
state(runr('
	pvalues %>% filter(dis == "All sites") %>% arrange(desc(nptm)) %>% pull(ptm) -> ptm_order
', 1), 1);

if (!switch('allsites'))
{
	state(runr('

	# Filter out "All sites"
	evoratedata %<>% filter(dis != "All sites")
	pvalues %<>% filter(dis != "All sites")
	percentages %<>% filter(dis != "All sites")

	# Reset factor levels
	evoratedata$dis <- as.factor(evoratedata$dis)
	pvalues$dis <- as.factor(pvalues$dis)
	percentages$dis <- as.factor(percentages$dis)

	', 1), 1);
}

starttime();

state(runr('

# evoratedata$dis <- as.character(evoratedata$dis)
print(head(evoratedata))
print(pvalues)
print(summary(evoratedata))
# print(levels(evoratedata$dis))
print(evoratedata$dis %>% unique %>% sort)

# Draw density plots
p <- ggplot(evoratedata, aes(rate)) +
theme_bw() '.$scaledens.' +
geom_density(aes(color = type, fill = type, alpha = 0.2), adjust='.$densityadjust.') +
geom_vline(data=pvalues, aes(xintercept='.$lineplace.'ptm), color = "#4f81bd", linetype = "'.$linetype.'"'.$linesize.') +
geom_vline(data=pvalues, aes(xintercept='.$lineplace.'control), color = "#bec1c0", linetype = "'.$linetype.'"'.$linesize.') +
geom_text(data=pvalues, aes(x='.$pvalpos.', y=Inf, label=paste("\n", " ", sig, sep=""), hjust=0, vjust=0.5), size = 3) +
geom_text(data=percentages, aes(x=Inf, y=Inf, label=paste("\n", pct, " ", sep=""), hjust=1, vjust=0.5), size = 2.5) +
scale_colour_manual(values = c('.$colours.')) +
scale_fill_manual(values = c('.$colours.')) +
# scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +			# make the plot a little higher
facet_grid(factor(ptm, levels = ptm_order) ~ '.$dis_order.') +
xlab("'.$xlabel.'") +
ylab("Probability density") +
labs(title = "'."$evorate $predictor all".'") +
theme_minimal() +
# theme(legend.position = "none", strip.background = element_rect(fill = NA, colour = NA), axis.title.y = element_text(angle=90), strip.text.y = element_text(angle=0))
theme(legend.position = "none", strip.background = element_rect(fill = NA, colour = NA), axis.text.y = element_text(size=7), axis.title.y = element_text(angle=90), strip.text.y = element_text(angle=0))

# # y limit
# ylim1 <- ggplot_build(p)
# # ylim1 <- ylim1$panel$ranges[[1]]$y.range
# ylim1 <- layer_scales(p)$y$range$range
# ylim2 <- ylim1[[2]]
# ylim1 <- ylim1[[1]]
# print(paste("Y LIMITS", ylim1, ylim2))
# old2 = ylim2
# ylim2 = min(c(ylim2, '.$ymax.'))
# ylim1 = (ylim2 / old2) * ylim1
# ylim2 = ylim2 * 1.2		# make it a little higher
# p <- p + coord_cartesian(xlim = xlim1*1.05, ylim = c(ylim1, ylim2))

# Save figure
# Note: ggsave dimensions in millimeters arent fully accurate. The actual dimensions e.g. Illustrator shows usually differ very slightly. The actual unit used by the PDF writer seems to be inches-based, introducing rounding errors.
# Note: geom_text with position set to Inf or -Inf, in combination with the log2 scale, is what causes the "Warning: Transformation introduced infinite values in continuous x-axis" warnings, which can safely be ignored.

# ggsave(p, file = "'.$outdens.'", width=234.48, height=154.28, units="mm")
# ggsave(p, file = "'.$outdens.'", width=180, height=170, units="mm")		# "Figures are best prepared at a width of 90 mm (single column) and 180 mm (double column) with a maximum height of 170mm. At this size, the font size should be 5-7pt." (https://www.nature.com/nature/for-authors/initial-submission)
# ggsave(p, file = "'.$outdens.'", width=183, height=247, units="mm")		# Absolute max Nature dimensions 183 mm x 247 mm ("For guidance, Natures standard figure sizes are 89 mm wide (single column) and 183 mm wide (double column). The full depth of a Nature page is 247 mm. Figures can also be a column-and-a-half where necessary (120–136 mm).", https://www.nature.com/nature/for-authors/final-submission)
# ggsave(p, file = "'.$outdens.'", width=183, height=length(unique(evoratedata$ptm)) * 14.53, units="mm")				# 17 PTMs = 14.53 * 17 = 247.01 mm (i.e. 183 mm x 247 mm, the Nature maximum)
# ggsave(p, file = "'.$outdens.'", width=183, height=24.385 + (length(unique(evoratedata$ptm)) * 13.095), units="mm")		# 17 PTMs = 24.385 + 17 * 13.095 = 247 mm
ggsave(p, file = "'.$outdens.'", width = 61 * (evoratedata$dis %>% unique %>% length), height = 24.385 + (length(unique(evoratedata$ptm)) * 13.095), units="mm")		# 17 PTMs = 24.385 + 17 * 13.095 = 247 mm, and 61 mm width per column

'), 1);

stoptime();

showmeall(1);

# state("Wrote boxplots to '$outbox'", 1);
state("Wrote density plots to '$outdens'", 1);
state("Wrote p-values to '$outtable'", 1);

done();
