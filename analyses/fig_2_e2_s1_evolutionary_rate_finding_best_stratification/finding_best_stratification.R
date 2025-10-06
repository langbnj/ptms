blang_init()

library("coin")


# "Figures are best prepared at a width of 90 mm (single column) and 180 mm (double column) with a maximum height of 170mm. At this size, the font size should be 5-7pt." (https://www.nature.com/nature/for-authors/initial-submission)
# Absolute max Nature dimensions 183 mm x 247 mm ("For guidance, Natures standard figure sizes are 89 mm wide (single column) and 183 mm wide (double column). The full depth of a Nature page is 247 mm. Figures can also be a column-and-a-half where necessary (120–136 mm).", https://www.nature.com/nature/for-authors/final-submission)

# To get source files:
# cd ~/pipeline/evolutionary_rate_analysis/tmp; mkdir -p ../tmp_copy; \cp -vpf tmp-dataframe-all-human-Ochoa-lichtarge_einsi_tree_1para-AlphaFold-coresurf.txt tmp-dataframe-all-human-Ochoa-rate4site_einsi_tree_1para-AlphaFold-coresurf.txt tmp-dataframe-all-human-Ochoa-rate4site_einsi_tree_1para-AlphaFold-nosurf.txt.rds tmp-dataframe-all-human-PhosphoSitePlus-lichtarge_einsi_tree_1para-AlphaFold-coresurf.txt tmp-dataframe-all-human-PhosphoSitePlus-rate4site_einsi_tree_1para-AlphaFold-coresurf.txt tmp-dataframe-all-human-UniProt-lichtarge_einsi_tree_1para-AlphaFold-coresurf.txt tmp-dataframe-all-human-UniProt-rate4site_einsi_tree_1para-AlphaFold-coresurf.txt tmp-dataframe-all-human-all-capra0_einsi_tree_1para-AlphaFold-coresurf.txt tmp-dataframe-all-human-all-lichtarge_einsi-AlphaFold-coresurf.txt tmp-dataframe-all-human-all-lichtarge_einsi-AlphaFold-nosurf.txt tmp-dataframe-all-human-all-lichtarge_einsi_tree-AlphaFold-coresurf.txt tmp-dataframe-all-human-all-lichtarge_einsi_tree-AlphaFold-nosurf.txt tmp-dataframe-all-human-all-lichtarge_einsi_tree_1para-AlphaFold-coresurf.txt tmp-dataframe-all-human-all-lichtarge_einsi_tree_1para-AlphaFold-nosurf.txt tmp-dataframe-all-human-all-lichtarge_ginsi-AlphaFold-coresurf.txt tmp-dataframe-all-human-all-lichtarge_ginsi-AlphaFold-nosurf.txt tmp-dataframe-all-human-all-lichtarge_ginsi_tree-AlphaFold-coresurf.txt tmp-dataframe-all-human-all-lichtarge_ginsi_tree-AlphaFold-nosurf.txt tmp-dataframe-all-human-all-lichtarge_linsi-AlphaFold-coresurf.txt tmp-dataframe-all-human-all-lichtarge_linsi-AlphaFold-nosurf.txt tmp-dataframe-all-human-all-lichtarge_linsi_tree-AlphaFold-coresurf.txt tmp-dataframe-all-human-all-lichtarge_linsi_tree-AlphaFold-nosurf.txt tmp-dataframe-all-human-all-rate4site_einsi-AlphaFold-coresurf.txt tmp-dataframe-all-human-all-rate4site_einsi_tree-AlphaFold-coresurf.txt tmp-dataframe-all-human-all-rate4site_einsi_tree_1para-AlphaFold-coresurf-accs.txt tmp-dataframe-all-human-all-rate4site_einsi_tree_1para-AlphaFold-coresurf.accsites.txt tmp-dataframe-all-human-all-rate4site_einsi_tree_1para-AlphaFold-coresurf.txt tmp-dataframe-all-human-all-rate4site_ginsi-AlphaFold-coresurf.txt tmp-dataframe-all-human-all-rate4site_ginsi_tree-AlphaFold-coresurf.txt tmp-dataframe-all-human-all-rate4site_linsi-AlphaFold-coresurf.txt tmp-dataframe-all-human-all-rate4site_linsi_tree-AlphaFold-coresurf.txt tmp-dataframe-all-human-dbPTM-lichtarge_einsi_tree_1para-AlphaFold-coresurf.txt tmp-dataframe-all-human-dbPTM-rate4site_einsi_tree_1para-AlphaFold-coresurf.txt tmp-dataframe-all-human-all-norm_capra0_einsi_tree_1para-AlphaFold-coresurf.txt tmp-dataframe-all-human-all-norm_lichtarge_einsi_tree_1para-AlphaFold-coresurf.txt tmp-dataframe-all-human-all-norm_rate4site_einsi_tree_1para-AlphaFold-coresurf.txt ../tmp_copy/
# Then copy all tmp-dataframe-... files in tmp_copy to this directory here.

# Function to assign significance label based on alternative="two.sided" and alternative="greater" (for directionality) p-values
# e.g. "***" (less), "(***)" (greater)
significance_label <- function(pvalue_twosided, pvalue_greater) {
  mysig <- ""
  if (pvalue_twosided < 0.05) { mysig <- "*" }
  if (pvalue_twosided < 0.01) { mysig <- "**" }
  if (pvalue_twosided < 0.001) { mysig <- "***" }
  if ((pvalue_greater < 0.05) & (mysig != "")) {
    mysig <- f("({mysig})")
  }
  return(mysig)
}

# Two-tailed, with "()" for inverse median:
filter <- "filter_none"
get_pvalues <- function(q, source, evorate, pred = "AlphaFold", filter = "filter_none", strat_minsize = 2) {
  tmppred <- ""
  if (pred != "AlphaFold") {
    tmppred <- f("-{pred}")
  }
  if (filter %in% c("filter_no_control_residues_strat", "filter_too_few_control_residues_strat", "filter_strat", "filter_strat_same_proteins")) {
    tmprds <- f("tmp-pvalues-{source}-{evorate}{tmppred}-{filter}-strat_minsize{strat_minsize}.rds")
  } else {
    tmprds <- f("tmp-pvalues-{source}-{evorate}{tmppred}-{filter}.rds")
  }
  if (file.exists(tmprds)) {
    # Load rds
    print(f("Reading pvalues from '{tmprds}'..."))
    pvalues <- read_rds(tmprds)
    # Remove "Surface ", turning "Surface Structured" and "Surface Disordered" into "Structured" and "Disordered" as in the other plots
    # pvalues %<>% mutate(dis = str_replace(dis, "^Surface ", ""))
    # Rename levels instead, since this is a factor
    levels(pvalues$dis) <- list("Buried"="Buried", "Structured"="Surface Structured", "Disordered"="Surface Disordered")
  } else {
    pvalues <- tibble()
    # myptm <- "C-nit"
    # mydis <- "dissurf"
    # myptm <- "Y-p"
    # mydis <- "dissurf"
    # myptm <- "S-p"
    # mydis <- "strcore"
    # mydis <- "Buried"
    myptm <- "C-glt"
    mydis <- "Disordered"
    # for (mydis in sort(unique(q$dis))) {
    for (mydis in q %>% pull(dis) %>% unique %>% sort) {
      # for (myptm in sort(unique(q$ptm))) {
      for (myptm in q %>% filter(dis == mydis) %>% pull(ptm) %>% unique %>% sort) {
        
        qtmp <- q %>% filter(ptm == myptm & dis == mydis)
        # print("QTMP=")
        # print(qtmp)
        
        # if (qtmp %>% nrow == 0) {
        # Make sure we have both PTM and Control sites (e.g. for C-glt with pLDDT70 filtering, there weren't any for dissurf)
        if (qtmp %>% pull(type) %>% unique %>% length != 2) {
          next
        }
        
        # Get counts before filtering out accs without controls
        nptm_pre <- qtmp %>% filter(ptm == myptm & dis == mydis & type == "P") %>% unique %>% nrow
        # ncontrol should not be affected
        ncontrol_pre <- qtmp %>% filter(ptm == myptm & dis == mydis & type == "C") %>% unique %>% nrow
        
        if (filter %in% c("filter_no_control_residues", "filter_no_control_residues_strat", "filter_both")) {
          # Limit to proteins that have both PTM and control residues
          # i.e. skip any proteins that do not have control residues
          # qtmp <- subset(qtmp, acc %in% intersect(subset(qtmp, type == "P")$acc, subset(qtmp, type == "C")$acc))
          qtmp$acc %>% unique %>% length
          qtmp %>% filter(type == "P") %>% select(acc, site) %>% unique %>% nrow
          sort(subset(qtmp, type == "P")$acc)
          sort(subset(qtmp, type == "C")$acc)
          sort(intersect(subset(qtmp, type == "P")$acc, subset(qtmp, type == "C")$acc))
          qtmp <- subset(qtmp, acc %in% intersect(subset(qtmp, type == "P")$acc, subset(qtmp, type == "C")$acc))
          qtmp$acc %>% unique %>% length
          qtmp %>% filter(type == "P") %>% select(acc, site) %>% unique %>% nrow
        }
        
        if (filter %in% c("filter_too_few_control_residues", "filter_too_few_control_residues_strat", "filter_both")) {
          # Filter out accs that do not have at least as many control residues as PTM residues for this ptm|dis
          # qtmp
          # qtmp %>% select(acc, site, type, rate)
          # qtmp %>% select(acc, site, type, rate) %>% unique
          # >> OK, acc|site|type|rate are the only informative columns (because qtmp is already filtered for ptm|dis), and there are no duplicate rows
          # qtmp %>% group_by(acc, type) %>% tally %>% pivot_wider(names_from = type, values_from = n)
          # qtmp %>% group_by(acc, type) %>% tally %>% pivot_wider(names_from = type, values_from = n) %>% filter(is.na(P))
          # qtmp %>% group_by(acc, type) %>% tally %>% pivot_wider(names_from = type, values_from = n) %>% filter(is.na(C))
          qtmp %>% group_by(acc, type) %>% tally %>% pivot_wider(names_from = type, values_from = n) %>% group_by(C >= P) %>% tally
          # >> For Y-p|dissurf, over half of accs have too few control residues! (Out of 5776 accs, 1941 have none and 1061 have too few vs. 2774 that have enough, i.e. only 48% of accs would be considered if I were to use this filter!
          # >> For S-p|dissurf as well, almost half of accs would be lost.
          # qtmp %>% group_by(acc, type) %>% tally %>% pivot_wider(names_from = type, values_from = n) %>% filter(C >= P) %>% pull(acc)
          qtmp <- subset(qtmp, acc %in% (qtmp %>% group_by(acc, type) %>% tally %>% pivot_wider(names_from = type, values_from = n) %>% filter(C >= P) %>% pull(acc)))
          qtmp$acc %>% unique %>% length
          qtmp %>% filter(type == "P") %>% select(acc, site) %>% unique %>% nrow
          qtmp %>% filter(type == "P") %>% pull(acc) %>% unique %>% length
        }
        
        # Get median rates
        medianptm <- median(subset(qtmp, type =="P")$rate)
        mediancontrol <- median(subset(qtmp, type =="C")$rate)
        
        # Get mean rates
        meanptm <- mean(subset(qtmp, type =="P")$rate)
        meancontrol <- mean(subset(qtmp, type =="C")$rate)

        # Get median rates (stratified by acc)
        if (filter %in% c("filter_no_control_residues_strat", "filter_too_few_control_residues_strat", "filter_strat", "filter_strat_same_proteins")) {
          # qtmp %>% filter(type =="P") %>% group_by(acc) %>% summarise(median_rate = median(rate))
          # median(qtmp %>% filter(type =="P") %>% group_by(acc) %>% summarise(median_rate = median(rate)) %>% pull(median_rate))
          medianptmstrat <- median(qtmp %>% filter(type =="P") %>% group_by(acc) %>% summarise(median_rate = median(rate)) %>% pull(median_rate))
          mediancontrolstrat <- median(qtmp %>% filter(type =="C") %>% group_by(acc) %>% summarise(median_rate = median(rate)) %>% pull(median_rate))
          
          # Get mean rates (stratified by acc)
          meanptmstrat <- mean(qtmp %>% filter(type =="P") %>% group_by(acc) %>% summarise(mean_rate = mean(rate)) %>% pull(mean_rate))
          meancontrolstrat <- mean(qtmp %>% filter(type =="C") %>% group_by(acc) %>% summarise(mean_rate = mean(rate)) %>% pull(mean_rate))
        }
        
        # Calculate p-value
        # # One-tailed (PTM sites < Control sites) using coin wilcox_test
        # pval <- pvalue(wilcox_test(rate ~ type, alternative="less", data=qtmp))
        # # One-tailed (PTM sites < Control sites) using built-in wilcox.test
        # pval <- wilcox.test(rate ~ type, alternative="less", data=qtmp)$p.value
        # Two-tailed p-value using built-in wilcox.test (default)
        pmw <- wilcox.test(rate ~ type, data=qtmp)$p.value
        pmwless <- wilcox.test(rate ~ type, alternative="less", data=qtmp)$p.value
        pmwgreater <- wilcox.test(rate ~ type, alternative="greater", data=qtmp)$p.value
        
        # Two-tailed p-value using coin wilcox_test, stratified by acc
        # qtmp
        # qtmp %>% select(acc, site, type, rate) %>% unique
        # Filter out accs with fewer than 2 observations
        # qtmp %>% group_by(acc) %>% tally %>% filter(n < 2) %>% pull(acc)
        # subset(qtmp, !(acc %in% (qtmp %>% group_by(acc) %>% tally %>% filter(n < 2) %>% pull(acc))))
        # qtmp <- subset(qtmp, !(acc %in% (qtmp %>% group_by(acc) %>% tally %>% filter(n < 2) %>% pull(acc))))
        # wilcox_test(rate ~ type | acc, data=qtmp)
        # wilcox_test(rate ~ type, data=qtmp) %>% pvalue %>% as.numeric
        # wilcox_test(rate ~ type | acc, data=subset(qtmp, !(acc %in% (qtmp %>% group_by(acc) %>% tally %>% filter(n < 2) %>% pull(acc)))) %>% mutate(acc = as.factor(acc))) %>% pvalue %>% as.numeric
        if (filter %in% c("filter_no_control_residues_strat", "filter_too_few_control_residues_strat", "filter_strat", "filter_strat_same_proteins")) {
          
          # Get fraction of PTM sites retained for current strat_minsize:
          strat_minsize_retained <- qtmp %>% group_by(acc) %>% tally %>% filter(n >= strat_minsize) %>% pull(acc) %>% unique %>% length
          strat_minsize_total <-    qtmp %>% pull(acc) %>% unique %>% length
          strat_minsize_fraction_retained <- strat_minsize_retained / strat_minsize_total
          
          if (strat_minsize_retained > 0) {
            pmwstrat <- wilcox_test(rate ~ type | acc, data=subset(qtmp, !(acc %in% (qtmp %>% group_by(acc) %>% tally %>% filter(n < strat_minsize) %>% pull(acc)))) %>% mutate(acc = as.factor(acc))) %>% pvalue %>% as.numeric
            pmwstratless <- wilcox_test(rate ~ type | acc, alternative="less", data=subset(qtmp, !(acc %in% (qtmp %>% group_by(acc) %>% tally %>% filter(n < strat_minsize) %>% pull(acc)))) %>% mutate(acc = as.factor(acc))) %>% pvalue %>% as.numeric
            pmwgreaterstrat <- wilcox_test(rate ~ type | acc, alternative="greater", data=subset(qtmp, !(acc %in% (qtmp %>% group_by(acc) %>% tally %>% filter(n < strat_minsize) %>% pull(acc)))) %>% mutate(acc = as.factor(acc))) %>% pvalue %>% as.numeric
          } else {
            # Insignificant
            pmwstrat <- 1
            pmwstratless <- 1
            pmwgreaterstrat <- 1
          }
        }
        # Note: Using sigstrat means losing the following number of accs from the wilcox_test analysis (since their blocks/strata would be of size 1):
        # qtmp %>% group_by(acc) %>% tally %>% filter(n < 2) %>% pull(acc) %>% unique %>% length
        # qtmp %>% pull(acc) %>% unique %>% length
        # >> 489 out of 9974 (4.9%) for S-p Buried
        # >> That's even with n < 2 (strat_minsize = 2).
        
        # Assign *** for pmw using medians
        sig <- ""
        if (pmw < 0.05) { sig <- "*" }
        if (pmw < 0.01) { sig <- "**" }
        if (pmw < 0.001) { sig <- "***" }
        
        # Directionality
        # # Using median
        # if ((medianptm > mediancontrol) & (sig != "")) {
        #   sig <- f("({sig})")
        # }
        # Using one-tailed test
        # pmwless
        # pmwgreater
        # if ((pmwless > 0.95) & (pmwgreater > 0.95)) { stop(f("Error: Both pmwless and pmwgreater are close to 1: {myptm} {mydis} medianptm {medianptm} mediancontrol {mediancontrol} pmwless {pmwless} pmwgreater {pmwgreater}")) }
        # if ((pmwless < 0.95) & (pmwgreater < 0.95)) { stop(f("Error: Neither pmwless nor pmwgreater is close to 1: {myptm} {mydis} medianptm {medianptm} mediancontrol {mediancontrol} pmwless {pmwless} pmwgreater {pmwgreater}")) }
        # if ((pmwless > 0.95) & (pmwgreater < 0.95)) {
        if ((pmwgreater < 0.05) & (sig != "")) {
          sig <- f("({sig})")
        }
        
        # Assign *** for sigstrat using medianstrats
        if (filter %in% c("filter_no_control_residues_strat", "filter_too_few_control_residues_strat", "filter_strat", "filter_strat_same_proteins")) {
          sigstrat <- ""
          if (pmwstrat < 0.05) { sigstrat <- "*" }
          if (pmwstrat < 0.01) { sigstrat <- "**" }
          if (pmwstrat < 0.001) { sigstrat <- "***" }

          # Directionality
          # Using median
          # # Using median
          # if ((medianptmstrat > mediancontrolstrat) & (sigstrat != "")) {
          #   sigstrat <- f("({sigstrat})")
          # }
          # Using one-tailed test
          # pmwstratless
          # pmwgreaterstrat
          # if ((pmwstratless == 1) & (pmwgreaterstrat == 1)) { stop(f("Error: Both pmwstratless and pmwgreaterstrat are 1: {myptm} {mydis} medianptmstrat {medianptmstrat} mediancontrolstrat {mediancontrolstrat} pmwstratless {pmwstratless} pmwgreaterstrat {pmwgreaterstrat}")) }
          # if ((pmwstratless != 1) & (pmwgreaterstrat != 1)) { stop(f("Error: Neither pmwless nor pmwgreater is 1: {myptm} {mydis} medianptmstrat {medianptmstrat} mediancontrolstrat {mediancontrolstrat} pmwstratless {pmwstratless} pmwgreaterstrat {pmwgreaterstrat}")) }
          # if ((pmwstratless == 1) & (pmwgreaterstrat < 1)) {
          if ((pmwgreaterstrat < 0.05) & (sigstrat != "")) {
            sigstrat <- f("({sigstrat})")
          }
        }
        
        # Print
        if (filter %in% c("filter_no_control_residues_strat", "filter_too_few_control_residues_strat", "filter_strat", "filter_strat_same_proteins")) {
          print(f(" >> {source} >> {evorate} >> {filter} >> {mydis} >> {myptm}\t >> {pmw}\t >> {sig} >> strat {pmwstrat}\t >> strat {sigstrat}\t >> strat_minsize {strat_minsize} >> fraction retained {strat_minsize_fraction_retained} ({strat_minsize_retained} / {strat_minsize_total})"))
        } else {
          print(f(" >> {source} >> {evorate} >> {filter} >> {mydis} >> {myptm}\t >> {pmw}\t >> {sig}"))
        }
        
        nptm <- qtmp %>% filter(ptm == myptm & dis == mydis & type == "P") %>% unique %>% nrow
        ncontrol <- qtmp %>% filter(ptm == myptm & dis == mydis & type == "C") %>% unique %>% nrow
        
        if (filter %in% c("filter_no_control_residues_strat", "filter_too_few_control_residues_strat", "filter_strat", "filter_strat_same_proteins")) {
          pvalues <- bind_rows(pvalues, tibble_row(ptm=myptm, 
                                                 dis=mydis, 
                                                 medianptm=medianptm, 
                                                 mediancontrol=mediancontrol, 
                                                 meanptm=meanptm, 
                                                 meancontrol=meancontrol, 
                                                 medianptmstrat=medianptmstrat, 
                                                 mediancontrolstrat=mediancontrolstrat, 
                                                 meanptmstrat=meanptmstrat, 
                                                 meancontrolstrat=meancontrolstrat, 
                                                 nptm=nptm, 
                                                 nptm_pre=nptm_pre, 
                                                 ncontrol=ncontrol, 
                                                 ncontrol_pre=ncontrol_pre, 
                                                 pmw=pmw, 
                                                 pmwstrat=pmwstrat,
                                                 pmwgreater=pmwgreater,
                                                 pmwgreaterstrat=pmwgreaterstrat,
                                                 sig=sig,
                                                 sigstrat=sigstrat,
                                                 strat_minsize=strat_minsize,
                                                 strat_minsize_retained=strat_minsize_retained,
                                                 strat_minsize_total=strat_minsize_total,
                                                 strat_minsize_fraction_retained=strat_minsize_fraction_retained,
          ))
        } else {
          pvalues <- bind_rows(pvalues, tibble_row(ptm=myptm, 
                                                   dis=mydis, 
                                                   medianptm=medianptm, 
                                                   mediancontrol=mediancontrol, 
                                                   meanptm=meanptm, 
                                                   meancontrol=meancontrol, 
                                                   # medianptmstrat=medianptmstrat, 
                                                   # mediancontrolstrat=mediancontrolstrat, 
                                                   # meanptmstrat=meanptmstrat, 
                                                   # meancontrolstrat=meancontrolstrat, 
                                                   nptm=nptm, 
                                                   nptm_pre=nptm_pre, 
                                                   ncontrol=ncontrol, 
                                                   ncontrol_pre=ncontrol_pre, 
                                                   pmw=pmw, 
                                                   # pmwstrat=pmwstrat,
                                                   pmwgreater=pmwgreater,
                                                   # pmwgreaterstrat=pmwgreaterstrat,
                                                   sig=sig,
                                                   # sigstrat=sigstrat,
          ))
        }
      }
    }
    pvalues %<>% arrange(desc(nptm))

    # Process pvalues
    pvalues$dis <- factor(pvalues$dis, levels=c("strcore", "strsurf", "dissurf"))
    # levels(pvalues$dis) <- list("All sites"="A", "Buried"="strcore", "Structured"="strsurf", "Disordered"="dissurf")
    levels(pvalues$dis) <- list("Buried"="strcore", "Structured"="strsurf", "Disordered"="dissurf")

    # Write to RDS
    print(f("Writing pvalues to '{tmprds}'..."))
    write_rds(pvalues, tmprds)
  }

  return(pvalues)
}




# Main analysis function
# analyse("all", f("{evorate}_einsi_tree_1para"), "filter_none", "none")
source <- "all"
# evorate <- "capra0_einsi_tree_1para"
# evorate <- "lichtarge_einsi_tree_1para"
evorate <- "rate4site_einsi_tree_1para"
pred <- "AlphaFold_pLDDT70"
# pvalue_type <- "load_resampled"
# pvalue_type <- "filter_none"
# pvalue_type <- "filter_no_control_residues"
# pvalue_type <- "acc_strat_mean_medianline"
pvalue_type <- "acc_strat_mean_medianline_same_proteins"
# multtest <- "none"    # No multiple testing correction
multtest <- "BH"    # Benjamini-Hochberg = FDR
# multtest <- "holm"  # More powerful alternative to Bonferroni (which is always underpowered)
strat_minsize <- 2

analyse <- function(source = "all", evorate = "rate4site_einsi_tree_1para", pvalue_type = "filter_none", multtest = "none", strat_minsize = 2, pred = "AlphaFold") {

  # Load complete df (all sources, rate4site_einsi_tree_1para, coresurf):
  # q <- tibble(read_tsv("tmp-dataframe-all-human-all-rate4site_einsi_tree_1para-AlphaFold-coresurf.txt"))
  # q <- tibble(read_tsv(f("tmp-dataframe-all-human-{source}-{evorate}-AlphaFold-coresurf.txt")))
  q <- tibble(read_tsv(f("tmp-dataframe-all-human-{source}-{evorate}-{pred}-coresurf.txt")))
  # Get order of PTM sites (by number of sites, descending)
  q
  # via dis=="A":
  q %>% filter(type == "P" & dis == "A") %>% group_by(ptm) %>% tally %>% arrange(desc(n))
  q %>% filter(type == "P" & dis == "A") %>% group_by(ptm) %>% tally %>% arrange(desc(n)) %>% pull(ptm) -> ptm_order
  ptm_order
  # via unique acc|site count:
  q %>% filter(type == "P") %>% select(ptm, acc, site) %>% unique %>% group_by(ptm) %>% tally %>% arrange(desc(n))
  q %>% filter(type == "P") %>% select(ptm, acc, site) %>% unique %>% group_by(ptm) %>% tally %>% arrange(desc(n)) %>% pull(ptm) -> ptm_order
  # >> Same outcome
  ptm_order
  q %<>% filter(dis != "A")
  # Replace "discore" with "strcore" (since residues clearly can't be both buried and disordered, and there are vanishingly few of these. they're definitely not "molten globule" residues or there would be more of them.)
  q %<>% mutate(dis = replace(dis, dis == "discore", "strcore"))
  q$type <- factor(q$type, levels=c("P", "C"))
  q$dis <- factor(q$dis, levels=c("strcore", "strsurf", "dissurf"))
  q
  q %>% summary
  q$dis %>% unique %>% sort
  q %>% filter(dis == "discore")
  q
  
  tmppred <- ""
  if (pred != "AlphaFold") {
    # Filename suffix (e.g. "AlphaFold_pLDDT70")
    tmppred <- f("-{pred}")
  }
  
  # Re-calculate p-values using the specified p-value type
  pvalue_type
  tmp_strat_minsize <- ""
  strat <- 0
  if (pvalue_type %in% c("", "filter_none")) {
    # Calculate p-values without leaving out any proteins (accs)
    mypvalues <- get_pvalues(q, source, evorate, pred = pred, filter = "filter_none")
  } else if (pvalue_type == "filter_no_control_residues") {
    # Calculate p-values and leave out any proteins (accs) without any control residues
    mypvalues <- get_pvalues(q, source, evorate, pred = pred, filter = "filter_no_control_residues")
  } else if (pvalue_type == "filter_no_control_residues_strat") {
    # Calculate p-values and leave out any proteins (accs) without any control residues, and get acc-stratified values below
    mypvalues <- get_pvalues(q, source, evorate, pred = pred, filter = "filter_no_control_residues_strat", strat_minsize = strat_minsize)
    tmp_strat_minsize <- f("-filter_no_control_residues_strat-strat_minsize{strat_minsize}")
    strat <- 1
  } else if (pvalue_type == "filter_too_few_control_residues_strat") {
    # Calculate p-values and leave out any proteins (accs) with too few control residues (fewer than PTM residues), and get acc-stratified values below
    mypvalues <- get_pvalues(q, source, evorate, pred = pred, filter = "filter_too_few_control_residues_strat", strat_minsize = strat_minsize)
    tmp_strat_minsize <- f("-filter_too_few_control_residues_strat-strat_minsize{strat_minsize}")
    strat <- 1
  } else if (pvalue_type == "filter_strat") {
    # Calculate p-values and get acc-stratified values below
    mypvalues <- get_pvalues(q, source, evorate, pred = pred, filter = "filter_strat", strat_minsize = strat_minsize)
    tmp_strat_minsize <- f("-filter_strat-strat_minsize{strat_minsize}")
    strat <- 1
  } else if (pvalue_type == "filter_too_few_control_residues") {
    # Calculate p-values and leave out any proteins (accs) with too few control residues (fewer than PTM residues)
    mypvalues <- get_pvalues(q, source, evorate, pred = pred, filter = "filter_too_few_control_residues")
  } else if (pvalue_type == "load_resampled") {
    # Load resampled p-values from ~/pipeline/evolutionary_rate_analysis/draw.pl (copied from ~/pipeline/evolutionary_rate_analysis/output)
    mypvalues <- dget(f("tmp-R-pvalues-10000-human-{source}-{evorate}-{pred}-coresurf.txt"))
    # Remove "All sites"
    mypvalues %<>% filter(dis != "All sites")
    # Remove "Surface ", turning "Surface Structured" and "Surface Disordered" into "Structured" and "Disordered" as in the other plots
    # mypvalues %<>% mutate(dis = str_replace(dis, "^Surface ", ""))
    # Rename levels instead, since this is a factor
    levels(mypvalues$dis) <- list("Buried"="Buried", "Structured"="Surface Structured", "Disordered"="Surface Disordered")
  } else if (pvalue_type %in% c("acc_strat_mean", "acc_strat_mean_medianline")) {
    # Plot everything at the acc level (averaging across sites)
    # strat_minsize will ensure that each acc has PTM and Control residues, so no need for explicit filter_no_control_residues
    q %<>% group_by(acc, ptm, dis, type) %>% summarise(rate = mean(rate))
    # Use filter_strat to get acc-stratified p-values
    # Calculate p-values and get acc-stratified values below
    mypvalues <- get_pvalues(q, source, evorate, pred = pred, filter = "filter_strat", strat_minsize = strat_minsize)
    tmp_strat_minsize <- f("-filter_strat-strat_minsize{strat_minsize}")
    strat <- 1
  } else if (pvalue_type %in% c("acc_strat_mean_same_proteins", "acc_strat_mean_medianline_same_proteins")) {
    # Plot everything at the acc level (averaging across sites)
    # strat_minsize will ensure that each acc has PTM and Control residues, so no need for explicit filter_no_control_residues
    q %<>% group_by(acc, ptm, dis, type) %>% summarise(rate = mean(rate))
    q
    qp_ptmdisaccs <- q %>% filter(type == "P") %>% select(ptm, dis, acc) %>% distinct
    qc_ptmdisaccs <- q %>% filter(type == "C") %>% select(ptm, dis, acc) %>% distinct
    # qp_ptmdisaccs %>% nrow
    # qc_ptmdisaccs %>% nrow
    # qp_ptmdisaccs %>% inner_join(qc_ptmdisaccs, by = c("ptm", "dis", "acc")) %>% nrow
    # qc_ptmdisaccs %>% inner_join(qp_ptmdisaccs, by = c("ptm", "dis", "acc")) %>% nrow
    # qp_ptmdisaccs %>% inner_join(qc_ptmdisaccs) %>% nrow
    # qc_ptmdisaccs %>% inner_join(qp_ptmdisaccs) %>% nrow
    # qp_ptmdisaccs %>% inner_join(qc_ptmdisaccs) %>% distinct %>% arrange(ptm, dis, acc)
    # qc_ptmdisaccs %>% inner_join(qp_ptmdisaccs) %>% distinct %>% arrange(ptm, dis, acc)
    # identical(qp_ptmdisaccs %>% inner_join(qc_ptmdisaccs) %>% distinct %>% arrange(ptm, dis, acc),
    #           qc_ptmdisaccs %>% inner_join(qp_ptmdisaccs) %>% distinct %>% arrange(ptm, dis, acc))
    # q
    # Filter q
    q %<>% inner_join(qp_ptmdisaccs %>% inner_join(qc_ptmdisaccs) %>% distinct, by = c("ptm", "dis", "acc"))
    q
    # Use filter_strat_same_proteins (just a new name for the same filter_strat to create new pvalues files) to get acc-stratified p-values
    # Calculate p-values and get acc-stratified values below
    mypvalues <- get_pvalues(q, source, evorate, pred = pred, filter = "filter_strat_same_proteins", strat_minsize = strat_minsize)
    tmp_strat_minsize <- f("-filter_strat_same_proteins-strat_minsize{strat_minsize}")
    strat <- 1
  } else if (pvalue_type == "acc_strat_median") {
    # Plot everything at the acc level (averaging across sites)
    q %<>% group_by(acc, ptm, dis, type) %>% summarise(rate = median(rate))
    # Use filter_strat to get acc-stratified p-values
    # Calculate p-values and get acc-stratified values below
    mypvalues <- get_pvalues(q, source, evorate, pred = pred, filter = "filter_strat", strat_minsize = strat_minsize)
    tmp_strat_minsize <- f("-filter_strat-strat_minsize{strat_minsize}")
    strat <- 1
  }
  
  # pLDDT
  my_plddt_thresh <- 0
  if (str_detect(pred, "pLDDT([0-9]+)")) {
    my_plddt_thresh <- str_match(pred, "pLDDT([0-9]+)")[2]
    # pLDDT-filtered: Don't plot disordered
    q %<>% filter(dis != "dissurf")
    mypvalues %<>% filter(dis != "Disordered")
  }
  # PAE
  my_pae_thresh <- 0
  if (str_detect(pred, "PAE([0-9]+)")) {
    my_pae_thresh <- str_match(pred, "PAE([0-9]+)")[2]
    # PAE-filtered: Don't plot disordered
    q %<>% filter(dis != "dissurf")
    mypvalues %<>% filter(dis != "Disordered")
  }
  
  # # Diagnostics:
  # "In scale_x_log10(expand = c(0, 0)) :
  # log-10 transformation introduced infinite values."
  # >> This is actually from geom_text with Inf, not from evorate (rate). 
  # print("rate == 0:")
  # print(q %>% filter(rate == 0))
  # print("rate < 0:")
  # print(q %>% filter(rate < 0))
  # print("log10(rate) == Inf or -Inf:")
  # print(q %>% filter(log10(rate) == Inf | log10(rate) == -Inf))
  print("summary(q$rate):")
  print(q %>% pull(rate) %>% summary)

  # mypvalues <- get_pvalues()
  # mypvalues_filtered_no_control_residues <- get_pvalues(filter = "filter_no_control_residues")
  # mypvalues_filtered_too_few_control_residues <- get_pvalues(filter = "filter_too_few_control_residues")
  # mypvalues_filtered_both <- get_pvalues(filter = "filter_both")
  # identical(mypvalues_filtered_too_few_control_residues, mypvalues_filtered_both)
  # # >> TRUE, i.e. no point in running "filter_both"
  
  # mypvalues_filtered_too_few_control_residues
  # mypvalues_filtered_no_control_residues
  # mypvalues
  
  # # Table showing the effect of filtering out proteins with no or with too few control sites:
  # bind_rows(
  #   mypvalues %>% select(ptm, dis, sig, nptm) %>% mutate(filter = "none"),
  #   mypvalues_filtered_no_control_residues %>% select(ptm, dis, sig, nptm) %>% mutate(filter = "no_controls"),
  #   mypvalues_filtered_too_few_control_residues %>% select(ptm, dis, sig, nptm) %>% mutate(filter = "too_few_controls"),
  # ) %>%
  #   pivot_wider(names_from = filter, values_from = c(nptm, sig)) %>%
  #   arrange(desc(dis == "Disordered"), dis, desc(nptm_none), ptm) %>%
  #   print(n = 100)
  
  # Multiple testing correction of choice
  if (pvalue_type == "load_resampled") {
     # Resampled p-values: don't expect pmwgreater (use meandif instead)
    mypvalues %<>% mutate(pmwadj = p.adjust(pmean, method = multtest), pmwgreater = ifelse(meandif < 0, 1, 0), pmwgreateradj = ifelse(meandif < 0, 1, 0), sigadj = mapply(significance_label, pmwadj, pmwgreateradj))
  } else if (pvalue_type %in% c("filter_no_control_residues_strat", "filter_too_few_control_residues_strat", "filter_strat", "acc_strat_mean", "acc_strat_mean_medianline", "acc_strat_median", "acc_strat_mean_same_proteins", "acc_strat_mean_medianline_same_proteins", "acc_strat_median_same_proteins")) {
    # stratified by acc
    mypvalues %<>% mutate(pmwstratadj = p.adjust(pmwstrat, method = multtest), pmwgreateradj = p.adjust(pmwgreaterstrat, method = multtest), sigadj = mapply(significance_label, pmwstratadj, pmwgreateradj))
  } else {
    # default
    mypvalues %<>% mutate(pmwadj = p.adjust(pmw, method = multtest), pmwgreateradj = p.adjust(pmwgreater, method = multtest), sigadj = mapply(significance_label, pmwadj, pmwgreateradj))
  }

  # # Independent hypothesis weighting (IHW), i.e. a weighted Benjamini-Hochberg correction:
  # # https://www.nature.com/articles/nmeth.3885
  # # All examples given by Wolfgang Huber use alpha = 0.1 rather than 0.05.
  # # https://bioconductor.org/packages/release/bioc/vignettes/IHW/inst/doc/introduction_to_ihw.html
  # # "For rank-based tests (e.g., Wilcoxon) we can use any function that does not depend on the order of arguments (Bourgon, Gentleman, and Huber 2010: https://www.pnas.org/doi/full/10.1073/pnas.0914005107).
  # # n does not depend on the order of arguments.
  # mypvalues %>% mutate(pmwadj = adj_pvalues(ihw(pmw ~ nptm,  data = mypvalues, alpha = 0.1)), pmwgreateradj = adj_pvalues(ihw(pmwgreater ~ nptm,  data = mypvalues, alpha = 0.1)), sigadj = mapply(significance_label, pmwadj, pmwgreateradj))
  # mypvalues
  # # myihw <- ihw(pmw ~ nptm,  data = mypvalues, alpha = 0.1)
  # # >> leads to nbins(myihw) == 1, which means ihw is equivalent to BH
  # # Instead, the covariate needs to be a factor (categorical).
  # myihw <- ihw(pmw ~ dis,  data = mypvalues, alpha = 0.1)
  # # This says: "We recommend that you supply (many) more than 1000 p-values for meaningful data-driven hypothesis weighting results."
  # # >> Presumably IHW is over the top to use here.
  
  
  
  
  
  
  
  # Draw density plot
  # levels(q$dis) <- list("All sites"="A", "Buried"="strcore", "Structured"="strsurf", "Disordered"="dissurf")
  # levels(q$dis) <- list("Buried"="strcore", "Structured"="strsurf", "Disordered"="dissurf")
  # levels(q$dis)
  mycoord <- coord_cartesian()
  mybounds <- c(-Inf, Inf)
  if (my_plddt_thresh != 0 & my_pae_thresh != 0) {
    tmpmyxlab <- f(" (pLDDT ≥ {my_plddt_thresh}, PAE ≤ {my_pae_thresh})")
  } else if (my_plddt_thresh != 0) {
    tmpmyxlab <- f(" (pLDDT ≥ {my_plddt_thresh})")
  } else {
    tmpmyxlab <- ""
  }
  # Get x-axis based on evorate type
  if (str_detect(evorate, "^capra")) {
    myx <- scale_x_log10(expand = c(0, 0), breaks = pretty_breaks(5))
    myy <- scale_y_continuous(expand = c(0, 0.2), breaks = c(0, 1))
    if (my_plddt_thresh == 0) {
      mycoord <- coord_cartesian(xlim = c(NA, 1.19))
    }
    # myxlab <- xlab("Variation score (inverse Jensen-Shannon divergence)")
    # myxlab <- xlab("Variation score (inverse JSD)")
    myxlab <- xlab(f("Variation score (inverse JSD){tmpmyxlab}"))
    my_sig_left <- 0
  } else if (str_detect(evorate, "^lichtarge")) {
    myx <- scale_x_continuous(expand = c(0, 0), trans = "log2", breaks = trans_breaks("log2", function(x) 2^x))
    myy <- scale_y_continuous(expand = c(0, 0.2), breaks = c(0, 1))
    # myxlab <- xlab("Variation score (rvET)")
    myxlab <- xlab(f("Variation score (rvET){tmpmyxlab}"))
    if (my_plddt_thresh == 0) {
      mycoord <- coord_cartesian(xlim = c(NA, 127))
    }
    my_sig_left <- 0
  } else if (str_detect(evorate, "^rate4site")) {
    # myx <- scale_x_continuous(limits = c(2e-3, 9e1), expand = c(0, 0), trans = "log10", breaks = trans_breaks("log10", function(x) 10^x, n=4), labels = trans_format("log10", math_format(10^.x)))
    # myx <- scale_x_continuous(limits = c(4e-3, 3e2), expand = c(0, 0), trans = "log10", breaks = trans_breaks("log10", function(x) 10^x, n=4), labels = c(0.001, 0.01, 0.1, 1, 10, 100, 1000))
    # myx <- scale_x_continuous(limits = c(2e-2, 8e1), expand = c(0, 0), trans = "log10", breaks = trans_breaks("log10", function(x) 10^x, n=4), labels = c(0.01, 0.1, 1, 10, 100))
    # myx <- scale_x_continuous(expand = c(0, 0), trans = "log10", breaks = trans_breaks("log10", function(x) 10^x, n=4), labels = c(0.01, 0.1, 1, 10, 100))
    myx <- scale_x_log10(expand = c(0, 0), breaks = c(0.01, 0.1, 1, 10), labels = c("0.01", "0.1", "1", "10"))
    # myx <- scale_x_log10(expand = c(0, 0))
    # myx <- scale_x_continuous(expand = c(0, 0), trans = "log10", breaks = trans_breaks("log10", function(x) 10^x, n=4), labels = c(0.01, 0.1, 1, 10, 100))
    # myx <- scale_x_continuous(limits = c(1e-4, 9e1), expand = c(0, 0), trans = "log10", breaks = trans_breaks("log10", function(x) 10^x, n=5), labels = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100))
    # myx <- scale_x_continuous(limits = c(1e-4, 14e1), expand = c(0, 0), trans = "log10", breaks = trans_breaks("log10", function(x) 10^x, n=6), labels = trans_format("log10", math_format(10^.x)))
    # myx <- scale_x_continuous(expand = c(0, 0), trans = "log10", breaks = trans_breaks("log10", function(x) 10^x))
    # myx <- scale_x_continuous(expand = c(0, 0), breaks = pretty_breaks(5))
    # myy <- scale_y_continuous(expand = c(0, 0.2), breaks = c(0, 1))
    myy <- scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), breaks = c(0, 1))
    # myxlab <- xlab("Evolutionary rate (Rate4Site)")
    # IMPORTANT Adjusting data so <2e-2 is set to 2e-2, so that we plot all the data:
    q %<>% mutate(rate = ifelse(rate < 0.02, 0.02, rate))
    # q %<>% mutate(rate = ifelse(rate < 1e-2, 1e-2, rate))
    # q %<>% mutate(rate = ifelse(rate > 4e1, 4e1, rate))
    # mycoord <- coord_cartesian(xlim = c(2e-2, 3e1))
    # mycoord <- coord_cartesian(xlim = c(0.02, 4))
    # mybounds <- c(0.02, 4)
    myxlab <- xlab(f("Evolutionary rate{tmpmyxlab}"))
    if (my_plddt_thresh == 0) {
      mycoord <- coord_cartesian(xlim = c(0.02, 40))
    }
    my_sig_left <- 0
  } else if (str_detect(evorate, "^norm_capra")) {
    myx <- scale_x_continuous(expand = c(0, 0), breaks = pretty_breaks(5))
    myy <- scale_y_continuous(expand = c(0, 0.2), breaks = c(0, 1))
    # myxlab <- xlab("Normalised variation score (inverse JSD)")
    myxlab <- xlab(f("Normalised inverse JSD{tmpmyxlab}"))
    mycoord <- coord_cartesian(xlim = c(-2.5, 3))
    my_sig_left <- -Inf
  } else if (str_detect(evorate, "^norm_lichtarge")) {
    myx <- scale_x_continuous(expand = c(0, 0), breaks = pretty_breaks(5))
    myy <- scale_y_continuous(expand = c(0, 0.2), breaks = c(0, 1))
    # myxlab <- xlab("Normalised variation score (rvET)")
    myxlab <- xlab(f("Normalised rvET{tmpmyxlab}"))
    mycoord <- coord_cartesian(xlim = c(-2.5, 3))
    my_sig_left <- -Inf
  } else if (str_detect(evorate, "^norm_rate4site")) {
    myx <- scale_x_continuous(expand = c(0, 0), breaks = pretty_breaks(5))
    myy <- scale_y_continuous(expand = c(0, 0.2), breaks = c(0, 1))
    myxlab <- xlab(f("Normalised evolutionary rate{tmpmyxlab}"))
    mycoord <- coord_cartesian(xlim = c(-2, 2.5))
    my_sig_left <- -Inf
  }
  
  
  if (strat == 0) {
    # Show number of PTM sites (default)
    my_n_ptm <- geom_text(data=mypvalues, aes(x=Inf, y=Inf, label=paste("\nn=", comma(nptm), " ", sep=""), hjust=1, vjust=0.5), size = 5/.pt)
    my_sig <- geom_text(data=mypvalues, aes(x=my_sig_left, y=Inf, label=paste("\n", " ", sigadj, sep=""), hjust=0, vjust=0.5), size = 6/.pt)
    my_vline_ptm <- geom_vline(data=mypvalues, aes(xintercept=medianptm), color = ptmb, linewidth = 1/.pt, linetype = "solid")
    my_vline_control <- geom_vline(data=mypvalues, aes(xintercept=mediancontrol), color = ptmg, linewidth = 1/.pt, linetype = "solid")
  } else {
    # Stratifying by acc: Show number of retained accs (proteins)
    my_n_ptm <- geom_text(data=mypvalues, aes(x=Inf, y=Inf, label=paste("\nn=", comma(strat_minsize_retained), " ", sep=""), hjust=1, vjust=0.5), size = 5/.pt)
    # my_sig <- geom_text(data=mypvalues, aes(x=0, y=Inf, label=paste("\n", " ", sigadj, sep=""), hjust=0, vjust=0.5), size = 6/.pt)
    my_sig <- geom_text(data=mypvalues, aes(x=my_sig_left, y=Inf, label=paste("\n", " ", sigadj, sep=""), hjust=0, vjust=0.5), size = 5/.pt)
    # my_vline_ptm <- geom_vline(data=mypvalues, aes(xintercept=medianptmstrat), color = ptmb, linewidth = 1/.pt, linetype = "solid")
    # my_vline_control <- geom_vline(data=mypvalues, aes(xintercept=mediancontrolstrat), color = ptmg, linewidth = 1/.pt, linetype = "solid")
    if (pvalue_type == "acc_strat_mean") {
      # Show mean lines
      my_vline_ptm <- geom_vline(data=mypvalues, aes(xintercept=meanptmstrat), color = ptmb, linewidth = 1/.pt, linetype = "solid")
      my_vline_control <- geom_vline(data=mypvalues, aes(xintercept=meancontrolstrat), color = ptmg, linewidth = 1/.pt, linetype = "solid")
    } else {
      # Show median lines
      my_vline_ptm <- geom_vline(data=mypvalues, aes(xintercept=medianptmstrat), color = ptmb, linewidth = 1/.pt, linetype = "solid")
      my_vline_control <- geom_vline(data=mypvalues, aes(xintercept=mediancontrolstrat), color = ptmg, linewidth = 1/.pt, linetype = "solid")
    }
  }
  # Plot
  q %>% 
    mutate(type = fct_recode(type, "Modified"="P", "Control"="C")) %>%
    mutate(dis = fct_recode(dis, "Buried"="strcore", "Structured"="strsurf", "Disordered"="dissurf")) %>%
    ggplot(aes(rate)) +
    # geom_density(aes(colour = fct_rev(type), fill = type, alpha = 0.2), adjust = 1, linewidth = 1/.pt) +
    # geom_density(aes(colour = fct_rev(type), fill = fct_rev(type)), alpha = 0.3, adjust = 1, linewidth = 1/.pt) +
    geom_density(aes(colour = fct_rev(type), fill = fct_rev(type)), alpha = 0.3, bounds = mybounds) +
    # geom_density(aes(colour = fct_rev(type), alpha = 0.2), adjust = 1, linewidth = 1/.pt) +
    
    # # Draw lines at mean
    # geom_vline(data=mypvalues, aes(xintercept=meanptm), color = ptmb, linetype = "dashed") +
    # geom_vline(data=mypvalues, aes(xintercept=meancontrol), color = ptmg, linetype = "dashed") +
    
    # Draw lines at median
    # geom_vline(data=mypvalues, aes(xintercept=medianptm), color = ptmb, linewidth = 1/.pt, linetype = "solid") +
    # geom_vline(data=mypvalues, aes(xintercept=mediancontrol), color = ptmg, linewidth = 1/.pt, linetype = "solid") +
    
    my_vline_ptm +
    my_vline_control +
    
    # # Draw lines at median (orange / blue)
    # geom_vline(data=mypvalues, aes(xintercept=mediancontrol), color = ptmbd, linewidth = 1/.pt, linetype = "solid") +
    # geom_vline(data=mypvalues, aes(xintercept=medianptm), color = ptmod, linewidth = 1/.pt, linetype = "solid") +
    
    # Significance: *** etc.
    # geom_text(data=mypvalues, aes(x=0, y=Inf, label=paste("\n", " ", sigadj, sep=""), hjust=0, vjust=0.5), size = 6/.pt) +
    # geom_text(data=mypvalues, aes(x=0, y=Inf, label=paste("\n", " ", sigadj, sep=""), hjust=0, vjust=0.5), size = 5/.pt) +
    # Significance: *** etc. (stratified by acc)
    # geom_text(data=mypvalues, aes(x=0, y=Inf, label=paste("\n", " ", sigadjstrat, sep=""), hjust=0, vjust=0.5), size = 6/.pt) +
    
    my_sig +
    
    my_n_ptm +
    
    # Showing nptm_pre and nptm:
    # geom_text(data=mypvalues, aes(x=Inf, y=Inf, label=paste("\n\npre=", nptm_pre, " ", "\nn=", nptm, " ", sep=""), hjust=1, vjust=0.5), size = 2.5) +
    # >> This shows that filtering out accs without any control residues reduces the number of PTM sites by around 10%.
    
    # Showing ncontrol_pre and ncontrol:
    # geom_text(data=mypvalues, aes(x=Inf, y=Inf, label=paste("\n\nnc_pre=", ncontrol_pre, " ", "\nnc=", ncontrol, " ", sep=""), hjust=1, vjust=0.5), size = 2) +
    # >> Shows that control residues are not affected by the filtering done above (where we remove accs that don't have any control residues for a given ptm|dis|surf), as expected.
    
    # Showing nptm and ncontrol:
    # geom_text(data=mypvalues, aes(x=Inf, y=Inf, label=paste("\n\nn=", nptm, " ", "\nnc=", ncontrol, " ", sep=""), hjust=1, vjust=0.5), size = 2) +
    # >> Shows that the number of control residues is generally a good amount higher than the number of PTM sites (by a factor of 1.5 - 10, normally around 2).
    
    # log10 x-axis ticks (would appear below each facet)
    # annotation_logticks(sides = "b", size = 0.5/.weight, short = unit(1, "points"), mid = unit(1, "points"), long = unit(1, "points"),) +
    
    # scale_x_continuous(limits = c(2e-3, 9e1), expand = c(0, 0), trans = "log10", breaks = trans_breaks("log10", function(x) 10^x, n=4), labels = trans_format("log10", math_format(10^.x))) +
    # scale_x_continuous(limits = c(8e-3, 3e2), expand = c(0, 0), trans = "log10", breaks = trans_breaks("log10", function(x) 10^x, n=4), labels = trans_format("log10", math_format(10^.x))) +
    # scale_x_continuous(limits = c(2e-3, 3e2), expand = c(0, 0), trans = "log10", breaks = trans_breaks("log10", function(x) 10^x, n=4), labels = trans_format("log10", math_format(10^.x))) +
    # scale_x_continuous(limits = c(2e-3, 4e2), expand = c(0, 0), trans = "log10", breaks = trans_breaks("log10", function(x) 10^x, n=4), labels = trans_format("log10", math_format(10^.x))) +
    # scale_x_continuous(limits = c(4e-3, 3e2), expand = c(0, 0), trans = "log10", breaks = trans_breaks("log10", function(x) 10^x, n=4), labels = trans_format("log10", math_format(10^.x))) +
    # Manually labelled x-axis
    # scale_x_continuous(limits = c(4e-3, 3e2), expand = c(0, 0), trans = "log10", breaks = trans_breaks("log10", function(x) 10^x, n=4), labels = c(0.001, 0.01, 0.1, 1, 10, 100, 1000)) +
    myx +
    # scale_y_continuous(expand = c(0, 0.2)) +
    # scale_y_continuous(expand = c(0, 0.2), breaks = pretty_breaks(2)) +
    myy +
    mycoord +
    # scale_colour_manual(values = c("P" = "#4f81bd", "C" = "#bec1c0")) +
    # scale_fill_manual(values = c("P" = "#4f81bd", "C" = "#bec1c0")) +
    scale_colour_manual(values = c("Modified" = ptmb, "Control" = ptmg), aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
    # scale_colour_manual(values = c("Modified" = ptmod, "Control" = ptmbd), aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
    # facet_grid(factor(ptm, levels = ptm_order) ~ fct_relevel(dis, "Buried", "Structured", "Disordered")) +
    # facet_grid(factor(ptm, levels = ptm_order) ~ fct_relevel(dis, "Buried", "Structured", "Disordered"), scales = "free_x") +
    facet_grid(factor(ptm, levels = ptm_order) ~ fct_relevel(dis, "Buried", "Structured", "Disordered"), scales = "free_y") +
    # xlab("Evolutionary rate (Rate4Site)") +
    myxlab +
    ylab("Probability density") +
    # labs(tag = "a") +
    # labs(title = "rate4site_einsi_tree_1para AlphaFold all") +
    theme_nature() +
    # theme(axis.text = element_text(size = 5)) +
    # # theme(legend.position = "none", strip.background = element_rect(fill = NA, colour = NA), axis.text.y = element_text(size=7), axis.title.y = element_text(angle=90), strip.text.y = element_text(angle=0))
    # # theme(legend.position = "none", strip.background = element_rect(fill = NA, colour = NA), axis.title.y = element_text(angle=90), strip.text.y = element_text(angle=0))
    # theme(legend.position = c(1, -0.06), legend.direction = "horizontal", legend.background = element_blank())
    theme(panel.spacing.x = unit(1, "mm")) +
    theme(panel.spacing.y = unit(1, "mm")) +
    theme(axis.title.x = element_text(hjust = 0.1))
  # qsave(file = "output-density.pdf", width = 61 * (q$dis %>% unique %>% length), height = 24.385 + (length(unique(q$ptm)) * 13.095), units="mm")		# 17 PTMs = 24.385 + 17 * 13.095 = 247 mm, and 61 mm width per column
  # qsave(file = "output-density.pdf", width = 90, height = 17 * 10, units="mm")
  # qsave(file = "output-density.pdf", width = 90, height = 120, units="mm")
  # qsave(file = "output-density.pdf", width = 90, height = 100, units="mm")
  qsave(f("output-density-{source}-{evorate}-{pvalue_type}-{multtest}{tmp_strat_minsize}{tmppred}.pdf"), width = 90, height = 120, units="mm")
  
  
  # Draw density plot (select subset of 11 PTMs)
  # Initialise plot
  ptm_order_select <- c("S-p", "T-p", "Y-p", "K-ub", "K-sum", "K-ac", "K-mal", "K-suc", "K-me", "R-me", "N-gly")
  mypvalues_select <- mypvalues %>% filter(ptm %in% ptm_order_select)
  # # BH FDR correction (affects output here: K-suc Structured * -> '')
  # mypvalues_select %<>% mutate(pmwadj = p.adjust(pmw, method = "BH"), pmwgreateradj = p.adjust(pmwgreater, method = "BH"), sigadj = mapply(significance_label, pmwadj, pmwgreateradj))
  # Plot
  if (strat == 0) {
    # Show number of PTM sites (default)
    my_n_ptm <- geom_text(data=mypvalues_select, aes(x=Inf, y=Inf, label=paste("\nn=", comma(nptm), " ", sep=""), hjust=1, vjust=0.5), size = 5/.pt)
    my_sig <- geom_text(data=mypvalues_select, aes(x=my_sig_left, y=Inf, label=paste("\n", " ", sigadj, sep=""), hjust=0, vjust=0.5), size = 6/.pt)
    my_vline_ptm <- geom_vline(data=mypvalues_select, aes(xintercept=medianptm), color = ptmb, linewidth = 1/.pt, linetype = "solid")
    my_vline_control <- geom_vline(data=mypvalues_select, aes(xintercept=mediancontrol), color = ptmg, linewidth = 1/.pt, linetype = "solid")
  } else {
    # Stratifying by acc: Show number of retained accs (proteins)
    my_n_ptm <- geom_text(data=mypvalues_select, aes(x=Inf, y=Inf, label=paste("\nn=", comma(strat_minsize_retained), " ", sep=""), hjust=1, vjust=0.5), size = 5/.pt)
    # my_sig <- geom_text(data=mypvalues_select, aes(x=0, y=Inf, label=paste("\n", " ", sigadj, sep=""), hjust=0, vjust=0.5), size = 6/.pt)
    my_sig <- geom_text(data=mypvalues_select, aes(x=my_sig_left, y=Inf, label=paste("\n", " ", sigadj, sep=""), hjust=0, vjust=0.5), size = 5/.pt)
    # my_vline_ptm <- geom_vline(data=mypvalues_select, aes(xintercept=medianptmstrat), color = ptmb, linewidth = 1/.pt, linetype = "solid")
    # my_vline_control <- geom_vline(data=mypvalues_select, aes(xintercept=mediancontrolstrat), color = ptmg, linewidth = 1/.pt, linetype = "solid")
    if (pvalue_type == "acc_strat_mean") {
      # Show mean lines
      my_vline_ptm <- geom_vline(data=mypvalues_select, aes(xintercept=meanptmstrat), color = ptmb, linewidth = 1/.pt, linetype = "solid")
      my_vline_control <- geom_vline(data=mypvalues_select, aes(xintercept=meancontrolstrat), color = ptmg, linewidth = 1/.pt, linetype = "solid")
    } else {
      # Show median lines
      my_vline_ptm <- geom_vline(data=mypvalues_select, aes(xintercept=medianptmstrat), color = ptmb, linewidth = 1/.pt, linetype = "solid")
      my_vline_control <- geom_vline(data=mypvalues_select, aes(xintercept=mediancontrolstrat), color = ptmg, linewidth = 1/.pt, linetype = "solid")
    }
  }
  q %>% 
    filter(ptm %in% ptm_order_select) %>%
    mutate(type = fct_recode(type, "Modified"="P", "Control"="C")) %>%
    mutate(dis = fct_recode(dis, "Buried"="strcore", "Structured"="strsurf", "Disordered"="dissurf")) %>%
    ggplot(aes(rate)) +
    # geom_density(aes(colour = fct_rev(type), fill = fct_rev(type)), alpha = 0.3, adjust = 1, linewidth = 1/.pt) +
    geom_density(aes(colour = fct_rev(type), fill = fct_rev(type)), alpha = 0.3, bounds = mybounds) +
    # geom_vline(data=mypvalues_select, aes(xintercept=medianptm), color = ptmb, linewidth = 1/.pt, linetype = "solid") +
    # geom_vline(data=mypvalues_select, aes(xintercept=mediancontrol), color = ptmg, linewidth = 1/.pt, linetype = "solid") +
    # # geom_vline(data=mypvalues_select, aes(xintercept=mediancontrol), color = ptmbd, linewidth = 1/.pt, linetype = "solid") +
    # # geom_vline(data=mypvalues_select, aes(xintercept=medianptm), color = ptmod, linewidth = 1/.pt, linetype = "solid") +
    my_vline_ptm +
    my_vline_control +
    # geom_text(data=mypvalues_select, aes(x=0, y=Inf, label=paste("\n", " ", sigadj, sep=""), hjust=0, vjust=0.5), size = 6/.pt) +
    # # geom_text(data=mypvalues_select, aes(x=Inf, y=Inf, label=paste("\nn=", nptm, " ", sep=""), hjust=1, vjust=0.5), size = 5/.pt) +
    my_sig +
    my_n_ptm +
    # scale_x_continuous(limits = c(4e-3, 3e2), expand = c(0, 0), trans = "log10", breaks = trans_breaks("log10", function(x) 10^x, n=4), labels = c(0.001, 0.01, 0.1, 1, 10, 100, 1000)) +
    # scale_x_continuous(limits = c(2e-2, 9e1), expand = c(0, 0), trans = "log10", breaks = trans_breaks("log10", function(x) 10^x, n=4), labels = c(0.01, 0.1, 1, 10, 100)) +
    myx +
    # scale_y_continuous(expand = c(0, 0.2), breaks = pretty_breaks(2)) +
    myy +
    mycoord +
    scale_colour_manual(values = c("Modified" = ptmb, "Control" = ptmg), aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
    # scale_colour_manual(values = c("Modified" = ptmod, "Control" = ptmbd), aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
    # facet_grid(factor(ptm, levels = ptm_order_select) ~ fct_relevel(dis, "Buried", "Structured", "Disordered")) +
    # facet_grid(factor(ptm, levels = ptm_order_select) ~ fct_relevel(dis, "Buried", "Structured", "Disordered"), scales = "free_x") +
    facet_grid(factor(ptm, levels = ptm_order_select) ~ fct_relevel(dis, "Buried", "Structured", "Disordered"), scales = "free_y") +
    # xlab("Evolutionary rate (Rate4Site)") +
    myxlab +
    ylab("Probability density") +
    theme_nature() +
    # theme(axis.text = element_text(size = 5)) +
    # theme(legend.position = c(1, -0.08), legend.direction = "horizontal", legend.background = element_blank())
    theme(panel.spacing.x = unit(1, "mm")) +
    theme(panel.spacing.y = unit(1, "mm")) +
    theme(axis.title.x = element_text(hjust = 0.1))
  # qsave(file = "output-density-select11.pdf", width = 90, height = 70)
  qsave(f("output-density-select11-{source}-{evorate}-{pvalue_type}-{multtest}{tmp_strat_minsize}{tmppred}.pdf"), width = 90, height = 100)
  
  
  # Draw density plot (select subset of 9 PTMs)
  # Initialise plot
  ptm_order_select <- c("S-p", "T-p", "Y-p", "K-ub", "K-ac", "K-mal", "K-suc", "K-me", "R-me")
  mypvalues_select <- mypvalues %>% filter(ptm %in% ptm_order_select)
  # # BH FDR correction (does not affect output here)
  # mypvalues_select %<>% mutate(pmwadj = p.adjust(pmw, method = "BH"), pmwgreateradj = p.adjust(pmwgreater, method = "BH"), sigadj = mapply(significance_label, pmwadj, pmwgreateradj))
  if (strat == 0) {
    # Show number of PTM sites (default)
    # my_n_ptm <- geom_text(data=mypvalues_select, aes(x=Inf, y=Inf, label=paste("\nn=", comma(nptm), " ", sep=""), hjust=1, vjust=0.5), size = 5/.pt)
    my_n_ptm <- geom_text(data=mypvalues_select, aes(x=Inf, y=Inf, label=paste("\nn=", comma(nptm), " ", sep=""), hjust=0.5, vjust=0.5), size = 5/.pt)
    my_sig <- geom_text(data=mypvalues_select, aes(x=my_sig_left, y=Inf, label=paste("\n", " ", sigadj, sep=""), hjust=0, vjust=0.5), size = 6/.pt)
    # my_sig <- geom_text(data=mypvalues_select, aes(x=0, y=Inf, label=paste("\n", " ", sigadj, sep=""), hjust=0, vjust=0.5), size = 5/.pt)
    my_vline_ptm <- geom_vline(data=mypvalues_select, aes(xintercept=medianptm), color = ptmb, linewidth = 1/.pt, linetype = "solid")
    my_vline_control <- geom_vline(data=mypvalues_select, aes(xintercept=mediancontrol), color = ptmg, linewidth = 1/.pt, linetype = "solid")
  } else {
    # Stratifying by acc: Show number of retained accs (proteins)
    my_n_ptm <- geom_text(data=mypvalues_select, aes(x=Inf, y=Inf, label=paste("\nn=", comma(strat_minsize_retained), " ", sep=""), hjust=1, vjust=0.5), size = 5/.pt)
    my_sig <- geom_text(data=mypvalues_select, aes(x=my_sig_left, y=Inf, label=paste("\n", " ", sigadj, sep=""), hjust=0, vjust=0.5), size = 6/.pt)
    # my_vline_ptm <- geom_vline(data=mypvalues_select, aes(xintercept=medianptmstrat), color = ptmb, linewidth = 1/.pt, linetype = "solid")
    # my_vline_control <- geom_vline(data=mypvalues_select, aes(xintercept=mediancontrolstrat), color = ptmg, linewidth = 1/.pt, linetype = "solid")
    if (pvalue_type == "acc_strat_mean") {
      # Show mean lines
      my_vline_ptm <- geom_vline(data=mypvalues_select, aes(xintercept=meanptmstrat), color = ptmb, linewidth = 1/.pt, linetype = "solid")
      my_vline_control <- geom_vline(data=mypvalues_select, aes(xintercept=meancontrolstrat), color = ptmg, linewidth = 1/.pt, linetype = "solid")
    } else {
      # # Show median lines
      # my_vline_ptm <- geom_vline(data=mypvalues_select, aes(xintercept=medianptmstrat), color = ptmb, linewidth = 1/.pt, linetype = "solid")
      # my_vline_control <- geom_vline(data=mypvalues_select, aes(xintercept=mediancontrolstrat), color = ptmg, linewidth = 1/.pt, linetype = "solid")
      # Show median lines
      my_vline_ptm <- geom_vline(data=mypvalues_select, aes(xintercept=medianptmstrat), color = ptmb, linetype = "solid")
      my_vline_control <- geom_vline(data=mypvalues_select, aes(xintercept=mediancontrolstrat), color = ptmg, linetype = "solid")
    }
  }
  # Plot
  (
    q %>% 
      filter(ptm %in% ptm_order_select) %>%
      mutate(type = fct_recode(type, "Modified"="P", "Control"="C")) %>%
      mutate(dis = fct_recode(dis, "Buried"="strcore", "Structured"="strsurf", "Disordered"="dissurf")) %>%
      ggplot(aes(rate)) +
      # geom_density(aes(colour = fct_rev(type), fill = fct_rev(type)), alpha = 0.3, adjust = 1, linewidth = 0.75/.weight) +
      geom_density(aes(colour = fct_rev(type), fill = fct_rev(type)), alpha = 0.3, bounds = mybounds) +
      # # geom_vline(data=mypvalues_select, aes(xintercept=medianptm), color = ptmb, linewidth = 0.75/.weight, linetype = "solid") +
      # # geom_vline(data=mypvalues_select, aes(xintercept=mediancontrol), color = ptmg, linewidth = 0.75/.weight, linetype = "solid") +
      # # geom_vline(data=mypvalues_select, aes(xintercept=mediancontrol), color = ptmbd, linewidth = 0.75/.weight, linetype = "solid") +
      # # geom_vline(data=mypvalues_select, aes(xintercept=medianptm), color = ptmod, linewidth = 0.75/.weight, linetype = "solid") +
      # # geom_vline(data=mypvalues_select, aes(xintercept=mediancontrol), color = scales::viridis_pal()(2)[1], linewidth = 0.75/.weight, linetype = "solid") +
      # # geom_vline(data=mypvalues_select, aes(xintercept=medianptm), color = scales::viridis_pal()(2)[2], linewidth = 0.75/.weight, linetype = "solid") +
      # # geom_vline(data=mypvalues_select, aes(xintercept=mediancontrol), color = scales::viridis_pal()(9)[4], linewidth = 0.75/.weight, linetype = "solid") +
      # # geom_vline(data=mypvalues_select, aes(xintercept=medianptm), color = scales::viridis_pal()(9)[8], linewidth = 0.75/.weight, linetype = "solid") +
      # geom_vline(data=mypvalues_select, aes(xintercept=mediancontrol), color = ptmvir0, linewidth = 0.75/.weight, linetype = "solid") +
      # geom_vline(data=mypvalues_select, aes(xintercept=medianptm), color = ptmvir1, linewidth = 0.75/.weight, linetype = "solid") +
      my_vline_ptm +
      my_vline_control +
      # geom_text(data=mypvalues_select, aes(x=0, y=Inf, label=paste("\n", " ", sigadj, sep=""), hjust=0, vjust=0.5), size = 6/.pt) +
      # # geom_text(data=mypvalues_select, aes(x=Inf, y=Inf, label=paste("\nn=", nptm, " ", sep=""), hjust=1, vjust=0.5), size = 5/.pt) +
      # # geom_text(data=mypvalues_select, aes(x=Inf, y=Inf, label=paste("\nn=", nptm, " ", sep=""), hjust=1, vjust=0.5), size = 5/.pt) +
      my_sig +
      my_n_ptm +
      # scale_x_continuous(limits = c(4e-3, 3e2), expand = c(0, 0), trans = "log10", breaks = trans_breaks("log10", function(x) 10^x, n=4), labels = c(0.001, 0.01, 0.1, 1, 10, 100, 1000)) +
      # scale_x_continuous(limits = c(2e-2, 9e1), expand = c(0, 0), trans = "log10", breaks = trans_breaks("log10", function(x) 10^x, n=4), labels = c(0.01, 0.1, 1, 10, 100)) +
      myx +
      # scale_y_continuous(expand = c(0, 0.2), breaks = pretty_breaks(2)) +
      myy +
      mycoord +
      # scale_colour_manual(values = c("Modified" = ptmb, "Control" = ptmg), aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
      # scale_colour_manual(values = c("Modified" = ptmod, "Control" = ptmbd), aesthetics = c("colour", "fill"), name = NULL , guide = guide_legend(reverse = T)) +
      # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
      # scale_colour_manual(values = c("Modified" = scales::viridis_pal()(9)[8], "Control" = scales::viridis_pal()(9)[4]), aesthetics = c("colour", "fill"), name = NULL , guide = guide_legend(reverse = T)) +
      scale_colour_manual(values = c("Modified" = ptmvir1, "Control" = ptmvir0), aesthetics = c("colour", "fill"), name = NULL , guide = guide_legend(reverse = T)) +
      # facet_grid(factor(ptm, levels = ptm_order_select) ~ fct_relevel(dis, "Buried", "Structured", "Disordered")) +
      # facet_grid(factor(ptm, levels = ptm_order_select) ~ fct_relevel(dis, "Buried", "Structured", "Disordered"), scales = "free_x") +
      facet_grid(factor(ptm, levels = ptm_order_select) ~ fct_relevel(dis, "Buried", "Structured", "Disordered"), scales = "free_y") +
      # xlab("Evolutionary rate (Rate4Site)") +
      myxlab +
      ylab("Probability density") +
      theme_nature() +
      theme(panel.spacing.x = unit(1, "mm")) +
      theme(panel.spacing.y = unit(1, "mm")) +
      theme(axis.title.x = element_text(hjust = 0.1))
  ) %>%
    # qsave("output-density-select9.pdf", width = 90, height = 60)
    qsave(f("output-density-select9-{source}-{evorate}-{pvalue_type}-{multtest}{tmp_strat_minsize}{tmppred}.pdf"), width = 90, height = 80)
  # q
  
  # # Draw simplified evorate plot (not faceting by PTM type)
  # Doing this in alphasa_relasa_vs_ptms.R instead because I have weighted controls there (aaweight -> q$freq)
  # (
  #   q %>% 
  #     mutate(type = fct_recode(type, "Modified"="P", "Control"="C")) %>%
  #     mutate(dis = fct_recode(dis, "Buried"="strcore", "Structured"="strsurf", "Disordered"="dissurf")) %>%
  #     ggplot(aes(rate)) +
  #     geom_density(aes(colour = fct_rev(type), fill = fct_rev(type)), alpha = 0.3, adjust = 1, linewidth = 0.75/.weight) +
  #     # geom_vline(data=mypvalues, aes(xintercept=mediancontrol), color = ptmvir0, linewidth = 0.75/.weight, linetype = "solid") +
  #     # geom_vline(data=mypvalues, aes(xintercept=medianptm), color = ptmvir1, linewidth = 0.75/.weight, linetype = "solid") +
  #     # geom_text(data=mypvalues, aes(x=0, y=Inf, label=paste("\n", " ", sigadj, sep=""), hjust=0, vjust=0.5), size = 6/.pt) +
  #     # geom_text(data=mypvalues, aes(x=Inf, y=Inf, label=paste("\nn=", nptm, " ", sep=""), hjust=1, vjust=0.5), size = 5/.pt) +
  #     myx +
  #     scale_y_continuous(expand = c(0, 0.2), breaks = pretty_breaks(2)) +
  #     scale_colour_manual(values = c("Modified" = ptmvir1, "Control" = ptmvir0), aesthetics = c("colour", "fill"), name = NULL , guide = guide_legend(reverse = T)) +
  #     facet_grid(. ~ fct_relevel(dis, "Buried", "Structured", "Disordered")) +
  #     myxlab +
  #     ylab("Probability density") +
  #     theme_nature()
  # ) %>%
  #   qsave(f("output-density-noptmfacet-{source}-{evorate}-{pvalue_type}-{multtest}.pdf"), width = 90, height = 60)
}
# End main analysis function

# Analyses:
# Cycle through evorate types for mafftmode einsi_tree_1para
for (evorate in c("rate4site", "lichtarge", "capra0")) {
  analyse("all", f("{evorate}_einsi_tree_1para"), "filter_no_control_residues", "BH")
}

# rate4site

# Cycle through sources for evorate and mafftmode rate4site_einsi_tree_1para
for (source in c("UniProt", "Ochoa", "PhosphoSitePlus", "dbPTM")) {
  analyse(source, "rate4site_einsi_tree_1para", "filter_no_control_residues", "BH")
}
# Cycle through p-value types for rate4site_einsi_tree_1para
# for (pvalue_type in c("load_resampled", "filter_too_few_control_residues", "filter_no_control_residues", "filter_none")) {
for (pvalue_type in c("filter_too_few_control_residues", "filter_no_control_residues", "filter_none")) {
  analyse("all", "rate4site_einsi_tree_1para", pvalue_type, "BH")
}
# Cycle through multiple testing corrections for rate4site_einsi_tree_1para
for (multtest in c("none", "BH", "holm")) {
  analyse("all", "rate4site_einsi_tree_1para", "filter_no_control_residues", multtest)
}
# Cycle through mafftmodes for evorate rate4site
for (mafftmode in c("einsi_tree_1para", "einsi_tree", "einsi")) {
  analyse("all", f("rate4site_{mafftmode}"), "filter_no_control_residues", "BH")
}
# Cycle through mafftmodes for evorate rate4site
for (mafftmode in c("einsi_tree", "linsi_tree", "ginsi_tree", "einsi", "linsi", "ginsi")) {
  analyse("all", f("rate4site_{mafftmode}"), "filter_no_control_residues", "BH")
}

# lichtarge

# Cycle through sources for evorate and mafftmode lichtarge_einsi_tree_1para
for (source in c("UniProt", "Ochoa", "PhosphoSitePlus", "dbPTM")) {
  analyse(source, "lichtarge_einsi_tree_1para", "filter_no_control_residues", "BH")
}
# Cycle through p-value types for lichtarge_einsi_tree_1para
# for (pvalue_type in c("load_resampled", "filter_too_few_control_residues", "filter_no_control_residues", "filter_none")) {
for (pvalue_type in c("filter_too_few_control_residues", "filter_no_control_residues", "filter_none")) {
  analyse("all", "lichtarge_einsi_tree_1para", pvalue_type, "BH")
}
# Cycle through multiple testing corrections for lichtarge_einsi_tree_1para
for (multtest in c("none", "BH", "holm")) {
  analyse("all", "lichtarge_einsi_tree_1para", "filter_no_control_residues", multtest)
}
# Cycle through mafftmodes for evorate lichtarge
for (mafftmode in c("einsi_tree_1para", "einsi_tree", "einsi")) {
  analyse("all", f("lichtarge_{mafftmode}"), "filter_no_control_residues", "BH")
}
# Cycle through mafftmodes for evorate lichtarge
for (mafftmode in c("einsi_tree", "linsi_tree", "ginsi_tree", "einsi", "linsi", "ginsi")) {
  analyse("all", f("lichtarge_{mafftmode}"), "filter_no_control_residues", "BH")
}

# capra0

# Cycle through sources for evorate and mafftmode capra0_einsi_tree_1para
for (source in c("UniProt", "Ochoa", "PhosphoSitePlus", "dbPTM")) {
  analyse(source, "capra0_einsi_tree_1para", "filter_no_control_residues", "BH")
}
# Cycle through p-value types for capra0_einsi_tree_1para
# for (pvalue_type in c("load_resampled", "filter_too_few_control_residues", "filter_no_control_residues", "filter_none")) {
for (pvalue_type in c("filter_too_few_control_residues", "filter_no_control_residues", "filter_none")) {
  analyse("all", "capra0_einsi_tree_1para", pvalue_type, "BH")
}
# Cycle through multiple testing corrections for capra0_einsi_tree_1para
for (multtest in c("none", "BH", "holm")) {
  analyse("all", "capra0_einsi_tree_1para", "filter_no_control_residues", multtest)
}
# Cycle through mafftmodes for evorate capra0
for (mafftmode in c("einsi_tree_1para", "einsi_tree", "einsi")) {
  analyse("all", f("capra0_{mafftmode}"), "filter_no_control_residues", "BH")
}
# Cycle through mafftmodes for evorate capra0
for (mafftmode in c("einsi_tree", "linsi_tree", "ginsi_tree", "einsi", "linsi", "ginsi")) {
  analyse("all", f("capra0_{mafftmode}"), "filter_no_control_residues", "BH")
}

# main (already run above)
analyse("all", "rate4site_einsi_tree_1para", "filter_no_control_residues", "BH")
analyse("all", "lichtarge_einsi_tree_1para", "filter_no_control_residues", "BH")
analyse("all", "capra0_einsi_tree_1para", "filter_no_control_residues", "BH")



# stratified by acc using resampling-based wilcox_test from coin package (much slower)
# strat_minsize 2:
analyse("all", "rate4site_einsi_tree_1para", "filter_strat", "BH", strat_minsize=2)
analyse("all", "lichtarge_einsi_tree_1para", "filter_strat", "BH", strat_minsize=2)
analyse("all", "capra0_einsi_tree_1para", "filter_strat", "BH", strat_minsize=2)
# strat_minsize 2:
analyse("all", "rate4site_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=2)
analyse("all", "lichtarge_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=2)
analyse("all", "capra0_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=2)
# strat_minsize 3:
analyse("all", "rate4site_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=3)
analyse("all", "lichtarge_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=3)
analyse("all", "capra0_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=3)
# strat_minsize 4:
analyse("all", "rate4site_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=4)
analyse("all", "lichtarge_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=4)
analyse("all", "capra0_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=4)
# strat_minsize 5:
analyse("all", "rate4site_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=5)
analyse("all", "lichtarge_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=5)
analyse("all", "capra0_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=5)
# strat_minsize 7:
analyse("all", "rate4site_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=7)
analyse("all", "lichtarge_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=7)
analyse("all", "capra0_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=7)
# strat_minsize 10:
analyse("all", "rate4site_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=10)
analyse("all", "lichtarge_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=10)
analyse("all", "capra0_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=10)
# strat_minsize 15:
analyse("all", "rate4site_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=15)
analyse("all", "lichtarge_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=15)
analyse("all", "capra0_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=15)
# strat_minsize 20:
analyse("all", "rate4site_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=20)
analyse("all", "lichtarge_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=20)
analyse("all", "capra0_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=20)
# strat_minsize 50:
analyse("all", "rate4site_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=50)
analyse("all", "lichtarge_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=50)
analyse("all", "capra0_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=50)
# strat_minsize 100:
analyse("all", "rate4site_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=100)
analyse("all", "lichtarge_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=100)
analyse("all", "capra0_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=100)
# >> almost nothing is retained

# filter_too_few_control strat_minsize 2:
analyse("all", "rate4site_einsi_tree_1para", "filter_too_few_control_residues_strat", "BH", strat_minsize=2)
analyse("all", "lichtarge_einsi_tree_1para", "filter_too_few_control_residues_strat", "BH", strat_minsize=2)
analyse("all", "capra0_einsi_tree_1para", "filter_too_few_control_residues_strat", "BH", strat_minsize=2)

# filter_too_few_control strat_minsize 5:
analyse("all", "rate4site_einsi_tree_1para", "filter_too_few_control_residues_strat", "BH", strat_minsize=5)
analyse("all", "lichtarge_einsi_tree_1para", "filter_too_few_control_residues_strat", "BH", strat_minsize=5)
analyse("all", "capra0_einsi_tree_1para", "filter_too_few_control_residues_strat", "BH", strat_minsize=5)

# filter_too_few_control strat_minsize 10:
analyse("all", "rate4site_einsi_tree_1para", "filter_too_few_control_residues_strat", "BH", strat_minsize=10)
analyse("all", "lichtarge_einsi_tree_1para", "filter_too_few_control_residues_strat", "BH", strat_minsize=10)
analyse("all", "capra0_einsi_tree_1para", "filter_too_few_control_residues_strat", "BH", strat_minsize=10)

# normalised strat_minsize 2:
analyse("all", "norm_rate4site_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=2)
analyse("all", "norm_lichtarge_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=2)
analyse("all", "norm_capra0_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=2)

# normalised strat_minsize 5:
analyse("all", "norm_rate4site_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=5)
analyse("all", "norm_lichtarge_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=5)
analyse("all", "norm_capra0_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=5)

# normalised strat_minsize 10:
analyse("all", "norm_rate4site_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=10)
analyse("all", "norm_lichtarge_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=10)
analyse("all", "norm_capra0_einsi_tree_1para", "filter_no_control_residues_strat", "BH", strat_minsize=10)

# normalised
analyse("all", "norm_rate4site_einsi_tree_1para", "filter_no_control_residues", "BH")
analyse("all", "norm_lichtarge_einsi_tree_1para", "filter_no_control_residues", "BH")
analyse("all", "norm_capra0_einsi_tree_1para", "filter_no_control_residues", "BH")

# normalised
analyse("all", "norm_rate4site_einsi_tree_1para", "filter_too_few_control_residues", "BH")
analyse("all", "norm_lichtarge_einsi_tree_1para", "filter_too_few_control_residues", "BH")
analyse("all", "norm_capra0_einsi_tree_1para", "filter_too_few_control_residues", "BH")

# Plot at the acc level (mean across each protein)
analyse("all", "rate4site_einsi_tree_1para", "acc_strat_mean", "BH", strat_minsize=2)
analyse("all", "lichtarge_einsi_tree_1para", "acc_strat_mean", "BH", strat_minsize=2)
analyse("all", "capra0_einsi_tree_1para", "acc_strat_mean", "BH", strat_minsize=2)

# Plot at the acc level (median across each protein)
analyse("all", "rate4site_einsi_tree_1para", "acc_strat_median", "BH", strat_minsize=2)
analyse("all", "lichtarge_einsi_tree_1para", "acc_strat_median", "BH", strat_minsize=2)
analyse("all", "capra0_einsi_tree_1para", "acc_strat_median", "BH", strat_minsize=2)

# Plot at the acc level (mean across each protein), but draw lines at the median
analyse("all", "rate4site_einsi_tree_1para", "acc_strat_mean_medianline", "BH", strat_minsize=2)
analyse("all", "lichtarge_einsi_tree_1para", "acc_strat_mean_medianline", "BH", strat_minsize=2)
analyse("all", "capra0_einsi_tree_1para", "acc_strat_mean_medianline", "BH", strat_minsize=2)

# Plot at the acc level (mean across each protein), but draw lines at the median
analyse("all", "rate4site_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2)
analyse("all", "lichtarge_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2)
analyse("all", "capra0_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2)

# Plot at the acc level (mean across each protein), but draw lines at the median, normalised
analyse("all", "norm_rate4site_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2)
analyse("all", "norm_lichtarge_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2)
analyse("all", "norm_capra0_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2)

# Plot at the acc level (mean across each protein), but draw lines at the median, and filter for pLDDT ≥ 70
analyse("all", "rate4site_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2, "AlphaFold_pLDDT70")
analyse("all", "lichtarge_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2, "AlphaFold_pLDDT70")
analyse("all", "capra0_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2, "AlphaFold_pLDDT70")

# Plot at the acc level (mean across each protein), but draw lines at the median, and filter for pLDDT ≥ 90
analyse("all", "rate4site_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2, "AlphaFold_pLDDT90")
analyse("all", "lichtarge_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2, "AlphaFold_pLDDT90")
analyse("all", "capra0_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2, "AlphaFold_pLDDT90")

# Plot at the acc level (mean across each protein), but draw lines at the median, and filter for pLDDT ≥ 70, normalised
analyse("all", "norm_rate4site_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2, "AlphaFold_pLDDT70")
analyse("all", "norm_lichtarge_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2, "AlphaFold_pLDDT70")
analyse("all", "norm_capra0_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2, "AlphaFold_pLDDT70")

# Plot at the acc level (mean across each protein), but draw lines at the median, and filter for pLDDT ≥ 90
analyse("all", "rate4site_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2, "AlphaFold_pLDDT90")
analyse("all", "lichtarge_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2, "AlphaFold_pLDDT90")
analyse("all", "capra0_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2, "AlphaFold_pLDDT90")

# Plot at the acc level (mean across each protein), but draw lines at the median, and filter for pLDDT ≥ 70 and PAE ≥ 2 (>> MAIN EXPORT FIGURE 2 IS RATE4SITE)
analyse("all", "rate4site_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2, "AlphaFold_pLDDT70_PAE2")
analyse("all", "lichtarge_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2, "AlphaFold_pLDDT70_PAE2")
analyse("all", "capra0_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2, "AlphaFold_pLDDT70_PAE2")

# Plot at the acc level (mean across each protein), but draw lines at the median, and filter for pLDDT ≥ 90 and PAE ≤ 1
analyse("all", "rate4site_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2, "AlphaFold_pLDDT90_PAE1")
analyse("all", "lichtarge_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2, "AlphaFold_pLDDT90_PAE1")
analyse("all", "capra0_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2, "AlphaFold_pLDDT90_PAE1")

# Plot at the acc level (mean across each protein), but draw lines at the median, and filter for pLDDT ≥ 70 and PAE ≥ 2, normalised
analyse("all", "norm_rate4site_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2, "AlphaFold_pLDDT70_PAE2")
analyse("all", "norm_lichtarge_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2, "AlphaFold_pLDDT70_PAE2")
analyse("all", "norm_capra0_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2, "AlphaFold_pLDDT70_PAE2")

# Plot at the acc level (mean across each protein), but draw lines at the median, and filter for pLDDT ≥ 90 and PAE ≥ 1, normalised
analyse("all", "norm_rate4site_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2, "AlphaFold_pLDDT90_PAE1")
analyse("all", "norm_lichtarge_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2, "AlphaFold_pLDDT90_PAE1")
analyse("all", "norm_capra0_einsi_tree_1para", "acc_strat_mean_medianline_same_proteins", "BH", strat_minsize=2, "AlphaFold_pLDDT90_PAE1")
