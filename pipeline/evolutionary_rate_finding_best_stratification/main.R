#!/usr/bin/env Rscript --vanilla

# initialize
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("doBy"))
suppressPackageStartupMessages(library("coin"))
suppressPackageStartupMessages(library("tidyverse"))

# Clear workspace
# rm(list=ls())

# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
nosurf <- args[1]
test <- args[2]
alternative <- args[3]
source <- args[4]
evorate <- args[5]
mafftmode <- args[6]
pred <- args[7]
minsize <- args[8]
disfilt <- args[9]

# Main function
GetResults <- function (test, file) {
  
#   test <- "wilcoxsign_test"
# #   evorate <- "capra0"
# #   evorate <- "lichtarge"
#   evorate <- "rate4site"
#   file <- paste("../../evolutionary_rate_analysis/tmp/tmp-dataframe-all-human-", evorate, "_linsi_tree-MULTICOM.txt", sep="")

  table <- read.delim(file, header=T, quote="", stringsAsFactors = TRUE)
#   str(table)
  head(table)

  # type: Make sure P comes before C (factor order)
  table$type <- factor(table$type, levels=c("P", "C"))
  
  table <- table[table$dis!="A",]
  table$dis <- factor(as.character(table$dis))
  table <- table
  
  # Apply disfilt to e.g. analyse only buried residues if it's "strcore"
  if (disfilt != "all")
  {
    if (disfilt == "strcore") {
      # also include discore in tmpdisfilt
      tmpdisfilt <- c("strcore", "discore")
    } else {
      tmpdisfilt <- c(disfilt)
    }
    table <- table[table$dis %in% tmpdisfilt,]
  }
  
  # Dis/Str
  table$disstr <- NA
#   str(table$disstr)
  if (nrow(table[table$dis %in% c("dissurf", "discore", "strsurf", "strcore"),]) > 0)
  {
    # if not -nosurf:
    try(table[table$dis=="dissurf",]$disstr <- "dis", silent=T)
    # table[table$dis=="discore",]$disstr <- "dis"
    try(table[table$dis=="discore",]$disstr <- "str", silent=T)
    try(table[table$dis=="strsurf",]$disstr <- "str", silent=T)
    try(table[table$dis=="strcore",]$disstr <- "str", silent=T)
  }
  else
  {
    # if -nosurf:
    try(table[table$dis=="dis",]$disstr <- "dis", silent=T)
    try(table[table$dis=="str",]$disstr <- "str", silent=T)
  }
  table$disstr <- factor(as.character(table$disstr))
#   str(table$disstr)
  
  # Core/Surf
  table$coresurf <- NA
#   str(table$coresurf)
  if (nrow(table[table$dis %in% c("dissurf", "discore", "strsurf", "strcore"),]) > 0)
  {
    # if not -nosurf:
    # table[table$dis=="discore",]$coresurf <- "core"
    try(table[table$dis=="discore",]$coresurf <- "core", silent=T)
    try(table[table$dis=="strcore",]$coresurf <- "core", silent=T)
    try(table[table$dis=="dissurf",]$coresurf <- "surf", silent=T)
    try(table[table$dis=="strsurf",]$coresurf <- "surf", silent=T)
  }
  else
  {
    # if -nosurf:
    # table[table$dis=="dis",]$coresurf <- NA
    # table[table$dis=="str",]$coresurf <- NA
    # Using NA leads to a "do not contain data" warning, so use "n/a" as a string instead
    try(table[table$dis=="dis",]$coresurf <- "n/a", silent=T)
    try(table[table$dis=="str",]$coresurf <- "n/a", silent=T)
  }
  table$coresurf <- factor(as.character(table$coresurf))
#   str(table$coresurf)
  
  # ptm-disstr
  table$ptmdisstr <- as.factor(paste(table$ptm, table$disstr))
  # ptm-coresurf
  table$ptmcoresurf <- as.factor(paste(table$ptm, table$coresurf))
  # ptm-acc
  table$ptmacc <- as.factor(paste(table$ptm, table$acc))
  # disstr-coresurf
  table$disstrcoresurf <- as.factor(paste(table$disstr, table$coresurf))
  # disstr-acc
  table$disstracc <- as.factor(paste(table$disstr, table$acc))
  # coresurf-acc
  table$coresurfacc <- as.factor(paste(table$coresurf, table$acc))
  # ptm-disstr-coresurf
  table$ptmdisstrcoresurf <- as.factor(paste(table$ptm, table$disstr, table$coresurf))
  # ptm-disstr-acc
  table$ptmdisstracc <- as.factor(paste(table$ptm, table$disstr, table$acc))
  # ptm-coresurf-acc
  table$ptmcoresurfacc <- as.factor(paste(table$ptm, table$coresurf, table$acc))
  # disstr-coresurf-acc
  table$disstrcoresurfacc <- as.factor(paste(table$disstr, table$coresurf, table$acc))
  # ptm-disstr-coresurf-acc
  table$ptmdisstrcoresurfacc <- as.factor(paste(table$ptm, table$disstr, table$coresurf, table$acc))
  
  # Acc
  table$acc <- factor(as.character(table$acc))
  
  # Save processed table to rds file for debugging
  # print(file)
  # options("max.print" = 100000000)
  # print(table)
  # write_rds(table, file = "tmp-table.rds")
  # # write_rds(table, file = "tmp-dataframe-all-human-Ochoa-rate4site_einsi_tree_1para-AlphaFold-nosurf.txt.rds")
  
  # Progress bar
#   bar <- txtProgressBar(min = 0, max = 11, style = 3)
  bar <- txtProgressBar(min = 0, max = 16, style = 3)
  
#   # Get complete block design
#   controls <- 1
#   str(table)
#   numtable <- summaryBy(rate ~ acc + type, data = table, FUN = length, keep.accs=TRUE)
#   numtable <- numtable[numtable$type=="C",]
#   numtable$type <- NULL
#   numtable$controls <- numtable$rate
#   numtable$rate <- NULL
#   str(numtable)
#   head(numtable)
#   table <- table[table$acc %in% numtable[numtable$controls==controls,]$acc,]
#   wilcox.test(table[table$type=="P",]$rate, table[table$type=="C",]$rate)
#   wilcoxsign_test(table[table$type=="P",]$rate ~ table[table$type=="C",]$rate)


  # # Tests, with increasingly elaborate stratification including by protein:
  # # To filter out cases where there aren't enough data for a stratification, leading to statistical test crashes if in a stratum n=1:
  # # Replace "table" with "subset(table, !(acc %in% (table %>% group_by(acc) %>% tally %>% filter(n < minsize) %>% pull(acc))))" (and then replace "acc" with the actual stratification in use)
  # results <- data.frame(strat=character(), pval=numeric())
  # # 1
  # print("none")
  # results <- rbind(results, data.frame(strat="none", pval=pvalue(get(test)(rate ~ type, alternative=alternative, data=table))))
  # setTxtProgressBar(bar, 1)
  # print("ptm")
  # # results <- rbind(results, data.frame(strat="ptm", pval=pvalue(get(test)(rate ~ type | ptm, alternative=alternative, data=table))))
  # results <- rbind(results, data.frame(strat="ptm", pval=pvalue(get(test)(rate ~ type | ptm, alternative=alternative, data=subset(table, !(ptm %in% (table %>% group_by(ptm) %>% tally %>% filter(n < minsize) %>% pull(ptm))))))))
  # setTxtProgressBar(bar, 2)
  # print("disstr")
  # # results <- rbind(results, data.frame(strat="disstr", pval=pvalue(get(test)(rate ~ type | disstr, alternative=alternative, data=table))))
  # results <- rbind(results, data.frame(strat="disstr", pval=pvalue(get(test)(rate ~ type | disstr, alternative=alternative, data=subset(table, !(disstr %in% (table %>% group_by(disstr) %>% tally %>% filter(n < minsize) %>% pull(disstr))))))))
  # setTxtProgressBar(bar, 3)
  # print("coresurf")
  # # results <- rbind(results, data.frame(strat="coresurf", pval=pvalue(get(test)(rate ~ type | coresurf, alternative=alternative, data=table))))
  # results <- rbind(results, data.frame(strat="coresurf", pval=pvalue(get(test)(rate ~ type | coresurf, alternative=alternative, data=subset(table, !(coresurf %in% (table %>% group_by(coresurf) %>% tally %>% filter(n < minsize) %>% pull(coresurf))))))))
  # setTxtProgressBar(bar, 4)
  # print("acc")
  # # results <- rbind(results, data.frame(strat="acc", pval=pvalue(get(test)(rate ~ type | acc, alternative=alternative, data=table))))
  # results <- rbind(results, data.frame(strat="acc", pval=pvalue(get(test)(rate ~ type | acc, alternative=alternative, data=subset(table, !(acc %in% (table %>% group_by(acc) %>% tally %>% filter(n < minsize) %>% pull(acc))))))))
  # setTxtProgressBar(bar, 5)
  # # 2
  # print("ptm-disstr")
  # # results <- rbind(results, data.frame(strat="ptm-disstr", pval=pvalue(get(test)(rate ~ type | ptmdisstr, alternative=alternative, data=table))))
  # results <- rbind(results, data.frame(strat="ptm-disstr", pval=pvalue(get(test)(rate ~ type | ptmdisstr, alternative=alternative, data=subset(table, !(ptmdisstr %in% (table %>% group_by(ptmdisstr) %>% tally %>% filter(n < minsize) %>% pull(ptmdisstr))))))))
  # setTxtProgressBar(bar, 6)
  # print("ptm-coresurf")
  # # results <- rbind(results, data.frame(strat="ptm-coresurf", pval=pvalue(get(test)(rate ~ type | ptmcoresurf, alternative=alternative, data=table))))
  # results <- rbind(results, data.frame(strat="ptm-coresurf", pval=pvalue(get(test)(rate ~ type | ptmcoresurf, alternative=alternative, data=subset(table, !(ptmcoresurf %in% (table %>% group_by(ptmcoresurf) %>% tally %>% filter(n < minsize) %>% pull(ptmcoresurf))))))))
  # setTxtProgressBar(bar, 7)
  # print("ptm-acc")
  # # results <- rbind(results, data.frame(strat="ptm-acc", pval=pvalue(get(test)(rate ~ type | ptmacc, alternative=alternative, data=table))))
  # results <- rbind(results, data.frame(strat="ptm-acc", pval=pvalue(get(test)(rate ~ type | ptmacc, alternative=alternative, data=subset(table, !(ptmacc %in% (table %>% group_by(ptmacc) %>% tally %>% filter(n < minsize) %>% pull(ptmacc))))))))
  # setTxtProgressBar(bar, 8)
  # print("disstr-coresurf")
  # # results <- rbind(results, data.frame(strat="disstr-coresurf", pval=pvalue(get(test)(rate ~ type | disstrcoresurf, alternative=alternative, data=table))))
  # results <- rbind(results, data.frame(strat="disstr-coresurf", pval=pvalue(get(test)(rate ~ type | disstrcoresurf, alternative=alternative, data=subset(table, !(disstrcoresurf %in% (table %>% group_by(disstrcoresurf) %>% tally %>% filter(n < minsize) %>% pull(disstrcoresurf))))))))
  # setTxtProgressBar(bar, 9)
  # print("disstr-acc")
  # # results <- rbind(results, data.frame(strat="disstr-acc", pval=pvalue(get(test)(rate ~ type | disstracc, alternative=alternative, data=table))))
  # results <- rbind(results, data.frame(strat="disstr-acc", pval=pvalue(get(test)(rate ~ type | disstracc, alternative=alternative, data=subset(table, !(disstracc %in% (table %>% group_by(disstracc) %>% tally %>% filter(n < minsize) %>% pull(disstracc))))))))
  # setTxtProgressBar(bar, 10)
  # print("coresurf-acc")
  # # results <- rbind(results, data.frame(strat="coresurf-acc", pval=pvalue(get(test)(rate ~ type | coresurfacc, alternative=alternative, data=table))))
  # results <- rbind(results, data.frame(strat="coresurf-acc", pval=pvalue(get(test)(rate ~ type | coresurfacc, alternative=alternative, data=subset(table, !(coresurfacc %in% (table %>% group_by(coresurfacc) %>% tally %>% filter(n < minsize) %>% pull(coresurfacc))))))))
  # setTxtProgressBar(bar, 11)
  # # 3
  # print("ptm-disstr-coresurf")
  # # results <- rbind(results, data.frame(strat="ptm-disstr-coresurf", pval=pvalue(get(test)(rate ~ type | ptmdisstrcoresurf, alternative=alternative, data=table))))
  # results <- rbind(results, data.frame(strat="ptm-disstr-coresurf", pval=pvalue(get(test)(rate ~ type | ptmdisstrcoresurf, alternative=alternative, data=subset(table, !(ptmdisstrcoresurf %in% (table %>% group_by(ptmdisstrcoresurf) %>% tally %>% filter(n < minsize) %>% pull(ptmdisstrcoresurf))))))))
  # setTxtProgressBar(bar, 12)
  # print("ptm-disstr-acc")
  # # results <- rbind(results, data.frame(strat="ptm-disstr-acc", pval=pvalue(get(test)(rate ~ type | ptmdisstracc, alternative=alternative, data=table))))
  # results <- rbind(results, data.frame(strat="ptm-disstr-acc", pval=pvalue(get(test)(rate ~ type | ptmdisstracc, alternative=alternative, data=subset(table, !(ptmdisstracc %in% (table %>% group_by(ptmdisstracc) %>% tally %>% filter(n < minsize) %>% pull(ptmdisstracc))))))))
  # setTxtProgressBar(bar, 13)
  # print("ptm-coresurf-acc")
  # # results <- rbind(results, data.frame(strat="ptm-coresurf-acc", pval=pvalue(get(test)(rate ~ type | ptmcoresurfacc, alternative=alternative, data=table))))
  # results <- rbind(results, data.frame(strat="ptm-coresurf-acc", pval=pvalue(get(test)(rate ~ type | ptmcoresurfacc, alternative=alternative, data=subset(table, !(ptmcoresurfacc %in% (table %>% group_by(ptmcoresurfacc) %>% tally %>% filter(n < minsize) %>% pull(ptmcoresurfacc))))))))
  # setTxtProgressBar(bar, 14)
  # print("disstr-coresurf-acc")
  # # results <- rbind(results, data.frame(strat="disstr-coresurf-acc", pval=pvalue(get(test)(rate ~ type | disstrcoresurfacc, alternative=alternative, data=table))))
  # results <- rbind(results, data.frame(strat="disstr-coresurf-acc", pval=pvalue(get(test)(rate ~ type | disstrcoresurfacc, alternative=alternative, data=subset(table, !(disstrcoresurfacc %in% (table %>% group_by(disstrcoresurfacc) %>% tally %>% filter(n < minsize) %>% pull(disstrcoresurfacc))))))))
  # setTxtProgressBar(bar, 15)
  # # 4
  # print("ptm-disstr-coresurf-acc")
  # # results <- rbind(results, data.frame(strat="ptm-disstr-coresurf-acc", pval=pvalue(get(test)(rate ~ type | ptmdisstrcoresurfacc, alternative=alternative, data=table))))
  # results <- rbind(results, data.frame(strat="ptm-disstr-coresurf-acc", pval=pvalue(get(test)(rate ~ type | ptmdisstrcoresurfacc, alternative=alternative, data=subset(table, !(ptmdisstrcoresurfacc %in% (table %>% group_by(ptmdisstrcoresurfacc) %>% tally %>% filter(n < minsize) %>% pull(ptmdisstrcoresurfacc))))))))
  # setTxtProgressBar(bar, 16)
  




  # With "try":
  # Tests, with increasingly elaborate stratification including by protein:
  # To filter out cases where there aren't enough data for a stratification, leading to statistical test crashes if in a stratum n=1:
  # Replace "table" with "subset(table, !(acc %in% (table %>% group_by(acc) %>% tally %>% filter(n < minsize) %>% pull(acc))))" (and then replace "acc" with the actual stratification in use)
  results <- data.frame(strat=character(), pval=numeric())
  # 1
  print("none")
  mypvalue <- NA
  try(mypvalue <- pvalue(get(test)(rate ~ type, alternative=alternative, data=table)))
  results <- rbind(results, data.frame(strat="none", pval=mypvalue))
  setTxtProgressBar(bar, 1)
  print("ptm")
  # results <- rbind(results, data.frame(strat="ptm", pval=pvalue(get(test)(rate ~ type | ptm, alternative=alternative, data=table))))
  mypvalue <- NA
  try(mypvalue <- pvalue(get(test)(rate ~ type | ptm, alternative=alternative, data=subset(table, !(ptm %in% (table %>% group_by(ptm) %>% tally %>% filter(n < minsize) %>% pull(ptm)))))))
  results <- rbind(results, data.frame(strat="ptm", pval=mypvalue))
  setTxtProgressBar(bar, 2)
  print("disstr")
  # results <- rbind(results, data.frame(strat="disstr", pval=pvalue(get(test)(rate ~ type | disstr, alternative=alternative, data=table))))
  mypvalue <- NA
  try(mypvalue <- pvalue(get(test)(rate ~ type | disstr, alternative=alternative, data=subset(table, !(disstr %in% (table %>% group_by(disstr) %>% tally %>% filter(n < minsize) %>% pull(disstr)))))))
  results <- rbind(results, data.frame(strat="disstr", pval=mypvalue))
  setTxtProgressBar(bar, 3)
  print("coresurf")
  # results <- rbind(results, data.frame(strat="coresurf", pval=pvalue(get(test)(rate ~ type | coresurf, alternative=alternative, data=table))))
  mypvalue <- NA
  try(mypvalue <- pvalue(get(test)(rate ~ type | coresurf, alternative=alternative, data=subset(table, !(coresurf %in% (table %>% group_by(coresurf) %>% tally %>% filter(n < minsize) %>% pull(coresurf)))))))
  results <- rbind(results, data.frame(strat="coresurf", pval=mypvalue))
  setTxtProgressBar(bar, 4)
  print("acc")
  # results <- rbind(results, data.frame(strat="acc", pval=pvalue(get(test)(rate ~ type | acc, alternative=alternative, data=table))))
  mypvalue <- NA
  try(mypvalue <- pvalue(get(test)(rate ~ type | acc, alternative=alternative, data=subset(table, !(acc %in% (table %>% group_by(acc) %>% tally %>% filter(n < minsize) %>% pull(acc)))))))
  results <- rbind(results, data.frame(strat="acc", pval=mypvalue))
  setTxtProgressBar(bar, 5)
  # 2
  print("ptm-disstr")
  # results <- rbind(results, data.frame(strat="ptm-disstr", pval=pvalue(get(test)(rate ~ type | ptmdisstr, alternative=alternative, data=table))))
  mypvalue <- NA
  try(mypvalue <- pvalue(get(test)(rate ~ type | ptmdisstr, alternative=alternative, data=subset(table, !(ptmdisstr %in% (table %>% group_by(ptmdisstr) %>% tally %>% filter(n < minsize) %>% pull(ptmdisstr)))))))
  results <- rbind(results, data.frame(strat="ptm-disstr", pval=mypvalue))
  setTxtProgressBar(bar, 6)
  print("ptm-coresurf")
  # results <- rbind(results, data.frame(strat="ptm-coresurf", pval=pvalue(get(test)(rate ~ type | ptmcoresurf, alternative=alternative, data=table))))
  mypvalue <- NA
  try(mypvalue <- pvalue(get(test)(rate ~ type | ptmcoresurf, alternative=alternative, data=subset(table, !(ptmcoresurf %in% (table %>% group_by(ptmcoresurf) %>% tally %>% filter(n < minsize) %>% pull(ptmcoresurf)))))))
  results <- rbind(results, data.frame(strat="ptm-coresurf", pval=mypvalue))
  setTxtProgressBar(bar, 7)
  print("ptm-acc")
  # results <- rbind(results, data.frame(strat="ptm-acc", pval=pvalue(get(test)(rate ~ type | ptmacc, alternative=alternative, data=table))))
  mypvalue <- NA
  try(mypvalue <- pvalue(get(test)(rate ~ type | ptmacc, alternative=alternative, data=subset(table, !(ptmacc %in% (table %>% group_by(ptmacc) %>% tally %>% filter(n < minsize) %>% pull(ptmacc)))))))
  results <- rbind(results, data.frame(strat="ptm-acc", pval=mypvalue))
  setTxtProgressBar(bar, 8)
  print("disstr-coresurf")
  # results <- rbind(results, data.frame(strat="disstr-coresurf", pval=pvalue(get(test)(rate ~ type | disstrcoresurf, alternative=alternative, data=table))))
  mypvalue <- NA
  try(mypvalue <- pvalue(get(test)(rate ~ type | disstrcoresurf, alternative=alternative, data=subset(table, !(disstrcoresurf %in% (table %>% group_by(disstrcoresurf) %>% tally %>% filter(n < minsize) %>% pull(disstrcoresurf)))))))
  results <- rbind(results, data.frame(strat="disstr-coresurf", pval=mypvalue))
  setTxtProgressBar(bar, 9)
  print("disstr-acc")
  # results <- rbind(results, data.frame(strat="disstr-acc", pval=pvalue(get(test)(rate ~ type | disstracc, alternative=alternative, data=table))))
  mypvalue <- NA
  try(mypvalue <- pvalue(get(test)(rate ~ type | disstracc, alternative=alternative, data=subset(table, !(disstracc %in% (table %>% group_by(disstracc) %>% tally %>% filter(n < minsize) %>% pull(disstracc)))))))
  results <- rbind(results, data.frame(strat="disstr-acc", pval=mypvalue))
  setTxtProgressBar(bar, 10)
  print("coresurf-acc")
  # results <- rbind(results, data.frame(strat="coresurf-acc", pval=pvalue(get(test)(rate ~ type | coresurfacc, alternative=alternative, data=table))))
  mypvalue <- NA
  try(mypvalue <- pvalue(get(test)(rate ~ type | coresurfacc, alternative=alternative, data=subset(table, !(coresurfacc %in% (table %>% group_by(coresurfacc) %>% tally %>% filter(n < minsize) %>% pull(coresurfacc)))))))
  results <- rbind(results, data.frame(strat="coresurf-acc", pval=mypvalue))
  setTxtProgressBar(bar, 11)
  # 3
  print("ptm-disstr-coresurf")
  # results <- rbind(results, data.frame(strat="ptm-disstr-coresurf", pval=pvalue(get(test)(rate ~ type | ptmdisstrcoresurf, alternative=alternative, data=table))))
  mypvalue <- NA
  try(mypvalue <- pvalue(get(test)(rate ~ type | ptmdisstrcoresurf, alternative=alternative, data=subset(table, !(ptmdisstrcoresurf %in% (table %>% group_by(ptmdisstrcoresurf) %>% tally %>% filter(n < minsize) %>% pull(ptmdisstrcoresurf)))))))
  results <- rbind(results, data.frame(strat="ptm-disstr-coresurf", pval=mypvalue))
  setTxtProgressBar(bar, 12)
  print("ptm-disstr-acc")
  # results <- rbind(results, data.frame(strat="ptm-disstr-acc", pval=pvalue(get(test)(rate ~ type | ptmdisstracc, alternative=alternative, data=table))))
  mypvalue <- NA
  try(mypvalue <- pvalue(get(test)(rate ~ type | ptmdisstracc, alternative=alternative, data=subset(table, !(ptmdisstracc %in% (table %>% group_by(ptmdisstracc) %>% tally %>% filter(n < minsize) %>% pull(ptmdisstracc)))))))
  results <- rbind(results, data.frame(strat="ptm-disstr-acc", pval=mypvalue))
  setTxtProgressBar(bar, 13)
  print("ptm-coresurf-acc")
  # results <- rbind(results, data.frame(strat="ptm-coresurf-acc", pval=pvalue(get(test)(rate ~ type | ptmcoresurfacc, alternative=alternative, data=table))))
  mypvalue <- NA
  try(mypvalue <- pvalue(get(test)(rate ~ type | ptmcoresurfacc, alternative=alternative, data=subset(table, !(ptmcoresurfacc %in% (table %>% group_by(ptmcoresurfacc) %>% tally %>% filter(n < minsize) %>% pull(ptmcoresurfacc)))))))
  results <- rbind(results, data.frame(strat="ptm-coresurf-acc", pval=mypvalue))
  setTxtProgressBar(bar, 14)
  print("disstr-coresurf-acc")
  # results <- rbind(results, data.frame(strat="disstr-coresurf-acc", pval=pvalue(get(test)(rate ~ type | disstrcoresurfacc, alternative=alternative, data=table))))
  mypvalue <- NA
  try(mypvalue <- pvalue(get(test)(rate ~ type | disstrcoresurfacc, alternative=alternative, data=subset(table, !(disstrcoresurfacc %in% (table %>% group_by(disstrcoresurfacc) %>% tally %>% filter(n < minsize) %>% pull(disstrcoresurfacc)))))))
  results <- rbind(results, data.frame(strat="disstr-coresurf-acc", pval=mypvalue))
  setTxtProgressBar(bar, 15)
  # 4
  print("ptm-disstr-coresurf-acc")
  # results <- rbind(results, data.frame(strat="ptm-disstr-coresurf-acc", pval=pvalue(get(test)(rate ~ type | ptmdisstrcoresurfacc, alternative=alternative, data=table))))
  mypvalue <- NA
  try(mypvalue <- pvalue(get(test)(rate ~ type | ptmdisstrcoresurfacc, alternative=alternative, data=subset(table, !(ptmdisstrcoresurfacc %in% (table %>% group_by(ptmdisstrcoresurfacc) %>% tally %>% filter(n < minsize) %>% pull(ptmdisstrcoresurfacc)))))))
  results <- rbind(results, data.frame(strat="ptm-disstr-coresurf-acc", pval=mypvalue))
  setTxtProgressBar(bar, 16)

  
  
  results
  results[results$pval==min(results$pval),]
  # str(paste(sort(as.character(results[results$pval==min(results$pval),]$strat)), collapse=","))
  
  # Store results for this file (evorate type)
  # allresults <- data.frame()
  # allresults <- rbind(allresults, data.frame(test="wilcox", evorate=evorate, min=paste(sort(as.character(results[results$pval==min(results$pval),]$strat)), collapse=","), results=results))

  #   allresults <- rbind(allresults, data.frame(test=test, evorate=evorate, results=results, min=results$pval==min(results$pval)))
  #   allresults
  
  close(bar)

  return(results)
}
# debug(GetResults)


# start

if (!exists("allresults"))
{
  # run this to clear the results dataframe:
  allresults <- data.frame()
}

# tests = c("wilcox_test", "oneway_test", "normal_test", "median_test")
# tests = c("maxstat_test")
# tests = c("fligner_test", "ansari_test", "lbl_test")
# nosurfs = c("-nosurf", "")
# nosurfs = c("", "-nosurf")
# nosurfs = c("-coresurf", "-nosurf")
# nosurfs = c("")
# tests = c("wilcox_test")
# The test functions (e.g. wilcox_test) are provided by the coin package (they have an underscore instead of a period, i.e. wilcox_test instead of wilcox.test).
# These support stratification, unlike the standard wilcox.test.
# Description from ?wilcox_test in R: "oneway_test, wilcox_test, kruskal_test, normal_test, median_test and savage_test provide the Fisher-Pitman permutation test, the Wilcoxon-Mann-Whitney test, the Kruskal-Wallis test, the van der Waerden test, the Brown-Mood median test and the Savage test".
# i.e. oneway_test is a permutation test.

# combinations <- list()
# # All
# nosurfs = c("-coresurf")
# tests = c("wilcox_test", "oneway_test")
# sources = c("Ochoa", "UniProt", "PhosphoSitePlus", "dbPTM", "all")
# evorates = c("rate4site", "lichtarge", "capra0")
# mafftmodes = c("einsi_tree_1para", "einsi_tree", "einsi", "linsi_tree", "linsi", "ginsi_tree", "ginsi")
# for (nosurf in nosurfs) {
#   for (test in tests) {
#     for (source in sources) {
#       for (evorate in evorates) {
#         for (mafftmode in mafftmodes) {
#           combinations <- c(combinations, list(c(nosurf, test, source, evorate, mafftmode)))
#         }
#       }
#     }
#   }
# }
# 
# # Finding out best evorate_mafftmode:
# # wilcox_test and Ochoa only (faster)
# nosurfs = c("-coresurf")
# tests = c("wilcox_test")
# sources = c("Ochoa")
# evorates = c("rate4site", "lichtarge", "capra0")
# mafftmodes = c("einsi_tree_1para", "einsi_tree", "einsi", "linsi_tree", "linsi", "ginsi_tree", "ginsi")
# for (nosurf in nosurfs) {
#   for (test in tests) {
#     for (source in sources) {
#       for (evorate in evorates) {
#         for (mafftmode in mafftmodes) {
#           combinations <- c(combinations, list(c(nosurf, test, source, evorate, mafftmode)))
#         }
#       }
#     }
#   }
# }
# 
# # Finding out best source and test:
# # rate4site_einsi_tree_1para only (faster)
# nosurfs = c("-coresurf")
# tests = c("wilcox_test", "oneway_test")
# sources = c("Ochoa", "UniProt", "PhosphoSitePlus", "dbPTM", "all")
# evorates = c("rate4site")
# mafftmodes = c("einsi_tree_1para")
# for (nosurf in nosurfs) {
#   for (test in tests) {
#     for (source in sources) {
#       for (evorate in evorates) {
#         for (mafftmode in mafftmodes) {
#           combinations <- c(combinations, list(c(nosurf, test, source, evorate, mafftmode)))
#         }
#       }
#     }
#   }
# }
# combinations
# combinations[1]
# combinations[[1]]

myalpha = 0.05

start1 <- proc.time()
# for (combination in combinations) {
#   nosurf <- combination[1]
#   test <- combination[2]
#   source <- combination[3]
#   evorate <- combination[4]
#   mafftmode <- combination[5]

  start2 <- proc.time()
  file <- paste("../../evolutionary_rate_analysis/tmp/tmp-dataframe-all-human-", source, "-", evorate, "_", mafftmode, "-", pred, nosurf, ".txt", sep="")

  # cat(paste("\n\nNow running '", nosurf, "' (", match(nosurf, nosurfs), "/", length(nosurfs), ") test '", test, "' (", match(test, tests), "/", length(tests), ") on source '", source, "' (", match(source, sources), "/", length(sources), ") and evorate '", evorate, "' (", match(evorate, evorates), "/", length(evorates), ") and mafftmode '", mafftmode, "' (", match(mafftmode, mafftmodes), "/", length(mafftmodes), ")\n\n", sep=""))
  cat(paste("\n\nRunning '", nosurf, "' '", test, "' (alternative '", alternative, "') on source '", source, "' and evorate '", evorate, "' and mafftmode '", mafftmode, "' and pred '", pred, "' and minsize ", minsize, " and disfilt '", disfilt, "':\n\n", sep=""))
  results <- GetResults(test, file)

  new <- data.frame(nosurf=nosurf, test=test, alternative=alternative, source=source, evorate=evorate, mafftmode=mafftmode, pred=pred, minsize=minsize, disfilt=disfilt, strat=results$strat, pval=results$pval, sig=results$pval<myalpha, min=results$pval==min(results$pval))
  allresults <- rbind(allresults, new)
  cat("\nAll results:\n\n")
  print(allresults)

#   cat("\nJust added:\n\n")
#   print(new)
#   cat("\n")

  outfile <- paste0("../output/output-", nosurf, "-", test, "-", alternative, "-", source, "-", evorate, "-", mafftmode, "-", pred, "-", minsize, "-", disfilt, ".tsv")
  write_tsv(allresults, file = outfile)
  cat(paste0("\nWrote to '", outfile, "'\n\n"))

  cat(paste("\nTime elapsed: ", round((proc.time() - start2)[[3]]), " seconds\n\n", sep=""))
# }

# print(allresults)
cat(paste("\nTotal time elapsed: ", round((proc.time() - start1)[[3]]), " seconds\n\n", sep=""))

# # Plot stratifications to find the best one overall
# tmpresults <- allresults
# # All tests: disstr, dis, ptmdisstr, ptmdis are best, and roughly equal.
# # tmpresults <- allresults[allresults$test %in% c("wilcox_test"),]
# # tmpresults <- allresults[allresults$test %in% c("normal_test"),]
# # tmpresults <- allresults[allresults$test %in% c("wilcox_test", "normal_test),]
# # tmpresults <- allresults[allresults$test %in% c("wilcox_test", "normal_test", "oneway_test", "median_test"),]
# 
# str(tmpresults)
# tmpresults
# min(tmpresults[tmpresults$pval!=0,]$pval)
# # tmpresults[tmpresults$pval == 0,]$pval <- 2.2e-16
# min(tmpresults$pval)
# head(tmpresults$pval)
# head(log(tmpresults$pval))
# 
# # Comparing p-values for tests:
# ggplot(tmpresults, aes(x=test, y=pval)) + geom_boxplot() + scale_y_log10()
# # >> wilcox_test and normal_test give "best" results
# 
# # Comparing stratifications:
# ggplot(tmpresults, aes(x=strat, y=pval)) + geom_boxplot() + scale_y_log10() + geom_point()
# ggplot(tmpresults, aes(x=strat, y=pval)) + geom_boxplot() + scale_y_log10()

# # All tests: disstr, dis, ptmdisstr, ptmdis are best, and roughly equal.
# # wilcox_test and normal_test: ptmdisstr, ptmdis, disstr, dis
# # wilcox_test: disstr, dis
# # normal_test: disstr, dis


# # Reformat allresults table
# library(reshape2)
# tmp <- melt(allresults, c("nosurf", "test", "source", "evorate", "mafftmode", "strat"), c("pval", "sig", "min"))
# tmp
# output <- dcast(tmp, nosurf + source + evorate + mafftmode + strat ~ test)
# output

# # Write all results to file
# # write.table(allresults, file = paste("output-allresults.txt", sep = ""), row.accs = FALSE, col.accs = TRUE, quote = FALSE, sep = "\t")
# outfile <- paste0("../output/output-", nosurf, "-", test, "-", source, "-", evorate, "-", mafftmode, ".tsv")
# write_tsv(allresults, file = outfile)
