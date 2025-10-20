# install.packages("Matching")
evoratedata <- read.delim("~/Desktop/tmp/tmp-dataframe-all-human-lichtarge_einsi_tree-MULTICOM.txt", header=T, quote="")

# evoratedata$name <- NULL
evoratedata$site <- NULL

# evoratedata <- evoratedata[evoratedata$dis != "A", ]
levels(evoratedata$type) <- list("Control"="C", "PTM"="P")
# evoratedata$rate <- log2(evoratedata$rate)
# levels(evoratedata$dis) <- list("Core / Structured"="strcore", "Core / Disordered"="discore", "Surface / Structured"="strsurf", "Surface / Disordered"="dissurf", "Combined"="A")
# levels(evoratedata$dis) <- list("All sites"="A", "Core / Structured"="strcore", "Core / Disordered"="discore", "Surface / Structured"="strsurf", "Surface / Disordered"="dissurf")
levels(evoratedata$dis) <- list("All sites"="A", "Core / Disordered"="discore", "Core / Structured"="strcore", "Surface / Disordered"="dissurf", "Surface / Structured"="strsurf")

# Remove Core / Disordered since it doesnt occur often enough?
# evoratedata <- evoratedata[evoratedata$dis != "discore", ]
# levels(evoratedata$dis) <- list("Core / Structured"="strcore", "Surface / Structured"="strsurf", "Surface / Disordered"="dissurf", "Combined"="A")
# levels(evoratedata$dis) <- list("All sites"="A", "Core / Structured"="strcore", "Surface / Structured"="strsurf", "Surface / Disordered"="dissurf")

library(Matching)		# for ks.boot test (Kolmogorov-Smirnov with 1000x bootstrapping, allows for ties in contrast to regular KS test)

pvalues <- data.frame()

str(evoratedata$ptm)
head(evoratedata$ptm)

str(evoratedata$dis)
head(evoratedata$dis)

str(evoratedata$type)
head(evoratedata$type)


teststatistic <- data.frame(ptm=character(), dis=character(), i=numeric(), p=numeric(), nptm=numeric(), ncontrol=numeric(), w=numeric())



for (x in unique(evoratedata$ptm)) {
  print(paste("x is now:", x))
  for (y in sort(unique(evoratedata$dis))) {
    print(paste("y is now:", y))
    
    # ptmr <- subset(evoratedata, ptm == x & dis == y & evoratedata$type == "PTM")$rate
    # controlr <- subset(evoratedata, ptm == x & dis == y & evoratedata$type == "Control")$rate
    
    ptmr <- subset(evoratedata, ptm == x & dis == y & evoratedata$type == "PTM")
    controlr <- subset(evoratedata, ptm == x & dis == y & evoratedata$type == "Control")
    
    # ptmr <- subset(evoratedata, ptm == x & evoratedata$type == "PTM")
    # controlr <- subset(evoratedata, ptm == x & evoratedata$type == "Control")
    
    wtmp <- wilcox.test(ptmr$rate, controlr$rate, alternative="two.sided", paired=FALSE)
    ttmp <- t.test(ptmr$rate, controlr$rate, alternative="two.sided", paired=FALSE)
    
    teststatistic <- rbind(teststatistic, data.frame(ptm=x, dis=y, i=1, p=wtmp$p.value, nptm=length(ptmr$rate), ncontrol=length(controlr$rate), w=wtmp$statistic))
    
    pvalues <- rbind(pvalues, data.frame(ptm=x, dis=y, nptm=length(ptmr$rate), ncontrol=length(controlr$rate), medianptm=median(ptmr$rate), mediancontrol=median(controlr$rate), meanptm=mean(ptmr$rate), meancontrol=mean(controlr$rate), pmw=wtmp$p.value, ptt=ttmp$p.value, mediandif=round(median(ptmr$rate)-median(controlr$rate), 3), meandif=round(mean(ptmr$rate)-mean(controlr$rate), 3)))
  }
}

write.table(pvalues, "~/Desktop/tmp/output-pvalues-dataframe-all-human-lichtarge_einsi_tree-MULTICOM.txt", sep="\t", quote=F, row.names=F)
