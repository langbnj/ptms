blang_init()

library(viridisLite)

# Minimum number of sites for a PTM type to be included
minsites <- 1000





# # Analyses
# 
# # # Allele frequency distribution in gnomAD
# # q_gnomad <- Query("SELECT ac, COUNT(*) AS c FROM snps_gnomad WHERE source IN ('v4_wes', 'v4_wgs') GROUP BY ac ORDER BY ac")
# # q_gnomad %>%
# #   mutate(acmax = ifelse(ac <= 10, ac, 10)) %>%
# #   group_by(acmax) %>%
# #   summarise(n = sum(c)) %>%
# #   ggplot(aes(x = acmax, y = n)) +
# #   geom_bar(stat = "identity") +
# #   scale_x_continuous(breaks = seq(0, 10, 1)) +
# #   theme_nature() +
# #   ggtitle("COSMIC")
# # qsave("output-af-cosmic.pdf")
# 
# # Sample count distribution in COSMIC
# qcosmic <- Query("SELECT ac, COUNT(*) AS c FROM snps_cosmic WHERE cell_line=0 GROUP BY ac ORDER BY ac")
# qcosmic %>%
#   mutate(acmax = ifelse(ac <= 10, ac, 10)) %>%
#   group_by(acmax) %>%
#   summarise(n = sum(c)) %>%
#   ggplot(aes(x = acmax, y = n)) +
#   geom_bar(stat = "identity") +
#   scale_x_continuous(breaks = seq(0, 10, 1)) +
#   theme_nature() +
#   ggtitle("COSMIC")
# qsave("output-ac-cosmic.pdf")
# 
# # Papers count distribution in Mastermind
# qmas <- Query("SELECT papers AS ac, COUNT(*) AS c FROM snps_mastermind GROUP BY papers ORDER BY papers")
# qmas %>%
#   mutate(acmax = ifelse(ac <= 10, ac, 10)) %>%
#   group_by(acmax) %>%
#   summarise(n = sum(c)) %>%
#   ggplot(aes(x = acmax, y = n)) +
#   geom_bar(stat = "identity") +
#   scale_x_continuous(breaks = seq(0, 10, 1)) +
#   theme_nature() +
#   ggtitle("Genomenon Mastermind")
# qsave("output-ac-mastermind.pdf")
# 
# # Best papers threshold for Genomenon Mastermind?
# snps_mastermind
# snps_mastermind %>% nrow
# snps_mastermind %>% group_by(acc, site, original, variant) %>% tally %>% nrow
# snps_mastermind %>% select(acc, site, variant) %>% unique %>% nrow
# snps_mastermind %>% select(acc, site, original, variant) %>% unique %>% nrow
# # >> 10,457,743
# snps_mastermind %>% filter(ac >= 1) %>% tally / nrow(snps_mastermind)
# snps_mastermind %>% filter(ac >= 2) %>% tally / nrow(snps_mastermind)
# snps_mastermind %>% filter(ac >= 3) %>% tally / nrow(snps_mastermind)
# snps_mastermind %>% filter(ac >= 4) %>% tally / nrow(snps_mastermind)
# snps_mastermind %>% filter(ac >= 5) %>% tally / nrow(snps_mastermind)












# plddt_filtering <- 0 # Disable
plddt_filtering <- "70_PAE2"
analyse_contingency_coresurf <- function(plddt_filtering = 0) {

  # "Figures are best prepared at a width of 90 mm (single column) and 180 mm (double column) with a maximum height of 170mm. At this size, the font size should be 5-7pt." (https://www.nature.com/nature/for-authors/initial-submission)
  # Absolute max Nature dimensions 183 mm x 247 mm ("For guidance, Natures standard figure sizes are 89 mm wide (single column) and 183 mm wide (double column). The full depth of a Nature page is 247 mm. Figures can also be a column-and-a-half where necessary (120–136 mm).", https://www.nature.com/nature/for-authors/final-submission)
  
  
  
  
  # Load data
  
  # # As in ~/pipeline/evolutionary_rate_analysis (acc|sites with evolutionary rates only):
  # # Get PTM sites
  # qpr <- Query("SELECT m.ptm, m.acc, r.ensp, m.site, SUBSTRING(m.ptm, 1, 1) AS aa, 'ptm' AS type, SUBSTRING(s.seq, m.site, 1) AS aatest, SUBSTRING(d.seq, m.site, 1) AS dis, SUBSTRING(cs.seq, m.site, 1) AS coresurf, MIN(r.rate) AS minrate, COUNT(DISTINCT r.rate) AS rates FROM unimod m, uniseq s, uniseq cs, uniseq d, uniens ue, evorate_rate4site_einsi_tree_1para r WHERE m.species='human' AND m.ptm IS NOT NULL AND m.acc=s.acc AND m.acc=cs.acc AND m.acc=d.acc AND s.type IN ('UniProt', 'UniIso') AND cs.type='CoreSurf' AND SUBSTRING(cs.seq, m.site, 1) IN ('C', 'S') AND d.type='AlphaFold' AND ue.acc=m.acc AND ue.ensp=r.ensp AND m.site=r.site GROUP BY m.acc, m.site, m.ptm ORDER BY m.acc, m.site, m.ptm")
  # # Get Control sites
  # qcr <- Query("SELECT t.ptm, m.acc, r.ensp, m.site, m.aa, 'control' AS type, SUBSTRING(s.seq, m.site, 1) AS aatest, SUBSTRING(d.seq, m.site, 1) AS dis, SUBSTRING(cs.seq, m.site, 1) AS coresurf, MIN(r.rate) AS minrate, COUNT(DISTINCT r.rate) AS rates FROM unimod_control m, uniseq s, uniseq cs, uniseq d, uniens ue, evorate_rate4site_einsi_tree_1para r, (SELECT m.ptm, m.acc FROM unimod m, uniseq s, uniseq cs, uniseq d, uniens ue, evorate_rate4site_einsi_tree_1para r WHERE m.species='human' AND m.ptm IS NOT NULL AND m.acc=s.acc AND m.acc=cs.acc AND m.acc=d.acc AND s.type IN ('UniProt', 'UniIso') AND cs.type='CoreSurf' AND SUBSTRING(cs.seq, m.site, 1) IN ('C', 'S') AND d.type='AlphaFold' AND ue.acc=m.acc AND ue.ensp=r.ensp AND m.site=r.site GROUP BY m.ptm, m.acc ORDER BY m.ptm, m.acc) t WHERE m.acc=t.acc AND m.aa=SUBSTRING(t.ptm, 1, 1) AND m.species='human' AND m.acc=s.acc AND m.acc=cs.acc AND m.acc=d.acc AND s.type IN ('UniProt', 'UniIso') AND cs.type='CoreSurf' AND SUBSTRING(cs.seq, m.site, 1) IN ('C', 'S') AND d.type='AlphaFold' AND ue.acc=m.acc AND ue.ensp=r.ensp AND m.site=r.site GROUP BY m.acc, m.site, t.ptm, m.aa ORDER BY m.acc, m.site, t.ptm, m.aa")
  # >> ...but why care about evorates here? They aren't relevant.
  # >> Can use the above qpr and qcr queries for the evorate plots, though, and for the combined evorate/af "sail" plots.
  
  # # Verify the number of evolutionary rates per residue is as expected (== 1)
  # # stopifnot(q %>% filter(rates != 1) %>% nrow == 0)
  # q %>% filter(rates != 1) %>% pull(rates) %>% as.factor %>% summary
  # q %>% pull(rates) %>% as.factor %>% summary
  # # What percentage of sites have more than one rate?
  # q %>% filter(rates != 1) %>% nrow / nrow(q) * 100
  # # >> Only 0.26% of sites (one in 400) have more than one rate. This is negligible, so I can ignore it. I'll use the minimum rate below.
  # # Check for duplicate rows by ENSP
  # stopifnot(q %>% group_by(ptm, acc, site) %>% nrow == q %>% group_by(ptm, ensp, site) %>% nrow)
  # # Clean up
  # q %<>% mutate(rate = minrate) %>% select(-rates, -minrate) %>% arrange(acc, site, ptm, type)
  
  # Without requiring evorates (we don't need them in this analysis here):
  
  # Run pLDDT filtering?
  tmp_plddt <- ""
  tmp_uniseq_suffix = ""
  if (plddt_filtering != 0) {
    tmp_plddt <- f("-pLDDT{plddt_filtering}")
    tmp_uniseq_suffix = f("_pLDDT{plddt_filtering}")
  }
  
  
  
  # Get PTM sites
  tmpfile <- f("tmp-q{tmp_plddt}.rds")
  if (!file.exists(tmpfile)) {
    # with uniiso
    # qp <- Query("SELECT m.ptm, m.acc, m.site, SUBSTRING(m.ptm, 1, 1) AS aa, 'ptm' AS type, SUBSTRING(s.seq, m.site, 1) AS aatest, SUBSTRING(d.seq, m.site, 1) AS disstr, SUBSTRING(cs.seq, m.site, 1) AS coresurf FROM unimod m, uniseq s, uniseq cs, uniseq d WHERE m.species='human' AND m.ptm IS NOT NULL AND m.acc=s.acc AND m.acc=cs.acc AND m.acc=d.acc AND s.type IN ('UniProt', 'UniIso') AND cs.type='CoreSurf' AND SUBSTRING(cs.seq, m.site, 1) IN ('C', 'S') AND d.type='AlphaFold' GROUP BY m.acc, m.site, m.ptm ORDER BY m.acc, m.site, m.ptm")
    # no UniIso
    # qp <- Query("SELECT m.ptm, m.acc, m.site, SUBSTRING(m.ptm, 1, 1) AS aa, 'ptm' AS type, SUBSTRING(s.seq, m.site, 1) AS aatest, SUBSTRING(d.seq, m.site, 1) AS disstr, SUBSTRING(cs.seq, m.site, 1) AS coresurf FROM unimod m, uniseq s, uniseq cs, uniseq d WHERE m.species='human' AND m.ptm IS NOT NULL AND m.acc=s.acc AND m.acc=cs.acc AND m.acc=d.acc AND s.type IN ('UniProt') AND LENGTH(s.seq)>=16 AND cs.type='CoreSurf' AND SUBSTRING(cs.seq, m.site, 1) IN ('C', 'S') AND d.type='AlphaFold' GROUP BY m.acc, m.site, m.ptm ORDER BY m.acc, m.site, m.ptm")
    qp <- Query(f("SELECT m.ptm, m.acc, m.site, SUBSTRING(m.ptm, 1, 1) AS aa, 'ptm' AS type, SUBSTRING(s.seq, m.site, 1) AS aatest, SUBSTRING(d.seq, m.site, 1) AS disstr, SUBSTRING(cs.seq, m.site, 1) AS coresurf FROM unimod m, uniseq s, uniseq cs, uniseq d WHERE m.species='human' AND m.ptm IS NOT NULL AND m.acc=s.acc AND m.acc=cs.acc AND m.acc=d.acc AND s.type IN ('UniProt') AND LENGTH(s.seq)>=16 AND cs.type='CoreSurf{tmp_uniseq_suffix}' AND SUBSTRING(cs.seq, m.site, 1) IN ('C', 'S') AND d.type='AlphaFold{tmp_uniseq_suffix}' GROUP BY m.acc, m.site, m.ptm ORDER BY m.acc, m.site, m.ptm"))
    
    # Get Control sites
    # with uniiso
    # qc <- Query("SELECT t.ptm, m.acc, m.site, m.aa, 'control' AS type, SUBSTRING(s.seq, m.site, 1) AS aatest, SUBSTRING(d.seq, m.site, 1) AS disstr, SUBSTRING(cs.seq, m.site, 1) AS coresurf FROM unimod_control m, uniseq s, uniseq cs, uniseq d, (SELECT m.ptm, m.acc FROM unimod m, uniseq s, uniseq cs, uniseq d WHERE m.species='human' AND m.ptm IS NOT NULL AND m.acc=s.acc AND m.acc=cs.acc AND m.acc=d.acc AND s.type IN ('UniProt', 'UniIso') AND cs.type='CoreSurf' AND SUBSTRING(cs.seq, m.site, 1) IN ('C', 'S') AND d.type='AlphaFold' GROUP BY m.ptm, m.acc ORDER BY m.ptm, m.acc) t WHERE m.acc=t.acc AND m.aa=SUBSTRING(t.ptm, 1, 1) AND m.species='human' AND m.acc=s.acc AND m.acc=cs.acc AND m.acc=d.acc AND s.type IN ('UniProt', 'UniIso') AND cs.type='CoreSurf' AND SUBSTRING(cs.seq, m.site, 1) IN ('C', 'S') AND d.type='AlphaFold' GROUP BY m.acc, m.site, t.ptm, m.aa ORDER BY m.acc, m.site, t.ptm, m.aa")
    # >> Around 20% more sites to analyse this way [but they might be redundant]
    # without uniiso
    # qc <- Query("SELECT t.ptm, m.acc, m.site, m.aa, 'control' AS type, SUBSTRING(s.seq, m.site, 1) AS aatest, SUBSTRING(d.seq, m.site, 1) AS disstr, SUBSTRING(cs.seq, m.site, 1) AS coresurf FROM unimod_control m, uniseq s, uniseq cs, uniseq d, (SELECT m.ptm, m.acc FROM unimod m, uniseq s, uniseq cs, uniseq d WHERE m.species='human' AND m.ptm IS NOT NULL AND m.acc=s.acc AND m.acc=cs.acc AND m.acc=d.acc AND s.type IN ('UniProt') AND LENGTH(s.seq)>=16 AND cs.type='CoreSurf' AND SUBSTRING(cs.seq, m.site, 1) IN ('C', 'S') AND d.type='AlphaFold' GROUP BY m.ptm, m.acc ORDER BY m.ptm, m.acc) t WHERE m.acc=t.acc AND m.aa=SUBSTRING(t.ptm, 1, 1) AND m.species='human' AND m.acc=s.acc AND m.acc=cs.acc AND m.acc=d.acc AND s.type IN ('UniProt') AND cs.type='CoreSurf' AND SUBSTRING(cs.seq, m.site, 1) IN ('C', 'S') AND d.type='AlphaFold' GROUP BY m.acc, m.site, t.ptm, m.aa ORDER BY m.acc, m.site, t.ptm, m.aa")
    qc <- Query(f("SELECT t.ptm, m.acc, m.site, m.aa, 'control' AS type, SUBSTRING(s.seq, m.site, 1) AS aatest, SUBSTRING(d.seq, m.site, 1) AS disstr, SUBSTRING(cs.seq, m.site, 1) AS coresurf FROM unimod_control m, uniseq s, uniseq cs, uniseq d, (SELECT m.ptm, m.acc FROM unimod m, uniseq s, uniseq cs, uniseq d WHERE m.species='human' AND m.ptm IS NOT NULL AND m.acc=s.acc AND m.acc=cs.acc AND m.acc=d.acc AND s.type IN ('UniProt') AND LENGTH(s.seq)>=16 AND cs.type='CoreSurf{tmp_uniseq_suffix}' AND SUBSTRING(cs.seq, m.site, 1) IN ('C', 'S') AND d.type='AlphaFold{tmp_uniseq_suffix}' GROUP BY m.ptm, m.acc ORDER BY m.ptm, m.acc) t WHERE m.acc=t.acc AND m.aa=SUBSTRING(t.ptm, 1, 1) AND m.species='human' AND m.acc=s.acc AND m.acc=cs.acc AND m.acc=d.acc AND s.type IN ('UniProt') AND cs.type='CoreSurf{tmp_uniseq_suffix}' AND SUBSTRING(cs.seq, m.site, 1) IN ('C', 'S') AND d.type='AlphaFold{tmp_uniseq_suffix}' GROUP BY m.acc, m.site, t.ptm, m.aa ORDER BY m.acc, m.site, t.ptm, m.aa"))
    
    # #DEBUG test updated uniseq
    # qp <- Query("SELECT m.ptm, m.acc, m.site, SUBSTRING(m.ptm, 1, 1) AS aa, 'ptm' AS type, SUBSTRING(s.seq, m.site, 1) AS aatest, SUBSTRING(d.seq, m.site, 1) AS disstr, SUBSTRING(cs.seq, m.site, 1) AS coresurf FROM unimod m, uniseq_backup_updated_with_alphasync_mappings_but_not_used s, uniseq_backup_updated_with_alphasync_mappings_but_not_used cs, uniseq_backup_updated_with_alphasync_mappings_but_not_used d WHERE m.species='human' AND m.ptm IS NOT NULL AND m.acc=s.acc AND m.acc=cs.acc AND m.acc=d.acc AND s.type IN ('UniProt') AND cs.type='CoreSurf' AND SUBSTRING(cs.seq, m.site, 1) IN ('C', 'S') AND d.type='AlphaFold' GROUP BY m.acc, m.site, m.ptm ORDER BY m.acc, m.site, m.ptm")
    # qc <- Query("SELECT t.ptm, m.acc, m.site, m.aa, 'control' AS type, SUBSTRING(s.seq, m.site, 1) AS aatest, SUBSTRING(d.seq, m.site, 1) AS disstr, SUBSTRING(cs.seq, m.site, 1) AS coresurf FROM unimod_control m, uniseq_backup_updated_with_alphasync_mappings_but_not_used s, uniseq_backup_updated_with_alphasync_mappings_but_not_used cs, uniseq_backup_updated_with_alphasync_mappings_but_not_used d, (SELECT m.ptm, m.acc FROM unimod m, uniseq_backup_updated_with_alphasync_mappings_but_not_used s, uniseq_backup_updated_with_alphasync_mappings_but_not_used cs, uniseq_backup_updated_with_alphasync_mappings_but_not_used d WHERE m.species='human' AND m.ptm IS NOT NULL AND m.acc=s.acc AND m.acc=cs.acc AND m.acc=d.acc AND s.type IN ('UniProt') AND cs.type='CoreSurf' AND SUBSTRING(cs.seq, m.site, 1) IN ('C', 'S') AND d.type='AlphaFold' GROUP BY m.ptm, m.acc ORDER BY m.ptm, m.acc) t WHERE m.acc=t.acc AND m.aa=SUBSTRING(t.ptm, 1, 1) AND m.species='human' AND m.acc=s.acc AND m.acc=cs.acc AND m.acc=d.acc AND s.type IN ('UniProt') AND cs.type='CoreSurf' AND SUBSTRING(cs.seq, m.site, 1) IN ('C', 'S') AND d.type='AlphaFold' GROUP BY m.acc, m.site, t.ptm, m.aa ORDER BY m.acc, m.site, t.ptm, m.aa")
    # #END DEBUG
    
    # qpr
    qp
    # qcr
    qc
    
    # Analysis:
    qc %>% group_by(aa) %>% summarise(accs = n_distinct(acc))
    qp %>% group_by(aa) %>% summarise(accs = n_distinct(acc))
    
    qc %>% group_by(aa) %>% summarise(accs = n_distinct(acc)) %>% summarise(aa_accs = sum(accs))
    qp %>% group_by(aa) %>% summarise(accs = n_distinct(acc)) %>% summarise(aa_accs = sum(accs))
    # >> qp has quite a few more aa|acc combinations

    # Limit to ptm|accs that have both PTM and control residues
    # This ensures proteins are modifiable with a given PTM
    # i.e. skip any ptm|acc combinations that do not have control residues
    qp_accs <- qp %>% select(acc) %>% unique
    qp_accs %>% nrow
    qc_accs <- qc %>% select(acc) %>% unique
    qc_accs %>% nrow
    # How many mod accs are in control?
    qp_accs %>% mutate(in_qc = acc %in% qc_accs$acc) %>% group_by(in_qc) %>% tally
    # >> 395 accs only have mods, but no controls
    # How many control accs are in mod?
    qc_accs %>% mutate(in_qp = acc %in% qp_accs$acc) %>% group_by(in_qp) %>% tally
    # >> All control accs are in mod (good!)
    qp$acc %>% unique %>% length
    qp %>% select(acc, site) %>% unique %>% nrow
    qp %>% select(acc)
    qc %>% select(acc)
    sort(intersect(qp$acc, qc$acc)) %>% head
    qc %>% select(ptm, acc) %>% unique %>% nrow
    qp %>% select(ptm, acc) %>% unique %>% nrow
    qp %>% select(ptm, acc) %>% unique
    qc
    qp
    qp0 <- qp
    qp %<>% inner_join(qc %>% select(ptm, acc) %>% distinct, by = c("ptm", "acc"))
    qc %>% select(ptm, acc) %>% unique %>% nrow
    qp %>% select(ptm, acc) %>% unique %>% nrow
    qc %>% select(ptm, acc) %>% distinct %>% arrange(ptm, acc)
    qp %>% select(ptm, acc) %>% distinct %>% arrange(ptm, acc)
    identical(qc %>% select(ptm, acc) %>% distinct %>% arrange(ptm, acc), qp %>% select(ptm, acc) %>% distinct %>% arrange(ptm, acc))
    qc %>% select(acc) %>% unique %>% nrow
    qp %>% select(acc) %>% unique %>% nrow
    identical(qc %>% select(acc) %>% distinct, qp %>% select(acc) %>% distinct)
    
    # Combine into one df
    q <- bind_rows(qp, qc)
    q
    # Verify the amino acid is as expected
    stopifnot(q %>% filter(aa != aatest) %>% nrow == 0)
    q$aatest <- NULL
    q
    # Check for duplicate rows
    stopifnot(q %>% group_by(ptm, acc, site) %>% nrow == q %>% nrow)
    # Define structural categories (dis)
    q$dis <- NA
    q %<>% mutate(dis = ifelse(coresurf == "C", "Buried", dis))
    # q %<>% mutate(dis = ifelse(coresurf == "S", "Surface", dis))
    q %<>% mutate(dis = ifelse(coresurf == "S", "Structured", dis))
    q %<>% mutate(dis = ifelse(coresurf == "S" & disstr == "*", "Disordered", dis))
    # q$dis <- fct_relevel(q$dis, "Buried", "Surface", "Disordered")
    q$dis <- fct_relevel(q$dis, "Buried", "Structured", "Disordered")
    q %>% group_by(dis, disstr, coresurf) %>% tally
    q %>% group_by(dis) %>% tally
    q
    # Convert type to factor
    q %<>% mutate(type = as.factor(type))
    # Clean up
    q %<>% select(-disstr, -coresurf)
    q %<>% arrange(acc, site, ptm, type)
    q
    
    write_rds(q, tmpfile)
  } else {
    q <- read_rds(tmpfile)
  }
  
  # pLDDT filtering: remove disordered sites from analysis (too few)
  if (plddt_filtering != 0) {
    q$dis <- as.character(q$dis)
    q %<>% filter(dis != "Disordered")
    q$dis <- factor(q$dis, levels = c("Buried", "Structured"))
  }

  # Store original q in q0
  q0 <- q
  
  
  
  
  # Get accs from uniens to make sure they are mapped to Ensembl 108 (for VEP, which is used for all of these variant sets)
  # uniens_accs <- Query("SELECT DISTINCT ue.acc FROM uniens ue, uniseq s WHERE ue.species='human' AND ue.acc=s.acc AND s.type='UniProt' AND LENGTH(s.seq)>=16") %>% pull(acc)
  uniens_accs <- read_tsv("uniens_accs.tsv.gz") %>% pull(acc)
  uniens_accs %>% length
  
  # Get COSMIC
  tmpfile <- f("tmp-cosmic.rds")
  if (!file.exists(tmpfile)) {
    # Note: no cell lines (cell_line=0)
    # snps_cosmic <- Query("SELECT chr, pos, originalbase, variantbase, source, acc, site, original, variant, ac, COUNT(DISTINCT ac) AS acs FROM snps_cosmic WHERE cell_line=0 GROUP BY chr, pos, originalbase, variantbase, source, acc, site, original, variant")
    # Condense genomic to protein variants (using sum of sample counts (ac) as a very rough approximation of allele frequency (not plotting this, only plotting presence/absence of variants for an acc|site)
    snps_cosmic <- Query("SELECT acc, site, original, variant, SUM(ac) AS ac, cgc FROM snps_cosmic WHERE acc=canon AND cell_line=0 GROUP BY acc, site, original, variant")
    # # Check uniqueness of rows (no other purpose)
    # stopifnot(snps_cosmic %>% filter(acs != 1) %>% nrow == 0)
    # Drop sources & unique
    # snps_cosmic %<>% select(-source, -acs) %>% unique
    write_rds(snps_cosmic, tmpfile)
  } else {
    snps_cosmic <- read_rds(tmpfile)
  }
  snps_cosmic
  
  # Get ClinVar
  tmpfile <- f("tmp-clinvar.rds")
  if (!file.exists(tmpfile)) {
    # Condense genomic to protein variants (using sum of submitters (ac) as a very rough approximation of allele frequency (not plotting this, only plotting presence/absence of variants for an acc|site)
    snps_clinvar <- Query("SELECT acc, site, original, variant, MAX(clinsig) AS clinsig, SUM(submitters) AS ac FROM snps_clinvar WHERE acc=canon GROUP BY acc, site, original, variant")
    # stopifnot(snps_clinvar %>% filter(clinsigs != 1) %>% nrow == 0)
    # This happens sometimes, that there is one site that has the same variant as clinsig=0 and as clinsig=1. Obviously I should consider it clinsig=1 in that case, hence MAX() above.
    # SELECT acc, site, original, variant, clinsig, COUNT(DISTINCT clinsig) AS clinsigs FROM snps_clinvar GROUP BY acc, site, original, variant HAVING clinsigs!=1;
    # SELECT * FROM snps_clinvar WHERE acc='A8MTJ6' AND site=234 AND variant='L';
    # snps_clinvar %<>% select(-clinsigs) %>% unique
    write_rds(snps_clinvar, tmpfile)
  } else {
    snps_clinvar <- read_rds(tmpfile)
  }
  snps_clinvar
  
  # Get Genomenon Mastermind
  tmpfile <- f("tmp-mastermind.rds")
  if (!file.exists(tmpfile)) {
    # snps_mastermind <- Query("SELECT chr, pos, originalbase, variantbase, acc, site, original, variant, papers AS ac, COUNT(DISTINCT papers) AS papers FROM snps_mastermind GROUP BY chr, pos, originalbase, variantbase, acc, site, original, variant")
    # Condense genomic to protein variants (using sum of papers (ac) as a very rough approximation of allele frequency (not plotting this, only plotting presence/absence of variants for an acc|site)
    snps_mastermind <- Query("SELECT acc, site, original, variant, SUM(papers) AS ac FROM snps_mastermind WHERE acc=canon GROUP BY acc, site, original, variant")
    # stopifnot(snps_mastermind %>% filter(papers != 1) %>% nrow == 0)
    # Drop sources & unique
    # snps_mastermind %<>% select(-papers) %>% unique
    write_rds(snps_mastermind, tmpfile)
  } else {
    snps_mastermind <- read_rds(tmpfile)
  }
  snps_mastermind
  
  # # Get PCGP
  # # snps_pcgp <- Query("SELECT acc, site, original, variant, patients AS ac, COUNT(DISTINCT patients) AS acs, samples, status FROM snps_pcgp GROUP BY acc, site, original, variant")
  # snps_pcgp <- Query("SELECT acc, site, original, variant, SUM(patients) AS ac, COUNT(DISTINCT patients) AS acs, samples, status FROM snps_pcgp GROUP BY acc, site, original, variant")
  # # snps_pcgp %>% filter(acs != 1)
  # # snps_pcgp %>% filter(acc == "P01111" & site == 12 & original == "G" & variant == "S")
  # # stopifnot(snps_pcgp %>% filter(acs != 1) %>% nrow == 0)
  # # This happens sometimes, that there is one site that has the same amino acid variant present multiple times (but only for 27 rows). Using SUM(patients) to address this (different base-level mutants leading to same aa change, I guess).
  # # Drop sources & unique
  # snps_pcgp %<>% select(-acs) %>% unique
  # snps_pcgp
  # 
  # # Get PCGP status='valid' (not 'putative')
  # snps_pcgp_valid <- Query("SELECT acc, site, original, variant, patients AS ac, COUNT(DISTINCT patients) AS acs, samples, status FROM snps_pcgp WHERE status='valid' GROUP BY acc, site, original, variant")
  # snps_pcgp_valid %>% filter(acs != 1)
  # snps_pcgp_valid %>% filter(acc == "P01111" & site == 12 & original == "G" & variant == "S")
  # stopifnot(snps_pcgp_valid %>% filter(acs != 1) %>% nrow == 0)
  # # Drop sources & unique
  # snps_pcgp_valid %<>% select(-acs) %>% unique
  # snps_pcgp_valid
  
  
  # Get gnomAD MAFs
  tmpfile <- f("tmp-gnomad.rds")
  if (!file.exists(tmpfile)) {
    # Query("SELECT 2+2")
    # Query("SELECT 2+a")
    # snps_gnomad <- Query("SELECT acc, site, original, variant, af FROM snps_gnomad")
    # snps_gnomad <- Query("SELECT chr, pos, originalbase, variantbase, source, acc, site, original, variant, af FROM snps_gnomad GROUP BY chr, pos, originalbase, variantbase, source, acc, site, original, variant")
    # GROUP BY and MAX(af) AS af to collapse duplicates introduced by v2 liftover from GRCh37 to GRCh38 (only affects v2):
    # After liftover, some GRCh37 repeat regions seem to have gotten condensed into single GRCh38 loci.
    # >> Using MAX(af) to collapse these (summing would lead to af>1 in one case).
    # >> Actually simply using v4, as this includes the v2 and v3 individuals.
    # Still using GROUP BY to avoid ENST duplicates:
    snps_gnomad <- Query("SELECT chr, pos, originalbase, variantbase, source, acc, site, original, variant, ac, af, COUNT(DISTINCT af) AS afs FROM snps_gnomad WHERE source IN ('v4_wes', 'v4_wgs') AND acc=canon GROUP BY chr, pos, originalbase, variantbase, source, acc, site, original, variant")
    
    # Verify that v4_wgs and v4_wes have identical AFs (AF_joint)
    snps_gnomad
    snps_gnomad %>% group_by(source) %>% tally
    qs <- full_join(snps_gnomad %>% filter(source=="v4_wgs"), snps_gnomad %>% filter(source=="v4_wes"), by = join_by(chr, pos, originalbase, variantbase, acc, site, original, variant), suffix = c("_wgs", "_wes"))
    qs
    stopifnot(qs %>% filter(afs_wgs!=1 | afs_wes!=1) %>% nrow == 0)
    stopifnot(qs %>% filter(af_wgs != af_wes) %>% nrow == 0)
    summary(qs %>% select(af_wgs, af_wes))
    # Drop sources & unique
    snps_gnomad %<>% select(-source, -afs) %>% unique
    write_rds(snps_gnomad, tmpfile)
  } else {
    snps_gnomad <- read_rds(tmpfile)
  }
  snps_gnomad
  
  
  
  
  
  
  
  
  
  
  
  
  
  # gnomAD: Condense genomic to protein variants (using sum of MAFs since these variants are very likely to be independent, so SUM(af) should be a good approximation of 1-majorAF)
  snps_gnomad
  snps_gnomad_accsite <- snps_gnomad %>% group_by(acc, site) %>% summarise(sum_af = sum(af), min_af = min(af), max_af = max(af), afs = n()) %>% mutate(af = sum_af)
  snps_gnomad_accsite
  snps_gnomad_accsite %<>% select(-min_af, -max_af, -sum_af, -afs)
  snps_gnomad_accsite
  # Replace >1 with 1
  # snps_gnomad_accsite %>% group_by(af > 1) %>% summary
  snps_gnomad_accsite %>% filter(af > 1) %>% arrange(desc(af))
  # >> Only 42 acc|sites have SUM(af) > 1. That's negligible.
  # snps_gnomad_accsite %>% filter(acc == "P20591" & site == 379)
  # SELECT * FROM snps_gnomad WHERE acc='Q7Z3E5' AND site='180' AND source IN ('v4_wgs', 'v4_wes');
  # >> Two SNVs with AF ~ 1, one from v4_wgs, one from v4_wes (hence not merged into a larger variant).
  # >> Can simply set these to AF = 1.
  snps_gnomad_accsite %<>% mutate(af = ifelse(af > 1, 1, af))
  snps_gnomad_accsite
  snps_gnomad_accsite %<>% mutate(source = "gnomAD")
  
  # The other SNP sources are less fancy (no AF):
  snps_clinvar
  snps_cosmic
  snps_mastermind
  # >> Will plot these separately below.
  
  # Combine PTM data with gnomAD accsite (amino acid level)
  q <- left_join(q, snps_gnomad_accsite)
  q %>% pull(af) %>% summary
  q %>% pull(af) %>% is.na %>% summary
  # Percentage that has af=NA (i.e. no variants exist for these acc|sites)
  q %>% filter(is.na(af)) %>% nrow / nrow(q) * 100
  # >> 52.6% of residues don't have any variants.
  # Replace NA with 0
  q %<>% mutate(af = ifelse(is.na(af), 0, af))
  q$af %>% summary
  q
  
  # # Analyse
  # # Group by ptm and type and summarise af into lists
  # q %>% group_by(ptm, dis)
  # q %>% group_by(ptm) %>% tally
  # q %>% group_by(ptm, dis) %>% tally
  # q %>% group_by(ptm, dis, type) %>% tally
  # q %>% group_by(ptm) %>% do(wilcox_test = wilcox.test(af ~ type, data = .)) %>% summarise(ptm, p_value = wilcox_test$p.value)
  # # Identify PTMs where not all 3 dis types are present
  # q %>% group_by(ptm) %>% filter(n_distinct(dis) != 3) %>% select(ptm, dis) %>% unique
  # q %>% group_by(ptm) %>% filter(n_distinct(dis) == 3) %>% select(ptm) %>% unique
  # q %>% filter(ptm %in% (q %>% group_by(ptm) %>% filter(n_distinct(dis) == 3) %>% pull(ptm) %>% unique)) %>% group_by(ptm) %>% filter(n_distinct(dis) != 3) %>% select(ptm, dis) %>% unique
  # q %>% filter(ptm %in% (q %>% group_by(ptm) %>% filter(n_distinct(dis) == 3) %>% pull(ptm) %>% unique)) %>% group_by(ptm) %>% filter(n_distinct(dis) == 3) %>% select(ptm, dis) %>% unique
  # q %>% filter(ptm=="S-p") %>% group_by(ptm, dis) %>% do(wilcox_test = wilcox.test(af ~ type, data = .)) %>% summarise(ptm, dis, p_value = wilcox_test$p.value)
  # Get PTMs where all 3 dis types and ptm and control are present
  myptms <- q %>% group_by(ptm) %>% summarise(combinations = n_distinct(dis, type)) %>% filter(combinations == 6) %>% pull(ptm) %>% unique
  myptms
  # qsum AF
  qsum <- q %>% 
    filter(ptm %in% myptms) %>% 
    nest_by(ptm, dis) %>% 
    mutate(
      af_mean_ptm = data %>% filter(type == "ptm") %>% pull(af) %>% mean,
      af_mean_control = data %>% filter(type == "control") %>% pull(af) %>% mean,
      af_mean_dif = af_mean_ptm - af_mean_control,
      af_median_ptm = data %>% filter(type == "ptm") %>% pull(af) %>% median,
      af_median_control = data %>% filter(type == "control") %>% pull(af) %>% median,
      af_median_dif = af_median_ptm - af_median_control,
      # wilcox_test = wilcox.test(af ~ type, data = data, exact = F)$p.value
      wilcox_test = tryCatch(
        expr = wilcox.test(af ~ type, data = data, exact = F)$p.value,
        error = function(e) NA_real_
      )
    ) %>% 
    summarise(ptm, dis, af_mean_ptm, af_mean_control, af_mean_dif, af_median_ptm, af_median_control, af_median_dif, p_value = wilcox_test)
  qsum
  qsum %>% summary
  # qsum AF and frac0
  qsum <- q %>% 
    filter(ptm %in% myptms) %>% 
    nest_by(ptm, dis) %>% 
    mutate(
      af_mean_ptm = data %>% filter(type == "ptm") %>% pull(af) %>% mean,
      af_mean_control = data %>% filter(type == "control") %>% pull(af) %>% mean,
      af_mean_dif = af_mean_ptm - af_mean_control,
      af_median_ptm = data %>% filter(type == "ptm") %>% pull(af) %>% median,
      af_median_control = data %>% filter(type == "control") %>% pull(af) %>% median,
      af_median_dif = af_median_ptm - af_median_control,
      # wilcox_test = wilcox.test(af ~ type, data = data, exact = F)$p.value,
      wilcox_test = tryCatch(
        expr = wilcox.test(af ~ type, data = data, exact = F)$p.value,
        error = function(e) NA_real_
      ),
      n0_ptm = data %>% filter(type == "ptm" & af == 0) %>% nrow,
      n0_control = data %>% filter(type == "control" & af == 0) %>% nrow,
      n_ptm = data %>% filter(type == "ptm") %>% nrow,
      n_control = data %>% filter(type == "control") %>% nrow,
      frac0_ptm = n0_ptm / n_ptm,
      frac0_control = n0_control / n_control,
      frac0_dif = frac0_ptm - frac0_control,
      # frac0_fisher = fisher.test(matrix(c(n0_ptm, n0_control, n_ptm - n0_ptm, n_control - n0_control), nrow = 2))$p.value
      frac0_fisher = tryCatch(
        expr = fisher.test(matrix(c(n0_ptm, n0_control, n_ptm - n0_ptm, n_control - n0_control), nrow = 2))$p.value,
        error = function(e) NA_real_
      )
    ) %>%
    summarise(ptm, dis, af_mean_ptm, af_mean_control, af_mean_dif, af_median_ptm, af_median_control, af_median_dif, p_value = wilcox_test, frac0_ptm, frac0_control, frac0_dif, frac0_p_value = frac0_fisher)
  qsum
  qsum %>% summary
  
  # How many significant AF Wilcoxon p-values?
  qsum %>% filter(p_value < 0.05) %>% nrow
  # >> 29
  # How many significant frac0 Fisher p-values?
  qsum %>% filter(frac0_p_value < 0.05) %>% nrow
  # >> 30
  # How many significant for both?
  qsum %>% filter(p_value < 0.05 & frac0_p_value < 0.05) %>% nrow
  # >> 25 out of 30
  
  # Apply FDR correction to p-values
  qsum$p_value_fdr <- p.adjust(qsum$p_value, method = "BH")
  qsum$frac0_p_value_fdr <- p.adjust(qsum$frac0_p_value, method = "BH")
  qsum
  
  # How many significant AF Wilcoxon p-values after FDR?
  qsum %>% filter(p_value_fdr < 0.05) %>% nrow
  # >> 23
  # How many significant frac0 Fisher p-values after FDR?
  qsum %>% filter(frac0_p_value_fdr < 0.05) %>% nrow
  # >> 23
  # How many significant for both after FDR?
  qsum %>% filter(p_value_fdr < 0.05 & frac0_p_value_fdr < 0.05) %>% nrow
  # >> 19 out of 23
  
  # # Compare PTM vs. Control af_mean and af_median (simply across ptm|dis, here)
  # wilcox.test(qsum$af_mean_ptm, qsum$af_mean_control)
  # wilcox.test(qsum$af_median_ptm, qsum$af_median_control)
  # 
  # # Compare PTM vs. Control frac0 (simply across ptm|dis, here)
  # wilcox.test(qsum$frac0_ptm, qsum$frac0_control)
  
  # # Get PTM types with at least 1000 sites
  # # ptms1000 <- Query("SELECT ptm FROM (SELECT ptm, COUNT(DISTINCT acc, site) AS c FROM unimod WHERE species='HUMAN' AND ptm IS NOT NULL GROUP BY ptm ORDER BY ptm) AS t WHERE c > 1000");
  # ptms1000 <- Query("SELECT ptm FROM (SELECT ptm, COUNT(DISTINCT acc, site) AS c FROM unimod WHERE species='HUMAN' AND ptm IS NOT NULL GROUP BY ptm ORDER BY c DESC) AS t WHERE c > 1000");
  # ptms1000 %<>% pull(ptm)
  # ptms1000
  # Get PTM types with at least 1000 sites
  ptms1000 <- c("R-me", "K-sum", "S-gly", "K-me", "T-gly", "S-p", "K-ac", "T-p", "K-ub", "K-mal", "K-suc", "N-gly", "Y-p", "M-ox", "C-pal", "C-glt", "C-nit") # Order as in alphasa "all ptms" figure (all >1000, ptms1000 I suppose)
  
  
  # # Get PTM types with at least 2000 sites
  # # ptms2000 <- Query("SELECT ptm FROM (SELECT ptm, COUNT(DISTINCT acc, site) AS c FROM unimod WHERE species='HUMAN' AND ptm IS NOT NULL GROUP BY ptm ORDER BY ptm) AS t WHERE c > 2000");
  # ptms2000 <- Query("SELECT ptm FROM (SELECT ptm, COUNT(DISTINCT acc, site) AS c FROM unimod WHERE species='HUMAN' AND ptm IS NOT NULL GROUP BY ptm ORDER BY c DESC) AS t WHERE c > 2000");
  # ptms2000 %<>% pull(ptm)
  # ptms2000
  
  # Set standard PTM types
  stdptms <- c('S-p', 'T-p', 'Y-p', 'K-ac', 'K-mal', 'K-me', 'K-suc', 'K-sum', 'K-ub', 'N-gly', 'R-me', 'S-gly', 'T-gly')
  # favptms <- c('S-p', 'T-p', 'Y-p', 'K-ac', 'K-mal', 'K-me', 'K-suc', 'K-sum', 'K-ub', 'N-gly', 'R-me')
  ptms17 <- c("R-me", "K-sum", "S-gly", "K-me", "T-gly", "S-p", "K-ac", "T-p", "K-ub", "K-mal", "K-suc", "N-gly", "Y-p", "M-ox", "C-pal", "C-glt", "C-nit") # as in alphasa "all ptms" figure (all >1000, ptms1000 I suppose)
  ptms11 <- c("S-p", "T-p", "Y-p", "K-ub", "K-sum", "K-ac", "K-mal", "K-suc", "K-me", "R-me", "N-gly") # select11
  ptms9 <- c("S-p", "T-p", "Y-p", "K-ub", "K-ac", "K-mal", "K-suc", "K-me", "R-me") # select9
  favptms <- ptms11
  mainptms <- c('S-p', 'T-p', 'Y-p', 'K-ac', 'K-ub')
  
  # What fraction of PTM and Control has af=0?
  q %>% filter(ptm %in% mainptms) %>% group_by(ptm, type) %>% summarise(mean_af = mean(af), median_af = median(af), median_af_non0 = median(ifelse(af == 0, NA, af), na.rm=T), n0 = sum(af == 0), n = n(), frac0 = n0 / n) %>% print(n = 1000)
  q %>% filter(ptm %in% mainptms) %>% group_by(ptm, dis, type) %>% summarise(mean_af = mean(af), median_af = median(af), median_af_non0 = median(ifelse(af == 0, NA, af), na.rm=T), n0 = sum(af == 0), n = n(), frac0 = n0 / n) %>% print(n = 1000)
  # Summarise to calculate frac0 PTM - frac0 Control
  q %>% filter(ptm %in% mainptms) %>% group_by(ptm, dis, type) %>% summarise(mean_af = mean(af), median_af = median(af), median_af_non0 = median(ifelse(af == 0, NA, af), na.rm=T), n0 = sum(af == 0), n = n(), frac0 = n0 / n) %>% group_by(ptm, dis) %>% summarise(frac0_ptm = frac0[type == "ptm"], frac0_control = frac0[type == "control"], frac0_dif = frac0[type == "ptm"] - frac0[type == "control"]) %>% print(n = 1000)
  q %>% filter(ptm %in% mainptms) %>% group_by(ptm, dis, type) %>% summarise(mean_af = mean(af), median_af = median(af), median_af_non0 = median(ifelse(af == 0, NA, af), na.rm=T), n0 = sum(af == 0), n = n(), frac0 = n0 / n) %>% group_by(ptm, dis) %>% summarise(frac0_ptm = frac0[type == "ptm"], frac0_control = frac0[type == "control"], frac0_dif = frac0[type == "ptm"] - frac0[type == "control"]) %>% arrange(desc(frac0_dif)) %>% print(n = 1000)
  # Summarise to calculate frac0 PTM - frac0 Control, as well as mean_af and median_af PTM - Control
  q %>% filter(ptm %in% mainptms) %>% group_by(ptm, dis, type) %>% summarise(mean_af = mean(af), median_af = median(af), median_af_non0 = median(ifelse(af == 0, NA, af), na.rm=T), n0 = sum(af == 0), n = n(), frac0 = n0 / n) %>% group_by(ptm, dis) %>% summarise(frac0_ptm = frac0[type == "ptm"], frac0_control = frac0[type == "control"], frac0_dif = frac0[type == "ptm"] - frac0[type == "control"], mean_af_dif = mean_af[type == "ptm"] - mean_af[type == "control"], median_af_dif = median_af[type == "ptm"] - median_af[type == "control"], median_af_non0_dif = median_af_non0[type == "ptm"] - median_af_non0[type == "control"]) %>% arrange(desc(frac0_dif)) %>% print(n = 1000)
  q %>% filter(ptm %in% mainptms) %>% group_by(ptm, dis, type) %>% summarise(mean_af = mean(af), median_af = median(af), median_af_non0 = median(ifelse(af == 0, NA, af), na.rm=T), n0 = sum(af == 0), n = n(), frac0 = n0 / n) %>% group_by(ptm, dis) %>% summarise(frac0_ptm = frac0[type == "ptm"], frac0_control = frac0[type == "control"], frac0_dif = frac0[type == "ptm"] - frac0[type == "control"], mean_af_dif = mean_af[type == "ptm"] - mean_af[type == "control"], median_af_dif = median_af[type == "ptm"] - median_af[type == "control"], median_af_non0_dif = median_af_non0[type == "ptm"] - median_af_non0[type == "control"]) %>% arrange(desc(mean_af_dif)) %>% print(n = 1000)
  q %>% filter(ptm %in% mainptms) %>% group_by(ptm, dis, type) %>% summarise(mean_af = mean(af), median_af = median(af), median_af_non0 = median(ifelse(af == 0, NA, af), na.rm=T), n0 = sum(af == 0), n = n(), frac0 = n0 / n) %>% group_by(ptm, dis) %>% summarise(frac0_ptm = frac0[type == "ptm"], frac0_control = frac0[type == "control"], frac0_dif = frac0[type == "ptm"] - frac0[type == "control"], mean_af_dif = mean_af[type == "ptm"] - mean_af[type == "control"], median_af_dif = median_af[type == "ptm"] - median_af[type == "control"], median_af_non0_dif = median_af_non0[type == "ptm"] - median_af_non0[type == "control"]) %>% arrange(desc(median_af_non0_dif)) %>% print(n = 1000)
  # >> Buried PTM sites are more often invariant, while disordered PTM sites more often have variants.
  
  # AF geom_histogram
  { q %>%
      mutate(af = ifelse(af == 0, 1e-10, af)) %>%
      ggplot(aes(x = af, fill = type, colour = type)) +
      geom_density(alpha = 0.3) +
      scale_colour_viridis_d(aesthetics = c("fill", "colour")) +
      scale_x_log10() +
      theme_nature() } %>%
    qsave(f("output-snps_gnomad-af-density{tmp_plddt}.pdf"))
  
  # AF geom_tile
  qsum$af_mean_dif %>% summary
  min(qsum$af_mean_dif)
  max(qsum$af_mean_dif)
  qsum %>%
    # mutate(af_mean_dif = ifelse(af_mean_dif > 0, NA, af_mean_dif)) %>%
    mutate(af_mean_dif = ifelse(p_value >= 0.05, NA, af_mean_dif)) %>%
    ggplot(aes(x = ptm, y = fct_rev(dis), fill = af_mean_dif)) +
    geom_tile() +
    # geom_text(aes(label = ifelse(is.na(af_mean_dif), "", round(af_mean_dif, 4))), size = 1) +
    geom_text(aes(label = ifelse(af_mean_dif > 0, "+", "")), size = 2) +      # "+" means PTM has more variants than Control
    # scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    # scale_fill_viridis_c(na.value = "white", breaks = pretty_breaks()) +
    # scale_fill_viridis_c(limits = c(-max(abs(qsum$af_mean_dif), na.rm = T), max(abs(qsum$af_mean_dif), na.rm = T)), breaks = pretty_breaks(8), na.value = "white") +
    scale_fill_viridis_c(limits = c(-max(abs(qsum$af_mean_dif), na.rm = T), max(abs(qsum$af_mean_dif), na.rm = T)), na.value = "white") +
    # scale_fill_viridis_c(limits = c(-max(qsum$af_mean_dif, na.rm = T), max(qsum$af_mean_dif, na.rm = T)), breaks = pretty_breaks(8), na.value = "white") +
    # scale_fill_viridis_c(na.value = "white") +
    # scale_fill_viridis_c(option = "E") + # cividis
    # scale_fill_gradient2(midpoint = 0) +
    # scale_fill_viridis_c(option = "turbo", na.value = "white", limits = c(min(qsum$af_mean_dif), max(qsum$af_mean_dif))) + # turbo (diverging)
    # scale_fill_viridis_c(option = "turbo", na.value = "white", limits = c(-max(abs(qsum$af_mean_dif)), max(abs(qsum$af_mean_dif)))) + # turbo (diverging)
    # scale_fill_gradientn(colours=c("blue", "black", "yellow"), na.value="white", limits = c(-max(abs(qsum$af_mean_dif), na.rm = T), max(abs(qsum$af_mean_dif), na.rm = T)), breaks = pretty_breaks()) +
    # scale_fill_gradientn(colours=c(viridis(2)[1], "white", viridis(2)[2]), na.value="white", limits = c(-max(abs(qsum$af_mean_dif), na.rm = T), max(abs(qsum$af_mean_dif), na.rm = T)), breaks = pretty_breaks()) +
    # scale_fill_gradientn(colours=c(viridis(2)[1], "white", viridis(2)[2]), na.value="white", limits = c(min(qsum$af_mean_dif, na.rm = T), max(qsum$af_mean_dif, na.rm = T)), breaks = pretty_breaks()) +
    # scale_fill_gradientn(colours=c(viridis(2)[1], "black", viridis(2)[2]), na.value="white", limits = c(min(qsum$af_mean_dif, na.rm = T), max(qsum$af_mean_dif, na.rm = T)), breaks = pretty_breaks()) +
    # scale_fill_gradientn(colours=c(ptmvir0, "black", ptmvir1), na.value="white", limits = c(min(qsum$af_mean_dif, na.rm = T), max(qsum$af_mean_dif, na.rm = T)), breaks = pretty_breaks()) +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # theme(legend.position = "right") +
    # theme(legend.key.size = unit(1, "cm"))
    theme(legend.position = "none")
  qsave(f("output-snps_gnomad-af_mean_dif{tmp_plddt}.pdf"), width = 120, height = 20)
  
  # AF geom_tile vertical
  tmpptm <- qsum %>% mutate(af_mean_dif = ifelse(p_value >= 0.05, 0, af_mean_dif)) %>% group_by(ptm) %>% summarise(mean_af_mean_dif = mean(af_mean_dif)) %>% arrange(desc(mean_af_mean_dif)) %>% pull(ptm)
  qsum %>%
    # mutate(af_mean_dif = ifelse(af_mean_dif > 0, NA, af_mean_dif)) %>%
    arrange(desc(af_mean_dif)) %>%
    # mutate(af_mean_dif = ifelse(p_value >= 0.05, NA, af_mean_dif)) %>%
    ggplot(aes(x = dis, y = factor(ptm, levels = tmpptm), fill = af_mean_dif)) +
    geom_tile() +
    geom_text(aes(label = ifelse(af_mean_dif > 0, "+", "")), size = 2) +      # "+" means PTM has more variants than Control
    # scale_fill_viridis_c(limits = c(-max(abs(qsum$af_mean_dif), na.rm = T), max(abs(qsum$af_mean_dif), na.rm = T)), na.value = "white") +
    scale_fill_viridis_c(na.value = "white") +
    theme_nature() +
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "none") +
    xlab(NULL) +
    ylab(NULL)
  qsave(f("output-snps_gnomad-af_mean_dif-vertical{tmp_plddt}.pdf"), height = 80)
  
  # AF geom_tile vertical filtered
  tmpptm <- qsum %>% mutate(af_mean_dif = ifelse(p_value >= 0.05, 0, af_mean_dif)) %>% group_by(ptm) %>% summarise(mean_af_mean_dif = mean(af_mean_dif)) %>% arrange(desc(mean_af_mean_dif)) %>% pull(ptm)
  qsum %>%
    # mutate(af_mean_dif = ifelse(af_mean_dif > 0, NA, af_mean_dif)) %>%
    arrange(desc(af_mean_dif)) %>%
    mutate(af_mean_dif = ifelse(p_value >= 0.05, NA, af_mean_dif)) %>%
    ggplot(aes(x = dis, y = factor(ptm, levels = tmpptm), fill = af_mean_dif)) +
    geom_tile() +
    geom_text(aes(label = ifelse(af_mean_dif > 0, "+", "")), size = 2) +      # "+" means PTM has more variants than Control
    # scale_fill_viridis_c(limits = c(-max(abs(qsum$af_mean_dif), na.rm = T), max(abs(qsum$af_mean_dif), na.rm = T)), na.value = "white") +
    scale_fill_viridis_c(na.value = "white") +
    theme_nature() +
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "none") +
    xlab(NULL) +
    ylab(NULL)
  qsave(f("output-snps_gnomad-af_mean_dif-vertical-filtered{tmp_plddt}.pdf"), height = 80)
  
  # AF geom_tile vertical filtered (FDR)
  tmpptm <- qsum %>% mutate(af_mean_dif = ifelse(p_value >= 0.05, 0, af_mean_dif)) %>% group_by(ptm) %>% summarise(mean_af_mean_dif = mean(af_mean_dif)) %>% arrange(desc(mean_af_mean_dif)) %>% pull(ptm)
  qsum %>%
    # mutate(af_mean_dif = ifelse(af_mean_dif > 0, NA, af_mean_dif)) %>%
    arrange(desc(af_mean_dif)) %>%
    mutate(af_mean_dif = ifelse(p_value_fdr >= 0.05, NA, af_mean_dif)) %>%
    ggplot(aes(x = dis, y = factor(ptm, levels = tmpptm), fill = af_mean_dif)) +
    geom_tile() +
    geom_text(aes(label = ifelse(af_mean_dif > 0, "+", "")), size = 2) +      # "+" means PTM has more variants than Control
    # scale_fill_viridis_c(limits = c(-max(abs(qsum$af_mean_dif), na.rm = T), max(abs(qsum$af_mean_dif), na.rm = T)), na.value = "white") +
    scale_fill_viridis_c(na.value = "white") +
    theme_nature() +
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "none") +
    xlab(NULL) +
    ylab(NULL)
  qsave(f("output-snps_gnomad-af_mean_dif-vertical-filtered_fdr{tmp_plddt}.pdf"), height = 80)
  
  # AF geom_tile vertical filtered (FDR) mainptms
  tmpptm <- qsum %>% arrange(desc(ptm)) %>% pull(ptm) %>% unique
  qsum %>%
    # mutate(af_mean_dif = ifelse(af_mean_dif > 0, NA, af_mean_dif)) %>%
    arrange(desc(af_mean_dif)) %>%
    mutate(af_mean_dif = ifelse(p_value_fdr >= 0.05, NA, af_mean_dif)) %>%
    filter(ptm %in% mainptms) %>%
    ggplot(aes(x = dis, y = factor(ptm, levels = tmpptm), fill = af_mean_dif)) +
    geom_tile() +
    geom_text(aes(label = ifelse(af_mean_dif > 0, "+", "")), size = 2) +      # "+" means PTM has more variants than Control
    geom_text(aes(label = ifelse(is.na(af_mean_dif), "", round(af_mean_dif, 6))), size = 1) +
    # scale_fill_viridis_c(limits = c(-max(abs(qsum$af_mean_dif), na.rm = T), max(abs(qsum$af_mean_dif), na.rm = T)), na.value = "white") +
    scale_fill_viridis_c(na.value = "white") +
    # scale_fill_viridis_c(direction = -1, na.value = "white") +
    theme_nature() +
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "none") +
    xlab(NULL) +
    ylab(NULL)
  qsave(f("output-snps_gnomad-af_mean_dif-vertical-filtered_fdr-mainptms{tmp_plddt}.pdf"))
  # >> High af_mean_dif means the PTM sites have higher minor allele frequencies than Control.
  
  
  
  # # AC0 geom_tile vertical mainptms
  # qsum is now identical to qfrac0, no need to retain qfrac0
  # qfrac0 <- q %>% filter(ptm %in% mainptms) %>% group_by(ptm, dis, type) %>% summarise(mean_af = mean(af), median_af = median(af), median_af_non0 = median(ifelse(af == 0, NA, af), na.rm=T), n0 = sum(af == 0), n = n(), frac0 = n0 / n) %>% group_by(ptm, dis) %>% summarise(frac0_ptm = frac0[type == "ptm"], frac0_control = frac0[type == "control"], frac0_dif = frac0[type == "ptm"] - frac0[type == "control"], mean_af_dif = mean_af[type == "ptm"] - mean_af[type == "control"], median_af_dif = median_af[type == "ptm"] - median_af[type == "control"], median_af_non0_dif = median_af_non0[type == "ptm"] - median_af_non0[type == "control"]) %>% arrange(desc(frac0_dif))
  # qfrac0 %>% print(n = 1000)
  # q %>% filter(ptm %in% mainptms) %>% group_by(ptm, dis, type) %>% summarise(mean_af = mean(af), median_af = median(af), median_af_non0 = median(ifelse(af == 0, NA, af), na.rm=T), n0 = sum(af == 0), n = n(), frac0 = n0 / n) %>% group_by(ptm, dis) %>% summarise(frac0_ptm = frac0[type == "ptm"], frac0_control = frac0[type == "control"], frac0_dif = frac0[type == "ptm"] - frac0[type == "control"], mean_af_dif = mean_af[type == "ptm"] - mean_af[type == "control"], median_af_dif = median_af[type == "ptm"] - median_af[type == "control"], median_af_non0_dif = median_af_non0[type == "ptm"] - median_af_non0[type == "control"]) %>% arrange(desc(frac0_dif)) %>%
  #   ggplot(aes(x = dis, y = factor(ptm, levels = tmpptm), fill = frac0_dif)) +
  #   geom_tile() +
  #   geom_text(aes(label = ifelse(frac0_dif > 0, "+", "")), size = 2) +
  #   scale_fill_viridis_c(na.value = "white") +
  #   theme_nature() +
  #   theme(legend.position = "none") +
  #   xlab(NULL) +
  #   ylab(NULL)
  # 
  # # AC0 geom_tile vertical
  # qfrac0 <- q %>% filter(ptm %in% myptms) %>% group_by(ptm, dis, type) %>% summarise(mean_af = mean(af), median_af = median(af), median_af_non0 = median(ifelse(af == 0, NA, af), na.rm=T), n0 = sum(af == 0), n = n(), frac0 = n0 / n) %>% group_by(ptm, dis) %>% summarise(frac0_ptm = frac0[type == "ptm"], frac0_control = frac0[type == "control"], frac0_dif = frac0[type == "ptm"] - frac0[type == "control"], mean_af_dif = mean_af[type == "ptm"] - mean_af[type == "control"], median_af_dif = median_af[type == "ptm"] - median_af[type == "control"], median_af_non0_dif = median_af_non0[type == "ptm"] - median_af_non0[type == "control"]) %>% arrange(desc(frac0_dif))
  # # tmpptm <- qfrac0 %>% group_by(ptm) %>% summarise(mean_frac0_dif = mean(frac0_dif)) %>% arrange(desc(mean_frac0_dif)) %>% pull(ptm)
  # tmpptm <- qfrac0 %>% group_by(ptm) %>% summarise(mean_frac0_dif = mean(frac0_dif)) %>% arrange(mean_frac0_dif) %>% pull(ptm)
  # qfrac0 %>% 
  #   ggplot(aes(x = dis, y = factor(ptm, levels = tmpptm), fill = frac0_dif)) +
  #   geom_tile() +
  #   geom_text(aes(label = ifelse(frac0_dif > 0, "+", "")), size = 2) +
  #   scale_fill_viridis_c(na.value = "white") +
  #   theme_nature() +
  #   theme(legend.position = "none") +
  #   xlab(NULL) +
  #   ylab(NULL)
  # qsave("output-snps_gnomad-frac0_dif-vertical.pdf", height = 80)
  
  # AC0 geom_tile vertical
  tmpptm <- qsum %>% group_by(ptm) %>% summarise(mean_frac0_dif = mean(frac0_dif)) %>% arrange(mean_frac0_dif) %>% pull(ptm)
  qsum %>% 
    # mutate(frac0_dif = ifelse(p_value >= 0.05, NA, frac0_dif)) %>%
    ggplot(aes(x = dis, y = factor(ptm, levels = tmpptm), fill = frac0_dif)) +
    geom_tile() +
    geom_text(aes(label = ifelse(frac0_dif < 0, "+", "")), size = 2) +      # "+" means PTM has more variants than Control
    scale_fill_viridis_c(na.value = "white") +
    theme_nature() +
    theme(legend.position = "none") +
    xlab(NULL) +
    ylab(NULL)
  qsave(f("output-snps_gnomad-frac0_dif-vertical{tmp_plddt}.pdf"), height = 80)
  
  # AC0 geom_tile vertical filtered
  tmpptm <- qsum %>% group_by(ptm) %>% summarise(mean_frac0_dif = mean(frac0_dif)) %>% arrange(mean_frac0_dif) %>% pull(ptm)
  qsum %>%
    mutate(frac0_dif = ifelse(frac0_p_value_fdr >= 0.05, NA, frac0_dif)) %>%
    filter(ptm %in% tmpptm) %>% pull(frac0_dif) %>% abs %>% min(na.rm=T) -> minfrac0_dif
  # Get maxfrac0_dif (highest absolute value for the colour scale, which will be symmetric around 0)
  qsum %>%
    mutate(frac0_dif = ifelse(frac0_p_value_fdr >= 0.05, NA, frac0_dif)) %>%
    filter(ptm %in% tmpptm) %>% pull(frac0_dif) %>% abs %>% max(na.rm=T) -> maxfrac0_dif
  maxfrac0_dif
  qsum %>% 
    mutate(frac0_dif = ifelse(p_value >= 0.05, NA, frac0_dif)) %>%
    ggplot(aes(x = dis, y = factor(ptm, levels = tmpptm), fill = frac0_dif)) +
    geom_tile() +
    # geom_text(aes(label = ifelse(frac0_dif < 0, "+", "")), size = 2) +      # "+" means PTM has more variants than Control
    scale_fill_gradientn(colours = c(scales::viridis_pal()(5)[1], scales::viridis_pal()(5)[2], "#FFFFFFFF", scales::viridis_pal()(5)[4], scales::viridis_pal()(5)[5]), values = scales::rescale(c(-maxfrac0_dif, -minfrac0_dif, 0, minfrac0_dif, maxfrac0_dif)), na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = NULL) +
    theme_nature() +
    theme(legend.position = "none") +
    xlab(NULL) +
    ylab(NULL)
  qsave(f("output-snps_gnomad-frac0_dif-vertical-filtered{tmp_plddt}.pdf"), height = 80)
  
  # AC0 geom_tile vertical filtered (FDR)
  # tmpptm <- qsum %>% group_by(ptm) %>% summarise(mean_frac0_dif = mean(frac0_dif)) %>% arrange(mean_frac0_dif) %>% pull(ptm)
  tmpptm <- ptms17
  qsum %>%
    mutate(frac0_dif = ifelse(frac0_p_value_fdr >= 0.05, NA, frac0_dif)) %>%
    filter(ptm %in% tmpptm) %>% pull(frac0_dif) %>% abs %>% min(na.rm=T) -> minfrac0_dif
  # Get maxfrac0_dif (highest absolute value for the colour scale, which will be symmetric around 0)
  qsum %>%
    mutate(frac0_dif = ifelse(frac0_p_value_fdr >= 0.05, NA, frac0_dif)) %>%
    filter(ptm %in% tmpptm) %>% pull(frac0_dif) %>% abs %>% max(na.rm=T) -> maxfrac0_dif
  maxfrac0_dif
  qsum %>% 
    mutate(frac0_dif = ifelse(frac0_p_value_fdr >= 0.05, NA, frac0_dif)) %>%
    filter(ptm %in% tmpptm) %>%
    ggplot(aes(x = dis, y = fct_rev(factor(ptm, levels = tmpptm)), fill = -frac0_dif)) +
    geom_tile() +
    # geom_text(aes(label = ifelse(frac0_dif < 0, "+", "")), size = 2) +      # "+" means PTM has more variants than Control
    # scale_fill_gradientn(colours = c(scales::viridis_pal()(5)[1], scales::viridis_pal()(5)[2], "#FFFFFFFF", scales::viridis_pal()(5)[4], scales::viridis_pal()(5)[5]), values = scales::rescale(c(-maxfrac0_dif, -minfrac0_dif, 0, minfrac0_dif, maxfrac0_dif)), na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = NULL) +
    scale_fill_gradientn(colours = c(scales::viridis_pal()(5)[1], scales::viridis_pal()(5)[2], "#FFFFFFFF", scales::viridis_pal()(5)[4], scales::viridis_pal()(5)[5]), values = scales::rescale(c(-maxfrac0_dif, -minfrac0_dif, 0, minfrac0_dif, maxfrac0_dif)), na.value = "white", breaks = c(-0.25, 0, 0.25), labels = c("-0.25", "0", "0.25"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = "Variable residues") +
    theme_nature(legend_position = "bottom", legend_nudge_top=-8.25, legend_nudge_right=-1.5) +
    theme(legend.text = element_text(size = 5, margin = margin(t = 1))) +
    theme(axis.text.y = element_text(size = 6)) +
    theme(legend.title.position = "left", legend.title = element_text(size = 5, vjust = 1, margin = margin(t = unit(0.75, "mm"), r = unit(1, "mm")))) +
    xlab(NULL) +
    ylab(NULL)
  qsave(f("output-snps_gnomad-frac0_dif-vertical-filtered_fdr{tmp_plddt}.pdf"), height = 80)
  
  # AC0 geom_tile horizontal filtered (FDR) in same order as alphasa_relasa_vs_ptms/output-density-canonical-allptms-weighted.pdf (17 ptms)
  tmpptm <- c("R-me", "K-sum", "S-gly", "K-me", "T-gly", "S-p", "K-ac", "T-p", "K-ub", "K-mal", "K-suc", "N-gly", "Y-p", "M-ox", "C-pal", "C-glt", "C-nit")
  # Get minfrac0_dif (lowest absolute frac0_dif where frac0_p_value_fdr is significant (to make sure it falls into what we plot as white on the fill scale))
  qsum %>%
    mutate(frac0_dif = ifelse(frac0_p_value_fdr >= 0.05, NA, frac0_dif)) %>%
    filter(ptm %in% tmpptm) %>% pull(frac0_dif) %>% abs %>% min(na.rm=T) -> minfrac0_dif
  # Get maxfrac0_dif (highest absolute value for the colour scale, which will be symmetric around 0)
  qsum %>%
    mutate(frac0_dif = ifelse(frac0_p_value_fdr >= 0.05, NA, frac0_dif)) %>%
    filter(ptm %in% tmpptm) %>% pull(frac0_dif) %>% abs %>% max(na.rm=T) -> maxfrac0_dif
  maxfrac0_dif
  qsum %>% 
    mutate(frac0_dif = ifelse(frac0_p_value_fdr >= 0.05, NA, frac0_dif)) %>%
    filter(ptm %in% tmpptm) %>%
    ggplot(aes(x = factor(ptm, levels = tmpptm), y = dis, fill = -frac0_dif)) +
    geom_tile() +
    scale_fill_gradientn(colours = c(scales::viridis_pal()(5)[1], scales::viridis_pal()(5)[2], "#FFFFFFFF", scales::viridis_pal()(5)[4], scales::viridis_pal()(5)[5]), values = scales::rescale(c(-maxfrac0_dif, -minfrac0_dif, 0, minfrac0_dif, maxfrac0_dif)), na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = NULL) +
    theme_nature(legend_position = "top") +
    theme(text = element_text(family = "Helvetica Neue", size = 5)) + # Base font size
    theme(legend.text = element_text(size = 5, margin = margin(t = 1))) +
    # legend.key.spacing.y = unit(0, "mm"),
    # theme(legend.position = "none") +
    # ggtitle("") +
    xlab(NULL) +
    # xlab("") +
    ylab(NULL)
  qsave(f("output-snps_gnomad-frac0_dif-horizontal-filtered_fdr-alphasa_allptms{tmp_plddt}.pdf"), width = 80, height = 30)
  
  # AC0 geom_tile vertical filtered (FDR) mainptms
  tmpptm <- qsum %>% arrange(desc(ptm)) %>% pull(ptm) %>% unique
  qsum %>% 
    mutate(frac0_dif = ifelse(frac0_p_value_fdr >= 0.05, NA, frac0_dif)) %>%
    filter(ptm %in% mainptms) %>%
    ggplot(aes(x = dis, y = factor(ptm, levels = tmpptm), fill = frac0_dif)) +
    geom_tile() +
    geom_text(aes(label = ifelse(frac0_dif < 0, "+", "")), size = 2) +      # "+" means PTM has more variants than Control
    scale_fill_viridis_c(na.value = "white") +
    theme_nature() +
    theme(legend.position = "none") +
    xlab(NULL) +
    ylab(NULL)
  qsave(f("output-snps_gnomad-frac0_dif-vertical-filtered_fdr-mainptms{tmp_plddt}.pdf"))
  # >> High frac0_dif means more invariant sites in PTM than Control. An example is M-ac: the initiator methionine is of course usually invariant.
  
  # Useful guide to legends:
  # https://www.tidyverse.org/blog/2024/02/ggplot2-3-5-0-legends/
  
  # AC0 geom_tile vertical filtered (FDR) favptms (select11)
  tmpptm <- ptms11
  # Get minfrac0_dif (lowest absolute frac0_dif where frac0_p_value_fdr is significant (to make sure it falls into what we plot as white on the fill scale))
  qsum %>%
    mutate(frac0_dif = ifelse(frac0_p_value_fdr >= 0.05, NA, frac0_dif)) %>%
    filter(ptm %in% tmpptm) %>% pull(frac0_dif) %>% abs %>% min(na.rm=T) -> minfrac0_dif
  # Get maxfrac0_dif (highest absolute value for the colour scale, which will be symmetric around 0)
  qsum %>%
    mutate(frac0_dif = ifelse(frac0_p_value_fdr >= 0.05, NA, frac0_dif)) %>%
    filter(ptm %in% tmpptm) %>% pull(frac0_dif) %>% abs %>% max(na.rm=T) -> maxfrac0_dif
  maxfrac0_dif
  qsum %>% 
    mutate(frac0_dif = ifelse(frac0_p_value_fdr >= 0.05, NA, frac0_dif)) %>%
    filter(ptm %in% tmpptm) %>%
    ggplot(aes(x = dis, y = fct_rev(factor(ptm, levels = tmpptm)), fill = -frac0_dif)) +
    geom_tile() +
    scale_fill_gradientn(colours = c(scales::viridis_pal()(5)[1], scales::viridis_pal()(5)[2], "#FFFFFFFF", scales::viridis_pal()(5)[4], scales::viridis_pal()(5)[5]), values = scales::rescale(c(-maxfrac0_dif, -minfrac0_dif, 0, minfrac0_dif, maxfrac0_dif)), na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = "Variable residues") +
    theme_nature(legend_position = "bottom", legend_nudge_top=-8.25, legend_nudge_right=-0.25) +
    # theme(text = element_text(family = "Helvetica Neue", size = 5)) + # Base font size
    theme(legend.text = element_text(size = 5, margin = margin(t = 1))) +
    # theme(legend.text = element_text(size = 6, margin = margin(t = 1))) +
    theme(axis.text.y = element_text(size = 6)) +
    theme(legend.title.position = "left", legend.title = element_text(size = 5, vjust = 1, margin = margin(t = unit(0.75, "mm"), r = unit(1, "mm")))) +
    xlab(NULL) +
    ylab(NULL)
  qsave(f("output-snps_gnomad-frac0_dif-vertical-filtered_fdr-select11{tmp_plddt}.pdf"), height = 70)
  
  # AC0 geom_tile vertical filtered (FDR) select9
  tmpptm <- ptms9
  # Get minfrac0_dif (lowest absolute frac0_dif where frac0_p_value_fdr is significant (to make sure it falls into what we plot as white on the fill scale))
  qsum %>%
    mutate(frac0_dif = ifelse(frac0_p_value_fdr >= 0.05, NA, frac0_dif)) %>%
    filter(ptm %in% ptms9) %>% pull(frac0_dif) %>% abs %>% min(na.rm=T) -> minfrac0_dif
  # Get maxfrac0_dif (highest absolute value for the colour scale, which will be symmetric around 0)
  qsum %>%
    mutate(frac0_dif = ifelse(frac0_p_value_fdr >= 0.05, NA, frac0_dif)) %>%
    filter(ptm %in% ptms9) %>% pull(frac0_dif) %>% abs %>% max(na.rm=T) -> maxfrac0_dif
  maxfrac0_dif
  qsum %>% 
    mutate(frac0_dif = ifelse(frac0_p_value_fdr >= 0.05, NA, frac0_dif)) %>%
    filter(ptm %in% ptms9) %>%
    ggplot(aes(x = dis, y = fct_rev(factor(ptm, levels = tmpptm)), fill = -frac0_dif)) +
    geom_tile() +
    scale_fill_gradientn(colours = c(scales::viridis_pal()(5)[1], scales::viridis_pal()(5)[2], "#FFFFFFFF", scales::viridis_pal()(5)[4], scales::viridis_pal()(5)[5]), values = scales::rescale(c(-maxfrac0_dif, -minfrac0_dif, 0, minfrac0_dif, maxfrac0_dif)), na.value = "white", breaks = c(-0.03, 0, 0.03), labels = c("-0.03", "0", "0.03"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = "Variable residues") +
    theme_nature(legend_position = "bottom", legend_nudge_top=-14, legend_nudge_right=-0.25) +
    # theme(text = element_text(family = "Helvetica Neue", size = 5)) + # Base font size
    theme(legend.text = element_text(size = 5, margin = margin(t = 1))) +
    # theme(legend.text = element_text(size = 6, margin = margin(t = 1))) +
    theme(axis.text.y = element_text(size = 6)) +
    theme(legend.title.position = "left", legend.title = element_text(size = 5, vjust = 1, margin = margin(t = unit(0.75, "mm"), r = unit(1, "mm")))) +
    xlab(NULL) +
    ylab(NULL)
  qsave(f("output-snps_gnomad-frac0_dif-vertical-filtered_fdr-select9{tmp_plddt}.pdf"), height = 50)
  
  # AC0 geom_tile horizontal filtered (FDR) favptms (select11)
  # tmpptm <- qsum %>% arrange(ptm) %>% pull(ptm) %>% unique
  # tmpptm <- sort(favptms)
  tmpptm <- favptms
  # tmpptm
  # Get minfrac0_dif (lowest absolute frac0_dif where frac0_p_value_fdr is significant (to make sure it falls into what we plot as white on the fill scale))
  qsum %>%
    mutate(frac0_dif = ifelse(frac0_p_value_fdr >= 0.05, NA, frac0_dif)) %>%
    filter(ptm %in% favptms) %>% pull(frac0_dif) %>% abs %>% min(na.rm=T) -> minfrac0_dif
  # Get maxfrac0_dif (highest absolute value for the colour scale, which will be symmetric around 0)
  qsum %>%
    mutate(frac0_dif = ifelse(frac0_p_value_fdr >= 0.05, NA, frac0_dif)) %>%
    filter(ptm %in% favptms) %>% pull(frac0_dif) %>% abs %>% max(na.rm=T) -> maxfrac0_dif
  maxfrac0_dif
  qsum %>% 
    mutate(frac0_dif = ifelse(frac0_p_value_fdr >= 0.05, NA, frac0_dif)) %>%
    filter(ptm %in% favptms) %>%
    ggplot(aes(x = factor(ptm, levels = tmpptm), y = dis, fill = -frac0_dif)) +
    geom_tile() +
    # geom_text(aes(label = ifelse(frac0_dif < 0, "+", "")), size = 2) +      # "+" means PTM has more variants than Control
    # scale_fill_viridis_c(na.value = "white", breaks = pretty_breaks(2), name = NULL) +
    # scale_fill_viridis_c(na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = NULL) +
    # scale_fill_gradient2(low = scales::viridis_pal()(5)[2], mid = "white", high = scales::viridis_pal()(5)[4], na.value = "white") +
    # scale_fill_gradient2(low = scales::viridis_pal()(5)[2], mid = "white", high = scales::viridis_pal()(5)[4], na.value = "white") +
    # scale_fill_gradient2(low = scales::viridis_pal()(9)[4], mid = "white", high = scales::viridis_pal()(9)[8], na.value = "white") +
    # scale_fill_gradient2(low = scales::viridis_pal()(9)[4], mid = "white", high = scales::viridis_pal()(9)[8], na.value = "white") +
    # scale_fill_gradientn(colours=c("blue", "black", "yellow"), na.value="white", limits = c(-maxfrac0_dif, maxfrac0_dif), breaks=pretty_breaks()) +
    # scale_fill_gradientn(colours = c("#440154", "#414487", "#FFFFFF", "#2A788E", "#3CBB75"), values = scales::rescale(c(-1, -0.5, 0, 0.5, 1)), na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = NULL) +
    # scale_fill_gradient2(low = "#440154", mid = "white", high = "#2A788E", midpoint = 0, na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = NULL) +
    # scale_fill_gradientn(colours = c("#440154", "#414487", "#FFFFFF", "#2A788E", "#3CBB75", "#20908C"), values = scales::rescale(c(-1, -0.5, 0, 0.25, 0.5, 1)), na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = NULL) +
    # scale_fill_gradientn(colours = c("#440154", "#414487", "#482878", "#3E4989", "#31688E", "#FFFFFF", "#35B779", "#6DCD59", "#B4DE2C", "#FDE725"), values = scales::rescale(c(-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.7, 1)), na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = NULL) +
    # scale_fill_gradientn(colours = c("#440154", "#414487", "#2A788E", "#2A788E", "#2A788E", "#7AD151", "#FDE725"), values = scales::rescale(c(-1, -0.6, -0.1, 0, 0.1, 0.6, 1)), na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = NULL) +
    # scale_fill_gradientn(colours = c("#440154", "#414487", "#2A788E", "#FFFFFF", "#2A788E", "#7AD151", "#FDE725"), values = scales::rescale(c(-1, -0.6, -0.1, 0, 0.1, 0.6, 1)), na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = NULL) +
    # scale_fill_gradientn(colours = c("#440154", "#414487", "#2A788E", "#FFFFFF", "#FFFFFF", "#2A788E", "#7AD151", "#FDE725"), values = scales::rescale(c(-1, -0.6, -0.15, -0.05, 0.05, 0.15, 0.6, 1)), na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = NULL) +
    # scale_fill_gradientn(colours = c("#440154", "#414487", "#482878", "#3E4989", "#31688E", "#FFFFFF", "#35B779", "#6DCD59", "#B4DE2C", "#FDE725"), values = scales::rescale(c(-1, -0.7, -0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5, 0.7)), na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = NULL) +
    # scale_fill_gradientn(colours = c("#440154", "#414487", "#482878", "#3E4989", "#FFFFFF", "#35B779", "#6DCD59", "#B4DE2C", "#FDE725"), values = scales::rescale(c(-1, -0.7, -0.5, -0.2, 0, 0.2, 0.5, 0.7, 1)), na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = NULL) +
    # scale_fill_gradientn(colours = c("#440154", "#414487", "#482878", "#3E4989", "#31688E", "#FFFFFF", "#35B779", "#6DCD59", "#B4DE2C", "#FDE725"), values = scales::rescale(c(-1, -0.8, -0.6, -0.3, -0.1, 0, 0.1, 0.3, 0.6, 0.8)), na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = NULL)
    # scale_fill_gradientn(colours = c("#440154", "#414487", "#482878", "#3E4989", "#31688E", "#FFFFFF", "#35B779", "#6DCD59", "#B4DE2C", "#FDE725", "#FDE725"), values = scales::rescale(c(-1, -0.8, -0.6, -0.3, -0.1, 0, 0.1, 0.3, 0.6, 0.8, 1)), na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = NULL)
    # scale_fill_gradientn(colours = c("#440154", "#443A83", "#31688E", "#21918C", "#35B779", "#FFFFFF", "#35B779", "#21918C", "#31688E", "#443A83", "#440154"), values = scales::rescale(c(-1, -0.8, -0.6, -0.3, -0.1, 0, 0.1, 0.3, 0.6, 0.8, 1)), na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = NULL)
    # scale_fill_gradientn(colours = c("#440154", "#443983", "#3B528B", "#21908C", "#2A788E", "#FFFFFF", "#35B779", "#5DC963", "#8FD744", "#B8DE29", "#FDE725"), values = scales::rescale(c(-1, -0.8, -0.6, -0.3, -0.1, 0, 0.1, 0.3, 0.6, 0.8, 1)), na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = NULL)
    # scale_fill_gradientn(colours = c("#440154", "#443983", "#3B528B", "#2A788E", "#21908C", "#FFFFFF", "#35B779", "#5DC963", "#8FD744", "#B8DE29", "#FDE725"), values = scales::rescale(c(-1, -0.8, -0.6, -0.3, -0.1, 0, 0.1, 0.3, 0.6, 0.8, 1)), na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = NULL) +
    # scale_fill_gradientn(colours = c(scales::viridis_pal()(5)[1], scales::viridis_pal()(5)[2], "#FFFFFFFF", scales::viridis_pal()(5)[4], scales::viridis_pal()(5)[5]), values = scales::rescale(c(-1, -0.1, 0, 0.1, 1)), na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = NULL) +
    # scale_fill_gradientn(colours = c(scales::viridis_pal()(5)[1], scales::viridis_pal()(5)[2], "#FFFFFFFF", "#FFFFFFFF", "#FFFFFFFF", scales::viridis_pal()(5)[4], scales::viridis_pal()(5)[5]), values = scales::rescale(c(-maxfrac0_dif, -minfrac0_dif, -minfrac0_dif, 0, minfrac0_dif, minfrac0_dif, maxfrac0_dif)), na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = NULL) +
    # scale_fill_gradientn(colours = c(scales::viridis_pal()(5)[1], scales::viridis_pal()(5)[2], scales::viridis_pal()(5)[2], "#FFFFFFFF", "#FFFFFFFF", "#FFFFFFFF", scales::viridis_pal()(5)[4], scales::viridis_pal()(5)[4], scales::viridis_pal()(5)[5]), values = scales::rescale(c(-maxfrac0_dif, -minfrac0_dif, -minfrac0_dif, -minfrac0_dif, 0, minfrac0_dif, minfrac0_dif, minfrac0_dif, maxfrac0_dif)), na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = NULL) +
    # scale_fill_gradientn(colours = c(scales::viridis_pal()(7)[1], scales::viridis_pal()(7)[2], scales::viridis_pal()(7)[3], "#FFFFFFFF", scales::viridis_pal()(7)[5], scales::viridis_pal()(7)[6], scales::viridis_pal()(7)[7]), values = scales::rescale(c(-maxfrac0_dif, -minfrac0_dif, 0, minfrac0_dif, maxfrac0_dif)), na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = NULL) +
    # scale_fill_gradientn(colours = c(scales::viridis_pal()(7)[1], scales::viridis_pal()(7)[2], scales::viridis_pal()(7)[3], "#FFFFFFFF", scales::viridis_pal()(7)[5], scales::viridis_pal()(7)[6], scales::viridis_pal()(7)[7]), values = scales::rescale(c(-1, -0.5, -0.1, 0, 0.1, 0.5, 1)), na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = NULL) +
    scale_fill_gradientn(colours = c(scales::viridis_pal()(5)[1], scales::viridis_pal()(5)[2], "#FFFFFFFF", scales::viridis_pal()(5)[4], scales::viridis_pal()(5)[5]), values = scales::rescale(c(-maxfrac0_dif, -minfrac0_dif, 0, minfrac0_dif, maxfrac0_dif)), na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = NULL) +
    # theme_nature() +
    theme_nature(legend_position = "top") +
    theme(text = element_text(family = "Helvetica Neue", size = 5)) + # Base font size
    theme(legend.text = element_text(size = 5, margin = margin(t = 1))) +
    # legend.key.spacing.y = unit(0, "mm"),
    # theme(legend.position = "none") +
    # ggtitle("") +
    xlab(NULL) +
    # xlab("") +
    ylab(NULL)
  qsave(f("output-snps_gnomad-frac0_dif-horizontal-filtered_fdr-favptms{tmp_plddt}.pdf"), width = 80, height = 30)
  # >> High frac0_dif means more invariant sites in PTM than Control. An example is M-ac: the initiator methionine is of course usually invariant.
  
  # # Viridis colours (uncomment for RStudio to display them):
  # 2 "#440154FF" "#FDE725FF"
  # 3 "#440154FF" "#21908CFF" "#FDE725FF"
  # 4 "#440154FF" "#31688EFF" "#35B779FF" "#FDE725FF"
  # 4             "#31688EFF"             "#FDE725FF" 2 and 4
  # 5 "#440154FF" "#3B528BFF" "#21908CFF" "#5DC863FF" "#FDE725FF"
  # 5             "#3B528BFF"             "#5DC863FF" 2 and 4
  # 6 "#440154FF" "#414487FF" "#2A788EFF" "#22A884FF" "#7AD151FF" "#FDE725FF"
  # 7 "#440154FF" "#443A83FF" "#31688EFF" "#21908CFF" "#35B779FF" "#8FD744FF" "#FDE725FF"
  # 8 "#440154FF" "#46337EFF" "#365C8DFF" "#277F8EFF" "#1FA187FF" "#4AC16DFF" "#9FDA3AFF" "#FDE725FF"
  # 9 "#440154FF" "#472D7BFF" "#3B528BFF" "#2C728EFF" "#21908CFF" "#27AD81FF" "#5DC863FF" "#AADC32FF" "#FDE725FF"
  # 9             "#472D7BFF"                                                             "#AADC32FF" 2 and 8
  # 9                                     "#2C728EFF"                                     "#AADC32FF" 4 and 8
  # 9                                                 "#21908CFF"                         "#AADC32FF" 5 and 8
  # See also: blang.R for ptmcol1 and ptmcol0, which are scales::viridis_pal()(5) [2] and [4]
  
  
  
  
  
  
  
  
  # Variant plots
  # plot_variants_bin2d("a", ptms11, 36, "Common natural variants", "3a-snps_gnomad-common_variants", snps_gnomad %>% filter(af > 0.0001))
  mytag <- "a"
  mpptm <- ptms11
  my_width <- 36
  mytitle <- "Common natural variants"
  myfilename <- "3a-snps_gnomad-common_variants"
  snps_tmp <- snps_gnomad %>% filter(af > 0.0001)
  show_insig <- F
  plot_variants_bin2d <- function(mytag, tmpptm, my_width, mytitle, myfilename, snps_tmp, show_insig = F) {
    
    # Combine
    qsnp <- left_join(q0, snps_tmp, relationship = "many-to-many", by = join_by(acc, site))
    
    # accs_with_ptms <- q %>% filter(type == "ptm") %>% pull(acc) %>% unique
    # accs_with_controls <- q %>% filter(type == "control") %>% pull(acc) %>% unique
    # # Get intersection
    # accs_with_ptms_and_controls <- intersect(accs_with_ptms, accs_with_controls)
    # qsnp %<>% filter(acc %in% accs_with_ptms_and_controls)
    # qsnp
    
    # To ensure that these residues can have variants: only include uniens accs (Ensembl 108, and I used VEP 108)
    print(f("Before uniens filter: {comma(qsnp %>% pull(acc) %>% unique %>% length)} accs"))
    qsnp %<>% filter(acc %in% uniens_accs)
    print(f("After uniens filter:  {comma(qsnp %>% pull(acc) %>% unique %>% length)} accs"))
    
    # Set ac to 0 if NA
    qsnp %<>% mutate(ac = ifelse(is.na(ac), 0, ac))
    
    # Using ptms1000 instead for a consistent set of PTMs with a reasonable number of known sites
    # # Ensure we have all 6 combinations of modified/control and buried/structured/surface for each PTM type
    # myptms <- q %>% group_by(ptm) %>% summarise(combinations = n_distinct(dis, type)) %>% filter(combinations == 6) %>% pull(ptm) %>% unique
    # myptms
    # length(myptms)
    
    # qsum frac0
    qsnp %<>% 
      # filter(ptm %in% myptms) %>%
      filter(ptm %in% ptms1000) %>%
      nest_by(ptm, dis) %>% 
      mutate(
        n0_ptm = data %>% filter(type == "ptm" & ac == 0) %>% nrow,
        n0_control = data %>% filter(type == "control" & ac == 0) %>% nrow,
        n_ptm = data %>% filter(type == "ptm") %>% nrow,
        n_control = data %>% filter(type == "control") %>% nrow,
        frac0_ptm = n0_ptm / n_ptm,
        frac0_control = n0_control / n_control,
        frac0_dif = frac0_ptm - frac0_control,
        frac0_fisher = fisher.test(matrix(c(n0_ptm, n0_control, n_ptm - n0_ptm, n_control - n0_control), nrow = 2))$p.value
      ) %>%
      summarise(ptm, dis, frac0_ptm, frac0_control, frac0_dif, frac0_p_value = frac0_fisher, .groups = "drop")
    
    # Replace NaN with NA (e.g. for C-pal disordered)
    qsnp %<>% filter(!is.na(frac0_dif))
    # qsnp %<>% mutate(frac0_ptm = ifelse(is.nan(frac0_ptm), NA, frac0_ptm))
    # qsnp %<>% mutate(frac0_control = ifelse(is.nan(frac0_control), NA, frac0_control))
    # qsnp %<>% mutate(frac0_dif = ifelse(is.nan(frac0_dif), NA, frac0_dif))
    
    # Apply FDR correction to p-values
    qsnp$frac0_p_value_fdr <- p.adjust(qsnp$frac0_p_value, method = "BH")
    qsnp
    
    # Adjust legend position for number of PTMs (ptms11, ptms17=ptms1000)
    if (identical(tmpptm, ptms9)) {
      my_ptms <- "select9"
      # my_legend_nudge_top <- -8.25
      my_legend_nudge_top <- 0
      my_height <- 60
      # my_theme <- theme()
      levels(qsnp$dis) <- c("Buried", "Structured  ", "  Disordered")
    } else if (identical(tmpptm, ptms11)) {
      my_ptms <- "select11"
      # my_legend_nudge_top <- -8.25
      my_legend_nudge_top <- 0
      my_height <- 70
      # my_theme <- theme()
      levels(qsnp$dis) <- c("Buried", "Structured  ", "  Disordered")
    } else if (identical(tmpptm, ptms1000)) {
      my_ptms <- "ptms1000"
      # my_legend_nudge_top <- -8.25
      my_legend_nudge_top <- 0
      my_height <- 80
      # my_theme <- theme(legend.position = "bottom")
      # my_theme <- theme()
      levels(qsnp$dis) <- c("Buried", "Structured  ", "  Disordered")
    }
    
    tmp_insig <- ""
    if (show_insig == F) {
      # Filter out insignificant comparisons
      qsnp %<>% mutate(frac0_dif = ifelse(frac0_p_value_fdr >= 0.05, NA, frac0_dif))
    } else {
      tmp_insig <- "-insig"
    }
    
    print(qsnp %>% filter(frac0_p_value_fdr < 0.05), n=200)
    
    # Check if there is anything significant before plotting
    if (qsnp %>% filter(frac0_p_value_fdr < 0.05) %>% nrow > 0) {
      # AC0 geom_tile vertical filtered (FDR) favptms (select11)
      # Get minfrac0_dif (lowest absolute frac0_dif where p_value_fdr is significant (to make sure it falls into what we plot as white on the fill scale))
      qsnp %>% 
        mutate(frac0_dif = ifelse(frac0_p_value_fdr >= 0.05, NA, frac0_dif)) %>%
        filter(ptm %in% tmpptm) %>% pull(frac0_dif) %>% abs %>% min(na.rm=T) -> minfrac0_dif
      # Get maxfrac0_dif (highest absolute value for the colour scale, which will be symmetric around 0)
      qsnp %>% 
        mutate(frac0_dif = ifelse(frac0_p_value_fdr >= 0.05, NA, frac0_dif)) %>%
        filter(ptm %in% tmpptm) %>% pull(frac0_dif) %>% abs %>% max(na.rm=T) -> maxfrac0_dif
      maxfrac0_dif
      # Define break (e.g. maxfrac0_dif 0.056 will lead to 0.05 -> -0.05, 0, 0.05)
      # mybreak <- round(maxfrac0_dif - 0.005, 2)
      if (maxfrac0_dif >= 0.01) {
        mybreak <- round(maxfrac0_dif - 0.005, 2)
      } else {
        mybreak <- round(maxfrac0_dif - 0.0005, 3)
      }
      # print(maxfrac0_dif)
      # Add spacing between Structured and Disordered by renaming level "Structured" to "Structured " and Disordered to " Disordered"" 
      qsnp %>% 
        # mutate(frac0_dif = ifelse(frac0_p_value_fdr >= 0.05, NA, frac0_dif)) %>%
        filter(ptm %in% tmpptm) %>%
        ggplot(aes(x = dis, y = fct_rev(factor(ptm, levels = tmpptm)), fill = -frac0_dif)) +
        geom_tile() +
        # scale_fill_gradientn(colours = c(scales::viridis_pal()(5)[1], scales::viridis_pal()(5)[2], "#FFFFFFFF", scales::viridis_pal()(5)[4], scales::viridis_pal()(5)[5]), values = scales::rescale(c(-maxfrac0_dif, -minfrac0_dif, 0, minfrac0_dif, maxfrac0_dif)), na.value = "white", breaks = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"), limits = c(-maxfrac0_dif, maxfrac0_dif), name = "Variable residues") +
        scale_fill_gradientn(colours = c(scales::viridis_pal()(5)[1], scales::viridis_pal()(5)[2], "#FFFFFFFF", scales::viridis_pal()(5)[4], scales::viridis_pal()(5)[5]), values = scales::rescale(c(-maxfrac0_dif, -minfrac0_dif, 0, minfrac0_dif, maxfrac0_dif)), na.value = "white", breaks = c(-mybreak, 0, mybreak), labels = c(f("-{mybreak}"), "0", f("{mybreak}")), limits = c(-maxfrac0_dif, maxfrac0_dif), name = "Variable residues") +
        theme_nature(legend_position = "bottom", legend_nudge_top=my_legend_nudge_top, legend_nudge_right=-2, extra_margin_bottom=2) +
        # theme(text = element_text(family = "Helvetica Neue", size = 5)) + # Base font size
        theme(legend.text = element_text(size = 5, margin = margin(t = 1))) +
        # theme(legend.text = element_text(size = 6, margin = margin(t = 1))) +
        theme(axis.text.y = element_text(size = 6)) +
        theme(legend.title.position = "left", legend.title = element_text(size = 5, vjust = 1, margin = margin(t = unit(0.75, "mm"), r = unit(1, "mm")))) +
        # theme(legend.title.position = "top", legend.title = element_text(size = 5, hjust = 1, vjust = 1, margin = margin(b = unit(0.75, "mm"), t = unit(1.5, "mm"), r = unit(1, "mm")))) +
        # my_theme +
        labs(title = mytitle) +
        labs(tag = mytag) +
        xlab("") +
        ylab(NULL)
      qsave(f("output-variants{tmp_insig}-{myfilename}-frac0_dif-vertical-filtered_fdr-{my_ptms}{tmp_plddt}.pdf"), width = my_width, height = my_height)
      
      # Save qsnp to TSV
      write_tsv(qsnp, f("output-tsv-variants{tmp_insig}-{myfilename}-frac0_dif-vertical-filtered_fdr-{my_ptms}{tmp_plddt}.tsv"))
    }
  }
  
  if (plddt_filtering == 0) {
    mywidth <- 36
  } else {
    # pLDDT filtering on: Plot without disorder (narrower)
    mywidth <- 36
  }
    
  # # # Figure 3
  # # plot_variants_bin2d("a", ptms11, 40, "Common natural variants", "3a-snps_gnomad-common_variants", snps_gnomad %>% filter(af > 0.0001))
  # # plot_variants_bin2d("b", ptms11, 40, "Rare natural variants", "3b-snps_gnomad-rare_variants", snps_gnomad %>% filter(af <= 0.0001))
  # # plot_variants_bin2d("c", ptms11, 40, "Clinically significant variants", "3c-snps_clinvar-clinsig", snps_clinvar %>% filter(clinsig == 1))
  # # plot_variants_bin2d("d", ptms11, 40, "Cancer Gene Census mutations", "3d-snps_cosmic-cgc1", snps_cosmic %>% filter(cgc == 1))
  # 
  # # Figure 3 (slightly narrower to fit 180 = 36 * 5 plots)
  # plot_variants_bin2d("a", ptms11, 36, "Common natural variants", "3a-snps_gnomad-common_variants", snps_gnomad %>% filter(af > 0.0001))
  # plot_variants_bin2d("b", ptms11, 36, "Rare natural variants", "3b-snps_gnomad-rare_variants", snps_gnomad %>% filter(af <= 0.0001))
  # plot_variants_bin2d("c", ptms11, 36, "Recurrent cancer mutations", "3c-snps_cosmic-recurrent-ac-gte10", snps_cosmic %>% filter(ac >= 10))
  # plot_variants_bin2d("d", ptms11, 36, "Cancer Gene Census", "3d-snps_cosmic-cgc1", snps_cosmic %>% filter(cgc == 1))
  # plot_variants_bin2d("e", ptms11, 36, "Clinically significant variants", "3e-snps_clinvar-clinsig", snps_clinvar %>% filter(clinsig == 1))
  # 
  # # Figure 3 (ptms9) (slightly narrower to fit 160 = 32 * 5 plots)
  # plot_variants_bin2d("a", ptms9, 32, "Common variants", "ptms9-3a-snps_gnomad-common_variants", snps_gnomad %>% filter(af > 0.0001))
  # plot_variants_bin2d("b", ptms9, 32, "Rare variants", "ptms9-3b-snps_gnomad-rare_variants", snps_gnomad %>% filter(af <= 0.0001))
  # plot_variants_bin2d("c", ptms9, 32, "Recurrent in cancer", "ptms9-3c-snps_cosmic-recurrent-ac-gte10", snps_cosmic %>% filter(ac >= 10))
  # plot_variants_bin2d("d", ptms9, 32, "Cancer Gene Census", "ptms9-3d-snps_cosmic-cgc1", snps_cosmic %>% filter(cgc == 1))
  # plot_variants_bin2d("e", ptms9, 32, "Clinically significant", "ptms9-3e-snps_clinvar-clinsig", snps_clinvar %>% filter(clinsig == 1))
  # 
  # # Figure E3 (Extended Figure 3 with ptms17 instead of ptms11)
  # plot_variants_bin2d("a", ptms1000, 36, "Common natural variants", "E3a-snps_gnomad-common_variants", snps_gnomad %>% filter(af > 0.0001))
  # plot_variants_bin2d("b", ptms1000, 36, "Rare natural variants", "E3b-snps_gnomad-rare_variants", snps_gnomad %>% filter(af <= 0.0001))
  # plot_variants_bin2d("c", ptms1000, 36, "Recurrent cancer mutations", "E3c-snps_cosmic-recurrent-ac-gte10", snps_cosmic %>% filter(ac >= 10))
  # plot_variants_bin2d("d", ptms1000, 36, "Cancer Gene Census", "E3d-snps_cosmic-cgc1", snps_cosmic %>% filter(cgc == 1))
  # plot_variants_bin2d("e", ptms1000, 36, "Clinically significant variants", "E3e-snps_clinvar-clinsig", snps_clinvar %>% filter(clinsig == 1))
  # 
  # # Figure E3 (Extended Figure 3 with ptms17 instead of ptms11), and showing insignificant comparisons (not used, too confusing/useless)
  # plot_variants_bin2d("a", ptms1000, 36, "Common natural variants", "E3a-snps_gnomad-common_variants", snps_gnomad %>% filter(af > 0.0001), show_insig = T)
  # plot_variants_bin2d("b", ptms1000, 36, "Rare natural variants", "E3b-snps_gnomad-rare_variants", snps_gnomad %>% filter(af <= 0.0001), show_insig = T)
  # plot_variants_bin2d("c", ptms1000, 36, "Recurrent cancer mutations", "E3c-snps_cosmic-recurrent-ac-gte10", snps_cosmic %>% filter(ac >= 10), show_insig = T)
  # plot_variants_bin2d("d", ptms1000, 36, "Cancer Gene Census", "E3d-snps_cosmic-cgc1", snps_cosmic %>% filter(cgc == 1), show_insig = T)
  # plot_variants_bin2d("e", ptms1000, 36, "Clinically significant variants", "E3e-snps_clinvar-clinsig", snps_clinvar %>% filter(clinsig == 1), show_insig = T)
  # 
  # # Figure S3 (Larger sets of variants: "all" etc.)
  # plot_variants_bin2d("a", ptms1000, 40, "gnomAD (all)", "S3a-snps_gnomad-all_variants", snps_gnomad)
  # plot_variants_bin2d("b", ptms1000, 40, "ClinVar (all)", "S3b-snps_clinvar-all_variants", snps_clinvar)
  # plot_variants_bin2d("c", ptms1000, 40, "COSMIC (all)", "S3c-snps_cosmic-all_variants", snps_cosmic)
  # plot_variants_bin2d("d", ptms1000, 40, "COSMIC CGC (all)", "S3d-snps_cosmic-cgc", snps_cosmic %>% filter(cgc != 0))
  # plot_variants_bin2d("e", ptms1000, 40, "COSMIC CGC tier 2 (candidate)", "S3e-snps_cosmic-cgc2", snps_cosmic %>% filter(cgc == 2))
  # plot_variants_bin2d("f", ptms1000, 40, "COSMIC ≥10 samples", "S3f-snps_cosmic-recurrent-ac-gte10", snps_cosmic %>% filter(ac >= 10))
  # plot_variants_bin2d("g", ptms1000, 40, "Genomenon (all)", "S3g-snps_mastermind-all_variants", snps_mastermind)
  # plot_variants_bin2d("h", ptms1000, 40, "Genomenon ≥10 articles", "S3h-snps_mastermind-papers-gte10", snps_mastermind %>% filter(ac >= 10))
  # # plot_variants_bin2d("h", ptms1000, 40, "Genomenon ≥2 articles", "snps_mastermind-papers-gte2", snps_mastermind %>% filter(ac >= 2))
  # # plot_variants_bin2d("h", ptms1000, 40, "Genomenon ≥3 articles", "snps_mastermind-papers-gte3", snps_mastermind %>% filter(ac >= 3))
  # # plot_variants_bin2d("h", ptms1000, 40, "Genomenon ≥4 articles", "snps_mastermind-papers-gte4", snps_mastermind %>% filter(ac >= 4))
  # # plot_variants_bin2d("h", ptms1000, 40, "Genomenon ≥5 articles", "snps_mastermind-papers-gte5", snps_mastermind %>% filter(ac >= 5))
  # # plot_variants_bin2d(snps_pcgp, "snps_pcgp-all_variants")
  # # plot_variants_bin2d(snps_pcgp_valid, "snps_pcgp-valid")
  
  
  
  
  
  
  
  
  
  # Mimic SNPs
  q0
  snps_gnomad
  # Combine
  qsnp <- left_join(q0, snps_gnomad, relationship = "many-to-many")
  qsnp
  qsnp %>% pull(af) %>% summary
  qsnp %>% pull(af) %>% is.na %>% summary
  # Percentage that has af=NA (i.e. no variants exist for these acc|sites)
  qsnp %>% filter(is.na(af)) %>% nrow / nrow(q) * 100
  # >> 52.6% of residues don't have any variants.
  # Set original and variant residues for invariant positions (e.g. to K and K for K-ac)
  qsnp %<>% mutate(original = ifelse(is.na(af), aa, original), variant = ifelse(is.na(af), aa, variant))
  qsnp
  # Replace AF=NA with 0
  qsnp %<>% mutate(af = ifelse(is.na(af), 0, af))
  # # Filter out AF=NA
  # qsnp %<>% filter(!is.na(af))
  # # qsnp %<>% filter(af != 0)
  qsnp$af %>% summary
  qsnp
  qsnp %>% filter(ptm == "Y-p" & variant == "D")
  # qsnp %>% filter(ptm %in% mainptms) %>% group_by(ptm, dis, type, original, variant) %>% tally
  # qsnp %>% filter(ptm %in% mainptms) %>% group_by(ptm, dis, type, original, variant) %>% summarise(n = n())
  # identical(qsnp %>% filter(ptm %in% mainptms) %>% group_by(ptm, dis, type, original, variant) %>% tally, qsnp %>% filter(ptm %in% mainptms) %>% group_by(ptm, dis, type, original, variant) %>% summarise(n = n()))
  qsnp %>% filter(ptm %in% mainptms) %>% group_by(ptm, dis, type, original, variant) %>% tally
  qsnp %>% filter(ptm %in% mainptms) %>% group_by(ptm, dis, type, original, variant) %>% tally %>% mutate(frac = n / sum(n))
  { qsnp %>% 
      filter(ptm %in% mainptms) %>% 
      group_by(ptm, dis, type, original, variant) %>% 
      tally %>%
      mutate(frac = n / sum(n)) %>%
      ggplot(aes(x = variant, y = frac, fill = fct_rev(type))) +
      geom_col(position = "dodge") +
      # facet_wrap(~ ptm, scales = "free") +
      facet_grid(cols = vars(ptm), rows = vars(fct_rev(type)), scales = "free") +
      # scale_fill_viridis_d() +
      # scale_fill_manual(values = c("ptm" = ptmvir1, "control" = ptmvir0)) +
      scale_fill_manual(values = c("ptm" = scales::viridis_pal()(5)[4], "control" = scales::viridis_pal()(5)[2])) +
      theme_nature()
  } %>% qsave(f("output-mimic-variants-letters-snps_gnomad{tmp_plddt}.pdf"), width = 120)
  
  # Enrichment PTM vs Control
  qsnp %>% filter(ptm %in% mainptms) %>% group_by(ptm, dis, type, original, variant) %>% tally %>% mutate(frac = n / sum(n))
  # qsnp %>% filter(ptm %in% mainptms) %>% group_by(ptm, dis, type, original, variant) %>% tally %>% mutate(frac = n / sum(n)) %>% filter(is.na(n))
  # qsnp %>% filter(ptm %in% mainptms) %>% group_by(ptm, dis, type, original, variant) %>% tally %>% mutate(frac = n / sum(n)) %>% filter(is.na(frac))
  # qsnp %>% filter(ptm %in% mainptms) %>% group_by(ptm, dis, type, original, variant) %>% tally %>% mutate(frac = n / sum(n)) %>% select(-n) %>% pivot_wider(names_from = type, values_from = frac, names_prefix = "frac_") %>% print(n=1000)
  # qsnp %>% filter(ptm %in% mainptms) %>% group_by(ptm, dis, type, original, variant) %>% tally %>% mutate(frac = n / sum(n)) %>% select(-n) %>% pivot_wider(names_from = type, values_from = frac, names_prefix = "frac_") %>% pivot_longer(cols = starts_with("frac_"), names_to = "type", values_to = "frac", names_prefix = "frac_")
  qsnp %>% filter(ptm %in% mainptms) %>% group_by(ptm, dis, type, original, variant) %>% tally %>% mutate(frac = n / sum(n)) %>% select(-n) %>% pivot_wider(names_from = type, values_from = frac, names_prefix = "frac_") %>% mutate(enrichment = frac_ptm / frac_control)
  { qsnp %>% 
      filter(ptm %in% mainptms) %>% 
      group_by(ptm, dis, type, original, variant) %>% 
      tally %>% 
      mutate(frac = n / sum(n)) %>% 
      select(-n) %>% 
      pivot_wider(names_from = type, values_from = frac, names_prefix = "frac_") %>%
      # pivot_longer(cols = starts_with("frac_"), names_to = "type", values_to = "frac", names_prefix = "frac_") %>%
      mutate(enrichment = frac_ptm / frac_control) %>%
      ggplot(aes(x = variant, y = enrichment)) +
      # geom_col(position = "dodge") +
      # geom_point() +
      geom_text(aes(label = variant), size = 1) +
      geom_hline(yintercept = 1) +
      # facet_wrap(~ ptm, scales = "free") +
      # facet_grid(cols = vars(ptm), rows = vars(fct_rev(type)), scales = "free") +
      facet_grid(rows = vars(ptm), cols = vars(dis)) +
      # scale_fill_viridis_d() +
      # scale_fill_manual(values = c("ptm" = ptmvir1, "control" = ptmvir0)) +
      scale_fill_manual(values = c("ptm" = scales::viridis_pal()(5)[4], "control" = scales::viridis_pal()(5)[2])) +
      theme_nature()
  } %>% qsave(f("output-snps_gnomad-mimic-variants{tmp_plddt}.pdf"), width = 120, height = 60)
  

  # Mimic bar plot
  # # plot_mimic_variants_bars("c", ptms11, "Clinically significant variants", "4c-snps_clinvar-clinsig", snps_clinvar %>% filter(clinsig == 1))
  # mytag <- "c"
  # tmpptm <- ptms11
  # mytitle <- "Clinically significant variants"
  # myfilename <- "4c-snps_clinvar-clinsig"
  # snps_tmp <- snps_clinvar %>% filter(clinsig == 1)
  # # plot_mimic_variants_bars("b", ptms11, "Rare natural variants", "4b-snps_gnomad-rare_variants", snps_gnomad %>% filter(af <= 0.0001))
  # mytag <- "b"
  # tmpptm <- ptms11
  # mytitle <- "Rare natural variants"
  # myfilename <- "4b-snps_gnomad-rare_variants"
  # snps_tmp <- snps_gnomad %>% filter(af <= 0.0001)
  # disfilt <- "nodis"
  plot_mimic_variants_bars <- function(mytag, tmpptm, my_width, my_height, minsites, mytitle, myfilename, snps_tmp, disfilt = "", show_insig = F) {
    
    # Combine
    qsnp <- left_join(q0, snps_tmp, relationship = "many-to-many", by = join_by(acc, site))
    
    # To ensure that these residues can have variants: only include uniens accs (Ensembl 108, and I used VEP 108)
    print(f("Before uniens filter: {comma(qsnp %>% pull(acc) %>% unique %>% length)} accs"))
    qsnp %<>% filter(acc %in% uniens_accs)
    print(f("After uniens filter:  {comma(qsnp %>% pull(acc) %>% unique %>% length)} accs"))
    
    # Set original and variant residues for invariant positions (e.g. to K and K for K-ac)
    # Without setting these, the question would be: Where there is a variant, is there an enrichment for variant X at PTM sites vs. Control?
    # With setting these, the question is: Is there an enrichment for residue X at PTM sites vs. Control? Which is what we want.
    qsnp %<>% mutate(original = ifelse(is.na(ac), aa, original), variant = ifelse(is.na(ac), aa, variant))
    # # Replace AC=NA or AF=NA with 0
    # if ("ac" %in% colnames(qsnp))
    # {
    #   qsnp %<>% mutate(ac = ifelse(is.na(ac), 0, ac))
    # }
    # else
    # {
    #   qsnp %<>% mutate(af = ifelse(is.na(af), 0, af))
    # }
    # # qsnp %>% filter(original != variant)
    # qsnp %>% filter(original != variant) %>% summarise(nrow = nrow / nrow(q))
    # (qsnp %>% filter(original != variant) %>% nrow) / (nrow(q))
    # print(f("Percentage of residues with variants: {qsnp %>% filter(original != variant) %>% nrow / nrow(q) * 100}%"))
    # 
    qfilt <- qsnp %>% 
      # filter(ptm %in% tmpptm) %>% 
      filter(ptm %in% ptms1000) %>% 
      # # At least 3 individuals/samples/submitters/papers reporting a variant
      # filter(ac >= 3) %>% 
      group_by(ptm, dis, type, original, variant) %>% 
      tally %>% 
      mutate(frac = n / sum(n), total = sum(n)) %>% 
      pivot_wider(names_from = type, values_from = c(n, total, frac)) %>% 
      filter(!is.na(n_ptm)) %>% # Avoid some NAs
      filter(!is.na(n_control)) %>% # Avoid some NAs
      # Require at least 10 PTM and control sites per structural category
      filter(n_ptm >= minsites & n_control >= minsites) %>% # Avoid some NAs
      mutate(enrichment = frac_ptm / frac_control) %>%
      rowwise %>% 
      mutate(fisher = fisher.test(matrix(c(n_ptm, n_control, total_ptm - n_ptm, total_control - n_control), nrow = 2))$p.value) %>% 
      ungroup %>% # Turn off rowwise again
      mutate(log2enrich = log2(enrichment)) %>%
      # filter(fisher < 0.05) %>% # Significant only (raw p-values)
      mutate(fisher_fdr = p.adjust(fisher, method = "fdr")) %>%
      # filter(enrichment > 1) %>% # Positive only (enriched, not depleted)
      filter(original != variant) %>% # Variants only (don't report enrichment of original amino acid)
      # filter(fisher_fdr < 0.05) %>% # Significant only (FDR)
      filter(ptm %in% tmpptm) %>% # PTMs being plotted only
      arrange(fisher_fdr)
    if (show_insig == F) {
      qfilt <- qfilt %>% filter(fisher_fdr < 0.05)
    }
    qfilt
    
    # disfilt
    tmpdisfilt <- ""
    if (disfilt == "nodis") {
      qfilt %<>% filter(dis != "Disordered")
      qfilt$dis %<>% droplevels
      qfilt$dis
      # # Always show Buried and Structured
      # # qfilt$dis <- factor(qfilt$dis, levels = c("Buried", "Structured"))
      # levels(qfilt$dis) <- c("Buried", "Structured")
      tmpdisfilt <- "-nodis"
    }
    
    # Filter for e.g. Phosphorylation only, based on filename
    qfilt$ptm <- factor(qfilt$ptm, levels = tmpptm)
    if (str_detect(myfilename, "STY-p")) {
      # Phosphorylation only
      qfilt %<>% filter(ptm %in% c("S-p", "T-p", "Y-p"))
      qfilt$ptm %<>% droplevels
    }
    if (str_detect(myfilename, "K-ub-sum-mal")) {
      # K-ub, K-sum and K-mal only (these show interesting trends)
      qfilt %<>% filter(ptm %in% c("K-ub", "K-sum", "K-mal"))
      qfilt$ptm %<>% droplevels
    }

    # Get minenrich (lowest absolute enrichment_dif where p_value_fdr is significant (to make sure it falls into what we plot as white on the fill scale))
    qfilt %>% pull(log2enrich) %>% abs %>% min(na.rm=T) -> minenrich
    minenrich <- 0.1
    minenrich
    # Get maxenrich (highest absolute value for the colour scale, which will be symmetric around 0)
    qfilt %>% pull(log2enrich) %>% abs %>% max(na.rm=T) -> maxenrich
    maxenrich
    # Define break (e.g. maxenrich 0.056 will lead to 0.05 -> -0.05, 0, 0.05)
    # mybreak <- round(maxenrich - 0.005, 2)
    # mybreak <- round(maxenrich - 0.05, 1)
    mybreak <- round(maxenrich - 0.1, 1)
    
    # Adjust legend position for number of PTMs (ptms11, ptms17=ptms1000)
    my_legend_nudge_top <- -6
    if (identical(tmpptm, ptms11)) {
      my_legend_nudge_top <- -6
      # my_height <- 80
    } else if (identical(tmpptm, ptms1000)) {
      my_legend_nudge_top <- -6
      # my_height <- 100
    } else if (identical(tmpptm, c("S-p", "T-p", "Y-p"))) {
      my_legend_nudge_top <- -6
      my_height <- 45
    }
    
    print(qfilt)
    print(f("{qfilt %>% nrow} rows"))
    
    # Only plot if any enrichments are available (otherwise plotting would crash)  
    if (qfilt %>% nrow > 0) {
      {
        qfilt %>%
          # mutate(enriched = ifelse(enrichment > 0, "enriched", "depleted")) %>%
          ggplot(aes(x = variant, y = log2enrich)) +
          # geom_col(aes(fill = enriched), position = "dodge") +
          geom_col(aes(fill = log2enrich), position = "dodge") +
          # geom_point() +
          # geom_text(aes(label = variant), size = 5/.pt, nudge_y = 1) +
          # geom_text(aes(label = variant, vjust = ifelse(log2enrich >= 0, -0.5, 1.5)), size = 5/.pt) +
          geom_text(aes(label = variant, vjust = ifelse(log2enrich >= 0, -0.1, 1.1)), size = 5/.pt) +
          # geom_text(aes(label = variant, vjust = ifelse(log2enrich >= 0, 0, 1)), size = 5/.pt) +
          geom_hline(yintercept = 0, linetype = "solid") +
          # coord_cartesian(ylim = c(-1, 1)) +
          # scale_y_continuous(expand = c(0.4, 0.4)) + # Expand more to make sure letters get plotted fully
          scale_y_continuous(expand = c(0.6, 0.6)) + # Expand more to make sure letters get plotted fully
          # facet_wrap(~ ptm, scales = "free") +
          # facet_grid(cols = vars(ptm), rows = vars(fct_rev(type)), scales = "free") +
          # facet_grid(rows = vars(as.factor(ptm)), cols = vars(as.factor(dis)), scales = "free", space = "free_x") +
          # facet_grid(rows = vars(factor(ptm, levels = tmpptm)), cols = vars(as.factor(dis)), scales = "free", space = "free_x") +
          # facet_grid(rows = vars(factor(ptm, levels = tmpptm)), cols = vars(as.factor(dis)), drop = T, scales = "free", space = "free_x", switch = "y") +
          facet_grid(rows = vars(factor(ptm, levels = tmpptm)), cols = vars(as.factor(dis)), drop = F, scales = "free", space = "free_x", switch = "y") +
          # scale_fill_viridis_d() +
          # scale_fill_manual(values = c("ptm" = ptmvir1, "control" = ptmvir0)) +
          # scale_fill_manual(values = c("enriched" = scales::viridis_pal()(5)[4], "depleted" = scales::viridis_pal()(5)[2])) +
          # scale_fill_viridis_c() +
          # scale_fill_gradientn(name = "Variable residues", colours = c(scales::viridis_pal()(5)[1], scales::viridis_pal()(5)[2], "#FFFFFFFF", scales::viridis_pal()(5)[4], scales::viridis_pal()(5)[5]), values = scales::rescale(c(-maxenrich, -minenrich, 0, minenrich, maxenrich)), na.value = "white", breaks = c(-mybreak, 0, mybreak), labels = c(f("-{mybreak}"), "0", f("{mybreak}")), limits = c(-maxenrich, maxenrich)) +
          scale_fill_gradientn(name = "Variant enrichment (log2)", colours = c(scales::viridis_pal()(5)[1], scales::viridis_pal()(5)[2], "#FFFFFFFF", scales::viridis_pal()(5)[4], scales::viridis_pal()(5)[5]), values = scales::rescale(c(-maxenrich, -minenrich, 0, minenrich, maxenrich)), na.value = "white", breaks = c(-mybreak, 0, mybreak), labels = c(f("-{mybreak}"), "0", f("{mybreak}")), limits = c(-maxenrich, maxenrich)) +
          theme_nature(legend_position = "bottom", legend_nudge_top=my_legend_nudge_top, legend_nudge_right=0, extra_margin_right=2) +
          theme(legend.text = element_text(size = 5, margin = margin(t = 1))) +
          # theme(legend.text = element_text(size = 6, margin = margin(t = 1))) +
          theme(axis.text.y = element_text(size = 6)) +
          # theme(legend.title.position = "left", legend.title = element_text(size = 5, vjust = 1, margin = margin(t = unit(0.75, "mm"), r = unit(1, "mm")))) +
          theme(legend.title.position = "top", legend.title = element_text(size = 5, hjust = 1, vjust = 1, margin = margin(b = unit(0.75, "mm"), t = unit(1.5, "mm"), r = unit(1, "mm")))) +
          theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            strip.placement = "outside",
            strip.text.y.left = element_text(angle = 0, hjust = 1),
            # strip.text.y.left = element_text(angle = 0, hjust = 1, margin = margin(r = 5)),
            # strip.switch.pad.grid = unit(0, "mm"),
            panel.spacing.y = unit(0, "mm"),
            strip.clip = "off",
          ) +
          labs(title = mytitle) +
          labs(tag = mytag) +
          xlab("") +
          ylab(NULL)
      } %>% qsave(f("output-mimic-variants-bars-filtered-fdr{tmpdisfilt}-{myfilename}-minsites{minsites}{tmp_plddt}.pdf"), width = my_width, height = my_height)
      
      # Save qfilt to TSV
      write_tsv(qfilt, f("output-mimic-variants-bars-tsv-filtered-fdr{tmpdisfilt}-{myfilename}-minsites{minsites}{tmp_plddt}.tsv"))
    }
  }
  
  # Figure 3 fghi (without showing "Disordered", to make these plots more compact)
  # plot_mimic_variants_bars("a", ptms11, 45, 10, "Common natural variants", "4a-snps_gnomad-common_variants", snps_gnomad %>% filter(af > 0.0001))
  plot_mimic_variants_bars("f", ptms11, 35, 60, 10, "Common variants", "3f-snps_gnomad-common_variants", snps_gnomad %>% filter(af > 0.0001), "nodis")
  plot_mimic_variants_bars("g", ptms11, 50, 60, 10, "Rare variants (gnomAD)", "3g-snps_gnomad-rare_variants", snps_gnomad %>% filter(af <= 0.0001), "nodis")
  plot_mimic_variants_bars("h", ptms11, 35, 60, 10, "Clinically significant", "3h-snps_clinvar-clinsig", snps_clinvar %>% filter(clinsig == 1), "nodis")
  plot_mimic_variants_bars("i", ptms11, 60, 60, 10, "Cancer Gene Census (COSMIC)", "3i-snps_cosmic-cgc1", snps_cosmic %>% filter(cgc == 1), "nodis")
  
  # # Supplementary Figure 4 efgh (without showing "Disordered", to make these plots more compact) (not used)
  # plot_mimic_variants_bars("e", ptms11, 35, 60, 10, "Common variants", "s4e-snps_gnomad-common_variants", snps_gnomad %>% filter(af > 0.0001), "nodis")
  # plot_mimic_variants_bars("f", ptms11, 50, 60, 10, "Rare variants (gnomAD)", "s4f-snps_gnomad-rare_variants", snps_gnomad %>% filter(af <= 0.0001), "nodis")
  # plot_mimic_variants_bars("g", ptms11, 35, 60, 10, "Clinically significant", "s4g-snps_clinvar-clinsig", snps_clinvar %>% filter(clinsig == 1), "nodis")
  # plot_mimic_variants_bars("h", ptms11, 60, 60, 10, "Cancer Gene Census (COSMIC)", "s4h-snps_cosmic-cgc1", snps_cosmic %>% filter(cgc == 1), "nodis")
  
  # Figure 3 gh (phosphorylation only, showing gnomad-rare_variants and cosmic-cgc1)
  # plot_mimic_variants_bars("f", ptms11, 35, 60, 10, "Common variants", "3f-snps_gnomad-common_variants", snps_gnomad %>% filter(af > 0.0001), "nodis")
  # plot_mimic_variants_bars("h", ptms11, 35, 60, 10, "Clinically significant", "3h-snps_clinvar-clinsig", snps_clinvar %>% filter(clinsig == 1), "nodis")
  plot_mimic_variants_bars("g", c("S-p", "T-p", "Y-p"), 30, 36, 10, "Rare natural variants", "3g-STY-p-snps_gnomad-rare_variants", snps_gnomad %>% filter(af <= 0.0001), "nodis")
  plot_mimic_variants_bars("h", c("S-p", "T-p", "Y-p"), 45, 36, 10, "Cancer Gene Census", "3h-STY-p-snps_cosmic-cgc1", snps_cosmic %>% filter(cgc == 1), "nodis")
  
  # Figure 4 (without showing "Disordered", to make these plots more compact)
  plot_mimic_variants_bars("a", ptms11, 45, 80, 10, "Common natural variants", "4a-snps_gnomad-common_variants", snps_gnomad %>% filter(af > 0.0001), "nodis")
  plot_mimic_variants_bars("b", ptms11, 85, 80, 10, "Rare natural variants", "4b-snps_gnomad-rare_variants", snps_gnomad %>% filter(af <= 0.0001), "nodis")
  plot_mimic_variants_bars("c", ptms11, 45, 80, 10, "Clinically significant variants", "4c-snps_clinvar-clinsig", snps_clinvar %>% filter(clinsig == 1), "nodis")
  plot_mimic_variants_bars("d", ptms11, 85, 80, 10, "Cancer Gene Census", "4d-snps_cosmic-cgc1", snps_cosmic %>% filter(cgc == 1), "nodis")
  
  # Figure 4
  plot_mimic_variants_bars("a", ptms11, 45, 80, 10, "Common natural variants", "4a-snps_gnomad-common_variants", snps_gnomad %>% filter(af > 0.0001))
  plot_mimic_variants_bars("b", ptms11, 85, 80, 10, "Rare natural variants", "4b-snps_gnomad-rare_variants", snps_gnomad %>% filter(af <= 0.0001))
  plot_mimic_variants_bars("c", ptms11, 45, 80, 10, "Clinically significant variants", "4c-snps_clinvar-clinsig", snps_clinvar %>% filter(clinsig == 1))
  plot_mimic_variants_bars("d", ptms11, 85, 80, 10, "Cancer Gene Census", "4d-snps_cosmic-cgc1", snps_cosmic %>% filter(cgc == 1))
  
  # Figure E4 (Extended Figure 4 with ptms17 instead of ptms11)
  plot_mimic_variants_bars("a", ptms1000, 45, 100, 10, "Common natural variants", "E4a-snps_gnomad-common_variants", snps_gnomad %>% filter(af > 0.0001))
  plot_mimic_variants_bars("b", ptms1000, 85, 100, 10, "Rare natural variants", "E4b-snps_gnomad-rare_variants", snps_gnomad %>% filter(af <= 0.0001))
  plot_mimic_variants_bars("c", ptms1000, 45, 100, 10, "Clinically significant variants", "E4c-snps_clinvar-clinsig", snps_clinvar %>% filter(clinsig == 1))
  plot_mimic_variants_bars("d", ptms1000, 85, 100, 10, "Cancer Gene Census", "E4d-snps_cosmic-cgc1", snps_cosmic %>% filter(cgc == 1))
  
  # Figure S4 (Larger sets of variants: "all" etc.)
  plot_mimic_variants_bars("a", ptms1000, 80, 100, 10, "gnomAD (all)", "S4a-snps_gnomad-all_variants", snps_gnomad)
  plot_mimic_variants_bars("b", ptms1000, 80, 100, 10, "ClinVar (all)", "S4b-snps_clinvar-all_variants", snps_clinvar)
  plot_mimic_variants_bars("c", ptms1000, 80, 100, 10, "COSMIC (all)", "S4c-snps_cosmic-all_variants", snps_cosmic)
  plot_mimic_variants_bars("d", ptms1000, 80, 100, 10, "COSMIC CGC (all)", "S4d-snps_cosmic-cgc", snps_cosmic %>% filter(cgc != 0))
  plot_mimic_variants_bars("e", ptms1000, 60, 100, 10, "COSMIC CGC tier 2 (candidate)", "S4e-snps_cosmic-cgc2", snps_cosmic %>% filter(cgc == 2))
  plot_mimic_variants_bars("f", ptms1000, 60, 100, 10, "COSMIC ≥10 samples", "S4f-snps_cosmic-recurrent-ac-gte10", snps_cosmic %>% filter(ac >= 10))
  plot_mimic_variants_bars("g", ptms1000, 100, 100, 10, "Genomenon ≥10 articles", "S4g-snps_mastermind-papers-gte10", snps_mastermind %>% filter(ac >= 10))
  plot_mimic_variants_bars("h", ptms1000, 100, 100, 10, "Genomenon (all)", "S4h-snps_mastermind-all_variants", snps_mastermind)
}


system.time(analyse_contingency_coresurf("70_PAE2"))
# analyse_contingency_coresurf(0)
# analyse_contingency_coresurf("90_PAE1")
# analyse_contingency_coresurf("70")
# analyse_contingency_coresurf("90")
