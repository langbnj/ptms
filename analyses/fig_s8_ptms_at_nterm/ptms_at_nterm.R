blang_init()
library(viridisLite)
library(coin)
library(Hmisc)


# Are buried PTM sites more N-terminal (which would support cotranslational folding)?

# Enable pLDDT filtering:
# plddt_filtering <- 0 # Disable
plddt_filtering <- 70
# plddt_filtering <- 90

# Run pLDDT filtering?
tmp_plddt <- ""
tmp_uniseq_suffix = ""
tmp_xlab_suffix = ""
disfilt <- 0
mywidth <- 80
tag_a <- "a"
tag_b <- "b"
if (plddt_filtering != 0) {
  tmp_plddt <- f("-pLDDT{plddt_filtering}")
  tmp_uniseq_suffix = f("_pLDDT{plddt_filtering}")
  tmp_xlab_suffix = f(" (pLDDT ≥ {plddt_filtering})")
  disfilt <- 1
  mywidth <- 60
  tag_a <- "c"
  tag_b <- "d"
}

qp <- Query(f("SELECT m.ptm, m.acc, m.site, SUBSTRING(m.ptm, 1, 1) AS aa, 'ptm' AS type, SUBSTRING(s.seq, m.site, 1) AS aatest, SUBSTRING(d.seq, m.site, 1) AS disstr, SUBSTRING(cs.seq, m.site, 1) AS coresurf, LENGTH(s.seq) AS len, (m.site - 1) / (LENGTH(s.seq) - 1) AS relsite FROM unimod m, uniseq s, uniseq cs, uniseq d WHERE m.species='human' AND m.ptm IS NOT NULL AND m.acc=s.acc AND m.acc=cs.acc AND m.acc=d.acc AND s.type IN ('UniProt') AND LENGTH(s.seq)>=16 AND cs.type='CoreSurf{tmp_uniseq_suffix}' AND SUBSTRING(cs.seq, m.site, 1) IN ('C', 'S') AND d.type='AlphaFold{tmp_uniseq_suffix}' GROUP BY m.acc, m.site, m.ptm ORDER BY m.acc, m.site, m.ptm"))
qc <- Query(f("SELECT t.ptm, m.acc, m.site, m.aa, 'control' AS type, SUBSTRING(s.seq, m.site, 1) AS aatest, SUBSTRING(d.seq, m.site, 1) AS disstr, SUBSTRING(cs.seq, m.site, 1) AS coresurf, LENGTH(s.seq) AS len, (m.site - 1) / (LENGTH(s.seq) - 1) AS relsite FROM unimod_control m, uniseq s, uniseq cs, uniseq d, (SELECT m.ptm, m.acc FROM unimod m, uniseq s, uniseq cs, uniseq d WHERE m.species='human' AND m.ptm IS NOT NULL AND m.acc=s.acc AND m.acc=cs.acc AND m.acc=d.acc AND s.type IN ('UniProt') AND LENGTH(s.seq)>=16 AND cs.type='CoreSurf{tmp_uniseq_suffix}' AND SUBSTRING(cs.seq, m.site, 1) IN ('C', 'S') AND d.type='AlphaFold{tmp_uniseq_suffix}' GROUP BY m.ptm, m.acc ORDER BY m.ptm, m.acc) t WHERE m.acc=t.acc AND m.aa=SUBSTRING(t.ptm, 1, 1) AND m.species='human' AND m.acc=s.acc AND m.acc=cs.acc AND m.acc=d.acc AND s.type IN ('UniProt') AND cs.type='CoreSurf{tmp_uniseq_suffix}' AND SUBSTRING(cs.seq, m.site, 1) IN ('C', 'S') AND d.type='AlphaFold{tmp_uniseq_suffix}' GROUP BY m.acc, m.site, t.ptm, m.aa ORDER BY m.acc, m.site, t.ptm, m.aa"))
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
q %<>% mutate(type = factor(type, levels = c("ptm", "control")))
# Convert acc to factor
q %<>% mutate(acc = as.factor(acc))
# Clean up
q %<>% select(-disstr, -coresurf)
q %<>% arrange(acc, site, ptm, type)
q


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
# How did this happen? How did I get quite a few mods without controls?
# >> I think it's normal. I had to implement a filter in evorate_finding_best_stratification as well.
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


# Keep only PTMs that have ≥1000 sites
# q %>% filter(type == "ptm") %>% group_by(ptm) %>% summarise(sites = n_distinct(acc, site)) %>% arrange(desc(sites)) %>% filter(sites >= 1000)
# q %>% filter(type == "ptm") %>% group_by(ptm) %>% summarise(sites = n_distinct(acc, site)) %>% arrange(desc(sites)) %>% filter(sites >= 1000) %>% pull(ptm) %>% sort
q %>% filter(type == "ptm") %>% group_by(ptm) %>% tally %>% arrange(desc(n))
q %>% select(acc, site, ptm, type) %>% unique %>% filter(type == "ptm") %>% group_by(ptm) %>% tally %>% arrange(desc(n))
ptm1000 <- q %>% select(acc, site, ptm, type) %>% unique %>% filter(type == "ptm") %>% group_by(ptm) %>% summarise(sites = n_distinct(acc, site)) %>% arrange(desc(sites)) %>% filter(sites >= 1000) %>% pull(ptm) %>% sort
ptm1000
length(ptm1000)
q %<>% filter(ptm %in% ptm1000)

# Add "freq" column (for weighting of control residues to the same proportion as PTM sites, by AA type)
# aaweight <- q %>% filter(type == "ptm") %>% group_by(aa) %>% tally %>% mutate(freq = n / sum(n)) %>% arrange(desc(n)) %>% select(aa, freq)
# aaweight <- q %>% select(acc, site, aa, type) %>% unique %>% filter(type == "ptm") %>% group_by(aa) %>% tally %>% mutate(freq = n / sum(n)) %>% arrange(desc(n)) %>% select(aa, freq)
# aaweight
aadisweight <- q %>% select(acc, site, aa, type, dis) %>% unique %>% filter(type == "ptm") %>% group_by(aa, dis) %>% tally %>% mutate(freq = n / sum(n)) %>% arrange(desc(n)) %>% select(aa, dis, freq)
aadisweight
# Add weights to q
# q %<>% left_join(aaweight, join_by(aa))
q %<>% left_join(aadisweight, join_by(aa, dis))
# Set PTM site weight to 1 (only Control should be weighted)
q %<>% mutate(freq = ifelse(type == "ptm", 1, freq))
q




# Relative site (0-1)
# tmpptm <- q %>% filter(type == "ptm") %>% group_by(ptm) %>% summarise(mean_relsite = mean(relsite)) %>% arrange(mean_relsite) %>% pull(ptm)
# tmpptm <- q %>% filter(type == "ptm") %>% group_by(ptm) %>% summarise(median_relsite = median(relsite)) %>% arrange(median_relsite) %>% pull(ptm)
tmpptm <- q %>% filter(type == "ptm" & dis == "Buried") %>% group_by(ptm) %>% summarise(median_relsite = median(relsite)) %>% arrange(median_relsite) %>% pull(ptm)
tmpptm
(q %>%
    filter(type == "ptm") %>%
    # select(acc, site, relsite, ptm) %>% unique %>%
    # mutate(ptm = fct_rev(factor(ptm, levels = tmpptm))) %>%
    mutate(ptm = factor(ptm, levels = tmpptm)) %>%
    ggplot(aes(x = relsite, colour = ptm, fill = ptm)) +
    geom_density(alpha = 0.3) +
    scale_fill_manual(aesthetics = c("colour", "fill"), values = c(viridis_pal()(length(tmpptm))), name = NULL, guide = guide_legend(reverse = T)) +
    # scale_x_continuous(expand = c(0, 0), breaks = pretty_breaks(4)) +
    # scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), breaks = pretty_breaks(5)) +
    # scale_x_continuous(breaks = pretty_breaks(3)) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = NULL) +
    facet_grid(rows = vars(ptm), cols = vars(dis), scales = "free_y") +
    theme_minimal() +
    xlab("Relative PTM site position") +
    ylab("Probability density") +
    theme_nature() +
    guides(colour = "none", fill = "none") +
    theme(legend.text = element_text(size = 5), strip.text.y = element_text(size = 5), legend.key.size = unit(1, "mm"), axis.line.y = element_blank())
) %>% qsave(f("output-density-relsite-allptms{tmp_plddt}.pdf"), width = 80, height = 80)











# FIGURES

# Relative site (0-1) with aa-frequency-weighted control
tmpptm <- q %>% filter(type == "ptm" & dis == "Buried") %>% group_by(ptm) %>% summarise(median_relsite = median(relsite)) %>% arrange(desc(median_relsite)) %>% pull(ptm)
# tmpptm <- q %>% filter(type == "ptm" & dis == "Buried") %>% group_by(ptm) %>% summarise(mean_relsite = mean(relsite)) %>% arrange(desc(mean_relsite)) %>% pull(ptm)
# Use order as for 1-300 weighted plot below
# tmpptm <- q %>% filter(type == "ptm" & dis == "Buried" & site <= maxsite & len >= 300) %>% group_by(ptm) %>% summarise(median_site = median(site)) %>% arrange(desc(median_site)) %>% pull(ptm)
# tmpptm <- q %>% filter(type == "ptm" & dis == "Buried" & site <= maxsite & len >= 300) %>% group_by(ptm) %>% summarise(mean_site = mean(site)) %>% arrange(desc(mean_site)) %>% pull(ptm)
tmpptm
qtmp <- q
if (disfilt == 1) {
  qtmp %<>% filter(dis != "Disordered")
}
(qtmp %>%
    mutate(ptm = ifelse(type == "control", "Control", ptm)) %>%
    mutate(ptm = fct_rev(factor(ptm, levels = c("Control", tmpptm)))) %>%
    select(ptm, acc, site, relsite, dis, freq) %>% unique %>%
    ggplot(aes(x = relsite, colour = ptm, fill = ptm, weight = freq)) +
    geom_density(alpha = 0.3, bounds = c(0, 1)) +
    geom_vline(data = function(x) x %>% group_by(ptm, dis) %>% summarise(mediansite = median(relsite)), aes(xintercept = mediansite), linetype = "solid") +
    geom_vline(data = function(x) x %>% filter(ptm == "Control") %>% group_by(dis) %>% summarise(control_median = Hmisc::wtd.quantile(relsite, weights=freq, probs=0.5)), aes(xintercept = control_median), color = "grey50", linetype = "solid") +
    # geom_vline(data = function(x) x %>% group_by(ptm, dis) %>% summarise(meansite = mean(site)), aes(xintercept = meansite), linetype = "dashed") +
    # geom_vline(xintercept = q %>% filter(len >= maxsite & site <= maxsite & type == "control") %>% pull(site) %>% median, colour = "grey", linetype = "dashed") +
    # geom_vline(xintercept = q %>% filter(len >= maxsite & site <= maxsite & type == "control") %>% pull(site) %>% mean, colour = "grey", linetype = "dashed") +
    # geom_vline(xintercept = q %>% filter(len >= maxsite & site <= maxsite & type == "control") %>% pull(site, freq) %>% weighted.mean(na.rm=T), colour = "grey", linetype = "dashed") +
    scale_fill_manual(aesthetics = c("colour", "fill"), values = c(viridis_pal()(length(tmpptm) + 1)[-1], "grey"), name = NULL, guide = guide_legend(reverse = T)) +
    # scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
    # scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), limits = c(0, 1)) +
    coord_cartesian(xlim = c(0, 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = NULL) +
    facet_grid(rows = vars(ptm), cols = vars(dis), scales = "free_y") +
    theme_minimal() +
    xlab(f("Relative modification position{tmp_xlab_suffix}")) +
    ylab("Probability density") +
    labs(tag = tag_a) +
    theme_nature() +
    theme(panel.spacing.x = unit(3, "mm")) +
    guides(colour = "none", fill = "none") +
    theme(legend.text = element_text(size = 5), strip.text.y = element_text(size = 5), legend.key.size = unit(1, "mm"), axis.line.y = element_blank())
) %>% qsave(f("output-density-relsite-allptms-weighted{tmp_plddt}.pdf"), width = mywidth, height = 80)




# 
# 
# # Absolute site (up to 2500)
# # Use PTM order from relative site (0-1) weighted
# # tmpptm <- q %>% filter(type == "ptm" & dis == "Buried") %>% group_by(ptm) %>% summarise(median_relsite = median(relsite)) %>% arrange(desc(median_relsite)) %>% pull(ptm)
# tmpptm <- q %>% filter(type == "ptm" & dis == "Buried") %>% group_by(ptm) %>% summarise(mean_relsite = mean(relsite)) %>% arrange(desc(mean_relsite)) %>% pull(ptm)
# tmpptm
# quantile(q$site, probs = c(0.95,0.975,0.99))
# # >> 95th percentile is 2404, so I'll use 2500 as the maximum shown
# maxsite <- 2500
# (q %>%
#     mutate(site = ifelse(site > maxsite, maxsite, site)) %>%
#     mutate(ptm = ifelse(type == "control", "Control", ptm)) %>%
#     mutate(ptm = fct_rev(factor(ptm, levels = c("Control", tmpptm)))) %>%
#     select(ptm, acc, site, relsite, dis, freq) %>% unique %>%
#     ggplot(aes(x = site, colour = ptm, fill = ptm, weight = freq)) +
#     geom_density(alpha = 0.3, bounds = c(1, maxsite)) +
#     scale_fill_manual(aesthetics = c("colour", "fill"), values = c(viridis_pal()(length(tmpptm) + 1)[-1], "grey"), name = NULL, guide = guide_legend(reverse = T)) +
#     # scale_x_continuous(expand = c(0, 0), breaks = pretty_breaks(5)) +
#     # scale_x_log10(expand = c(0, 0), breaks = pretty_breaks(5)) +
#     scale_x_continuous(expand = c(0, 0), breaks = c(1, 500, 1000, 1500, 2000, 2500), labels = c("1", "500", "1000", "1500", "2000", "≥ 2500")) +
#     scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = NULL) +
#     coord_cartesian(xlim = c(1, maxsite)) +
#     facet_grid(rows = vars(ptm), cols = vars(dis), scales = "free") +
#     xlab("PTM site position") +
#     ylab("Probability density") +
#     # labs(tag = "a") +
#     theme_nature() +
#     theme(panel.spacing.x = unit(3, "mm")) +
#     guides(colour = "none", fill = "none") +
#     theme(legend.text = element_text(size = 5), strip.text.y = element_text(size = 5), legend.key.size = unit(1, "mm"), axis.line.y = element_blank())
# ) %>% qsave(f("output-density-site2500-allptms-weighted{tmp_plddt}.pdf"), width = 80, height = 80)
# 
# 
# 
# # Absolute site (up to 2500)
# # Use PTM order from relative site (0-1) weighted
# # tmpptm <- q %>% filter(type == "ptm" & dis == "Buried") %>% group_by(ptm) %>% summarise(median_relsite = median(relsite)) %>% arrange(desc(median_relsite)) %>% pull(ptm)
# tmpptm <- q %>% filter(type == "ptm" & dis == "Buried") %>% group_by(ptm) %>% summarise(mean_relsite = mean(relsite)) %>% arrange(desc(mean_relsite)) %>% pull(ptm)
# tmpptm
# quantile(q$site, probs = c(0.5))
# median(q$site)
# # >> median is 395, which will give us half of the proteins if we use ≥395 (good)
# maxsite <- median(q$site)
# (q %>%
#     mutate(site = ifelse(site > maxsite, maxsite, site)) %>%
#     mutate(ptm = ifelse(type == "control", "Control", ptm)) %>%
#     mutate(ptm = fct_rev(factor(ptm, levels = c("Control", tmpptm)))) %>%
#     select(ptm, acc, site, relsite, dis, freq) %>% unique %>%
#     ggplot(aes(x = site, colour = ptm, fill = ptm, weight = freq)) +
#     geom_density(alpha = 0.3, bounds = c(1, maxsite)) +
#     scale_fill_manual(aesthetics = c("colour", "fill"), values = c(viridis_pal()(length(tmpptm) + 1)[-1], "grey"), name = NULL, guide = guide_legend(reverse = T)) +
#     # scale_x_continuous(expand = c(0, 0), breaks = pretty_breaks(5)) +
#     # scale_x_log10(expand = c(0, 0), breaks = pretty_breaks(5)) +
#     scale_x_continuous(expand = c(0, 0), breaks = c(1, 500, 1000, 1500, 2000, 2500), labels = c("1", "500", "1000", "1500", "2000", "≥ 2500")) +
#     scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = NULL) +
#     coord_cartesian(xlim = c(1, maxsite)) +
#     facet_grid(rows = vars(ptm), cols = vars(dis), scales = "free") +
#     xlab("PTM site position") +
#     ylab("Probability density") +
#     # labs(tag = "a") +
#     theme_nature() +
#     theme(panel.spacing.x = unit(3, "mm")) +
#     guides(colour = "none", fill = "none") +
#     theme(legend.text = element_text(size = 5), strip.text.y = element_text(size = 5), legend.key.size = unit(1, "mm"), axis.line.y = element_blank())
# ) %>% qsave(f("output-density-site-allptms-weighted{tmp_plddt}.pdf"), width = 80, height = 80)



# # Absolute site (1-300)
# # tmpptm <- q %>% filter(type == "ptm") %>% group_by(ptm) %>% summarise(mean_relsite = mean(relsite)) %>% arrange(mean_relsite) %>% pull(ptm)
# # tmpptm <- q %>% filter(type == "ptm") %>% group_by(ptm) %>% summarise(median_relsite = median(relsite)) %>% arrange(median_relsite) %>% pull(ptm)
# maxsite <- 300
# tmpptm <- q %>% filter(type == "ptm" & dis == "Buried" & site <= maxsite) %>% group_by(ptm) %>% summarise(mean_site = mean(site)) %>% arrange(mean_site) %>% pull(ptm)
# # tmpptm <- q %>% filter(type == "ptm" & dis == "Buried" & site <= maxsite) %>% group_by(ptm) %>% summarise(median_site = median(site)) %>% arrange(median_site) %>% pull(ptm)
# tmpptm
# (q %>%
#     filter(type == "ptm") %>%
#     # filter(site <= maxsite) %>%
#     # select(acc, site, relsite, ptm) %>% unique %>%
#     # mutate(ptm = fct_rev(factor(ptm, levels = tmpptm))) %>%
#     mutate(ptm = factor(ptm, levels = tmpptm)) %>%
#     ggplot(aes(x = site, colour = ptm, fill = ptm)) +
#     geom_density(alpha = 0.3) +
#     scale_fill_manual(aesthetics = c("colour", "fill"), values = c(viridis_pal()(length(tmpptm))), name = NULL, guide = guide_legend(reverse = T)) +
#     # scale_x_continuous(expand = c(0, 0), breaks = pretty_breaks(4)) +
#     # scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), breaks = pretty_breaks(5)) +
#     # scale_x_continuous(breaks = pretty_breaks(3)) +
#     # scale_x_continuous(expand = c(0, 0), breaks = c(1, 150, 300)) +
#     # scale_x_continuous(breaks = c(1, 150, 300)) +
#     scale_x_continuous(breaks = pretty_breaks(3), limits = c(1, maxsite)) +
#     scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = NULL) +
#     # coord_cartesian(xlim = c(1, maxsite)) +
#     facet_grid(rows = vars(ptm), cols = vars(dis), scales = "free_y") +
#     theme_minimal() +
#     xlab("Relative PTM site position") +
#     ylab("Probability density") +
#     theme_nature() +
#     guides(colour = "none", fill = "none") +
#     theme(legend.text = element_text(size = 5), strip.text.y = element_text(size = 5), legend.key.size = unit(1, "mm"), axis.line.y = element_blank())
# ) %>% qsave(f("output-density-site{maxsite}-allptms{tmp_plddt}.pdf"), width = 80, height = 80)




# Absolute site (1-300) with aa-frequency-weighted control
# tmpptm <- q %>% filter(type == "ptm") %>% group_by(ptm) %>% summarise(mean_relsite = mean(relsite)) %>% arrange(mean_relsite) %>% pull(ptm)
# tmpptm <- q %>% filter(type == "ptm") %>% group_by(ptm) %>% summarise(median_relsite = median(relsite)) %>% arrange(median_relsite) %>% pull(ptm)
maxsite <- 300
# tmpptm <- q %>% filter(type == "ptm" & dis == "Buried" & site <= maxsite) %>% group_by(ptm) %>% summarise(median_site = median(site)) %>% arrange(median_site) %>% pull(ptm)
# tmpptm <- q %>% filter(type == "ptm" & dis == "Buried" & site <= maxsite) %>% group_by(ptm) %>% summarise(mean_site = mean(site)) %>% arrange(desc(mean_site)) %>% pull(ptm)
# tmpptm <- q %>% filter(type == "ptm" & dis == "Buried" & site <= maxsite & len >= 300) %>% group_by(ptm) %>% summarise(median_site = median(site)) %>% arrange(desc(median_site)) %>% pull(ptm)
# tmpptm <- q %>% filter(type == "ptm" & dis == "Buried" & site <= maxsite & len >= 300) %>% group_by(ptm) %>% summarise(mean_site = mean(site)) %>% arrange(desc(mean_site)) %>% pull(ptm)
# Use PTM order from relative site (0-1) weighted
tmpptm <- q %>% filter(type == "ptm" & dis == "Buried") %>% group_by(ptm) %>% summarise(median_relsite = median(relsite)) %>% arrange(desc(median_relsite)) %>% pull(ptm)
# tmpptm <- q %>% filter(type == "ptm" & dis == "Buried") %>% group_by(ptm) %>% summarise(mean_relsite = mean(relsite)) %>% arrange(desc(mean_relsite)) %>% pull(ptm)
tmpptm
qtmp <- q
if (disfilt == 1) {
  qtmp %<>% filter(dis != "Disordered")
}
(qtmp %>%
    # Protein length ≥ 300
    filter(len >= maxsite) %>%
    # filter(type == "ptm") %>%
    mutate(ptm = ifelse(type == "control", "Control", ptm)) %>%
    mutate(ptm = fct_rev(factor(ptm, levels = c("Control", tmpptm)))) %>%
    select(ptm, acc, site, relsite, dis, freq) %>% unique %>%
    filter(site <= maxsite) %>%
    # select(acc, site, relsite, ptm) %>% unique %>%
    # mutate(ptm = fct_rev(factor(ptm, levels = tmpptm))) %>%
    # mutate(ptm = factor(ptm, levels = tmpptm)) %>%
    ggplot(aes(x = site, colour = ptm, fill = ptm, weight = freq)) +
    geom_density(alpha = 0.3, bounds = c(1, maxsite)) +
    geom_vline(data = function(x) x %>% group_by(ptm, dis) %>% summarise(mediansite = median(site)), aes(xintercept = mediansite), linetype = "solid") +
    geom_vline(data = function(x) x %>% filter(ptm == "Control") %>% group_by(dis) %>% summarise(control_median = Hmisc::wtd.quantile(site, weights=freq, probs=0.5)), aes(xintercept = control_median), color = "grey50", linetype = "solid") +
    # geom_vline(data = function(x) x %>% group_by(ptm, dis) %>% summarise(meansite = mean(site)), aes(xintercept = meansite), linetype = "dashed") +
    # geom_vline(xintercept = q %>% filter(len >= maxsite & site <= maxsite & type == "control") %>% pull(site) %>% median, colour = "grey", linetype = "dashed") +
    # geom_vline(xintercept = q %>% filter(len >= maxsite & site <= maxsite & type == "control") %>% pull(site) %>% mean, colour = "grey", linetype = "dashed") +
    # geom_vline(xintercept = q %>% filter(len >= maxsite & site <= maxsite & type == "control") %>% pull(site, freq) %>% weighted.mean(na.rm=T), colour = "grey", linetype = "dashed") +
    scale_fill_manual(aesthetics = c("colour", "fill"), values = c(viridis_pal()(length(tmpptm) + 1)[-1], "grey"), name = NULL, guide = guide_legend(reverse = T)) +
    # scale_fill_manual(aesthetics = c("colour", "fill"), values = c(viridis_pal()(length(tmpptm))), name = NULL, guide = guide_legend(reverse = T)) +
    # scale_x_continuous(expand = c(0, 0), breaks = pretty_breaks(4)) +
    # scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), breaks = pretty_breaks(5)) +
    # scale_x_continuous(breaks = pretty_breaks(3)) +
    # scale_x_continuous(expand = c(0, 0), breaks = c(1, 150, 300)) +
    # scale_x_continuous(breaks = c(1, 150, 300)) +
    # For the kernel, add 10% of non-plotted data
    # scale_x_continuous(breaks = pretty_breaks(3), limits = c(1, maxsite * 1.1)) +
    scale_x_continuous(expand = c(0, 0), breaks = c(1, 100, 200, 300)) +
    # maxsite * x: not necessary since I can define bounds in geom_density.
    # scale_x_continuous(expand = c(0, 0), breaks = c(1, 100, 200, 300), limits = c(1, maxsite * 1.2)) +
    # scale_x_continuous(expand = c(0, 0), breaks = c(1, 100, 200, 300), limits = c(1, maxsite * 2)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = NULL) +
    # scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = pretty_breaks(2)) +
    coord_cartesian(xlim = c(1, maxsite)) +
    # facet_grid(rows = vars(ptm), cols = vars(dis), scales = "free_y") +
    facet_grid(rows = vars(ptm), cols = vars(dis), scales = "free") +
    # facet_wrap(facets = vars(ptm, dis), ncol = 3, scales = "free") +
    # theme_minimal() +
    xlab(f("Mod. position{tmp_xlab_suffix} (first 300 aa, prot. ≥ 300 aa)")) +
    ylab("Probability density") +
    labs(tag = tag_b) +
    theme_nature() +
    theme(panel.spacing.x = unit(3, "mm")) +
    guides(colour = "none", fill = "none") +
    theme(legend.text = element_text(size = 5), strip.text.y = element_text(size = 5), legend.key.size = unit(1, "mm"), axis.line.y = element_blank())
) %>% qsave(f("output-density-site{maxsite}-allptms-weighted{tmp_plddt}.pdf"), width = mywidth, height = 80)



# # Absolute site (1-395, median) with aa-frequency-weighted control
# maxsite <- median(q$site)
# # >> 395
# # Use PTM order from relative site (0-1) weighted
# # tmpptm <- q %>% filter(type == "ptm" & dis == "Buried") %>% group_by(ptm) %>% summarise(median_relsite = median(relsite)) %>% arrange(desc(median_relsite)) %>% pull(ptm)
# tmpptm <- q %>% filter(type == "ptm" & dis == "Buried") %>% group_by(ptm) %>% summarise(mean_relsite = mean(relsite)) %>% arrange(desc(mean_relsite)) %>% pull(ptm)
# tmpptm
# (q %>%
#     # Protein length ≥ 395
#     filter(len >= maxsite) %>%
#     # filter(type == "ptm") %>%
#     mutate(ptm = ifelse(type == "control", "Control", ptm)) %>%
#     mutate(ptm = fct_rev(factor(ptm, levels = c("Control", tmpptm)))) %>%
#     select(ptm, acc, site, relsite, dis, freq) %>% unique %>%
#     filter(site <= maxsite) %>%
#     # select(acc, site, relsite, ptm) %>% unique %>%
#     # mutate(ptm = fct_rev(factor(ptm, levels = tmpptm))) %>%
#     # mutate(ptm = factor(ptm, levels = tmpptm)) %>%
#     ggplot(aes(x = site, colour = ptm, fill = ptm, weight = freq)) +
#     geom_density(alpha = 0.3, bounds = c(1, maxsite)) +
#     scale_fill_manual(aesthetics = c("colour", "fill"), values = c(viridis_pal()(length(tmpptm) + 1)[-1], "grey"), name = NULL, guide = guide_legend(reverse = T)) +
#     # scale_fill_manual(aesthetics = c("colour", "fill"), values = c(viridis_pal()(length(tmpptm))), name = NULL, guide = guide_legend(reverse = T)) +
#     # scale_x_continuous(expand = c(0, 0), breaks = pretty_breaks(4)) +
#     # scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), breaks = pretty_breaks(5)) +
#     # scale_x_continuous(breaks = pretty_breaks(3)) +
#     # scale_x_continuous(expand = c(0, 0), breaks = c(1, 150, 300)) +
#     # scale_x_continuous(breaks = c(1, 150, 300)) +
#     # For the kernel, add 10% of non-plotted data
#     # scale_x_continuous(breaks = pretty_breaks(3), limits = c(1, maxsite * 1.1)) +
#     scale_x_continuous(expand = c(0, 0), breaks = c(1, 100, 200, 300, maxsite)) +
#     # maxsite * x: not necessary since I can define bounds in geom_density.
#     # scale_x_continuous(expand = c(0, 0), breaks = c(1, 100, 200, 300), limits = c(1, maxsite * 1.2)) +
#     # scale_x_continuous(expand = c(0, 0), breaks = c(1, 100, 200, 300), limits = c(1, maxsite * 2)) +
#     scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = NULL) +
#     # scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = pretty_breaks(2)) +
#     coord_cartesian(xlim = c(1, maxsite)) +
#     # facet_grid(rows = vars(ptm), cols = vars(dis), scales = "free_y") +
#     facet_grid(rows = vars(ptm), cols = vars(dis), scales = "free") +
#     # facet_wrap(facets = vars(ptm, dis), ncol = 3, scales = "free") +
#     # theme_minimal() +
#     xlab(f("PTM site position (first {maxsite} residues, in proteins ≥ median length, i.e. 395 aa)")) +
#     ylab("Probability density") +
#     labs(tag = "b") +
#     theme_nature() +
#     theme(panel.spacing.x = unit(3, "mm")) +
#     guides(colour = "none", fill = "none") +
#     theme(legend.text = element_text(size = 5), strip.text.y = element_text(size = 5), legend.key.size = unit(1, "mm"), axis.line.y = element_blank())
# ) %>% qsave(f("output-density-site{maxsite}-allptms-weighted{tmp_plddt}.pdf"), width = 80, height = 80)








# RELATIVE (RELSITE) WITH UNIQUING FOR ACC|SITES, PLDDT FILTERING ≥ 70 (same code as below, just different comments)
# Adding these everywhere:
#  %>% select(acc, site, relsite, type) %>% unique
#  %>% select(ptmacc, site, relsite, type) %>% unique
# Compare relsite for PTM vs Control
minsize <- 2

# Disordered
wilcox.test(relsite ~ type, conf.int = T, q %>% filter(dis == "Disordered") %>% select(acc, site, relsite, type) %>% unique)
wilcox.test(relsite ~ type, conf.int = T, q %>% filter(dis == "Disordered") %>% select(acc, site, relsite, type) %>% unique)$p.value
wilcox_test(relsite ~ type | acc, q %>% filter(dis == "Disordered") %>% select(acc, site, relsite, type) %>% unique %>% filter(!(acc %in% (q %>% filter(dis == "Disordered") %>% select(acc, site, relsite, type) %>% unique %>% group_by(acc, type) %>% tally %>% filter(n < minsize) %>% pull(acc)))))
wilcox_test(relsite ~ type | ptmacc, q %>% mutate(ptmacc = as.factor(paste(ptm, acc))) %>% filter(dis == "Disordered") %>% select(ptmacc, site, relsite, type) %>% unique %>% filter(!(ptmacc %in% (q %>% mutate(ptmacc = as.factor(paste(ptm, acc))) %>% filter(dis == "Disordered") %>% select(ptmacc, site, relsite, type) %>% unique %>% group_by(ptmacc, type) %>% tally %>% filter(n < minsize) %>% pull(ptmacc)))))
# Two-tailed:
# >> W = 154586776, p-value = 0.1328
# >> loc -0.004016978 (-0.4%)
# Asymptotic Wilcoxon-Mann-Whitney Test
# >> acc stratified: Z = -2.3858, p-value = 0.01704
# >> ptm|acc stratified: -1.7199, p-value = 0.08545
q %>% filter(type == "ptm" & dis == "Disordered") %>% select(acc, site, relsite, type) %>% unique %>% nrow
q %>% filter(type == "control" & dis == "Disordered") %>% select(acc, site, relsite, type) %>% unique %>% nrow
# >> n(ptm)     >> 13,273
# >> n(control) >> 23,515

# Structured
wilcox.test(relsite ~ type, conf.int = T, q %>% filter(dis == "Structured") %>% select(acc, site, relsite, type) %>% unique)
wilcox.test(relsite ~ type, conf.int = T, q %>% filter(dis == "Structured") %>% select(acc, site, relsite, type) %>% unique)$p.value
wilcox_test(relsite ~ type | acc, q %>% filter(dis == "Structured") %>% select(acc, site, relsite, type) %>% unique %>% filter(!(acc %in% (q %>% filter(dis == "Structured") %>% select(acc, site, relsite, type) %>% unique %>% group_by(acc, type) %>% tally %>% filter(n < minsize) %>% pull(acc)))))
wilcox_test(relsite ~ type | ptmacc, q %>% mutate(ptmacc = as.factor(paste(ptm, acc))) %>% filter(dis == "Structured") %>% select(ptmacc, site, relsite, type) %>% unique %>% filter(!(ptmacc %in% (q %>% mutate(ptmacc = as.factor(paste(ptm, acc))) %>% filter(dis == "Structured") %>% select(ptmacc, site, relsite, type) %>% unique %>% group_by(ptmacc, type) %>% tally %>% filter(n < minsize) %>% pull(ptmacc)))))
# Two-tailed:
# >> W = 5.1876e+10, p-value = 0.3825
# >> loc 0.0006406653 (+0.06%)
# >> acc stratified: Z = -3.0662, p-value = 0.002168
# >> ptm|acc stratified: Z = -2.9726, p-value = 0.002953
q %>% filter(type == "ptm" & dis == "Structured") %>% select(acc, site, relsite, type) %>% unique %>% nrow
q %>% filter(type == "control" & dis == "Structured") %>% select(acc, site, relsite, type) %>% unique %>% nrow
# >> n(ptm)     >> 178,910
# >> n(control) >> 579,120

# Buried
wilcox.test(relsite ~ type, conf.int = T, q %>% filter(dis == "Buried") %>% select(acc, site, relsite, type) %>% unique)
wilcox.test(relsite ~ type, conf.int = T, q %>% filter(dis == "Buried") %>% select(acc, site, relsite, type) %>% unique)$p.value
wilcox_test(relsite ~ type | acc, q %>% filter(dis == "Buried") %>% select(acc, site, relsite, type) %>% unique %>% filter(!(acc %in% (q %>% filter(dis == "Buried") %>% select(acc, site, relsite, type) %>% unique %>% group_by(acc, type) %>% tally %>% filter(n < minsize) %>% pull(acc)))))
wilcox_test(relsite ~ type | ptmacc, q %>% mutate(ptmacc = as.factor(paste(ptm, acc))) %>% filter(dis == "Buried") %>% select(ptmacc, site, relsite, type) %>% unique %>% filter(!(ptmacc %in% (q %>% mutate(ptmacc = as.factor(paste(ptm, acc))) %>% filter(dis == "Buried") %>% select(ptmacc, site, relsite, type) %>% unique %>% group_by(ptmacc, type) %>% tally %>% filter(n < minsize) %>% pull(ptmacc)))))
# Two-tailed:
# W = 2.2529e+10, p-value 2.293402e-21
# loc  -0.009163429 (-0.9%)
# >> acc stratified: Z = -8.8586, p-value < 2.2e-16
# >> ptm|acc stratified: Z = -5.8452, p-value = 5.06e-09
q %>% filter(type == "ptm" & dis == "Buried") %>% select(acc, site, relsite, type) %>% unique %>% nrow
q %>% filter(type == "control" & dis == "Buried") %>% select(acc, site, relsite, type) %>% unique %>% nrow
# >> n(ptm)     >> 105,183
# >> n(control) >> 436,596





# ABSOLUTE (SITE) 1-300 CLEANED UP WITH UNIQUING FOR ACC|SITES, PLDDT FILTERING ≥ 70 (same code as below, just different comments)
# Adding these everywhere:
#  %>% select(acc, site, type) %>% unique
#  %>% select(ptmacc, site, type) %>% unique
# Compare relsite for PTM vs Control
minsize <- 2
maxsite <- 300
# In absolute terms (truncated at 300 and requiring protein to be ≥300):
q %>% filter(type == "ptm" & site <= maxsite & len >= maxsite) %>% select(ptm, acc, site) %>% unique %>% nrow
q %>% filter(type == "ptm") %>% select(ptm, acc, site) %>% unique %>% nrow
q %>% filter(type == "ptm" & site <= maxsite & len >= maxsite) %>% select(ptm, acc, site) %>% unique %>% nrow / q %>% filter(type == "ptm") %>% select(ptm, acc, site) %>% unique %>% nrow
# >> 37.7% of ptm|acc|sites
q %>% filter(type == "ptm" & site <= maxsite & len >= maxsite) %>% select(acc, site) %>% unique %>% nrow / q %>% filter(type == "ptm") %>% select(acc, site) %>% unique %>% nrow
# >> 37.7% (35%) of acc|sites
# >> EXPORT
q %>% filter(type == "ptm" & site <= maxsite & len >= maxsite) %>% select(acc) %>% unique %>% nrow / q %>% filter(type == "ptm") %>% select(acc) %>% unique %>% nrow
# >> 67.4% (67%) of proteins
# >> EXPORT

# Disordered
wilcox.test(site ~ type, conf.int = T, q %>% filter(dis == "Disordered" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique)
wilcox.test(site ~ type, conf.int = T, q %>% filter(dis == "Disordered" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique)$p.value
wilcox_test(site ~ type, conf.int = T, q %>% filter(dis == "Disordered" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique)
wilcox_test(site ~ type | acc, q %>% filter(dis == "Disordered" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique %>% filter(!(acc %in% (q %>% filter(dis == "Disordered" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique %>% group_by(acc, type) %>% tally %>% filter(n < minsize) %>% pull(acc)))))
wilcox_test(site ~ type | ptmacc, q %>% mutate(ptmacc = as.factor(paste(ptm, acc))) %>% filter(dis == "Disordered" & site <= maxsite & len >= maxsite) %>% select(ptmacc, site, type) %>% unique %>% filter(!(ptmacc %in% (q %>% mutate(ptmacc = as.factor(paste(ptm, acc))) %>% filter(dis == "Disordered" & site <= maxsite & len >= maxsite) %>% select(ptmacc, site, type) %>% unique %>% group_by(ptmacc, type) %>% tally %>% filter(n < minsize) %>% pull(ptmacc)))))
# Two-tailed:
# >> W = 18515698, p-value = 0.0003512
# >> loc 5.000028
# >> acc stratified: Z = 3.5743, p-value = 0.0003512
# >> ptm|acc stratified: Z = 0.27936, p-value = 0.78
q %>% filter(type == "ptm" & dis == "Disordered" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique %>% nrow
q %>% filter(type == "control" & dis == "Disordered" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique %>% nrow
# >> n(ptm)     >> 4,281
# >> n(control) >> 8,327

# Structured
wilcox.test(site ~ type, conf.int = T, q %>% filter(dis == "Structured" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique)
wilcox.test(site ~ type, conf.int = T, q %>% filter(dis == "Structured" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique)$p.value
wilcox_test(site ~ type | acc, q %>% filter(dis == "Structured" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique %>% filter(!(acc %in% (q %>% filter(dis == "Structured" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique %>% group_by(acc, type) %>% tally %>% filter(n < minsize) %>% pull(acc)))))
wilcox_test(site ~ type | ptmacc, q %>% mutate(ptmacc = as.factor(paste(ptm, acc))) %>% filter(dis == "Structured" & site <= maxsite & len >= maxsite) %>% select(ptmacc, site, type) %>% unique %>% filter(!(ptmacc %in% (q %>% mutate(ptmacc = as.factor(paste(ptm, acc))) %>% filter(dis == "Structured" & site <= maxsite & len >= maxsite) %>% select(ptmacc, site, type) %>% unique %>% group_by(ptmacc, type) %>% tally %>% filter(n < minsize) %>% pull(ptmacc)))))
# W = 6596263675, p-value 7.248868e-17
# loc -3.000033
# >> acc stratified: Z = -5.7799, p-value = 7.474e-09
# >> ptm|acc stratified: Z = -6.4449, p-value = 1.157e-10
q %>% filter(type == "ptm" & dis == "Structured" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique %>% nrow
q %>% filter(type == "control" & dis == "Structured" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique %>% nrow
# >> n(ptm)     >> 66,349
# >> n(control) >> 203,212

# Buried
wilcox.test(site ~ type, conf.int = T, q %>% filter(dis == "Buried" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique)
wilcox.test(site ~ type, conf.int = T, q %>% filter(dis == "Buried" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique)$p.value
wilcox_test(site ~ type, conf.int = T, q %>% filter(dis == "Buried" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique)
wilcox_test(site ~ type | acc, q %>% filter(dis == "Buried" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique %>% filter(!(acc %in% (q %>% filter(dis == "Buried" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique %>% group_by(acc, type) %>% tally %>% filter(n < minsize) %>% pull(acc)))))
wilcox_test(site ~ type | ptmacc, q %>% mutate(ptmacc = as.factor(paste(ptm, acc))) %>% filter(dis == "Buried" & site <= maxsite & len >= maxsite) %>% select(ptmacc, site, type) %>% unique %>% filter(!(ptmacc %in% (q %>% mutate(ptmacc = as.factor(paste(ptm, acc))) %>% filter(dis == "Buried" & site <= maxsite & len >= maxsite) %>% select(ptmacc, site, type) %>% unique %>% group_by(ptmacc, type) %>% tally %>% filter(n < minsize) %>% pull(ptmacc)))))
# W = 3.12e+09, p-value 1.39541e-58
# loc -7.000044
# >> acc stratified: Z = -12.266, p-value < 2.2e-16
# >> ptm|acc stratified: Z = -9.37, p-value < 2.2e-16
q %>% filter(type == "ptm" & dis == "Buried" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique %>% nrow
q %>% filter(type == "control" & dis == "Buried" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique %>% nrow
# >> n(ptm)     >> 41,582
# >> n(control) >> 158,189






# NO PLDDT FILTERING BELOW

# RELATIVE (RELSITE) WITH UNIQUING FOR ACC|SITES, NO PLDDT FILTERING
# Adding these everywhere:
#  %>% select(acc, site, relsite, type) %>% unique
#  %>% select(ptmacc, site, relsite, type) %>% unique
# Compare relsite for PTM vs Control
minsize <- 2

# Disordered
wilcox.test(relsite ~ type, conf.int = T, q %>% filter(dis == "Disordered") %>% select(acc, site, relsite, type) %>% unique)
wilcox.test(relsite ~ type, conf.int = T, q %>% filter(dis == "Disordered") %>% select(acc, site, relsite, type) %>% unique)$p.value
wilcox_test(relsite ~ type | acc, q %>% filter(dis == "Disordered") %>% select(acc, site, relsite, type) %>% unique %>% filter(!(acc %in% (q %>% filter(dis == "Disordered") %>% select(acc, site, relsite, type) %>% unique %>% group_by(acc, type) %>% tally %>% filter(n < minsize) %>% pull(acc)))))
wilcox_test(relsite ~ type | ptmacc, q %>% mutate(ptmacc = as.factor(paste(ptm, acc))) %>% filter(dis == "Disordered") %>% select(ptmacc, site, relsite, type) %>% unique %>% filter(!(ptmacc %in% (q %>% mutate(ptmacc = as.factor(paste(ptm, acc))) %>% filter(dis == "Disordered") %>% select(ptmacc, site, relsite, type) %>% unique %>% group_by(ptmacc, type) %>% tally %>% filter(n < minsize) %>% pull(ptmacc)))))
# Two-tailed:
# >> W = 9.607e+10, p-value 2.236629e-20
# >> loc 0.005707067 (+0.57%)
# >> acc stratified: Z = 2.2016, p-value = 0.02769
# >> ptm|acc stratified: Z = -2.0055, p-value = 0.04491
q %>% filter(type == "ptm" & dis == "Disordered") %>% select(acc, site, relsite, type) %>% unique %>% nrow
q %>% filter(type == "control" & dis == "Disordered") %>% select(acc, site, relsite, type) %>% unique %>% nrow
# >> n(ptm)     >> 311,545
# >> n(control) >> 609,563

# Structured
wilcox.test(relsite ~ type, conf.int = T, q %>% filter(dis == "Structured") %>% select(acc, site, relsite, type) %>% unique)
wilcox.test(relsite ~ type, conf.int = T, q %>% filter(dis == "Structured") %>% select(acc, site, relsite, type) %>% unique)$p.value
wilcox_test(relsite ~ type | acc, q %>% filter(dis == "Structured") %>% select(acc, site, relsite, type) %>% unique %>% filter(!(acc %in% (q %>% filter(dis == "Structured") %>% select(acc, site, relsite, type) %>% unique %>% group_by(acc, type) %>% tally %>% filter(n < minsize) %>% pull(acc)))))
wilcox_test(relsite ~ type | ptmacc, q %>% mutate(ptmacc = as.factor(paste(ptm, acc))) %>% filter(dis == "Structured") %>% select(ptmacc, site, relsite, type) %>% unique %>% filter(!(ptmacc %in% (q %>% mutate(ptmacc = as.factor(paste(ptm, acc))) %>% filter(dis == "Structured") %>% select(ptmacc, site, relsite, type) %>% unique %>% group_by(ptmacc, type) %>% tally %>% filter(n < minsize) %>% pull(ptmacc)))))
# Two-tailed:
# >> W = 8.3056e+10, p-value = 0.002595
# >> loc 0.002126313 (+0.21%)
# >> acc stratified: Z = -1.2849, p-value = 0.1988
# >> ptm|acc stratified: Z = -1.5851, p-value = 0.1129
q %>% filter(type == "ptm" & dis == "Structured") %>% select(acc, site, relsite, type) %>% unique %>% nrow
q %>% filter(type == "control" & dis == "Structured") %>% select(acc, site, relsite, type) %>% unique %>% nrow
# >> n(ptm)     >> 215,106
# >> n(control) >> 768,975

# Buried
wilcox.test(relsite ~ type, conf.int = T, q %>% filter(dis == "Buried") %>% select(acc, site, relsite, type) %>% unique)
wilcox.test(relsite ~ type, conf.int = T, q %>% filter(dis == "Buried") %>% select(acc, site, relsite, type) %>% unique)$p.value
wilcox_test(relsite ~ type | acc, q %>% filter(dis == "Buried") %>% select(acc, site, relsite, type) %>% unique %>% filter(!(acc %in% (q %>% filter(dis == "Buried") %>% select(acc, site, relsite, type) %>% unique %>% group_by(acc, type) %>% tally %>% filter(n < minsize) %>% pull(acc)))))
wilcox_test(relsite ~ type | ptmacc, q %>% mutate(ptmacc = as.factor(paste(ptm, acc))) %>% filter(dis == "Buried") %>% select(ptmacc, site, relsite, type) %>% unique %>% filter(!(ptmacc %in% (q %>% mutate(ptmacc = as.factor(paste(ptm, acc))) %>% filter(dis == "Buried") %>% select(ptmacc, site, relsite, type) %>% unique %>% group_by(ptmacc, type) %>% tally %>% filter(n < minsize) %>% pull(ptmacc)))))
# Two-tailed:
# >> W = 2.7548e+10, p-value 4.421865e-17
# >> loc -0.007808839 (-0.8%)
# >> acc stratified: Z = -8.0336, p-value = 9.462e-16
# >> ptm|acc stratified: Z = -5.0019, p-value = 5.676e-07
q %>% filter(type == "ptm" & dis == "Buried") %>% select(acc, site, relsite, type) %>% unique %>% nrow
q %>% filter(type == "control" & dis == "Buried") %>% select(acc, site, relsite, type) %>% unique %>% nrow
# >> n(ptm)     >> 112,451
# >> n(control) >> 497,924





# ABSOLUTE (SITE) 1-300 CLEANED UP WITH UNIQUING FOR ACC|SITES, NO PLDDT FILTERING
# Adding these everywhere:
#  %>% select(acc, site, type) %>% unique
#  %>% select(ptmacc, site, type) %>% unique
# Compare relsite for PTM vs Control
minsize <- 2
maxsite <- 300
# In absolute terms (truncated at 300 and requiring protein to be ≥300):
q %>% filter(type == "ptm" & site <= maxsite & len >= maxsite) %>% select(ptm, acc, site) %>% unique %>% nrow
q %>% filter(type == "ptm") %>% select(ptm, acc, site) %>% unique %>% nrow
q %>% filter(type == "ptm" & site <= maxsite & len >= maxsite) %>% select(ptm, acc, site) %>% unique %>% nrow / q %>% filter(type == "ptm") %>% select(ptm, acc, site) %>% unique %>% nrow
# >> 34.6% of ptm|acc|sites
q %>% filter(type == "ptm" & site <= maxsite & len >= maxsite) %>% select(acc, site) %>% unique %>% nrow / q %>% filter(type == "ptm") %>% select(acc, site) %>% unique %>% nrow
# >> 34.5% (35%) of acc|sites
# >> EXPORT
q %>% filter(type == "ptm" & site <= maxsite & len >= maxsite) %>% select(acc) %>% unique %>% nrow / q %>% filter(type == "ptm") %>% select(acc) %>% unique %>% nrow
# >> 69.8% (70%) of proteins
# >> EXPORT

# Disordered
wilcox.test(site ~ type, conf.int = T, q %>% filter(dis == "Disordered" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique)
wilcox.test(site ~ type, conf.int = T, q %>% filter(dis == "Disordered" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique)$p.value
wilcox_test(site ~ type, conf.int = T, q %>% filter(dis == "Disordered" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique)
wilcox_test(site ~ type | acc, q %>% filter(dis == "Disordered" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique %>% filter(!(acc %in% (q %>% filter(dis == "Disordered" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique %>% group_by(acc, type) %>% tally %>% filter(n < minsize) %>% pull(acc)))))
wilcox_test(site ~ type | ptmacc, q %>% mutate(ptmacc = as.factor(paste(ptm, acc))) %>% filter(dis == "Disordered" & site <= maxsite & len >= maxsite) %>% select(ptmacc, site, type) %>% unique %>% filter(!(ptmacc %in% (q %>% mutate(ptmacc = as.factor(paste(ptm, acc))) %>% filter(dis == "Disordered" & site <= maxsite & len >= maxsite) %>% select(ptmacc, site, type) %>% unique %>% group_by(ptmacc, type) %>% tally %>% filter(n < minsize) %>% pull(ptmacc)))))
# Two-tailed, Wilcoxon rank sum test with continuity correction:
# >> W = 9955191152, p-value 1.532819e-54
# >> loc 4.999978 
# >> acc stratified: Z = 4.9651, p-value = 6.867e-07
# >> ptm|acc stratified: Z = 4.8913, p-value = 1.002e-06
q %>% filter(type == "ptm" & dis == "Disordered" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique %>% nrow
q %>% filter(type == "control" & dis == "Disordered" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique %>% nrow
# >> n(ptm)     >> 98,860
# >> n(control) >> 194,576

# Structured
wilcox.test(site ~ type, conf.int = T, q %>% filter(dis == "Structured" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique)
wilcox.test(site ~ type, conf.int = T, q %>% filter(dis == "Structured" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique)$p.value
wilcox_test(site ~ type, conf.int = T, q %>% filter(dis == "Structured" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique)
wilcox_test(site ~ type | acc, q %>% filter(dis == "Structured" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique %>% filter(!(acc %in% (q %>% filter(dis == "Structured" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique %>% group_by(acc, type) %>% tally %>% filter(n < minsize) %>% pull(acc)))))
wilcox_test(site ~ type | ptmacc, q %>% mutate(ptmacc = as.factor(paste(ptm, acc))) %>% filter(dis == "Structured" & site <= maxsite & len >= maxsite) %>% select(ptmacc, site, type) %>% unique %>% filter(!(ptmacc %in% (q %>% mutate(ptmacc = as.factor(paste(ptm, acc))) %>% filter(dis == "Structured" & site <= maxsite & len >= maxsite) %>% select(ptmacc, site, type) %>% unique %>% group_by(ptmacc, type) %>% tally %>% filter(n < minsize) %>% pull(ptmacc)))))
# W = 1.0233e+10, p-value = 5.609e-15
# loc -2.999999
# >> acc stratified: Z = -6.753, p-value = 1.449e-11
# >> ptm|acc stratified: Z = -8.1949, p-value = 2.508e-16
q %>% filter(type == "ptm" & dis == "Structured" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique %>% nrow
q %>% filter(type == "control" & dis == "Structured" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique %>% nrow
# >> n(ptm)     >> 77,927
# >> n(control) >> 267,534

# Buried
wilcox.test(site ~ type, conf.int = T, q %>% filter(dis == "Buried" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique)
wilcox.test(site ~ type, conf.int = T, q %>% filter(dis == "Buried" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique)$p.value
wilcox_test(site ~ type | acc, q %>% filter(dis == "Buried" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique %>% filter(!(acc %in% (q %>% filter(dis == "Buried" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique %>% group_by(acc, type) %>% tally %>% filter(n < minsize) %>% pull(acc)))))
wilcox_test(site ~ type | ptmacc, q %>% mutate(ptmacc = as.factor(paste(ptm, acc))) %>% filter(dis == "Buried" & site <= maxsite & len >= maxsite) %>% select(ptmacc, site, type) %>% unique %>% filter(!(ptmacc %in% (q %>% mutate(ptmacc = as.factor(paste(ptm, acc))) %>% filter(dis == "Buried" & site <= maxsite & len >= maxsite) %>% select(ptmacc, site, type) %>% unique %>% group_by(ptmacc, type) %>% tally %>% filter(n < minsize) %>% pull(ptmacc)))))
# W = 3770603049, p-value 8.40149e-59
# loc -7.000024
# >> acc stratified: Z = -12.941, p-value < 2.2e-16
# >> ptm|acc stratified: Z = -9.8482, p-value < 2.2e-16
q %>% filter(type == "ptm" & dis == "Buried" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique %>% nrow
q %>% filter(type == "control" & dis == "Buried" & site <= maxsite & len >= maxsite) %>% select(acc, site, type) %>% unique %>% nrow
# >> n(ptm)     >> 43,926
# >> n(control) >> 180,651
