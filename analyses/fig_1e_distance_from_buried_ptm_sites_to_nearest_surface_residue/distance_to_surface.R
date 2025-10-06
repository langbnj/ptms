blang_init()

tmpfile <- "q.rds"

# Enable pLDDT filtering:
# plddt_filtering <- 0 # Disable
plddt_filtering <- 70
# plddt_filtering <- 90

if (!file.exists(tmpfile)) {
  # Takes around 15 minutes to run (possibly more)
  # Full query (gives all unimod columns, but this means introducing some duplicate acc|sites)
  # Note: alphaseq match: OK! The query below ensures the alphasa and uniseq sequences match (since it only uses type='CoreSurf' entries from uniseq).
  # Note: residue uniqueness: OK! This script ensures each residue is only counted once. It also only uses canonical isoform residues (due to type='CoreSurf').
  # q <- Query("SELECT m.*, SUBSTRING(s.seq, m.site, 1), MIN(ABS(a.site - m.site)) AS surfdist FROM unimod m, uniseq s, alphasa a WHERE m.species='human' AND m.ptm IS NOT NULL AND m.acc=s.acc AND s.type='CoreSurf' AND SUBSTRING(s.seq, m.site, 1)='C' AND m.acc=a.acc AND a.surf='S' GROUP BY m.id")
  
  # Using alphamap:
  # q <- Query("SELECT m.*, a.plddt, a.plddt10, SUBSTRING(cs.seq, m.site, 1), MIN(ABS(a.site - m.site)) AS surfdist FROM unimod m, uniseq cs, alphamap am, alphasa a WHERE m.species='human' AND m.ptm IS NOT NULL AND m.acc=cs.acc AND cs.type='CoreSurf' AND SUBSTRING(cs.seq, m.site, 1)='C' AND m.acc=am.value AND am.type='uniprot' AND am.version='2022_04' AND am.best=1 AND am.map=a.acc AND am.afdb=a.afdb AND a.surf='S' GROUP BY m.id")
  q <- Query("SELECT m.*, ac.plddt AS plddt_c, ac.plddt10 AS plddt10_c, a.plddt AS plddt_s, a.plddt10 AS plddt10_s, SUBSTRING(cs.seq, m.site, 1), MIN(ABS(a.site - m.site)) AS surfdist FROM unimod m, uniseq cs, alphamap am, alphasa a, alphasa ac WHERE m.species='human' AND m.ptm IS NOT NULL AND m.acc=cs.acc AND cs.type='CoreSurf' AND SUBSTRING(cs.seq, m.site, 1)='C' AND m.acc=am.value AND am.type='uniprot' AND am.version='2022_04' AND am.best=1 AND am.map=a.acc AND am.afdb=a.afdb AND a.surf='S' AND am.map=ac.acc AND am.afdb=ac.afdb AND ac.site=m.site AND ac.surf='C' GROUP BY m.id")
  # # DISTINCT query (no duplicates, and only the columns needed)
  # q <- Query("SELECT DISTINCT m.acc, m.site, m.ptm, m.source, MIN(ABS(a.site - m.site)) AS surfdist FROM unimod m, uniseq s, alphasa a WHERE m.species='human' AND m.ptm IS NOT NULL AND m.acc=s.acc AND s.type='CoreSurf' AND SUBSTRING(s.seq, m.site, 1)='C' AND m.acc=a.acc AND a.surf='S'")
  saveRDS(q, tmpfile)
} else {
  # Milliseconds
  q <- readRDS(tmpfile)
  q0 <- q
  # Currently 156,988 rows, which is correct according to:
  # SELECT COUNT(*) FROM unimod m, uniseq s WHERE m.species='human' AND m.ptm IS NOT NULL AND m.acc=s.acc AND s.type='CoreSurf' AND SUBSTRING(s.seq, m.site, 1)='C';
  # >> Expecting 156,988 rows
}



# Run pLDDT filtering?
tmp_plddt <- ""
uniseq_type = "CoreSurf"
if (plddt_filtering != 0) {
  tmp_plddt <- f("-pLDDT{plddt_filtering}")
  uniseq_type = f("CoreSurf_pLDDT{plddt_filtering}")
  
  q %>% nrow
  # 158539
  q %>% filter(plddt_c >= plddt_filtering) %>% nrow
  # 147319
  q %>% filter(plddt_s >= plddt_filtering) %>% nrow
  # 11138 >> too few
  q %>% filter(plddt_c >= plddt_filtering & plddt_s >= plddt_filtering) %>% nrow
  # 10732 >> too few
  # >> Can't filter by plddt_s
  q %>% filter(plddt_c >= plddt_filtering) %>% nrow
  q %>% filter(plddt10_c >= plddt_filtering) %>% nrow
  # >> Only slightly fewer when filtering by plddt10_c than plddt_c, great
  q %>% filter(plddt_c >= plddt_filtering & plddt10_c >= plddt_filtering) %>% nrow
  # >> Only slightly fewer when filtering by both plddt_c and plddt10_c, great
  
  # Filter by core residue pLDDT ≥ plddt_filtering and pLDDT10 ≥ plddt_filtering. Don't filter by surface residue pLDDT (would only retain ~10% of cases).
  q %<>% filter(plddt_c >= plddt_filtering & plddt10_c >= plddt_filtering)
}





# Summarise by PTM
q %>% 
  group_by(ptm) %>%
  summarise(n = n_distinct(acc, site), min_surfdist = min(surfdist), mean_surfdist = mean(surfdist), median_surfdist = median(surfdist), max_surfdist = max(surfdist)) %>%
  arrange(desc(n))

# Summarise by PTM
q %>% 
  group_by(source) %>%
  summarise(n = n_distinct(acc, site), min_surfdist = min(surfdist), mean_surfdist = mean(surfdist), median_surfdist = median(surfdist), max_surfdist = max(surfdist)) %>%
  arrange(desc(n))

# Summarise by PTM and source
q %>% 
  group_by(ptm, source) %>%
  summarise(n = n_distinct(acc, site), min_surfdist = min(surfdist), mean_surfdist = mean(surfdist), median_surfdist = median(surfdist), max_surfdist = max(surfdist)) %>%
  arrange(desc(n))

# Summarise Ochoa
q %>% 
  filter(source %in% c("Ochoa")) %>%
  group_by(ptm, source) %>%
  summarise(n = n_distinct(acc, site), min_surfdist = min(surfdist), mean_surfdist = mean(surfdist), median_surfdist = median(surfdist), max_surfdist = max(surfdist)) %>%
  arrange(desc(n))

# How many Ochoa PTM sites have surfdist ≥ 2?
q %>% 
  filter(source %in% c("Ochoa")) %>%
  group_by(source, surfdist >= 2) %>%
  summarise(n = n_distinct(acc, site), min_surfdist = min(surfdist), mean_surfdist = mean(surfdist), median_surfdist = median(surfdist), max_surfdist = max(surfdist)) %>%
  arrange(desc(n))
# >> 2136/(2136+5104) = 29.5% are buried at least 2 aa from the nearest surface residue. That's still quite substantial (2000 fairly buried Ochoa sites).

# How many Ochoa PTM sites have surfdist ≥ 3?
q %>% 
  filter(source %in% c("Ochoa")) %>%
  group_by(source, surfdist >= 3) %>%
  summarise(n = n_distinct(acc, site), min_surfdist = min(surfdist), mean_surfdist = mean(surfdist), median_surfdist = median(surfdist), max_surfdist = max(surfdist)) %>%
  arrange(desc(n))
# >> 1058/(1058+6182) = 14.61% are buried at least 3 aa from the nearest surface residue. That's still quite substantial (1000 deeply buried Ochoa sites).

# How many Ochoa PTM sites have surfdist ≥ 2 (by PTM type)?
q %>% 
  filter(source %in% c("Ochoa")) %>%
  group_by(ptm, source, surfdist >= 2) %>%
  summarise(n = n_distinct(acc, site), min_surfdist = min(surfdist), mean_surfdist = mean(surfdist), median_surfdist = median(surfdist), max_surfdist = max(surfdist)) %>%
  arrange(desc(n))

# How many Ochoa PTM sites have surfdist ≥ 3 (by PTM type)?
q %>% 
  filter(source %in% c("Ochoa")) %>%
  group_by(ptm, source, surfdist >= 3) %>%
  summarise(n = n_distinct(acc, site), min_surfdist = min(surfdist), mean_surfdist = mean(surfdist), median_surfdist = median(surfdist), max_surfdist = max(surfdist)) %>%
  arrange(desc(n))

# How many Ochoa PTM sites have surfdist ≥ 5 (by PTM type)?
q %>% 
  filter(source %in% c("Ochoa")) %>%
  group_by(ptm, source, surfdist >= 5) %>%
  summarise(n = n_distinct(acc, site), min_surfdist = min(surfdist), mean_surfdist = mean(surfdist), median_surfdist = median(surfdist), max_surfdist = max(surfdist)) %>%
  arrange(desc(n))

# How many Ochoa PTM sites have surfdist ≥ 10 (by PTM type)?
q %>% 
  filter(source %in% c("Ochoa")) %>%
  group_by(ptm, source, surfdist >= 10) %>%
  summarise(n = n_distinct(acc, site), min_surfdist = min(surfdist), mean_surfdist = mean(surfdist), median_surfdist = median(surfdist), max_surfdist = max(surfdist)) %>%
  arrange(desc(n))

# How many PTM sites have surfdist ≥ 2?
q %>% 
  group_by(surfdist >= 2) %>%
  summarise(n = n_distinct(acc, site), min_surfdist = min(surfdist), mean_surfdist = mean(surfdist), median_surfdist = median(surfdist), max_surfdist = max(surfdist)) %>%
  arrange(desc(n))
# >> 43286/(43286+68088) = 38.87% are buried at least 2 aa from the nearest surface residue

# How many PTM sites have surfdist ≥ 3?
q %>% 
  group_by(surfdist >= 3) %>%
  summarise(n = n_distinct(acc, site), min_surfdist = min(surfdist), mean_surfdist = mean(surfdist), median_surfdist = median(surfdist), max_surfdist = max(surfdist)) %>%
  arrange(desc(n))
# >> 24249/(24249+87125) = 21.77% are buried at least 3 aa from the nearest surface residue

# Plot by PTM and source
q %>%
  distinct(ptm, source, acc, site, surfdist) %>%
  ggplot(aes(x = surfdist)) +
  geom_bar(stat = "count") +
  facet_grid(rows = vars(ptm), cols = vars(source), scales = "free_y") +
  theme_minimal()

# Plot Ochoa
q %>%
  filter(source %in% c("Ochoa")) %>%
  distinct(ptm, source, acc, site, surfdist) %>%
  ggplot(aes(x = surfdist)) +
  geom_bar(stat = "count") +
  facet_grid(rows = vars(ptm), cols = vars(source)) +
  theme_minimal()

# Plot
q %>%
  distinct(acc, site, surfdist) %>%
  ggplot(aes(x = surfdist)) +
  geom_bar(stat = "count") +
  # geom_histogram() +
  # scale_x_continuous(limits = c(0, max(q$surfdist)), breaks = pretty_breaks(max(q$surfdist)), expand = c(0, 0)) +
  scale_x_continuous(breaks = pretty_breaks(8)) +
  # geom_vline(xintercept = 3.5, linetype = "dashed") +
  # geom_vline(xintercept = 3.5, linetype = "dotted") +
  theme_minimal() +
  xlab("Distance to nearest surface residue (amino acids)") +
  ylab("PTM sites")
ggsave("output-surfdist.pdf", width=5, height=4)
  
# Plot with max 10
q %>%
  distinct(acc, site, surfdist) %>%
  mutate(surfdist = ifelse(surfdist > 10, ">10", surfdist)) %>%
  # pull(surfdist) %>% as.factor %>% summary
  ggplot(aes(x = fct_relevel(surfdist, c(1:10, ">10")))) +
  geom_bar(stat = "count") +
  # scale_x_continuous(limits = c(0, max(q$surfdist)), breaks = pretty_breaks(max(q$surfdist)), expand = c(0, 0)) +
  # scale_x_discrete(breaks = c(as.character(1:9), "≥10"), limits = c(as.character(1:9), "≥10")) +
  # geom_vline(xintercept = 3.5, linetype = "dashed") +
  # geom_vline(xintercept = 3.5, linetype = "dotted") +
  theme_minimal() +
  xlab("Distance to nearest surface residue (amino acids)") +
  ylab("Buried human PTM sites")
# ggsave("output-surfdist.pdf", width=89, height=89, units = "mm")
ggsave("output-surfdist-max10.pdf", width = 5, height = 4)

       


# ALL RESIDUES, rather than only PTM sites:
# # Alternative, based on pipeline-uniseq_coresurf_corelen:
# corelen <- tibble(read_tsv("tmp-corelengths-human.txt"))
# corelen %<>% mutate(len = as.integer(len))
# corelen
# Now replaced by query below:
# # All proteins
# qtmp <- Query("SELECT seq FROM uniseq WHERE species='human' AND type='CoreSurf'")
# qtmp
# Only PTM-containing proteins
qtmp <- Query(f("SELECT s.acc, s.seq FROM unimod m, uniseq s WHERE m.species='human' AND m.ptm IS NOT NULL AND m.acc=s.acc AND s.type='{uniseq_type}' GROUP BY s.acc"))
qtmp
qtmp %<>%
  mutate(c = str_extract_all(seq, "C+")) %>%
  select(-seq) %>%
  unnest_longer(c) %>%
  mutate(len = str_length(c)) %>%
  select(-c)
# identical(corelen, qtmp)
# identical(sort(corelen$len), sort(qtmp$len))
qtmp
corelen <- qtmp

summary(corelen)
corelen %>% tally
corelen %>%
  ggplot(aes(x = len)) +
  geom_bar(stat = "count") +
  theme_minimal() +
  xlab("Length of buried residue stretch (amino acids)") +
  ylab("Residues")
ggsave("output-corelen-including_non_ptm.pdf", width = 5, height = 4)

corelen %>%
  mutate(len = ifelse(len > 10, ">10", len)) %>%
  ggplot(aes(x = fct_relevel(len, c(1:10, ">10")))) +
  geom_bar(stat = "count") +
  theme_minimal() +
  xlab("Length of buried residue stretch (amino acids)") +
  ylab("Residues")
ggsave("output-corelen-including_non_ptm-max10.pdf", width = 5, height = 4)

corelen %>%
  mutate(len = ifelse(len > 20, ">20", len)) %>%
  ggplot(aes(x = fct_relevel(len, c(1:20, ">20")))) +
  geom_bar(stat = "count") +
  theme_minimal() +
  xlab("Length of buried residue stretch (amino acids)") +
  ylab("Residues")
ggsave("output-corelen-including_non_ptm-max20.pdf", width = 5, height = 4)

# Calculate sum of the surfdists in each stretch
corelen
# 1 C -> 1 = 1
# 2 CC -> 1+1 = 2
# 3 CCC -> 1+2+1 = 4
# 4 CCCC -> 1+2+2+1 = 6
# 5 CCCCC -> 1+2+3+2+1 = 9
# 6 CCCCCC -> 1+2+3+3+2+1 = 12
# 7 CCCCCCC -> 1+2+3+4+3+2+1 = 16
# 8 CCCCCCCC -> 1+2+3+4+4+3+2+1 = 20
# 9 CCCCCCCCC -> 1+2+3+4+5+4+3+2+1 = 25

# # Python function for calculating sum:
# def sum_c(n):
#   sum = 0
# for i in range(1, n+1):
#   if i % 2 == 0:
#   # even: add half
#   sum += i//2
# else:
#   # odd: add half, rounded up
#   sum += (i+1)//2
# return sum

surfdists_sum <- function(n) {
  sum <- 0
  for (i in 1:n) {  # Iterate from 1 to n (inclusive)
    if (i %% 2 == 0) {  # Check for even numbers
      sum <- sum + i %/% 2  # Add half (integer division)
    } else {  # Odd numbers
      sum <- sum + (i + 1) %/% 2  # Add half, rounded up
    }
  }
  return(sum)  # Return the final sum
}

surfdists_sum(1)
surfdists_sum(2)
surfdists_sum(3)
surfdists_sum(4)
surfdists_sum(5)
surfdists_sum(6)
surfdists_sum(7)
surfdists_sum(8)
surfdists_sum(9)

# Calculate surfdists_sum column
corelen$surfdists_sum <- NULL
corelen
corelen %>% head(n = 1000) %>% mutate(surfdists_sum = map_int(len, surfdists_sum))
corelen %<>% mutate(surfdists_sum = map_int(len, surfdists_sum))
corelen

# Calculate average surfdist for each region
corelen %<>% mutate(surfdists_avg = surfdists_sum / len)
corelen

corelen %>%
  ggplot(aes(x = surfdists_avg)) +
  geom_histogram(binwidth = 1) +
  theme_minimal() +
  xlab("Average distance to nearest surface residue (amino acids)") +
  ylab("Buried regions")
ggsave("output-surfdists_avg-including_non_ptm-float.pdf", width = 5, height = 4)

corelen %>%
  mutate(surfdists_avg = as.integer(round(surfdists_avg))) %>%
  ggplot(aes(x = surfdists_avg)) +
  geom_bar(stat = "count") +
  theme_minimal() +
  xlab("Average distance to nearest surface residue (amino acids)") +
  ylab("Buried regions")
ggsave("output-surfdists_avg-including_non_ptm.pdf", width = 5, height = 4)

corelen %>%
  mutate(surfdists_avg = as.integer(round(surfdists_avg))) %>%
  mutate(surfdists_avg = ifelse(surfdists_avg > 10, ">10", surfdists_avg)) %>%
  ggplot(aes(x = fct_relevel(surfdists_avg, c(1:10, ">10")))) +
  geom_bar(stat = "count") +
  theme_minimal() +
  xlab("Average distance to nearest surface residue (amino acids)") +
  ylab("Buried regions")
ggsave("output-surfdists_avg-including_non_ptm-max10.pdf", width = 5, height = 4)

# Instead of using averages per Core region (which skews towards 1), get per-residue values so I can really compare the PTM distribution vs. all residues
surfdists_list <- function(n) {
  res <- list()
  for (i in 1:n) {  # Iterate from 1 to n (inclusive)
    if (i %% 2 == 0) {  # Check for even numbers
      res <- append(res, i %/% 2)  # Add half (integer division)
    } else {  # Odd numbers
      res <- append(res, (i + 1) %/% 2)  # Add half, rounded up
    }
  }
  return(res)  # Return the final sum
}
surfdists_list(1)
surfdists_list(2)
surfdists_list(3)
surfdists_list(4)
surfdists_list(5)
surfdists_list(6)
surfdists_list(7)
surfdists_list(8)
surfdists_list(9)

# Calculate surfdists_list column
corelen$surfdists_list <- NULL
corelen
corelen %<>% mutate(surfdists_list = map(len, surfdists_list))
corelen %<>% unnest_longer(surfdists_list)
corelen %>% head(n = 20)

corelen %>%
  mutate(surfdists_list = as.integer(round(surfdists_list))) %>%
  ggplot(aes(x = surfdists_list)) +
  geom_bar(stat = "count") +
  theme_minimal() +
  xlab("Distance to nearest surface residue (amino acids)") +
  ylab("Buried residues")
ggsave("output-surfdists-including_non_ptm.pdf", width = 5, height = 4)

corelen %>%
  mutate(surfdists_list = as.integer(round(surfdists_list))) %>%
  mutate(surfdists_list = ifelse(surfdists_list > 10, ">10", surfdists_list)) %>%
  ggplot(aes(x = fct_relevel(surfdists_list, c(1:10, ">10")))) +
  geom_bar(stat = "count") +
  theme_minimal() +
  xlab("Distance to nearest surface residue (amino acids)") +
  ylab("Buried residues")
ggsave("output-surfdists-including_non_ptm-max10.pdf", width = 5, height = 4)



# Make box plot comparing PTM distance distribution to background
# Merge q (PTM sites) and corelen (background) tibbles
q
q %>% distinct(acc, site, ptm, source, surfdist)
q %>% distinct(acc, site, surfdist)
corelen
wilcox.test(q %>% pull(surfdist), corelen$surfdists_list)$p.value
# >> p = 3.67e-20 (grouped by unimod row)
wilcox.test(q %>% distinct(acc, site, ptm, source, surfdist) %>% pull(surfdist), corelen$surfdists_list)$p.value
# >> p = 1.13e-21 (grouped by acc|site|ptm|source)
wilcox.test(q %>% distinct(acc, site, surfdist) %>% pull(surfdist), corelen$surfdists_list)$p.value
# >> p = 0.03892094 (grouped by acc|site only, i.e. no residue redundancy)
# qptm <- tibble(type = "ptm", surfdist = q$surfdist)
qptm <- tibble(type = "ptm", surfdist = q %>% distinct(acc, site, surfdist) %>% pull(surfdist))
qptm
qptm_plddt70 <- tibble(type = "ptm", surfdist = q %>% filter(plddt_c >= 70 & plddt_s >= 70) %>% distinct(acc, site, surfdist) %>% pull(surfdist))
qptm_plddt70
qbackground <- tibble(type = "background", surfdist = corelen$surfdists_list)
qbackground
qdists <- bind_rows(qptm, qbackground)
qdists

# Add zeroes (surface residues) for PTM sites
qtmp_surface_ptm <- Query(f("SELECT 'ptm' AS type, 0 AS surfdist FROM unimod m, uniseq s WHERE m.species='human' AND m.ptm IS NOT NULL AND m.acc=s.acc AND s.type='{uniseq_type}' AND SUBSTRING(s.seq, m.site, 1)='S' GROUP BY m.acc, m.site"))
qtmp_surface_ptm
qtmp_surface_ptm_plddt70 <- Query(f("SELECT 'ptm' AS type, 0 AS surfdist FROM unimod m, uniseq s, alphamap am, alphasa a WHERE m.species='human' AND m.ptm IS NOT NULL AND m.acc=s.acc AND s.type='{uniseq_type}' AND SUBSTRING(s.seq, m.site, 1)='S' AND m.acc=am.value AND am.type='uniprot' AND am.version='2022_04' AND am.best=1 AND am.map=a.acc AND am.afdb=a.afdb AND m.site=a.site AND a.plddt >= 70 GROUP BY m.acc, m.site"))
qtmp_surface_ptm_plddt70
# Add zeroes (surface residues) for background (PTM-containing proteins only)
qtmp_surface_background_seqs <- Query(f("SELECT s.seq FROM unimod m, uniseq s WHERE m.species='human' AND m.ptm IS NOT NULL AND m.acc=s.acc AND s.type='{uniseq_type}' GROUP BY s.acc"))
qtmp_surface_background_seqs
qtmp_surface_background_seqs %>%
  mutate(s = str_extract_all(seq, "S")) %>%
  select(s) %>%
  unnest_longer(s) %>%
  transmute(type = "background", surfdist = 0) -> qtmp_surface_background
qtmp_surface_background

qdists <- bind_rows(qtmp_surface_ptm, qptm, qbackground, qtmp_surface_background)
qdists
qdists_plddt70 <- bind_rows(qptm_plddt70, qbackground)
qdists_plddt70

wilcox.test(qdists %>% filter(type == "ptm") %>% filter(surfdist > 0) %>% pull(surfdist),
            qdists %>% filter(type == "background") %>% filter(surfdist > 0) %>% pull(surfdist))$p.value -> pvalue_buried_only
pvalue_buried_only
# >> p = 0.03892094 (when filtering out surfdist = 0, i.e. Surface residues)
wilcox.test(qdists %>% filter(type == "ptm") %>% pull(surfdist),
            qdists %>% filter(type == "background") %>% pull(surfdist))$p.value -> pvalue_including_surface
pvalue_including_surface
# >> p = 0 (when including surface residues)


# qdists %>%
#   ggplot(aes(x = surfdist, fill = type)) +
#   geom_bar(stat = "count") +
#   facet_grid(rows = vars(type), scales = "free_y") +
#   theme_minimal()
# 
# qdists %>%
#   ggplot(aes(x = surfdist, fill = type)) +
#   geom_bar(stat = "count", position = "dodge") +
#   theme_minimal()
# 
# qdists %>% 
#   group_by(type) %>% 
#   tally
#   
# # qdists$relsurfdist <- NULL
# qdists %<>% 
#   group_by(type) %>% 
#   mutate(relsurfdist = surfdist / n())
#   # filter(type == "ptm") %>% mutate(relsurfdist_test = relsurfdist * 156988)
#   # >> OK, n() works correctly within this group
#   # filter(type == "background") %>% mutate(relsurfdist_test = relsurfdist * 3348120)
#   # >> OK, n() works correctly within this group
# 
# qdists
# qdists %>%
#   ggplot(aes(x = relsurfdist, fill = type)) +
#   geom_histogram(position = "dodge") +
#   theme_minimal()
# 
# qdists %>% 
#   group_by(type, surfdist) %>% 
#   tally
# 
# qdists %<>% 
#   group_by(surfdist) %>% 
#   mutate(relsurfdist = surfdist / n())
# qdists
# qdists %>%
#   ggplot(aes(x = relsurfdist, fill = type)) +
#   geom_histogram(position = "dodge") +
#   theme_minimal()
# 
# qdists %>% 
#   group_by(type, surfdist) %>% 
#   tally %>%
#   ggplot(aes(x = surfdist, y = n, fill = type)) +
#   geom_col(position = "dodge") +
#   theme_minimal()

qdists %>% 
  filter(surfdist > 0) %>%
  group_by(type) %>% 
  tally
# >> background 3348120
# >> ptm         111374

qdists %>% 
  # filter(surfdist > 0) %>%
  group_by(type, surfdist) %>%
  tally %>%
  arrange(surfdist, type)
# >> background 3348120
# >> ptm         111374

# Bar plots
qdists %>% 
  filter(surfdist > 0) %>%
  group_by(type, surfdist) %>%
  tally %>%
  group_by(type) %>%
  # mutate(rel_n = n / n()) %>%
  mutate(rel_n = n / sum(n)) %>%
  arrange(surfdist, type) %>%
  ggplot(aes(x = surfdist, y = rel_n, fill = type, colour = type, alpha = 0.2)) +
  geom_col(position = "dodge") +
  # geom_line() +
  # geom_area(outline.type = "full") +
  # geom_area(stat = "identity") +
  scale_x_continuous(breaks = pretty_breaks(10)) +
  scale_colour_manual(values = ptmcol, aesthetics = c("colour", "fill")) +
  theme_minimal() +
  ggtitle(f("Buried residues only: p = ", wilcox.test(qdists %>% filter(type == "ptm") %>% filter(surfdist > 0) %>% pull(surfdist),
                                                       qdists %>% filter(type == "background") %>% filter(surfdist > 0) %>% pull(surfdist))$p.value))
ggsave("output-distance_to_surface-bar-buried_only.pdf", width = 5, height = 4)

qdists %>% 
  # filter(surfdist > 0) %>%
  group_by(type, surfdist) %>%
  tally %>%
  group_by(type) %>%
  # mutate(rel_n = n / n()) %>%
  mutate(rel_n = n / sum(n)) %>%
  arrange(surfdist, type) %>%
  ggplot(aes(x = surfdist, y = rel_n, fill = type, colour = type, alpha = 0.2)) +
  geom_col(position = "dodge") +
  # geom_line() +
  # geom_area(outline.type = "full") +
  # geom_area(stat = "identity") +
  scale_x_continuous(breaks = pretty_breaks(10)) +
  scale_colour_manual(values = ptmcol, aesthetics = c("colour", "fill")) +
  theme_minimal() +
  ggtitle(f("Including surface residues: p = ", wilcox.test(qdists %>% filter(type == "ptm") %>% pull(surfdist),
                                                           qdists %>% filter(type == "background") %>% pull(surfdist))$p.value))
ggsave("output-distance_to_surface-bar-including_surface.pdf", width = 5, height = 4)


qdists %>% 
  # filter(surfdist > 0) %>%
  group_by(type, surfdist) %>%
  ggplot(aes(x = surfdist, fill = type, colour = type, alpha = 0.2)) +
  geom_density(bw = 1/3) +
  # geom_histogram(position = "dodge") +
  # geom_area(outline.type = "full") +
  # geom_area(stat = "identity") +
  scale_x_continuous(breaks = pretty_breaks(10)) +
  # facet_grid(rows = vars(type), scales = "free_y") +
  scale_colour_manual(values = ptmcol, aesthetics = c("colour", "fill")) +
  coord_cartesian(xlim = c(0, 10)) +
  theme_minimal() +
  ggtitle(f("Including surface residues: p = {pvalue_including_surface}"))
ggsave("output-distance_to_surface-density-including_surface-zoom10.pdf", width = 5, height = 4)

# Density plots
qdists %>% 
  filter(surfdist > 0) %>%
  group_by(type, surfdist) %>%
  ggplot(aes(x = surfdist, fill = type, colour = type, alpha = 0.2)) +
  geom_density(bw = 1/3) +
  # geom_histogram(position = "dodge") +
  # geom_area(outline.type = "full") +
  # geom_area(stat = "identity") +
  scale_x_continuous(breaks = pretty_breaks(10)) +
  # facet_grid(rows = vars(type), scales = "free_y") +
  scale_colour_manual(values = ptmcol, aesthetics = c("colour", "fill")) +
  coord_cartesian(xlim = c(0, 10)) +
  theme_minimal() +
  ggtitle(f("Buried residues only: p = {pvalue_buried_only}"))
ggsave("output-distance_to_surface-density-buried_only-zoom10.pdf", width = 5, height = 4)

# Faceted bar plot
# tmpmax <- 5
tmpmax <- 7
(qdists %>% 
  filter(surfdist > 0) %>%
  group_by(type, surfdist) %>%
  mutate(surfdist = as.character(ifelse(surfdist >= tmpmax, f("≥{tmpmax}"), surfdist))) %>%
  ggplot(aes(x = fct_relevel(surfdist, c(1:(tmpmax - 1), f("≥{tmpmax}"))), fill = type, colour = type, alpha = 0.2)) +
  # geom_histogram(position = "dodge") +
  # geom_area(outline.type = "full") +
  # geom_area(stat = "identity") +
  # geom_density(bw = 1/3) +
  geom_bar(stat = "count") +
  # scale_x_continuous(breaks = pretty_breaks(10)) +
  # facet_grid(rows = vars(type), scales = "free_y") +
  scale_colour_manual(values = ptmcol, aesthetics = c("colour", "fill")) +
  # coord_cartesian(xlim = c(0, 10)) +
  facet_grid(rows = vars(type), scales = "free_y") +
  theme_nature() +
  ggtitle(f("Buried residues only: p = {pvalue_buried_only}")) +
  xlab("surfdist")
) %>% qsave(f("output-distance_to_surface-facetbar-buried_only-zoom{tmpmax}{tmp_plddt}.pdf"), height = 80, width = 80)

# Dodged bar plot
# tmpmax <- 5
tmpmax <- 7
(
  qdists %>% 
  filter(surfdist > 0) %>%
  group_by(type, surfdist) %>%
  mutate(surfdist = as.character(ifelse(surfdist >= tmpmax, f("≥{tmpmax}"), surfdist))) %>%
  group_by(type, surfdist) %>%
  tally %>%
  mutate(freq = n / sum(n)) %>%
  mutate(type = replace(type, type == "ptm", "PTM")) %>%
  mutate(type = replace(type, type == "background", "Control")) %>%
  mutate(type = factor(type, levels = c("PTM", "Control"))) %>%
  # ggplot(aes(x = fct_relevel(surfdist, c(1:(tmpmax - 1), f("≥{tmpmax}"))), y = freq, fill = type, colour = type, alpha = 0.2)) +
  ggplot(aes(x = fct_relevel(surfdist, c(1:(tmpmax - 1), f("≥{tmpmax}"))), y = freq, fill = type, colour = type)) +
  # geom_histogram(position = "dodge") +
  # geom_area(outline.type = "full") +
  # geom_area(stat = "identity") +
  # geom_density(bw = 1/3) +
  geom_col(position = "dodge", alpha = 0.5) +
  # scale_x_continuous(breaks = pretty_breaks(10)) +
  scale_y_continuous(expand = c(0, 0)) +
  # facet_grid(rows = vars(type), scales = "free_y") +
  scale_colour_manual(values = ptmcol, aesthetics = c("colour", "fill"), name = NULL) +
  # coord_cartesian(xlim = c(0, 10)) +
  # facet_grid(rows = vars(type), scales = "free_y") +
  guides(alpha = "none") +
  theme_nature(legend_position = "top") +
  ggtitle(f("Buried residues only: p = {pvalue_buried_only}")) +
  xlab("surfdist")
) %>% qsave(f("output-distance_to_surface-dodgebar-buried_only-zoom{tmpmax}{tmp_plddt}.pdf"))

# Line plot

# # Statistics (no stratification)
# wilcox.test(surfdist ~ type, qdists %>% filter(surfdist > 0))$p.value
# # >> p = 0.005594778 (identical to pvalue_buried_only above)
# tmppval <- wilcox.test(surfdist ~ type, qdists %>% filter(surfdist > 0))$p.value
# wilcox.test(surfdist ~ type, qdists)$p.value
# # >> p = ~0
# qdists %>% group_by(type) %>% filter(surfdist > 0) %>% summarise(median(surfdist), mean(surfdist))
# # >> Buried only:       ptm 1.93  control 1.98
# qdists %>% group_by(type) %>% summarise(median(surfdist), mean(surfdist))
# # >> Including surface: ptm 0.339 control 0.620
# ptmdist <- qdists %>% group_by(type) %>% filter(surfdist > 0 & type == "ptm") %>% pull(surfdist) %>% mean
# controldist <- qdists %>% group_by(type) %>% filter(surfdist > 0 & type == "background") %>% pull(surfdist) %>% mean
# tmpdelta <- ptmdist - controldist

# Statistics (can't stratify, but trying poisson.test)
qdists
qdists %>% group_by(type) %>% tally
wilcox.test(surfdist ~ type, conf.int=T, qdists %>% filter(surfdist > 0))
wilcox.test(surfdist ~ type, conf.int=T, qdists %>% filter(surfdist > 0))$p.value
# poisson.test(surfdist ~ type, qdists %>% filter(surfdist > 0))
t.test(surfdist ~ type, conf.int=T, qdists %>% filter(surfdist > 0))
library(coin)
wilcox_test(surfdist ~ as.factor(type), qdists %>% filter(surfdist > 0))
# ks.test(surfdist ~ type, qdists %>% filter(surfdist > 0))
# wilcox_test(surfdist ~ type | acc, conf.int=T, qdists %>% filter(surfdist > 0))
# >> p = 0.005594778 (identical to pvalue_buried_only above)
wilcox.test(surfdist ~ type, qdists %>% filter(surfdist > 0))
tmppval <- wilcox.test(surfdist ~ type, qdists %>% filter(surfdist > 0))$p.value
tmppval
wilcox.test(surfdist ~ type, qdists)$p.value
# >> p = ~0
qdists %>% group_by(type) %>% filter(surfdist > 0) %>% summarise(median(surfdist), mean(surfdist))
# >> Buried only:       ptm 1.93  control 1.98
qdists %>% group_by(type) %>% summarise(median(surfdist), mean(surfdist))
# >> Including surface: ptm 0.339 control 0.620
ptmdist <- qdists %>% group_by(type) %>% filter(surfdist > 0 & type == "ptm") %>% pull(surfdist) %>% mean
controldist <- qdists %>% group_by(type) %>% filter(surfdist > 0 & type == "background") %>% pull(surfdist) %>% mean
tmpdelta <- ptmdist - controldist

# tmpmax <- 5
tmpmax <- 7
(
  qdists %>% 
  filter(surfdist > 0) %>%
  mutate(surfdist = replace(surfdist, surfdist > tmpmax, tmpmax)) %>%
  group_by(type, surfdist) %>%
  # mutate(surfdist = as.character(ifelse(surfdist >= tmpmax, f("≥{tmpmax}"), surfdist))) %>%
  group_by(type, surfdist) %>%
  tally %>%
  mutate(freq = n / sum(n)) %>%
  mutate(type = replace(type, type == "ptm", "Modified")) %>%
  mutate(type = replace(type, type == "background", "Control")) %>%
  mutate(type = factor(type, levels = c("Modified", "Control"))) %>%
  # ggplot(aes(x = fct_relevel(surfdist, c(1:(tmpmax - 1), f("≥{tmpmax}"))), y = freq, fill = type, colour = type, alpha = 0.2)) +
  ggplot(aes(x = surfdist, ymin = 0, ymax = freq, fill = fct_rev(type), colour = fct_rev(type))) +
  # geom_histogram(position = "dodge") +
  # geom_area(outline.type = "full") +
  # geom_area(stat = "identity") +
  # geom_density(bw = 1/3) +
  # geom_col(position = "dodge", alpha = 0.5) +
  # geom_line() +
  geom_ribbon(alpha = 0.3) +
  # # annotate(geom = "text", x = 7, y = 0.45, hjust = 1, vjust = 1, label = f("p = {round(tmppval, 3)} Δµ[D]{round(ptmdif, 2)}"), parse = T, size = 5/.pt) +
  # # annotate(geom = "text", x = 7, y = 0.45, hjust = 1, vjust = 1, label = f("paste(Delta, mu[D]) == plain({round(tmpdelta, 3)})"), parse = T, size = 5/.pt) +
  # # annotate(geom = "text", x = 7, y = 0.45, hjust = 1, vjust = 1, label = f("µ(PTM) = {round(ptmdist, 2)}\nµ(Control) = {round(controldist, 2)}"), size = 5/.pt) +
  # annotate(geom = "text", x = 7, y = 0.45, hjust = 1, vjust = 1, label = f("µ(Control) = {round(controldist, 2)}\nµ(PTM) = {round(ptmdist, 2)}\np = {round(tmppval, 3)}"), size = 5/.pt) +
  # scale_x_continuous(expand = c(0, 0), breaks = pretty_breaks(10)) +
  scale_x_continuous(expand = c(0, 0), breaks = 1:tmpmax, labels = c(1:(tmpmax - 1), f("≥{tmpmax}"))) +
  scale_y_continuous(expand = c(0, 0), breaks = pretty_breaks(3)) +
  # facet_grid(rows = vars(type), scales = "free_y") +
  # scale_colour_manual(values = c("PTM" = ptmb, "Control" = ptmg), aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
  # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
  scale_colour_manual(values = ptmcol, aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
  coord_cartesian(xlim = c(1, tmpmax)) +
  # facet_grid(rows = vars(type), scales = "free_y") +
  guides(alpha = "none") +
  theme_nature(legend_position = "top", extra_margin_right = 2) +
  theme(axis.title.y = element_text(hjust = 0)) +
  # theme(legend.key.size = unit(1, "mm")) +
  # ggtitle(f("Buried residues only: p = {pvalue_buried_only}")) +
  labs(tag = "e") +
  xlab("Distance to nearest surface residue") +
  ylab("Fraction")
) %>% qsave(f("output-distance_to_surface-line-buried_only-zoom{tmpmax}{tmp_plddt}.pdf"), height = 20)



# tmpmax <- 5
tmpmax <- 7
(
  qdists %>% 
  filter(surfdist > 0) %>%
  mutate(surfdist = replace(surfdist, surfdist > tmpmax, tmpmax)) %>%
  group_by(type, surfdist) %>%
  # mutate(surfdist = as.character(ifelse(surfdist >= tmpmax, f("≥{tmpmax}"), surfdist))) %>%
  group_by(type, surfdist) %>%
  tally %>%
  mutate(freq = n / sum(n)) %>%
  mutate(type = replace(type, type == "ptm", "Modified")) %>%
  mutate(type = replace(type, type == "background", "Control")) %>%
  mutate(type = factor(type, levels = c("Modified", "Control"))) %>%
  # ggplot(aes(x = fct_relevel(surfdist, c(1:(tmpmax - 1), f("≥{tmpmax}"))), y = freq, fill = type, colour = type, alpha = 0.2)) +
  ggplot(aes(x = surfdist, ymin = 0, ymax = freq, fill = fct_rev(type), colour = fct_rev(type))) +
  # geom_histogram(position = "dodge") +
  # geom_area(outline.type = "full") +
  # geom_area(stat = "identity") +
  # geom_density(bw = 1/3) +
  # geom_col(position = "dodge", alpha = 0.5) +
  # geom_line() +
  geom_ribbon(alpha = 0.3) +
  # annotate(geom = "text", x = 7, y = 0.45, hjust = 1, vjust = 1, label = f("p = {round(tmppval, 3)} Δµ[D]{round(ptmdif, 2)}"), parse = T, size = 5/.pt) +
  # annotate(geom = "text", x = 7, y = 0.45, hjust = 1, vjust = 1, label = f("paste(Delta, mu[D]) == plain({round(tmpdelta, 3)})"), parse = T, size = 5/.pt) +
  # annotate(geom = "text", x = 7, y = 0.45, hjust = 1, vjust = 1, label = f("µ(PTM) = {round(ptmdist, 2)}\nµ(Control) = {round(controldist, 2)}"), size = 5/.pt) +
  # annotate(geom = "text", x = 7, y = 0.45, hjust = 1, vjust = 1, label = f("µ(Control) = {round(controldist, 2)}\nµ(PTM) = {round(ptmdist, 2)}\np = {round(tmppval, 3)}"), size = 5/.pt) +
  # annotate(geom = "text", x = 7, y = 0.45, hjust = 1, vjust = 1, label = f("µ(Control) = {round(controldist, 2)}\nµ(PTM) = {round(ptmdist, 2)}\np = ", formatC(tmppval, format = "g", digits = 1)), size = 5/.pt) +
  annotate(geom = "text", x = 7, y = 0.45, hjust = 1, vjust = 1, label = f("µ(Control) = {round(controldist, 2)}\nµ(PTM) = {round(ptmdist, 2)}\np = 1e-7"), size = 5/.pt) +
  # scale_x_continuous(expand = c(0, 0), breaks = pretty_breaks(10)) +
  scale_x_continuous(expand = c(0, 0), breaks = 1:tmpmax, labels = c(1:(tmpmax - 1), f("≥{tmpmax}"))) +
  scale_y_continuous(expand = c(0, 0), breaks = pretty_breaks(3)) +
  # facet_grid(rows = vars(type), scales = "free_y") +
  # scale_colour_manual(values = c("PTM" = ptmb, "Control" = ptmg), aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
  # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
  scale_colour_manual(values = ptmcol, aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
  coord_cartesian(xlim = c(1, tmpmax)) +
  # facet_grid(rows = vars(type), scales = "free_y") +
  guides(alpha = "none") +
  theme_nature(legend_position = "top", extra_margin_right = 2) +
  theme(axis.title.y = element_text(hjust = 0)) +
  # theme(legend.key.size = unit(1, "mm")) +
  # ggtitle(f("Buried residues only: p = {pvalue_buried_only}")) +
  labs(tag = "e") +
  xlab("Distance to nearest surface residue") +
  ylab("Fraction")
) %>% qsave(f("output-distance_to_surface-line-buried_only-zoom{tmpmax}{tmp_plddt}-with_text.pdf"), height = 20)

# Box plot
(
  qdists %>% 
  filter(surfdist > 0) %>%
  mutate(type = replace(type, type == "ptm", "PTM")) %>%
  mutate(type = replace(type, type == "background", "Control")) %>%
  mutate(type = factor(type, levels = c("PTM", "Control"))) %>%
  # ggplot(aes(x = fct_relevel(surfdist, c(1:(tmpmax - 1), f("≥{tmpmax}"))), y = freq, fill = type, colour = type, alpha = 0.2)) +
  ggplot(aes(x = type, y = surfdist, fill = type, colour = type, alpha = 0.2)) +
  # geom_histogram(position = "dodge") +
  # geom_area(outline.type = "full") +
  # geom_area(stat = "identity") +
  # geom_density(bw = 1/3) +
  # geom_col(position = "dodge", alpha = 0.5) +
  geom_boxplot(notch = T, outlier.shape = NA) +
  # geom_violin() +
  # scale_x_continuous(breaks = pretty_breaks(10)) +
  scale_y_continuous(expand = c(0, 0)) +
  # facet_grid(rows = vars(type), scales = "free_y") +
  scale_colour_manual(values = ptmcol, aesthetics = c("colour", "fill"), name = NULL) +
  coord_cartesian(ylim = c(0, 10)) +
  # facet_grid(rows = vars(type), scales = "free_y") +
  # guides(alpha = "none") +
  theme_nature(legend_position = "top") +
  ggtitle(f("Buried residues only: p = {pvalue_buried_only}")) +
  xlab("surfdist")
) %>% qsave(f("output-distance_to_surface-box-buried_only-zoom{tmpmax}{tmp_plddt}.pdf"))
