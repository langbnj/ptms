blang_init()

# To reset:
# file.remove("tmp-qm.rds")
# # file.remove("tmp-qm1.rds") # no longer used
# file.remove("tmp-qc.rds")
# file.remove("tmp-q.rds")

# q <- Query("SELECT m.acc, m.site, a.relasa10 FROM unimod m JOIN alphasa a ON m.acc=a.acc AND m.site=a.site AND m.species='human' AND m.source='UniProt'")
# q <- Query("SELECT a.acc, a.site, m.ptm, a.relasa10 FROM alphasa a LEFT JOIN unimod m ON m.acc=a.acc AND m.site=a.site AND m.ptm IS NOT NULL WHERE a.species='human' LIMIT 1000000")
# q <- Query("SELECT a.acc, a.site, m.ptm, a.relasa10 FROM alphasa a JOIN uniens ue ON ue.acc=a.acc LEFT JOIN unimod m ON m.acc=a.acc AND m.site=a.site AND m.ptm IS NOT NULL WHERE a.species='human' LIMIT 1000000")

# q <- Query("SELECT a.acc, a.site, m.ptm IS NOT NULL AS ptm, a.relasa10 FROM alphasa a JOIN uniens ue ON ue.acc=a.acc LEFT JOIN unimod m ON m.acc=a.acc AND m.site=a.site WHERE a.species='human' AND ue.canonical=1 LIMIT 10000")

# q <- Query("SELECT a.acc, a.site, m.ptm IS NOT NULL AS ptm, a.relasa10 FROM alphasa a JOIN uniens ue ON ue.acc=a.acc LEFT JOIN unimod m ON m.acc=a.acc AND m.site=a.site WHERE a.species='human' AND ue.canonical=1")

# q <- Query("SELECT a.acc, a.site, 1 AS ptm, a.relasa, a.relasa10, a.dis, a.dis10 FROM alphasa a JOIN uniens ue ON ue.acc=a.acc JOIN unimod         m ON m.acc=a.acc AND m.site=a.site WHERE a.species='human' AND ue.canonical=1 UNION
#             SELECT a.acc, a.site, 0 AS ptm, a.relasa, a.relasa10, a.dis, a.dis10 FROM alphasa a JOIN uniens ue ON ue.acc=a.acc JOIN unimod_control m ON m.acc=a.acc AND m.site=a.site WHERE a.species='human' AND ue.canonical=1")

# qm <- Query("SELECT a.acc, a.site, 1 AS ptm, a.relasa, a.relasa10, a.dis, a.dis10 FROM unimod         m JOIN uniens ue ON ue.acc=m.acc JOIN alphasa a ON m.acc=a.acc AND m.site=a.site WHERE m.species='human' AND ue.canonical=1 LIMIT 1000")
# qc <- Query("SELECT a.acc, a.site, 0 AS ptm, a.relasa, a.relasa10, a.dis, a.dis10 FROM unimod_control m JOIN uniens ue ON ue.acc=m.acc JOIN alphasa a ON m.acc=a.acc AND m.site=a.site WHERE m.species='human' AND ue.canonical=1 LIMIT 1000")

# Note: The uniens join looks unnecessary, but joining with uniens with canonical=1 speeds up the query since it happens prior to the join with alphasa, so unimod gets filtered to a subset before joining with alphasa.

# file.remove("tmp-q.rds")
if (!file.exists("tmp-q.rds")) {
  time0 <- proc.time()
  # file.remove("tmp-qm.rds")
  # file.remove("tmp-qm1.rds")
  if (!file.exists("tmp-qm.rds")) {
    # Run query (~30 mins)
    # with ptm IS NOT NULL (use this one):
    # qm <- Query("SELECT a.acc, a.site, 1 AS ptm, a.relasa, a.relasa10, a.dis, a.dis10 FROM unimod         m JOIN uniens ue ON ue.acc=m.acc JOIN alphasa a ON m.acc=a.acc AND m.site=a.site WHERE m.species='human' AND m.ptm IS NOT NULL AND ue.canonical=1")
    # qm <- Query("SELECT a.acc, a.site, a.aa, m.ptm, m.source, m.subset, m.scale, a.relasa, a.relasa10, a.dis, a.dis10 FROM unimod         m JOIN uniens ue ON ue.acc=m.acc JOIN alphasa a ON m.acc=a.acc AND m.site=a.site WHERE m.species='human' AND m.ptm IS NOT NULL AND ue.canonical=1")
    # qm <- Query("SELECT a.acc, a.site, a.aa, m.ptm, m.source, m.subset, m.scale, a.relasa, a.relasa10, a.dis, a.dis10 FROM unimod m, uniseq s, alphaseq aseq, alphasa a WHERE m.species='human' AND m.acc=s.acc AND s.type='UniProt' AND m.acc=aseq.acc AND s.seq=aseq.seq AND m.acc=a.acc AND m.site=a.site AND a.species='human' AND m.ptm IS NOT NULL;")
    
    # With alphamap, and without intra- or inter-membrane residues (membrane is either 1 or NULL, so we're querying for membrane IS NULL):
    # system.time(qm <- Query("SELECT a.acc, a.site, a.aa, m.ptm, m.source, m.subset, m.scale, a.relasa, a.relasa10, a.dis, a.dis10, a.plddt, a.plddt10, a.membrane FROM unimod m, uniseq s, alphamap am, alphasa a WHERE m.species='human' AND m.ptm IS NOT NULL AND m.acc=s.acc AND s.type='UniProt' AND LENGTH(s.seq)>=16 AND m.acc=am.value AND am.type='uniprot' AND am.version='2022_04' AND am.best=1 AND am.map=a.acc AND am.afdb=a.afdb AND m.site=a.site AND a.membrane IS NULL"))
    # With PAE values:
    system.time(qm <- Query("SELECT 'Modified' AS ptmbin, m.acc, m.site, m.aa, m.ptm, m.source, m.subset, m.scale, a.relasa, a.relasa10, a.dis, a.dis10, a.plddt, a.plddt10, a.membrane, IFNULL(COUNT(DISTINCT c.site1, c.site2), 0) AS contacted_residues, MIN(c.pae) AS min_pae, AVG(c.pae) AS avg_pae, MAX(c.pae) AS max_pae FROM unimod m, uniseq s, alphamap am, alphasa a LEFT OUTER JOIN alphacon c ON a.acc=c.acc AND a.afdb=c.afdb AND (a.site=c.site1 OR a.site=c.site2) WHERE m.species='human' AND m.ptm IS NOT NULL AND m.acc=s.acc AND s.type='UniProt' AND LENGTH(s.seq)>=16 AND m.acc=am.value AND am.type='uniprot' AND am.version='2022_04' AND am.best=1 AND am.map=a.acc AND am.afdb=a.afdb AND m.site=a.site AND a.membrane IS NULL GROUP BY m.acc, m.site, m.aa, m.ptm, m.source, m.subset, m.scale"))
    write_rds(qm, "tmp-qm.rds")

    # # Run query (~30 mins)
    # # includes ptm IS NULL sites (just for diagnostics, not used below - the difference is extremely marginal since there are almost no ptm=NULL sites)
    # # qm1 <- Query("SELECT a.acc, a.site, 1 AS ptm, a.relasa, a.relasa10, a.dis, a.dis10 FROM unimod         m JOIN uniens ue ON ue.acc=m.acc JOIN alphasa a ON m.acc=a.acc AND m.site=a.site WHERE m.species='human' AND ue.canonical=1")
    # qm1 <- Query("SELECT a.acc, a.site, a.aa, m.ptm, a.relasa, a.relasa10, a.dis, a.dis10 FROM unimod         m JOIN uniens ue ON ue.acc=m.acc JOIN alphasa a ON m.acc=a.acc AND m.site=a.site WHERE m.species='human' AND ue.canonical=1")
    # write_rds(qm1, "tmp-qm1.rds")
  } else {
    # Load rds
    qm <- read_rds("tmp-qm.rds")
  }
  proc.time() - time0
  
  # qm1
  # qm
  # # >> The size difference between these two is actually very small. This is because there are actually very few ptm=NULL sites.
  
  time0 <- proc.time()
  # file.remove("tmp-qc.rds")
  if (!file.exists("tmp-qc.rds")) {
    # Run query (~30 mins)
    # With PAE values:
    system.time(qc <- Query("SELECT 'Control' AS ptmbin, m.acc, m.site, m.aa, '0' AS ptm, '0' AS source, '0' AS subset, '0' AS scale, a.relasa, a.relasa10, a.dis, a.dis10, a.plddt, a.plddt10, a.membrane, IFNULL(COUNT(DISTINCT c.site1, c.site2), 0) AS contacted_residues, MIN(c.pae) AS min_pae, AVG(c.pae) AS avg_pae, MAX(c.pae) AS max_pae FROM unimod_control m, uniseq s, alphamap am, alphasa a LEFT OUTER JOIN alphacon c ON a.acc=c.acc AND a.afdb=c.afdb AND (a.site=c.site1 OR a.site=c.site2) WHERE m.species='human' AND m.acc=s.acc AND s.type='UniProt' AND LENGTH(s.seq)>=16 AND m.acc=am.value AND am.type='uniprot' AND am.version='2022_04' AND am.best=1 AND am.map=a.acc AND am.afdb=a.afdb AND m.site=a.site AND a.membrane IS NULL GROUP BY m.acc, m.site, m.aa"))
    write_rds(qc, "tmp-qc.rds")
  } else {
    # Load rds
    qc <- read_rds("tmp-qc.rds")
  }
  proc.time() - time0
  
  qm
  qc
  
  # tmp qm:
  # Add aa column retroactively to save time (from first character of ptm column)
  # qc %<>% mutate(source = '0', subset = '0', scale = '0', .after = ptm)
  # qc %<>% mutate(source = '0', subset = '0', scale = '0', .after = ptm)
  # qm %>% mutate(aa = str_sub(ptm, 1, 1), .after = site) %>% mutate(source = '0', subset = '0', scale = '0', .after = ptm)
  # qm %<>% mutate(aa = str_sub(ptm, 1, 1), .after = site) %>% mutate(source = '0', subset = '0', scale = '0', .after = ptm)
  
  # acc-level filtering (not used)
  # # Filter out accs from the control that don't have PTMs
  # qm_accs <- qm %>% select(acc) %>% unique
  # qm_accs %>% nrow
  # qc_accs <- qc %>% select(acc) %>% unique
  # qc_accs %>% nrow
  # # How many mod accs are in control?
  # qm_accs %>% mutate(in_qc = acc %in% qc_accs$acc) %>% group_by(in_qc) %>% tally
  # # >> All of them (good)
  # # How many control accs are in mod?
  # qc_accs %>% mutate(in_qm = acc %in% qm_accs$acc) %>% group_by(in_qm) %>% tally
  # # >> 114 are not
  # # Filtering these out from qc:
  # qc %<>% filter(acc %in% qm_accs$acc)
  # qc %>% select(acc) %>% unique %>% nrow
  
  # ptm|acc-level filtering (used)
  # Limit to ptm|accs that have both PTM and control residues
  # i.e. skip any ptm|acc combinations that do not have control residues
  qc %>% select(acc, site) %>% unique %>% nrow
  qm %>% select(acc, site) %>% unique %>% nrow
  qc %>% select(ptm, acc) %>% unique %>% nrow
  qm %>% select(ptm, acc) %>% unique %>% nrow
  # How many mod aa|accs are in control?
  # Intersection:
  qm %>% inner_join(qc %>% select(aa, acc) %>% distinct, by = c("aa", "acc")) %>% select(aa, acc) %>% distinct %>% arrange(aa, acc) %>% nrow
  qc %>% inner_join(qm %>% select(aa, acc) %>% distinct, by = c("aa", "acc")) %>% select(aa, acc) %>% distinct %>% arrange(aa, acc) %>% nrow
  identical(qm %>% inner_join(qc %>% select(aa, acc) %>% distinct, by = c("aa", "acc")) %>% select(aa, acc) %>% distinct %>% arrange(aa, acc) %>% nrow,
            qc %>% inner_join(qm %>% select(aa, acc) %>% distinct, by = c("aa", "acc")) %>% select(aa, acc) %>% distinct %>% arrange(aa, acc) %>% nrow)
  # >> The two directions are of course identical (doesn't matter if we start with qm or qc for the intersection)

  # mod aa|accs:
  qm %>% left_join(qc %>% select(aa, acc) %>% distinct, by = c("aa", "acc")) %>% select(aa, acc) %>% distinct %>% arrange(aa, acc) %>% nrow
  # mod aa|accs that have controls:
  qm %>% inner_join(qc %>% select(aa, acc) %>% distinct, by = c("aa", "acc")) %>% select(aa, acc) %>% distinct %>% arrange(aa, acc) %>% nrow
  # mod aa|accs that are missing controls:
  qm %>% left_join(qc %>% select(aa, acc) %>% distinct, by = c("aa", "acc")) %>% select(aa, acc) %>% distinct %>% arrange(aa, acc) %>% nrow -
    qm %>% inner_join(qc %>% select(aa, acc) %>% distinct, by = c("aa", "acc")) %>% select(aa, acc) %>% distinct %>% arrange(aa, acc) %>% nrow
  
  # control aa|accs:
  qc %>% left_join(qm %>% select(aa, acc) %>% distinct, by = c("aa", "acc")) %>% select(aa, acc) %>% distinct %>% arrange(aa, acc) %>% nrow
  # control aa|accs that have mods:
  qc %>% inner_join(qm %>% select(aa, acc) %>% distinct, by = c("aa", "acc")) %>% select(aa, acc) %>% distinct %>% arrange(aa, acc) %>% nrow
  # control aa|accs that are missing mods:
  qc %>% left_join(qm %>% select(aa, acc) %>% distinct, by = c("aa", "acc")) %>% select(aa, acc) %>% distinct %>% arrange(aa, acc) %>% nrow -
  qc %>% inner_join(qm %>% select(aa, acc) %>% distinct, by = c("aa", "acc")) %>% select(aa, acc) %>% distinct %>% arrange(aa, acc) %>% nrow
  
  # # Make backups
  # qm0 <- qm
  # qc0 <- qc

  # Filter qm
  qm %<>% inner_join(qc %>% select(aa, acc) %>% distinct, by = c("aa", "acc"))
  # Filter qc
  qc %<>% inner_join(qm %>% select(aa, acc) %>% distinct, by = c("aa", "acc"))
  
  # mod aa|accs:
  qm %>% left_join(qc %>% select(aa, acc) %>% distinct, by = c("aa", "acc")) %>% select(aa, acc) %>% distinct %>% arrange(aa, acc) %>% nrow
  # mod aa|accs that have controls:
  qm %>% inner_join(qc %>% select(aa, acc) %>% distinct, by = c("aa", "acc")) %>% select(aa, acc) %>% distinct %>% arrange(aa, acc) %>% nrow
  # mod aa|accs that are missing controls:
  qm %>% left_join(qc %>% select(aa, acc) %>% distinct, by = c("aa", "acc")) %>% select(aa, acc) %>% distinct %>% arrange(aa, acc) %>% nrow -
  qm %>% inner_join(qc %>% select(aa, acc) %>% distinct, by = c("aa", "acc")) %>% select(aa, acc) %>% distinct %>% arrange(aa, acc) %>% nrow
  
  # control aa|accs:
  qc %>% left_join(qm %>% select(aa, acc) %>% distinct, by = c("aa", "acc")) %>% select(aa, acc) %>% distinct %>% arrange(aa, acc) %>% nrow
  # control aa|accs that have mods:
  qc %>% inner_join(qm %>% select(aa, acc) %>% distinct, by = c("aa", "acc")) %>% select(aa, acc) %>% distinct %>% arrange(aa, acc) %>% nrow
  # control aa|accs that are missing mods:
  qc %>% left_join(qm %>% select(aa, acc) %>% distinct, by = c("aa", "acc")) %>% select(aa, acc) %>% distinct %>% arrange(aa, acc) %>% nrow -
  qc %>% inner_join(qm %>% select(aa, acc) %>% distinct, by = c("aa", "acc")) %>% select(aa, acc) %>% distinct %>% arrange(aa, acc) %>% nrow
  
  qc %>% select(aa, acc) %>% unique %>% nrow
  qm %>% select(aa, acc) %>% unique %>% nrow
  # qc %>% select(aa, acc) %>% distinct %>% arrange(aa, acc)
  # qm %>% select(aa, acc) %>% distinct %>% arrange(aa, acc)
  identical(qc %>% select(aa, acc) %>% distinct %>% arrange(aa, acc), qm %>% select(aa, acc) %>% distinct %>% arrange(aa, acc))
  qm %>% select(acc) %>% unique %>% nrow
  qc %>% select(acc) %>% unique %>% nrow
  identical(qc %>% select(acc) %>% distinct, qm %>% select(acc) %>% distinct)
  
  # Make a binary factor PTM column
  qm %<>% mutate(ptmbin = "Modified")
  qc %<>% mutate(ptmbin = "Control")

  # Combine qm and qc to form q
  q <- bind_rows(qm, qc)
  q %<>% mutate(ptmbin = as.factor(ptmbin))
  rm(qm, qc)

  # # Unique q
  # Unnecessary
  # q <- unique(q)
  
  write_rds(q, "tmp-q.rds")
} else {
  # qm <- read_rds("tmp-qm.rds")
  # # qm1 <- read_rds("tmp-qm1.rds")
  # qc <- read_rds("tmp-qm.rds")
  q <- read_rds("tmp-q.rds")
}
q

# Diagnostics
q %>% select(acc) %>% unique
# >> 19,236 accs
# Check if there are any isoforms (containing "-")
q %>% filter(str_detect(acc, "-")) %>% select(acc) %>% unique
# >> no isoforms (good!)
q %>% group_by(ptm) %>% tally %>% arrange(desc(n)) %>% print(n=50)
q %>% filter(ptmbin == 'Modified') %>% select(acc) %>% unique
# >> 19,236 modified accs
q %>% filter(ptmbin == 'Control') %>% select(acc) %>% unique
# >> 19,236 control accs

# q %>% nrow
# q %>% unique %>% nrow
# q %>% select(acc, site, ptm) %>% unique %>% nrow
# >> OK, acc|site|ptm is unique!

# # Remove accs where the uniseq and alphaseq sequences don't match
# # First, get all accs with matching sequences
# qa <- Query("SELECT DISTINCT s.acc FROM uniseq s, alphaseq a WHERE s.type='UniProt' AND s.species='human' AND s.acc=a.acc AND s.seq=a.seq ORDER BY s.acc")
# qa %>% nrow
# q
# q %<>% filter(acc %in% qa$acc)
# q
# # >> This doesn't actually remove anything anymore (the above queries already verify sequences between uniseq and alphaseq):
# q %>% filter(ptm != "0") %>% select(acc) %>% unique
# # >> 19,236 modified accs
# q %>% filter(ptm == "0") %>% select(acc) %>% unique
# # >> 19,236 control accs

# Make a binary factor PTM column
# # q %<>% mutate(ptmbin = if_else(ptm == 0, "No", "Yes"))
# q %<>% mutate(ptmbin = if_else(ptm == 0, "Control", "Modified"))
# q %<>% mutate(ptmbin = as.factor(ptmbin))

# Set relasa and relasa10 values larger than 1 to 1 (since a residue can't be more than 100% exposed)
q[q$relasa > 1,]$relasa <- 1
q[q$relasa10 > 1,]$relasa10 <- 1

# Keep only PTMs that have ≥1000 sites
# q %>% filter(ptmbin == "Modified") %>% group_by(ptm) %>% summarise(sites = n_distinct(acc, site)) %>% arrange(desc(sites)) %>% filter(sites >= 1000)
# q %>% filter(ptmbin == "Modified") %>% group_by(ptm) %>% summarise(sites = n_distinct(acc, site)) %>% arrange(desc(sites)) %>% filter(sites >= 1000) %>% pull(ptm) %>% sort
q %>% filter(ptmbin == "Modified") %>% group_by(ptm) %>% tally %>% arrange(desc(n))
q %>% select(acc, site, ptm, ptmbin) %>% unique %>% filter(ptmbin == "Modified") %>% group_by(ptm) %>% tally %>% arrange(desc(n))
# ptm1000 <- q %>% select(acc, site, ptm, ptmbin) %>% unique %>% filter(ptmbin == "Modified") %>% group_by(ptm) %>% summarise(sites = n_distinct(acc, site)) %>% arrange(desc(sites)) %>% filter(sites >= 1000) %>% pull(ptm) %>% sort
ptm1000 <- q %>% select(acc, site, ptm, ptmbin) %>% unique %>% filter(ptmbin == "Modified") %>% group_by(ptm) %>% summarise(sites = n_distinct(acc, site)) %>% arrange(desc(sites)) %>% filter(sites >= 950) %>% pull(ptm) %>% sort
ptm1000
length(ptm1000)
q
q %<>% filter(ptm %in% c(ptm1000, 0))
q %>% filter(ptmbin == "Modified") %>% group_by(ptm) %>% tally %>% arrange(desc(n))
q
# Get residues modified by these PTM types
q %>% filter(ptmbin == "Modified") %>% group_by(aa) %>% tally %>% arrange(desc(n))
aa1000 <- q %>% select(acc, site, aa, ptmbin) %>% unique %>% filter(ptmbin == "Modified") %>% group_by(aa) %>% tally %>% arrange(desc(n)) %>% pull(aa) %>% sort
aa1000
aa1000 %>% length
# Filter to only keep these residues
q %>% filter(ptmbin == "Control") %>% group_by(aa) %>% tally
q
q %<>% filter(aa %in% aa1000)
q

# Calculate aa weight to be applied to control residues (so their AA frequencies match those of the PTMs)
q %>% filter(ptmbin == "Modified") %>% group_by(aa) %>% tally %>% mutate(freq = n / sum(n)) %>% arrange(desc(n))
q %>% filter(ptmbin == "Control") %>% group_by(aa) %>% tally %>% mutate(freq = n / sum(n)) %>% arrange(desc(n))
q %>% filter(ptmbin == "Modified") %>% group_by(aa) %>% tally %>% mutate(freq = n / sum(n)) %>% arrange(desc(n)) %>% 
  left_join(q %>% filter(ptmbin == "Control") %>% group_by(aa) %>% tally %>% mutate(freq = n / sum(n)) %>% arrange(desc(n)), by = join_by(aa), suffix = c("_ptm", "_control")) %>%
  mutate(enrich = freq_ptm / freq_control, log2enrich = log2(enrich))
q %>% filter(ptmbin == "Modified") %>% group_by(aa, buried = relasa == 0) %>% tally %>% mutate(freq = n / sum(n)) %>% arrange(desc(buried), desc(n)) %>% 
  left_join(q %>% filter(ptmbin == "Control") %>% group_by(aa, buried = relasa == 0) %>% tally %>% mutate(freq = n / sum(n)) %>% arrange(desc(buried), desc(n)), by = join_by(aa, buried), suffix = c("_ptm", "_control")) %>%
  mutate(enrich = freq_ptm / freq_control, log2enrich = log2(enrich))
q %>% filter(ptmbin == "Modified") %>% group_by(aa, surface = relasa == 1) %>% tally %>% mutate(freq = n / sum(n)) %>% arrange(desc(surface), desc(n)) %>% 
  left_join(q %>% filter(ptmbin == "Control") %>% group_by(aa, surface = relasa == 1) %>% tally %>% mutate(freq = n / sum(n)) %>% arrange(desc(surface), desc(n)), by = join_by(aa, surface), suffix = c("_ptm", "_control")) %>%
  mutate(enrich = freq_ptm / freq_control, log2enrich = log2(enrich))
q %>% filter(ptmbin == "Modified" & aa == "S" & relasa == 0) %>% nrow
q %>% filter(ptmbin == "Control" & aa == "S" & relasa == 0) %>% nrow
q %>% filter(ptmbin == "Modified" & aa == "S" & relasa == 1) %>% nrow
q %>% filter(ptmbin == "Control" & aa == "S" & relasa == 1) %>% nrow
q %>% filter(ptmbin == "Modified") %>% group_by(aa) %>% tally %>% arrange(desc(n))
q %>% filter(ptmbin == "Modified") %>% group_by(aa) %>% tally %>% mutate(freq = n / sum(n)) %>% arrange(desc(n))
q %>% select(acc, site, aa, ptmbin) %>% unique %>% filter(ptmbin == "Modified") %>% group_by(aa) %>% tally %>% mutate(freq = n / sum(n)) %>% arrange(desc(n))
# Before uniquing acc|site (multiple occurrences from multiple sources)
# Add "freq" column (for weighting of control residues to the same proportion as PTM sites, by AA type)
# aaweight <- q %>% filter(ptmbin == "Modified") %>% group_by(aa) %>% tally %>% mutate(freq = n / sum(n)) %>% arrange(desc(n)) %>% select(aa, freq)
aaweight <- q %>% select(acc, site, aa, ptmbin) %>% unique %>% filter(ptmbin == "Modified") %>% group_by(aa) %>% tally %>% mutate(freq = n / sum(n)) %>% arrange(desc(n)) %>% select(aa, freq)
aaweight
# Add weights to q
# q$freq <- NULL
# q
identical(q %>% left_join(aaweight, by = "aa"), q %>% left_join(aaweight, by = join_by(aa)))
q
q %<>% left_join(aaweight, join_by(aa))
# Set PTM site weight to 1 (only Control should be weighted)
q %<>% mutate(freq = ifelse(ptmbin == "Modified", 1, freq))
q
# q %>% filter(ptmbin == "Control")
# q %>% filter(ptmbin == "Control" & aa == "S")

# PTM types with ≥100 buried sites:
# By fracburied
q %>% group_by(ptm, relasa <= 0.25) %>% tally %>% pivot_wider(names_from = `relasa <= 0.25`, values_from=n) %>% transmute(nburied = `TRUE`, fracburied = `TRUE` / (`FALSE` + `TRUE`)) %>% arrange(desc(fracburied)) %>% filter(nburied >= 100) %>% filter(ptm != '0')
# By nburied
q %>% group_by(ptm, relasa <= 0.25) %>% tally %>% pivot_wider(names_from = `relasa <= 0.25`, values_from=n) %>% transmute(nburied = `TRUE`, fracburied = `TRUE` / (`FALSE` + `TRUE`)) %>% arrange(desc(nburied)) %>% filter(nburied >= 100) %>% filter(ptm != '0')

# Set "dis" column using dis/str call from dis10 and core/surf call from relasa
q %<>% mutate(dis = paste0(ifelse(dis10 == "*", "dis", "str"), ifelse(relasa <= 0.25, "core", "surf")))
# Replace discore
q %<>% mutate(dis = replace(dis, dis == "discore", "strcore"))

# Calculate weight not just by aa, but also by dissurf classification
# q %>% select(acc, site, dis, aa, ptmbin) %>% unique
# q %>% select(acc, site, dis, aa, ptmbin) %>% unique %>% filter(ptmbin == "Modified") %>% group_by(dis, aa) %>% tally
# q %>% select(acc, site, dis, aa, ptmbin) %>% unique %>% filter(ptmbin == "Modified") %>% group_by(dis, aa) %>% tally %>% mutate(freq = n / sum(n)) %>% arrange(desc(n)) %>% select(dis, aa, freq)
dissurf_aaweight <- q %>% select(acc, site, dis, aa, ptmbin) %>% unique %>% filter(ptmbin == "Modified") %>% group_by(dis, aa) %>% tally %>% mutate(dissurf_aafreq = n / sum(n)) %>% arrange(desc(n)) %>% select(dis, aa, dissurf_aafreq)
dissurf_aaweight
# Add to q
q %<>% left_join(dissurf_aaweight, join_by(dis, aa))
# Set PTM site weight to 1 (only Control should be weighted)
q %<>% mutate(freq = ifelse(ptmbin == "Modified", 1, dissurf_aafreq))

# Load evorates from ../evolutionary_rate_finding_best_stratification/tmp-dataframe-….txt
source <- "all"
evorate <- "rate4site_einsi_tree_1para"
qevo <- tibble(read_tsv(f("../evolutionary_rate_finding_best_stratification/tmp-dataframe-all-human-{source}-{evorate}-AlphaFold-coresurf.txt")))
qevo_plddt70 <- tibble(read_tsv(f("../evolutionary_rate_finding_best_stratification/tmp-dataframe-all-human-{source}-{evorate}-AlphaFold_pLDDT70-coresurf.txt")))
# There should be 2,068,762 residues (acc|sites) in here: cat tmp-dataframe-all-human-all-rate4site_einsi_tree_1para-AlphaFold-coresurf.txt | tail -n+2 | cut -f2,3 | sort | uniq | wc -l
# qevo %<>% filter(dis != "A") # This removes ~200,000 residues for which there doesn't seem to be an AlphaFold structure
qevo %<>% mutate(dis = replace(dis, dis == "discore", "strcore"))
qevo$type <- factor(qevo$type, levels=c("P", "C"))
qevo$dis <- factor(qevo$dis, levels=c("strcore", "strsurf", "dissurf"))
# qevo %>% mutate(dis = fct_recode(dis, "Buried"="strcore", "Structured"="strsurf", "Disordered"="dissurf"))
levels(qevo$dis) <- list("Buried"="strcore", "Structured"="strsurf", "Disordered"="dissurf")
qevo
# qevo_plddt70
qevo_plddt70 %<>% mutate(dis = replace(dis, dis == "discore", "strcore"))
qevo_plddt70$type <- factor(qevo_plddt70$type, levels=c("P", "C"))
qevo_plddt70$dis <- factor(qevo_plddt70$dis, levels=c("strcore", "strsurf", "dissurf"))
# qevo_plddt70 %>% mutate(dis = fct_recode(dis, "Buried"="strcore", "Structured"="strsurf", "Disordered"="dissurf"))
levels(qevo_plddt70$dis) <- list("Buried"="strcore", "Structured"="strsurf", "Disordered"="dissurf")
qevo_plddt70

# Collapse to per-residue level (qevo currently still has a PTM column)
qevo %>% select(acc, site) %>% unique %>% nrow
qevo %>% select(acc, site, rate) %>% unique %>% nrow
qevo %>% select(acc, site, type, rate) %>% unique %>% nrow
# >> All 2,068,762 >> OK, each residue has a single evorate value (as expected) (even when including type=C/PTM)
qevo %>% select(acc, site, rate) %>% unique -> qevo_per_residue
qevo_per_residue %>% nrow
# >> We have evorates for 2,068,762 residues (acc|sites).
q %<>% left_join(qevo_per_residue, join_by(acc, site))
q
summary(q$rate)
qevo_per_residue %>% select(acc, site) %>% unique %>% nrow
# >> 2,068,762
q %>% filter(!is.na(rate)) %>% select(acc, site) %>% unique %>% nrow
# >> 2,068,762
# >> All qevo sites are in q, so I should be able to make an equivalent plot to finding_best_stratification.R in here!







# Plots

# relasa10 (smoothed)
q %>%
  ggplot(aes(x = relasa10, colour = ptmbin, fill = ptmbin)) +
  geom_density(alpha = 0.3) +
  # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
  # scale_colour_manual(values = c("Yes" = myo, "No" = mybd), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
  scale_colour_manual(values = c("Yes" = ptmb, "No" = ptmg), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
  # scale_x_continuous(limits = c(0, 1)) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_minimal() +
  xlab("AlphaFold relative accessible surface area (smoothed, ±10 aa)") +
  ylab("Density")
ggsave("output-density-canonical-smoothed.pdf", width = 6, height = 4)

# relasa (unsmoothed)
q
identical(q %>% filter(ptm == 0), q %>% filter(ptmbin == "Control"))
q %>% filter(ptm == 0)
q %>% filter(ptm == 0) %>% summary
q %>% filter(ptm == 0) %>% group_by(ptm) %>% tally
q %>% filter(ptm == 0) %>% group_by(acc) %>% tally %>% nrow
# >> 19,236 accs now (after making sure uniseq sequence matches alphaseq)

q %>% filter(ptm == 0) %>% select(acc) %>% unique %>% nrow
# >> 19,236 accs
q %>% filter(ptm == 0) %>% select(acc, site) %>% nrow
# >> 3,294,049 total control acc|sites
q %>% filter(ptm == 0) %>% select(acc, site) %>% unique %>% nrow
# >> 3,294,049 total unique control acc|sites
# Buried sites:
# q %>% filter(ptmbin == "Modified" & relasa <= 0.25) %>% select(acc, site) %>% unique %>% nrow
# q %>% filter(ptmbin == "Modified" & relasa <= 0.25) %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_total_controls <- q %>% filter(ptmbin == "Control") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_buried_controls <- q %>% filter(ptmbin == "Control" & dis == "strcore") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_fully_buried_controls <- q %>% filter(ptmbin == "Control" & dis == "strcore" & relasa == 0) %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_fully_buried_sites <- q %>% filter(ptmbin == "Modified" & dis == "strcore" & relasa == 0) %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_buried_sites <- q %>% filter(ptmbin == "Modified" & dis == "strcore") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_structured_sites <- q %>% filter(ptmbin == "Modified" & dis == "strsurf") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_disordered_sites <- q %>% filter(ptmbin == "Modified" & dis == "dissurf") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_total_sites <- q %>% filter(ptmbin == "Modified") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_buried_controls / my_total_controls
my_buried_controls_fraction <- my_buried_controls / my_total_controls
my_buried_controls_fraction
round(my_buried_controls_fraction, 3)
# >> 26.6% of Control residues are buried
# >> EXPORT
my_total_controls
my_fully_buried_controls
my_fully_buried_sites
# >> 12,509 sites are fully buried (SASA = 0) (!)
# >> EXPORT
q %>% filter(ptmbin == "Modified" & dis == "strcore" & relasa == 0) %>% select(ptm, acc, site) %>% unique
# >> 12,717 when counting different PTM types separately (but acc|sites is the less confusing number to report)
# Breakdown by PTM type:
q %>% filter(ptmbin == "Modified" & dis == "strcore" & relasa == 0) %>% select(ptm, acc, site) %>% unique %>% group_by(ptm) %>% tally %>% arrange(n)
q %>% filter(ptmbin == "Modified" & dis == "strcore" & relasa == 0) %>% select(ptm, acc, site) %>% unique %>% filter(ptm %in% c("S-p", "T-p", "Y-p")) %>% select(acc, site) %>% unique %>% nrow
# >> 11,046
# >> EXPORT
my_buried_sites
my_structured_sites
my_disordered_sites
my_fully_buried_control_fraction <- my_fully_buried_controls / my_total_controls
my_fully_buried_fraction <- my_fully_buried_sites / my_total_sites
my_buried_fraction <- my_buried_sites / my_total_sites
my_structured_fraction <- my_structured_sites / my_total_sites
my_disordered_fraction <- my_disordered_sites / my_total_sites
my_fully_buried_control_fraction
round(my_fully_buried_control_fraction, 3)
# >> 3.3% of control residues are fully buried (!)
# >> EXPORT
my_fully_buried_fraction
# >> 2% of PTM sites are fully buried (!) That's a really substantial number compared to control!
# >> EXPORT
my_buried_fraction
my_structured_fraction
my_disordered_fraction
# Format as percentage with one decimal
my_buried_fraction <- scales::percent(my_buried_fraction, accuracy = 0.1)
my_structured_fraction <- scales::percent(my_structured_fraction, accuracy = 0.1)
my_disordered_fraction <- scales::percent(my_disordered_fraction, accuracy = 0.1)
# q %>% filter(ptm != 0 & acc == "P04637") %>% arrange(acc, site, aa, ptm, source, subset)
(q %>%
  select(acc, site, relasa, ptmbin, freq) %>% unique %>%
  # ggplot(aes(x = relasa, colour = ptmbin, fill = ptmbin)) +
  ggplot(aes(x = relasa, colour = ptmbin, fill = ptmbin, weight = freq)) +
    # geom_density(alpha = 0.3) +
    geom_density(alpha = 0.3, bounds = c(0, 1)) +
    # geom_density(alpha = 0.2) +
    # geom_density(alpha = 0.15) +
    # geom_density(alpha = 0.3, bw= 0.3) +
  # geom_vline(xintercept = 0.25, linetype = "dashed") +
  # geom_vline(xintercept = 0.55, linetype = "dashed") +
  # annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 2.25), linetype = "dotted", linewidth = 0.5/.weight) +
  # annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 2.25), linetype = "dotted", linewidth = 0.5/.weight) +
  # Before using bounds:
  # annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 2.5), linetype = "dashed", linewidth = 0.5/.weight) +
  # annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 2.5), linetype = "dashed", linewidth = 0.5/.weight) +
  # annotate(geom = "text", x = 0.125, y = 0.15, label = "Buried", size = 5/.pt) +
  # annotate(geom = "text", x = 0.40, y = 0.15, label = "Structured", size = 5/.pt) +
  # annotate(geom = "text", x = 0.775, y = 0.15, label = "Disordered", size = 5/.pt) +
  # annotate(geom = "text", x = 0.125, y = 0.375, label = my_buried_fraction, size = 5/.pt) +
  # annotate(geom = "text", x = 0.40, y = 0.375, label = my_structured_fraction, size = 5/.pt) +
  # annotate(geom = "text", x = 0.775, y = 0.375, label = my_disordered_fraction, size = 5/.pt) +
  # Using bounds:
  annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 3.25), linetype = "dashed", linewidth = 0.5/.weight) +
  annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 3.25), linetype = "dashed", linewidth = 0.5/.weight) +
  annotate(geom = "text", x = 0.125, y = 2.95, label = "Buried", size = 5/.pt) +
  annotate(geom = "text", x = 0.40, y = 2.95, label = "Structured", size = 5/.pt) +
  annotate(geom = "text", x = 0.775, y = 2.95, label = "Disordered", size = 5/.pt) +
  annotate(geom = "text", x = 0.125, y = 2.6, label = my_buried_fraction, size = 5/.pt) +
  annotate(geom = "text", x = 0.40, y = 2.6, label = my_structured_fraction, size = 5/.pt) +
  annotate(geom = "text", x = 0.775, y = 2.6, label = my_disordered_fraction, size = 5/.pt) +
  # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
  # scale_colour_manual(values = c("Yes" = myo, "No" = mybd), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
  # scale_colour_manual(values = c("Yes" = ptmb, "No" = ptmg), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
  # scale_colour_manual(values = c("Modified" = ptmb, "Control" = ptmg), aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
  # scale_colour_manual(values = c("Modified" = scales::viridis_pal()(9)[8], "Control" = scales::viridis_pal()(9)[4]), aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
  # scale_colour_manual(values = c("Modified" = scales::viridis_pal()(5)[4], "Control" = scales::viridis_pal()(5)[2]), aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
  scale_colour_manual(values = ptmcol, aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
  # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = pretty_breaks(2)) +
  # scale_x_continuous(limits = c(0, 1)) +
  # coord_cartesian(xlim = c(0, 1)) +
  # coord_cartesian(xlim = c(0, 1), ylim = c(0, 3)) +
  theme_minimal() +
  # ggtitle("") +
  labs(tag = "a") +
  xlab("Relative accessible surface area") +
  ylab("Probability density") +
  # theme_nature()
  theme_nature(legend_position = "top") +
  theme(axis.title.y = element_text(hjust = 0))
  # theme_nature(legend_position = "top", extra_margin_right = 2, legend_nudge_right = 1.5)
) %>% qsave("output-density-canonical-weighted.pdf", height = 30)

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

# Filtered for pLDDT ≥ 70
# (q %>%
#     filter(plddt >= 70) %>%
#     select(acc, site, relasa, ptmbin, freq) %>% unique %>%
#     ggplot(aes(x = relasa, colour = ptmbin, fill = ptmbin, weight = freq)) +
#     geom_density(alpha = 0.3) +
#     annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 2.5), linetype = "dashed", linewidth = 0.5/.weight) +
#     annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 2.5), linetype = "dashed", linewidth = 0.5/.weight) +
#     annotate(geom = "text", x = 0.125, y = 0.15, label = "Buried", size = 5/.pt) +
#     annotate(geom = "text", x = 0.40, y = 0.15, label = "Structured", size = 5/.pt) +
#     annotate(geom = "text", x = 0.775, y = 0.15, label = "Disordered", size = 5/.pt) +
#     scale_colour_manual(values = ptmcol, aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
#     scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
#     scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = pretty_breaks(2)) +
#     theme_minimal() +
#     labs(tag = "a") +
#     xlab("RSA (pLDDT ≥ 70, confident)") +
#     ylab("Probability density") +
#     theme_nature(legend_position = "top") +
#     theme(axis.title.y = element_text(hjust = 0))
# ) %>% qsave("output-density-canonical-weighted-plddt70.pdf", height = 30)
my_plddt70_total_controls <- q %>% filter(plddt >= 70 & ptmbin == "Control") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt70_fully_buried_controls <- q %>% filter(plddt >= 70 & ptmbin == "Control" & dis == "strcore" & relasa == 0) %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt70_fully_buried_sites <- q %>% filter(plddt >= 70 & ptmbin == "Modified" & dis == "strcore" & relasa == 0) %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt70_buried_sites <- q %>% filter(plddt >= 70 & ptmbin == "Modified" & dis == "strcore") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt70_structured_sites <- q %>% filter(plddt >= 70 & ptmbin == "Modified" & dis == "strsurf") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt70_disordered_sites <- q %>% filter(plddt >= 70 & ptmbin == "Modified" & dis == "dissurf") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt70_total_sites <- q %>% filter(plddt >= 70 & ptmbin == "Modified") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt70_fully_buried_control_fraction <- my_plddt70_fully_buried_controls / my_plddt70_total_controls
my_plddt70_fully_buried_fraction <- my_plddt70_fully_buried_sites / my_plddt70_total_sites
my_plddt70_buried_fraction <- my_plddt70_buried_sites / my_plddt70_total_sites
my_plddt70_structured_fraction <- my_plddt70_structured_sites / my_plddt70_total_sites
my_plddt70_disordered_fraction <- my_plddt70_disordered_sites / my_plddt70_total_sites
my_plddt70_fully_buried_control_fraction
my_plddt70_fully_buried_fraction
my_plddt70_buried_fraction
my_plddt70_structured_fraction
my_plddt70_disordered_fraction
my_plddt70_buried_fraction <- scales::percent(my_plddt70_buried_fraction, accuracy = 0.1)
my_plddt70_structured_fraction <- scales::percent(my_plddt70_structured_fraction, accuracy = 0.1)
my_plddt70_disordered_fraction <- scales::percent(my_plddt70_disordered_fraction, accuracy = 0.1)
(q %>%
    filter(plddt >= 70) %>%
    select(acc, site, relasa, ptmbin, freq) %>% unique %>%
    ggplot(aes(x = relasa, colour = ptmbin, fill = ptmbin, weight = freq)) +
    geom_density(alpha = 0.3, bounds = c(0, 1)) +
    annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 6), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 6), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "text", x = 0.125, y = 5.5, label = "Buried", size = 5/.pt) +
    annotate(geom = "text", x = 0.40, y = 5.5, label = "Structured", size = 5/.pt) +
    annotate(geom = "text", x = 0.775, y = 5.5, label = "Disordered", size = 5/.pt) +
    annotate(geom = "text", x = 0.125, y = 4.75, label = my_plddt70_buried_fraction, size = 5/.pt) +
    annotate(geom = "text", x = 0.40, y = 4.75, label = my_plddt70_structured_fraction, size = 5/.pt) +
    annotate(geom = "text", x = 0.775, y = 4.75, label = my_plddt70_disordered_fraction, size = 5/.pt) +
    scale_colour_manual(values = ptmcol, aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = pretty_breaks(2)) +
    theme_minimal() +
    labs(tag = "a") +
    xlab("RSA (pLDDT ≥ 70, confident)") +
    ylab("Probability density") +
    theme_nature(legend_position = "top") +
    theme(axis.title.y = element_text(hjust = 0))
) %>% qsave("output-density-canonical-weighted-plddt70.pdf", height = 30)

# Filtered for pLDDT ≥ 70 & PAE ≤ 2
my_plddt70_pae2_total_controls <- q %>% filter(plddt >= 70 & min_pae <= 2 & ptmbin == "Control") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt70_pae2_fully_buried_controls <- q %>% filter(plddt >= 70 & min_pae <= 2 & ptmbin == "Control" & dis == "strcore" & relasa == 0) %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt70_pae2_fully_buried_sites <- q %>% filter(plddt >= 70 & min_pae <= 2 & ptmbin == "Modified" & dis == "strcore" & relasa == 0) %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt70_pae2_buried_sites <- q %>% filter(plddt >= 70 & min_pae <= 2 & ptmbin == "Modified" & dis == "strcore") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt70_pae2_structured_sites <- q %>% filter(plddt >= 70 & min_pae <= 2 & ptmbin == "Modified" & dis == "strsurf") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt70_pae2_disordered_sites <- q %>% filter(plddt >= 70 & min_pae <= 2 & ptmbin == "Modified" & dis == "dissurf") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt70_pae2_total_sites <- q %>% filter(plddt >= 70 & min_pae <= 2 & ptmbin == "Modified") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt70_pae2_fully_buried_control_fraction <- my_plddt70_pae2_fully_buried_controls / my_plddt70_pae2_total_controls
my_plddt70_pae2_fully_buried_fraction <- my_plddt70_pae2_fully_buried_sites / my_plddt70_pae2_total_sites
my_plddt70_pae2_buried_fraction <- my_plddt70_pae2_buried_sites / my_plddt70_pae2_total_sites
my_plddt70_pae2_structured_fraction <- my_plddt70_pae2_structured_sites / my_plddt70_pae2_total_sites
my_plddt70_pae2_disordered_fraction <- my_plddt70_pae2_disordered_sites / my_plddt70_pae2_total_sites
my_plddt70_pae2_fully_buried_control_fraction
my_plddt70_pae2_fully_buried_fraction
my_plddt70_pae2_buried_fraction
my_plddt70_pae2_structured_fraction
my_plddt70_pae2_disordered_fraction
my_plddt70_pae2_buried_fraction <- scales::percent(my_plddt70_pae2_buried_fraction, accuracy = 0.1)
my_plddt70_pae2_structured_fraction <- scales::percent(my_plddt70_pae2_structured_fraction, accuracy = 0.1)
my_plddt70_pae2_disordered_fraction <- scales::percent(my_plddt70_pae2_disordered_fraction, accuracy = 0.1)
(q %>% 
    filter(plddt >= 70 & min_pae <= 2) %>%
    select(acc, site, relasa, ptmbin, freq) %>% unique %>%
    ggplot(aes(x = relasa, colour = ptmbin, fill = ptmbin, weight = freq)) +
    geom_density(alpha = 0.3, bounds = c(0, 1)) +
    annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 6.75), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 6.75), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "text", x = 0.125, y = 5.5, label = "Buried", size = 5/.pt) +
    annotate(geom = "text", x = 0.40, y = 5.5, label = "Structured", size = 5/.pt) +
    annotate(geom = "text", x = 0.775, y = 5.5, label = "Disordered", size = 5/.pt) +
    annotate(geom = "text", x = 0.125, y = 4.75, label = my_plddt70_pae2_buried_fraction, size = 5/.pt) +
    annotate(geom = "text", x = 0.40, y = 4.75, label = my_plddt70_pae2_structured_fraction, size = 5/.pt) +
    annotate(geom = "text", x = 0.775, y = 4.75, label = my_plddt70_pae2_disordered_fraction, size = 5/.pt) +
    scale_colour_manual(values = ptmcol, aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = pretty_breaks(2)) +
    theme_minimal() +
    labs(tag = "b") +
    xlab("RSA (pLDDT ≥ 70, PAE ≤ 2)") +
    ylab("Probability density") +
    theme_nature(legend_position = "top") +
    theme(axis.title.y = element_text(hjust = 0))
) %>% qsave("output-density-canonical-weighted-plddt70_pae2.pdf", height = 30)

# Filtered for pLDDT ≥ 90
# (q %>%
#     filter(plddt >= 90) %>%
#     select(acc, site, relasa, ptmbin, freq) %>% unique %>%
#     ggplot(aes(x = relasa, colour = ptmbin, fill = ptmbin, weight = freq)) +
#     geom_density(alpha = 0.3) +
#     annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 2.5), linetype = "dashed", linewidth = 0.5/.weight) +
#     annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 2.5), linetype = "dashed", linewidth = 0.5/.weight) +
#     annotate(geom = "text", x = 0.125, y = 0.15, label = "Buried", size = 5/.pt) +
#     annotate(geom = "text", x = 0.40, y = 0.15, label = "Structured", size = 5/.pt) +
#     annotate(geom = "text", x = 0.775, y = 0.15, label = "Disordered", size = 5/.pt) +
#     scale_colour_manual(values = ptmcol, aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
#     scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
#     scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = pretty_breaks(2)) +
#     theme_minimal() +
#     labs(tag = "b") +
#     xlab("RSA (pLDDT ≥ 90, highly confident)") +
#     ylab("Probability density") +
#     theme_nature(legend_position = "top") +
#     theme(axis.title.y = element_text(hjust = 0))
# ) %>% qsave("output-density-canonical-weighted-plddt90.pdf", height = 30)
my_plddt90_total_controls <- q %>% filter(plddt >= 90 & ptmbin == "Control") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt90_fully_buried_controls <- q %>% filter(plddt >= 90 & ptmbin == "Control" & dis == "strcore" & relasa == 0) %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt90_fully_buried_sites <- q %>% filter(plddt >= 90 & ptmbin == "Modified" & dis == "strcore" & relasa == 0) %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt90_buried_sites <- q %>% filter(plddt >= 90 & ptmbin == "Modified" & dis == "strcore") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt90_structured_sites <- q %>% filter(plddt >= 90 & ptmbin == "Modified" & dis == "strsurf") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt90_disordered_sites <- q %>% filter(plddt >= 90 & ptmbin == "Modified" & dis == "dissurf") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt90_total_sites <- q %>% filter(plddt >= 90 & ptmbin == "Modified") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt90_fully_buried_control_fraction <- my_plddt90_fully_buried_controls / my_plddt90_total_controls
my_plddt90_fully_buried_fraction <- my_plddt90_fully_buried_sites / my_plddt90_total_sites
my_plddt90_buried_fraction <- my_plddt90_buried_sites / my_plddt90_total_sites
my_plddt90_structured_fraction <- my_plddt90_structured_sites / my_plddt90_total_sites
my_plddt90_disordered_fraction <- my_plddt90_disordered_sites / my_plddt90_total_sites
my_plddt90_fully_buried_control_fraction
my_plddt90_fully_buried_fraction
my_plddt90_buried_fraction
my_plddt90_structured_fraction
my_plddt90_disordered_fraction
my_plddt90_buried_fraction <- scales::percent(my_plddt90_buried_fraction, accuracy = 0.1)
my_plddt90_structured_fraction <- scales::percent(my_plddt90_structured_fraction, accuracy = 0.1)
my_plddt90_disordered_fraction <- scales::percent(my_plddt90_disordered_fraction, accuracy = 0.1)
(q %>%
    filter(plddt >= 90) %>%
    select(acc, site, relasa, ptmbin, freq) %>% unique %>%
    ggplot(aes(x = relasa, colour = ptmbin, fill = ptmbin, weight = freq)) +
    geom_density(alpha = 0.3, bounds = c(0, 1)) +
    annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 8), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 8), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "text", x = 0.125, y = 7.5, label = "Buried", size = 5/.pt) +
    annotate(geom = "text", x = 0.40, y = 7.5, label = "Structured", size = 5/.pt) +
    annotate(geom = "text", x = 0.775, y = 7.5, label = "Disordered", size = 5/.pt) +
    annotate(geom = "text", x = 0.125, y = 6.5, label = my_plddt90_buried_fraction, size = 5/.pt) +
    annotate(geom = "text", x = 0.40, y = 6.5, label = my_plddt90_structured_fraction, size = 5/.pt) +
    annotate(geom = "text", x = 0.775, y = 6.5, label = my_plddt90_disordered_fraction, size = 5/.pt) +
    scale_colour_manual(values = ptmcol, aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = pretty_breaks(2)) +
    theme_minimal() +
    labs(tag = "c") +
    xlab("RSA (pLDDT ≥ 90, highly confident)") +
    ylab("Probability density") +
    theme_nature(legend_position = "top") +
    theme(axis.title.y = element_text(hjust = 0))
) %>% qsave("output-density-canonical-weighted-plddt90.pdf", height = 30)

# Filtered for pLDDT ≥ 90 & PAE ≤ 1
my_plddt90_pae1_total_controls <- q %>% filter(plddt >= 90 & min_pae <= 1 & ptmbin == "Control") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt90_pae1_fully_buried_controls <- q %>% filter(plddt >= 90 & min_pae <= 1 & ptmbin == "Control" & dis == "strcore" & relasa == 0) %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt90_pae1_fully_buried_sites <- q %>% filter(plddt >= 90 & min_pae <= 1 & ptmbin == "Modified" & dis == "strcore" & relasa == 0) %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt90_pae1_buried_sites <- q %>% filter(plddt >= 90 & min_pae <= 1 & ptmbin == "Modified" & dis == "strcore") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt90_pae1_structured_sites <- q %>% filter(plddt >= 90 & min_pae <= 1 & ptmbin == "Modified" & dis == "strsurf") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt90_pae1_disordered_sites <- q %>% filter(plddt >= 90 & min_pae <= 1 & ptmbin == "Modified" & dis == "dissurf") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt90_pae1_total_sites <- q %>% filter(plddt >= 90 & min_pae <= 1 & ptmbin == "Modified") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_plddt90_pae1_fully_buried_control_fraction <- my_plddt90_pae1_fully_buried_controls / my_plddt90_pae1_total_controls
my_plddt90_pae1_fully_buried_fraction <- my_plddt90_pae1_fully_buried_sites / my_plddt90_pae1_total_sites
my_plddt90_pae1_buried_fraction <- my_plddt90_pae1_buried_sites / my_plddt90_pae1_total_sites
my_plddt90_pae1_structured_fraction <- my_plddt90_pae1_structured_sites / my_plddt90_pae1_total_sites
my_plddt90_pae1_disordered_fraction <- my_plddt90_pae1_disordered_sites / my_plddt90_pae1_total_sites
my_plddt90_pae1_fully_buried_control_fraction
my_plddt90_pae1_fully_buried_fraction
my_plddt90_pae1_buried_fraction
my_plddt90_pae1_structured_fraction
my_plddt90_pae1_disordered_fraction
my_plddt90_pae1_buried_fraction <- scales::percent(my_plddt90_pae1_buried_fraction, accuracy = 0.1)
my_plddt90_pae1_structured_fraction <- scales::percent(my_plddt90_pae1_structured_fraction, accuracy = 0.1)
my_plddt90_pae1_disordered_fraction <- scales::percent(my_plddt90_pae1_disordered_fraction, accuracy = 0.1)
(q %>%
    filter(plddt >= 90 & min_pae <= 1) %>%
    select(acc, site, relasa, ptmbin, freq) %>% unique %>%
    ggplot(aes(x = relasa, colour = ptmbin, fill = ptmbin, weight = freq)) +
    geom_density(alpha = 0.3, bounds = c(0, 1)) +
    annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 9), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 9), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "text", x = 0.125, y = 7.5, label = "Buried", size = 5/.pt) +
    annotate(geom = "text", x = 0.40, y = 7.5, label = "Structured", size = 5/.pt) +
    annotate(geom = "text", x = 0.775, y = 7.5, label = "Disordered", size = 5/.pt) +
    annotate(geom = "text", x = 0.125, y = 6.5, label = my_plddt90_pae1_buried_fraction, size = 5/.pt) +
    annotate(geom = "text", x = 0.40, y = 6.5, label = my_plddt90_pae1_structured_fraction, size = 5/.pt) +
    annotate(geom = "text", x = 0.775, y = 6.5, label = my_plddt90_pae1_disordered_fraction, size = 5/.pt) +
    scale_colour_manual(values = ptmcol, aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = pretty_breaks(2)) +
    theme_minimal() +
    labs(tag = "c") +
    xlab("RSA (pLDDT ≥ 90, PAE ≤ 1)") +
    ylab("Probability density") +
    theme_nature(legend_position = "top") +
    theme(axis.title.y = element_text(hjust = 0))
) %>% qsave("output-density-canonical-weighted-plddt90_pae1.pdf", height = 30)

# # Unweighted (old, don't use)
# (q %>%
#   ggplot(aes(x = relasa, colour = ptmbin, fill = ptmbin)) +
#   geom_density(alpha = 0.3) +
#   annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 2.5), linetype = "dashed", linewidth = 0.5/.weight) +
#   annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 2.5), linetype = "dashed", linewidth = 0.5/.weight) +
#   annotate(geom = "text", x = 0.125, y = 0.15, label = "Buried", size = 5/.pt) +
#   annotate(geom = "text", x = 0.40, y = 0.15, label = "Structured", size = 5/.pt) +
#   annotate(geom = "text", x = 0.775, y = 0.15, label = "Disordered", size = 5/.pt) +
#   scale_colour_manual(values = c("Modified" = ptmb, "Control" = ptmg), aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
#   scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = pretty_breaks(2)) +
#   theme_minimal() +
#   xlab("Relative accessible surface area") +
#   ylab("Probability density") +
#   theme_nature(legend_position = "top")
# ) %>% qsave("output-density-canonical.pdf", height = 30)

# All PTMs ≥1000 sites relasa (unsmoothed) with ptm breakdown
# q
# q %>% filter(ptm != '0') %>% group_by(ptm) %>% tally %>% arrange(desc(n)) %>% filter(n >= 950) %>% pull(ptm) -> tmpptms
# q %>% filter(ptm != '0') %>% group_by(ptm) %>% tally %>% arrange(desc(n)) %>% filter(n >= 950) %>% pull(ptm)
# q %>% filter(ptm != '0') %>% group_by(ptm) %>% tally %>% arrange(desc(n)) %>% filter(n >= 950)
# q %>% filter(ptm != '0') %>% group_by(ptm, buried = relasa <= 0.25) %>% tally %>% pivot_wider(names_from = buried, values_from = n) %>% rename(buried = `TRUE`, surface = `FALSE`) %>% mutate(buried = replace_na(buried, 0), surface = replace_na(surface, 0)) %>% mutate(n = buried + surface) %>% arrange(desc(buried)) %>% filter(n >= 950) %>% pull(ptm)
# q %>% filter(ptm != '0') %>% group_by(ptm, buried = relasa <= 0.25) %>% tally %>% pivot_wider(names_from = buried, values_from = n) %>% rename(buried = `TRUE`, surface = `FALSE`) %>% mutate(buried = replace_na(buried, 0), surface = replace_na(surface, 0)) %>% mutate(n = buried + surface, buriedfrac = buried / n) %>% arrange(desc(buriedfrac)) %>% filter(n >= 950)
# # # Arrange PTMs by their buried-ness
# tmpptm <- q %>% filter(ptm != '0') %>% group_by(ptm, buried = relasa <= 0.25) %>% tally %>% pivot_wider(names_from = buried, values_from = n) %>% rename(buried = `TRUE`, surface = `FALSE`) %>% mutate(buried = replace_na(buried, 0), surface = replace_na(surface, 0)) %>% mutate(n = buried + surface, buriedfrac = buried / n) %>% arrange(desc(buriedfrac)) %>% filter(n >= 950) %>% pull(ptm)
# tmpptm
# # Arrange PTMs by their average RSA (simpler and looks better)
# tmpptm <- q %>% filter(ptmbin == "Modified") %>% select(acc, site, relasa, ptm) %>% group_by(ptm) %>% summarise(mean_relasa = mean(relasa)) %>% arrange(mean_relasa) %>% pull(ptm)
# Arrange PTMs by their median RSA
tmpptm <- q %>% filter(ptmbin == "Modified") %>% select(acc, site, relasa, ptm) %>% group_by(ptm) %>% summarise(median_relasa = median(relasa)) %>% arrange(median_relasa) %>% pull(ptm)
tmpptm
tmpptm %>% length
(q %>%
  select(acc, site, relasa, ptm, freq) %>% unique %>%
  # filter(ptm %in% (q %>% filter(ptm != '0') %>% group_by(ptm) %>% tally %>% arrange(desc(n)) %>% filter(n >= 950) %>% pull(ptm))) %>%
  # filter(ptm %in% (q %>% filter(ptm != '0') %>% group_by(ptm, buried = relasa <= 0.25) %>% tally %>% pivot_wider(names_from = buried, values_from = n) %>% rename(buried = `TRUE`, surface = `FALSE`) %>% mutate(buried = replace_na(buried, 0), surface = replace_na(surface, 0)) %>% mutate(n = buried + surface, buriedfrac = buried / n) %>% arrange(desc(buriedfrac)) %>% filter(n >= 950) %>% pull(ptm))) %>%
  # filter(ptm %in% (tmpptm)) %>%
  # mutate(ptm = fct_rev(factor(ptm, levels = tmpptm))) %>%
  mutate(ptm = ifelse(ptm == 0, "Control", ptm)) %>%
  mutate(ptm = fct_rev(factor(ptm, levels = c("Control", tmpptm)))) %>%
    # ggplot(aes(x = relasa, colour = ptm, fill = ptm)) +
    ggplot(aes(x = relasa, colour = ptm, fill = ptm, weight = freq)) +
    # geom_density(alpha = 0.3) +
    geom_density(alpha = 0.3, bounds = c(0, 1)) +
    geom_vline(xintercept = 0.25, linetype = "dashed", linewidth = 0.5/.weight) +
    geom_vline(xintercept = 0.55, linetype = "dashed", linewidth = 0.5/.weight) +
    # annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 2.75), linetype = "dashed", linewidth = 0.5/.weight) +
    # annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 2.75), linetype = "dashed", linewidth = 0.5/.weight) +
    # annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 2.25), linetype = "dotted", linewidth = 0.5/.weight) +
    # annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 2.25), linetype = "dotted", linewidth = 0.5/.weight) +
    # annotate(geom = "text", x = 0.125, y = 0.15, label = "Buried", size = 5/.pt) +
    # annotate(geom = "text", x = 0.40, y = 0.15, label = "Structured", size = 5/.pt) +
    # annotate(geom = "text", x = 0.775, y = 0.15, label = "Disordered", size = 5/.pt) +
    # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
    # scale_colour_manual(values = c("Yes" = myo, "No" = mybd), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
    # scale_colour_manual(values = c("Yes" = ptmb, "No" = ptmg), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
    # scale_colour_manual(values = c("Modified" = ptmb, "Control" = ptmg), aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
    # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
    scale_fill_manual(aesthetics = c("colour", "fill"), values = c(viridis_pal()(length(tmpptm) + 1)[-1], "grey"), name = NULL, guide = guide_legend(reverse = T)) +
    # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL) +
    # scale_colour_viridis_d(aesthetics = c("colour"), name = NULL, guide = guide_legend(reverse = T)) +
    # scale_colour_viridis_d(name = NULL) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = NULL) +
    # scale_y_continuous(expand = c(0, 0), breaks = pretty_breaks(1)) +
    # scale_x_continuous(limits = c(0, 1)) +
    # coord_cartesian(xlim = c(0, 1)) +
    # coord_cartesian(xlim = c(0, 1), ylim = c(0, 3)) +
    facet_grid(rows = vars(ptm), scales = "free_y") +
    guides(colour = "none", fill = "none") +
    # ggtitle("") +
    labs(tag = "c") +
    xlab("Relative accessible surface area") +
    ylab("Probability density") +
    theme_nature() +
    theme(legend.text = element_text(size = 5), strip.text.y = element_text(size = 5), legend.key.size = unit(1, "mm"), axis.line.y = element_blank())
) %>% qsave("output-density-canonical-allptms-weighted.pdf", height = 80)

# Filtered for pLDDT ≥ 70
# tmpptm <- q %>% filter(plddt >= 70) %>% filter(ptmbin == "Modified") %>% select(acc, site, relasa, ptm) %>% group_by(ptm) %>% summarise(mean_relasa = mean(relasa)) %>% arrange(mean_relasa) %>% pull(ptm)
# Arrange PTMs by their median RSA
tmpptm <- q %>% filter(plddt >= 70) %>% filter(ptmbin == "Modified") %>% select(acc, site, relasa, ptm) %>% group_by(ptm) %>% summarise(median_relasa = median(relasa)) %>% arrange(median_relasa) %>% pull(ptm)
tmpptm
(q %>%
    filter(plddt >= 70) %>%
    select(acc, site, relasa, ptm, freq) %>% unique %>%
    mutate(ptm = ifelse(ptm == 0, "Control", ptm)) %>%
    mutate(ptm = fct_rev(factor(ptm, levels = c("Control", tmpptm)))) %>%
      ggplot(aes(x = relasa, colour = ptm, fill = ptm, weight = freq)) +
      # geom_density(alpha = 0.3) +
      geom_density(alpha = 0.3, bounds = c(0, 1)) +
      geom_vline(xintercept = 0.25, linetype = "dashed", linewidth = 0.5/.weight) +
      geom_vline(xintercept = 0.55, linetype = "dashed", linewidth = 0.5/.weight) +
      scale_fill_manual(aesthetics = c("colour", "fill"), values = c(viridis_pal()(length(tmpptm) + 1)[-1], "grey"), name = NULL, guide = guide_legend(reverse = T)) +
      scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = NULL) +
      facet_grid(rows = vars(ptm), scales = "free_y") +
      guides(colour = "none", fill = "none") +
      labs(tag = "i") +
      xlab("RSA (pLDDT ≥ 70, confident)") +
      ylab("Probability density") +
      theme_nature() +
      theme(legend.text = element_text(size = 5), strip.text.y = element_text(size = 5), legend.key.size = unit(1, "mm"), axis.line.y = element_blank())
) %>% qsave("output-density-canonical-allptms-weighted-plddt70.pdf", height = 80)

# Filtered for pLDDT ≥ 90
# tmpptm <- q %>% filter(plddt >= 90) %>% filter(ptmbin == "Modified") %>% select(acc, site, relasa, ptm) %>% group_by(ptm) %>% summarise(mean_relasa = mean(relasa)) %>% arrange(mean_relasa) %>% pull(ptm)
# Arrange PTMs by their median RSA
tmpptm <- q %>% filter(plddt >= 90) %>% filter(ptmbin == "Modified") %>% select(acc, site, relasa, ptm) %>% group_by(ptm) %>% summarise(median_relasa = median(relasa)) %>% arrange(median_relasa) %>% pull(ptm)
tmpptm
(q %>%
    filter(plddt >= 90) %>%
    select(acc, site, relasa, ptm, freq) %>% unique %>%
    mutate(ptm = ifelse(ptm == 0, "Control", ptm)) %>%
    mutate(ptm = fct_rev(factor(ptm, levels = c("Control", tmpptm)))) %>%
      ggplot(aes(x = relasa, colour = ptm, fill = ptm, weight = freq)) +
      # geom_density(alpha = 0.3) +
      geom_density(alpha = 0.3, bounds = c(0, 1)) +
      geom_vline(xintercept = 0.25, linetype = "dashed", linewidth = 0.5/.weight) +
      geom_vline(xintercept = 0.55, linetype = "dashed", linewidth = 0.5/.weight) +
      scale_fill_manual(aesthetics = c("colour", "fill"), values = c(viridis_pal()(length(tmpptm) + 1)[-1], "grey"), name = NULL, guide = guide_legend(reverse = T)) +
      scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = NULL) +
      facet_grid(rows = vars(ptm), scales = "free_y") +
      guides(colour = "none", fill = "none") +
      labs(tag = "j") +
      xlab("RSA (pLDDT ≥ 90, highly confident)") +
      ylab("Probability density") +
      theme_nature() +
      theme(legend.text = element_text(size = 5), strip.text.y = element_text(size = 5), legend.key.size = unit(1, "mm"), axis.line.y = element_blank())
) %>% qsave("output-density-canonical-allptms-weighted-plddt90.pdf", height = 80)

# q %>% 
#   select(acc, site, relasa, ptm, freq) %>% unique %>%
#   mutate(ptm = ifelse(ptm == 0, "Control", ptm)) %>%
#   mutate(ptm = fct_rev(factor(ptm, levels = c("Control", tmpptm)))) %>%
#   filter(ptm == "K-ac")
# # Verification:
# q %>%
#   select(acc, site, relasa, ptm, freq) %>% unique %>%
#   mutate(ptm = ifelse(ptm == 0, "Control", ptm)) %>%
#   mutate(ptm = fct_rev(factor(ptm, levels = c("Control", tmpptm)))) %>%
#   filter(ptm == "K-ac") %>% select(acc, site) %>% unique
# # >> OK, acc|sites are unique!

# # Unweighted (old)
# tmpptm <- q %>% filter(ptmbin == "Modified") %>% group_by(ptm) %>% summarise(mean_relasa = mean(relasa)) %>% arrange(mean_relasa) %>% pull(ptm)
# (q %>%
#   select(acc, site, relasa, ptm, freq) %>% unique %>%
#   mutate(ptm = ifelse(ptm == 0, "Control", ptm)) %>%
#   mutate(ptm = fct_rev(factor(ptm, levels = c("Control", tmpptm)))) %>%
#     ggplot(aes(x = relasa, colour = ptm, fill = ptm)) +
#     geom_density(alpha = 0.3) +
#     geom_vline(xintercept = 0.25, linetype = "dashed", linewidth = 0.5/.weight) +
#     geom_vline(xintercept = 0.55, linetype = "dashed", linewidth = 0.5/.weight) +
#     scale_fill_manual(aesthetics = c("colour", "fill"), values = c(viridis_pal()(length(tmpptm) + 1)[-1], "grey"), name = NULL, guide = guide_legend(reverse = T)) +
#     scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
#     scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = NULL) +
#     facet_grid(rows = vars(ptm), scales = "free_y") +
#     theme_minimal() +
#     xlab("Relative accessible surface area") +
#     ylab("Probability density") +
#     theme_nature() +
#     guides(colour = "none", fill = "none") +
#     theme(legend.text = element_text(size = 5), strip.text.y = element_text(size = 5), legend.key.size = unit(1, "mm"), axis.line.y = element_blank())
# ) %>% qsave("output-density-canonical-allptms.pdf", height = 80)

# Source breakdown for relasa (unsmoothed)
(
  q %>%
    filter(ptmbin == "Modified") %>%
    # mutate(source = ifelse(source == "PhosphoSitePlus", "PSP", source)) %>%
    mutate(source = fct_relevel(source, "dbPTM", "PhosphoSitePlus", "Ochoa", "UniProt")) %>%
    select(acc, site, source, relasa) %>% unique %>%
    ggplot(aes(x = relasa, colour = source, fill = source)) +
    # ggplot(aes(x = relasa, colour = source)) +
    geom_density(alpha = 0.3) +
    annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 3.2), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 3.2), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "text", x = 0.125, y = 0.2, label = "Buried", size = 5/.pt) +
    annotate(geom = "text", x = 0.40, y = 0.2, label = "Structured", size = 5/.pt) +
    annotate(geom = "text", x = 0.775, y = 0.2, label = "Disordered", size = 5/.pt) +
    scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = pretty_breaks(2)) +
    theme_minimal() +
    xlab("Relative accessible surface area") +
    ylab("Probability density") +
    theme_nature(legend_position = "topleft") +
    # guides(colour = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))
    guides(colour = guide_legend(nrow = 2, byrow = T), fill = guide_legend(nrow = 2, byrow = T)) +
    # theme(legend.spacing.x = unit(0.5, "mm")) +
    theme(legend.spacing.y = unit(0.5, "mm"))
    # theme(legend.direction = "vertical")
    # theme(legend.position = "bottom", legend.direction = "vertical") +
    # theme(legend.text = element_text(size = 5), legend.key.size = unit(1, "mm"))
) %>% qsave("output-density-canonical-sources.pdf", height = 30)

# Source breakdown for relasa (unsmoothed) (common PTMs only, i.e. S-p T-p Y-p)
# q %>% filter(ptmbin == "Modified") %>% group_by(ptm) %>% summarise(sources = n_distinct(source)) %>% arrange(desc(sources)) %>% filter(sources == max(sources)) %>% pull(ptm)
tmpptm <- q %>% select(acc, site, ptmbin, ptm, source) %>% unique %>% filter(ptmbin == "Modified") %>% group_by(ptm) %>% summarise(sources = n_distinct(source)) %>% arrange(desc(sources)) %>% filter(sources == max(sources)) %>% pull(ptm)
tmpptm
(
  q %>%
    filter(ptmbin == "Modified") %>%
    filter(ptm %in% tmpptm) %>%
    select(acc, site, source, relasa) %>% unique %>%
  # mutate(source = ifelse(source == "PhosphoSitePlus", "PSP", source)) %>%
    mutate(source = fct_relevel(source, "dbPTM", "PhosphoSitePlus", "Ochoa", "UniProt")) %>%
    select(acc, site, source, relasa) %>% unique %>%
    ggplot(aes(x = relasa, colour = source, fill = source)) +
    # ggplot(aes(x = relasa, colour = source)) +
    geom_density(alpha = 0.3) +
    annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 3.2), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 3.2), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "text", x = 0.125, y = 0.2, label = "Buried", size = 5/.pt) +
    annotate(geom = "text", x = 0.40, y = 0.2, label = "Structured", size = 5/.pt) +
    annotate(geom = "text", x = 0.775, y = 0.2, label = "Disordered", size = 5/.pt) +
    scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = pretty_breaks(2)) +
    theme_minimal() +
    xlab("Relative accessible surface area") +
    ylab("Probability density") +
    theme_nature(legend_position = "topleft") +
    # guides(colour = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))
    guides(colour = guide_legend(nrow = 2, byrow = T), fill = guide_legend(nrow = 2, byrow = T)) +
    # theme(legend.spacing.x = unit(0.5, "mm")) +
    theme(legend.spacing.y = unit(0.5, "mm"))
    # theme(legend.direction = "vertical")
    # theme(legend.position = "bottom", legend.direction = "vertical") +
    # theme(legend.text = element_text(size = 5), legend.key.size = unit(1, "mm"))
) %>% qsave(f("output-density-canonical-sources-{str_flatten(tmpptm, collapse='_')}-only.pdf"), height = 30)

# Source subset breakdown for relasa (unsmoothed)
(
  q %>%
    mutate(subset = str_replace(subset, "_(small|large)$", " (\\1-scale)")) %>%
    mutate(subset = ifelse(ptm == 0, "Control", subset)) %>%
    mutate(subset = fct_relevel(subset, "Control", after = Inf)) %>%
    select(acc, site, subset, relasa, freq) %>% unique %>%
    # select(acc, site, subset, relasa) %>% unique %>%
    # filter(ptm %in% (q %>% filter(ptm != '0') %>% group_by(ptm) %>% tally %>% arrange(desc(n)) %>% filter(n >= 950) %>% pull(ptm))) %>%
    # filter(ptm %in% (q %>% filter(ptm != '0') %>% group_by(ptm, buried = relasa <= 0.25) %>% tally %>% pivot_wider(names_from = buried, values_from = n) %>% rename(buried = `TRUE`, surface = `FALSE`) %>% mutate(buried = replace_na(buried, 0), surface = replace_na(surface, 0)) %>% mutate(n = buried + surface, buriedfrac = buried / n) %>% arrange(desc(buriedfrac)) %>% filter(n >= 950) %>% pull(ptm))) %>%
    # filter(ptm %in% (tmpptm)) %>%
    # mutate(ptm = fct_rev(factor(ptm, levels = tmpptm))) %>%
    # ggplot(aes(x = relasa, colour = subset, fill = subset)) +
    ggplot(aes(x = relasa, colour = subset, fill = subset, weight = freq)) +
    geom_density(alpha = 0.3) +
    geom_vline(xintercept = 0.25, linetype = "dashed", linewidth = 0.5/.weight) +
    geom_vline(xintercept = 0.55, linetype = "dashed", linewidth = 0.5/.weight) +
    # annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 2.75), linetype = "dashed", linewidth = 0.5/.weight) +
    # annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 2.75), linetype = "dashed", linewidth = 0.5/.weight) +
    # annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 2.25), linetype = "dotted", linewidth = 0.5/.weight) +
    # annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 2.25), linetype = "dotted", linewidth = 0.5/.weight) +
    # annotate(geom = "text", x = 0.125, y = 0.15, label = "Buried", size = 5/.pt) +
    # annotate(geom = "text", x = 0.40, y = 0.15, label = "Structured", size = 5/.pt) +
    # annotate(geom = "text", x = 0.775, y = 0.15, label = "Disordered", size = 5/.pt) +
    # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
    # scale_colour_manual(values = c("Yes" = myo, "No" = mybd), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
    # scale_colour_manual(values = c("Yes" = ptmb, "No" = ptmg), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
    # scale_colour_manual(values = c("Modified" = ptmb, "Control" = ptmg), aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
    # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
    scale_fill_manual(aesthetics = c("colour", "fill"), values = c(viridis_pal()(length(q %>% pull(subset) %>% unique))[-1], "grey"), name = NULL, guide = guide_legend(reverse = T)) +
    # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL) +
    # scale_colour_viridis_d(aesthetics = c("colour"), name = NULL, guide = guide_legend(reverse = T)) +
    # scale_colour_viridis_d(name = NULL) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = NULL) +
    # scale_y_continuous(expand = c(0, 0), breaks = pretty_breaks(1)) +
    # scale_x_continuous(limits = c(0, 1)) +
    # coord_cartesian(xlim = c(0, 1)) +
    # coord_cartesian(xlim = c(0, 1), ylim = c(0, 3)) +
    facet_grid(rows = vars(subset), scales = "free_y") +
    theme_minimal() +
    # ggtitle("") +
    # labs(tag = "a") +
    xlab("Relative accessible surface area") +
    ylab("Probability density") +
    guides(colour = "none", fill = "none") +
    # theme_nature()
    theme_nature() +
    theme(legend.text = element_text(size = 5), strip.text.y = element_text(size = 5), legend.key.size = unit(1, "mm"), axis.line.y = element_blank()) +
    theme(panel.spacing = unit(0.5, "mm")) +
    theme(axis.title.x = element_text(hjust = 0))
) %>% qsave("output-density-canonical-subsets.pdf", width = 50, height = 30)

# (Not actually used in plot below, just for using the numbers in the manuscript)
# relasa (unsmoothed) with ptm breakdown (buried PTMs ≥ 1000, manually excluding C-glt, C-nit, C-pal, M-ox)
# Note: This doesn't need to be weighted since it doesn't include control residues. Weight is included here to avoid confusion, but it makes no difference.
my_fully_buried_controls_y <- q %>% filter(ptmbin == "Control" & aa == "Y" & dis == "strcore" & relasa == 0) %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_fully_buried_sites_yp <- q %>% filter(ptmbin == "Modified" & ptm == "Y-p" & dis == "strcore" & relasa == 0) %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_buried_controls_y <- q %>% filter(ptmbin == "Control" & aa == "Y" & dis == "strcore" & relasa <= 0.25) %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_buried_sites_yp <- q %>% filter(ptmbin == "Modified" & ptm == "Y-p" & dis == "strcore" & relasa <= 0.25) %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_total_sites_yp <- q %>% filter(ptmbin == "Modified" & ptm == "Y-p") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_total_controls_y <- q %>% filter(ptmbin == "Control" & aa == "Y") %>% select(acc, site, relasa, ptmbin, freq) %>% unique %>% nrow
my_fully_buried_controls_y
my_fully_buried_sites_yp
my_buried_controls_y
my_buried_sites_yp
my_total_sites_yp
my_total_controls_y
my_fully_buried_control_fraction_y <- my_fully_buried_controls_y / my_total_controls_y
my_fully_buried_fraction_yp <- my_fully_buried_sites_yp / my_total_sites_yp
my_buried_fraction_p <- my_buried_sites_yp / my_total_sites_yp
my_buried_control_fraction_y <- my_buried_controls_y / my_total_controls_y
my_buried_control_fraction_y
round(my_buried_control_fraction_y, 3)
# >> EXPORT
my_fully_buried_control_fraction_y
round(my_fully_buried_control_fraction_y, 3)
# >> 5.4% of control Y are fully buried (!)
my_fully_buried_fraction_yp
round(my_fully_buried_fraction_yp, 3)
# >> 3.8% of Y-p sites are fully buried (!) That's a really substantial number compared to control!
# >> EXPORT
my_buried_fraction_p
# >> 43.3% of Y-p sites are buried (!) That's a really substantial number compared to control!
# >> EXPORT

buriedptms <- q %>% select(acc, site, ptm, ptmbin, relasa) %>% unique %>% filter(relasa <= 0.25) %>% filter(ptmbin == 'Modified') %>% group_by(ptm) %>% tally %>% arrange(desc(n)) %>% filter(n >= 950) %>% pull(ptm)

# q
# q %>% group_by(ptm, relasa <= 0.25) %>% tally %>% pivot_wider(names_from = `relasa <= 0.25`, values_from=n) %>% transmute(nburied = `TRUE`, fracburied = `TRUE` / (`FALSE` + `TRUE`)) %>% arrange(desc(fracburied))
# q %>% group_by(ptm, relasa <= 0.25) %>% tally %>% pivot_wider(names_from = `relasa <= 0.25`, values_from=n) %>% transmute(nburied = `TRUE`, fracburied = `TRUE` / (`FALSE` + `TRUE`)) %>% arrange(desc(fracburied)) %>% filter(nburied >= 1000)
# q %>% group_by(ptm, relasa <= 0.25) %>% tally %>% pivot_wider(names_from = `relasa <= 0.25`, values_from=n) %>% transmute(nburied = `TRUE`, fracburied = `TRUE` / (`FALSE` + `TRUE`)) %>% arrange(desc(fracburied)) %>% filter(nburied >= 1000) %>% filter(ptm != '0')
# q %>% group_by(ptm, relasa <= 0.25) %>% tally %>% pivot_wider(names_from = `relasa <= 0.25`, values_from=n) %>% transmute(nburied = `TRUE`, fracburied = `TRUE` / (`FALSE` + `TRUE`)) %>% arrange(desc(fracburied)) %>% filter(nburied >= 100) %>% filter(ptm != '0')
# q %>% filter(relasa <= 0.25) %>% filter(ptm != '0') %>% group_by(ptm) %>% tally %>% arrange(desc(n)) %>% filter(n >= 950) %>% pull(ptm) -> tmpptms
# tmpptm <- buriedptms
tmpptm <- buriedptms
(
  q %>%
  select(acc, site, ptm, relasa, freq) %>% unique %>%
  filter(ptm %in% tmpptm) %>%
  filter(!(ptm %in% c('C-glt', 'C-nit', 'C-pal', 'M-ox'))) %>%
  # ggplot(aes(x = relasa, colour = ptm, fill = ptm)) +
  ggplot(aes(x = relasa, colour = ptm, fill = ptm, weight = freq)) +
  # ggplot(aes(x = relasa, colour = ptm)) +
  # geom_density(alpha = 0.3) +
  geom_density(alpha = 0.3, bounds = c(0, 1)) +
  # geom_vline(xintercept = 0.25, linetype = "dashed") +
  # geom_vline(xintercept = 0.55, linetype = "dashed") +
  annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 2.75), linetype = "dashed", linewidth = 0.5/.weight) +
  annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 2.75), linetype = "dashed", linewidth = 0.5/.weight) +
  # annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 2.25), linetype = "dotted", linewidth = 0.5/.weight) +
  # annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 2.25), linetype = "dotted", linewidth = 0.5/.weight) +
  # annotate(geom = "text", x = 0.125, y = 0.15, label = "Buried", size = 5/.pt) +
  # annotate(geom = "text", x = 0.40, y = 0.15, label = "Structured", size = 5/.pt) +
  # annotate(geom = "text", x = 0.775, y = 0.15, label = "Disordered", size = 5/.pt) +
  annotate(geom = "text", x = 0.125, y = 0.2, label = "Buried", size = 5/.pt) +
  annotate(geom = "text", x = 0.40, y = 0.2, label = "Structured", size = 5/.pt) +
  annotate(geom = "text", x = 0.775, y = 0.2, label = "Disordered", size = 5/.pt) +
  # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
  # scale_colour_manual(values = c("Yes" = myo, "No" = mybd), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
  # scale_colour_manual(values = c("Yes" = ptmb, "No" = ptmg), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
  # scale_colour_manual(values = c("Modified" = ptmb, "Control" = ptmg), aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
  scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
  # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL) +
  # scale_colour_viridis_d(aesthetics = c("colour"), name = NULL, guide = guide_legend(reverse = T)) +
  # scale_colour_viridis_d(name = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks = pretty_breaks(2)) +
  # scale_x_continuous(limits = c(0, 1)) +
  # coord_cartesian(xlim = c(0, 1)) +
  # coord_cartesian(xlim = c(0, 1), ylim = c(0, 3)) +
  # ggtitle("") +
  labs(tag = "b") +
  xlab("Relative accessible surface area") +
  ylab("Probability density") +
  # theme_nature()
  theme_nature(legend_position = "top") +
  # theme(legend.text = element_text(size = 5), legend.key.size = unit(1, "mm")) +
  theme(legend.text = element_text(size = 6), legend.key.size = unit(1, "mm")) +
  theme(axis.title.y = element_text(hjust = 0))
) %>% qsave("output-density-canonical-buriedptms.pdf", height = 30)

# Filter for pLDDT ≥ 70
# Use same set of PTMs as for non-pLDDT-filtered
tmpptm <- buriedptms
(
  q %>%
    filter(plddt >= 70) %>%
    select(acc, site, ptm, relasa, freq) %>% unique %>%
    filter(ptm %in% tmpptm) %>%
    filter(!(ptm %in% c('C-glt', 'C-nit', 'C-pal', 'M-ox'))) %>%
    ggplot(aes(x = relasa, colour = ptm, fill = ptm, weight = freq)) +
    # geom_density(alpha = 0.3) +
    geom_density(alpha = 0.3, bounds = c(0, 1)) +
    annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 6), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 6), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "text", x = 0.125, y = 5, label = "Buried", size = 5/.pt) +
    annotate(geom = "text", x = 0.40, y = 5, label = "Structured", size = 5/.pt) +
    annotate(geom = "text", x = 0.775, y = 5, label = "Disordered", size = 5/.pt) +
    scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks = pretty_breaks(2)) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(tag = "d") +
    xlab("RSA (pLDDT ≥ 70, confident)") +
    ylab("Probability density") +
    theme_nature(legend_position = "top") +
    theme(legend.text = element_text(size = 6), legend.key.size = unit(1, "mm")) +
    theme(axis.title.y = element_text(hjust = 0))
) %>% qsave("output-density-canonical-buriedptms-plddt70.pdf", height = 30)

# Filter for pLDDT ≥ 70 & PAE ≤ 2
# Use same set of PTMs as for non-pLDDT-filtered
tmpptm <- buriedptms
(
  q %>%
    filter(plddt >= 70 & min_pae <= 2) %>%
    select(acc, site, ptm, relasa, freq) %>% unique %>%
    filter(ptm %in% tmpptm) %>%
    filter(!(ptm %in% c('C-glt', 'C-nit', 'C-pal', 'M-ox'))) %>%
    ggplot(aes(x = relasa, colour = ptm, fill = ptm, weight = freq)) +
    # geom_density(alpha = 0.3) +
    geom_density(alpha = 0.3, bounds = c(0, 1)) +
    annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 6), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 6), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "text", x = 0.125, y = 5, label = "Buried", size = 5/.pt) +
    annotate(geom = "text", x = 0.40, y = 5, label = "Structured", size = 5/.pt) +
    annotate(geom = "text", x = 0.775, y = 5, label = "Disordered", size = 5/.pt) +
    scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks = pretty_breaks(2)) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(tag = "f") +
    xlab("RSA (pLDDT ≥ 70, PAE ≤ 2)") +
    ylab("Probability density") +
    theme_nature(legend_position = "top") +
    theme(legend.text = element_text(size = 6), legend.key.size = unit(1, "mm")) +
    theme(axis.title.y = element_text(hjust = 0))
) %>% qsave("output-density-canonical-buriedptms-plddt70_pae2.pdf", height = 30)

# Filter for pLDDT ≥ 90
# Use same set of PTMs as for non-pLDDT-filtered
tmpptm <- buriedptms
(
  q %>%
    filter(plddt >= 90) %>%
    select(acc, site, ptm, relasa, freq) %>% unique %>%
    filter(ptm %in% tmpptm) %>%
    filter(!(ptm %in% c('C-glt', 'C-nit', 'C-pal', 'M-ox'))) %>%
    ggplot(aes(x = relasa, colour = ptm, fill = ptm, weight = freq)) +
    # geom_density(alpha = 0.3) +
    geom_density(alpha = 0.3, bounds = c(0, 1)) +
    annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 7.5), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 7.5), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "text", x = 0.125, y = 6.5, label = "Buried", size = 5/.pt) +
    annotate(geom = "text", x = 0.40, y = 6.5, label = "Structured", size = 5/.pt) +
    annotate(geom = "text", x = 0.775, y = 6.5, label = "Disordered", size = 5/.pt) +
    scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks = pretty_breaks(2)) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(tag = "g") +
    xlab("RSA (pLDDT ≥ 90, highly confident)") +
    ylab("Probability density") +
    theme_nature(legend_position = "top") +
    theme(legend.text = element_text(size = 6), legend.key.size = unit(1, "mm")) +
    theme(axis.title.y = element_text(hjust = 0))
) %>% qsave("output-density-canonical-buriedptms-plddt90.pdf", height = 30)

# Filter for pLDDT ≥ 90 & PAE ≤ 1
# Use same set of PTMs as for non-pLDDT-filtered
tmpptm <- buriedptms
(
  q %>%
    filter(plddt >= 90 & min_pae <= 1) %>%
    select(acc, site, ptm, relasa, freq) %>% unique %>%
    filter(ptm %in% tmpptm) %>%
    filter(!(ptm %in% c('C-glt', 'C-nit', 'C-pal', 'M-ox'))) %>%
    ggplot(aes(x = relasa, colour = ptm, fill = ptm, weight = freq)) +
    # geom_density(alpha = 0.3) +
    geom_density(alpha = 0.3, bounds = c(0, 1)) +
    annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 7.5), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 7.5), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "text", x = 0.125, y = 6.5, label = "Buried", size = 5/.pt) +
    annotate(geom = "text", x = 0.40, y = 6.5, label = "Structured", size = 5/.pt) +
    annotate(geom = "text", x = 0.775, y = 6.5, label = "Disordered", size = 5/.pt) +
    scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks = pretty_breaks(2)) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(tag = "h") +
    xlab("RSA (pLDDT ≥ 90, PAE ≤ 1)") +
    ylab("Probability density") +
    theme_nature(legend_position = "top") +
    theme(legend.text = element_text(size = 6), legend.key.size = unit(1, "mm")) +
    theme(axis.title.y = element_text(hjust = 0))
) %>% qsave("output-density-canonical-buriedptms-plddt90_pae1.pdf", height = 30)

# Source breakdown of relasa (unsmoothed) with ptm breakdown (buried PTMs ≥ 1000, manually excluding C-glt, C-nit, C-pal, M-ox)
# Note: This doesn't need to be weighted since it doesn't include control residues. Weight is included here to avoid confusion, but it makes no difference.
# q
# q %>% group_by(ptm, relasa <= 0.25) %>% tally %>% pivot_wider(names_from = `relasa <= 0.25`, values_from=n) %>% transmute(nburied = `TRUE`, fracburied = `TRUE` / (`FALSE` + `TRUE`)) %>% arrange(desc(fracburied))
# q %>% group_by(ptm, relasa <= 0.25) %>% tally %>% pivot_wider(names_from = `relasa <= 0.25`, values_from=n) %>% transmute(nburied = `TRUE`, fracburied = `TRUE` / (`FALSE` + `TRUE`)) %>% arrange(desc(fracburied)) %>% filter(nburied >= 1000)
# q %>% group_by(ptm, relasa <= 0.25) %>% tally %>% pivot_wider(names_from = `relasa <= 0.25`, values_from=n) %>% transmute(nburied = `TRUE`, fracburied = `TRUE` / (`FALSE` + `TRUE`)) %>% arrange(desc(fracburied)) %>% filter(nburied >= 1000) %>% filter(ptm != '0')
# q %>% group_by(ptm, relasa <= 0.25) %>% tally %>% pivot_wider(names_from = `relasa <= 0.25`, values_from=n) %>% transmute(nburied = `TRUE`, fracburied = `TRUE` / (`FALSE` + `TRUE`)) %>% arrange(desc(fracburied)) %>% filter(nburied >= 100) %>% filter(ptm != '0')
# q %>% filter(relasa <= 0.25) %>% filter(ptm != '0') %>% group_by(ptm) %>% tally %>% arrange(desc(n)) %>% filter(n >= 950) %>% pull(ptm) -> tmpptms
tmpptm <- q %>% select(acc, site, ptm, relasa) %>% unique %>% filter(relasa <= 0.25) %>% filter(ptm != '0') %>% group_by(ptm) %>% tally %>% arrange(desc(n)) %>% filter(n >= 950) %>% pull(ptm)
(
  q %>%
  select(acc, site, ptm, source, relasa, freq) %>% unique %>%
  filter(ptm %in% tmpptm) %>%
  filter(!(ptm %in% c('C-glt', 'C-nit', 'C-pal', 'M-ox'))) %>%
  # ggplot(aes(x = relasa, colour = ptm, fill = ptm)) +
  ggplot(aes(x = relasa, colour = ptm, fill = ptm, weight = freq)) +
  # ggplot(aes(x = relasa, colour = ptm)) +
  geom_density(alpha = 0.3) +
  # geom_vline(xintercept = 0.25, linetype = "dashed") +
  # geom_vline(xintercept = 0.55, linetype = "dashed") +
  annotate(geom = "line", x = c(0.25, 0.25), y = c(0, Inf), linetype = "dashed", linewidth = 0.5/.weight) +
  annotate(geom = "line", x = c(0.55, 0.55), y = c(0, Inf), linetype = "dashed", linewidth = 0.5/.weight) +
  # annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 2.25), linetype = "dotted", linewidth = 0.5/.weight) +
  # annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 2.25), linetype = "dotted", linewidth = 0.5/.weight) +
  annotate(geom = "text", x = 0.125, y = Inf, label = "Buried", size = 5/.pt, vjust = 1) +
  annotate(geom = "text", x = 0.40, y = Inf, label = "Structured", size = 5/.pt, vjust = 1) +
  annotate(geom = "text", x = 0.775, y = Inf, label = "Disordered", size = 5/.pt, vjust = 1) +
  # annotation_custom(grid::textGrob("Special Point"), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
  # scale_colour_manual(values = c("Yes" = myo, "No" = mybd), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
  # scale_colour_manual(values = c("Yes" = ptmb, "No" = ptmg), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
  # scale_colour_manual(values = c("Modified" = ptmb, "Control" = ptmg), aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
  # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
  scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL) +
  # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL) +
  # scale_colour_viridis_d(aesthetics = c("colour"), name = NULL, guide = guide_legend(reverse = T)) +
  # scale_colour_viridis_d(name = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = pretty_breaks(2)) +
  # scale_x_continuous(limits = c(0, 1)) +
  # coord_cartesian(xlim = c(0, 1)) +
  # coord_cartesian(xlim = c(0, 1), ylim = c(0, 3)) +
  facet_wrap(facets = vars(source), ncol = 2, scales = "free") +
  # theme_minimal() +
  # ggtitle("") +
  # labs(tag = "a") +
  xlab("Relative accessible surface area") +
  ylab("Probability density") +
  # theme_nature()
  # theme_nature(legend_position = "topcentre") +
  theme_nature(legend_position = "bottom") +
  theme(axis.title.x = element_text(hjust = 0))
  # theme(legend.spacing.x = unit(4, "mm"))
  # theme(legend.text = element_text(size = 5), legend.key.size = unit(1, "mm"))
) %>% qsave("output-density-canonical-buriedptms-by-source.pdf", width = 80, height = 50)

# Subset breakdown of relasa (unsmoothed) with ptm breakdown (buried PTMs ≥ 1000, manually excluding C-glt, C-nit, C-pal, M-ox)
# Note: This doesn't need to be weighted since it doesn't include control residues. Weight is included here to avoid confusion, but it makes no difference.
# q
# q %>% group_by(ptm, relasa <= 0.25) %>% tally %>% pivot_wider(names_from = `relasa <= 0.25`, values_from=n) %>% transmute(nburied = `TRUE`, fracburied = `TRUE` / (`FALSE` + `TRUE`)) %>% arrange(desc(fracburied))
# q %>% group_by(ptm, relasa <= 0.25) %>% tally %>% pivot_wider(names_from = `relasa <= 0.25`, values_from=n) %>% transmute(nburied = `TRUE`, fracburied = `TRUE` / (`FALSE` + `TRUE`)) %>% arrange(desc(fracburied)) %>% filter(nburied >= 1000)
# q %>% group_by(ptm, relasa <= 0.25) %>% tally %>% pivot_wider(names_from = `relasa <= 0.25`, values_from=n) %>% transmute(nburied = `TRUE`, fracburied = `TRUE` / (`FALSE` + `TRUE`)) %>% arrange(desc(fracburied)) %>% filter(nburied >= 1000) %>% filter(ptm != '0')
# q %>% group_by(ptm, relasa <= 0.25) %>% tally %>% pivot_wider(names_from = `relasa <= 0.25`, values_from=n) %>% transmute(nburied = `TRUE`, fracburied = `TRUE` / (`FALSE` + `TRUE`)) %>% arrange(desc(fracburied)) %>% filter(nburied >= 100) %>% filter(ptm != '0')
# q %>% filter(relasa <= 0.25) %>% filter(ptm != '0') %>% group_by(ptm) %>% tally %>% arrange(desc(n)) %>% filter(n >= 950) %>% pull(ptm) -> tmpptms
# q %>% select(acc, site, ptm, subset) %>% unique %>% filter(ptm %in% tmpptm) %>% group_by(subset, ptm) %>% tally %>% arrange(subset, ptm) %>% print(n = 50)
tmpptm <- q %>% select(acc, site, ptm, ptmbin, relasa) %>% unique %>% filter(relasa <= 0.25) %>% filter(ptmbin == 'Modified') %>% group_by(ptm) %>% tally %>% arrange(desc(n)) %>% filter(n >= 950) %>% filter(!(ptm %in% c('C-glt', 'C-nit', 'C-pal', 'M-ox'))) %>% pull(ptm)
tmpptm
(
  q %>%
    select(acc, site, ptm, subset, relasa, freq) %>% unique %>%
    filter(ptm %in% tmpptm) %>%
    # Too few sites in UniProt_large (but _small is okay)
    filter(subset != "UniProt_large") %>%
    mutate(subset = str_replace(subset, "_(small|large)$", " (\\1-scale)")) %>%
    select(acc, site, ptm, subset, relasa, freq) %>% unique %>%
    # filter(!(ptm %in% c('C-glt', 'C-nit', 'C-pal', 'M-ox'))) %>%
    # ggplot(aes(x = relasa, colour = ptm, fill = ptm)) +
    ggplot(aes(x = relasa, colour = ptm, fill = ptm, weight = freq)) +
    # ggplot(aes(x = relasa, colour = ptm)) +
    geom_density(alpha = 0.3) +
    # geom_vline(xintercept = 0.25, linetype = "dashed") +
    # geom_vline(xintercept = 0.55, linetype = "dashed") +
    annotate(geom = "line", x = c(0.25, 0.25), y = c(0, Inf), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "line", x = c(0.55, 0.55), y = c(0, Inf), linetype = "dashed", linewidth = 0.5/.weight) +
    # annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 2.25), linetype = "dotted", linewidth = 0.5/.weight) +
    # annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 2.25), linetype = "dotted", linewidth = 0.5/.weight) +
    annotate(geom = "text", x = 0.125, y = Inf, label = "Buried", size = 5/.pt, vjust = 1) +
    annotate(geom = "text", x = 0.40, y = Inf, label = "Structured", size = 5/.pt, vjust = 1) +
    annotate(geom = "text", x = 0.775, y = Inf, label = "Disordered", size = 5/.pt, vjust = 1) +
    # annotation_custom(grid::textGrob("Special Point"), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
    # scale_colour_manual(values = c("Yes" = myo, "No" = mybd), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
    # scale_colour_manual(values = c("Yes" = ptmb, "No" = ptmg), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
    # scale_colour_manual(values = c("Modified" = ptmb, "Control" = ptmg), aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
    # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
    scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL) +
    # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL) +
    # scale_colour_viridis_d(aesthetics = c("colour"), name = NULL, guide = guide_legend(reverse = T)) +
    # scale_colour_viridis_d(name = NULL) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = pretty_breaks(2)) +
    # scale_x_continuous(limits = c(0, 1)) +
    # coord_cartesian(xlim = c(0, 1)) +
    # coord_cartesian(xlim = c(0, 1), ylim = c(0, 3)) +
    facet_wrap(facets = vars(subset), ncol = 3, scales = "free", dir = "v") +
    theme_minimal() +
    # ggtitle("") +
    # labs(tag = "a") +
    xlab("Relative accessible surface area") +
    ylab("Probability density") +
    # theme_nature()
    # theme_nature(legend_position = "topcentre") +
    theme_nature(legend_position = "bottom")
    # theme(axis.title.x = element_text(hjust = 0))
    # theme(legend.spacing.x = unit(4, "mm"))
    # theme(legend.text = element_text(size = 5), legend.key.size = unit(1, "mm"))
) %>% qsave("output-density-canonical-buriedptms-by-subset.pdf", width = 120, height = 50)

# Vertical subset breakdown of relasa (unsmoothed) with ptm breakdown (buried PTMs ≥ 1000, manually excluding C-glt, C-nit, C-pal, M-ox)
# Note: This doesn't need to be weighted since it doesn't include control residues. Weight is included here to avoid confusion, but it makes no difference.
# q
# q %>% group_by(ptm, relasa <= 0.25) %>% tally %>% pivot_wider(names_from = `relasa <= 0.25`, values_from=n) %>% transmute(nburied = `TRUE`, fracburied = `TRUE` / (`FALSE` + `TRUE`)) %>% arrange(desc(fracburied))
# q %>% group_by(ptm, relasa <= 0.25) %>% tally %>% pivot_wider(names_from = `relasa <= 0.25`, values_from=n) %>% transmute(nburied = `TRUE`, fracburied = `TRUE` / (`FALSE` + `TRUE`)) %>% arrange(desc(fracburied)) %>% filter(nburied >= 1000)
# q %>% group_by(ptm, relasa <= 0.25) %>% tally %>% pivot_wider(names_from = `relasa <= 0.25`, values_from=n) %>% transmute(nburied = `TRUE`, fracburied = `TRUE` / (`FALSE` + `TRUE`)) %>% arrange(desc(fracburied)) %>% filter(nburied >= 1000) %>% filter(ptm != '0')
# q %>% group_by(ptm, relasa <= 0.25) %>% tally %>% pivot_wider(names_from = `relasa <= 0.25`, values_from=n) %>% transmute(nburied = `TRUE`, fracburied = `TRUE` / (`FALSE` + `TRUE`)) %>% arrange(desc(fracburied)) %>% filter(nburied >= 100) %>% filter(ptm != '0')
# q %>% filter(relasa <= 0.25) %>% filter(ptm != '0') %>% group_by(ptm) %>% tally %>% arrange(desc(n)) %>% filter(n >= 950) %>% pull(ptm) -> tmpptms
tmpptm <- q %>% select(acc, site, ptm, ptmbin, relasa) %>% unique %>% filter(relasa <= 0.25) %>% filter(ptmbin == 'Modified') %>% group_by(ptm) %>% tally %>% arrange(desc(n)) %>% filter(n >= 950) %>% filter(!(ptm %in% c('C-glt', 'C-nit', 'C-pal', 'M-ox'))) %>% pull(ptm)
tmpptm
q %>% filter(ptm %in% tmpptm) %>% group_by(subset, ptm) %>% tally %>% arrange(subset, ptm) %>% print(n = 50)
(
  q %>%
    filter(ptm %in% tmpptm) %>%
    # Too few sites in UniProt_large (but _small is okay)
    filter(subset != "UniProt_large") %>%
    mutate(subset = str_replace(subset, "_(small|large)$", " (\\1-scale)")) %>%
    select(acc, site, ptm, subset, relasa, freq) %>% unique %>%
    # filter(!(ptm %in% c('C-glt', 'C-nit', 'C-pal', 'M-ox'))) %>%
    # ggplot(aes(x = relasa, colour = ptm, fill = ptm)) +
    ggplot(aes(x = relasa, colour = ptm, fill = ptm, weight = freq)) +
    # ggplot(aes(x = relasa, colour = ptm)) +
    # geom_density(alpha = 0.3) +
    geom_density(alpha = 0.3, bounds = c(0, 1)) +
    # geom_vline(xintercept = 0.25, linetype = "dashed") +
    # geom_vline(xintercept = 0.55, linetype = "dashed") +
    annotate(geom = "line", x = c(0.25, 0.25), y = c(0, Inf), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "line", x = c(0.55, 0.55), y = c(0, Inf), linetype = "dashed", linewidth = 0.5/.weight) +
    # annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 2.25), linetype = "dotted", linewidth = 0.5/.weight) +
    # annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 2.25), linetype = "dotted", linewidth = 0.5/.weight) +
    annotate(geom = "text", x = 0.125, y = Inf, label = "Buried", size = 5/.pt, vjust = 1) +
    annotate(geom = "text", x = 0.40, y = Inf, label = "Structured", size = 5/.pt, vjust = 1) +
    annotate(geom = "text", x = 0.775, y = Inf, label = "Disordered", size = 5/.pt, vjust = 1) +
    # annotation_custom(grid::textGrob("Special Point"), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
    # scale_colour_manual(values = c("Yes" = myo, "No" = mybd), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
    # scale_colour_manual(values = c("Yes" = ptmb, "No" = ptmg), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
    # scale_colour_manual(values = c("Modified" = ptmb, "Control" = ptmg), aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
    # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
    scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL) +
    # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL) +
    # scale_colour_viridis_d(aesthetics = c("colour"), name = NULL, guide = guide_legend(reverse = T)) +
    # scale_colour_viridis_d(name = NULL) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = pretty_breaks(2)) +
    # scale_x_continuous(limits = c(0, 1)) +
    # coord_cartesian(xlim = c(0, 1)) +
    # coord_cartesian(xlim = c(0, 1), ylim = c(0, 3)) +
    facet_wrap(facets = vars(subset), ncol = 2, scales = "free", dir = "h") +
    theme_minimal() +
    # ggtitle("") +
    labs(tag = "d") +
    xlab("Relative accessible surface area") +
    ylab("Probability density") +
    guides(colour = guide_legend(reverse = T), fill = guide_legend(reverse = T)) +
    # theme_nature()
    # theme_nature(legend_position = "topcentre") +
    theme_nature(legend_position = "bottom") +
    theme(axis.title.x = element_text(hjust = 0))
    # theme(legend.spacing.x = unit(4, "mm")) 
    # theme(legend.text = element_text(size = 5), legend.key.size = unit(1, "mm"))
) %>% qsave("output-density-canonical-buriedptms-by-subset-vertical.pdf", width = 80, height = 80)

# Filtered for pLDDT ≥ 70
# Use same set of PTMs as for non-pLDDT-filtered
tmpptm <- q %>% select(acc, site, ptm, ptmbin, relasa) %>% unique %>% filter(relasa <= 0.25) %>% filter(ptmbin == 'Modified') %>% group_by(ptm) %>% tally %>% arrange(desc(n)) %>% filter(n >= 950) %>% filter(!(ptm %in% c('C-glt', 'C-nit', 'C-pal', 'M-ox'))) %>% pull(ptm)
tmpptm
q %>% filter(ptm %in% tmpptm) %>% group_by(subset, ptm) %>% tally %>% arrange(subset, ptm) %>% print(n = 50)
(
  q %>%
    filter(plddt >= 70) %>%
    filter(ptm %in% tmpptm) %>%
    filter(subset != "UniProt_large") %>%
    mutate(subset = str_replace(subset, "_(small|large)$", " (\\1-scale)")) %>%
    select(acc, site, ptm, subset, relasa, freq) %>% unique %>%
    ggplot(aes(x = relasa, colour = ptm, fill = ptm, weight = freq)) +
    # geom_density(alpha = 0.3) +
    geom_density(alpha = 0.3, bounds = c(0, 1)) +
    annotate(geom = "line", x = c(0.25, 0.25), y = c(0, Inf), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "line", x = c(0.55, 0.55), y = c(0, Inf), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "text", x = 0.125, y = Inf, label = "Buried", size = 5/.pt, vjust = 1) +
    annotate(geom = "text", x = 0.40, y = Inf, label = "Structured", size = 5/.pt, vjust = 1) +
    annotate(geom = "text", x = 0.775, y = Inf, label = "Disordered", size = 5/.pt, vjust = 1) +
    scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = pretty_breaks(2)) +
    facet_wrap(facets = vars(subset), ncol = 2, scales = "free", dir = "h") +
    theme_minimal() +
    labs(tag = "k") +
    xlab("RSA (pLDDT ≥ 70, confident)") +
    ylab("Probability density") +
    theme_nature(legend_position = "bottom") +
    theme(axis.title.x = element_text(hjust = 0))
) %>% qsave("output-density-canonical-buriedptms-by-subset-vertical-plddt70.pdf", width = 80, height = 80)

# Filtered for pLDDT ≥ 90
# Use same set of PTMs as for non-pLDDT-filtered
tmpptm <- q %>% select(acc, site, ptm, ptmbin, relasa) %>% unique %>% filter(relasa <= 0.25) %>% filter(ptmbin == 'Modified') %>% group_by(ptm) %>% tally %>% arrange(desc(n)) %>% filter(n >= 950) %>% filter(!(ptm %in% c('C-glt', 'C-nit', 'C-pal', 'M-ox'))) %>% pull(ptm)
tmpptm
q %>% filter(ptm %in% tmpptm) %>% group_by(subset, ptm) %>% tally %>% arrange(subset, ptm) %>% print(n = 50)
(
  q %>%
    filter(plddt >= 90) %>%
    filter(ptm %in% tmpptm) %>%
    filter(subset != "UniProt_large") %>%
    mutate(subset = str_replace(subset, "_(small|large)$", " (\\1-scale)")) %>%
    select(acc, site, ptm, subset, relasa, freq) %>% unique %>%
    ggplot(aes(x = relasa, colour = ptm, fill = ptm, weight = freq)) +
    # geom_density(alpha = 0.3) +
    geom_density(alpha = 0.3, bounds = c(0, 1)) +
    annotate(geom = "line", x = c(0.25, 0.25), y = c(0, Inf), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "line", x = c(0.55, 0.55), y = c(0, Inf), linetype = "dashed", linewidth = 0.5/.weight) +
    annotate(geom = "text", x = 0.125, y = Inf, label = "Buried", size = 5/.pt, vjust = 1) +
    annotate(geom = "text", x = 0.40, y = Inf, label = "Structured", size = 5/.pt, vjust = 1) +
    annotate(geom = "text", x = 0.775, y = Inf, label = "Disordered", size = 5/.pt, vjust = 1) +
    scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = NULL) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = pretty_breaks(2)) +
    facet_wrap(facets = vars(subset), ncol = 2, scales = "free", dir = "h") +
    theme_minimal() +
    labs(tag = "l") +
    xlab("RSA (pLDDT ≥ 90, highly confident)") +
    ylab("Probability density") +
    theme_nature(legend_position = "bottom") +
    theme(axis.title.x = element_text(hjust = 0))
) %>% qsave("output-density-canonical-buriedptms-by-subset-vertical-plddt90.pdf", width = 80, height = 80)

# relasa (unsmoothed) with ptm breakdown (major PTMs)
q %>%
  filter(ptm %in% c('S-p', 'T-p', 'Y-p', 'K-ac', 'K-ub')) %>%
  # ggplot(aes(x = relasa, colour = ptm, fill = ptm)) +
  ggplot(aes(x = relasa, colour = ptm)) +
  geom_density(alpha = 0.3) +
  # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
  # scale_colour_manual(values = c("Yes" = myo, "No" = mybd), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
  # scale_colour_manual(values = c("Yes" = ptmb, "No" = ptmg), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
  # scale_x_continuous(limits = c(0, 1)) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_minimal() +
  # theme_nature() +
  xlab("AlphaFold relative accessible surface area") +
  ylab("Density")
ggsave("output-density-canonical-STYp-Kac-Kub.pdf", width = 6, height = 4)

# relasa (unsmoothed) with ptm breakdown (my PTMs of interest)
q %>%
  filter(ptm %in% c('S-p', 'T-p', 'Y-p', 'K-ac', 'K-mal', 'K-me', 'K-suc', 'K-sum', 'K-ub', 'N-gly', 'R-me', 'S-gly', 'T-gly')) %>%
  # ggplot(aes(x = relasa, colour = ptm, fill = ptm)) +
  ggplot(aes(x = relasa, colour = ptm)) +
  geom_density(alpha = 0.3) +
  # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
  # scale_colour_manual(values = c("Yes" = myo, "No" = mybd), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
  # scale_colour_manual(values = c("Yes" = ptmb, "No" = ptmg), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
  # scale_x_continuous(limits = c(0, 1)) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_minimal() +
  # theme_nature() +
  xlab("AlphaFold relative accessible surface area") +
  ylab("Density")
ggsave("output-density-canonical-STYp-Kac-Kmal-KRme-Ksuc-Ksum-Kub-NSTgly.pdf", width = 6, height = 4)

# relasa (unsmoothed) with ptm breakdown (minor PTMs)
q %>%
  filter(ptm %in% c('K-sum', 'R-me', 'K-mal', 'K-suc', 'K-me')) %>%
  # ggplot(aes(x = relasa, colour = ptm, fill = ptm)) +
  ggplot(aes(x = relasa, colour = ptm)) +
  geom_density(alpha = 0.3) +
  # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
  # scale_colour_manual(values = c("Yes" = myo, "No" = mybd), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
  # scale_colour_manual(values = c("Yes" = ptmb, "No" = ptmg), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
  # scale_x_continuous(limits = c(0, 1)) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_minimal() +
  # theme_nature() +
  xlab("AlphaFold relative accessible surface area") +
  ylab("Density")
ggsave("output-density-canonical-KRme-Kmal-Ksuc-Ksum.pdf", width = 6, height = 4)




# Histograms

# relasa (unsmoothed) histogram
(
  q %>%
    select(acc, site, ptmbin, relasa) %>% unique %>%
    ggplot(aes(x = relasa, colour = ptmbin, fill = ptmbin)) +
    # Buried (blue) / Structured (green) / Disordered (orange)
    annotate(geom = "rect", xmin = -Inf, xmax = 0.25, ymin = 0, ymax = Inf, fill = ptmbl, alpha = 0.2) +
    annotate(geom = "rect", xmin = 0.25, xmax = 0.55, ymin = 0, ymax = Inf, fill = ptmgr, alpha = 0.2) +
    annotate(geom = "rect", xmin = 0.55, xmax = Inf, ymin = 0, ymax = Inf, fill = ptmod, alpha = 0.2) +
    geom_histogram(binwidth = 0.05, alpha = 0.3, boundary = 0) +
    # geom_vline(xintercept = 0.25, linetype = "dotted") +
    # geom_vline(xintercept = 0.55, linetype = "dotted") +
    geom_vline(xintercept = 0.25, linetype = "dashed") +
    geom_vline(xintercept = 0.55, linetype = "dashed") +
    # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
    # scale_colour_manual(values = c("Yes" = myo, "No" = mybd), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
    # scale_colour_manual(values = c("Yes" = ptmb, "No" = ptmg), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
    # scale_colour_manual(values = c("Yes" = ptmb, "No" = ptmg), aesthetics = c("colour", "fill"), name = "Modified", guide = "none") +
    scale_colour_manual(values = c("Yes" = ptmbd, "No" = ptmg), aesthetics = c("colour", "fill"), name = "Modified", guide = "none") +
    scale_x_continuous(expand = c(0, 0)) +
    # scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(labels = label_comma()) +
    # coord_cartesian(xlim = c(0, 1)) +
    facet_grid(rows = vars(fct_rev(ptmbin))) +
    # theme_minimal() +
    theme_nature() +
    xlab("AlphaFold relative accessible surface area") +
    ylab("Residues")
) %>% qsave("output-histogram-canonical.pdf", width = 60)

# relasa (unsmoothed) histogram small
q %>%
  ggplot(aes(x = relasa, colour = ptmbin, fill = ptmbin)) +
  geom_histogram(binwidth = 0.05, alpha = 0.3) +
  # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
  # scale_colour_manual(values = c("Yes" = myo, "No" = mybd), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
  scale_colour_manual(values = c("Yes" = ptmb, "No" = ptmg), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
  # scale_x_continuous(limits = c(0, 1)) +
  coord_cartesian(xlim = c(0, 1)) +
  facet_grid(rows = vars(fct_rev(ptmbin))) +
  theme_minimal() +
  xlab("AlphaFold-based relative\nsolvent-accessible surface area") +
  ylab("Residues")
ggsave("output-histogram-canonical-small.pdf", width = 2, height = 2)


# relasa (unsmoothed) histogram with ptm breakdown (my PTMs of interest)
(
  q %>%
    select(acc, site, ptm, relasa) %>% unique %>%
    filter(ptm %in% c('S-p', 'T-p', 'Y-p', 'K-ac', 'K-mal', 'K-me', 'K-suc', 'K-sum', 'K-ub', 'N-gly', 'R-me', 'S-gly', 'T-gly')) %>%
    # ggplot(aes(x = relasa, colour = ptm, fill = ptm)) +
    ggplot(aes(x = relasa, fill = ptm)) +
    # geom_histogram(binwidth = 0.05, alpha = 0.3) +
    # geom_histogram(binwidth = 0.05, alpha = 0.3, boundary = 0) +
    geom_histogram(binwidth = 0.05, alpha = 0.7, boundary = 0) +
    # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
    # scale_colour_manual(values = c("Yes" = myo, "No" = mybd), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
    # scale_colour_manual(values = c("Yes" = ptmb, "No" = ptmg), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
    # scale_x_continuous(limits = c(0, 1)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(breaks = pretty_breaks(2)) +
    coord_cartesian(xlim = c(0, 1)) +
    # facet_grid(rows = vars(ptm)) +
    facet_grid(rows = vars(ptm), scales = "free_y") +
    guides(colour = "none", fill = "none") +
    theme_nature() +
    xlab("AlphaFold relative accessible surface area") +
    ylab("Residues")
) %>% qsave("output-histogram-canonical-STYp-Kac-Kmal-KRme-Ksuc-Ksum-Kub-NSTgly.pdf", width = 60, height = 70)

# relasa (unsmoothed) histogram with ptm breakdown (my PTMs of interest) (buried only)
(
  q %>%
    select(acc, site, ptm, relasa) %>% unique %>%
    filter(ptm %in% c('S-p', 'T-p', 'Y-p', 'K-ac', 'K-mal', 'K-me', 'K-suc', 'K-sum', 'K-ub', 'N-gly', 'R-me', 'S-gly', 'T-gly')) %>%
    filter(relasa <= 0.25) %>%
    ggplot(aes(x = relasa, colour = ptm, fill = ptm)) +
    geom_histogram(binwidth = 0.05, alpha = 0.3) +
    # scale_colour_viridis_d(aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
    # scale_colour_manual(values = c("Yes" = myo, "No" = mybd), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
    # scale_colour_manual(values = c("Yes" = ptmb, "No" = ptmg), aesthetics = c("colour", "fill"), name = "Modified", guide = guide_legend(reverse = T)) +
    # scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(breaks = pretty_breaks(2)) +
    coord_cartesian(xlim = c(0, 1)) +
    # facet_grid(rows = vars(ptm), scales) +
    facet_grid(rows = vars(ptm), scales = "free_y") +
    guides(colour = "none", fill = "none") +
    theme_nature() +
    xlab("AlphaFold relative accessible surface area") +
    ylab("Residues")
) %>% qsave("output-histogram-canonical-STYp-Kac-Kmal-KRme-Ksuc-Ksum-Kub-NSTgly-buried.pdf", width = 60, height = 70)










# Example output-density-canonical-weighted.pdf plot from above
# ptm_order_select <- c("S-p", "T-p", "Y-p", "K-ub", "K-ac", "K-mal", "K-suc", "K-me", "R-me")
# mypvalues_select <- mypvalues %>% filter(ptm %in% ptm_order_select)
(q %>%
    # filter(ptm %in% c(ptm_order_select, "0")) %>%
    # filter(ptm %in% ptm_order_select) %>%
    mutate(dis = fct_recode(dis, "Buried"="strcore", "Structured"="strsurf", "Disordered"="dissurf")) %>%
    select(acc, site, relasa, ptmbin, freq, dis, rate) %>% unique %>%
    # ggplot(aes(x = relasa, colour = ptmbin, fill = ptmbin, weight = freq)) +
    ggplot(aes(x = rate, colour = ptmbin, fill = ptmbin)) +
    # ggplot(aes(x = rate, colour = ptmbin, fill = ptmbin, weight = freq)) +
    geom_density(alpha = 0.3, adjust = 1, linewidth = 0.75/.weight) +
    # annotate(geom = "line", x = c(0.25, 0.25), y = c(0, 2.5), linetype = "dashed", linewidth = 0.5/.weight) +
    # annotate(geom = "line", x = c(0.55, 0.55), y = c(0, 2.5), linetype = "dashed", linewidth = 0.5/.weight) +
    # annotate(geom = "text", x = 0.125, y = 0.15, label = "Buried", size = 5/.pt) +
    # annotate(geom = "text", x = 0.40, y = 0.15, label = "Structured", size = 5/.pt) +
    # annotate(geom = "text", x = 0.775, y = 0.15, label = "Disordered", size = 5/.pt) +
    scale_colour_manual(values = ptmcol, aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = T)) +
    # scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.55, 1), labels = c("0", "0.25", "0.55", "1")) +
    # scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = pretty_breaks(2)) +
    scale_x_continuous(limits = c(2e-2, 9e1), expand = c(0, 0), trans = "log10", breaks = trans_breaks("log10", function(x) 10^x, n=4), labels = c(0.01, 0.1, 1, 10, 100)) +
    scale_y_continuous(expand = c(0, 0.01), breaks = pretty_breaks(1)) +
    facet_grid(. ~ fct_relevel(dis, "Buried", "Structured", "Disordered")) +
    # facet_grid(factor(ptm, levels = c(ptm_order_select, "0")) ~ fct_relevel(dis, "Buried", "Structured", "Disordered")) +
    # facet_grid(factor(ptm, levels = ptm_order_select) ~ fct_relevel(dis, "Buried", "Structured", "Disordered")) +
    theme_minimal() +
    labs(tag = "a") +
    # xlab("Relative accessible surface area") +
    xlab("Evolutionary rate (Rate4Site)") +
    ylab("Probability density") +
    theme_nature() +
    theme(axis.title.y = element_text(hjust = 0))
) %>% qsave("output-density-canonical-evorate.pdf", width = 90, height = 30)
# ) %>% qsave("output-density-canonical-weighted2.pdf", height = 30)


# # Draw simplified evorate plot (not faceting by PTM type)
# # Copied the plotting code from alphasa_relasa_vs_ptms.R because I have weighted controls here (aaweight -> q$freq)
# (
#   q %>%
#     # mutate(type = fct_recode(type, "Modified"="P", "Control"="C")) %>%
#     mutate(dis = fct_recode(dis, "Buried"="strcore", "Structured"="strsurf", "Disordered"="dissurf")) %>%
#     ggplot(aes(rate)) +
#     # geom_density(aes(colour = fct_rev(type), fill = fct_rev(type)), alpha = 0.3, adjust = 1, linewidth = 0.75/.weight) +
#     geom_density(aes(colour = fct_rev(ptmbin), fill = fct_rev(ptmbin)), alpha = 0.3, adjust = 1, linewidth = 0.75/.weight) +
#     # geom_vline(data=mypvalues, aes(xintercept=mediancontrol), color = ptmvir0, linewidth = 0.75/.weight, linetype = "solid") +
#     # geom_vline(data=mypvalues, aes(xintercept=medianptm), color = ptmvir1, linewidth = 0.75/.weight, linetype = "solid") +
#     # geom_text(data=mypvalues, aes(x=0, y=Inf, label=paste("\n", " ", sigadj, sep=""), hjust=0, vjust=0.5), size = 6/.pt) +
#     # geom_text(data=mypvalues, aes(x=Inf, y=Inf, label=paste("\nn=", nptm, " ", sep=""), hjust=1, vjust=0.5), size = 5/.pt) +
#     scale_x_continuous(limits = c(2e-2, 9e1), expand = c(0, 0), trans = "log10", breaks = trans_breaks("log10", function(x) 10^x, n=4), labels = c(0.01, 0.1, 1, 10, 100)) +
#     scale_y_continuous(expand = c(0, 0.2), breaks = pretty_breaks(2)) +
#     scale_colour_manual(values = c("Modified" = ptmvir1, "Control" = ptmvir0), aesthetics = c("colour", "fill"), name = NULL , guide = guide_legend(reverse = T)) +
#     facet_grid(. ~ fct_relevel(dis, "Buried", "Structured", "Disordered")) +
#     xlab("Evolutionary rate (Rate4Site)") +
#     ylab("Probability density") +
#     theme_nature()
# ) %>%
#   qsave(f("output-density-noptmfacet-{source}-{evorate}.pdf"), width = 90, height = 60)


# "Figures are best prepared at a width of 90 mm (single column) and 180 mm (double column) with a maximum height of 170mm. At this size, the font size should be 5-7pt." (https://www.nature.com/nature/for-authors/initial-submission)
# Absolute max Nature dimensions 183 mm x 247 mm ("For guidance, Natures standard figure sizes are 89 mm wide (single column) and 183 mm wide (double column). The full depth of a Nature page is 247 mm. Figures can also be a column-and-a-half where necessary (120–136 mm).", https://www.nature.com/nature/for-authors/final-submission)
