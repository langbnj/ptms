source("~/Documents/R/blang.R")
blang_init()
setwd("~/Library/CloudStorage/Dropbox/Madan/SJ-Primary-Projects/PTMs/ptm_contact_aas")

# Get PTM types with at least 1000 sites
ptms1000 <- c("R-me", "K-sum", "S-gly", "K-me", "T-gly", "S-p", "K-ac", "T-p", "K-ub", "K-mal", "K-suc", "N-gly", "Y-p", "M-ox", "C-pal", "C-glt", "C-nit") # Order as in alphasa "all ptms" figure (all >1000, ptms1000 I suppose)
ptms11 <- c("S-p", "T-p", "Y-p", "K-ub", "K-sum", "K-ac", "K-mal", "K-suc", "K-me", "R-me", "N-gly") # Order as in Figure 3



plot_contacts_bars <- function(mytag, tmpptm, my_width, my_height, minsites, mytitle, myfilename, source = "all", disfilt = "", show_insig = F) {
  
  # # Load data
  # system.time(q <- read_tsv("ptm_contact_aas.txt.gz"))
  system.time(q <- read_tsv(f("output-contacts-human-{source}-1000.txt.gz"), col_types = "ccciicccccddcccccdd"))
  
  # PAE as integer, since it's always an integer anyway (even though it's printed as e.g. 2.0 in the file)
  q %<>% mutate(pae = as.integer(pae))
  
  if (str_detect(myfilename, "pLDDT70")) {
    q %<>% filter(plddt >= 70)
  }
  if (str_detect(myfilename, "pLDDT90")) {
    q %<>% filter(plddt >= 90)
  }
  if (str_detect(myfilename, "PAE2")) {
    # Filter out fragmented structures (these will have pae 'None')
    q %<>% filter(!is.na(pae))
    q %<>% filter(pae <= 2)
  }
  if (str_detect(myfilename, "PAE1")) {
    # Filter out fragmented structures (these will have pae 'None')
    q %<>% filter(!is.na(pae))
    q %<>% filter(pae <= 1)
  }
  
  q$ptm <- factor(q$ptm, levels = tmpptm)
  if (str_detect(myfilename, "STY-p")) {
    # Phosphorylation only
    q %<>% filter(ptm %in% c("S-p", "T-p", "Y-p"))
    q$ptm %<>% droplevels
  }
  if (str_detect(myfilename, "K-ub-sum-mal")) {
    # K-ub, K-sum and K-mal only (these show interesting trends)
    q %<>% filter(ptm %in% c("K-ub", "K-sum", "K-mal"))
    q$ptm %<>% droplevels
  }
  
  # # Source (all, UniProt etc.)
  # q %<>% filter(source == source)
  
  # Rename some columns etc. to align with other analysis scripts
  q %<>% mutate(disstr = dis, coresurf = surf, dis = class, class = NULL, surf = NULL)
  q %<>% mutate(type = ifelse(type == "mod", "ptm", type))
  q
  q %>% group_by(dis) %>% tally
  q %<>% mutate(dis = ifelse(dis == "core", "Buried", dis))
  q %<>% mutate(dis = ifelse(dis == "strsurf", "Structured", dis))
  q %<>% mutate(dis = ifelse(dis == "dissurf", "Disordered", dis))
  q %>% group_by(dis) %>% tally
  q %<>% mutate(dis = factor(dis, levels = c("Buried", "Structured", "Disordered")))
  q
  q %>% group_by(dis) %>% tally
  
  q %>%
    # filter(ptm %in% ptms1000) %>%
    # Unique so as not to double-count e.g. a Y with two D contacts, and different atoms etc.
    select(ptm, dis, type, acc, site, aa1, aa2) %>% unique %>%
    group_by(ptm, dis, type, aa1, aa2) %>% 
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
    # filter(fisher_fdr < 0.05) %>% # Significant only (FDR)
    filter(ptm %in% tmpptm) %>% # PTMs being plotted only
    arrange(fisher_fdr) -> 
    qfilt
  if (show_insig == F) {
    qfilt %<>% filter(fisher_fdr < 0.05)
  }
  qfilt
  # {
  #   # q %>% filter(source == "all") %>% filter(ptm %in% c("S-p", "T-p", "Y-p")) %>% group_by(ptm, class, type) %>% count(aa2) %>% pivot_wider(names_from = type, values_from = n) %>% group_by(ptm, class) %>% mutate(controlfreq = control / sum(control), controlsum = sum(control), modfreq = mod / sum(mod), modsum = sum(mod), enrich = modfreq / controlfreq, controlfreqsum = sum(controlfreq), modfreqsum = sum(modfreq), enrichsum = sum(enrich)) %>%
  #   # q %>% group_by(ptm, dis, type) %>% count(aa2) %>% pivot_wider(names_from = type, values_from = n) %>% group_by(ptm, class) %>% mutate(controlfreq = control / sum(control), controlsum = sum(control), modfreq = mod / sum(mod), modsum = sum(mod), enrich = modfreq / controlfreq, controlfreqsum = sum(controlfreq), modfreqsum = sum(modfreq), enrichsum = sum(enrich)) %>%
  #   qproc  %>%
  #     ggplot(aes(x = aa2, y = log2enrich, fill = factor(log2enrich >= 0))) +
  #     geom_bar(stat = "identity") +
  #     # scale_fill_viridis_d() +
  #     facet_grid(rows = vars(paste(ptm, dis))) +
  #     # facet_grid(cols = vars(dis), rows = vars(ptm)) +
  #     theme_nature() +
  #     theme(strip.text.y = element_text(angle = 0, hjust = 0)) +
  #     guides(fill = "none") +
  #     xlab("Contacted amino acid") +
  #     ylab("log2 Enrichment (PTM sites relative to control residues)")
  # } %>% qsave(f("output-ptm_contact_aas-new-enrichment-STY-p-only.pdf"), width=80, height=80)
  
  # disfilt
  tmpdisfilt <- ""
  if (disfilt == "nodis") {
    # qfilt %<>% filter(dis != "dissurf")
    qfilt %<>% filter(dis != "Disordered")
    qfilt$dis %<>% droplevels
    qfilt$dis
    tmpdisfilt <- "-nodis"
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
  my_legend_nudge_top <- -4
  
  print(qfilt)
  print(f("{qfilt %>% nrow} rows"))
  
  # Only plot if any enrichments are available (otherwise plotting would crash)  
  if (qfilt %>% nrow > 0) {
    {
      qfilt %>%
        # mutate(enriched = ifelse(enrichment > 0, "enriched", "depleted")) %>%
        ggplot(aes(x = aa2, y = log2enrich)) +
        # geom_col(aes(fill = enriched), position = "dodge") +
        geom_col(aes(fill = log2enrich), position = "dodge") +
        # geom_point() +
        # geom_text(aes(label = aa2), size = 5/.pt, nudge_y = 1) +
        # geom_text(aes(label = aa2, vjust = ifelse(log2enrich >= 0, -0.5, 1.5)), size = 5/.pt) +
        geom_text(aes(label = aa2, vjust = ifelse(log2enrich >= 0, -0.1, 1.1)), size = 5/.pt) +
        # geom_text(aes(label = aa2, vjust = ifelse(log2enrich >= 0, 0, 1)), size = 5/.pt) +
        geom_hline(yintercept = 0, linetype = "solid") +
        # coord_cartesian(ylim = c(-1, 1)) +
        # scale_y_continuous(expand = c(0.4, 0.4)) + # Expand more to make sure letters get plotted fully
        scale_y_continuous(expand = c(0.6, 0.6)) + # Expand more to make sure letters get plotted fully
        # facet_wrap(~ ptm, scales = "free") +
        # facet_grid(cols = vars(ptm), rows = vars(fct_rev(type)), scales = "free") +
        # facet_grid(rows = vars(as.factor(ptm)), cols = vars(as.factor(dis)), scales = "free", space = "free_x") +
        # facet_grid(rows = vars(factor(ptm, levels = tmpptm)), cols = vars(as.factor(dis)), scales = "free", space = "free_x") +
        # facet_grid(rows = vars(factor(ptm, levels = tmpptm)), cols = vars(as.factor(dis)), drop = T, scales = "free", space = "free_x", switch = "y") +
        facet_grid(rows = vars(ptm), cols = vars(as.factor(dis)), drop = F, scales = "free", space = "free_x", switch = "y") +
        # scale_fill_viridis_d() +
        # scale_fill_manual(values = c("ptm" = ptmvir1, "control" = ptmvir0)) +
        # scale_fill_manual(values = c("enriched" = scales::viridis_pal()(5)[4], "depleted" = scales::viridis_pal()(5)[2])) +
        # scale_fill_viridis_c() +
        # scale_fill_gradientn(name = "Variable residues", colours = c(scales::viridis_pal()(5)[1], scales::viridis_pal()(5)[2], "#FFFFFFFF", scales::viridis_pal()(5)[4], scales::viridis_pal()(5)[5]), values = scales::rescale(c(-maxenrich, -minenrich, 0, minenrich, maxenrich)), na.value = "white", breaks = c(-mybreak, 0, mybreak), labels = c(f("-{mybreak}"), "0", f("{mybreak}")), limits = c(-maxenrich, maxenrich)) +
        scale_fill_gradientn(name = "Contact enrichment (log2)", colours = c(scales::viridis_pal()(5)[1], scales::viridis_pal()(5)[2], "#FFFFFFFF", scales::viridis_pal()(5)[4], scales::viridis_pal()(5)[5]), values = scales::rescale(c(-maxenrich, -minenrich, 0, minenrich, maxenrich)), na.value = "white", breaks = c(-mybreak, 0, mybreak), labels = c(f("-{mybreak}"), "0", f("{mybreak}")), limits = c(-maxenrich, maxenrich)) +
        theme_nature(legend_position = "bottom", legend_nudge_top=my_legend_nudge_top, legend_nudge_right=0, extra_margin_right=2) +
        theme(legend.text = element_text(size = 5, margin = margin(t = 1))) +
        # theme(legend.text = element_text(size = 6, margin = margin(t = 1))) +
        theme(axis.text.y = element_text(size = 6)) +
        # theme(legend.title.position = "top", legend.title = element_text(size = 5, hjust = 1, vjust = 1, margin = margin(b = unit(0.75, "mm"), t = unit(1.5, "mm"), r = unit(1, "mm")))) +
        theme(legend.title.position = "left", legend.title = element_text(size = 5, vjust = 1, margin = margin(t = unit(0.75, "mm"), r = unit(1, "mm")))) +
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
      # } %>% qsave(f("output-ptm_contact_aas-new-enrichment-STY-p-only-pLDDT70_PAE2.pdf"), width=my_width, height=my_height)
    } %>% qsave(f("output-ptm_contact_aas-bars-filtered-fdr-{myfilename}-{source}-minsites{minsites}{tmpdisfilt}.pdf"), width = my_width, height = my_height)
    
    # Save qfilt to TSV
    write_tsv(qfilt, f("output-ptm_contact_aas-bars-tsv-filtered-fdr-{myfilename}-{source}-minsites{minsites}{tmpdisfilt}.tsv"))
  }
}

# # Load data
# system.time(q <- read_tsv("output-contacts-human-all-1000.txt.gz"))
# q %<>% filter(plddt >= 70)
# q %<>% filter(pae <= 2)
# q %<>% filter(plddt >= 90)
# q %<>% filter(pae <= 1)
# # # q %<>% filter(class != 'dissurf')
# # q %<>% filter(ptm %in% c("S-p", "T-p", "Y-p"))
# # Rename some columns etc. to align with other analysis scripts
# q %<>% mutate(disstr = dis, coresurf = surf, dis = class, class = NULL, surf = NULL)
# q %<>% mutate(type = ifelse(type == "mod", "ptm", type))




# # minsites 1000
# 
# # All sources, unfiltered
# plot_contacts_bars("a", ptms11, 80, 80, 1000, NULL, "a-all", source = "all", show_insig = F)
# plot_contacts_bars("a", ptms11, 80, 80, 1000, NULL, "a-all", source = "all", disfilt = "nodis", show_insig = F)
# # Different pLDDT & PAE filtering
# plot_contacts_bars("b", ptms11, 80, 80, 1000, NULL, "b-pLDDT70_PAE2", source = "all", show_insig = F)
# plot_contacts_bars("c", ptms11, 80, 80, 1000, NULL, "c-pLDDT70_PAE2", source = "all", disfilt = "nodis", show_insig = F)
# # plot_contacts_bars("d", ptms11, 80, 80, 1000, NULL, "d-pLDDT90_PAE2", source = "all", show_insig = F)
# # plot_contacts_bars("e", ptms11, 80, 80, 1000, NULL, "e-pLDDT90_PAE2", source = "all", disfilt = "nodis", show_insig = F)
# plot_contacts_bars("f", ptms11, 80, 80, 1000, NULL, "f-pLDDT90_PAE1", source = "all", show_insig = F)
# plot_contacts_bars("g", ptms11, 80, 80, 1000, NULL, "g-pLDDT90_PAE1", source = "all", disfilt = "nodis", show_insig = F)
# # Different sources, pLDDT70_PAE2
# # plot_contacts_bars("h", ptms11, 80, 80, 1000, NULL, "h-pLDDT70_PAE2", source = "uniprot", disfilt = "nodis", show_insig = F)
# # plot_contacts_bars("i", ptms11, 80, 80, 1000, NULL, "i-pLDDT70_PAE2", source = "ochoa", disfilt = "nodis", show_insig = F)
# plot_contacts_bars("j", ptms11, 80, 80, 1000, NULL, "j-pLDDT70_PAE2", source = "phosphositeplus", disfilt = "nodis", show_insig = F)
# plot_contacts_bars("k", ptms11, 80, 80, 1000, NULL, "k-pLDDT70_PAE2", source = "dbptm", disfilt = "nodis", show_insig = F)
# # Phosphorylation only, different pLDDT & PAE filtering
# plot_contacts_bars("l", ptms11, 80, 36, 1000, NULL, "l-STY-p", source = "all", show_insig = F)
# plot_contacts_bars("m", ptms11, 80, 36, 1000, NULL, "m-STY-p", source = "all", disfilt = "nodis", show_insig = F)
# plot_contacts_bars("n", ptms11, 80, 36, 1000, NULL, "n-STY-p-pLDDT70_PAE2", source = "all", show_insig = F)
# plot_contacts_bars("o", ptms11, 80, 36, 1000, NULL, "o-STY-p-pLDDT70_PAE2", source = "all", disfilt = "nodis", show_insig = F)
# plot_contacts_bars("p", ptms11, 80, 36, 1000, NULL, "p-STY-p-pLDDT90_PAE1", source = "all", show_insig = F)
# plot_contacts_bars("q", ptms11, 80, 36, 1000, NULL, "q-STY-p-pLDDT90_PAE1", source = "all", disfilt = "nodis", show_insig = F)
# 
# 
# 
# 
# # minsites 100
# 
# # All sources, unfiltered
# plot_contacts_bars("a", ptms11, 80, 80, 100, NULL, "a-all", source = "all", show_insig = F)
# plot_contacts_bars("a", ptms11, 80, 80, 100, NULL, "a-all", source = "all", disfilt = "nodis", show_insig = F)
# # Different pLDDT & PAE filtering
# plot_contacts_bars("b", ptms11, 80, 80, 100, NULL, "b-pLDDT70_PAE2", source = "all", show_insig = F)
# plot_contacts_bars("c", ptms11, 80, 80, 100, NULL, "c-pLDDT70_PAE2", source = "all", disfilt = "nodis", show_insig = F)
# # plot_contacts_bars("d", ptms11, 80, 80, 100, NULL, "d-pLDDT90_PAE2", source = "all", show_insig = F)
# # plot_contacts_bars("e", ptms11, 80, 80, 100, NULL, "e-pLDDT90_PAE2", source = "all", disfilt = "nodis", show_insig = F)
# plot_contacts_bars("f", ptms11, 80, 80, 100, NULL, "f-pLDDT90_PAE1", source = "all", show_insig = F)
# plot_contacts_bars("g", ptms11, 80, 80, 100, NULL, "g-pLDDT90_PAE1", source = "all", disfilt = "nodis", show_insig = F)
# # Different sources, pLDDT70_PAE2
# plot_contacts_bars("h", ptms11, 80, 80, 100, NULL, "h-pLDDT70_PAE2", source = "uniprot", disfilt = "nodis", show_insig = F)
# plot_contacts_bars("i", ptms11, 80, 80, 100, NULL, "i-pLDDT70_PAE2", source = "ochoa", disfilt = "nodis", show_insig = F)
# plot_contacts_bars("j", ptms11, 80, 80, 100, NULL, "j-pLDDT70_PAE2", source = "phosphositeplus", disfilt = "nodis", show_insig = F)
# plot_contacts_bars("k", ptms11, 80, 80, 100, NULL, "k-pLDDT70_PAE2", source = "dbptm", disfilt = "nodis", show_insig = F)
# # Phosphorylation only, different pLDDT & PAE filtering
# plot_contacts_bars("l", ptms11, 80, 36, 100, NULL, "l-STY-p", source = "all", show_insig = F)
# plot_contacts_bars("m", ptms11, 80, 36, 100, NULL, "m-STY-p", source = "all", disfilt = "nodis", show_insig = F)
# plot_contacts_bars("n", ptms11, 80, 36, 100, NULL, "n-STY-p-pLDDT70_PAE2", source = "all", show_insig = F)
# plot_contacts_bars("o", ptms11, 80, 36, 100, NULL, "o-STY-p-pLDDT70_PAE2", source = "all", disfilt = "nodis", show_insig = F)
# plot_contacts_bars("p", ptms11, 80, 36, 100, NULL, "p-STY-p-pLDDT90_PAE1", source = "all", show_insig = F)
# plot_contacts_bars("q", ptms11, 80, 36, 100, NULL, "q-STY-p-pLDDT90_PAE1", source = "all", disfilt = "nodis", show_insig = F)



# minsites 40

# All sources, unfiltered
plot_contacts_bars("a", ptms11, 80, 80, 40, NULL, "a-all", source = "all", show_insig = F)
plot_contacts_bars("a", ptms11, 80, 80, 40, NULL, "a-all", source = "all", disfilt = "nodis", show_insig = F)
# Different pLDDT & PAE filtering
plot_contacts_bars("b", ptms11, 80, 80, 40, NULL, "b-pLDDT70_PAE2", source = "all", show_insig = F)
plot_contacts_bars("c", ptms11, 80, 80, 40, NULL, "c-pLDDT70_PAE2", source = "all", disfilt = "nodis", show_insig = F)
# plot_contacts_bars("d", ptms11, 80, 80, 40, NULL, "d-pLDDT90_PAE2", source = "all", show_insig = F)
# plot_contacts_bars("e", ptms11, 80, 80, 40, NULL, "e-pLDDT90_PAE2", source = "all", disfilt = "nodis", show_insig = F)
plot_contacts_bars("f", ptms11, 80, 80, 40, NULL, "f-pLDDT90_PAE1", source = "all", show_insig = F)
plot_contacts_bars("g", ptms11, 80, 80, 40, NULL, "g-pLDDT90_PAE1", source = "all", disfilt = "nodis", show_insig = F)
# Different sources, pLDDT70_PAE2
plot_contacts_bars("h", ptms11, 80, 80, 40, NULL, "h-pLDDT70_PAE2", source = "uniprot", disfilt = "nodis", show_insig = F)
plot_contacts_bars("i", ptms11, 80, 80, 40, NULL, "i-pLDDT70_PAE2", source = "ochoa", disfilt = "nodis", show_insig = F)
plot_contacts_bars("j", ptms11, 80, 80, 40, NULL, "j-pLDDT70_PAE2", source = "phosphositeplus", disfilt = "nodis", show_insig = F)
plot_contacts_bars("k", ptms11, 80, 80, 40, NULL, "k-pLDDT70_PAE2", source = "dbptm", disfilt = "nodis", show_insig = F)
# Phosphorylation only, different pLDDT & PAE filtering
plot_contacts_bars("l", ptms11, 80, 36, 40, NULL, "l-STY-p", source = "all", show_insig = F)
plot_contacts_bars("m", ptms11, 80, 36, 40, NULL, "m-STY-p", source = "all", disfilt = "nodis", show_insig = F)
plot_contacts_bars("n", ptms11, 80, 36, 40, NULL, "n-STY-p-pLDDT70_PAE2", source = "all", show_insig = F)
plot_contacts_bars("o", ptms11, 80, 36, 40, NULL, "o-STY-p-pLDDT70_PAE2", source = "all", disfilt = "nodis", show_insig = F)
plot_contacts_bars("p", ptms11, 80, 36, 40, NULL, "p-STY-p-pLDDT90_PAE1", source = "all", show_insig = F)
plot_contacts_bars("q", ptms11, 80, 36, 40, NULL, "q-STY-p-pLDDT90_PAE1", source = "all", disfilt = "nodis", show_insig = F)
# labelled for Fig. 3f,g
plot_contacts_bars("f", ptms11, 105, 45, 10, "Contacts", "3f-STY-p-pLDDT70_PAE2", source = "all", disfilt = "nodis", show_insig = F)
plot_contacts_bars("g", ptms11, 105, 45, 10, "Contacts", "3g-K-ub-sum-mal-pLDDT70_PAE2", source = "all", disfilt = "nodis", show_insig = F)



# # minsites 10
# 
# # All sources, unfiltered
# plot_contacts_bars("a", ptms11, 80, 80, 10, NULL, "a-all", source = "all", show_insig = F)
# plot_contacts_bars("a", ptms11, 80, 80, 10, NULL, "a-all", source = "all", disfilt = "nodis", show_insig = F)
# # Different pLDDT & PAE filtering
# plot_contacts_bars("b", ptms11, 80, 80, 10, NULL, "b-pLDDT70_PAE2", source = "all", show_insig = F)
# plot_contacts_bars("c", ptms11, 80, 80, 10, NULL, "c-pLDDT70_PAE2", source = "all", disfilt = "nodis", show_insig = F)
# # plot_contacts_bars("d", ptms11, 80, 80, 10, NULL, "d-pLDDT90_PAE2", source = "all", show_insig = F)
# # plot_contacts_bars("e", ptms11, 80, 80, 10, NULL, "e-pLDDT90_PAE2", source = "all", disfilt = "nodis", show_insig = F)
# plot_contacts_bars("f", ptms11, 80, 80, 10, NULL, "f-pLDDT90_PAE1", source = "all", show_insig = F)
# plot_contacts_bars("g", ptms11, 80, 80, 10, NULL, "g-pLDDT90_PAE1", source = "all", disfilt = "nodis", show_insig = F)
# # Different sources, pLDDT70_PAE2
# plot_contacts_bars("h", ptms11, 80, 80, 10, NULL, "h-pLDDT70_PAE2", source = "uniprot", disfilt = "nodis", show_insig = F)
# plot_contacts_bars("i", ptms11, 80, 80, 10, NULL, "i-pLDDT70_PAE2", source = "ochoa", disfilt = "nodis", show_insig = F)
# plot_contacts_bars("j", ptms11, 80, 80, 10, NULL, "j-pLDDT70_PAE2", source = "phosphositeplus", disfilt = "nodis", show_insig = F)
# plot_contacts_bars("k", ptms11, 80, 80, 10, NULL, "k-pLDDT70_PAE2", source = "dbptm", disfilt = "nodis", show_insig = F)
# # Phosphorylation only, different pLDDT & PAE filtering
# plot_contacts_bars("l", ptms11, 80, 36, 10, NULL, "l-STY-p", source = "all", show_insig = F)
# plot_contacts_bars("m", ptms11, 80, 36, 10, NULL, "m-STY-p", source = "all", disfilt = "nodis", show_insig = F)
# plot_contacts_bars("n", ptms11, 80, 36, 10, NULL, "n-STY-p-pLDDT70_PAE2", source = "all", show_insig = F)
# plot_contacts_bars("o", ptms11, 80, 36, 10, NULL, "o-STY-p-pLDDT70_PAE2", source = "all", disfilt = "nodis", show_insig = F)
# plot_contacts_bars("p", ptms11, 80, 36, 10, NULL, "p-STY-p-pLDDT90_PAE1", source = "all", show_insig = F)
# plot_contacts_bars("q", ptms11, 80, 36, 10, NULL, "q-STY-p-pLDDT90_PAE1", source = "all", disfilt = "nodis", show_insig = F)
