blang_init()

# Get primate percentage etc. for all Ensembl Compara alignments (alns)
q_einsi_tree <- Query("SELECT *, 'einsi_tree' AS mafftmode FROM (SELECT f.aln, COUNT(DISTINCT e.species) AS species, SUM(s.primate) AS primates, COUNT(DISTINCT e.species)-SUM(s.primate) AS nonprimates, SUM(s.primate) / COUNT(DISTINCT e.species) AS primate_percentage, GROUP_CONCAT(s.unispec) FROM comparafasta_einsi_tree f, comparaenspspecies e, comparaspecies s WHERE e.ensp=f.ensp AND e.species=s.species GROUP BY f.aln ORDER BY primate_percentage DESC, species DESC) t")
q_einsi_tree_1para <- Query("SELECT *, 'einsi_tree_1para' AS mafftmode FROM (SELECT f.aln, COUNT(DISTINCT e.species) AS species, SUM(s.primate) AS primates, COUNT(DISTINCT e.species)-SUM(s.primate) AS nonprimates, SUM(s.primate) / COUNT(DISTINCT e.species) AS primate_percentage, GROUP_CONCAT(s.unispec) FROM comparafasta_einsi_tree_1para f, comparaenspspecies e, comparaspecies s WHERE e.ensp=f.ensp AND e.species=s.species GROUP BY f.aln ORDER BY primate_percentage DESC, species DESC) t")
q <- bind_rows(q_einsi_tree, q_einsi_tree_1para)
# q$primate_percentage <- q$primates / q$species

# label_einsi_tree_1para <- '"einsi_tree_1para"\n\nOne-to-many:\nUse best paralog'
# label_einsi_tree_1para <- 'One-to-many:\nuse best paralog\n"einsi_tree_1para"'
# label_einsi_tree_1para <- 'One-to-many:\nUse best paralog'
# label_einsi_tree_1para <- 'One-to-many homologies:\nUse best-available paralog\n(einsi_tree_1para)'
label_einsi_tree_1para <- "Using best available paralogue\nfor one-to-many homologies\n(einsi_tree_1para)"
# label_einsi_tree <- "einsi_tree\n\nOne-to-many:\nSkip\n(One-to-one only)"
# label_einsi_tree <- "einsi_tree\n\nOne-to-many:\nIgnore"
# label_einsi_tree <- '"einsi_tree"\n\nOne-to-one only'
# label_einsi_tree <- 'One-to-one only\n"einsi_tree"'
# label_einsi_tree <- 'One-to-many:\nSkip\n(One-to-one only)'
# label_einsi_tree <- 'One-to-many:\nSkip\n(One-to-one only)'
# label_einsi_tree <- 'One-to-many homologies:\nUse best-available paralog\n(einsi_tree_1para)'
label_einsi_tree <- "Using only one-to-one orthologs\n(einsi_tree)"

 # qplot(qall$species) + xlab("Number of species") + ylab("Number of homology clusters")
# qplot(qall$species) + xlab("Number of species covered") + ylab("One-to-one ortholog clusters")
(
  q %>%
    mutate(mafftmode = ifelse(mafftmode == "einsi_tree_1para", label_einsi_tree_1para, mafftmode)) %>%
    mutate(mafftmode = ifelse(mafftmode == "einsi_tree", label_einsi_tree, mafftmode)) %>%
    ggplot(aes(x = species, fill = mafftmode)) +
    # geom_histogram(binwidth = 5, alpha = 0.9) +
    geom_histogram(binwidth = 5, boundary = 2) +
    # geom_density() +
    # scale_x_continuous(expand = c(0, 0)) +
    # scale_x_continuous(expand = c(0, 0), limits = c(0, 200)) +
    scale_x_continuous(expand = c(0, 0), limits = c(2, 202), breaks = c(2, 50, 100, 150, 200)) +
    scale_y_continuous(expand = c(0, 0)) +
    # scale_colour_manual(aesthetics = c("colour", "fill"), values = c(label_einsi_tree = ptmbd, label_einsi_tree_1para = ptmod)) +
    # scale_colour_manual(aesthetics = c("colour", "fill"), values = c("One-to-many:\nSkip\n(One-to-one only)" = ptmbd, "One-to-many:\nUse best paralog" = ptmod)) +
    scale_colour_manual(aesthetics = c("colour", "fill"), values = c("Using only one-to-one orthologues\n(einsi_tree)" = "grey", "Using best available paralogue\nfor one-to-many homologies\n(einsi_tree_1para)" = "black")) +
    coord_cartesian(xlim = c(2, 200)) +
    # facet_grid(rows = vars(fct_rev(mafftmode))) +
    facet_wrap(facets = vars(mafftmode), ncol = 1) +
    guides(colour = "none", fill = "none", alpha = "none") +
    xlab("Number of species covered") +
    ylab("Retained homology clusters") +
    labs(tag = "d") +
    theme_nature(extra_margin_right = 2)
) %>% qsave("Alignment size.pdf", width = 45, height = 40)
  
# qplot(q$primate_percentage) + scale_x_continuous(labels=percent_format()) + xlab("Primate percentage") + ylab("One-to-one ortholog clusters")
(
  q %>%
    mutate(mafftmode = ifelse(mafftmode == "einsi_tree_1para", label_einsi_tree_1para, mafftmode)) %>%
    mutate(mafftmode = ifelse(mafftmode == "einsi_tree", label_einsi_tree, mafftmode)) %>%
    # ggplot(aes(x = primate_percentage, colour = mafftmode, fill = mafftmode, alpha = 0.2)) + 
    ggplot(aes(x = primate_percentage, fill = mafftmode)) + 
    # geom_histogram() +
    # geom_histogram(binwidth = 0.01) +
    geom_histogram(binwidth = 0.025, boundary = 0) +
    # geom_histogram(binwidth = 0.05) +
    # scale_x_continuous(labels=percent_format()) +
    scale_x_continuous(expand = c(0, 0), labels=percent_format()) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    # scale_colour_manual(aesthetics = c("colour", "fill"), values = c(label_einsi_tree = ptmbd, label_einsi_tree_1para = ptmod)) +
    scale_colour_manual(aesthetics = c("colour", "fill"), values = c("Using only one-to-one orthologues\n(einsi_tree)" = "grey", "Using best available paralogue\nfor one-to-many homologies\n(einsi_tree_1para)" = "black")) +
    # scale_colour_manual(aesthetics = c("colour", "fill"), values = c('One-to-one only\n"einsi_tree"' = ptmbd, 'One-to-many:\nuse best paralog\n"einsi_tree_1para"' = ptmod)) +
    # facet_grid(rows = vars(fct_rev(mafftmode))) +
    # facet_grid(rows = vars(fct_rev(mafftmode)), scales = "free_y") +
    # facet_grid(rows = vars(mafftmode), scales = "free_y") +
    # facet_grid(rows = vars(mafftmode)) +
    # facet_grid(rows = vars(fct_rev(mafftmode))) +
    facet_wrap(facets = vars(mafftmode), ncol = 1) +
    guides(colour = "none", fill = "none", alpha = "none") +
    xlab("Primate percentage") +
    ylab("Filtered homology clusters") +
    labs(tag = "e") +
    theme_nature(extra_margin_right = 2)
) %>% qsave("Primate percentage.pdf", width = 45, height = 40)
