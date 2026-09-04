#!/usr/bin/env Rscript
# plot_surfdist_3d.R
# ------------------
# Density plots for surf_dist_3d_ang from output-surfdist-3d.csv.
# Produces:
#   output-surfdist-3d-density-control_vs_modified.pdf
#   output-surfdist-3d-density-control_vs_modified-by-ptm.pdf
#
# Usage:
#   Rscript plot_surfdist_3d.R
#   Rscript plot_surfdist_3d.R path/to/output-surfdist-3d.csv

source("relevant_files_for_context/blang.R")

# ---------------------------------------------------------------------------
# Read data
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
csv_path <- if (length(args) >= 1) args[1] else "output-surfdist-3d.csv"

d <- read_csv(csv_path, show_col_types = FALSE)
d <- d %>% filter(!is.na(surf_dist_3d_ang))

# Note: Python script already applied pLDDT >= 70, min_pae <= 2, and relasa <= 0.25 (buried).
cat(sprintf("Pre-filter rows: %d\n", nrow(d)))
# d <- d %>% filter(plddt >= 70)
# d <- d %>% filter(min_pae <= 2)
# d <- d %>% filter(dis == "strcore")
cat(sprintf("Post-filter (using Python-prefiltered data) rows: %d\n", nrow(d)))

# Hardcoded list of 17 PTM types observed in output-density-canonical-allptms-weighted.pdf
# (Corresponds to >= 950 sites globally in alphasa R analysis)
keep_ptms <- c(
  "S-p", "T-p", "Y-p", "K-ub", "K-sum", "K-ac", "K-mal", "K-suc", "K-me", 
  "R-me", "N-gly", "T-gly", "S-gly", "M-ox", "C-glt", "C-nit", "C-pal"
)

d <- d %>% filter(ptmbin == "Control" | ptm %in% keep_ptms)
cat(sprintf("Final rows after keeping %d specified PTM types: %d\n", length(keep_ptms), nrow(d)))

out_stem <- tools::file_path_sans_ext(csv_path)

# ---------------------------------------------------------------------------
# Prepare weighting for global density plot
# (Weight Control residues by AA frequencies of Modified residues)
# ---------------------------------------------------------------------------
aaweight <- d %>% 
  filter(ptmbin == "Modified") %>% 
  group_by(aa) %>% 
  tally(name = "n_mod") %>% 
  mutate(freq = n_mod / sum(n_mod)) %>% 
  select(aa, freq)

d <- d %>% left_join(aaweight, by = "aa")
d <- d %>% mutate(freq = ifelse(ptmbin == "Modified", 1, freq))
d <- d %>% mutate(freq = ifelse(is.na(freq), 0, freq))

# Convert to factor
d <- d %>% mutate(ptmbin = factor(ptmbin, levels = c("Modified", "Control")))

# ---------------------------------------------------------------------------
# 1. Overall density plot (Modified vs Control)
# ---------------------------------------------------------------------------
# Note: Wilcoxon test does not use weights trivially, so we run the standard unweighted test
res_wilcox <- wilcox.test(surf_dist_3d_ang ~ ptmbin, data = d)

# p-value formatting helper (reports down to IEEE 754 double precision limit)
fmt_pval <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p == 0) return("p < 2.22e-308")
  if (p < 0.001) {
    s <- sprintf("%.2e", p)
    s <- sub("\\.0+e", "e", s) # turn 1.00e-170 into 1e-170
    s <- sub("\\.(\\d)0e", ".\\1e", s) # turn 1.10e-170 into 1.1e-170
    return(paste0("p = ", s))
  }
  return(sprintf("p = %.3f", p))
}

pval_str <- fmt_pval(res_wilcox$p.value)

# Means SHOULD be weighted to reflect the distributions accurately
means <- d %>% group_by(ptmbin) %>% summarise(mu = weighted.mean(surf_dist_3d_ang, w = freq, na.rm=TRUE))
mu_mod  <- means %>% filter(ptmbin == "Modified") %>% pull(mu)
mu_ctrl <- means %>% filter(ptmbin == "Control") %>% pull(mu)

label_text <- sprintf("µ(Control) = %.2f Å\nµ(Modified) = %.2f Å\n%s", mu_ctrl, mu_mod, pval_str)

(d %>%
  ggplot(aes(x = surf_dist_3d_ang, colour = fct_rev(ptmbin), fill = fct_rev(ptmbin), weight = freq)) +
    geom_density(alpha = 0.3) +
    annotate("text", x = max(d$surf_dist_3d_ang, na.rm = TRUE), y = Inf,
             label = label_text,
             hjust = 1.05, vjust = 1.7, size = 5/.pt, family = "Helvetica Neue") +
    scale_colour_manual(values = ptmcol, aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = TRUE)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = pretty_breaks(3)) +
    coord_cartesian(xlim = c(0, NA)) +
    xlab("Distance to nearest heavy (C/N/O/S) surface atom (Å)") +
    ylab("Probability density") +
    theme_nature(legend_position = "top", extra_margin_right = 10) +
    theme(plot.margin = unit(c(5.5, 15, 5.5, 20), "pt"))
) %>% qsave(sprintf("%s-density-control_vs_modified.pdf", out_stem), width = 85, height = 45)

cat(sprintf("Saved %s-density-control_vs_modified.pdf\n", out_stem))

# ---------------------------------------------------------------------------
# 2. Faceted density by PTM type (Modified vs Control)
# ---------------------------------------------------------------------------
# Order PTMs as specified (descending site count)
ptm_order <- keep_ptms

# Filter for the relevant PTMs
d_modified <- d %>% filter(ptmbin == "Modified", ptm %in% ptm_order)

# For each PTM, include ALL controls of the matching amino acid(s) in its facet
# PTM-to-AA mapping for expanding controls
ptm_aa_map <- d %>% 
  filter(ptmbin == "Modified", ptm %in% ptm_order) %>% 
  group_by(ptm, aa) %>% 
  summarise(.groups = "drop")

if (nrow(ptm_aa_map) == 0) {
  cat("No specified PTM types found in data — skipping faceted density plot.\n")
} else {
  d_control_list <- list()
  for (i in seq_len(nrow(ptm_aa_map))) {
    target_ptm <- ptm_aa_map$ptm[i]
    target_aa  <- ptm_aa_map$aa[i]
    tmp_ctrl   <- d %>% filter(ptmbin == "Control", aa == target_aa)
    if (nrow(tmp_ctrl) > 0) {
      tmp_ctrl$ptm  <- target_ptm  # assign to this facet
      tmp_ctrl$freq <- 1           # 100% of the Control is this AA; weight is 1
      d_control_list[[length(d_control_list) + 1]] <- tmp_ctrl
    }
  }
  d_control_expanded <- bind_rows(d_control_list)
  
  d_facet <- bind_rows(d_modified, d_control_expanded) %>%
    mutate(ptm = factor(ptm, levels = ptm_order))

  # Calculate stats per PTM facet
  facet_stats <- d_facet %>%
    group_by(ptm) %>%
    summarise(
      pval = wilcox.test(surf_dist_3d_ang ~ ptmbin, data = pick(everything()))$p.value,
      mu_mod = weighted.mean(surf_dist_3d_ang[ptmbin == "Modified"], w = freq[ptmbin == "Modified"], na.rm = TRUE),
      mu_ctrl = weighted.mean(surf_dist_3d_ang[ptmbin == "Control"], w = freq[ptmbin == "Control"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      label = sprintf("µ(Ctrl) = %.2f Å\nµ(PTM) = %.2f Å\n%s", mu_ctrl, mu_mod, fmt_pval(pval))
    ) %>%
    ungroup()

  n_ptms <- length(ptm_order)

  (d_facet %>%
    ggplot(aes(x = surf_dist_3d_ang, colour = fct_rev(ptmbin), fill = fct_rev(ptmbin), weight = freq)) +
      geom_density(alpha = 0.3) +
      geom_text(data = facet_stats, aes(label = label), x = Inf, y = Inf,
                inherit.aes = FALSE, hjust = 1.05, vjust = 1.3, 
                size = 5/.pt, family = "Helvetica Neue") +
      scale_colour_manual(values = ptmcol, aesthetics = c("colour", "fill"), name = NULL, guide = guide_legend(reverse = TRUE)) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = pretty_breaks(2)) +
      coord_cartesian(xlim = c(0, NA)) +
      facet_wrap(vars(ptm), ncol = 3, scales = "free_y", dir = "h", axes = "all", axis.labels = "all") +
      xlab("Distance to nearest heavy (C/N/O/S) surface atom (Å)") +
      ylab("Probability density") +
      theme_nature(legend_position = "bottomright", axis_fontsize = 6, extra_margin_right = 10) +
      theme(plot.margin = unit(c(5.5, 15, 5.5, 20), "pt"))
  ) %>% qsave(sprintf("%s-density-control_vs_modified-by-ptm.pdf", out_stem),
              # Width 120 / 3 = 40mm. Height 120 / 6 = 20mm. 40/20 = 2.0 aspect ratio.
              width = 120, height = max(40, ceiling(n_ptms/3) * 20))

  cat(sprintf("Saved %s-density-control_vs_modified-by-ptm.pdf (%d PTM types)\n", out_stem, n_ptms))
}

cat("Done!\n")
