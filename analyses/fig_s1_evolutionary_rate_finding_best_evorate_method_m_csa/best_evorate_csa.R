blang_init()

library("coin")
library("ggbeeswarm")

# Read small table showing normalised median delta and median p-values
q <- read_tsv("output-table.tsv")
q

# Read raw table with individual residue rates (CSA vs. control)
qr <- read_tsv("tmp-t.txt")

(
  q %>%
    # Show non-normalised only
    filter(!str_detect(evorate, "^norm_")) %>%
    # ggplot(aes(x = -median_normalised_rank, y = -log10(median_pvalue))) +
    # ggplot(aes(x = median_normalised_rank, y = median_pvalue)) +
    # ggplot(aes(x = median_normalised_rank, y = median_pvalue, colour = evorate)) +
    ggplot(aes(x = median_normalised_rank, y = median_pvalue, colour = evorate)) +
    # ggplot(aes(x = -median_normalised_rank, y = -log10(median_pvalue), colour = evorate)) +
    # geom_point(size = 1) +
    # geom_point(size = 0.5) +
    geom_point() +
    # geom_text(data = subset(q, evorate == "rate4site_einsi_tree_1para"), aes(label = evorate), size = 5/.pt, hjust = "inward", nudge_x = -0.003) +
    # geom_text(data = subset(q, median_pvalue < 0.02), aes(label = evorate), size = 5/.pt, hjust = "inward", nudge_x = -0.003) +
    # geom_text_repel(data = subset(q, evorate == "rate4site_einsi_tree_1para"), aes(label = evorate), size = 5/.pt, min.segment.length = 0.01) +
    # geom_text_repel(data = subset(q, evorate == "rate4site_einsi_tree_1para"), aes(label = evorate), size = 5/.pt) +
    # geom_text_repel(aes(label = evorate)) +
    # geom_text_repel(aes(label = evorate), max.overlaps = 1, max.iter = 1000000, size = 5/.pt) +
    # geom_text_repel(aes(label = evorate), max.iter = 1000000, size = 5/.pt) +
    # geom_text_repel(aes(label = evorate), max.iter = 10000000, size = 5/.pt, min.segment.length = 0.02, segment.size = 0.5/.weight, point.padding = 0.1) +
    # geom_text_repel(aes(label = evorate)) +
    # geom_text_repel(aes(label = evorate), max.iter = 10000000, size = 5/.pt, min.segment.length = 0.02, segment.size = 0.5/.weight, point.padding = 0.1) +
    geom_text_repel(aes(label = evorate), max.iter = 10000000, min.segment.length = 0.02, point.padding = 0.1) +
    # geom_text_repel(aes(label = evorate, fontface = ifelse(evorate == "rate4site_einsi_tree_1para", "bold", "plain")), max.iter = 10000000, size = 5/.pt, min.segment.length = 0.02, segment.size = 0.25) +
    # geom_label(aes(label = evorate)) +
    # geom_label_repel(aes(label = evorate)) +
    # geom_text_repel(aes(label = ifelse(median_pvalue < 0.02, evorate, NA))) +
    # geom_text_repel(data = subset(q, median_pvalue < 0.02), aes(label = evorate)) +
    # scale_colour_viridis_d() +
    # scale_colour_viridis_d(option = "cividis") +
    # scale_colour_viridis_d() +
    # scale_colour_manual(values = c("rate4site_einsi_tree_1para" = ptmo), na.value = ptmg) +
    # scale_colour_manual(values = c("rate4site_einsi_tree_1para" = ptmod), na.value = ptmg) +
    # scale_colour_manual(values = c("rate4site_einsi_tree_1para" = ptmo), na.value = ptmbl) +
    # scale_colour_manual(values = c("rate4site_einsi_tree_1para" = ptmod), na.value = ptmbl) +
    # scale_colour_manual(values = c("rate4site_einsi_tree_1para" = ptmb), na.value = ptmg) +
    # scale_colour_manual(values = c("rate4site_einsi_tree_1para" = ptmod), na.value = ptmg) +
    # scale_colour_manual(values = c("rate4site_einsi_tree_1para" = "black"), na.value = "#bec1c0") +
    scale_colour_manual(values = c("rate4site_einsi_tree_1para" = "black"), na.value = "#bec1c0") +
    # scale_colour_manual(values = c("rate4site_einsi_tree_1para" = ptmvir1), na.value = ptmvir0) +
    # coord_cartesian(ylim = c(min(q$median_pvalue), 0.05)) +
    # coord_cartesian(ylim = c(0.005, 0.05)) +
    # coord_cartesian(ylim = c(0, 0.05)) +
    # coord_cartesian(ylim = c(0.05, 0.005)) +
    coord_cartesian(ylim = c(0.05, min(q$median_pvalue))) +
    # coord_cartesian(xlim = c(max(q$median_normalised_rank), min(q$median_normalised_rank) - 0.005), ylim = c(0.05, min(q$median_pvalue))) +
    # coord_fixed() +
    scale_x_reverse() +
    scale_y_reverse() +
    guides(colour = "none") +
    # xlab("Median difference in normalized residue rank, as compared to structured control residues") +
    # ylab("Median Mann-Whitney U test p-value (per protein and catalytic residue, against control)") +
    labs(tag = "a") +
    ggtitle("Conservation score selection") +
    xlab("Median normalized rank difference") +
    ylab("Median p-value") +
    # ggtitle("Conservation score evaluation") +
    theme_nature(extra_margin_right = 2)
    # theme(plot.tag.position = c(0, 1))
) %>%
# dev.off()
# ggsave(filename = "output-best-evorate.pdf", width = 40, height = 40, units = "mm")
# ggsave(filename = "output-best-evorate.pdf", dpi = "retina", width = 40, height = 40, units = "mm")
# ggsave(filename = "output-best-evorate.pdf", dpi = "retina", width = 40, height = 40, units = "mm", family = "Helvetica Neue")
# ggsave(filename = "output-best-evorate.pdf", dpi = "retina", width = 40, height = 40, units = "mm", family = "Helvetica Neue", device = cairo_pdf)
# ggsave(filename = "output-best-evorate.pdf", dpi = "retina", width = 40, height = 40, units = "mm", family = "Helvetica Neue", device = quartz)
# quartz(title = "ggplot", file = "output-best-evorate.pdf", type = "pdf", dpi = 320, width = conv_unit(40, "mm", "inch"), height = conv_unit(40, "mm", "inch"), family = "Helvetica Neue")
# ggsave(title = "ggplot", file = "output-best-evorate.pdf", type = "pdf", dpi = "retina", width = 40, height = 40, units = "mm", family = "Helvetica Neue", device = quartz)
# ggsave(title = "ggplot", file = "output-best-evorate.pdf", type = "pdf", dpi = "retina", width = conv_unit(40, "mm", "inch"), height = conv_unit(40, "mm", "inch"), family = "Helvetica Neue", device = quartz)
# ggsave(title = "ggplot", file = "output-best-evorate.pdf", type = "pdf", dpi = 600, width = conv_unit(40, "mm", "inch"), height = conv_unit(40, "mm", "inch"), family = "Helvetica Neue", device = quartz)
  qsave("output-best-evorate.pdf")
  # qsave("output-best-evorate.pdf", width = 39, height = 39)

# Version showing all 3 einsi_tree_1para methods in black
(
  q %>%
    filter(!str_detect(evorate, "^norm_")) %>%
    ggplot(aes(x = median_normalised_rank, y = median_pvalue, colour = evorate)) +
    geom_point() +
    geom_text_repel(aes(label = evorate), max.iter = 10000000, min.segment.length = 0.02, point.padding = 0.1) +
    scale_colour_manual(values = c("rate4site_einsi_tree_1para" = "black", "lichtarge_einsi_tree_1para" = "black", "capra0_einsi_tree_1para" = "black"), na.value = "#bec1c0") +
    coord_cartesian(ylim = c(0.05, min(q$median_pvalue))) +
    scale_x_reverse() +
    scale_y_reverse() +
    guides(colour = "none") +
    labs(tag = "a") +
    ggtitle("Conservation score selection") +
    xlab("Median normalized rank difference") +
    ylab("Median p-value") +
    theme_nature(extra_margin_right = 2)
) %>%
  qsave("output-best-evorate-3.pdf")




# Analyse raw table with individual residue rates (CSA vs. control)
qr
qr$...1 <- NULL
qr %<>% unique
qr
qr %>%
  ggplot(aes(x = rate, colour = type)) +
  geom_density() +
  scale_x_log10() +
  facet_wrap(facets = vars(evorate), ncol = 3) +
  theme_nature()
qsave("output-evorate-distributions.pdf", width = 3 * 40, height = 7 * 20)
# Try stratified wilcox_test
qr %<>% mutate(evorate = as.factor(evorate))
qr %<>% mutate(type = factor(type, levels=c("csa", "control")))
qr %<>% mutate(acc = as.factor(acc))
qr %<>% mutate(aa = as.factor(aa))
levels(qr$type)
qr

# Visualising the number of CSA and control residues per acc:
qr %>% 
  group_by(evorate, acc, type) %>% 
  tally %>%
  pivot_wider(names_from = "type", values_from = "n") %>%
  filter(is.na(control)) %>% pull(acc) %>% unique %>% length
# >> Only one acc, P00488, has no control residues. That's few enough to ignore.
qr %>% pull(acc) %>% unique %>% length
# >> In total, there are 102 accs.

# Plot: How many control residues are there vs. CSA residues?
(qr %>% 
  filter(evorate == "rate4site_einsi_tree_1para") %>%
  group_by(evorate, acc, type) %>% 
  tally %>%
  pivot_wider(names_from = "type", values_from = "n") %>%
  mutate(control = replace_na(control, 0)) %>%
  ggplot(aes(x = csa, y = control)) +
  # geom_hex(bins = 10) +
  # geom_beeswarm() +
  # geom_point(position = "jitter") +
  geom_point() +
  # geom_hline(yintercept = 5, linetype = "dashed") +
  annotate(geom = "line", x = c(1, 14), y = c(5, 5), linetype = "dashed", linewidth = 0.5/.weight) +
  annotate(geom = "text", x = Inf, y = 5, label = "5", hjust = 1, size = 5/.pt) +
  stat_smooth(linewidth = 1/.pt, method = "lm", se = F) +
  coord_cartesian(ylim = c(0, 300)) +
  # facet_wrap(facets = vars(evorate), ncol = 3) +
  labs(tag = "b") +
  ggtitle("Number of catalytic residues\ncompared to control residues") +
  xlab("Catalytic residues (M-CSA)") +
  ylab("Control residues") +
  theme_nature()
) %>% qsave("output-control-csa-residue-count.pdf")
# qsave("output-control-csa-residue-count.pdf", width = 50, height = 50)

# How many accs have catalytic sites that are less conserved than control?
qr %>% 
  group_by(evorate, acc, type) %>% 
  summarise(median_rate = median(rate)) %>%
  pivot_wider(names_from = "type", values_from = "median_rate") %>%
  filter(!is.na(control)) %>%
  group_by(evorate, csa < control) %>%
  tally %>% 
  pivot_wider(names_from = `csa < control`, values_from = n) %>%
  arrange(desc(`FALSE`)) %>%
  print(n = 50)
# capra0_einsi_tree_1para has only 2 proteins where the median evorate for CSA residues is ≥ control.
# rate4site_einsi_tree_1para has only 3 proteins where the median evorate for CSA residues is ≥ control.
# lichtarge_einsi_tree_1para has only 3 proteins where the median evorate for CSA residues is ≥ control.
# >> i.e. rate4site_einsi_tree_1para looks good.
# Using the mean instead of median:
qr %>% 
  group_by(evorate, acc, type) %>% 
  summarise(mean_rate = mean(rate)) %>%
  pivot_wider(names_from = "type", values_from = "mean_rate") %>%
  filter(!is.na(control)) %>%
  group_by(evorate, csa < control) %>%
  tally %>% 
  pivot_wider(names_from = `csa < control`, values_from = n) %>%
  arrange(desc(`FALSE`)) %>%
  print(n = 50)
# >> Now rate4site comes out best, with only 1 protein where the mean evorate for CSA residues is ≥ control.

# Conclusion: rate4site_einsi_tree_1para is the best method to use, according to M-CSA, even at the individual protein level.

# Let's see what happens if I don't look at the protein level, but calculate overall median and mean evorates per evorate type instead:
qr %>% 
  group_by(evorate, type) %>% 
  summarise(median_rate = median(rate)) %>%
  pivot_wider(names_from = "type", values_from = "median_rate") %>%
  group_by(evorate, csa < control) %>%
  mutate(delta = csa - control) %>%
  arrange(desc(delta)) %>%
  print(n = 50)
# >> CSA is always < control for all evorate types.
# >> einsi_tree_1para has the strongest difference for all evorate types (values aren't comparable between types, but within each type einsi_tree_1para has the biggest delta).

# Looking at the biggest mean delta across proteins:
qr %>% 
  group_by(evorate, acc, type) %>% 
  summarise(mean_rate = mean(rate)) %>%
  pivot_wider(names_from = "type", values_from = "mean_rate") %>%
  filter(!is.na(control)) %>%
  mutate(delta = csa - control) %>%
  group_by(evorate) %>%
  summarise(mean_delta = mean(delta)) %>% 
  arrange(desc(mean_delta)) %>%
  print(n = 50)
# Exactly as above.
# >> Again: CSA is always < control for all evorate types.
# >> einsi_tree_1para has the strongest difference for all evorate types (values aren't comparable between types, but within each type einsi_tree_1para has the biggest delta).

# Looking at the median instead:
qr %>%
  group_by(evorate, acc, type) %>% 
  summarise(median_rate = median(rate)) %>%
  pivot_wider(names_from = "type", values_from = "median_rate") %>%
  filter(!is.na(control)) %>%
  mutate(delta = csa - control) %>%
  group_by(evorate) %>%
  summarise(median_delta = median(delta)) %>% 
  arrange(desc(median_delta)) %>%
  print(n = 50)
# Exactly as above.
# >> Again: CSA is always < control for all evorate types.
# >> einsi_tree_1para has the strongest difference for all evorate types (values aren't comparable between types, but within each type einsi_tree_1para has the biggest delta).

# Using log10(rate) instead:
# Looking at the biggest mean delta across proteins:
qr %>% 
  mutate(rate = log10(rate)) %>%
  group_by(evorate, acc, type) %>% 
  summarise(mean_rate = mean(rate)) %>%
  pivot_wider(names_from = "type", values_from = "mean_rate") %>%
  filter(!is.na(control)) %>%
  mutate(delta = csa - control) %>%
  group_by(evorate) %>%
  summarise(mean_delta = mean(delta)) %>% 
  arrange(desc(mean_delta)) %>%
  print(n = 50)
# >> Again: CSA is always < control for all evorate types.
# >> rate4site_einsi_tree_1para has the strongest median difference now!

# Using log10(rate) instead:
# Looking at the median instead:
qr %>%
  mutate(rate = log10(rate)) %>%
  group_by(evorate, acc, type) %>% 
  summarise(median_rate = median(rate)) %>%
  pivot_wider(names_from = "type", values_from = "median_rate") %>%
  filter(!is.na(control)) %>%
  mutate(delta = csa - control) %>%
  group_by(evorate) %>%
  summarise(median_delta = median(delta)) %>% 
  arrange(desc(median_delta)) %>%
  print(n = 50)
# Exactly as above.
# >> Again: CSA is always < control for all evorate types.
# >> rate4site_einsi_tree_1para has the strongest median difference now!

# Conclusion: rate4site_einsi_tree_1para leads to the strongest mean and median per-acc delta, once we use log10(rate). It also reports the fewest proteins with CSA rate ≥ Control rate (1 or 2 out of 102).
# >> rate4site_einsi_tree_1para is best, according to M-CSA.
