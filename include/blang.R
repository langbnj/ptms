# Shared settings and plotting helpers for the R scripts in this repository.
# Each script starts with blang_init(), which sources this file. To define it,
# add a line like this to ~/.Rprofile, with the path to your copy of the repository:
#   blang_init <- function() source("~/ptms/include/blang.R")

# Start from a clean workspace, keeping blang_init() and any MySQL connection opened by mysql.R
rm(list = setdiff(ls(), c("blang_init", "superreservedcon", "superreserveddrv")))

options(nwarnings = 10000)
options(tibble.print_min = 6)

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(magrittr)
  library(ggrepel)
  library(measurements)
  library(ggbeeswarm)
})

# String interpolation, like Python's f-strings
f <- stringr::str_glue

# Colours
ptmod <- "#F77200"  # dark orange
ptmo <- "#FF9300"   # orange
ptmbd <- "#003D51"  # dark blue
ptmbl <- "#00799E"  # light blue
ptmgr <- "#46ba1f"  # green
ptmco <- "#e7a13d"  # light orange
ptmvir1 <- scales::viridis_pal()(5)[4]
ptmvir0 <- scales::viridis_pal()(5)[2]
ptmb <- ptmvir1
ptmg <- ptmvir0
ptmcol <- c("Modified" = ptmvir1, "Control" = ptmvir0)

# In RStudio, work in the folder of the script being run
try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)), silent = TRUE)

# Save a plot (default: the last one) with the given size in mm. Uses the macOS
# quartz device where available, otherwise ggsave().
qsave <- function(plot_or_filename = NULL, filename = "qsave.pdf", type = "pdf", dpi = 600,
                  width = 40, height = 40, units = "mm", family = "Helvetica Neue") {
  if (inherits(plot_or_filename, "ggplot")) {
    myplot <- plot_or_filename
  } else {
    myplot <- last_plot()
    filename <- plot_or_filename
  }
  if (capabilities("aqua")) {
    quartz(title = "ggplot", file = filename, type = type, dpi = dpi,
           width = conv_unit(width, units, "inch"), height = conv_unit(height, units, "inch"),
           family = family)
    print(myplot)
    invisible(dev.off())
  } else {
    ggsave(filename, myplot, width = width, height = height, units = units, dpi = dpi,
           device = if (type == "pdf") cairo_pdf else type)
  }
}

# ggplot2 linewidth per point of stroke width (cf. .pt), so that 0.5 / .weight gives 0.5 pt
.weight <- 2.134

# Default point shape and line widths
update_geom_defaults("point",   list(shape = 16, size = 0.5))
update_geom_defaults("bar",     list(linewidth = 0.5 / .weight))
update_geom_defaults("col",     list(linewidth = 0.5 / .weight))
update_geom_defaults("density", list(linewidth = 0.5 / .weight))
update_geom_defaults("hline",   list(linewidth = 0.5 / .weight))
update_geom_defaults("line",    list(linewidth = 0.5 / .weight))
update_geom_defaults("ribbon",  list(linewidth = 0.5 / .weight))
update_geom_defaults("vline",   list(linewidth = 0.5 / .weight))

# Default text
update_geom_defaults("text",        list(family = "Helvetica Neue"))
update_geom_defaults("text_repel",  list(family = "Helvetica Neue", size = 5 / .pt, segment.size = 0.5 / .weight))
update_geom_defaults("label_repel", list(family = "Helvetica Neue", size = 5 / .pt, segment.size = 0.5 / .weight))

# Nature-style theme.
# legend_position: "bottom" or "bottomright" (below the plot, level with the x-axis title),
#   "top" or "topright", "topleft", "topcentre" (inside the plot)
# legend_nudge_top, legend_nudge_right: move the legend down or right by this many mm
# extra_margin_bottom, extra_margin_right: extra space in mm, e.g. for long axis labels
theme_nature <- function(axis_fontsize = 5, legend_position = "bottom", legend_nudge_top = 0,
                         legend_nudge_right = 0, extra_margin_bottom = 0, extra_margin_right = 0) {
  if (legend_position %in% c("bottom", "bottomright")) {
    legend_position <- c(1, 0)
    legend_justification <- c(1, 1)
    legend_margin <- margin(((6.5 + axis_fontsize) / .pt) + legend_nudge_top, -legend_nudge_right,
                            5.5 / .pt, 5.5 / .pt, "mm")
  } else if (legend_position %in% c("top", "topright")) {
    legend_position <- c(1, 1)
    legend_justification <- c(1, 1)
    legend_margin <- margin(legend_nudge_top, -legend_nudge_right, 5.5 / .pt, 5.5 / .pt, "mm")
  } else if (legend_position == "topleft") {
    legend_position <- c(0, 1)
    legend_justification <- c(0, 1)
    legend_margin <- margin(legend_nudge_top, 5.5 / .pt, 5.5 / .pt, legend_nudge_right, "mm")
  } else if (legend_position %in% c("topcenter", "topcentre")) {
    legend_position <- c(0.5, 1)
    legend_justification <- c(0.5, 1)
    legend_margin <- margin(legend_nudge_top, 5.5 / .pt, 5.5 / .pt, 5.5 / .pt + legend_nudge_right, "mm")
  }

  theme_classic() +
    theme(
      plot.margin = margin(1, 1 + extra_margin_right, 1 + extra_margin_bottom, 1, unit = "mm"),
      text = element_text(family = "Helvetica Neue", size = 6),
      plot.tag = element_text(size = 7, hjust = 0, vjust = 1),
      plot.tag.position = c(0, 1),
      plot.title = element_text(size = 6),
      axis.title = element_text(size = 6),
      axis.text = element_text(size = axis_fontsize, colour = "black"),
      axis.line = element_line(linewidth = 0.5 / .weight, lineend = "square"),
      axis.ticks = element_line(linewidth = 0.5 / .weight, colour = "black"),
      legend.key.size = unit(2, "mm"),
      legend.key.spacing.x = unit(1, "mm"),
      legend.text = element_text(size = 6, margin = margin(l = unit(2, "mm"), r = unit(0, "mm"))),
      legend.title = element_text(size = 6),
      legend.background = element_blank(),
      legend.position = "inside",
      legend.position.inside = legend_position,
      legend.justification = legend_justification,
      legend.direction = "horizontal",
      legend.box.spacing = unit(0, "mm"),
      legend.margin = legend_margin,
      strip.background = element_blank(),
      strip.text = element_text(size = 6, margin = margin(1, 1, 1, 4.4)),
      strip.text.y = element_text(angle = 0, hjust = 0)
    )
}
