#!/usr/bin/env Rscript

# Image3: absolute LQ distribution and individual positions of species shared
# between the top and bottom 1% PSS tails.
# Run from any directory with: Rscript /path/to/images/Image3.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
})

arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", arguments, value = TRUE)
script_path <- if (length(file_argument)) {
  normalizePath(
    gsub(
      "~+~", " ", sub("^--file=", "", file_argument[[1]]),
      fixed = TRUE
    ),
    mustWork = TRUE
  )
} else {
  normalizePath("Image3.R", mustWork = TRUE)
}

image_dir <- dirname(script_path)
project_dir <- dirname(image_dir)
score_file <- file.path(project_dir, "results", "LQ.score_results.tsv")
dataset_file <- file.path(project_dir, "input", "longevity_quotient.tsv")
output_file <- file.path(image_dir, "Image3.png")
tail_proportion <- 0.01

stopifnot(file.exists(score_file), file.exists(dataset_file))
scores <- read.delim(
  score_file, check.names = FALSE, stringsAsFactors = FALSE
)
source_data <- read.delim(
  dataset_file, check.names = FALSE, stringsAsFactors = FALSE
)
required_score_columns <- c(
  "Species1", "Species2", "TraitValue1", "TraitValue2", "FinalScore"
)
if (!all(required_score_columns %in% names(scores))) {
  stop("The PSS result table lacks required columns.")
}
if (!all(c("accepted_tree_tip", "family", "LQ_mammal") %in% names(source_data))) {
  stop("The source dataset lacks required taxonomy or untransformed LQ columns.")
}
if (any(!is.finite(scores$FinalScore))) {
  stop("FinalScore contains non-finite values.")
}

tail_size <- max(1L, ceiling(nrow(scores) * tail_proportion))
ranking <- order(scores$FinalScore, decreasing = TRUE, method = "radix")
top_pairs <- scores[ranking[seq_len(tail_size)], , drop = FALSE]
bottom_pairs <- scores[tail(ranking, tail_size), , drop = FALSE]
top_species <- unique(c(top_pairs$Species1, top_pairs$Species2))
bottom_species <- unique(c(bottom_pairs$Species1, bottom_pairs$Species2))
shared_species <- sort(intersect(top_species, bottom_species))

phenotype_map <- unique(rbind(
  data.frame(Species = scores$Species1, LQ = scores$TraitValue1),
  data.frame(Species = scores$Species2, LQ = scores$TraitValue2)
))
if (anyDuplicated(phenotype_map$Species) || any(!is.finite(phenotype_map$LQ))) {
  stop("Absolute LQ values do not map uniquely to species.")
}
family_map <- unique(data.frame(
  Species = source_data$accepted_tree_tip,
  Family = source_data$family,
  stringsAsFactors = FALSE
))
if (anyDuplicated(family_map$Species) || anyNA(family_map$Family)) {
  stop("Family annotations do not map uniquely to species.")
}

plot_data <- phenotype_map[
  match(shared_species, phenotype_map$Species), , drop = FALSE
]
plot_data <- merge(
  plot_data, family_map, by = "Species", all.x = TRUE, sort = FALSE
)
if (nrow(plot_data) != length(shared_species) ||
    anyNA(plot_data$LQ) || anyNA(plot_data$Family)) {
  stop("Shared species could not be aligned with LQ and family data.")
}
if (!all.equal(
  unname(plot_data$LQ),
  unname(source_data$LQ_mammal[
    match(plot_data$Species, source_data$accepted_tree_tip)
  ]),
  tolerance = 1e-12
)) {
  stop("PSS phenotype values do not match the untransformed LQ source column.")
}

plot_data <- plot_data[order(plot_data$LQ, plot_data$Species), , drop = FALSE]
fg_threshold <- 1.3
bg_threshold <- 0.8
plot_data$Status <- factor(
  ifelse(
    plot_data$LQ > fg_threshold, "FG",
    ifelse(plot_data$LQ < bg_threshold, "BG", "Unselected")
  ),
  levels = c("BG", "Unselected", "FG")
)
plain_species_label <- gsub("_", " ", plot_data$Species, fixed = TRUE)
plot_data$SpeciesLabel <- ifelse(
  plot_data$Status == "FG", paste0("FG  |  ", plain_species_label),
  ifelse(
    plot_data$Status == "BG", paste0("BG  |  ", plain_species_label),
    plain_species_label
  )
)
plot_data$SpeciesLabel <- factor(
  plot_data$SpeciesLabel, levels = plot_data$SpeciesLabel
)

status_palette <- c(
  BG = "#0072B2",
  Unselected = "#9CA3AF",
  FG = "#D55E00"
)
stopifnot(sum(plot_data$Status == "FG") == 7L)
stopifnot(sum(plot_data$Status == "BG") == 6L)

phenotype_range <- range(plot_data$LQ)
x_padding <- diff(phenotype_range) * 0.035
x_limits <- phenotype_range + c(-x_padding, x_padding)
median_lq <- median(plot_data$LQ)
density_bandwidth <- bw.nrd0(plot_data$LQ)

density_panel <- ggplot(plot_data, aes(x = LQ)) +
  annotate(
    "rect", xmin = -Inf, xmax = bg_threshold,
    ymin = -Inf, ymax = Inf, fill = status_palette[["BG"]], alpha = 0.08
  ) +
  annotate(
    "rect", xmin = fg_threshold, xmax = Inf,
    ymin = -Inf, ymax = Inf, fill = status_palette[["FG"]], alpha = 0.08
  ) +
  geom_density(
    bw = density_bandwidth, kernel = "gaussian", trim = TRUE,
    fill = "#FACC15", colour = "#8A6500", linewidth = 0.9,
    alpha = 0.72
  ) +
  geom_rug(
    aes(colour = Status, linewidth = Status), sides = "b",
    outside = FALSE, show.legend = FALSE
  ) +
  geom_vline(
    xintercept = c(bg_threshold, fg_threshold),
    linetype = "dotdash", linewidth = 0.65,
    colour = c(status_palette[["BG"]], status_palette[["FG"]])
  ) +
  geom_vline(
    xintercept = median_lq, linetype = "dashed",
    linewidth = 0.7, colour = "grey28"
  ) +
  annotate(
    "text", x = median_lq, y = Inf,
    label = paste0("Median = ", sprintf("%.3f", median_lq)),
    hjust = -0.08, vjust = 1.5, size = 3.5, colour = "grey20"
  ) +
  annotate(
    "text", x = bg_threshold, y = Inf,
    label = "BG: LQ < 0.8 (n = 6)",
    hjust = 1.05, vjust = 3.2, size = 3.3,
    fontface = "bold", colour = status_palette[["BG"]]
  ) +
  annotate(
    "text", x = Inf, y = Inf,
    label = "FG: LQ > 1.3 (n = 7)",
    hjust = 1.05, vjust = 3.2, size = 3.3,
    fontface = "bold", colour = status_palette[["FG"]]
  ) +
  scale_colour_manual(values = status_palette) +
  scale_linewidth_manual(
    values = c(BG = 1.4, Unselected = 0.45, FG = 1.4),
    guide = "none"
  ) +
  scale_x_continuous(expand = expansion(mult = 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  coord_cartesian(xlim = x_limits) +
  labs(
    title = "Density of absolute LQ among species shared by both PSS tails",
    subtitle = paste0(
      length(shared_species),
      " shared-tail species · 7 FG (LQ > 1.3) · 6 BG (LQ < 0.8)"
    ),
    x = NULL,
    y = "Probability density"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 10.5, colour = "grey25"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(8, 18, 0, 8)
  )

species_panel <- ggplot(
  plot_data, aes(x = LQ, y = SpeciesLabel, colour = Status)
) +
  geom_hline(
    data = plot_data[plot_data$Status != "Unselected", , drop = FALSE],
    aes(yintercept = as.numeric(SpeciesLabel), colour = Status),
    inherit.aes = FALSE, linewidth = 3.2, alpha = 0.10
  ) +
  geom_segment(
    aes(x = x_limits[[1]], xend = LQ, yend = SpeciesLabel),
    colour = "grey88", linewidth = 0.35, show.legend = FALSE
  ) +
  geom_point(aes(size = Status, shape = Status), alpha = 0.98) +
  geom_vline(
    xintercept = c(bg_threshold, fg_threshold),
    linetype = "dotdash", linewidth = 0.65,
    colour = c(status_palette[["BG"]], status_palette[["FG"]])
  ) +
  geom_vline(
    xintercept = median_lq, linetype = "dashed",
    linewidth = 0.7, colour = "grey28"
  ) +
  scale_colour_manual(
    values = status_palette,
    breaks = c("FG", "BG", "Unselected"),
    labels = c("Foreground: LQ > 1.3", "Background: LQ < 0.8", "Unselected"),
    name = "Selection"
  ) +
  scale_shape_manual(
    values = c(BG = 24, Unselected = 16, FG = 25),
    breaks = c("FG", "BG", "Unselected"),
    labels = c("Foreground: LQ > 1.3", "Background: LQ < 0.8", "Unselected"),
    name = "Selection"
  ) +
  scale_size_manual(
    values = c(BG = 3.2, Unselected = 1.8, FG = 3.2),
    breaks = c("FG", "BG", "Unselected"),
    labels = c("Foreground: LQ > 1.3", "Background: LQ < 0.8", "Unselected"),
    name = "Selection"
  ) +
  scale_x_continuous(
    limits = x_limits, expand = expansion(mult = 0),
    breaks = pretty(phenotype_range, n = 6)
  ) +
  labs(
    x = "Untransformed longevity quotient (LQ_mammal)",
    y = NULL,
    caption = paste0(
      "FG and BG rows carry explicit name prefixes. Points show individual ",
      "values; dot-dash lines mark selection thresholds and the dashed line ",
      "marks the median."
    )
  ) +
  theme_classic(base_size = 9.5) +
  theme(
    axis.title.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 7.6, face = "italic"),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.35),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 8.5),
    legend.key.height = unit(0.38, "cm"),
    plot.caption = element_text(size = 8.8, colour = "grey25", hjust = 0),
    plot.margin = margin(0, 18, 8, 8)
  )

figure <- plot_grid(
  density_panel, species_panel,
  ncol = 1, rel_heights = c(0.34, 0.66), align = "v", axis = "lr"
)

temporary_output <- tempfile(fileext = ".png")
ggsave(
  temporary_output, figure,
  width = 12, height = 13.5, units = "in", dpi = 300, bg = "white"
)
if (!file.exists(temporary_output) || file.info(temporary_output)$size == 0L) {
  stop("The temporary PNG was not created correctly.")
}
if (file.exists(output_file) && !file.remove(output_file)) {
  stop("The previous PNG could not be replaced.")
}
if (!file.copy(temporary_output, output_file, overwrite = FALSE)) {
  stop("The PNG could not be copied into the images directory.")
}
unlink(temporary_output)

cat(sprintf(
  paste0(
    "Created Image3.png for %d shared species; absolute LQ range %.6f–%.6f; ",
    "median %.6f; %d FG and %d BG species.\n"
  ),
  nrow(plot_data), min(plot_data$LQ), max(plot_data$LQ), median_lq,
  sum(plot_data$Status == "FG"), sum(plot_data$Status == "BG")
))
