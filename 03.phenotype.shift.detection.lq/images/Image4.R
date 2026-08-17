#!/usr/bin/env Rscript

# Image4: phylogenetic placement and absolute LQ of species shared between the
# top and bottom 1% PSS tails.
# Run from any directory with: Rscript /path/to/images/Image4.R

suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
  library(ggtree)
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
  normalizePath("Image4.R", mustWork = TRUE)
}

image_dir <- dirname(script_path)
project_dir <- dirname(image_dir)
score_file <- file.path(project_dir, "results", "LQ.score_results.tsv")
dataset_file <- file.path(project_dir, "input", "longevity_quotient.tsv")
tree_file <- Sys.getenv(
  "PSS_TREE",
  unset = paste0(
    "/Users/fabio/pCloud Drive/Bio/Projects/active.research/",
    "traits.evolution.omar/score.approach/inputs/tree/",
    "science.abn7829_data_s4.nex.tree"
  )
)
output_file <- file.path(image_dir, "Image4.png")
tail_proportion <- 0.01

stopifnot(
  file.exists(score_file), file.exists(dataset_file), file.exists(tree_file)
)
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
family_map <- unique(data.frame(
  Species = source_data$accepted_tree_tip,
  Family = source_data$family,
  stringsAsFactors = FALSE
))
if (anyDuplicated(phenotype_map$Species) || any(!is.finite(phenotype_map$LQ))) {
  stop("Absolute LQ values do not map uniquely to species.")
}
if (anyDuplicated(family_map$Species) || anyNA(family_map$Family)) {
  stop("Family annotations do not map uniquely to species.")
}

tree <- read.tree(tree_file)
missing_tips <- setdiff(shared_species, tree$tip.label)
if (length(missing_tips)) {
  stop("Shared species absent from the S4 tree: ", paste(missing_tips, collapse = ", "))
}
tree <- drop.tip(tree, setdiff(tree$tip.label, shared_species))
tree <- ladderize(tree, right = FALSE)
if (Ntip(tree) != length(shared_species) ||
    !setequal(tree$tip.label, shared_species)) {
  stop("Unexpected tip set after pruning the S4 tree.")
}

plot_data <- phenotype_map[
  match(tree$tip.label, phenotype_map$Species), , drop = FALSE
]
plot_data <- merge(
  plot_data, family_map, by = "Species", all.x = TRUE, sort = FALSE
)
if (nrow(plot_data) != Ntip(tree) || anyNA(plot_data$LQ) ||
    anyNA(plot_data$Family)) {
  stop("Tree tips could not be aligned with LQ and family data.")
}
source_lq <- source_data$LQ_mammal[
  match(plot_data$Species, source_data$accepted_tree_tip)
]
if (!isTRUE(all.equal(plot_data$LQ, source_lq, tolerance = 1e-12))) {
  stop("PSS phenotype values do not match the untransformed LQ source column.")
}
fg_threshold <- 1.3
bg_threshold <- 0.8
plot_data$Status <- factor(
  ifelse(
    plot_data$LQ > fg_threshold, "FG",
    ifelse(plot_data$LQ < bg_threshold, "BG", "Unselected")
  ),
  levels = c("BG", "Unselected", "FG")
)
status_palette <- c(
  BG = "#0072B2",
  Unselected = "#A7A7A7",
  FG = "#D55E00"
)
stopifnot(sum(plot_data$Status == "FG") == 7L)
stopifnot(sum(plot_data$Status == "BG") == 6L)

base_tree <- ggtree(
  tree, layout = "rectangular", colour = "grey58", linewidth = 0.55
)
tree_coordinates <- as.data.frame(base_tree$data)
tip_coordinates <- tree_coordinates[
  tree_coordinates$isTip, c("label", "x", "y"), drop = FALSE
]
tip_data <- merge(
  tip_coordinates, plot_data,
  by.x = "label", by.y = "Species", all.x = TRUE, sort = FALSE
)
if (nrow(tip_data) != Ntip(tree) || anyNA(tip_data$LQ)) {
  stop("Plotting coordinates could not be aligned with tip metadata.")
}

tree_min <- min(tree_coordinates$x, na.rm = TRUE)
tree_max <- max(tree_coordinates$x, na.rm = TRUE)
tree_span <- tree_max - tree_min
if (!is.finite(tree_span) || tree_span <= 0) {
  stop("The pruned tree has an invalid horizontal span.")
}
bar_start <- tree_max + tree_span * 0.025
bar_span <- tree_span * 0.22
tip_data$bar_start <- bar_start
tip_data$bar_end <- bar_start + bar_span * tip_data$LQ / max(tip_data$LQ)
label_x <- bar_start + bar_span + tree_span * 0.025
tip_data$display_label <- paste0(
  ifelse(
    tip_data$Status == "FG", "FG  |  ",
    ifelse(tip_data$Status == "BG", "BG  |  ", "")
  ),
  gsub("_", " ", tip_data$label, fixed = TRUE),
  "   ", sprintf("%.3f", tip_data$LQ)
)

figure <- base_tree +
  geom_point(
    data = tip_data,
    aes(x = x, y = y, colour = Status, shape = Status, size = Status),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = tip_data,
    aes(
      x = bar_start, xend = bar_end, y = y, yend = y,
      colour = Status, linewidth = Status
    ),
    inherit.aes = FALSE, lineend = "butt"
  ) +
  geom_text(
    data = tip_data[tip_data$Status == "Unselected", , drop = FALSE],
    aes(x = label_x, y = y, label = display_label),
    inherit.aes = FALSE, hjust = 0, size = 2.65,
    fontface = "italic", colour = "grey25"
  ) +
  geom_text(
    data = tip_data[tip_data$Status == "FG", , drop = FALSE],
    aes(x = label_x, y = y, label = display_label),
    inherit.aes = FALSE, hjust = 0, size = 2.9,
    fontface = "bold.italic", colour = status_palette[["FG"]]
  ) +
  geom_text(
    data = tip_data[tip_data$Status == "BG", , drop = FALSE],
    aes(x = label_x, y = y, label = display_label),
    inherit.aes = FALSE, hjust = 0, size = 2.9,
    fontface = "bold.italic", colour = status_palette[["BG"]]
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
    values = c(BG = 3.0, Unselected = 1.6, FG = 3.0),
    breaks = c("FG", "BG", "Unselected"),
    labels = c("Foreground: LQ > 1.3", "Background: LQ < 0.8", "Unselected"),
    name = "Selection"
  ) +
  scale_linewidth_manual(
    values = c(BG = 3.6, Unselected = 1.6, FG = 3.6),
    guide = "none"
  ) +
  xlim(tree_min, label_x + tree_span * 0.62) +
  labs(
    title = "Phylogenetic placement of foreground and background LQ species",
    subtitle = paste0(
      Ntip(tree),
      " shared-tail species · 7 FG (LQ > 1.3) · 6 BG (LQ < 0.8)"
    ),
    caption = paste0(
      "The Kuderna et al. S4 tree is pruned to species represented in both ",
      "PSS tails. Bar length and labels report absolute LQ; selected names ",
      "are bold and carry explicit FG or BG prefixes."
    )
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 10.5, colour = "grey25"),
    plot.caption = element_text(size = 8.8, colour = "grey25", hjust = 0),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(12, 25, 12, 15)
  ) +
  guides(
    colour = guide_legend(order = 1),
    shape = guide_legend(order = 1),
    size = guide_legend(order = 1)
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
    "Created Image4.png: S4 phylogeny pruned to %d shared-tail species; ",
    "%d FG and %d BG species.\n"
  ),
  Ntip(tree), sum(tip_data$Status == "FG"), sum(tip_data$Status == "BG")
))
