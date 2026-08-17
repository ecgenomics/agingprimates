#!/usr/bin/env Rscript

# Nodal localization of the top 1% of untransformed-LQ PSS pairs.
# This follows Omar's shifts.per.node figures: every pair is assigned once to
# its MRCA on the trait-pruned S4 tree; node size is count and colour is rate.

suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(scales)
  library(cowplot)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath("plot_lq_shifts_per_node.R", mustWork = TRUE)
}
project_dir <- dirname(script_path)

pair_file <- file.path(project_dir, "results", "LQ.score_results.tsv")
dataset_file <- file.path(project_dir, "input", "longevity_quotient.tsv")
tree_file <- Sys.getenv(
  "PSS_TREE",
  unset = paste0(
    "/Users/fabio/pCloud Drive/Bio/Projects/active.research/",
    "traits.evolution.omar/score.approach/inputs/tree/",
    "science.abn7829_data_s4.nex.tree"
  )
)
figure_dir <- file.path(project_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(pair_file), file.exists(dataset_file), file.exists(tree_file))
pairs <- read.delim(pair_file, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c(
  "Species1", "Species2", "TraitValue1", "TraitValue2", "FinalScore"
)
if (!all(required_columns %in% names(pairs))) {
  stop("PSS table lacks required columns.")
}
if (any(!is.finite(pairs$FinalScore))) {
  stop("FinalScore contains non-finite values.")
}

source_data <- read.delim(
  dataset_file, check.names = FALSE, stringsAsFactors = FALSE
)
family_map <- unique(data.frame(
  species = source_data$accepted_tree_tip,
  family = source_data$family,
  stringsAsFactors = FALSE
))
if (anyNA(family_map$family) || anyDuplicated(family_map$species)) {
  stop("Family annotations do not map uniquely to species.")
}

species <- sort(unique(c(pairs$Species1, pairs$Species2)))
expected_pairs <- choose(length(species), 2)
if (nrow(pairs) != expected_pairs) {
  stop("PSS table is not the complete trait-specific pair universe.")
}

tree <- read.tree(tree_file)
missing_tips <- setdiff(species, tree$tip.label)
if (length(missing_tips)) {
  stop("LQ species absent from S4 tree: ", paste(missing_tips, collapse = ", "))
}
tree <- drop.tip(tree, setdiff(tree$tip.label, species))
family_map <- family_map[match(tree$tip.label, family_map$species), , drop = FALSE]
if (anyNA(family_map$family)) {
  stop("Family annotations are incomplete after tree pruning.")
}

family_palette <- c(
  Aotidae = "#E69F00", Atelidae = "#56B4E9",
  Callitrichidae = "#009E73", Cebidae = "#D55E00",
  Cercopithecidae = "#0072B2", Cheirogaleidae = "#A88F00",
  Daubentoniidae = "#5F6368", Galagidae = "#E76F51",
  Hominidae = "#00A6D6", Hylobatidae = "#5A8F29",
  Indriidae = "#C56A00", Lemuridae = "#2A9D8F",
  Lepilemuridae = "#A63D40", Lorisidae = "#3A6EA5",
  Pitheciidae = "#7A8500", Tarsiidae = "#E4572E"
)
unknown_families <- setdiff(unique(family_map$family), names(family_palette))
if (length(unknown_families)) {
  stop("Families missing from palette: ", paste(unknown_families, collapse = ", "))
}

phenotype_map <- unique(rbind(
  data.frame(species = pairs$Species1, phenotype = pairs$TraitValue1),
  data.frame(species = pairs$Species2, phenotype = pairs$TraitValue2)
))
if (anyNA(phenotype_map$phenotype) ||
    any(!is.finite(phenotype_map$phenotype)) ||
    any(phenotype_map$phenotype < 0) ||
    anyDuplicated(phenotype_map$species) ||
    !setequal(phenotype_map$species, tree$tip.label)) {
  stop("Untransformed LQ values do not map uniquely to pruned-tree tips.")
}

pair_mrca <- vapply(
  seq_len(nrow(pairs)),
  function(i) getMRCA(tree, c(pairs$Species1[[i]], pairs$Species2[[i]])),
  integer(1)
)
internal_nodes <- seq.int(Ntip(tree) + 1L, Ntip(tree) + Nnode(tree))
if (any(!pair_mrca %in% internal_nodes)) {
  stop("At least one pair did not map to an internal MRCA node.")
}
opportunity <- tabulate(
  match(pair_mrca, internal_nodes), nbins = length(internal_nodes)
)
if (sum(opportunity) != nrow(pairs) || any(opportunity < 1L)) {
  stop("Node opportunities do not partition the complete pair universe.")
}

tail_proportion <- 0.01
n_tail <- max(1L, ceiling(nrow(pairs) * tail_proportion))
ranking <- order(pairs$FinalScore, decreasing = TRUE, method = "radix")
tail_indices <- ranking[seq_len(n_tail)]
tail_count <- tabulate(
  match(pair_mrca[tail_indices], internal_nodes),
  nbins = length(internal_nodes)
)
tail_rate <- tail_count / opportunity
if (sum(tail_count) != n_tail) {
  stop("Nodal counts do not sum to the selected tail size.")
}

node_metrics <- data.frame(
  node = internal_nodes,
  opportunity = opportunity,
  top_1pct_pair_count = tail_count,
  top_1pct_pair_rate = tail_rate,
  stringsAsFactors = FALSE
)
write.table(
  node_metrics,
  file.path(figure_dir, "LQ.top_1pct.node_metrics.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

base_tree <- ggtree(
  tree, layout = "circular", colour = "grey68", linewidth = 0.34
)
tree_coordinates <- as.data.frame(base_tree$data)
node_metrics <- merge(
  node_metrics,
  tree_coordinates[, c("node", "x", "y")],
  by = "node", all.x = TRUE, sort = FALSE
)
if (anyNA(node_metrics$x) || anyNA(node_metrics$y)) {
  stop("Could not recover plotting coordinates for all internal nodes.")
}

tip_radius <- max(tree_coordinates$x[tree_coordinates$isTip])
tip_phenotypes <- merge(
  tree_coordinates[tree_coordinates$isTip, c("node", "label", "x", "y")],
  phenotype_map, by.x = "label", by.y = "species", all.x = TRUE,
  sort = FALSE
)
tip_phenotypes <- merge(
  tip_phenotypes, family_map, by.x = "label", by.y = "species",
  all.x = TRUE, sort = FALSE
)
if (nrow(tip_phenotypes) != Ntip(tree) ||
    anyNA(tip_phenotypes$phenotype) || anyNA(tip_phenotypes$family)) {
  stop("Could not align LQ values and families with plotted tips.")
}

phenotype_max <- max(tip_phenotypes$phenotype)
bar_start <- tip_radius * 1.015
bar_span <- tip_radius * 0.24
tip_phenotypes$bar_start <- bar_start
tip_phenotypes$bar_end <- bar_start +
  bar_span * tip_phenotypes$phenotype / phenotype_max

represented_families <- names(family_palette)[
  names(family_palette) %in% unique(tip_phenotypes$family)
]
count_breaks <- unique(pretty(c(0, max(node_metrics$top_1pct_pair_count)), n = 3))
count_breaks <- count_breaks[
  count_breaks >= 0 & count_breaks <= max(node_metrics$top_1pct_pair_count)
]
rate_limit <- max(node_metrics$top_1pct_pair_rate)
rate_breaks <- unique(c(0, rate_limit / 2, rate_limit))

tree_figure <- base_tree +
  geom_segment(
    data = tip_phenotypes,
    aes(x = bar_start, y = y, xend = bar_end, yend = y, colour = family),
    inherit.aes = FALSE, linewidth = 1.35, lineend = "butt"
  ) +
  geom_point(
    data = tip_phenotypes,
    aes(x = x, y = y, colour = family),
    inherit.aes = FALSE, shape = 16, size = 1.35
  ) +
  geom_point(
    data = node_metrics,
    aes(x = x, y = y, size = top_1pct_pair_count, fill = top_1pct_pair_rate),
    inherit.aes = FALSE, shape = 21, colour = "grey18", stroke = 0.28
  ) +
  scale_fill_gradient(
    low = "#F3E8FF", high = "#7200FF", trans = "sqrt",
    limits = c(0, rate_limit), oob = squish,
    breaks = rate_breaks, labels = label_percent(accuracy = 0.1),
    name = "Nodal tail rate\n(C / opportunity)"
  ) +
  scale_size_continuous(
    range = c(1.0, 8.0), breaks = count_breaks,
    limits = c(0, max(node_metrics$top_1pct_pair_count)),
    name = "Top-tail pair count"
  ) +
  scale_colour_manual(values = family_palette, guide = "none") +
  xlim(0, tip_radius * 1.30) +
  labs(
    title = "LQ (untransformed): top 1% PSS pairs",
    subtitle = paste0(
      Ntip(tree), " species · ", nrow(pairs), " available pairs · ",
      n_tail, " pairs in top tail"
    ),
    caption = paste0(
      "Each pair is assigned once to its MRCA. Colour is the observed top-tail ",
      "fraction among pairs available at that node; size is the raw count.\n",
      "Tip bars show untransformed LQ (LQ_mammal) on a linear zero-to-maximum ",
      "scale (maximum = ", format(phenotype_max, digits = 4, trim = TRUE),
      "); bars and terminal points share family colour."
    )
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10.5, hjust = 0.5),
    plot.caption = element_text(size = 9, colour = "grey25", hjust = 0.5),
    legend.position = "bottom", legend.box = "horizontal",
    legend.title = element_text(size = 9.5, face = "bold"),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 18, 10, 18)
  ) +
  guides(
    fill = guide_colourbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = unit(6.2, "cm"), barheight = unit(0.35, "cm"), order = 1
    ),
    size = guide_legend(
      title.position = "top", title.hjust = 0.5,
      override.aes = list(fill = "#A64DFF"), order = 2
    )
  )

family_legend_data <- data.frame(
  family = factor(represented_families, levels = represented_families),
  x = 1, y = seq_along(represented_families)
)
family_legend <- get_legend(
  ggplot(family_legend_data, aes(x = x, y = y, colour = family)) +
    geom_point(size = 3) +
    scale_colour_manual(
      values = family_palette, breaks = represented_families,
      drop = FALSE, name = "Family"
    ) +
    guides(colour = guide_legend(
      ncol = 1, byrow = TRUE,
      override.aes = list(shape = 15, size = 4)
    )) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.key.height = unit(0.44, "cm"),
      legend.margin = margin(0, 4, 0, 4)
    )
)

figure <- plot_grid(
  tree_figure, family_legend,
  nrow = 1, rel_widths = c(0.82, 0.18), align = "h", axis = "tb"
)

ggsave(
  file.path(figure_dir, "LQ.top_1pct.shifts_per_node.png"),
  figure, width = 12, height = 10.2, dpi = 300, bg = "white"
)
ggsave(
  file.path(figure_dir, "LQ.top_1pct.shifts_per_node.pdf"),
  figure, width = 12, height = 10.2, bg = "white"
)
