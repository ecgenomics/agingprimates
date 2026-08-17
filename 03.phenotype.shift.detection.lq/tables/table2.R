#!/usr/bin/env Rscript

# Table 2: foreground and background species selected from the species shared
# between the top and bottom 1% LQ PSS tails.
# Run from any directory with: Rscript /path/to/tables/table2.R

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
  normalizePath("table2.R", mustWork = TRUE)
}

table_dir <- dirname(script_path)
project_dir <- dirname(table_dir)
score_file <- file.path(project_dir, "results", "LQ.score_results.tsv")
dataset_file <- file.path(project_dir, "input", "longevity_quotient.tsv")
dataset_output <- file.path(table_dir, "table2.tsv")
markdown_output <- file.path(table_dir, "table2.md")
tail_proportion <- 0.01
fg_threshold <- 1.3
bg_threshold <- 0.8

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
required_source_columns <- c(
  "accepted_tree_tip", "superfamily", "family", "LQ_mammal"
)
if (!all(required_score_columns %in% names(scores))) {
  stop("The PSS result table lacks required columns.")
}
if (!all(required_source_columns %in% names(source_data))) {
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

shared_data <- phenotype_map[
  match(shared_species, phenotype_map$Species), , drop = FALSE
]
taxonomy <- source_data[
  match(shared_data$Species, source_data$accepted_tree_tip),
  c("accepted_tree_tip", "superfamily", "family", "LQ_mammal"),
  drop = FALSE
]
if (anyNA(taxonomy$family) ||
    !isTRUE(all.equal(shared_data$LQ, taxonomy$LQ_mammal, tolerance = 1e-12))) {
  stop("Shared species could not be verified against the source dataset.")
}

shared_data$Group <- ifelse(
  shared_data$LQ > fg_threshold, "FG",
  ifelse(shared_data$LQ < bg_threshold, "BG", "Unselected")
)
selected <- shared_data[shared_data$Group %in% c("FG", "BG"), , drop = FALSE]
selected$Superfamily <- taxonomy$superfamily[
  match(selected$Species, taxonomy$accepted_tree_tip)
]
selected$Family <- taxonomy$family[
  match(selected$Species, taxonomy$accepted_tree_tip)
]

top_counts <- table(c(top_pairs$Species1, top_pairs$Species2))
bottom_counts <- table(c(bottom_pairs$Species1, bottom_pairs$Species2))
count_occurrences <- function(species, counts) {
  values <- unname(counts[species])
  values[is.na(values)] <- 0L
  as.integer(values)
}
selected$TopTailPairCount <- count_occurrences(selected$Species, top_counts)
selected$BottomTailPairCount <- count_occurrences(selected$Species, bottom_counts)
selected$SelectionCriterion <- ifelse(
  selected$Group == "FG", "LQ > 1.3", "LQ < 0.8"
)

selected$group_order <- match(selected$Group, c("FG", "BG"))
selected$value_order <- ifelse(selected$Group == "FG", -selected$LQ, selected$LQ)
selected <- selected[
  order(selected$group_order, selected$value_order, selected$Species),
  , drop = FALSE
]

dataset <- selected[, c(
  "Group", "Species", "Superfamily", "Family", "LQ",
  "SelectionCriterion", "TopTailPairCount", "BottomTailPairCount"
)]
names(dataset) <- c(
  "group", "species", "superfamily", "family", "LQ_mammal",
  "selection_criterion", "top_1pct_pair_count", "bottom_1pct_pair_count"
)
if (sum(dataset$group == "FG") != 7L || sum(dataset$group == "BG") != 6L) {
  stop("Unexpected FG or BG species count.")
}
if (any(dataset$top_1pct_pair_count < 1L) ||
    any(dataset$bottom_1pct_pair_count < 1L)) {
  stop("A selected species is not represented in both PSS tails.")
}

display <- data.frame(
  Group = dataset$group,
  Species = paste0("*", gsub("_", " ", dataset$species, fixed = TRUE), "*"),
  Superfamily = dataset$superfamily,
  Family = dataset$family,
  `Absolute LQ` = sprintf("%.3f", dataset$LQ_mammal),
  Criterion = dataset$selection_criterion,
  `Top-tail pairs` = dataset$top_1pct_pair_count,
  `Bottom-tail pairs` = dataset$bottom_1pct_pair_count,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

escape_markdown <- function(value) {
  gsub("|", "\\|", as.character(value), fixed = TRUE)
}
markdown_row <- function(values) {
  paste0("| ", paste(escape_markdown(values), collapse = " | "), " |")
}
markdown_lines <- c(
  "# Table 2. Foreground and background species selected by absolute LQ",
  "",
  paste0(
    "Species are selected from the 47-species top–bottom shared-tail set ",
    "using strict thresholds (7 FG and 6 BG species)."
  ),
  "",
  markdown_row(names(display)),
  markdown_row(rep("---", ncol(display))),
  apply(display, 1L, markdown_row),
  ""
)

temporary_dataset <- tempfile(fileext = ".tsv")
temporary_markdown <- tempfile(fileext = ".md")
write.table(
  dataset, temporary_dataset, sep = "\t", quote = FALSE,
  row.names = FALSE, na = "NA"
)
writeLines(markdown_lines, temporary_markdown, useBytes = TRUE)
if (file.info(temporary_dataset)$size == 0L ||
    file.info(temporary_markdown)$size == 0L) {
  stop("Temporary table outputs were not created correctly.")
}
for (output in c(dataset_output, markdown_output)) {
  if (file.exists(output) && !file.remove(output)) {
    stop("Previous output could not be replaced: ", output)
  }
}
if (!file.copy(temporary_dataset, dataset_output, overwrite = FALSE) ||
    !file.copy(temporary_markdown, markdown_output, overwrite = FALSE)) {
  stop("Table outputs could not be copied into the tables directory.")
}
unlink(c(temporary_dataset, temporary_markdown))

cat(sprintf(
  paste0(
    "Created table2.tsv and table2.md: %d FG species (LQ > %.1f), ",
    "%d BG species (LQ < %.1f).\n"
  ),
  sum(dataset$group == "FG"), fg_threshold,
  sum(dataset$group == "BG"), bg_threshold
))
