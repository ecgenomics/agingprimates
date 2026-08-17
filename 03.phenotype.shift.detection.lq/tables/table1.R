#!/usr/bin/env Rscript

# Table 1: top-1% LQ PSS pairs whose two endpoint species are both represented
# in the bottom 1% of LQ PSS pairs.
# Run from any directory with: Rscript /path/to/tables/table1.R

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
  normalizePath("table1.R", mustWork = TRUE)
}

table_dir <- dirname(script_path)
project_dir <- dirname(table_dir)
score_file <- file.path(project_dir, "results", "LQ.score_results.tsv")
output_file <- file.path(table_dir, "table1.md")
tail_proportion <- 0.01

stopifnot(file.exists(score_file))
scores <- read.delim(
  score_file, check.names = FALSE, stringsAsFactors = FALSE
)
required_columns <- c(
  "Species1", "Species2", "TraitValue1", "TraitValue2",
  "PatristicDistance", "FinalScore"
)
if (!all(required_columns %in% names(scores))) {
  stop("The PSS result table lacks required columns.")
}
if (any(!is.finite(scores$FinalScore))) {
  stop("FinalScore contains non-finite values.")
}

tail_size <- max(1L, ceiling(nrow(scores) * tail_proportion))
ranking <- order(scores$FinalScore, decreasing = TRUE, method = "radix")
top_indices <- ranking[seq_len(tail_size)]
bottom_indices <- tail(ranking, tail_size)
top_pairs <- scores[top_indices, , drop = FALSE]
bottom_pairs <- scores[bottom_indices, , drop = FALSE]

bottom_species <- sort(unique(c(
  bottom_pairs$Species1, bottom_pairs$Species2
)))
top_species <- sort(unique(c(
  top_pairs$Species1, top_pairs$Species2
)))
shared_species <- sort(intersect(top_species, bottom_species))
bottom_pair_counts <- table(c(
  bottom_pairs$Species1, bottom_pairs$Species2
))

# Follow Omar's strict shared-tail network definition: first identify species
# represented in both tails, then retain a top-tail pair only when both of its
# endpoints belong to that shared-species set.
keep <- top_pairs$Species1 %in% shared_species &
  top_pairs$Species2 %in% shared_species
selected <- top_pairs[keep, , drop = FALSE]
selected$TopRank <- match(rownames(selected), rownames(top_pairs))

count_bottom_occurrences <- function(species) {
  counts <- unname(bottom_pair_counts[species])
  counts[is.na(counts)] <- 0L
  as.integer(counts)
}

table_data <- data.frame(
  `Top rank` = selected$TopRank,
  `Species 1` = selected$Species1,
  `LQ 1` = selected$TraitValue1,
  `Species 2` = selected$Species2,
  `LQ 2` = selected$TraitValue2,
  `Absolute LQ difference` = abs(
    selected$TraitValue1 - selected$TraitValue2
  ),
  `Patristic distance` = selected$PatristicDistance,
  `PSS FinalScore` = selected$FinalScore,
  `Bottom-tail pairs involving species 1` =
    count_bottom_occurrences(selected$Species1),
  `Bottom-tail pairs involving species 2` =
    count_bottom_occurrences(selected$Species2),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!nrow(table_data)) {
  stop("No top-tail pair satisfies the requested shared-species criterion.")
}
if (any(
  table_data[["Bottom-tail pairs involving species 1"]] == 0L |
  table_data[["Bottom-tail pairs involving species 2"]] == 0L
)) {
  stop("Selection invariant failed: a retained endpoint lacks a bottom-tail pair.")
}
if (anyDuplicated(paste(table_data[["Species 1"]], table_data[["Species 2"]]))) {
  stop("Duplicated species pairs in the selected table.")
}

format_species <- function(species) {
  paste0("*", gsub("_", " ", species, fixed = TRUE), "*")
}
format_number <- function(value, digits) {
  formatC(value, format = "f", digits = digits)
}

display_data <- table_data
display_data[["Species 1"]] <- format_species(display_data[["Species 1"]])
display_data[["Species 2"]] <- format_species(display_data[["Species 2"]])
display_data[["LQ 1"]] <- format_number(display_data[["LQ 1"]], 3)
display_data[["LQ 2"]] <- format_number(display_data[["LQ 2"]], 3)
display_data[["Absolute LQ difference"]] <- format_number(
  display_data[["Absolute LQ difference"]], 3
)
display_data[["Patristic distance"]] <- format_number(
  display_data[["Patristic distance"]], 3
)
display_data[["PSS FinalScore"]] <- format_number(
  display_data[["PSS FinalScore"]], 6
)

escape_markdown <- function(value) {
  gsub("|", "\\|", as.character(value), fixed = TRUE)
}
markdown_row <- function(values) {
  paste0("| ", paste(escape_markdown(values), collapse = " | "), " |")
}

lines <- c(
  "# Table 1. Top-1% LQ PSS pairs linked to the bottom 1% by species",
  "",
  paste0(
    "Top-tail pairs retained when both endpoint species belong to the ",
    "top–bottom shared-species set (", nrow(table_data), " retained pairs)."
  ),
  "",
  markdown_row(names(display_data)),
  markdown_row(rep("---", ncol(display_data))),
  apply(display_data, 1L, markdown_row),
  ""
)
temporary_output <- tempfile(fileext = ".md")
writeLines(lines, temporary_output, useBytes = TRUE)
if (!file.exists(temporary_output) || file.info(temporary_output)$size == 0L) {
  stop("The temporary Markdown table was not created correctly.")
}
if (file.exists(output_file) && !file.remove(output_file)) {
  stop("The previous Markdown table could not be replaced.")
}
if (!file.copy(temporary_output, output_file, overwrite = FALSE)) {
  stop("The Markdown table could not be copied into the tables directory.")
}
unlink(temporary_output)

cat(sprintf(
  paste0(
    "Complete pairs: %d; tail size: %d; shared-tail species: %d; ",
    "retained top pairs: %d.\n"
  ),
  nrow(scores), tail_size, length(shared_species), nrow(table_data)
))
