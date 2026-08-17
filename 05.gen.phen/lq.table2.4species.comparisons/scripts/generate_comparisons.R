#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
n_comparisons <- if (length(args) >= 1L) as.integer(args[[1L]]) else 100L
seed <- if (length(args) >= 2L) as.integer(args[[2L]]) else 260811L

if (is.na(n_comparisons) || n_comparisons < 1L) {
  stop("n_comparisons must be a positive integer")
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
root <- dirname(script_dir)

source_table <- if (length(args) >= 3L) {
  normalizePath(args[[3L]], mustWork = TRUE)
} else {
  normalizePath(
    file.path(root, "..", "..", "03.phenotype.shift.detection.lq", "tables", "table2.tsv"),
    mustWork = TRUE
  )
}

input_dir <- file.path(root, "input")
caas_dir <- file.path(root, "configurations", "caas")
rerc_dir <- file.path(root, "configurations", "rerconverge")
table_dir <- file.path(root, "tables")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(caas_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rerc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

table2 <- read.delim(source_table, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c("group", "species")
if (!all(required_columns %in% names(table2))) {
  stop("Input table must contain columns: group and species")
}
if (anyDuplicated(table2$species)) stop("Species names must be unique")
if (!all(table2$group %in% c("FG", "BG"))) stop("group values must be FG or BG")

fg_pool <- sort(table2$species[table2$group == "FG"])
bg_pool <- sort(table2$species[table2$group == "BG"])
if (length(fg_pool) < 4L || length(bg_pool) < 4L) {
  stop("At least four FG and four BG species are required")
}

file.copy(source_table, file.path(input_dir, "table2.tsv"), overwrite = TRUE)

fg_quartets <- combn(fg_pool, 4L, simplify = FALSE)
bg_quartets <- combn(bg_pool, 4L, simplify = FALSE)
caas_candidates <- expand.grid(
  fg_index = seq_along(fg_quartets),
  bg_index = seq_along(bg_quartets),
  KEEP.OUT.ATTRS = FALSE
)
if (n_comparisons > nrow(caas_candidates)) {
  stop(sprintf(
    "Requested %d unique CAAStools comparisons, but only %d are possible",
    n_comparisons, nrow(caas_candidates)
  ))
}

set.seed(seed)
caas_selection <- caas_candidates[
  sample(seq_len(nrow(caas_candidates)), n_comparisons, replace = FALSE),
  ,
  drop = FALSE
]

# RERconverge uses every possible foreground quartet exactly once.
rerc_selection <- sample(seq_along(fg_quartets))

# Remove stale generated configurations before regeneration.
unlink(list.files(caas_dir, "[.]caas[.]cfg$", full.names = TRUE))
unlink(list.files(rerc_dir, "[.]rerc[.]cfg$", full.names = TRUE))

write_cfg <- function(path, fg, bg) {
  cfg <- data.frame(
    species = c(sort(fg), sort(bg)),
    state = c(rep(1L, length(fg)), rep(0L, length(bg)))
  )
  write.table(
    cfg, path, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE
  )
}

membership_rows <- vector("list", n_comparisons + length(rerc_selection))
comparison_rows <- vector("list", n_comparisons + length(rerc_selection))

for (i in seq_len(n_comparisons)) {
  comparison_id <- sprintf("%03d", i)

  caas_fg_index <- caas_selection$fg_index[[i]]
  caas_bg_index <- caas_selection$bg_index[[i]]
  caas_fg <- fg_quartets[[caas_fg_index]]
  caas_bg <- bg_quartets[[caas_bg_index]]
  write_cfg(
    file.path(caas_dir, paste0(comparison_id, ".caas.cfg")),
    caas_fg, caas_bg
  )
  caas_row <- i
  membership_rows[[caas_row]] <- data.frame(
    comparison_id = comparison_id,
    method = "caas",
    species = c(caas_fg, caas_bg),
    role = c(rep("FG", 4L), rep("BG", 4L)),
    stringsAsFactors = FALSE
  )
  comparison_rows[[caas_row]] <- data.frame(
    comparison_id = comparison_id,
    method = "caas",
    fg_combination_id = sprintf("FG%03d", caas_fg_index),
    fg_combination_occurrence = 1L,
    bg_strategy = "sampled_4_from_BG_pool",
    fg_species = paste(sort(caas_fg), collapse = ","),
    bg_species = paste(sort(caas_bg), collapse = ","),
    stringsAsFactors = FALSE
  )

}

for (i in seq_along(rerc_selection)) {
  comparison_id <- sprintf("%03d", i)
  rerc_fg_index <- rerc_selection[[i]]
  rerc_fg <- fg_quartets[[rerc_fg_index]]
  rerc_bg <- bg_pool
  write_cfg(
    file.path(rerc_dir, paste0(comparison_id, ".rerc.cfg")),
    rerc_fg, rerc_bg
  )
  rerc_row <- n_comparisons + i
  membership_rows[[rerc_row]] <- data.frame(
    comparison_id = comparison_id,
    method = "rerconverge",
    species = c(rerc_fg, rerc_bg),
    role = c(rep("FG", 4L), rep("BG", length(rerc_bg))),
    stringsAsFactors = FALSE
  )
  comparison_rows[[rerc_row]] <- data.frame(
    comparison_id = comparison_id,
    method = "rerconverge",
    fg_combination_id = sprintf("FG%03d", rerc_fg_index),
    fg_combination_occurrence = 1L,
    bg_strategy = "fixed_complete_BG_pool",
    fg_species = paste(sort(rerc_fg), collapse = ","),
    bg_species = paste(sort(rerc_bg), collapse = ","),
    stringsAsFactors = FALSE
  )
}

membership <- do.call(rbind, membership_rows)
comparison_manifest <- do.call(rbind, comparison_rows)
write.table(
  membership, file.path(table_dir, "comparison_species_membership.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  comparison_manifest, file.path(table_dir, "comparison_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

metadata <- c(
  paste0("source_table\t", source_table),
  paste0("seed\t", seed),
  paste0("n_caas_comparisons\t", n_comparisons),
  paste0("n_rerconverge_comparisons\t", length(rerc_selection)),
  paste0("fg_pool_n\t", length(fg_pool)),
  paste0("bg_pool_n\t", length(bg_pool)),
  paste0("possible_fg_quartets\t", length(fg_quartets)),
  paste0("possible_bg_quartets\t", length(bg_quartets)),
  paste0("possible_unique_caas_comparisons\t", nrow(caas_candidates)),
  paste0("unique_caas_comparisons_generated\t", n_comparisons),
  paste0("unique_rerconverge_fg_combinations_generated\t", length(unique(rerc_selection))),
  "caas_design\t4 sampled FG vs 4 sampled BG; all generated configurations unique",
  "rerconverge_design\tall 35 unique 4-species FG combinations vs fixed complete 6-species BG",
  "cfg_format\ttab-separated species and binary state, no header; FG=1, BG=0"
)
writeLines(metadata, file.path(table_dir, "generation_metadata.tsv"))

cat(sprintf(
  "Generated %d CAAStools files and %d RERconverge files with seed %d.\n",
  n_comparisons, length(rerc_selection), seed
))
cat(sprintf(
  "CAAStools: %d unique configurations. RERconverge: %d unique FG quartets, no duplicates.\n",
  n_comparisons, length(unique(rerc_selection))
))
