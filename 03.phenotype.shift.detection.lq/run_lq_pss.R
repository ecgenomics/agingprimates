#!/usr/bin/env Rscript

# Run the analytical PSS implementation from the Omar project on the
# untransformed mammal-wide longevity quotient.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath("run_lq_pss.R", mustWork = TRUE)
}
project_dir <- dirname(script_path)

dataset_file <- file.path(project_dir, "input", "longevity_quotient.tsv")
tree_file <- Sys.getenv(
  "PSS_TREE",
  unset = paste0(
    "/Users/fabio/pCloud Drive/Bio/Projects/active.research/",
    "traits.evolution.omar/score.approach/inputs/tree/",
    "science.abn7829_data_s4.nex.tree"
  )
)
pss_core <- Sys.getenv(
  "PSS_CORE",
  unset = paste0(
    "/Users/fabio/pCloud Drive/Bio/Projects/active.research/",
    "traits.evolution.omar/score.approach/01.method.development/",
    "R.package/pss.core.R"
  )
)
results_dir <- file.path(project_dir, "results")

stopifnot(file.exists(dataset_file), file.exists(tree_file), file.exists(pss_core))
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

source(pss_core)
source_data <- utils::read.delim(
  dataset_file, check.names = FALSE, stringsAsFactors = FALSE
)
required <- c("accepted_tree_tip", "LQ_mammal", "log_LQ_mammal")
if (!all(required %in% names(source_data))) {
  stop("Input dataset lacks required LQ columns.")
}

# Deliberately select LQ_mammal, never log_LQ_mammal.
analysis_data <- data.frame(
  SpeciesBROAD = source_data$accepted_tree_tip,
  LQ = suppressWarnings(as.numeric(source_data$LQ_mammal)),
  stringsAsFactors = FALSE
)
if (anyDuplicated(analysis_data$SpeciesBROAD)) {
  stop("Duplicated accepted_tree_tip values in the LQ dataset.")
}
if (anyNA(analysis_data$SpeciesBROAD) || any(!is.finite(analysis_data$LQ))) {
  stop("Species or untransformed LQ contains invalid values.")
}

result <- phylogenetic_shift_score(
  data = analysis_data,
  tree = tree_file,
  trait = "LQ",
  species_col = "SpeciesBROAD",
  outdir = results_dir,
  write_output = TRUE,
  verbose = TRUE
)

utils::write.table(
  result$model_fit,
  file.path(results_dir, "LQ.model_fit.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
utils::write.table(
  data.frame(
    trait = "LQ",
    source_column = "LQ_mammal",
    transformed = FALSE,
    matched_species = length(result$tree$tip.label),
    pair_count = nrow(result$scores),
    selected_model = result$selected_model,
    min_final_score = min(result$scores$FinalScore),
    max_final_score = max(result$scores$FinalScore),
    stringsAsFactors = FALSE
  ),
  file.path(results_dir, "LQ.analysis_summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
saveRDS(result, file.path(results_dir, "LQ.pss_result.rds"))
print(result)
