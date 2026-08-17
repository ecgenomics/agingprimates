#!/usr/bin/env Rscript

# Standalone, reproducible phylogenetic-signal analysis of untransformed LQ.
# The script uses only the input snapshots bundled beside it and can be run
# from any working directory.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (!length(file_arg)) {
  stop("Run this file with Rscript: Rscript run_phylogenetic_signal.R")
}
script_arg <- sub("^--file=", "", file_arg[[1]])
# Rscript encodes spaces as "~+~" in the --file argument on some systems.
script_arg <- gsub("~\\+~", " ", script_arg)
script_path <- normalizePath(script_arg, mustWork = TRUE)
analysis_dir <- dirname(script_path)

required_packages <- c("ape", "phytools")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing R package(s): ", paste(missing_packages, collapse = ", "),
    ". Install them with install.packages(c(",
    paste(sprintf("\"%s\"", missing_packages), collapse = ", "),
    "))."
  )
}

dataset_file <- file.path(analysis_dir, "input", "longevity_quotient.tsv")
tree_file <- file.path(analysis_dir, "input", "science.abn7829_data_s4.nex.tree")
output_file <- file.path(analysis_dir, "results.md")
stopifnot(file.exists(dataset_file), file.exists(tree_file))

# Analysis constants. The fixed seed makes the randomization test reproducible.
trait_column <- "LQ_mammal"
species_column <- "accepted_tree_tip"
n_permutations <- 9999L
random_seed <- 20260817L
alpha <- 0.05

data <- utils::read.delim(
  dataset_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "\""
)
required_columns <- c(species_column, trait_column, "taxonomic_mapping_status")
if (!all(required_columns %in% names(data))) {
  stop(
    "Dataset lacks required column(s): ",
    paste(setdiff(required_columns, names(data)), collapse = ", ")
  )
}

species <- as.character(data[[species_column]])
trait_values <- suppressWarnings(as.numeric(data[[trait_column]]))
if (anyNA(species) || any(!nzchar(species))) {
  stop("Species identifiers contain missing or empty values.")
}
if (anyDuplicated(species)) {
  stop("Species identifiers are duplicated in the LQ dataset.")
}
if (any(!is.finite(trait_values))) {
  stop("LQ_mammal contains missing or non-finite values.")
}
if (any(trait_values <= 0)) {
  stop("LQ_mammal must be strictly positive.")
}

tree <- ape::read.tree(tree_file)
if (is.null(tree$edge.length) || any(!is.finite(tree$edge.length))) {
  stop("The phylogeny must have finite branch lengths.")
}
if (any(tree$edge.length <= 0)) {
  stop("The phylogeny contains non-positive branch lengths.")
}
if (anyDuplicated(tree$tip.label)) {
  stop("The phylogeny contains duplicated tip labels.")
}
if (!ape::is.rooted(tree)) {
  stop("The phylogeny must be rooted.")
}
if (!ape::is.binary(tree)) {
  stop("The phylogeny must be fully bifurcating.")
}
if (!ape::is.ultrametric(tree)) {
  stop("The phylogeny must be ultrametric.")
}

trait <- stats::setNames(trait_values, species)
matched_species <- tree$tip.label[tree$tip.label %in% names(trait)]
excluded_species <- setdiff(names(trait), tree$tip.label)
tree_only_species <- setdiff(tree$tip.label, names(trait))
if (length(matched_species) < 3L) {
  stop("Fewer than three species overlap between the LQ dataset and tree.")
}

pruned_tree <- ape::keep.tip(tree, matched_species)
trait <- trait[pruned_tree$tip.label]
stopifnot(identical(names(trait), pruned_tree$tip.label))

set.seed(random_seed)
blomberg <- phytools::phylosig(
  pruned_tree,
  trait,
  method = "K",
  test = TRUE,
  nsim = n_permutations
)
pagel <- phytools::phylosig(
  pruned_tree,
  trait,
  method = "lambda",
  test = TRUE
)

k_value <- unname(blomberg$K)
k_p <- unname(blomberg$P)
lambda_value <- unname(pagel$lambda)
lambda_p <- unname(pagel$P)
lambda_lr <- 2 * (unname(pagel$logL) - unname(pagel$logL0))

fmt_estimate <- function(x) formatC(x, digits = 6, format = "fg", flag = "#")
fmt_p <- function(x) {
  if (x < 0.001) {
    formatC(x, digits = 4, format = "e")
  } else {
    formatC(x, digits = 4, format = "f")
  }
}
yes_no <- function(x) if (isTRUE(x)) "sì" else "no"

k_strength <- if (k_value < 1) {
  "inferiore a 1, quindi la somiglianza tra parenti è più debole di quella attesa sotto moto browniano"
} else if (k_value > 1) {
  "superiore a 1, quindi la somiglianza tra parenti è più forte di quella attesa sotto moto browniano"
} else {
  "uguale a 1, in linea con l'attesa sotto moto browniano"
}
k_evidence <- if (k_p < alpha) {
  "Il test di randomizzazione rifiuta l'assenza di segnale filogenetico."
} else {
  "Il test di randomizzazione non rifiuta l'assenza di segnale filogenetico."
}
lambda_evidence <- if (lambda_p < alpha) {
  "Il likelihood-ratio test rifiuta λ = 0."
} else {
  "Il likelihood-ratio test non rifiuta λ = 0."
}
overall_interpretation <- if (k_p < alpha && lambda_p < alpha) {
  "Entrambi i test rilevano un segnale filogenetico significativo nel LQ."
} else if (k_p < alpha || lambda_p < alpha) {
  paste(
    "I due test danno evidenza statistica diversa.",
    "Questo non è necessariamente contraddittorio perché K e λ misurano aspetti diversi della struttura filogenetica."
  )
} else {
  "Nessuno dei due test rileva un segnale filogenetico significativo nel LQ."
}

md5 <- tools::md5sum(c(dataset_file, tree_file))
package_versions <- vapply(
  required_packages,
  function(pkg) as.character(utils::packageVersion(pkg)),
  character(1)
)

report <- c(
  "# 04. Phylogenetic signal detection: LQ",
  "",
  paste0("Report generato il ", Sys.Date(), "."),
  "",
  "## Risultati",
  "",
  "| Statistica | Stima | Test dell'assenza di segnale | P-value | Significativo (α = 0.05) |",
  "|---|---:|---|---:|:---:|",
  sprintf(
    "| K di Blomberg | %s | Randomizzazione delle etichette (%d permutazioni) | %s | %s |",
    fmt_estimate(k_value), n_permutations, fmt_p(k_p), yes_no(k_p < alpha)
  ),
  sprintf(
    "| λ di Pagel | %s | Likelihood-ratio test, H0: λ = 0 (LR = %s) | %s | %s |",
    fmt_estimate(lambda_value), fmt_estimate(lambda_lr), fmt_p(lambda_p),
    yes_no(lambda_p < alpha)
  ),
  "",
  "## Interpretazione",
  "",
  paste0("- K è ", fmt_estimate(k_value), ": ", k_strength, ". ", k_evidence),
  paste0(
    "- λ è ", fmt_estimate(lambda_value),
    "; 0 indica assenza di covarianza filogenetica e valori vicini a 1 una struttura simile all'attesa browniana. ",
    lambda_evidence
  ),
  paste0("- **Conclusione:** ", overall_interpretation),
  "",
  "## Dati e controlli",
  "",
  sprintf("- Fenotipo: `%s` (LQ non trasformato).", trait_column),
  sprintf("- Specie nel dataset: %d.", nrow(data)),
  sprintf("- Specie nell'albero completo: %d.", ape::Ntip(tree)),
  sprintf("- Specie analizzate dopo l'intersezione: %d.", length(trait)),
  sprintf("- Specie del dataset escluse perché assenti dall'albero: %d.", length(excluded_species)),
  sprintf("- Punte dell'albero escluse perché prive di LQ: %d.", length(tree_only_species)),
  paste0(
    "- Albero analizzato: radicato = ", yes_no(ape::is.rooted(pruned_tree)),
    "; biforcante = ", yes_no(ape::is.binary(pruned_tree)),
    "; ultrametrico = ", yes_no(ape::is.ultrametric(pruned_tree)), "."
  ),
  "",
  "Le specie senza corrispondenza esatta nell'albero sono:",
  "",
  paste0("`", paste(excluded_species, collapse = "`, `"), "`."),
  "",
  "## Metodo",
  "",
  paste0(
    "I nomi delle specie sono stati abbinati esattamente tra `accepted_tree_tip` e le punte dell'albero. ",
    "L'albero è stato potato alle specie con LQ disponibile. K di Blomberg è stato testato ",
    "permutando i valori del tratto sulle punte; λ di Pagel è stato stimato per massima verosimiglianza ",
    "e confrontato con λ = 0 mediante likelihood-ratio test. Entrambe le statistiche sono state calcolate ",
    "con `phytools::phylosig`."
  ),
  "",
  "## Riproducibilità",
  "",
  "Dalla cartella dell'analisi:",
  "",
  "```bash",
  "Rscript run_phylogenetic_signal.R",
  "```",
  "",
  sprintf("- Seed del test di randomizzazione: `%d`.", random_seed),
  sprintf("- R: `%s`.", R.version.string),
  sprintf("- ape: `%s`; phytools: `%s`.", package_versions[["ape"]], package_versions[["phytools"]]),
  sprintf("- MD5 `input/longevity_quotient.tsv`: `%s`.", unname(md5[[1]])),
  sprintf("- MD5 `input/science.abn7829_data_s4.nex.tree`: `%s`.", unname(md5[[2]])),
  "",
  "Gli input sono copie locali incluse nella cartella, quindi lo script non dipende da file esterni al modulo di analisi."
)

writeLines(report, output_file, useBytes = TRUE)
message("Wrote ", output_file)
message(
  "Matched species: ", length(trait),
  "; Blomberg K = ", signif(k_value, 6),
  " (P = ", signif(k_p, 6), ")",
  "; Pagel lambda = ", signif(lambda_value, 6),
  " (P = ", signif(lambda_p, 6), ")"
)
