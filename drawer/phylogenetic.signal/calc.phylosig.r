# Load required libraries
library(ape)
library(phytools)
library(geiger)

# --- Input ---
# Read tree (replace with your own Newick if needed)
tree <- read.tree("phylogeny.nw")

# Read dataframe (tab-separated)
df <- read.table("lqdf.tab", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# --- Prepare trait vector ---
# Ensure LQ is numeric and Species are character
df$LQ <- as.numeric(df$LQ)
df$Species <- as.character(df$Species)

# Drop rows with missing LQ
df <- df[!is.na(df$LQ), ]

# Create named trait vector
trait <- df$LQ
names(trait) <- df$Species

# --- Make tree and data match ---
# Drop tree tips not in trait and align the trait to the (new) tree tip order
tree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[tree$tip.label]

# Optional sanity check
if (length(trait) < 3) {
  stop("Not enough species overlap between tree and data to fit models.")
}

# --- Phylogenetic signal ---
# Pagel’s lambda
lambda_res <- phylosig(tree, trait, method = "lambda", test = TRUE)
print(lambda_res)

# Blomberg’s K
K_res <- phylosig(tree, trait, method = "K", test = TRUE)
print(K_res)

# --- Geiger fitContinuous for different models ---
fitBM <- fitContinuous(tree, trait, model = "BM")   # Brownian Motion
fitOU <- fitContinuous(tree, trait, model = "OU")   # Ornstein-Uhlenbeck
fitEB <- fitContinuous(tree, trait, model = "EB")   # Early Burst

# --- Compare models: extract AICc manually (Method 1) ---
# Each fit object stores results in $opt (lnL, aic, aicc, etc.)
get_opt <- function(fit) {
  # Safely extract fields from geiger::fitContinuous result
  data.frame(
    logLik = if (!is.null(fit$opt$lnL)) fit$opt$lnL else NA_real_,
    AIC    = if (!is.null(fit$opt$aic)) fit$opt$aic else NA_real_,
    AICc   = if (!is.null(fit$opt$aicc)) fit$opt$aicc else NA_real_,
    stringsAsFactors = FALSE
  )
}

bm <- get_opt(fitBM); ou <- get_opt(fitOU); eb <- get_opt(fitEB)

aic_results <- data.frame(
  Model  = c("BM", "OU", "EB"),
  logLik = c(bm$logLik, ou$logLik, eb$logLik),
  AIC    = c(bm$AIC,    ou$AIC,    eb$AIC),
  AICc   = c(bm$AICc,   ou$AICc,   eb$AICc),
  stringsAsFactors = FALSE
)

# Compute ΔAICc and Akaike weights for readability
min_aicc <- min(aic_results$AICc, na.rm = TRUE)
aic_results$DeltaAICc <- aic_results$AICc - min_aicc

w <- exp(-0.5 * aic_results$DeltaAICc)
aic_results$AkaikeWeight <- w / sum(w, na.rm = TRUE)

# Order by AICc (best first) and print
aic_results <- aic_results[order(aic_results$AICc), ]
print(aic_results, row.names = FALSE)

# --- Export results ---

# Save AIC comparison table
write.table(
  aic_results,
  file = "aic_results.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Build phylogenetic signal summary
phylosig_results <- data.frame(
  Method = c("Pagel_lambda", "Blomberg_K"),
  Value  = c(lambda_res$lambda, K_res$K),
  P_value = c(lambda_res$P, K_res$P),
  stringsAsFactors = FALSE
)

# Save phylogenetic signal table
write.table(
  phylosig_results,
  file = "phylosig_results.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)