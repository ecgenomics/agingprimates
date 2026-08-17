# =============================
# Circular ultrametric tree + LQ bars
# Branches in plain color (no superfamily mapping)
# =============================

library(ape)
library(dplyr)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(scales)
library(phytools)

# ---- Input ----
tree <- read.tree("phylogeny.nw")
df   <- read.table("lqdf.tab", header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE)

# ---- Clean & standardize column names ----
clean_names <- function(x) {
  x <- gsub("[\r\n\t]", "", x)
  x <- trimws(x)
  x <- tolower(x)
  gsub("\\s+", "_", x)
}
names(df) <- clean_names(names(df))

# Canonicalize columns
df <- df %>%
  transmute(
    species     = species,
    lq          = suppressWarnings(as.numeric(lq)),
    superfamily = superfamily
  ) %>%
  filter(!is.na(lq), !is.na(species), species != "")

# ---- Prune tree & make ultrametric ----
species_in_both <- intersect(tree$tip.label, df$species)
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, species_in_both))
tree_ultra  <- phytools::force.ultrametric(tree_pruned, method = "extend")
tree_ultra$edge.length[tree_ultra$edge.length <= 0] <- 1e-6

# Align df to ultrametric tree
df_plot <- df %>%
  filter(species %in% tree_ultra$tip.label) %>%
  mutate(species = factor(species, levels = tree_ultra$tip.label))

# ---- Base circular tree (plain branches) ----
p <- ggtree(tree_ultra, layout = "circular", size = 0.4, color = "grey30") +
  theme(
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  )

# ---- Bars (LQ) around the tree ----
bars <- df_plot %>%
  select(species, lq, superfamily)

p <- p + ggnewscale::new_scale_fill()

p <- p +
  ggtreeExtra::geom_fruit(
    data        = bars,
    geom        = geom_bar,
    mapping     = aes(y = species, x = lq, fill = superfamily),
    orientation = "y",
    stat        = "identity",
    width       = 0.6,
    offset      = 0.02,
    inherit.aes = FALSE
  ) +
  scale_fill_discrete(name = "Superfamily (bars)") +
  labs(
    title    = "LQ along the phylogeny",
    subtitle = "Circular ultrametric tree with tip bars for LQ"
  )

# Optional tip labels
# p <- p + geom_tiplab(size = 2)

ggsave("LQ_circular_tree.png", p, width = 10, height = 10, dpi = 300)
print(p)

# =============================
# Densities of LQ (overall and by Superfamily)
# =============================

# Build a minimal dataframe for density plots
den_df <- df_plot %>%
  dplyr::select(lq, superfamily) %>%
  dplyr::filter(!is.na(lq))

# ---------- (1) Overall density in grey ----------
p_den_overall <- ggplot(den_df, aes(x = lq)) +
  geom_density(fill = "grey70", color = "grey30", alpha = 0.85, adjust = 1) +
  labs(
    title = "Density of LQ (overall)",
    x = "LQ",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.margin = margin(5, 5, 5, 5))

ggsave("LQ_density_overall.png", p_den_overall, width = 7, height = 5, dpi = 300)

# ---------- (2) Density by Superfamily (same hue palette as bars) ----------
# Define levels & palette consistent with bar plot ring
sup_levels <- sort(unique(den_df$superfamily))
den_df$superfamily <- factor(den_df$superfamily, levels = sup_levels)
sup_palette <- setNames(scales::hue_pal()(length(sup_levels)), sup_levels)

p_den_by_sup <- ggplot(den_df, aes(x = lq, fill = superfamily)) +
  geom_density(alpha = 0.4, color = NA, adjust = 1) +   # translucent fills, no outlines
  scale_fill_manual(values = sup_palette, name = "Superfamily") +
  labs(
    title = "Density of LQ by Superfamily",
    x = "LQ",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  )

ggsave("LQ_density_by_superfamily.png", p_den_by_sup, width = 7, height = 5, dpi = 300)