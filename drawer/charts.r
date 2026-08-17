# Load required libraries
library(ggplot2)
library(dplyr)
library(ggrepel)

# Definisci cutoff
cutoff_up <- 1.3
cutoff_down <- 0.7
# Calcola la mediana globale di LQ_zscore
mediana <- median(lqdf$LQ, na.rm = TRUE)

# Filtra solo le specie significativamente longeve o brevilongeve
outliers <- lqdf %>%
  filter(LQ >= cutoff_up & LQ >= LQ_mean_family | LQ <= cutoff_down & LQ <= LQ_mean_family)

# Crea il boxplot con etichette solo per gli outlier ±cutoff
ggplot(data = lqdf, aes(y = Family, x = LQ, fill = Family)) +
  geom_boxplot(outlier.shape = NA) +                       # disattiva i pallini outlier di default
  geom_jitter(width = 0.1, alpha = 0.4) +                  # mostra tutti i punti
  geom_text_repel(data = outliers,                        # etichette solo per outlier ±cutoff
                  aes(label = Species),
                  box.padding = 0.3,
                  max.overlaps = Inf,
                  min.segment.length = 0) +
  geom_vline(xintercept = mediana, linetype = "dashed", color = "red") +  # linea mediana verticale
  ylab("Family") +
  xlab("LQ Z-score") +
  theme_minimal()