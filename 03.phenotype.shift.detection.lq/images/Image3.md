# Image 3. Absolute LQ distribution of species shared by the top and bottom PSS tails

## Description

Distribution of the original, untransformed mammal-wide longevity quotient
(`LQ_mammal`) among the 47 species represented in both the independently
selected top 1% and bottom 1% of phylogenetic shift score (PSS) pairs. The
upper panel follows the Omar shared-tail phenotype-distribution logic but uses
a Gaussian kernel density estimate in place of the histogram. The bandwidth is
selected reproducibly with R's bw.nrd0 rule, and one rug mark is shown per
species. The lower panel adds an explicit species-level map: every point is an
individual species, ordered by absolute LQ.

Strict phenotype thresholds identify seven foreground species (FG;
LQ_mammal > 1.3) and six background species (BG; LQ_mammal < 0.8). FG is
shown in vermillion and BG in blue; threshold regions and lines are repeated
across both panels. Selected rows have larger shaped points, tinted row
guides, and explicit FG or BG prefixes attached to the species names. The 34
shared-tail species within the inclusive interval from 0.8 to 1.3 are shown in
grey. PhyloPic silhouettes were not added because the labelled, colour- and
shape-coded species map identifies every selected taxon directly without
introducing external graphical dependencies.

The dashed vertical line marks the sample median (0.989). Absolute LQ ranges
from 0.672 to 1.347 across the shared-tail species. These values are not
log-transformed. The display is descriptive and does not treat species as
independent statistical replicates.

## Reproduction

Run `Rscript Image3.R` from this directory, or invoke the script by its full
path from any working directory. The script reconstructs the top and bottom 1%
tails from `LQ.score_results.tsv`, intersects their species sets, verifies the
phenotype values against `LQ_mammal`, and writes `Image3.png` beside the script.
It then applies the strict FG and BG thresholds used in `table2.tsv`.

## Output specifications

- Format: PNG
- Resolution: 300 dpi
- Dimensions: 12 × 13.5 inches
- Pixel dimensions: 3,600 × 4,050
- Colour space: RGB
- Background: white
