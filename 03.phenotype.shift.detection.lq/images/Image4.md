# Image 4. Phylogenetic distribution of species shared by the top and bottom PSS tails

## Description

Phylogenetic placement of the 47 primate species represented in both the
independently selected top 1% and bottom 1% of pairwise phylogenetic shift
scores for the untransformed longevity quotient. The Kuderna et al. S4 primate
tree is pruned directly to the shared-tail species set without changing its
topology or branch lengths.

Strict phenotype thresholds identify seven foreground species (FG;
LQ_mammal > 1.3) and six background species (BG; LQ_mammal < 0.8). FG tips,
bars, and names are highlighted in vermillion; BG elements are highlighted in
blue. Selected names are bold and carry explicit FG or BG prefixes. The 34
shared-tail species inside the inclusive interval from 0.8 to 1.3 are shown in
grey.

Bar length represents absolute LQ on a common linear scale from zero to the
maximum among the displayed species, and every tip label reports its
untransformed LQ value. The plot therefore shows the phylogenetic distribution
of the threshold-defined groups; it does not infer discrete evolutionary
regimes or ancestral states.

## Reproduction

Run `Rscript Image4.R` from this directory, or invoke the script by its full
path from any working directory. The script reconstructs the two 1% tails,
intersects their species sets, verifies absolute LQ values against
`LQ_mammal`, prunes the S4 tree, and writes `Image4.png` beside the script.
It then applies the same strict FG and BG thresholds recorded in `table2.tsv`.

## Output specifications

- Format: PNG
- Resolution: 300 dpi
- Dimensions: 12 × 13.5 inches
- Pixel dimensions: 3,600 × 4,050
- Colour space: RGB
- Background: white
