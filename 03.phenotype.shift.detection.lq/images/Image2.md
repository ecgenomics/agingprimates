# Image 2. Nodal localization of the bottom 1% of LQ PSS pairs

## Description

Circular phylogeny showing where the lowest-scoring 1% of pairwise
phylogenetic shift score (PSS) comparisons occur for the untransformed
longevity quotient (`LQ_mammal`). The analysis includes 126 species and all
7,875 possible species pairs; the displayed tail contains 79 pairs.

Each species pair is assigned exactly once to its most recent common ancestor
(MRCA) on the trait-pruned Kuderna et al. S4 primate tree. Internal-node size
represents the number of bottom-tail pairs assigned to that node. Internal-node
colour represents the corresponding nodal rate: the assigned bottom-tail count
divided by the number of all possible pairs whose MRCA is that node.

Terminal bars show absolute, untransformed LQ values on a common linear scale
from zero to the maximum observed value. Bars and terminal points are coloured
by family. The evolutionary model selected by PSS was OU.

## Reproduction

Run `Rscript Image2.R` from this directory, or invoke the script by its full
path from any working directory. The script locates the project inputs and PSS
results relative to its own path and writes `Image2.png` beside the script.

## Output specifications

- Format: PNG
- Resolution: 300 dpi
- Dimensions: 12 × 10.2 inches
- Colour space: RGB
- Background: white
