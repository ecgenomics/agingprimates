# Foreground and background species selection from shared PSS tails

## Analytical scope

Foreground (FG) and background (BG) species were selected from the 47 species
that occur in both independently defined extremes of the pairwise
phylogenetic shift score (PSS) distribution for longevity quotient. PSS was
calculated for all 7,875 species pairs using the original, untransformed
mammal-wide longevity quotient (`LQ_mammal`) and the OU evolutionary model
selected by AIC. The upper and lower tails each contain
`ceiling(7,875 × 0.01) = 79` pairs. Species occurring in at least one upper-tail
pair were intersected with species occurring in at least one lower-tail pair,
producing the shared-tail candidate set before phenotype thresholds were
applied.

## Threshold definitions

The thresholds are strict and operate on absolute LQ rather than on a
transformed value, residual rank, PSS score, or model-derived posterior
probability:

- **Foreground (FG):** `LQ_mammal > 1.3`.
- **Background (BG):** `LQ_mammal < 0.8`.
- **Unselected shared-tail species:** `0.8 ≤ LQ_mammal ≤ 1.3`.

Consequently, a hypothetical species with LQ exactly equal to 1.3 or 0.8 would
not be assigned to FG or BG. These boundaries are fixed classification rules;
they were not estimated from the kernel density curve and do not represent
formal significance thresholds.

## Selected species

Seven species satisfy the FG criterion: *Ateles geoffroyi*, *Lemur catta*,
*Eulemur fulvus*, *Aotus trivirgatus*, *Eulemur macaco*,
*Cercocebus torquatus*, and *Nomascus concolor*. Their absolute LQ values range
from 1.313 to 1.347. Six species satisfy the BG criterion:
*Trachypithecus pileatus*, *Macaca thibetana*, *Eulemur rubriventer*,
*Mico humeralifer*, *Macaca tonkeana*, and *Alouatta seniculus*. Their absolute
LQ values range from 0.672 to 0.746. The remaining 34 shared-tail species fall
inside the inclusive interval from 0.8 to 1.3 and are not assigned to either
contrast group.

## Interpretation and intended use

An LQ above 1 indicates longevity greater than the mammalian expectation for
body mass, whereas an LQ below 1 indicates longevity below that expectation.
The FG set therefore represents a stringent high-LQ extreme, while the BG set
represents a stringent low-LQ extreme. Restricting both groups to shared-tail
species preserves the preceding PSS-based screening context: every selected
species participates in at least one top-1% pair and at least one bottom-1%
pair. This requirement distinguishes the contrast from a simple thresholding
of all 126 species in the LQ dataset.

The resulting groups are suitable as explicit foreground and background labels
for downstream comparative analyses, provided those analyses account for
phylogenetic relatedness and unequal group sizes. Group membership is not
evidence that a lineage experienced a discrete evolutionary regime shift, and
the labels do not identify causality. The selected species are phylogenetically
distributed across several primate families; observations are therefore not
independent, but neither group is confined to a single clade.

## Reproducibility

`table2.R` reconstructs both 1% PSS tails, derives their shared-species
intersection, verifies PSS phenotype values against the `LQ_mammal` source
column, applies the strict FG and BG thresholds, counts each species'
participation in the two tails, and writes both `table2.tsv` and `table2.md`.
The TSV is the machine-readable dataset; the Markdown file is its formatted
human-readable representation.
