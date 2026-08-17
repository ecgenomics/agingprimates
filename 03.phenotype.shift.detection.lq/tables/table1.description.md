# Table 1. Top-ranked LQ shift pairs connected to the bottom tail through shared species

Table 1 identifies high-scoring longevity-quotient comparisons that are linked
at the species level to low-scoring comparisons. Phylogenetic shift scores
(PSS) were calculated from the untransformed mammal-wide longevity quotient
(`LQ_mammal`) for all 7,875 pairwise combinations among 126 primate species.
The top and bottom tails were selected independently using the highest and
lowest `ceiling(7,875 × 0.01) = 79` `FinalScore` values, following the
percentile-tail convention used in the Omar PSS workflow.

Species represented in the top tail were intersected with species represented
in the bottom tail, yielding 47 shared-tail species. Following the strict
shared-species network definition used in the Omar workflow, a top-tail pair
was retained only when both of its endpoint species belonged to this
intersection. Under this definition, 34 of the 79 top-tail pairs were retained.

Rows remain ordered by decreasing PSS `FinalScore`, and “Top rank” refers to
the rank within the complete 79-pair upper tail. LQ columns report the original,
untransformed phenotype values. The two final columns report how many of the
79 bottom-tail pairs involve each endpoint species; both counts must therefore
be positive for every retained row. This table is descriptive: repeated species create
dependence among rows, percentile tails are not inferential thresholds, and
membership in both score extremes does not by itself identify a causal
evolutionary mechanism.
