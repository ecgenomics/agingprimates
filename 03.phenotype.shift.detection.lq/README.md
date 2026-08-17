# 03. Phenotype shift detection: LQ

This project runs Omar's analytical phylogenetic shift score implementation on
the **untransformed** longevity quotient.

- Phenotype exposed to `pss.core.R`: `LQ`
- Source column: `LQ_mammal`
- Explicitly unused column: `log_LQ_mammal`
- Phylogeny: Kuderna et al. S4 tree used by `omar/score.approach`
- Nodal figure: top 1% of PSS pairs, each assigned once to its MRCA, following
  `omar/score.approach/03.primate.traits/figures/shifts.per.node`

Run from this directory:

```bash
Rscript run_lq_pss.R
Rscript plot_lq_shifts_per_node.R
```

The pCloud filesystem returned `Function not implemented` for both symbolic
and hard links. Therefore `input/longevity_quotient.tsv` is a byte-identical
copy of `../02.allometric.transformation/longevity_quotient.tsv`; matching
SHA-256 checksums are recorded in `input/SHA256SUMS`.
