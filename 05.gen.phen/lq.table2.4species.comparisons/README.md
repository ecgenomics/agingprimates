# LQ table 2: four-species comparison configurations

This directory contains random comparison configurations derived from
`03.phenotype.shift.detection.lq/tables/table2.tsv`.

## Design

The input contains 7 foreground (`FG`) and 6 background (`BG`) species.

- **CAAStools:** 100 unique configurations, each containing 4 randomly selected
  FG species and 4 randomly selected BG species.
- **RERconverge:** 35 unique configuration files, each containing 4 FG species
  and the complete fixed set of 6 BG species.

There are only `choose(7, 4) = 35` distinct FG quartets. Therefore, the 100
requested RERconverge comparisons are reduced to the 35 genuinely distinct
possibilities. Every FG quartet occurs exactly once.

All files are headerless, tab-separated, two-column manifests with species in
the first column and binary state in the second (`1 = FG`, `0 = BG`).

For RERconverge, downstream analysis must be restricted to the species in each
manifest if the six listed BG species are to be the only background.

## Contents

- `input/table2.tsv`: snapshot of the source table used for generation.
- `configurations/caas/`: `001.caas.cfg` through `100.caas.cfg`.
- `configurations/rerconverge/`: `001.rerc.cfg` through `035.rerc.cfg`.
- `tables/comparison_manifest.tsv`: one row per generated configuration.
- `tables/comparison_species_membership.tsv`: long-form species membership.
- `tables/generation_metadata.tsv`: seed, pool sizes and design summary.
- `scripts/generate_comparisons.R`: reproducible generator.

## Regeneration

From the `agingprimates` repository root:

```bash
Rscript 05.gen.phen/lq.table2.4species.comparisons/scripts/generate_comparisons.R 100 260811
```
