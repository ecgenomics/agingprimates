# Aging in primates

Reconstruction and extension of the comparative genomics workflow developed by
Farré et al. (2021), using the dense primate sampling and phylogeny released by
Kuderna et al. (2023).

The project asks which genes, evolutionary-rate shifts, and convergent
amino-acid substitutions are associated with unusually long or short maximum
lifespan in primates after accounting for body mass and shared evolutionary
history.

This repository is being rebuilt from the inputs upward. Historical analyses
are retained in `drawer/` for provenance, but they are not treated as current
results. Each active stage has its own inputs, scripts, outputs, checksums, and
methodological documentation.

## Scientific rationale

Maximum lifespan (MLS) is strongly associated with body mass across mammals.
Raw MLS therefore does not distinguish absolute lifespan from lifespan that is
unusually long or short for an animal's size. Our primary replication trait is
the mammal-wide longevity quotient (LQ):

```text
expected_MLS_years = 6.47 × body_mass_g^0.189
LQ_mammal = observed_MLS_years / expected_MLS_years
```

An LQ of 1 is the mammalian allometric expectation. Values above or below 1
represent longer- or shorter-than-expected lifespan, respectively. We retain
the fixed equation for comparability with the earlier mammalian analysis while
planning a primate-specific phylogenetic allometric model as a sensitivity
analysis.

The reconstructed workflow currently covers:

1. curation and provenance of phenotype and phylogeny inputs;
2. exact reproduction of the historical mammal-wide LQ;
3. phylogenetic phenotype-shift scoring;
4. estimation of phylogenetic signal;
5. preparation of foreground/background species comparisons and a Nextflow
   framework for genome–phenome analyses.

CAAS and RERconverge analyses will be reintroduced only after their inputs,
contrasts, missing-data rules, and multiple-testing strategy have been reviewed.

## Phylogeny

The sole operational primary tree is:

```text
inputs/science.abn7829_data_s4.nex.tree
```

This is an independent, byte-identical copy of Supplementary Data S4 from
Kuderna et al. (2023). Despite its filename extension, it contains plain Newick
text. The tree has 236 tips and is rooted, fully bifurcating,
fossil-calibrated, and ultrametric. Its root-to-tip distance is approximately
69.06116 time units, interpreted approximately as millions of years.

The original filename is deliberately retained to prevent confusion with the
former project tree. The MD5 checksum is:

```text
e6793ca5602c702b81ca255d329ee01b
```

## Repository structure

```text
.
├── inputs/                              approved and candidate source inputs
├── theory/                              numbered methodological notes
├── 02.allometric.transformation/        reproducible LQ construction
├── 03.phenotype.shift.detection.lq/     pairwise phylogenetic shift scoring
├── 04.phylognetic.signal.detection/     Pagel's λ and Blomberg's K
├── 05.gen.phen/                         genomic comparison preparation
└── drawer/                              archived historical work
```

The absence of a directory named `01` is intentional at present: input
assembly and project logic are documented in `inputs/` and `theory/` while the
independent input collection is still being curated.

### `inputs/`

Contains the selected Kuderna S4 tree, the current primate phenomic dataset,
trait metadata, genome/sample availability tables, candidate trait-selection
files, checksums, and an input-specific README.

Files in this directory are not automatically considered analytically
approved. In particular, maximum-lifespan and body-mass observations are being
audited at the record level for source, units, taxonomy, and missingness.

### `theory/`

Numbered notes record the scientific and statistical logic of the project:

- `01.project.outline.md`: project aims, workflow, outputs, and principles;
- `02.allometric.transformation.logic.md`: historical LQ implementation,
  input provenance, limitations, and proposed replication and sensitivity
  analyses.

New notes use two-digit increasing prefixes and dot-separated names.

### `02.allometric.transformation/`

`build_longevity_quotient.py` reconstructs the fixed-formula mammalian LQ from
the source phenomic data. It keeps phenotype validity separate from genome
availability and writes both an analysis table and a row-level audit table.

Current build summary:

- 526 source records;
- 162 records with valid positive MLS and body mass;
- 126 exact matches to tips in the Kuderna S4 tree;
- 140 valid phenotype records with genomic data available;
- exact reproduction of historical LQ values for all 140 overlapping species
  (maximum absolute difference approximately `8.88 × 10^-16`).

Run from this directory with:

```bash
python3 build_longevity_quotient.py
```

Principal outputs are `longevity_quotient.tsv`,
`longevity_quotient_audit.tsv`, `build_report.json`, and integrity checksums.

### `03.phenotype.shift.detection.lq/`

Applies the phylogenetic shift score (PSS) approach to the untransformed
`LQ_mammal` phenotype. The analysis includes 126 exactly matched species and
all 7,875 pairwise species comparisons. The OU model is currently selected
over Brownian motion for score construction.

The module produces complete score tables, model summaries, a top-1% nodal
visualization, and documented foreground/background selection tables. The
current strict absolute-LQ selection contains seven foreground species
(`LQ > 1.3`) and six background species (`LQ < 0.7`) drawn from species shared
between the independently defined top and bottom PSS tails. These thresholds
are a downstream contrast definition, not part of the allometric correction.

Run from the module directory with:

```bash
Rscript run_lq_pss.R
Rscript plot_lq_shifts_per_node.R
```

### `04.phylognetic.signal.detection/`

Measures phylogenetic signal in untransformed `LQ_mammal` for the 126 matched
species using `phytools::phylosig`.

Current estimates are:

| Statistic | Estimate | Test | P value |
|---|---:|---|---:|
| Blomberg's K | 0.145606 | 9,999 label randomizations | 0.00030003 |
| Pagel's λ | 0.785893 | likelihood-ratio test against λ = 0 | 3.94222 × 10^-9 |

Both tests detect significant phylogenetic signal. K below 1 indicates weaker
similarity among relatives than expected under Brownian motion, whereas λ
indicates substantial phylogenetic covariance. Full methods, unmatched taxa,
software versions, seed, and input checksums are recorded in `results.md`.

Run from the module directory with:

```bash
Rscript run_phylogenetic_signal.R
```

### `05.gen.phen/`

Contains the current handoff from phenotype definition to genome–phenome
analysis:

- `lq.table2.4species.comparisons/` generates explicit species-comparison
  manifests from the foreground/background selection;
- `lq.table2.nextflow/` provides Nextflow and Slurm entry points for the
  corresponding genomic workflow.

Large execution products, logs, and work directories are excluded from Git.

### `drawer/`

Contains the former analysis, datasets, provisional CAAS runs, top/bottom
contrasts, plots, and scripts. These files are preserved as a historical record
and as potential material for controlled recovery. Active analyses must not
silently depend on files in `drawer/`.

## Reproducibility principles

- Raw inputs, curated inputs, intermediate products, and results are kept
  conceptually separate.
- Species renaming, exclusion, duplication, and reconciliation must be
  explicit and auditable.
- Units and trait sources are recorded rather than inferred downstream.
- Genome availability is not used to define whether a phenotype is valid.
- The allometric equation, fitting population, and phylogeny are treated as
  part of the phenotype definition.
- Input copies and key outputs carry cryptographic checksums.
- Historical thresholds and gene lists are not accepted automatically.
- Site- and gene-level claims will require an explicit multiple-testing plan.

## Planned work

The next analytical steps are:

1. complete record-level auditing of MLS and adult body mass;
2. reconcile unresolved phenotype names to Kuderna S4 with a documented
   taxonomy table;
3. fit a primate-specific log–log allometric model using PGLS;
4. compare fixed mammalian LQ, log LQ, OLS residuals, and PGLS residuals;
5. repeat phylogenetic-signal and shift analyses across justified trait
   definitions;
6. freeze genomic foreground/background contrasts;
7. reconstruct CAAS and RERconverge analyses;
8. integrate resulting genes with Farré et al. (2021), known aging genes, and
   pathway enrichment.

## References

- Farré X, Molina R, Barteri F, Timmers PRHJ, Joshi PK. 2021. Comparative
  Analysis of Mammal Genomes Unveils Key Genomic Variability for Human Life
  Span. *Molecular Biology and Evolution* 38:4948–4961.
  https://doi.org/10.1093/molbev/msab219
- Kuderna LFK et al. 2023. A global catalog of whole-genome diversity from 233
  primate species. *Science* 380:906–913.
  https://doi.org/10.1126/science.abn7829
- de Magalhães JP, Costa J, Church GM. 2007. An analysis of the relationship
  between metabolism, developmental schedules, and longevity using
  phylogenetic independent contrasts. *The Journals of Gerontology: Series A*
  62:149–160. https://doi.org/10.1093/gerona/62.2.149
