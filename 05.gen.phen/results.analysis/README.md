# CAAStools bootstrap result analysis

## Status and scope

This directory contains the first, script-level implementation of quality
control and filtering for the CAAStools bootstrap analysis of primate longevity
comparisons. The code currently operates on individual headerless CAAStools
bootstrap result files, plus a dataset-wide reference table constructed from a
directory of those files. A Nextflow workflow for applying the analysis to all
genes has **not** been created yet.

The implementation was developed and validated on 18 August 2026 using the
16,128 result files currently available under:

```text
../lq.table2.bootstrap.nextflow/results.fromcluster/results/
└── run-260817235905/
    └── caas-bootstrap/
```

The original Nextflow execution contained 16,133 gene tasks. Therefore, the
cycle-burden reference must be rebuilt if the five currently absent results are
added or if any input result is replaced.

No raw CAAStools result is modified in place. All filtering produces new files
and explicit audit reports.

## Directory contents

```text
results.analysis/
├── README.md
├── examples/
│   └── bootstrap.example.tsv
├── scripts/
│   ├── build_cycle_burden_reference.py
│   └── filter_bootstrap_caas.py
└── tests/
    ├── test_build_cycle_burden_reference.py
    └── test_filter_bootstrap_caas.py
```

The scripts depend only on the Python standard library. No pandas, NumPy, or
SciPy installation is required for this analysis stage.

## Biological and technical motivation

The bootstrap output contains one row for every alignment position that reached
the CAAStools bootstrap stage. Most positions have no positive bootstrap cycle
and are not useful downstream. Among the positive positions, two additional
patterns can be caused by poor alignment quality, repetitive sequence, or a
single species configuration repeatedly matching the same amino-acid pattern:

1. **Local dense clusters.** Many consecutive or nearly consecutive positions
   are CAAS-positive, and the same bootstrap cycle recurs through most of the
   cluster. Long runs of this kind are less credible as independent molecular
   convergence events and can indicate a problematic alignment segment.

2. **Gene-wide cycle dominance.** A single bootstrap cycle can recur at many
   dispersed positions throughout a gene even after local dense clusters are
   removed. `THADA` is the motivating example: cycle `b_74` originally occurred
   in 55 of the gene's 65 positive positions.

Local clustering and gene-wide dominance are distinct phenomena. The complete
filter therefore applies a dataset-wide cycle-burden rule first and a local
interval rule second.

## CAAStools bootstrap input format

Each input is a headerless, tab-separated file. The downloaded legacy run has
six columns; the updated CAAStools workflow writes a seventh positional p-value
column immediately after the cycle identifiers:

| Column | Meaning |
|---|---|
| 1 | Position identifier in `gene@position` form |
| 2 | Number of bootstrap cycles positive at the position |
| 3 | Total number of bootstrap cycles tested |
| 4 | Positive-cycle frequency, equal to column 2 / column 3 |
| 5 | Comma-separated positive cycle identifiers such as `b_47,b_66` |
| 6 | Hypergeometric p-value in new outputs; trait-template filename in legacy outputs |
| 7 | Trait-template filename in new outputs |

Example:

```text
THADA@220    1    100    0.01    b_74    0.0125    longevity.template.001.caas.cfg
```

The actual separator is a tab, not spaces.

### Important format details

- CAAStools alignment coordinates are **zero-based**. Position `gene@0` is
  valid.
- A position may be positive for more than one cycle. Consequently, the sum of
  cycle occurrences can be greater than the number of positive positions.
- The number of cycle identifiers in column 5 must equal column 2.
- Column 4 must equal column 2 divided by column 3.
- When present, the positional hypergeometric p-value must be numeric and
  between zero and one.
- Every input file must contain a single gene and no duplicate coordinate.
- The scripts accept only canonical cycle identifiers matching `b_NUMBER`.
- Blank lines are ignored; malformed nonblank lines terminate the analysis with
  an explanatory error.

"Total positions" in this analysis means the number of position rows present
in the CAAStools result. It is not necessarily the full unfiltered protein
alignment length. The upstream CAAStools run already used the hypergeometric
significance filter, so positions rejected before bootstrap are absent from the
file.

## Terminology and notation

For a gene \(g\) and cycle \(c\):

- \(T_g\): number of input position rows for the gene;
- \(P_g\): number of rows with at least one positive cycle;
- \(H_{g,c}\): number of positions positive for cycle \(c\);
- cycle dominance: \(H_{g,c} / P_g\);
- cycle burden: \(H_{g,c} / T_g\).

The distinction between dominance and burden is important. A cycle can account
for 100% of the positive positions in a gene with only one positive position;
that case should not be treated like a cycle occurring at dozens or hundreds
of positions.

## Analysis overview

The full analysis has two executable phases.

### Phase 1: build a genomic cycle-burden reference

Script:

```text
scripts/build_cycle_burden_reference.py
```

For every input gene, the script counts the positions associated with each
cycle and calculates the cycle burden \(H_{g,c}/T_g\). It then constructs a
separate empirical distribution for each cycle.

Only genes in which the cycle occurs contribute to that cycle's conditional
distribution. The default cutoff is the 97.5th percentile of this nonzero
distribution. Quantiles are calculated using linear interpolation with rank:

```text
(number_of_values - 1) × quantile
```

This is the commonly used type-7 empirical quantile definition.

The reference is cycle-specific because the 100 longevity comparisons are
correlated and do not have identical baseline behaviour. A universal burden
cutoff would over-penalize cycles that are intrinsically more frequent and
under-penalize cycles that are usually rare.

#### Reference command

```bash
python3 scripts/build_cycle_burden_reference.py \
  ../lq.table2.bootstrap.nextflow/results.fromcluster/results/run-260817235905/caas-bootstrap \
  --output cycle-burden-reference.tsv \
  --quantile 0.975
```

The default file pattern is:

```text
*.bootstrap.caas.tsv
```

It can be changed with `--pattern`. Patterns are interpreted relative to the
input directory.

#### Reference output columns

| Column | Meaning |
|---|---|
| `cycle_id` | Bootstrap cycle identifier |
| `quantile` | Requested conditional quantile |
| `burden_threshold` | Empirical burden cutoff for the cycle |
| `genes_with_cycle` | Genes contributing a nonzero burden |
| `median_burden` | Conditional median burden |
| `maximum_burden` | Largest observed burden |
| `input_files` | Number of result files used |
| `genes_with_positive_cycles` | Input genes with at least one positive cycle |

In the current 16,128-file reference, `b_74` has:

```text
quantile                0.975
burden_threshold        0.307428371068
genes_with_cycle        9282
median_burden           0.0177909660402
maximum_burden          1
input_files             16128
genes_with_positive     13060
```

### Phase 2: filter one gene

Script:

```text
scripts/filter_bootstrap_caas.py
```

When a cycle-burden reference is supplied, the operations occur in this exact
order:

1. validate the complete CAAStools result;
2. detect globally anomalous gene-cycle associations;
3. remove only the flagged cycle identifiers from affected rows;
4. recalculate the positive-cycle count and frequency in every changed row;
5. remove original zero-cycle rows and rows emptied by global cycle removal;
6. detect and remove local dense recurrent-cycle intervals;
7. write the filtered CAAStools rows and both audit reports.

### Global gene-cycle rule

A cycle is flagged in a gene only when all three default conditions hold:

```text
H(g,c) >= 10
H(g,c) / P(g) >= 0.80
H(g,c) / T(g) >= cycle-specific 97.5th-percentile burden cutoff
```

The defaults are configurable:

```text
--min-cycle-positions 10
--min-global-cycle-dominance 0.8
```

The quantile itself is determined when the reference is built.

Every cycle satisfying the rule is removed, not only the single cycle with the
largest count. This matters because correlated longevity configurations can be
positive together at nearly every position. A gene can therefore contribute
more than one flagged gene-cycle pair.

#### Cycle removal is position-preserving when possible

The global filter removes a flagged cycle identifier rather than blindly
discarding the complete row. For example:

```text
before: b_47,b_74,b_66    count=3    frequency=0.03
flag:   b_74
after:  b_47,b_66         count=2    frequency=0.02
```

The row is deleted only if no unflagged cycle remains. This preserves evidence
for independent configurations and is deliberately less destructive than
removing every position touched by a globally recurrent cycle.

### Local dense-cluster rule

After global pruning, all pairs of remaining positive positions are considered
as possible interval endpoints. For start coordinate \(s\) and end coordinate
\(e\):

```text
inclusive span L = e - s + 1
K = number of positive CAAS positions from s through e
N(c) = number of those K positions containing cycle c
```

The interval is suspicious when:

```text
L >= 10 amino acids
K / L >= 0.80
max_c N(c) / K >= 0.80
```

The defaults are configurable:

```text
--min-span 10
--min-density 0.8
--min-cycle-recurrence 0.8
```

The density denominator is the inclusive coordinate span, not the number of
rows present in the input. Thus eight positive positions distributed across
coordinates 100–109 have density 8/10 and pass exactly at the default cutoff.

Dense intervals supported by different cycles are retained. For example, ten
consecutive positive positions carrying ten different cycles have CAAS density
1.0 but no cycle recurrence of 0.8, so they are not classified as an alignment
artifact by this rule.

For every possible start position, the furthest qualifying endpoint is kept.
Intervals wholly contained in a previously reported interval are omitted from
the report because they add no removed position. Partially overlapping maximal
intervals are retained as separate report rows. Positions covered by multiple
reported intervals are removed and counted only once.

## Complete per-gene command

Using explicit output paths is recommended for a future workflow:

```bash
python3 scripts/filter_bootstrap_caas.py \
  /path/to/THADA.Homo_sapiens.filter2.bootstrap.caas.tsv \
  --cycle-burden-reference cycle-burden-reference.tsv \
  --output filtered/THADA.bootstrap.caas.filtered.tsv \
  --flagged-cycle-report reports/THADA.flagged-cycles.tsv \
  --cluster-report reports/THADA.clusters.tsv
```

If explicit paths are omitted, the script inserts tags before the `.tsv`
suffix:

```text
gene.bootstrap.caas.filtered.tsv
gene.bootstrap.caas.flagged-cycles.tsv
gene.bootstrap.caas.clusters.tsv
```

The input is never overwritten. Output and report paths are required to be
distinct.

### Local-only mode

If `--cycle-burden-reference` is omitted, the script preserves its earlier
behaviour and applies only zero-cycle and local-cluster filtering:

```bash
python3 scripts/filter_bootstrap_caas.py gene.bootstrap.caas.tsv
```

This mode is useful for testing the local rule independently, but the complete
analysis should use the genomic reference.

## Output files

### Filtered CAAStools result

The filtered result remains headerless and preserves its original six-column
(legacy) or seven-column (current) layout. When global cycles are removed,
columns 2, 4, and 5 are recalculated:

- column 2: number of retained cycles;
- column 4: retained cycles divided by the original total cycles;
- column 5: retained comma-separated cycle identifiers.

The total cycle count in column 3 remains 100 because the original analysis
tested 100 configurations.

### Flagged-cycle report

The global report contains one row per removed gene-cycle association:

| Column | Meaning |
|---|---|
| `gene` | Gene identifier |
| `cycle_id` | Removed cycle |
| `cycle_positions` | Positions associated with that cycle |
| `positive_positions_before_filter` | All positive positions in the raw gene |
| `total_positions` | All input rows for the gene |
| `cycle_dominance` | Cycle positions / positive positions |
| `cycle_burden` | Cycle positions / total positions |
| `reference_quantile` | Quantile used to build the cutoff |
| `burden_threshold` | Cycle-specific genomic cutoff |
| `positions_affected` | Comma-separated alignment coordinates |

An empty report still contains its header.

### Local-cluster report

The local report contains one row per maximal suspicious interval:

| Column | Meaning |
|---|---|
| `cluster_id` | File-local cluster identifier |
| `gene` | Gene identifier |
| `start_position` | First positive endpoint |
| `end_position` | Last positive endpoint |
| `span_aa` | Inclusive coordinate length |
| `caas_positions` | Positive positions in the interval |
| `caas_density` | Positive positions / span |
| `recurrent_cycles` | Cycles meeting the recurrence cutoff |
| `maximum_cycle_occurrences` | Largest cycle count in the interval |
| `maximum_cycle_fraction` | Largest cycle count / positive positions |
| `positions_removed` | Comma-separated coordinates removed |

An empty report still contains its header.

## Worked examples

### Synthetic example

Input:

```text
examples/bootstrap.example.tsv
```

The example contains 34 rows:

- positions 1–10 are all positive and `b_1` occurs in exactly 8/10;
- positions 11–15 have zero positive cycles;
- positions 20–29 are all positive but carry ten different cycles;
- positions 40–48 all carry `b_3`, but the interval is only 9 aa long.

With local defaults:

- positions 1–10 are removed as a recurrent dense cluster;
- positions 11–15 are removed because they have zero cycles;
- positions 20–29 are retained because no cycle dominates;
- positions 40–48 are retained because the interval is shorter than 10 aa.

Observed summary:

```text
Total positions: 34
Zero-cycle positions removed: 5
Positive positions before cluster filtering: 29
Suspicious intervals: 1
Positive positions removed in suspicious intervals: 10
Positions retained: 19
```

### THADA

Raw `THADA` result:

```text
total input rows                         153
zero-cycle rows                           88
positive positions                        65
distinct positive cycles                  40
cycle-position associations              141
```

Before filtering, `b_74` occurred in 55/65 positive positions:

```text
b_74 dominance     55 / 65  = 0.846154
b_74 burden        55 / 153 = 0.359477
b_74 cutoff                  = 0.307428
```

The `b_74` burden is at approximately the 98.6th percentile among the 9,282
genes in which `b_74` occurs, so it exceeds the selected 97.5th-percentile
cutoff.

Complete filtering gives:

```text
globally flagged cycles                         1  (b_74)
cycle-position associations removed            55
positions emptied by global filtering          45
positive positions after global filtering      20
local suspicious intervals                      0
final retained positions                       20
```

Ten positions originally carried `b_74` together with other cycles. Those
positions were retained after removing `b_74`, with their counts and
frequencies recalculated. In particular, the secondary signals remain:

```text
b_47    16 retained positions
b_66    14 retained positions
```

All 20 retained THADA rows were checked: their cycle count, comma-separated
cycle list, and recalculated frequency agree exactly.

The species assignment for `b_74` is:

```text
FG: Eulemur_fulvus,Eulemur_macaco,Lemur_catta,Nomascus_concolor
BG: Alouatta_seniculus,Macaca_tonkeana,Mico_humeralifer,Trachypithecus_pileatus
```

### C4A local-filter example

`C4A` contains 467 input rows:

```text
zero-cycle positions        299
positive positions          168
local suspicious intervals    0
positions retained          168
```

This demonstrates that a gene can have many positive positions without being
removed by the local filter when they do not form a sufficiently dense,
same-cycle interval.

## Threshold calibration

A global dominance cutoff alone is too broad for these correlated bootstrap
comparisons. Before implementing the cycle-specific burden reference, the
following raw-data behaviour was observed:

```text
dominant-cycle share >= 0.80 and hits >= 10: 2,603 genes
```

Adding a cycle-specific empirical burden cutoff gave the following approximate
impact when considering the single most frequent cycle per gene:

| Conditional burden cutoff | Genes flagged | THADA flagged? |
|---:|---:|:---:|
| 95th percentile | 1,164 | yes |
| 97.5th percentile | 580 | yes |
| 99th percentile | 218 | no |

The 97.5th percentile was selected as a compromise. The production
implementation evaluates every cycle satisfying the rule, not only one
preselected dominant cycle. Therefore, the final number of genes containing at
least one flagged cycle is 648 rather than 580, and some genes contribute
multiple flagged gene-cycle pairs.

These cutoffs are quality-control heuristics, not formal hypothesis-test
p-values. They should remain configurable and their reports should be retained
for sensitivity analyses.

## Validation on the current dataset

### Raw inputs

| Metric | Value |
|---|---:|
| Bootstrap result files | 16,128 |
| Input position rows | 2,275,918 |
| Positive positions | 273,201 |
| Zero-cycle positions | 2,002,717 |
| Cycle-position associations | 1,198,696 |
| Genes with at least one positive cycle | 13,060 |
| Cycles represented | 100 |

### Local-only filter

Before global pruning, the local rule identified:

| Metric | Value |
|---|---:|
| Genes with at least one local cluster | 596 |
| Maximal suspicious intervals | 1,390 |
| Unique positive positions removed | 14,470 |
| Largest positive-position count in one gene | 1,240 (`FSIP2`) |
| Largest interval count in one gene | 28 (`ZNF845`) |

The complete local-only validation required approximately 39 seconds on the
development machine and did not write bulk output files.

### Global plus local filter

Using the 97.5th-percentile reference and default per-gene thresholds:

| Metric | Value |
|---|---:|
| Files validated | 16,128 |
| Genes with at least one globally flagged cycle | 648 |
| Flagged gene-cycle pairs | 1,952 |
| Cycle-position associations removed | 148,535 |
| Positions emptied by global cycle removal | 24,556 |
| Local intervals remaining after global pruning | 854 |
| Positions removed by the subsequent local filter | 8,882 |
| Final retained positions | 239,763 |

The arithmetic is internally consistent:

```text
273,201 raw positive positions
-24,556 emptied by global cycle removal
- 8,882 removed by local intervals
=239,763 final retained positions
```

The complete in-memory global-plus-local validation required approximately 35
seconds after the reference had been built. Every retained row was checked for
agreement between its count, cycle list, and frequency.

## Automated tests

Run all tests from this directory:

```bash
python3 -m unittest discover -s tests -p 'test_*.py' -v
```

Current status:

```text
17 tests passed
```

The tests cover:

- parsing and strict validation of legacy six-column and current seven-column
  files;
- acceptance of zero-based coordinates;
- rejection of mismatched cycle counts;
- exact 80% local-density and recurrence boundaries;
- nearly consecutive positions with gaps;
- retention of dense clusters supported by different cycles;
- retention of recurrent intervals shorter than 10 aa;
- overlapping intervals and unique removal counts;
- configurable local thresholds;
- output naming and row preservation;
- empirical quantile interpolation;
- construction and reloading of cycle-specific references;
- minimum global occurrence and dominance conditions;
- removal of only the globally flagged cycle;
- preservation of alternative cycles on the same position;
- recalculation of counts and frequencies;
- flagged-cycle and local-cluster report generation;
- failure when a reference lacks an observed cycle.

Syntax compilation can be checked with:

```bash
python3 -m py_compile scripts/*.py tests/*.py
```

## Performance and implementation notes

- Reference construction is linear in the number of input rows plus listed
  cycle associations.
- Global per-gene filtering is linear in the number of rows and cycle
  associations for that gene.
- Local interval detection follows the requested pairwise-position logic and
  is quadratic in the number of positive positions remaining in a gene.
- Applying global pruning before local interval detection reduces the number of
  positions entering the quadratic step.
- On the current dataset, the largest raw positive-position count is 1,240, and
  the complete validation remains practical on a local machine.
- The scripts write one gene at a time and do not load the full set of alignment
  rows into memory simultaneously. Reference construction retains only the
  per-cycle burden distributions.

## Known limitations and interpretation cautions

1. **The reference is conditional on cycle occurrence.** Genes with zero hits
   for a cycle are not included in that cycle's empirical distribution.

2. **The reference includes the genes later classified.** It is not a
   leave-one-gene-out reference. With thousands of observations per common
   cycle, the influence of one gene on its own cutoff is expected to be small,
   but leave-one-out or held-out references could be implemented later.

3. **Bootstrap cycles are correlated.** The 100 configurations reuse species,
   so cycle IDs are not independent replicates. This is why all cycle removals
   are reported and why a cycle-specific empirical reference is preferred over
   a universal theoretical cutoff.

4. **Input rows are already prefiltered.** Burden is normalized by positions
   present in the bootstrap output, not by the unfiltered protein length.
   Comparisons involving results generated with different upstream filters
   would require additional care.

5. **The global and local rules are quality-control heuristics.** A genuinely
   convergent gene could contain multiple substitutions associated with one
   longevity configuration. Flagged events should remain auditable rather than
   being permanently erased from the raw data.

6. **Alignment quality is inferred indirectly.** Neither filter inspects the
   amino-acid alignment itself. They identify output patterns consistent with
   alignment problems but do not prove that an alignment is incorrect.

7. **Current reference completeness.** The validation reference used 16,128
   files, whereas the upstream Nextflow run contained 16,133 tasks. Rebuild the
   reference when the final result collection changes.

8. **Partially overlapping local intervals remain separate report entries.**
   Their shared positions are removed only once. This retains the evidence for
   each qualifying pairwise interval without inflating removal counts.

9. **The current validation results predate the minimum-coverage rule.** The
   16,128 downloaded bootstrap tables were produced without requiring three
   observed species per FG and BG group. The upstream bootstrap workflow now
   requires at least 3/4 present, non-gap species in both groups for a cycle to
   count at a position. These existing tables remain useful for developing the
   downstream filters, but their cycle-burden reference must not be reused for
   the new cluster run. Rebuild the reference from the complete rerun.

## Recommended workflow design for the next stage

A future Nextflow implementation should:

1. treat the original bootstrap directory as read-only;
2. validate that the expected result set is complete;
3. build one versioned cycle-burden reference from the complete input set;
4. fan out one filtering task per gene using that immutable reference;
5. publish filtered rows, flagged-cycle reports, and cluster reports separately;
6. aggregate per-gene and per-cycle summaries;
7. retain genes with empty final output as valid completed results;
8. record the exact thresholds and input-file count in run metadata;
9. allow `-resume` without regenerating completed per-gene outputs when the
   reference has not changed;
10. keep all derived result directories outside Git.

Useful aggregate tables would include:

- gene × cycle position counts before and after filtering;
- dominant-cycle burden and empirical percentile per gene;
- every globally flagged gene-cycle pair;
- every locally filtered interval;
- counts of raw, globally retained, locally retained, and final positions per
  gene;
- mapping from cycle IDs to their foreground and background species.

## Development history

### 19 August 2026

- replaced the interpretation of bootstrap-cycle counts as independent
  evidence with an optional species- and amino-acid-aware event summary in the
  bundled local CAAStools copy;
- retained `b_N` identifiers for traceability but refer to them as hypotheses
  in the event output;
- inferred the complete FG and BG discovery pools from the resampling table and
  required stable group membership when species aggregation is requested;
- recorded the non-gap amino acid of every observed discovery species at each
  positive position;
- collapsed identical and compatible amino-acid signatures into maximal
  compatible events without an order-dependent greedy merge;
- replaced the list of input hypothesis patterns with one `event_pattern`
  calculated from the final merged signature, including pattern 4 for a
  many-vs-many amino-acid result;
- kept incompatible signatures as distinct events at the same position and
  ranked a primary event using balanced discovery support, total support, and
  nominal p-value;
- counted every supporting species once per event and emitted FG/BG support
  over both observed and complete discovery denominators;
- reported all-support and event-support amino-acid counts plus species whose
  amino acids conflict with each event signature;
- added `fg_dominant_amino_acid` and `bg_dominant_amino_acid` immediately after
  the corresponding event amino-acid sets; proportions use unique supporting
  species and preserve all co-dominant residues in a tie;
- introduced `ct pooled-discovery`, which accepts one complete fixed FG/BG
  pool file, generates unique subset hypotheses internally, saves the realized
  resampling-format table, and invokes the event summarizer without requiring
  an externally built 100-row input;
- added `--comparisons max|N` with `max` as the default, four-species FG/BG
  subset-size options, a seed, and automatic hypotheses output; the longevity
  design has 525 possible four-vs-four pairings and will be run with 100;
- made smaller pooled selections reproducible through seeded SHA-256 candidate
  ranking and retained the complete 7-FG/6-BG pools as event denominators even
  if a selected hypothesis set omits a species;
- implemented an exact conditional event-separation p-value and labelled it
  nominal because the amino-acid signature is data-selected;
- preserved the legacy headerless bootstrap output and added the new analysis
  as an optional companion long-format table;
- added integration and unit tests for compatible merges, incompatible events,
  fixed-side validation, exact p-values, unique-species counting, and legacy
  compatibility; all 18 bundled CAAStools tests pass;
- created a local C4BPA demonstration under
  `../caastools.bootstrap.events.local-test/` using the 100 longevity
  hypotheses and the 3/4 observed-species thresholds;
- converted the production Nextflow pipeline to pooled discovery: one initial
  process now prepares the shared deterministic 100-hypothesis table and its
  provenance metadata, after which one SLURM process per alignment reuses that
  exact table;
- separated the per-gene legacy and event-level products into
  `caas-pooled/` and `caas-pooled-events/`, while publishing the shared table,
  seed, design size, and SHA-256 checksum under `metadata/`;
- validated the complete two-stage workflow locally on the real C4BPA
  alignment, including a successful resumed run in which both preparation and
  analysis were recovered from the Nextflow cache.

### 18 August 2026

- inspected the real legacy six-column CAAStools bootstrap format;
- implemented strict parsing and zero-cycle filtering;
- implemented pairwise local cluster filtering with 10-aa, 80%-density, and
  80%-cycle-recurrence defaults;
- added maximal interval reports and overlapping-interval handling;
- identified `THADA` as a gene-wide `b_74` dominance case;
- measured dominance and burden distributions across all available genes;
- selected a cycle-specific conditional 97.5th-percentile burden reference;
- implemented reference construction;
- implemented global cycle pruning with minimum 10-position and 80%-dominance
  requirements;
- preserved alternative cycles and recalculated modified CAAStools rows;
- validated the complete two-stage procedure on all 16,128 available results;
- expanded the automated downstream test suite to 17 passing tests;
- added minimum observed-species controls to the bundled CAAStools bootstrap
  implementation and configured the longevity workflow to require three
  observed FG and three observed BG species per cycle and position;
- added the positional hypergeometric p-value after the positive-cycle field
  in new CAAStools output and retained compatibility with legacy six-column
  results in the downstream parser;
- marked the present downstream validation as belonging to the earlier,
  less-restrictive bootstrap run.
