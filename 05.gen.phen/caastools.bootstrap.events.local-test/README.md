# Local CAAStools bootstrap event demonstration

This directory records a local, reproducible demonstration of the optional
species- and amino-acid-aware bootstrap summary added to the CAAStools copy in
`../lq.table2.bootstrap.nextflow/bin/caastools/`.

Nothing in this directory is used by the SLURM workflow. It is deliberately a
small development fixture for inspecting the new output before changing or
rerunning the cluster analysis.

## Biological motivation

The 100 longevity resampling rows overlap in species composition. Their
identifiers (`b_1` through `b_100`) are useful discovery hypotheses, but their
raw count should not be interpreted as 100 independent pieces of evidence.
The event summary instead asks:

1. which amino-acid contrast was detected;
2. which unique FG and BG species support that contrast across all compatible
   positive hypotheses;
3. how many observed and total discovery species those supporters represent;
4. whether other discovery species carry amino acids that conflict with the
   contrast;
5. whether distinct, incompatible contrasts must remain separate at the same
   alignment position.

## Implementation tested here

The local CAAStools version now performs the following optional operations:

1. infers the complete discovery pools from all resampling hypotheses;
2. verifies that no species changes from FG to BG across hypotheses;
3. retains each positive hypothesis's observed species and amino acids;
4. represents each hypothesis as a disjoint FG/BG amino-acid signature;
5. collapses identical signatures;
6. finds all maximal groups of mutually compatible signatures;
7. creates one event per maximal group and keeps incompatible groups separate;
8. assigns the CAAS pattern to the final merged signature: pattern 1 for
   one-vs-one, 2 for one-vs-many, 3 for many-vs-one, and 4 for many-vs-many;
9. deduplicates species within each event;
10. reports the dominant FG and BG amino acid and its proportion among unique
    supporting species, preserving all co-dominant amino acids in a tie;
11. counts supporting and all-observed amino acids using unique species;
12. reports support over observed and complete discovery denominators;
13. reports conflicting species and their amino acids;
14. ranks the primary event by balanced discovery support, total unique support,
    nominal event p-value, and deterministic signature order;
15. calculates an exact conditional, oriented signature-separation p-value;
16. preserves the historical bootstrap table unchanged as a separate output.

The exact event test conditions on the number of observed FG and BG species.
Its statistic is the smaller of the FG and BG signature-match fractions, so a
large value requires support from both sides. The calculation uses counts of
FG-signature, BG-signature, and other amino acids and is therefore exact but
much faster than enumerating labeled species one by one. The value is called
`event_nominal_pvalue` because the event signature is selected from the data;
dataset-wide multiple-testing correction is not performed here.

## Command used

From this directory:

```bash
python3 ../lq.table2.bootstrap.nextflow/bin/caastools/ct bootstrap \
  -a ../lq.table2.nextflow/inputs/alignments/C4BPA.Homo_sapiens.filter2.phy \
  -t ../lq.table2.bootstrap.nextflow/inputs/longevity.template.001.caas.cfg \
  -s ../lq.table2.bootstrap.nextflow/inputs/longevity.100-comparisons.resampling.tsv \
  -o real-data/C4BPA.bootstrap.caas.tsv \
  --event-output real-data/C4BPA.bootstrap.caas.events.tsv \
  --fmt phylip-relaxed \
  --filter_significant no \
  --patterns 1,2,3 \
  --max_bg_gaps NO --max_fg_gaps NO --max_gaps NO \
  --max_gaps_per_position 0.5 \
  --max_bg_miss NO --max_fg_miss NO --max_miss NO \
  --min_fg_observed 3 --min_bg_observed 3
```

Providing `--event-output` enables event summarization. Alternatively,
`--summarize-species yes` derives the event filename automatically from the
ordinary output filename.

The significance prefilter was disabled in this local demonstration so that
the event representation could be inspected independently. This does not
change the current cluster workflow configuration.

## Integrated pooled-discovery command

CAAStools can now generate the hypothesis table itself from the complete
seven-FG/six-BG pool file. The production input is
`../lq.table2.bootstrap.nextflow/inputs/longevity.full-pools.caas.cfg`; a copy
is retained in `synthetic/longevity.full-pools.cfg` to keep this demonstration
self-contained.

The 100-comparison demonstration was generated with:

```bash
python3 ../lq.table2.bootstrap.nextflow/bin/caastools/ct pooled-discovery \
  -a synthetic/compatible.phy \
  -t synthetic/longevity.full-pools.cfg \
  -o synthetic/compatible.pooled100.caas.tsv \
  --event-output synthetic/compatible.pooled100.caas.events.tsv \
  --hypotheses-output synthetic/compatible.pooled100.hypotheses.tsv \
  --fmt phylip-relaxed \
  --fg-size 4 --bg-size 4 \
  --comparisons 100 --seed 260811 \
  --patterns 1,2,3 \
  --filter_significant no \
  --min-fg-observed 3 --min-bg-observed 3
```

The complete design contains 525 unique four-vs-four pairings. Omitting
`--comparisons` uses the new default `max` and generates all 525, as recorded
in `synthetic/compatible.pooledmax.hypotheses.tsv`. A smaller selection is
ranked deterministically from the seed with SHA-256. The 100 hypotheses created
here are therefore reproducible, unique, and explicitly saved, but are not
expected to be identical to the earlier 100 selected with R's `sample()`.

With the synthetic column, the selected 100 hypotheses contain ten positive
pattern-1/2/3 subsets; all merge into one final pattern-4 event supported by
the complete 7/7 FG and 6/6 BG discovery pools. With all 525 pairings, 49
hypotheses are positive, but the final event, its 13 supporting species, amino
acid signature, dominant amino-acid proportions, and nominal event p-value are
unchanged. This illustrates why the event is the biological unit and the raw
hypothesis count is not.

A real C4BPA run with the new deterministic 100-hypothesis selection produced
334 legacy position rows and 19 event rows at 19 positions. Its files are
`real-data/C4BPA.pooled100.caas.tsv`,
`real-data/C4BPA.pooled100.caas.events.tsv`, and the shared realized input
`real-data/longevity.pooled100.hypotheses.tsv`.

## Demonstration outputs

### Synthetic compatible-hypothesis example

The `synthetic/` fixture makes the merge directly visible. Its historical
output is:

```text
compatible@0  3  4  0.75  b_1,b_2,b_3  0.002797202797202797  synthetic/trait.cfg
```

The three positive hypotheses have signatures `A/V`, `A/{I,V}`, and
`{A,G}/V`. The fourth hypothesis is not positive under patterns 1–3. The event
table correctly combines the three compatible positive signatures into one
event:

- merged signature: FG `{A,G}` versus BG `{I,V}`;
- dominant amino acids: FG `A=0.8`; BG `V=0.8`;
- final event pattern: `4` (many-vs-many; here two amino acids versus two);
- supporting hypotheses: `b_1,b_2,b_3`;
- unique support: 5/7 FG discovery species and 5/6 BG discovery species;
- amino-acid support: FG `A=4,G=1`; BG `I=1,V=4`;
- balanced discovery support: `5/7 = 0.714285714286`;
- no conflicting species;
- nominal event p-value: `0.000582750582751`.

The output uses the complete longevity species names—for example,
`Aotus_trivirgatus=A`, `Eulemur_macaco=G`, and
`Mico_humeralifer=I`—rather than anonymous labels such as `f1` or `b1`.
Thus, the result is no longer interpreted simply as “three positive cycles.”
It becomes one pattern-4 amino-acid event supported by ten unique species,
with explicit FG/BG denominators and a trace back to all three discovery
hypotheses.

### Real C4BPA example

`real-data/C4BPA.bootstrap.caas.tsv` is the backwards-compatible bootstrap
table. It contains 334 evaluated positions, including zero-hypothesis rows.

`real-data/C4BPA.bootstrap.caas.events.tsv` is the new headered companion
table. It contains 32 positive events at 32 positions. In this particular gene
there are no multi-event positions, the largest event contains six compatible
hypotheses, and the largest event contains ten unique supporting species.

All 32 events contain at least one conflicting discovery species. This is a
useful illustration of why the new representation is more conservative than a
raw hypothesis count: a four-versus-four subset can be CAAS-positive even when
other species in the fixed seven-FG/six-BG discovery pools carry the opposite
event amino acid.

For example, the first event is at C4BPA position 1:

- positive hypothesis: `b_74`;
- event signature and final pattern: FG `{A,L}` versus BG `{S}`, pattern 3;
- dominant amino acids: FG `L=0.75`; BG `S=1`;
- unique support: 4 FG and 4 BG species;
- denominators: 7 observed/discovery FG and 6 observed/discovery BG species;
- balanced discovery support: `min(4/7, 4/6) = 4/7`;
- conflicts: three FG species carrying `S` and one BG species carrying `L`;
- nominal event p-value: `0.179487179487`;
- positional hypergeometric p-value: `0.002514907815566993`.

The difference between the two p-values is intentional. The positional value
asks whether a nominal four-versus-four CAAS is rare given the whole alignment
column. The event value evaluates the selected oriented signature across the
complete observed discovery pools and therefore exposes contrary species that
were absent from the positive subset.

## Automated validation

Run from `../lq.table2.bootstrap.nextflow/bin/caastools/`:

```bash
python3 -m unittest discover -s test -p 'test_*.py' -v
```

The 18-test suite covers the pre-existing significance and coverage behavior,
plus compatible signature merging, incompatible event separation,
non-transitive compatibility with all maximal alternatives, unique-species
counting, all four final event-pattern classes, explicit species names,
dominant amino-acid proportions and ties, amino-acid counts, exact event
p-values, fixed-side validation, deterministic hypothesis ordering, all-525
generation, fixed-count uniqueness, seeded reproducibility, invalid pooling
requests, automatic hypothesis-file output, and preservation of the legacy
bootstrap output.

## Current scope

This implementation does not yet replace downstream alignment-train filtering,
proteome-wide FDR, gene-level burden summaries, or SLURM timeout handling. It
also does not modify the Nextflow workflow in this development step. Those
decisions should be made after inspecting event summaries across more genes.
