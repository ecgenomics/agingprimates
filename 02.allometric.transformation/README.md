# 02. Allometric transformation

This directory contains the reproducible construction of the primary
mammal-wide longevity quotient described in
`../theory/02.allometric.transformation.logic.md`.

Run from any directory with:

```bash
python3 build_longevity_quotient.py
```

The calculation is fixed (not fitted to these primates):

```text
expected_MLS_years = 6.47 * body_mass_g^0.189
LQ_mammal = MLS_years / expected_MLS_years
log_LQ_mammal = log(MLS_years) - log(expected_MLS_years)
```

Outputs:

- `longevity_quotient.tsv`: analysis-ready dataset with all valid phenotypic
  observations. Genome availability is annotated but is not an inclusion rule.
- `longevity_quotient_audit.tsv`: every source row, including excluded records
  and reasons.
- `build_report.json`: counts, input checksums, and historical comparison.
- `MD5SUMS` and `longevity_quotient.tsv.md5`: output integrity checksums.

The dataset does not apply long/short-lived thresholds and does not calculate
family summaries. Those are downstream analytical decisions. Taxonomic matches
are exact string matches only; unresolved names remain explicitly flagged.
