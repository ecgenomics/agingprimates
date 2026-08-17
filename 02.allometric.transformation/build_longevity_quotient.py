#!/usr/bin/env python3
"""Build the fixed mammal-wide longevity quotient dataset reproducibly."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import pandas as pd


INTERCEPT = 6.47
EXPONENT = 0.189
MLS_TRAIT = "PTD00065"
MASS_TRAIT = "PTD00014"


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_two_column_no_header(path: Path, value_name: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", header=None, names=["Species", value_name])
    frame[value_name] = pd.to_numeric(frame[value_name], errors="coerce")
    return frame.dropna(subset=["Species"]).drop_duplicates("Species", keep="last")


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=here.parent / "inputs")
    parser.add_argument("--historical", type=Path, default=here.parent / "dataset" / "lqdf.tab")
    parser.add_argument("--output-dir", type=Path, default=here)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    inputs = args.inputs.resolve()
    outdir = args.output_dir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    phenomic_path = inputs / "nhp.phenomic.dataset.csv"
    sources_path = inputs / "trait.data.sources.csv"
    tree_path = inputs / "science.abn7829_data_s4.nex.tree"
    genomic_path = inputs / "tmb.sample.count.tab"
    required = [phenomic_path, sources_path, tree_path, genomic_path]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing required input(s): " + ", ".join(missing))

    raw = pd.read_csv(phenomic_path, low_memory=False)
    needed = {"GroupName", MLS_TRAIT, MASS_TRAIT}
    absent = sorted(needed.difference(raw.columns))
    if absent:
        raise ValueError(f"Phenomic input lacks columns: {absent}")
    if raw["GroupName"].duplicated().any():
        raise ValueError("GroupName must be unique; duplicate records were found")

    split = raw["GroupName"].astype("string").str.split("_", n=3, expand=True)
    if split.shape[1] != 4 or split.isna().any(axis=None):
        raise ValueError("Every GroupName must contain superfamily, family, genus, and species")

    audit = pd.DataFrame({
        "original_group_name": raw["GroupName"],
        "superfamily": split[0],
        "family": split[1],
        "original_species_name": split[2] + "_" + split[3],
        "accepted_tree_tip": split[2] + "_" + split[3],
        "original_MLS_value": raw[MLS_TRAIT],
        "original_MLS_unit": "years",
        "MLS_years": pd.to_numeric(raw[MLS_TRAIT], errors="coerce"),
        "original_body_mass_value": raw[MASS_TRAIT],
        "original_body_mass_unit": "kg",
        "body_mass_g": pd.to_numeric(raw[MASS_TRAIT], errors="coerce") * 1000.0,
    })

    # Exact matching is deliberately explicit; name reconciliation belongs in a mapping table.
    tree_text = tree_path.read_text(encoding="utf-8")
    audit["taxonomic_mapping_status"] = audit["accepted_tree_tip"].map(
        lambda tip: "exact_tree_tip_match" if str(tip) in tree_text else "not_in_tree_exactly"
    )
    sources = pd.read_csv(sources_path).set_index("ID")
    for trait in (MLS_TRAIT, MASS_TRAIT):
        if trait not in sources.index:
            raise ValueError(f"No provenance record for {trait}")
    audit["MLS_trait_id"] = MLS_TRAIT
    audit["MLS_source"] = sources.at[MLS_TRAIT, "Source"]
    audit["MLS_source_link"] = sources.at[MLS_TRAIT, "Link"]
    audit["mass_trait_id"] = MASS_TRAIT
    audit["mass_source"] = sources.at[MASS_TRAIT, "Source"]
    audit["mass_source_link"] = sources.at[MASS_TRAIT, "Link"]
    audit["source_access_date"] = "not recorded in source input"
    audit["captive_wild_context"] = "not recorded in source input"
    audit["sex"] = "not recorded in source input"
    audit["sample_size"] = "not recorded in source input"
    audit["duplicate_resolution_rule"] = "not applicable: GroupName is unique"

    mls_ok = audit["MLS_years"].map(lambda x: pd.notna(x) and math.isfinite(x) and x > 0)
    mass_ok = audit["body_mass_g"].map(lambda x: pd.notna(x) and math.isfinite(x) and x > 0)
    audit["included"] = mls_ok & mass_ok
    audit["exclusion_reason"] = ""
    audit.loc[~mls_ok & mass_ok, "exclusion_reason"] = "missing_or_nonpositive_MLS"
    audit.loc[mls_ok & ~mass_ok, "exclusion_reason"] = "missing_or_nonpositive_body_mass"
    audit.loc[~mls_ok & ~mass_ok, "exclusion_reason"] = "missing_or_nonpositive_MLS_and_body_mass"

    valid = audit.loc[audit["included"]].copy()
    valid["log_expected_MLS"] = math.log(INTERCEPT) + EXPONENT * valid["body_mass_g"].map(math.log)
    valid["expected_MLS_years"] = valid["log_expected_MLS"].map(math.exp)
    valid["LQ_mammal"] = valid["MLS_years"] / valid["expected_MLS_years"]
    valid["log_LQ_mammal"] = valid["MLS_years"].map(math.log) - valid["log_expected_MLS"]
    if not ((valid["LQ_mammal"].map(math.log) - valid["log_LQ_mammal"]).abs() < 1e-12).all():
        raise AssertionError("LQ and log-LQ calculations are inconsistent")

    genomic = read_two_column_no_header(genomic_path, "genome_samples")
    valid = valid.merge(genomic, how="left", left_on="accepted_tree_tip", right_on="Species")
    valid = valid.drop(columns="Species")
    valid["genome_samples"] = valid["genome_samples"].fillna(0).astype(int)
    valid["genome_available"] = valid["genome_samples"] > 0
    valid = valid.sort_values("accepted_tree_tip", kind="stable").reset_index(drop=True)

    final_columns = [
        "accepted_tree_tip", "original_species_name", "superfamily", "family",
        "MLS_years", "body_mass_g", "expected_MLS_years", "LQ_mammal",
        "log_LQ_mammal", "MLS_trait_id", "MLS_source", "MLS_source_link",
        "mass_trait_id", "mass_source", "mass_source_link",
        "taxonomic_mapping_status", "genome_samples", "genome_available",
    ]
    dataset_path = outdir / "longevity_quotient.tsv"
    audit_path = outdir / "longevity_quotient_audit.tsv"
    report_path = outdir / "build_report.json"
    valid[final_columns].to_csv(dataset_path, sep="\t", index=False, float_format="%.12g")
    audit.sort_values("original_group_name").to_csv(audit_path, sep="\t", index=False, float_format="%.12g")

    comparison = {"historical_file_present": args.historical.is_file()}
    if args.historical.is_file():
        old = pd.read_csv(args.historical, sep="\t")
        overlap = valid[["accepted_tree_tip", "LQ_mammal"]].merge(
            old[["Species", "LQ"]], left_on="accepted_tree_tip", right_on="Species"
        )
        delta = (overlap["LQ_mammal"] - overlap["LQ"]).abs()
        comparison.update({
            "overlapping_species": int(len(overlap)),
            "maximum_absolute_LQ_difference": float(delta.max()) if len(delta) else None,
            "differences_above_1e-12": int((delta > 1e-12).sum()),
        })

    report = {
        "formula": "expected_MLS_years = 6.47 * body_mass_g^0.189",
        "input_md5": {path.name: md5(path) for path in required},
        "records_total": int(len(audit)),
        "records_included": int(audit["included"].sum()),
        "records_excluded": int((~audit["included"]).sum()),
        "exclusion_counts": audit.loc[~audit["included"], "exclusion_reason"].value_counts().to_dict(),
        "exact_tree_tip_matches_in_dataset": int((valid["taxonomic_mapping_status"] == "exact_tree_tip_match").sum()),
        "genome_available_in_dataset": int(valid["genome_available"].sum()),
        "historical_comparison": comparison,
    }
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checksum_paths = [dataset_path, audit_path, report_path]
    checksum_path = outdir / "MD5SUMS"
    checksum_path.write_text("".join(f"{md5(path)}  {path.name}\n" for path in checksum_paths), encoding="utf-8")
    (outdir / "longevity_quotient.tsv.md5").write_text(
        f"{md5(dataset_path)}  {dataset_path.name}\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
