#!/usr/bin/env python3

import csv
import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = ROOT.parents[1] / "03.phenotype.shift.detection.lq"
PSS_RESULTS = DATA_ROOT / "results/LQ.score_results.tsv"
PHENOTYPE_TABLE = DATA_ROOT / "input/longevity_quotient.tsv"
PSS_CONFIGS = {
    "max_delta_lq_within_top_1pct_pss": ROOT / "configs/01.max-delta-lq.caas.cfg",
    "max_pss_final_score": ROOT / "configs/02.max-pss-score.caas.cfg",
}
ACTIVE_4V4_CONFIGS = {
    "max_delta_lq_within_top_1pct_pss": ROOT
    / "configs/01.max-delta-lq.4v4.caas.cfg",
    "max_pss_final_score": ROOT / "configs/02.max-pss-score.4v4.caas.cfg",
}
ACTIVE_4V4_RANKS = {
    "max_delta_lq_within_top_1pct_pss": [1, 2, 3, 5],
    "max_pss_final_score": [1, 2, 3, 4],
}
EXTREMES_CONFIG = ROOT / "configs/03.absolute-lq-extremes.caas.cfg"


def read_config(path: Path) -> dict[str, int]:
    membership: dict[str, int] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            species, state = line.rstrip("\n").split("\t")
            if species in membership:
                raise AssertionError(f"Duplicate species in {path}:{line_number}: {species}")
            membership[species] = int(state)
    return membership


class HypothesisConfigTests(unittest.TestCase):
    def read_pair_manifest(self) -> list[dict[str, str]]:
        with (ROOT / "inputs/pss_pair_manifest.tsv").open(
            "r", encoding="utf-8", newline=""
        ) as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def read_active_4v4_manifest(self) -> list[dict[str, str]]:
        with (ROOT / "inputs/pss_pair_manifest.4v4.tsv").open(
            "r", encoding="utf-8", newline=""
        ) as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_active_4v4_configs_match_manifest_and_group_sizes(self) -> None:
        expected = {criterion: {} for criterion in ACTIVE_4V4_CONFIGS}
        observed_ranks = {criterion: [] for criterion in ACTIVE_4V4_CONFIGS}
        for row in self.read_active_4v4_manifest():
            criterion = row["criterion"]
            expected[criterion][row["high_lq_species"]] = 1
            expected[criterion][row["low_lq_species"]] = 0
            observed_ranks[criterion].append(int(row["criterion_rank"]))

        for criterion, config_path in ACTIVE_4V4_CONFIGS.items():
            with self.subTest(criterion=criterion):
                observed = read_config(config_path)
                self.assertEqual(observed, expected[criterion])
                self.assertEqual(sum(state == 1 for state in observed.values()), 4)
                self.assertEqual(sum(state == 0 for state in observed.values()), 4)
                self.assertEqual(observed_ranks[criterion], ACTIVE_4V4_RANKS[criterion])

    def test_active_glob_selects_only_two_4v4_configs(self) -> None:
        expected_paths = set(ACTIVE_4V4_CONFIGS.values())
        self.assertEqual(set((ROOT / "configs").glob("*.4v4.caas.cfg")), expected_paths)

        expected_setting = 'caas_config_glob = "${projectDir}/configs/*.4v4.caas.cfg"'
        for config_path in (ROOT / "nextflow.config", ROOT / "conf/cluster.config"):
            with self.subTest(config=config_path.name):
                self.assertIn(
                    expected_setting,
                    config_path.read_text(encoding="utf-8"),
                )

    def test_active_4v4_manifest_is_subset_of_validated_pair_manifest(self) -> None:
        historical_rows = {
            (row["criterion"], row["criterion_rank"]): row
            for row in self.read_pair_manifest()
        }
        active_rows = self.read_active_4v4_manifest()
        self.assertEqual(len(active_rows), 8)
        for row in active_rows:
            key = (row["criterion"], row["criterion_rank"])
            self.assertIn(key, historical_rows)
            self.assertEqual(row, historical_rows[key])

    def test_pss_configs_match_pair_manifest(self) -> None:
        expected = {criterion: {} for criterion in PSS_CONFIGS}
        for row in self.read_pair_manifest():
            expected[row["criterion"]][row["high_lq_species"]] = 1
            expected[row["criterion"]][row["low_lq_species"]] = 0

        for criterion, config_path in PSS_CONFIGS.items():
            with self.subTest(criterion=criterion):
                observed = read_config(config_path)
                self.assertEqual(observed, expected[criterion])
                self.assertEqual(set(observed.values()), {0, 1})

    def test_pair_manifest_reproduces_both_pss_selection_criteria(self) -> None:
        with PSS_RESULTS.open("r", encoding="utf-8", newline="") as handle:
            pss_rows = list(csv.DictReader(handle, delimiter="\t"))

        by_score = sorted(
            pss_rows, key=lambda row: float(row["FinalScore"]), reverse=True
        )
        top_one_percent = by_score[: math.ceil(len(by_score) * 0.01)]
        expected = {
            "max_delta_lq_within_top_1pct_pss": sorted(
                top_one_percent,
                key=lambda row: abs(
                    float(row["TraitValue1"]) - float(row["TraitValue2"])
                ),
                reverse=True,
            )[:5],
            "max_pss_final_score": by_score[:5],
        }

        manifest_by_criterion = {criterion: [] for criterion in PSS_CONFIGS}
        for row in self.read_pair_manifest():
            manifest_by_criterion[row["criterion"]].append(row)

        for criterion, selected_rows in expected.items():
            observed_rows = sorted(
                manifest_by_criterion[criterion],
                key=lambda row: int(row["criterion_rank"]),
            )
            self.assertEqual(len(observed_rows), 5)
            for observed, selected in zip(observed_rows, selected_rows):
                trait_1 = float(selected["TraitValue1"])
                trait_2 = float(selected["TraitValue2"])
                if trait_1 > trait_2:
                    high_species, high_lq = selected["Species1"], trait_1
                    low_species, low_lq = selected["Species2"], trait_2
                else:
                    high_species, high_lq = selected["Species2"], trait_2
                    low_species, low_lq = selected["Species1"], trait_1

                self.assertEqual(observed["high_lq_species"], high_species)
                self.assertEqual(observed["low_lq_species"], low_species)
                self.assertAlmostEqual(float(observed["high_lq"]), high_lq)
                self.assertAlmostEqual(float(observed["low_lq"]), low_lq)
                self.assertAlmostEqual(
                    float(observed["absolute_lq_difference"]), high_lq - low_lq
                )
                self.assertAlmostEqual(
                    float(observed["pss_final_score"]), float(selected["FinalScore"])
                )

    def test_absolute_lq_extremes_config_and_manifest(self) -> None:
        with PHENOTYPE_TABLE.open("r", encoding="utf-8", newline="") as handle:
            phenotype_rows = list(csv.DictReader(handle, delimiter="\t"))

        by_lq = sorted(
            phenotype_rows, key=lambda row: float(row["LQ_mammal"])
        )
        expected_bg = by_lq[:5]
        descending_lq = list(reversed(by_lq))
        expected_fg = []
        seen_genera = set()
        for row in descending_lq:
            genus = row["accepted_tree_tip"].split("_", maxsplit=1)[0]
            if genus in seen_genera:
                continue
            seen_genera.add(genus)
            expected_fg.append(row)
            if len(expected_fg) == 5:
                break
        expected_membership = {
            **{row["accepted_tree_tip"]: 1 for row in expected_fg},
            **{row["accepted_tree_tip"]: 0 for row in expected_bg},
        }
        self.assertEqual(read_config(EXTREMES_CONFIG), expected_membership)

        with (ROOT / "inputs/lq_extremes_manifest.tsv").open(
            "r", encoding="utf-8", newline=""
        ) as handle:
            manifest_rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(len(manifest_rows), 10)
        for group, expected_rows, source_order in (
            ("FG_high_lq", expected_fg, descending_lq),
            ("BG_low_lq", expected_bg, by_lq),
        ):
            observed_rows = [row for row in manifest_rows if row["group"] == group]
            observed_rows.sort(key=lambda row: int(row["group_rank"]))
            for observed, expected in zip(observed_rows, expected_rows):
                self.assertEqual(observed["species"], expected["accepted_tree_tip"])
                self.assertEqual(
                    int(observed["source_extreme_rank"]),
                    source_order.index(expected) + 1,
                )
                self.assertAlmostEqual(
                    float(observed["lq_mammal"]), float(expected["LQ_mammal"])
                )
                self.assertEqual(observed["genome_samples"], expected["genome_samples"])
                self.assertEqual(
                    observed["genome_available"], expected["genome_available"]
                )


if __name__ == "__main__":
    unittest.main()
