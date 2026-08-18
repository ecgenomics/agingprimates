#!/usr/bin/env python3

import csv
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from merge_and_filter_caas import EXPECTED_HEADER, merge_tables  # noqa: E402


def result_row(gene: str, pvalue: str) -> list[str]:
    return [
        gene,
        "test.caas.cfg",
        "10",
        "A/V",
        pvalue,
        "pattern1",
        "6",
        "5",
        "0",
        "0",
        "0",
        "0",
        "fg_species",
        "bg_species",
        "-",
    ]


class MergeAndFilterTests(unittest.TestCase):
    def test_one_header_and_inclusive_significance_filter(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            first = root / "first.caas.tsv"
            second = root / "second.caas.tsv"
            merged = root / "merged.tsv"
            significant = root / "significant.tsv"

            with first.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(EXPECTED_HEADER)
                writer.writerow(result_row("GENE1", "0.01"))
                writer.writerow(result_row("GENE2", "0.2"))

            with second.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(EXPECTED_HEADER)
                writer.writerow(EXPECTED_HEADER)
                writer.writerow(result_row("GENE3", "0.05"))

            counts = merge_tables(
                [first, second], merged, significant, alpha=0.05
            )
            self.assertEqual(counts, (3, 2))

            with merged.open("r", encoding="utf-8", newline="") as handle:
                merged_rows = list(csv.reader(handle, delimiter="\t"))
            with significant.open("r", encoding="utf-8", newline="") as handle:
                significant_rows = list(csv.reader(handle, delimiter="\t"))

            self.assertEqual(merged_rows.count(EXPECTED_HEADER), 1)
            self.assertEqual(len(merged_rows), 4)
            self.assertEqual(significant_rows.count(EXPECTED_HEADER), 1)
            self.assertEqual([row[0] for row in significant_rows[1:]], ["GENE1", "GENE3"])


if __name__ == "__main__":
    unittest.main()
