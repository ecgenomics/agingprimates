#!/usr/bin/env python3

"""Tests for construction of the genomic cycle-burden reference."""

from pathlib import Path
import sys
import tempfile
import unittest


SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

from build_cycle_burden_reference import (  # noqa: E402
    build_reference_rows,
    empirical_quantile,
    write_reference,
)
from filter_bootstrap_caas import load_cycle_burden_reference  # noqa: E402


def bootstrap_line(gene: str, position: int, cycles: tuple[str, ...] = ()) -> str:
    count = len(cycles)
    return "\t".join(
        [
            f"{gene}@{position}",
            str(count),
            "100",
            str(count / 100),
            ",".join(cycles),
            "template.cfg",
        ]
    )


class CycleBurdenReferenceTests(unittest.TestCase):
    def test_empirical_quantile_uses_linear_interpolation(self):
        self.assertAlmostEqual(0.3, empirical_quantile([0.1, 0.5], 0.5))
        self.assertAlmostEqual(0.5, empirical_quantile([0.1, 0.5], 1.0))

    def test_reference_is_cycle_specific_and_reloadable(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            gene_a = directory / "A.bootstrap.caas.tsv"
            gene_b = directory / "B.bootstrap.caas.tsv"

            gene_a.write_text(
                "\n".join(
                    bootstrap_line(
                        "A",
                        position,
                        tuple(
                            cycle
                            for cycle, last_position in (("b_1", 5), ("b_2", 2))
                            if position <= last_position
                        ),
                    )
                    for position in range(1, 11)
                )
                + "\n",
                encoding="utf-8",
            )
            gene_b.write_text(
                "\n".join(
                    bootstrap_line(
                        "B",
                        position,
                        tuple(
                            cycle
                            for cycle, last_position in (("b_1", 1), ("b_2", 4))
                            if position <= last_position
                        ),
                    )
                    for position in range(1, 11)
                )
                + "\n",
                encoding="utf-8",
            )

            rows = build_reference_rows([gene_a, gene_b], quantile=0.5)
            output_path = directory / "reference.tsv"
            write_reference(rows, output_path)
            loaded = load_cycle_burden_reference(output_path)

            self.assertEqual(["b_1", "b_2"], [row.cycle_id for row in rows])
            self.assertAlmostEqual(0.3, loaded["b_1"].burden_threshold)
            self.assertAlmostEqual(0.3, loaded["b_2"].burden_threshold)
            self.assertEqual(2, loaded["b_1"].genes_with_cycle)


if __name__ == "__main__":
    unittest.main()
