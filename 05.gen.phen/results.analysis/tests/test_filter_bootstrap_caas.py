#!/usr/bin/env python3

"""Tests for the CAAStools bootstrap result filter."""

from pathlib import Path
import sys
import tempfile
import unittest


SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

from filter_bootstrap_caas import (  # noqa: E402
    CycleBurdenThreshold,
    default_output_path,
    default_report_path,
    filter_records,
    filter_records_with_global_reference,
    load_cycle_burden_reference,
    parse_bootstrap_result,
    write_flagged_cycle_report,
    write_cluster_report,
    write_filtered_result,
)


def bootstrap_line(position: int, cycles: tuple[str, ...] = ()) -> str:
    count = len(cycles)
    return "\t".join(
        [
            f"TEST@{position}",
            str(count),
            "100",
            str(count / 100),
            ",".join(cycles),
            "template.cfg",
        ]
    )


class BootstrapFilterTests(unittest.TestCase):
    def parse_lines(self, lines: list[str]):
        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)
        input_path = Path(temporary_directory.name) / "input.tsv"
        input_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return parse_bootstrap_result(input_path)

    def test_zero_rows_and_exact_eighty_percent_recurrent_cluster_are_removed(self):
        lines = []
        for position in range(1, 11):
            cycle = "b_1" if position <= 8 else f"b_{position}"
            lines.append(bootstrap_line(position, (cycle,)))
        lines.append(bootstrap_line(11))

        records = self.parse_lines(lines)
        retained, intervals, summary = filter_records(records)

        self.assertEqual([], retained)
        self.assertEqual(1, len(intervals))
        self.assertEqual((1, 10), (intervals[0].start_position, intervals[0].end_position))
        self.assertEqual(("b_1",), intervals[0].recurrent_cycles)
        self.assertEqual(11, summary.total_positions)
        self.assertEqual(1, summary.zero_cycle_positions)
        self.assertEqual(10, summary.positive_positions)
        self.assertEqual(10, summary.cluster_positions_removed)
        self.assertEqual(0, summary.retained_positions)

    def test_nearly_consecutive_cluster_uses_coordinate_density(self):
        positions = [100, 101, 102, 103, 105, 106, 108, 109]
        records = self.parse_lines([bootstrap_line(position, ("b_7",)) for position in positions])

        retained, intervals, summary = filter_records(records)

        self.assertEqual([], retained)
        self.assertEqual(1, len(intervals))
        self.assertEqual(10, intervals[0].span)
        self.assertAlmostEqual(0.8, intervals[0].density)
        self.assertEqual(8, summary.cluster_positions_removed)

    def test_dense_cluster_with_different_cycles_is_retained(self):
        records = self.parse_lines(
            [bootstrap_line(position, (f"b_{position}",)) for position in range(20, 30)]
        )

        retained, intervals, summary = filter_records(records)

        self.assertEqual([], intervals)
        self.assertEqual(10, len(retained))
        self.assertEqual(0, summary.cluster_positions_removed)

    def test_recurrent_cluster_shorter_than_minimum_span_is_retained(self):
        records = self.parse_lines(
            [bootstrap_line(position, ("b_3",)) for position in range(40, 49)]
        )

        retained, intervals, _ = filter_records(records)

        self.assertEqual([], intervals)
        self.assertEqual(9, len(retained))

    def test_thresholds_are_configurable(self):
        positions = [1, 2, 3, 4, 5, 6, 8]
        records = self.parse_lines([bootstrap_line(position, ("b_4",)) for position in positions])

        retained_default, intervals_default, _ = filter_records(records)
        retained_relaxed, intervals_relaxed, _ = filter_records(
            records,
            min_span=8,
            min_density=0.875,
            min_cycle_recurrence=1.0,
        )

        self.assertEqual(7, len(retained_default))
        self.assertEqual([], intervals_default)
        self.assertEqual([], retained_relaxed)
        self.assertEqual(1, len(intervals_relaxed))

    def test_zero_based_caas_position_is_valid(self):
        records = self.parse_lines([bootstrap_line(0, ("b_1",))])

        self.assertEqual(0, records[0].position)

    def test_seven_column_pvalue_is_validated_and_preserved(self):
        line = "TEST@7\t1\t100\t0.01\tb_9\t0.0125\ttemplate.cfg"
        records = self.parse_lines([line])

        self.assertEqual(0.0125, records[0].hypergeometric_pvalue)

        with tempfile.TemporaryDirectory() as temporary_directory:
            output_path = Path(temporary_directory) / "filtered.tsv"
            write_filtered_result(records, output_path)
            self.assertEqual(line + "\n", output_path.read_text())

    def test_invalid_hypergeometric_pvalue_is_rejected(self):
        malformed = "TEST@7\t1\t100\t0.01\tb_9\t1.5\ttemplate.cfg"
        with self.assertRaisesRegex(ValueError, "p-value must be between 0 and 1"):
            self.parse_lines([malformed])

    def test_overlapping_intervals_remove_shared_positions_only_once(self):
        positions = [1, 2, 3, 4, 5, 6, 7, 10, 13, 14, 15, 16, 17, 18, 19]
        records = self.parse_lines([bootstrap_line(position, ("b_1",)) for position in positions])

        retained, intervals, summary = filter_records(records)

        self.assertEqual([], retained)
        self.assertEqual(2, len(intervals))
        self.assertEqual(15, summary.cluster_positions_removed)

    def test_default_output_names_insert_tag_before_tsv_suffix(self):
        input_path = Path("gene.bootstrap.caas.tsv")

        self.assertEqual(
            Path("gene.bootstrap.caas.filtered.tsv"),
            default_output_path(input_path),
        )
        self.assertEqual(
            Path("gene.bootstrap.caas.clusters.tsv"),
            default_report_path(input_path),
        )

    def test_global_filter_removes_only_flagged_cycle_and_recalculates_rows(self):
        lines = []
        for position in range(1, 21):
            if position <= 6:
                cycles = ("b_1",)
            elif position <= 8:
                cycles = ("b_1", "b_2")
            elif position <= 10:
                cycles = ("b_2",)
            else:
                cycles = ()
            lines.append(bootstrap_line(position, cycles))
        records = self.parse_lines(lines)
        thresholds = {
            "b_1": CycleBurdenThreshold("b_1", 0.975, 0.30, 100),
            "b_2": CycleBurdenThreshold("b_2", 0.975, 0.90, 100),
        }

        retained, intervals, flagged, summary = filter_records_with_global_reference(
            records,
            thresholds,
            min_cycle_positions=5,
        )

        self.assertEqual([], intervals)
        self.assertEqual(["b_1"], [item.cycle_id for item in flagged])
        self.assertEqual([7, 8, 9, 10], [record.position for record in retained])
        self.assertTrue(all(record.cycles == ("b_2",) for record in retained))
        self.assertTrue(all(record.positive_cycles == 1 for record in retained))
        self.assertTrue(all(record.frequency == 0.01 for record in retained))
        self.assertEqual(6, summary.positions_emptied_by_global_filter)
        self.assertEqual(4, summary.positive_positions_after_global_filter)
        self.assertEqual(4, summary.retained_positions)

    def test_global_filter_requires_reference_for_every_observed_cycle(self):
        records = self.parse_lines([bootstrap_line(1, ("b_1",))])
        with self.assertRaisesRegex(ValueError, "lacks threshold.*b_1"):
            filter_records_with_global_reference(
                records,
                {},
                min_cycle_positions=1,
            )

    def test_reference_loader_and_flagged_cycle_report(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            reference_path = Path(temporary_directory) / "reference.tsv"
            reference_path.write_text(
                "cycle_id\tquantile\tburden_threshold\tgenes_with_cycle\n"
                "b_1\t0.975\t0.3\t100\n",
                encoding="utf-8",
            )
            thresholds = load_cycle_burden_reference(reference_path)
            records = self.parse_lines(
                [bootstrap_line(position, ("b_1",)) for position in range(1, 11)]
            )
            _, _, flagged, _ = filter_records_with_global_reference(
                records,
                thresholds,
            )
            report_path = Path(temporary_directory) / "flagged.tsv"
            write_flagged_cycle_report(flagged, report_path)

            self.assertEqual(0.3, thresholds["b_1"].burden_threshold)
            report_lines = report_path.read_text().splitlines()
            self.assertEqual(2, len(report_lines))
            self.assertIn("cycle_burden", report_lines[0])
            self.assertIn("\tb_1\t", report_lines[1])

    def test_writers_preserve_caas_rows_and_create_cluster_report(self):
        records = self.parse_lines(
            [bootstrap_line(position, ("b_1",)) for position in range(1, 11)]
            + [bootstrap_line(20, ("b_20",))]
        )
        retained, intervals, _ = filter_records(records)

        with tempfile.TemporaryDirectory() as temporary_directory:
            output_path = Path(temporary_directory) / "filtered.tsv"
            report_path = Path(temporary_directory) / "clusters.tsv"
            write_filtered_result(retained, output_path)
            write_cluster_report(intervals, report_path)

            self.assertEqual(bootstrap_line(20, ("b_20",)) + "\n", output_path.read_text())
            report_lines = report_path.read_text().splitlines()
            self.assertEqual(2, len(report_lines))
            self.assertIn("recurrent_cycles", report_lines[0])
            self.assertIn("b_1", report_lines[1])

    def test_mismatched_cycle_count_is_rejected(self):
        malformed = "TEST@1\t2\t100\t0.02\tb_1\ttemplate.cfg"
        with self.assertRaisesRegex(ValueError, "2.*1 cycle IDs"):
            self.parse_lines([malformed])


if __name__ == "__main__":
    unittest.main()
