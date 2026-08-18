#!/usr/bin/env python3

"""Filter globally recurrent cycles and suspicious clusters from CAAStools output.

CAAStools bootstrap writes one headerless, tab-separated row per tested alignment
position.  The expected columns are::

    gene@position  positive_cycles  total_cycles  frequency  cycle_ids  [pvalue]  template

When a genomic cycle-burden reference is supplied, this script performs three
filters in the following order:

1. cycle IDs that are anomalously recurrent across the whole gene are removed;
2. rows left with no positive bootstrap cycle are removed;
3. CAAS-positive positions are removed when they belong to a dense interval in
   which the same bootstrap cycle recurs unusually often.

A global cycle is flagged only when it has at least ``--min-cycle-positions``
occurrences, covers at least ``--min-global-cycle-dominance`` of the gene's
positive positions, and its burden exceeds the cycle-specific genomic cutoff.
Burden is the fraction of all input positions positive for that cycle.  The
reference is created by ``build_cycle_burden_reference.py``.

An interval is suspicious when all of the following are true:

* its inclusive coordinate span is at least ``--min-span`` amino acids;
* CAAS-positive positions occupy at least ``--min-density`` of that span;
* at least one cycle occurs in at least ``--min-cycle-recurrence`` of the
  CAAS-positive positions in the interval.

The interval endpoints are drawn pairwise from the positive CAAS positions.
Overlapping qualifying intervals are allowed; every positive position covered
by at least one qualifying interval is removed.  The output preserves the
legacy six-column or current seven-column CAAStools format with counts and
frequencies recalculated. Separate TSV reports record globally removed cycles
and locally filtered intervals.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from dataclasses import dataclass, replace
from pathlib import Path
import re
import sys
from typing import Iterable, Sequence


CYCLE_PATTERN = re.compile(r"b_[0-9]+$")


@dataclass(frozen=True)
class BootstrapRecord:
    """One validated row from a CAAStools bootstrap result."""

    fields: tuple[str, ...]
    gene: str
    position: int
    positive_cycles: int
    total_cycles: int
    frequency: float
    cycles: tuple[str, ...]
    hypergeometric_pvalue: float | None
    line_number: int

    @property
    def is_positive(self) -> bool:
        return self.positive_cycles > 0


@dataclass(frozen=True)
class CycleBurdenThreshold:
    """One cycle-specific cutoff loaded from the genomic reference."""

    cycle_id: str
    quantile: float
    burden_threshold: float
    genes_with_cycle: int


@dataclass(frozen=True)
class FlaggedCycle:
    """A gene-cycle association selected by the global burden filter."""

    gene: str
    cycle_id: str
    occurrences: int
    positive_positions: int
    total_positions: int
    dominance: float
    burden: float
    reference_quantile: float
    burden_threshold: float
    positions: tuple[int, ...]


@dataclass(frozen=True)
class SuspiciousInterval:
    """A maximal-by-containment interval that satisfies both cluster tests."""

    gene: str
    start_position: int
    end_position: int
    caas_positions: tuple[int, ...]
    recurrent_cycles: tuple[str, ...]
    cycle_counts: tuple[tuple[str, int], ...]

    @property
    def span(self) -> int:
        return self.end_position - self.start_position + 1

    @property
    def density(self) -> float:
        return len(self.caas_positions) / self.span

    @property
    def maximum_cycle_count(self) -> int:
        return max(count for _, count in self.cycle_counts)

    @property
    def maximum_cycle_fraction(self) -> float:
        return self.maximum_cycle_count / len(self.caas_positions)


@dataclass(frozen=True)
class FilterSummary:
    """Counts produced while filtering one bootstrap result."""

    total_positions: int
    zero_cycle_positions: int
    positive_positions: int
    suspicious_intervals: int
    cluster_positions_removed: int
    retained_positions: int


@dataclass(frozen=True)
class CombinedFilterSummary:
    """Counts spanning global cycle pruning and local cluster filtering."""

    total_positions: int
    input_zero_cycle_positions: int
    input_positive_positions: int
    globally_flagged_cycles: int
    cycle_associations_removed: int
    positions_emptied_by_global_filter: int
    positive_positions_after_global_filter: int
    suspicious_intervals: int
    cluster_positions_removed: int
    retained_positions: int


def _validate_fraction(value: float, option_name: str) -> float:
    if not 0.0 <= value <= 1.0:
        raise ValueError(
            f"{option_name} must be between 0 and 1 (received {value})"
        )
    return value


def fraction(value: str) -> float:
    """Argparse converter for probabilities and proportions."""

    try:
        parsed = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(f"not a number: {value}") from error
    try:
        return _validate_fraction(parsed, "fraction")
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error)) from error


def parse_bootstrap_result(input_path: Path) -> list[BootstrapRecord]:
    """Read and strictly validate one headerless CAAStools bootstrap result."""

    records: list[BootstrapRecord] = []
    observed_positions: set[int] = set()
    observed_gene: str | None = None

    with input_path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue

            fields = line.split("\t")
            if len(fields) not in (6, 7):
                raise ValueError(
                    f"{input_path}:{line_number}: expected 6 or 7 tab-separated columns, "
                    f"found {len(fields)}"
                )

            position_id, positive_text, total_text, frequency_text, cycles_text = fields[:5]
            pvalue_text = fields[5] if len(fields) == 7 else None
            gene, separator, position_text = position_id.rpartition("@")
            if not separator or not gene:
                raise ValueError(
                    f"{input_path}:{line_number}: invalid position identifier {position_id!r}; "
                    "expected gene@position"
                )

            try:
                position = int(position_text)
                positive_cycles = int(positive_text)
                total_cycles = int(total_text)
                empirical_frequency = float(frequency_text)
            except ValueError as error:
                raise ValueError(
                    f"{input_path}:{line_number}: position, counts, or frequency is not numeric"
                ) from error

            # CAAStools reports zero-based alignment coordinates, so position 0
            # is valid even though most result files begin at a later residue.
            if position < 0:
                raise ValueError(f"{input_path}:{line_number}: position must be zero or greater")
            if positive_cycles < 0 or total_cycles < 1 or positive_cycles > total_cycles:
                raise ValueError(
                    f"{input_path}:{line_number}: invalid cycle counts "
                    f"{positive_cycles}/{total_cycles}"
                )
            if not 0.0 <= empirical_frequency <= 1.0:
                raise ValueError(f"{input_path}:{line_number}: frequency must be between 0 and 1")

            expected_frequency = positive_cycles / total_cycles
            if abs(empirical_frequency - expected_frequency) > 1e-9:
                raise ValueError(
                    f"{input_path}:{line_number}: frequency {empirical_frequency} does not match "
                    f"cycle counts {positive_cycles}/{total_cycles}"
                )

            hypergeometric_pvalue = None
            if pvalue_text is not None:
                try:
                    hypergeometric_pvalue = float(pvalue_text)
                except ValueError as error:
                    raise ValueError(
                        f"{input_path}:{line_number}: hypergeometric p-value is not numeric"
                    ) from error
                if not 0.0 <= hypergeometric_pvalue <= 1.0:
                    raise ValueError(
                        f"{input_path}:{line_number}: hypergeometric p-value must be between 0 and 1"
                    )

            cycles = tuple(cycle for cycle in cycles_text.split(",") if cycle)
            if len(cycles) != positive_cycles:
                raise ValueError(
                    f"{input_path}:{line_number}: positive-cycle count is {positive_cycles}, "
                    f"but {len(cycles)} cycle IDs were listed"
                )
            if len(cycles) != len(set(cycles)):
                raise ValueError(f"{input_path}:{line_number}: duplicate cycle IDs in one position")
            invalid_cycles = [cycle for cycle in cycles if not CYCLE_PATTERN.fullmatch(cycle)]
            if invalid_cycles:
                raise ValueError(
                    f"{input_path}:{line_number}: invalid bootstrap cycle ID(s): "
                    + ",".join(invalid_cycles)
                )

            if observed_gene is None:
                observed_gene = gene
            elif gene != observed_gene:
                raise ValueError(
                    f"{input_path}:{line_number}: multiple gene identifiers in one result: "
                    f"{observed_gene!r} and {gene!r}"
                )
            if position in observed_positions:
                raise ValueError(f"{input_path}:{line_number}: duplicate position {position}")
            observed_positions.add(position)

            records.append(
                BootstrapRecord(
                    fields=tuple(fields),
                    gene=gene,
                    position=position,
                    positive_cycles=positive_cycles,
                    total_cycles=total_cycles,
                    frequency=empirical_frequency,
                    cycles=cycles,
                    hypergeometric_pvalue=hypergeometric_pvalue,
                    line_number=line_number,
                )
            )

    return records


def load_cycle_burden_reference(
    reference_path: Path,
) -> dict[str, CycleBurdenThreshold]:
    """Load and validate the reference produced by the companion builder."""

    required_fields = {
        "cycle_id",
        "quantile",
        "burden_threshold",
        "genes_with_cycle",
    }
    thresholds: dict[str, CycleBurdenThreshold] = {}

    with reference_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing_fields = required_fields.difference(reader.fieldnames or [])
        if missing_fields:
            raise ValueError(
                f"{reference_path}: missing reference column(s): "
                + ",".join(sorted(missing_fields))
            )

        for line_number, row in enumerate(reader, start=2):
            cycle_id = row["cycle_id"]
            if not CYCLE_PATTERN.fullmatch(cycle_id):
                raise ValueError(
                    f"{reference_path}:{line_number}: invalid cycle ID {cycle_id!r}"
                )
            if cycle_id in thresholds:
                raise ValueError(
                    f"{reference_path}:{line_number}: duplicate cycle {cycle_id}"
                )
            try:
                reference_quantile = float(row["quantile"])
                burden_threshold = float(row["burden_threshold"])
                genes_with_cycle = int(row["genes_with_cycle"])
            except ValueError as error:
                raise ValueError(
                    f"{reference_path}:{line_number}: non-numeric reference value"
                ) from error

            _validate_fraction(reference_quantile, "reference quantile")
            _validate_fraction(burden_threshold, "burden threshold")
            if genes_with_cycle < 1:
                raise ValueError(
                    f"{reference_path}:{line_number}: genes_with_cycle must be positive"
                )

            thresholds[cycle_id] = CycleBurdenThreshold(
                cycle_id=cycle_id,
                quantile=reference_quantile,
                burden_threshold=burden_threshold,
                genes_with_cycle=genes_with_cycle,
            )

    if not thresholds:
        raise ValueError(f"{reference_path}: reference contains no cycle thresholds")
    return thresholds


def find_globally_flagged_cycles(
    records: Sequence[BootstrapRecord],
    thresholds: dict[str, CycleBurdenThreshold],
    min_cycle_positions: int = 10,
    min_cycle_dominance: float = 0.8,
) -> list[FlaggedCycle]:
    """Select recurrent gene-cycle associations using genomic cutoffs."""

    if min_cycle_positions < 1:
        raise ValueError("min_cycle_positions must be at least 1")
    _validate_fraction(min_cycle_dominance, "min_cycle_dominance")

    positive_records = [record for record in records if record.is_positive]
    if not positive_records:
        return []

    positions_by_cycle: dict[str, list[int]] = {}
    for record in positive_records:
        for cycle_id in record.cycles:
            positions_by_cycle.setdefault(cycle_id, []).append(record.position)

    missing_cycles = set(positions_by_cycle).difference(thresholds)
    if missing_cycles:
        raise ValueError(
            "cycle-burden reference lacks threshold(s) for: "
            + ",".join(sorted(missing_cycles))
        )

    flagged = []
    gene = positive_records[0].gene
    for cycle_id, positions in positions_by_cycle.items():
        occurrences = len(positions)
        dominance = occurrences / len(positive_records)
        burden = occurrences / len(records)
        threshold = thresholds[cycle_id]
        if (
            occurrences >= min_cycle_positions
            and _meets_fraction(occurrences, len(positive_records), min_cycle_dominance)
            and burden + 1e-12 >= threshold.burden_threshold
        ):
            flagged.append(
                FlaggedCycle(
                    gene=gene,
                    cycle_id=cycle_id,
                    occurrences=occurrences,
                    positive_positions=len(positive_records),
                    total_positions=len(records),
                    dominance=dominance,
                    burden=burden,
                    reference_quantile=threshold.quantile,
                    burden_threshold=threshold.burden_threshold,
                    positions=tuple(positions),
                )
            )

    return sorted(flagged, key=lambda item: int(item.cycle_id.split("_", 1)[1]))


def remove_flagged_cycles(
    records: Sequence[BootstrapRecord],
    flagged_cycles: Sequence[FlaggedCycle],
) -> list[BootstrapRecord]:
    """Remove flagged IDs and recalculate each affected CAAStools row."""

    flagged_ids = {flagged.cycle_id for flagged in flagged_cycles}
    if not flagged_ids:
        return list(records)

    adjusted_records = []
    for record in records:
        retained_cycles = tuple(
            cycle_id for cycle_id in record.cycles if cycle_id not in flagged_ids
        )
        if retained_cycles == record.cycles:
            adjusted_records.append(record)
            continue

        positive_cycles = len(retained_cycles)
        frequency = positive_cycles / record.total_cycles
        adjusted_fields = (
            record.fields[0],
            str(positive_cycles),
            record.fields[2],
            str(frequency),
            ",".join(retained_cycles),
            *record.fields[5:],
        )
        adjusted_records.append(
            replace(
                record,
                fields=adjusted_fields,
                positive_cycles=positive_cycles,
                frequency=frequency,
                cycles=retained_cycles,
            )
        )

    return adjusted_records


def _meets_fraction(numerator: int, denominator: int, threshold: float) -> bool:
    """Compare a fraction to a threshold while avoiding boundary round-off."""

    return numerator / denominator + 1e-12 >= threshold


def find_suspicious_intervals(
    positive_records: Sequence[BootstrapRecord],
    min_span: int = 10,
    min_density: float = 0.8,
    min_cycle_recurrence: float = 0.8,
) -> list[SuspiciousInterval]:
    """Find maximal pairwise intervals satisfying density and cycle recurrence.

    For every possible positive-position start, the furthest qualifying endpoint
    is retained.  Intervals wholly contained in a previously retained interval
    are then omitted from the report because they do not remove any additional
    positions.  Partially overlapping intervals remain separate report entries.
    """

    if min_span < 1:
        raise ValueError("min_span must be at least 1")
    _validate_fraction(min_density, "min_density")
    _validate_fraction(min_cycle_recurrence, "min_cycle_recurrence")

    ordered = sorted(positive_records, key=lambda record: record.position)
    if not ordered:
        return []

    candidates: list[SuspiciousInterval] = []

    for start_index, start_record in enumerate(ordered):
        cycle_counts: Counter[str] = Counter()
        furthest_interval: SuspiciousInterval | None = None

        for end_index in range(start_index, len(ordered)):
            end_record = ordered[end_index]
            cycle_counts.update(end_record.cycles)

            span = end_record.position - start_record.position + 1
            number_of_caas_positions = end_index - start_index + 1
            if span < min_span:
                continue
            if not _meets_fraction(number_of_caas_positions, span, min_density):
                continue

            recurrent_cycles = tuple(
                sorted(
                    cycle
                    for cycle, count in cycle_counts.items()
                    if _meets_fraction(count, number_of_caas_positions, min_cycle_recurrence)
                )
            )
            if not recurrent_cycles:
                continue

            furthest_interval = SuspiciousInterval(
                gene=start_record.gene,
                start_position=start_record.position,
                end_position=end_record.position,
                caas_positions=tuple(
                    record.position for record in ordered[start_index : end_index + 1]
                ),
                recurrent_cycles=recurrent_cycles,
                cycle_counts=tuple(sorted(cycle_counts.items())),
            )

        if furthest_interval is not None:
            candidates.append(furthest_interval)

    maximal_intervals: list[SuspiciousInterval] = []
    furthest_reported_end = -1
    for interval in candidates:
        if interval.end_position <= furthest_reported_end:
            continue
        maximal_intervals.append(interval)
        furthest_reported_end = interval.end_position

    return maximal_intervals


def filter_records(
    records: Sequence[BootstrapRecord],
    min_span: int = 10,
    min_density: float = 0.8,
    min_cycle_recurrence: float = 0.8,
) -> tuple[list[BootstrapRecord], list[SuspiciousInterval], FilterSummary]:
    """Apply the zero-cycle and suspicious-cluster filters."""

    positive_records = [record for record in records if record.is_positive]
    intervals = find_suspicious_intervals(
        positive_records,
        min_span=min_span,
        min_density=min_density,
        min_cycle_recurrence=min_cycle_recurrence,
    )
    cluster_positions = {
        position
        for interval in intervals
        for position in interval.caas_positions
    }
    retained_records = [
        record for record in positive_records if record.position not in cluster_positions
    ]

    summary = FilterSummary(
        total_positions=len(records),
        zero_cycle_positions=len(records) - len(positive_records),
        positive_positions=len(positive_records),
        suspicious_intervals=len(intervals),
        cluster_positions_removed=len(cluster_positions),
        retained_positions=len(retained_records),
    )
    return retained_records, intervals, summary


def filter_records_with_global_reference(
    records: Sequence[BootstrapRecord],
    thresholds: dict[str, CycleBurdenThreshold],
    min_cycle_positions: int = 10,
    min_cycle_dominance: float = 0.8,
    min_span: int = 10,
    min_density: float = 0.8,
    min_cycle_recurrence: float = 0.8,
) -> tuple[
    list[BootstrapRecord],
    list[SuspiciousInterval],
    list[FlaggedCycle],
    CombinedFilterSummary,
]:
    """Apply global cycle pruning, zero filtering, and local cluster filtering."""

    input_positive_positions = sum(record.is_positive for record in records)
    flagged_cycles = find_globally_flagged_cycles(
        records,
        thresholds,
        min_cycle_positions=min_cycle_positions,
        min_cycle_dominance=min_cycle_dominance,
    )
    adjusted_records = remove_flagged_cycles(records, flagged_cycles)
    adjusted_positive_positions = sum(record.is_positive for record in adjusted_records)
    retained_records, intervals, local_summary = filter_records(
        adjusted_records,
        min_span=min_span,
        min_density=min_density,
        min_cycle_recurrence=min_cycle_recurrence,
    )

    summary = CombinedFilterSummary(
        total_positions=len(records),
        input_zero_cycle_positions=len(records) - input_positive_positions,
        input_positive_positions=input_positive_positions,
        globally_flagged_cycles=len(flagged_cycles),
        cycle_associations_removed=sum(item.occurrences for item in flagged_cycles),
        positions_emptied_by_global_filter=(
            input_positive_positions - adjusted_positive_positions
        ),
        positive_positions_after_global_filter=adjusted_positive_positions,
        suspicious_intervals=local_summary.suspicious_intervals,
        cluster_positions_removed=local_summary.cluster_positions_removed,
        retained_positions=local_summary.retained_positions,
    )
    return retained_records, intervals, flagged_cycles, summary


def write_filtered_result(records: Iterable[BootstrapRecord], output_path: Path) -> None:
    """Write retained rows in the original headerless CAAStools format."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        for record in records:
            handle.write("\t".join(record.fields) + "\n")


def write_cluster_report(intervals: Sequence[SuspiciousInterval], report_path: Path) -> None:
    """Write an auditable report describing every maximal suspicious interval."""

    report_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "cluster_id",
        "gene",
        "start_position",
        "end_position",
        "span_aa",
        "caas_positions",
        "caas_density",
        "recurrent_cycles",
        "maximum_cycle_occurrences",
        "maximum_cycle_fraction",
        "positions_removed",
    ]

    with report_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for cluster_number, interval in enumerate(intervals, start=1):
            writer.writerow(
                {
                    "cluster_id": f"cluster_{cluster_number}",
                    "gene": interval.gene,
                    "start_position": interval.start_position,
                    "end_position": interval.end_position,
                    "span_aa": interval.span,
                    "caas_positions": len(interval.caas_positions),
                    "caas_density": f"{interval.density:.6f}",
                    "recurrent_cycles": ",".join(interval.recurrent_cycles),
                    "maximum_cycle_occurrences": interval.maximum_cycle_count,
                    "maximum_cycle_fraction": f"{interval.maximum_cycle_fraction:.6f}",
                    "positions_removed": ",".join(map(str, interval.caas_positions)),
                }
            )


def write_flagged_cycle_report(
    flagged_cycles: Sequence[FlaggedCycle],
    report_path: Path,
) -> None:
    """Write the gene-cycle associations removed by the global filter."""

    report_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "gene",
        "cycle_id",
        "cycle_positions",
        "positive_positions_before_filter",
        "total_positions",
        "cycle_dominance",
        "cycle_burden",
        "reference_quantile",
        "burden_threshold",
        "positions_affected",
    ]

    with report_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for flagged in flagged_cycles:
            writer.writerow(
                {
                    "gene": flagged.gene,
                    "cycle_id": flagged.cycle_id,
                    "cycle_positions": flagged.occurrences,
                    "positive_positions_before_filter": flagged.positive_positions,
                    "total_positions": flagged.total_positions,
                    "cycle_dominance": f"{flagged.dominance:.6f}",
                    "cycle_burden": f"{flagged.burden:.6f}",
                    "reference_quantile": f"{flagged.reference_quantile:.6f}",
                    "burden_threshold": f"{flagged.burden_threshold:.6f}",
                    "positions_affected": ",".join(map(str, flagged.positions)),
                }
            )


def default_output_path(input_path: Path) -> Path:
    return input_path.with_name(input_path.stem + ".filtered.tsv")


def default_report_path(input_path: Path) -> Path:
    return input_path.with_name(input_path.stem + ".clusters.tsv")


def default_flagged_cycle_report_path(input_path: Path) -> Path:
    return input_path.with_name(input_path.stem + ".flagged-cycles.tsv")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Remove anomalously recurrent cycles, zero-cycle rows, and dense local "
            "clusters from one CAAStools bootstrap result."
        )
    )
    parser.add_argument("input", type=Path, help="Headerless CAAStools bootstrap TSV")
    parser.add_argument("-o", "--output", type=Path, help="Filtered TSV output path")
    parser.add_argument(
        "--cluster-report",
        type=Path,
        help="TSV report of suspicious intervals (default: insert .clusters before .tsv)",
    )
    parser.add_argument(
        "--cycle-burden-reference",
        type=Path,
        help="Cycle-specific genomic thresholds from build_cycle_burden_reference.py",
    )
    parser.add_argument(
        "--flagged-cycle-report",
        type=Path,
        help="TSV report of globally removed cycles (default: insert .flagged-cycles before .tsv)",
    )
    parser.add_argument(
        "--min-cycle-positions",
        type=int,
        default=10,
        help="Minimum positions required to flag a globally recurrent cycle (default: 10)",
    )
    parser.add_argument(
        "--min-global-cycle-dominance",
        type=fraction,
        default=0.8,
        help="Minimum fraction of positive positions carrying a global cycle (default: 0.8)",
    )
    parser.add_argument(
        "--min-span",
        type=int,
        default=10,
        help="Minimum inclusive interval length in amino acids (default: 10)",
    )
    parser.add_argument(
        "--min-density",
        type=fraction,
        default=0.8,
        help="Minimum CAAS-positive coordinate density (default: 0.8)",
    )
    parser.add_argument(
        "--min-cycle-recurrence",
        type=fraction,
        default=0.8,
        help="Minimum fraction of CAAS positions sharing one cycle (default: 0.8)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    arguments = parser.parse_args(argv)

    if arguments.min_span < 1:
        parser.error("--min-span must be at least 1")
    if arguments.min_cycle_positions < 1:
        parser.error("--min-cycle-positions must be at least 1")
    if not arguments.input.is_file():
        parser.error(f"input file not found: {arguments.input}")
    if (
        arguments.cycle_burden_reference is not None
        and not arguments.cycle_burden_reference.is_file()
    ):
        parser.error(
            f"cycle-burden reference not found: {arguments.cycle_burden_reference}"
        )
    if arguments.flagged_cycle_report and not arguments.cycle_burden_reference:
        parser.error("--flagged-cycle-report requires --cycle-burden-reference")

    output_path = arguments.output or default_output_path(arguments.input)
    report_path = arguments.cluster_report or default_report_path(arguments.input)
    flagged_report_path = (
        arguments.flagged_cycle_report
        or default_flagged_cycle_report_path(arguments.input)
    )
    input_resolved = arguments.input.resolve()
    if output_path.resolve() == input_resolved or report_path.resolve() == input_resolved:
        parser.error("output and cluster report must not overwrite the input file")
    if output_path.resolve() == report_path.resolve():
        parser.error("output and cluster report must be different files")
    if arguments.cycle_burden_reference:
        output_paths = {
            output_path.resolve(),
            report_path.resolve(),
            flagged_report_path.resolve(),
        }
        if input_resolved in output_paths:
            parser.error("no output report may overwrite the input file")
        if len(output_paths) != 3:
            parser.error("filtered output, cluster report, and flagged-cycle report must differ")

    try:
        records = parse_bootstrap_result(arguments.input)
        if arguments.cycle_burden_reference:
            thresholds = load_cycle_burden_reference(arguments.cycle_burden_reference)
            retained_records, intervals, flagged_cycles, combined_summary = (
                filter_records_with_global_reference(
                    records,
                    thresholds,
                    min_cycle_positions=arguments.min_cycle_positions,
                    min_cycle_dominance=arguments.min_global_cycle_dominance,
                    min_span=arguments.min_span,
                    min_density=arguments.min_density,
                    min_cycle_recurrence=arguments.min_cycle_recurrence,
                )
            )
            write_flagged_cycle_report(flagged_cycles, flagged_report_path)
        else:
            retained_records, intervals, summary = filter_records(
                records,
                min_span=arguments.min_span,
                min_density=arguments.min_density,
                min_cycle_recurrence=arguments.min_cycle_recurrence,
            )
        write_filtered_result(retained_records, output_path)
        write_cluster_report(intervals, report_path)
    except (OSError, ValueError) as error:
        parser.exit(2, f"ERROR: {error}\n")

    print(f"Input: {arguments.input}")
    if arguments.cycle_burden_reference:
        print(f"Total positions: {combined_summary.total_positions}")
        print(f"Input zero-cycle positions: {combined_summary.input_zero_cycle_positions}")
        print(f"Input positive positions: {combined_summary.input_positive_positions}")
        print(f"Globally flagged cycles: {combined_summary.globally_flagged_cycles}")
        print(
            "Cycle-position associations removed: "
            f"{combined_summary.cycle_associations_removed}"
        )
        print(
            "Positions emptied by global cycle filtering: "
            f"{combined_summary.positions_emptied_by_global_filter}"
        )
        print(
            "Positive positions after global filtering: "
            f"{combined_summary.positive_positions_after_global_filter}"
        )
        print(f"Suspicious intervals: {combined_summary.suspicious_intervals}")
        print(
            "Positive positions removed in suspicious intervals: "
            f"{combined_summary.cluster_positions_removed}"
        )
        print(f"Positions retained: {combined_summary.retained_positions}")
        print(f"Flagged-cycle report: {flagged_report_path}")
    else:
        print(f"Total positions: {summary.total_positions}")
        print(f"Zero-cycle positions removed: {summary.zero_cycle_positions}")
        print(f"Positive positions before cluster filtering: {summary.positive_positions}")
        print(f"Suspicious intervals: {summary.suspicious_intervals}")
        print(
            "Positive positions removed in suspicious intervals: "
            f"{summary.cluster_positions_removed}"
        )
        print(f"Positions retained: {summary.retained_positions}")
    print(f"Filtered output: {output_path}")
    print(f"Cluster report: {report_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
