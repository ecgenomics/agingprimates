#!/usr/bin/env python3

"""Build cycle-specific genomic burden thresholds from CAAStools bootstrap results.

For each input gene and bootstrap cycle, burden is defined as::

    positions positive for the cycle / positions present in the gene result

The output contains one empirical burden threshold per cycle.  Only genes in
which a cycle occurs contribute to that cycle's conditional distribution.  The
default threshold is its 97.5th percentile (linear interpolation, equivalent
to the common type-7 sample quantile).

The resulting TSV is consumed by ``filter_bootstrap_caas.py``.  Separating
reference construction from per-gene filtering makes the reference fully
auditable and reusable in a later workflow.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
from math import floor
from pathlib import Path
import sys
from typing import Iterable, Sequence

from filter_bootstrap_caas import BootstrapRecord, fraction, parse_bootstrap_result


REFERENCE_FIELDS = [
    "cycle_id",
    "quantile",
    "burden_threshold",
    "genes_with_cycle",
    "median_burden",
    "maximum_burden",
    "input_files",
    "genes_with_positive_cycles",
]


@dataclass(frozen=True)
class CycleReferenceRow:
    """Summary and empirical cutoff for one bootstrap cycle."""

    cycle_id: str
    quantile: float
    burden_threshold: float
    genes_with_cycle: int
    median_burden: float
    maximum_burden: float
    input_files: int
    genes_with_positive_cycles: int


def cycle_sort_key(cycle_id: str) -> tuple[int, str]:
    """Sort canonical b_NUMBER identifiers numerically."""

    try:
        return int(cycle_id.split("_", 1)[1]), cycle_id
    except (IndexError, ValueError):
        return sys.maxsize, cycle_id


def empirical_quantile(values: Sequence[float], quantile: float) -> float:
    """Return a linearly interpolated empirical quantile."""

    if not values:
        raise ValueError("cannot calculate a quantile from an empty sequence")
    if not 0.0 <= quantile <= 1.0:
        raise ValueError("quantile must be between 0 and 1")

    ordered = sorted(values)
    rank = (len(ordered) - 1) * quantile
    lower_index = floor(rank)
    upper_index = min(lower_index + 1, len(ordered) - 1)
    fraction_above_lower = rank - lower_index
    return (
        ordered[lower_index] * (1.0 - fraction_above_lower)
        + ordered[upper_index] * fraction_above_lower
    )


def cycle_counts(records: Iterable[BootstrapRecord]) -> Counter[str]:
    """Count the alignment positions positive for each cycle."""

    return Counter(cycle for record in records for cycle in record.cycles)


def build_reference_rows(
    input_paths: Sequence[Path],
    quantile: float = 0.975,
) -> list[CycleReferenceRow]:
    """Read all genes and calculate one conditional burden cutoff per cycle."""

    if not input_paths:
        raise ValueError("no CAAStools bootstrap result files were provided")
    if not 0.0 <= quantile <= 1.0:
        raise ValueError("quantile must be between 0 and 1")

    burdens_by_cycle: dict[str, list[float]] = defaultdict(list)
    genes_with_positive_cycles = 0

    for input_path in input_paths:
        records = parse_bootstrap_result(input_path)
        if not records:
            continue
        counts = cycle_counts(records)
        if counts:
            genes_with_positive_cycles += 1
        for cycle_id, count in counts.items():
            burdens_by_cycle[cycle_id].append(count / len(records))

    reference_rows = []
    for cycle_id in sorted(burdens_by_cycle, key=cycle_sort_key):
        burdens = burdens_by_cycle[cycle_id]
        reference_rows.append(
            CycleReferenceRow(
                cycle_id=cycle_id,
                quantile=quantile,
                burden_threshold=empirical_quantile(burdens, quantile),
                genes_with_cycle=len(burdens),
                median_burden=empirical_quantile(burdens, 0.5),
                maximum_burden=max(burdens),
                input_files=len(input_paths),
                genes_with_positive_cycles=genes_with_positive_cycles,
            )
        )

    if not reference_rows:
        raise ValueError("none of the input files contains a positive bootstrap cycle")
    return reference_rows


def write_reference(rows: Sequence[CycleReferenceRow], output_path: Path) -> None:
    """Write the cycle reference as a tab-separated table."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REFERENCE_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "cycle_id": row.cycle_id,
                    "quantile": f"{row.quantile:.12g}",
                    "burden_threshold": f"{row.burden_threshold:.12g}",
                    "genes_with_cycle": row.genes_with_cycle,
                    "median_burden": f"{row.median_burden:.12g}",
                    "maximum_burden": f"{row.maximum_burden:.12g}",
                    "input_files": row.input_files,
                    "genes_with_positive_cycles": row.genes_with_positive_cycles,
                }
            )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build cycle-specific genomic burden thresholds from CAAStools bootstrap results."
    )
    parser.add_argument("input_directory", type=Path, help="Directory containing bootstrap TSV files")
    parser.add_argument("-o", "--output", required=True, type=Path, help="Output reference TSV")
    parser.add_argument(
        "--pattern",
        default="*.bootstrap.caas.tsv",
        help="Glob relative to INPUT_DIRECTORY (default: *.bootstrap.caas.tsv)",
    )
    parser.add_argument(
        "--quantile",
        type=fraction,
        default=0.975,
        help="Conditional burden quantile used as cutoff (default: 0.975)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    arguments = parser.parse_args(argv)

    if not arguments.input_directory.is_dir():
        parser.error(f"input directory not found: {arguments.input_directory}")
    input_paths = sorted(
        path for path in arguments.input_directory.glob(arguments.pattern) if path.is_file()
    )
    if not input_paths:
        parser.error(
            f"no files matching {arguments.pattern!r} in {arguments.input_directory}"
        )

    try:
        rows = build_reference_rows(input_paths, quantile=arguments.quantile)
        write_reference(rows, arguments.output)
    except (OSError, ValueError) as error:
        parser.exit(2, f"ERROR: {error}\n")

    print(f"Input files: {len(input_paths)}")
    print(f"Cycles represented: {len(rows)}")
    print(f"Reference quantile: {arguments.quantile}")
    print(f"Reference output: {arguments.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
