#!/usr/bin/env python3

"""Merge CAAStools discovery tables and optionally filter by raw p-value."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Iterable


EXPECTED_HEADER = [
    "Gene",
    "Trait",
    "Position",
    "Substitution",
    "Pvalue",
    "Pattern",
    "FFGN",
    "FBGN",
    "GFG",
    "GBG",
    "MFG",
    "MBG",
    "FFG",
    "FBG",
    "MS",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Merge headered CAAStools discovery tables, retain one header, "
            "and optionally write rows whose Pvalue is at most alpha."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing the input tables staged by Nextflow.",
    )
    parser.add_argument(
        "--suffix",
        required=True,
        help="Only regular files ending with this suffix are merged.",
    )
    parser.add_argument("--merged", type=Path, required=True)
    parser.add_argument("--significant", type=Path)
    parser.add_argument("--alpha", type=float, default=0.05)
    return parser.parse_args()


def discover_inputs(input_dir: Path, suffix: str, outputs: Iterable[Path]) -> list[Path]:
    if not input_dir.is_dir():
        raise ValueError(f"Input directory does not exist: {input_dir}")
    excluded = {path.resolve() for path in outputs}
    inputs = sorted(
        (
            path
            for path in input_dir.iterdir()
            if path.is_file()
            and path.name.endswith(suffix)
            and path.resolve() not in excluded
        ),
        key=lambda path: path.name,
    )
    if not inputs:
        raise ValueError(
            f"No input files ending with {suffix!r} found in {input_dir}"
        )
    return inputs


def rows_from_table(path: Path) -> Iterable[list[str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        first_nonempty_seen = False
        for line_number, row in enumerate(reader, start=1):
            if not row or all(cell == "" for cell in row):
                continue
            if not first_nonempty_seen:
                first_nonempty_seen = True
                if row != EXPECTED_HEADER:
                    raise ValueError(
                        f"Unexpected or missing header in {path} at line {line_number}: {row}"
                    )
                continue
            if row == EXPECTED_HEADER:
                # Be tolerant of previously concatenated files with repeated headers.
                continue
            if len(row) != len(EXPECTED_HEADER):
                raise ValueError(
                    f"Expected {len(EXPECTED_HEADER)} columns in {path} "
                    f"at line {line_number}, found {len(row)}"
                )
            yield row


def parse_pvalue(row: list[str], source: Path) -> float:
    try:
        pvalue = float(row[4])
    except ValueError as error:
        raise ValueError(
            f"Invalid Pvalue {row[4]!r} in {source} for {row[0]}@{row[2]}"
        ) from error
    if not math.isfinite(pvalue) or not 0.0 <= pvalue <= 1.0:
        raise ValueError(
            f"Pvalue outside [0, 1] in {source} for {row[0]}@{row[2]}: {pvalue}"
        )
    return pvalue


def merge_tables(
    inputs: list[Path], merged: Path, significant: Path | None, alpha: float
) -> tuple[int, int]:
    merged.parent.mkdir(parents=True, exist_ok=True)
    if significant is not None:
        significant.parent.mkdir(parents=True, exist_ok=True)

    total_rows = 0
    significant_rows = 0

    with merged.open("w", encoding="utf-8", newline="") as merged_handle:
        merged_writer = csv.writer(
            merged_handle, delimiter="\t", lineterminator="\n"
        )
        merged_writer.writerow(EXPECTED_HEADER)

        significant_handle = None
        significant_writer = None
        try:
            if significant is not None:
                significant_handle = significant.open(
                    "w", encoding="utf-8", newline=""
                )
                significant_writer = csv.writer(
                    significant_handle, delimiter="\t", lineterminator="\n"
                )
                significant_writer.writerow(EXPECTED_HEADER)

            for source in inputs:
                for row in rows_from_table(source):
                    merged_writer.writerow(row)
                    total_rows += 1
                    if significant_writer is not None and parse_pvalue(row, source) <= alpha:
                        significant_writer.writerow(row)
                        significant_rows += 1
        finally:
            if significant_handle is not None:
                significant_handle.close()

    return total_rows, significant_rows


def main() -> None:
    args = parse_args()
    if not 0.0 <= args.alpha <= 1.0:
        raise ValueError(f"--alpha must be in [0, 1], found {args.alpha}")

    output_paths = [args.merged]
    if args.significant is not None:
        output_paths.append(args.significant)
    inputs = discover_inputs(args.input_dir, args.suffix, output_paths)
    total_rows, significant_rows = merge_tables(
        inputs, args.merged, args.significant, args.alpha
    )

    print(f"Input files: {len(inputs)}")
    print(f"Merged CAAS rows: {total_rows}")
    print(f"Merged output: {args.merged}")
    if args.significant is not None:
        print(f"Significant rows (Pvalue <= {args.alpha:g}): {significant_rows}")
        print(f"Significant output: {args.significant}")


if __name__ == "__main__":
    main()
