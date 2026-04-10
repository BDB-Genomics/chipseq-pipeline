#!/usr/bin/env python3
"""Parse and validate ChIP-seq QC metrics in a format-tolerant way."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


FRIP_LABELS = {"frip"}
PHANTOM_METRIC_NAMES = {"NSC", "RSC"}
SAMTOOLS_METRIC_LABELS = {
    "total_reads": {
        "raw total sequences",
        "sequences",
    },
    "mapping_rate": {
        "percentage of properly paired reads",
        "percentage of reads mapped",
        "percentage of reads mapped and paired",
    },
    "duplicates": {
        "reads duplicated",
    },
}


def read_nonempty_lines(path: str | Path):
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped:
                yield stripped


def parse_frip(frip_path: str | Path) -> float | None:
    """Parse FRiP from a small tabular report."""
    try:
        for line in read_nonempty_lines(frip_path):
            parts = re.split(r"\s+", line)
            if len(parts) < 2:
                continue
            label = parts[0].strip().lower().rstrip(":")
            if label in FRIP_LABELS:
                return float(parts[1])
    except Exception as exc:  # noqa: BLE001
        print(f"Error parsing FRiP: {exc}", file=sys.stderr)
    return None


def parse_phantom(phantom_path: str | Path) -> tuple[float | None, float | None]:
    """Parse NSC and RSC from PhantomPeakQualTools output.

    PhantomPeakQualTools is typically an 11-column tab-delimited report. This parser
    first looks for an explicit header and then falls back to the documented columns.
    """
    try:
        lines = list(read_nonempty_lines(phantom_path))
        if not lines:
            return None, None

        header_tokens = re.split(r"\s+", lines[0].lstrip("#"))
        header_map = {token.upper(): idx for idx, token in enumerate(header_tokens)}
        if PHANTOM_METRIC_NAMES.issubset(header_map):
            data_row = re.split(r"\s+", lines[1]) if len(lines) > 1 else []
            nsc_idx = header_map["NSC"]
            rsc_idx = header_map["RSC"]
            if len(data_row) > max(nsc_idx, rsc_idx):
                return float(data_row[nsc_idx]), float(data_row[rsc_idx])

        for line in lines:
            parts = re.split(r"\t+", line)
            if len(parts) >= 10:
                try:
                    return float(parts[8]), float(parts[9])
                except ValueError:
                    continue
    except Exception as exc:  # noqa: BLE001
        print(f"Error parsing PhantomPeakQualTools: {exc}", file=sys.stderr)
    return None, None


def parse_metric_from_sn_line(line: str) -> float | int | None:
    """Extract the value from a samtools SN line."""
    parts = line.split("\t")
    if len(parts) < 3:
        return None

    raw_value = parts[2].strip().rstrip("%")
    try:
        if "." in raw_value:
            return float(raw_value)
        return int(raw_value)
    except ValueError:
        return None


def parse_samtools_stats(stats_path: str | Path) -> dict[str, float | int | None]:
    """Parse samtools stats by label, not by fixed row number."""
    metrics = {
        "total_reads": None,
        "mapping_rate": None,
        "duplicates": None,
    }

    try:
        for line in read_nonempty_lines(stats_path):
            if not line.startswith("SN\t"):
                continue

            label = line.split("\t", 2)[1].strip().lower().rstrip(":")
            value = parse_metric_from_sn_line(line)
            if value is None:
                continue

            if label in SAMTOOLS_METRIC_LABELS["total_reads"]:
                metrics["total_reads"] = int(value)
            elif label in SAMTOOLS_METRIC_LABELS["mapping_rate"]:
                metrics["mapping_rate"] = float(value)
            elif label in SAMTOOLS_METRIC_LABELS["duplicates"]:
                metrics["duplicates"] = int(value)
    except Exception as exc:  # noqa: BLE001
        print(f"Error parsing samtools stats: {exc}", file=sys.stderr)

    return metrics


def safe_float(value, label: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} is not numeric: {value!r}") from exc


def write_log(log_path: str | Path, lines: list[str]) -> None:
    with open(log_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Parse and validate ChIP-seq QC metrics.")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--frip-file", required=True)
    parser.add_argument("--phantom-file", required=True)
    parser.add_argument("--stats-file", required=True)
    parser.add_argument("--min-frip", type=float, required=True)
    parser.add_argument("--min-nsc", type=float, required=True)
    parser.add_argument("--min-rsc", type=float, required=True)
    parser.add_argument("--min-mapping-rate", type=float, required=True)
    parser.add_argument("--max-duplicate-rate", type=float, required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    log_content = [f"QC Report for {args.sample}", "-------------------------------"]

    frip = parse_frip(args.frip_file)
    nsc, rsc = parse_phantom(args.phantom_file)
    stats = parse_samtools_stats(args.stats_file)

    parse_failures = []
    if frip is None:
        parse_failures.append("FRiP")
    if nsc is None or rsc is None:
        parse_failures.append("NSC/RSC")
    if any(value is None for value in stats.values()):
        parse_failures.append("Samtools Stats")

    if parse_failures:
        log_content.append(f"[ERROR] Failed to parse metrics: {', '.join(parse_failures)}")
        write_log(args.log, log_content)
        sys.exit(1)

    frip = safe_float(frip, "FRiP")
    nsc = safe_float(nsc, "NSC")
    rsc = safe_float(rsc, "RSC")
    total_reads = int(stats["total_reads"])
    duplicate_reads = int(stats["duplicates"])
    mapping_rate = safe_float(stats["mapping_rate"], "Mapping rate")

    dup_rate = 100.0 if total_reads <= 0 else (duplicate_reads * 100.0 / total_reads)

    log_content.extend(
        [
            f"FRiP: {frip:.4f} (Target: >= {args.min_frip})",
            f"NSC: {nsc:.4f} (Target: >= {args.min_nsc})",
            f"RSC: {rsc:.4f} (Target: >= {args.min_rsc})",
            f"Mapping Rate (%): {mapping_rate:.2f} (Target: >= {args.min_mapping_rate})",
            f"Duplicate Rate (%): {dup_rate:.2f} (Target: <= {args.max_duplicate_rate})",
        ]
    )

    failures: list[str] = []
    if frip < args.min_frip:
        failures.append(f"FRiP {frip:.4f} < {args.min_frip}")
    if nsc < args.min_nsc:
        failures.append(f"NSC {nsc:.4f} < {args.min_nsc}")
    if rsc < args.min_rsc:
        failures.append(f"RSC {rsc:.4f} < {args.min_rsc}")
    if mapping_rate < args.min_mapping_rate:
        failures.append(f"Mapping Rate {mapping_rate:.2f} < {args.min_mapping_rate}")
    if dup_rate > args.max_duplicate_rate:
        failures.append(f"Duplicate Rate {dup_rate:.2f} > {args.max_duplicate_rate}")

    if failures:
        log_content.append("-------------------------------")
        log_content.append("RESULT: FAILED")
        log_content.extend(f"[QC FAILURE] {failure}" for failure in failures)
        write_log(args.log, log_content)
        print("\n".join(log_content), file=sys.stderr)
        sys.exit(2)

    log_content.append("-------------------------------")
    log_content.append("RESULT: PASSED")
    write_log(args.log, log_content)

    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write(f"{args.sample}\tPASSED\n")


if __name__ == "__main__":
    main()
