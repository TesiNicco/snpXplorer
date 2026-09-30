#!/usr/bin/env python3
"""Benchmark annotation-only runs with the multiprocessing script.

Run inside the snpXplorer container:
    python3 /Annotation/BIN/performance/benchmark_annotation.py

Each worker setting gets a separate output directory and log. No email is
provided, so enrichment, plotting, result compression, and mail are skipped.
"""
import argparse
import csv
import os
from pathlib import Path
import resource
import subprocess
import sys
import time
from datetime import datetime, timezone


PERFORMANCE_DIR = Path(__file__).resolve().parent
ANNOTATION_SCRIPT = PERFORMANCE_DIR.parent / "standalone_annotation_optimized_mp.py"
DEFAULT_INPUT = PERFORMANCE_DIR / "snplist_example.txt"
FIELDS = (
    "workers", "wall_s", "user_cpu_s", "system_cpu_s", "cpu_percent",
    "speedup_vs_1", "max_rss_kb_reported", "exit_code", "selected_variants",
    "output_rows", "gnomad_af_present", "clinvar_gene_present",
    "gnomad_responses_ok", "gnomad_failed_attempts", "alphagenome_errors", "ld_errors",
    "log_file", "output_dir",
)


def count_data_rows(path):
    if not path.exists():
        return 0
    with path.open(encoding="utf-8", errors="replace") as stream:
        return max(0, sum(1 for _ in stream) - 1)


def count_present(path, column):
    if not path.exists():
        return 0
    with path.open(newline="", encoding="utf-8", errors="replace") as stream:
        return sum(
            row.get(column) not in (None, "", "NA")
            for row in csv.DictReader(stream, delimiter="\t")
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--workers", type=int, nargs="+", default=[1, 2, 4])
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if not args.input.is_file():
        parser.error(f"Input file not found: {args.input}")
    if not ANNOTATION_SCRIPT.is_file():
        parser.error(f"Annotation script not found: {ANNOTATION_SCRIPT}")
    if not args.workers or any(not 1 <= n <= 4 for n in args.workers):
        parser.error("Worker counts must be between 1 and 4")

    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_root = PERFORMANCE_DIR / f"benchmark_runs_{timestamp}"
    summary = PERFORMANCE_DIR / f"benchmark_results_{timestamp}.tsv"
    commands = []
    for workers in args.workers:
        run_id = f"benchmark_w{workers}"
        output_dir = run_root / f"snpXplorer_annotation_{run_id}"
        log_file = run_root / f"workers_{workers}.log"
        command = [
            sys.executable, str(ANNOTATION_SCRIPT),
            "-q", str(args.input), "-t", "annotation", "-b", "hg38",
            "-o", str(run_root), "-r", run_id,
            "--workers", str(workers), "--gsea-workers", "1",
        ]
        commands.append((workers, output_dir, log_file, command))

    if args.dry_run:
        for workers, _, _, command in commands:
            print(f"{workers} workers: {' '.join(command)}")
        return

    run_root.mkdir(parents=True, exist_ok=False)
    baseline = None
    with summary.open("w", newline="", encoding="utf-8") as table:
        writer = csv.DictWriter(table, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        table.flush()
        for workers, output_dir, log_file, command in commands:
            print(f"Running {workers} workers; log: {log_file}", flush=True)
            with log_file.open("w", encoding="utf-8") as log:
                started = time.monotonic()
                process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
                _, status, usage = os.wait4(process.pid, 0)
                elapsed = time.monotonic() - started
                process.returncode = os.waitstatus_to_exitcode(status)

            if workers == 1 and process.returncode == 0:
                baseline = elapsed
            log_text = log_file.read_text(encoding="utf-8", errors="replace")
            cpu_s = usage.ru_utime + usage.ru_stime
            row = {
                "workers": workers,
                "wall_s": f"{elapsed:.3f}",
                "user_cpu_s": f"{usage.ru_utime:.3f}",
                "system_cpu_s": f"{usage.ru_stime:.3f}",
                "cpu_percent": f"{100 * cpu_s / elapsed:.1f}" if elapsed else "",
                "speedup_vs_1": f"{baseline / elapsed:.2f}" if baseline and process.returncode == 0 else "",
                "max_rss_kb_reported": usage.ru_maxrss,
                "exit_code": process.returncode,
                "selected_variants": count_data_rows(output_dir / "input_selection" / "selected_variants.tsv"),
                "output_rows": count_data_rows(output_dir / "variant_annotation_combined.tsv"),
                "gnomad_af_present": count_present(output_dir / "target_annotations" / "gnomad_annotation_target.tsv", "Total Af"),
                "clinvar_gene_present": count_present(output_dir / "target_annotations" / "clinvar_annotation_target.tsv", "Gene Symbol"),
                "gnomad_responses_ok": log_text.count("[GNOMAD_OK] "),
                "gnomad_failed_attempts": log_text.count("[GNOMAD] ") - log_text.count("[GNOMAD] Annotation stopped:"),
                "alphagenome_errors": log_text.count("Error in alphagenome_annotation:"),
                "ld_errors": log_text.count("Error in query_ld:"),
                "log_file": str(log_file),
                "output_dir": str(output_dir),
            }
            writer.writerow(row)
            table.flush()
            print(
                f"{workers} workers: {elapsed:.1f}s, {row['cpu_percent']}% CPU, "
                f"exit {process.returncode}, {row['selected_variants']} selected variants",
                flush=True,
            )
    print(f"Summary: {summary}", flush=True)
    print("CPU percentage can exceed 100% when multiple cores work simultaneously.", flush=True)


if __name__ == "__main__":
    main()
