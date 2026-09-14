#!/usr/bin/env python3
"""
Lightweight response-time benchmark for snpXplorer.

This script is intentionally separate from the web server. It performs
sequential requests against the production site and reports median/IQR
response time per query type.

Covered query groups
1. Single-variant annotation
2. Haplotype module, variant query
3. Haplotype module, gene query
4. Haplotype module, trait query

Notes
- Annotation timings measure request/submit latency for `/annotation/`.
  That endpoint launches backend work asynchronously, so this does NOT measure
  full job completion time.
- Haplotype timings measure synchronous response time for `/haplotype/`.
- Keep query lists modest and run sequentially to avoid stressing production.
"""

from __future__ import annotations
import argparse
import csv
import json
import statistics
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Tuple
from urllib.parse import urlencode
import requests

# Define the base URL and endpoint paths for the snpXplorer service.
BASE_URL = "https://snpxplorer.net"
ANNOTATION_PATH = "/snpbot/api/variant/query"
HAPLOTYPE_PATH = "/haplotype/"
SCRIPT_DIR = Path(__file__).resolve().parent

# Fill these lists before running.
RSIDS_ANNOTATION: List[str] = ["rs2526377","rs12590654","rs10947943","rs17125924","rs10933431","rs679515","rs1140239","rs7225151","rs113706587","rs13237518","rs17020490","rs62374257","rs7767350","rs7384878","rs2242595","rs3851179","rs2830489","rs10952097","rs602602","rs7908662"]
RSIDS_HAPLOTYPE: List[str] = ["rs2526377","rs12590654","rs10947943","rs17125924","rs10933431","rs679515","rs1140239","rs7225151","rs113706587","rs13237518","rs17020490","rs62374257","rs7767350","rs7384878","rs2242595","rs3851179","rs2830489","rs10952097","rs602602","rs7908662"]
GENES_HAPLOTYPE: List[str] = ["CLU", "PICALM", "CR1", "ABCA7", "BIN1", "CD33", "MS4A6A", "MS4A4E", "EPHA1", "CD2AP", "SORL1", "PTK2B", "NME8", "ZCWPW1", "CELF1", "FERMT2", "SLC24A4", "APOE", "APOC1", "ABCA1"]
TRAITS_HAPLOTYPE: List[str] = ["Alzheimer's disease", "Type 2 diabetes", "Coronary heart disease", "Breast cancer", "Prostate cancer", "Rheumatoid arthritis", "Inflammatory bowel disease", "Asthma", "Schizophrenia", "Bipolar disorder"]

# Use an address you control if you actually want annotation jobs to be accepted.
# The annotation endpoint validates the email before submitting a job.
ANNOTATION_EMAIL = ""

# Define a dataclass to hold timing results for each query.
@dataclass
class TimingResult:
    query_type: str
    query_value: str
    replicate: int
    elapsed_seconds: float
    status_code: int
    ok: bool
    response_verified: bool = False
    verification_detail: str = ""
    request_url: str = ""
    error_message: str = ""

# Function to compute the percentile of a sorted list of values.
def percentile(sorted_values: List[float], p: float) -> float:
    if not sorted_values:
        raise ValueError("Cannot compute percentile of empty list.")
    if len(sorted_values) == 1:
        return sorted_values[0]
    idx = (len(sorted_values) - 1) * p
    lower = int(idx)
    upper = min(lower + 1, len(sorted_values) - 1)
    weight = idx - lower
    return sorted_values[lower] * (1 - weight) + sorted_values[upper] * weight

# Function to summarize a list of timing values, returning key statistics.
def summarize(values: List[float]) -> dict:
    ordered = sorted(values)
    return {
        "n": len(ordered),
        "median": statistics.median(ordered),
        "q1": percentile(ordered, 0.25),
        "q3": percentile(ordered, 0.75),
        "min": ordered[0],
        "max": ordered[-1],
    }

# Function to perform a POST request and measure the elapsed time.
def post_timed(session: requests.Session, url: str, data: dict, timeout: float) -> Tuple[float, requests.Response]:
    start = time.perf_counter()
    response = session.post(url, data=data, timeout=timeout)
    elapsed = time.perf_counter() - start
    return elapsed, response

# Function to perform a GET request and measure the elapsed time.
def get_timed(session: requests.Session, url: str, params: dict, timeout: float) -> Tuple[float, requests.Response]:
    start = time.perf_counter()
    response = session.get(url, params=params, timeout=timeout)
    elapsed = time.perf_counter() - start
    return elapsed, response

# Function to verify that the SNPbot API response looks like a successful query evaluation.
def verify_snpbot_response(response: requests.Response) -> Tuple[bool, str]:
    if not response.ok:
        return False, f"http_{response.status_code}"
    try:
        payload = response.json()
    except json.JSONDecodeError:
        return False, "non_json_response"
    if not isinstance(payload, dict):
        return False, "json_not_dict"
    if payload.get("error"):
        return False, f"api_error:{payload.get('error')}"
    expected_keys = {"query", "build", "Chromosome", "Position_hg38", "rsid", "chr_hg38", "pos_hg38"}
    if any(key in payload for key in expected_keys):
        return True, "json_payload_present"
    return False, "json_missing_expected_keys"

# Function to verify that the haplotype HTML response looks like a successful query evaluation.
def verify_haplotype_response(response: requests.Response, browse_type: str) -> Tuple[bool, str]:
    if not response.ok:
        return False, f"http_{response.status_code}"
    text = response.text
    if browse_type == "trait":
        markers = [
            "Traits related to your query",
            "Haplotype Associations for Selected Traits",
            "GWAS Associations for Selected Traits",
        ]
    else:
        markers = [
            "Graphical overview of haplotypes",
            "Haplotypes in this region",
            "Average effect of haplotypes on phenotypes",
            "No haplotype found for the SNP of interest",
            "Haplotype of the SNP of interest is reported in the first row.",
        ]
    for marker in markers:
        if marker in text:
            return True, f"found_marker:{marker}"
    return False, "expected_marker_not_found"

# Function to benchmark single-variant annotation queries.
def benchmark_annotation_rsids(session: requests.Session, rsids: Iterable[str], repeats: int, timeout: float) -> List[TimingResult]:
    results: List[TimingResult] = []
    for rsid in rsids:
        params = {
            "query": rsid,
            "build": "hg38",
        }
        for rep in range(1, repeats + 1):
            request_url = f"{BASE_URL}{ANNOTATION_PATH}?{urlencode(params)}"
            try:
                elapsed, response = get_timed(session, f"{BASE_URL}{ANNOTATION_PATH}", params, timeout)
                verified, detail = verify_snpbot_response(response)
                results.append(
                    TimingResult(
                        query_type="annotation_single_variant",
                        query_value=rsid,
                        replicate=rep,
                        elapsed_seconds=elapsed,
                        status_code=response.status_code,
                        ok=response.ok,
                        response_verified=verified,
                        verification_detail=detail,
                        request_url=request_url,
                    )
                )
            except requests.RequestException as exc:
                results.append(
                    TimingResult(
                        query_type="annotation_single_variant",
                        query_value=rsid,
                        replicate=rep,
                        elapsed_seconds=float("nan"),
                        status_code=0,
                        ok=False,
                        response_verified=False,
                        verification_detail="request_exception",
                        request_url=request_url,
                        error_message=str(exc),
                    )
                )
    return results

# Function to benchmark haplotype queries for different browse types (rsid, gene, trait).
def benchmark_haplotype_queries(session: requests.Session, queries: Iterable[str], browse_type: str, repeats: int, timeout: float) -> List[TimingResult]:
    label_map = {
        "rsid": "haplotype_variant",
        "gene": "haplotype_gene",
        "trait": "haplotype_trait",
    }
    if browse_type not in label_map:
        raise ValueError(f"Unsupported browse_type: {browse_type}")
    results: List[TimingResult] = []
    for query in queries:
        payload = {
            "browse": query,
            "browse_type": browse_type,
        }
        for rep in range(1, repeats + 1):
            request_url = f"{BASE_URL}{HAPLOTYPE_PATH}"
            try:
                elapsed, response = post_timed(session, f"{BASE_URL}{HAPLOTYPE_PATH}", payload, timeout)
                verified, detail = verify_haplotype_response(response, browse_type)
                results.append(
                    TimingResult(
                        query_type=label_map[browse_type],
                        query_value=query,
                        replicate=rep,
                        elapsed_seconds=elapsed,
                        status_code=response.status_code,
                        ok=response.ok,
                        response_verified=verified,
                        verification_detail=detail,
                        request_url=request_url,
                    )
                )
            except requests.RequestException as exc:
                results.append(
                    TimingResult(
                        query_type=label_map[browse_type],
                        query_value=query,
                        replicate=rep,
                        elapsed_seconds=float("nan"),
                        status_code=0,
                        ok=False,
                        response_verified=False,
                        verification_detail="request_exception",
                        request_url=request_url,
                        error_message=str(exc),
                    )
                )
    return results

# Function to write raw timing results to a CSV file.
def write_raw_results(results: List[TimingResult], outpath: Path) -> None:
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with outpath.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["query_type", "query_value", "replicate", "elapsed_seconds", "status_code", "ok", "response_verified", "verification_detail", "request_url", "error_message"])
        for row in results:
            writer.writerow(
                [
                    row.query_type,
                    row.query_value,
                    row.replicate,
                    f"{row.elapsed_seconds:.6f}",
                    row.status_code,
                    row.ok,
                    row.response_verified,
                    row.verification_detail,
                    row.request_url,
                    row.error_message,
                ]
            )

# Function to print a summary of timing results, including median and IQR.
def print_summary(results: List[TimingResult], title: str = "Response time summary") -> None:
    grouped: dict[str, List[float]] = {}
    for row in results:
        if row.ok and not (isinstance(row.elapsed_seconds, float) and row.elapsed_seconds != row.elapsed_seconds):
            grouped.setdefault(row.query_type, []).append(row.elapsed_seconds)
    if not grouped:
        print(f"\n{title}")
        print("-" * len(title))
        print("No successful results collected.")
        return
    print(f"\n{title}")
    print("-" * len(title))
    for query_type, values in grouped.items():
        stats = summarize(values)
        verified_total = sum(1 for row in results if row.query_type == query_type and row.response_verified)
        total = sum(1 for row in results if row.query_type == query_type)
        print(
            f"{query_type}: "
            f"median={stats['median']:.3f}s, "
            f"IQR={stats['q1']:.3f}-{stats['q3']:.3f}s, "
            f"n={stats['n']}, "
            f"verified={verified_total}/{total}"
        )

# Function to run a benchmark block and report results, including writing to CSV and printing summary.
def run_and_report_block(block_name: str, results: List[TimingResult], output_dir: Path, n_queries: int, repeats: int) -> None:
    outpath = output_dir / f"{block_name}_{n_queries}_{repeats}x.csv"
    write_raw_results(results, outpath)
    print_summary(results, title=f"{block_name} summary")
    failures = sum(1 for row in results if not row.ok)
    if failures:
        print(f"Failed requests: {failures} / {len(results)}")
    print(f"Raw results written to: {outpath}")

# Main function to parse arguments, run benchmarks, and output results.
def main() -> None:
    parser = argparse.ArgumentParser(description="Sequential response-time benchmark for snpXplorer.")
    parser.add_argument(
        "mode",
        choices=["rsid_annotation", "rsid_haplo", "gene", "trait"],
        help="Which benchmark modality to run.",
    )
    parser.add_argument("--repeats", type=int, default=3, help="Number of repeats per query.")
    parser.add_argument("--timeout", type=float, default=120.0, help="Per-request timeout in seconds.")
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.0,
        help="Optional sleep in seconds between requests to reduce production load.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=SCRIPT_DIR / "results",
        help="Output directory for per-modality CSV files.",
    )
    args = parser.parse_args()

    mode_to_block = {
        "rsid_annotation": ("rsids_annotation", RSIDS_ANNOTATION, "annotation"),
        "rsid_haplo": ("rsids_haplotype", RSIDS_HAPLOTYPE, "rsid"),
        "gene": ("genes_haplotype", GENES_HAPLOTYPE, "gene"),
        "trait": ("traits_haplotype", TRAITS_HAPLOTYPE, "trait"),
    }
    block_name, queries, mode_kind = mode_to_block[args.mode]

    # Use a requests session to manage HTTP connections and headers.
    with requests.Session() as session:
        session.headers.update({"User-Agent": "snpXplorer-performance-test/1.0"})
        if mode_kind == "annotation":
            results = benchmark_annotation_rsids(
                session=session,
                rsids=queries,
                repeats=args.repeats,
                timeout=args.timeout,
            )
        else:
            results = benchmark_haplotype_queries(
                session=session,
                queries=queries,
                browse_type=mode_kind,
                repeats=args.repeats,
                timeout=args.timeout,
            )
        run_and_report_block(block_name, results, args.output, n_queries=len(list(queries)), repeats=args.repeats)


if __name__ == "__main__":
    main()
