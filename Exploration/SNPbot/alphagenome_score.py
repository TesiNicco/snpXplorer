from __future__ import annotations
import argparse
import datetime as dt
import json
import logging
import os
import pickle
import sys
import re
from pathlib import Path
import pandas as pd
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

LOGGER = logging.getLogger(__name__)

def _load_default_api_key() -> str | None:
    """Load the AlphaGenome API key from config.txt next to this script."""
    try:
        config_path = Path(__file__).resolve().with_name("config.txt")
        cfg = pd.read_table(config_path, sep=r"\s+", engine="python")
        api_key = str(cfg["alphagenome_api"].iloc[0]).strip()
        return api_key or None
    except Exception:
        return None

# Functions
# Helper function to convert DataFrames with object columns to a format that can be safely written to parquet.
def _to_parquet_safe(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of df with object columns converted to parquet-safe types."""
    safe_df = df.copy()
    for col in safe_df.columns:
        if safe_df[col].dtype != "object":
            continue
        # Keep pure string columns as-is. Convert mixed/custom Python objects.
        non_null = safe_df[col].dropna()
        if non_null.empty:
            continue
        if non_null.map(lambda v: isinstance(v, str)).all():
            continue
        safe_df[col] = safe_df[col].map(lambda v: None if pd.isna(v) else str(v))
    return safe_df

# Function to parse user variant input into (chrom, pos, ref, alt)
def _parse_variant_input(variant_str: str) -> tuple[str, int, str, str]:
    """Parse user variant input into (chrom, pos, ref, alt)."""
    value = variant_str.strip()
    patterns = [
        r"^(chr[\w]+|[\w]+):(\d+):([ACGTNacgtn]+)>([ACGTNacgtn]+)$",
        r"^(chr[\w]+|[\w]+):(\d+):([ACGTNacgtn]+):([ACGTNacgtn]+)$",
        r"^(chr[\w]+|[\w]+)-(\d+)-([ACGTNacgtn]+)-([ACGTNacgtn]+)$",
    ]
    match = None
    for pattern in patterns:
        match = re.fullmatch(pattern, value)
        if match is not None:
            break
    if match is None:
        raise ValueError(
            "Invalid variant format. Use one of: "
            "chr11:11111111:A>G, chr11:11111111:A:G, chr11-11111111-A-G"
        )
    chrom, pos_str, ref, alt = match.groups()
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    pos = int(pos_str)
    if pos < 1:
        raise ValueError(f"Position must be >= 1. Got: {pos}")
    ref = ref.upper()
    alt = alt.upper()
    if ref == alt:
        raise ValueError("REF and ALT cannot be the same.")
    return chrom, pos, ref, alt

# Function to build argument parser 
def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run AlphaGenome prediction-only scoring for one variant and save "
            "df_scores into an output directory."
        )
    )
    parser.add_argument(
        "variant",
        help=(
            "Variant string. Supported: chr11:11111111:A>G, "
            "chr11:11111111:A:G, chr11-11111111-A-G"
        ),
    )
    parser.add_argument(
        "--api-key",
        default=_load_default_api_key(),
        help="AlphaGenome API key (if different from default).",
    )
    parser.add_argument(
        "--sequence-length",
        default="1MB",
        choices=["2KB", "16KB", "100KB", "500KB", "1MB"],
        help="Sequence length around variant for prediction.",
    )
    parser.add_argument(
        "--organism",
        default="human",
        choices=["human"],
        help="Organism (currently human only).",
    )
    parser.add_argument(
        "--output-root",
        default="outputs",
        help="Root output directory.",
    )
    parser.add_argument(
        "--stdout",
        action="store_true",
        help="Write the main dataframe as a pickled Python object to stdout instead of creating output files.",
    )
    parser.add_argument(
        "--genes",
        action="store_true",
        help='Keep only entries with a non-null "gene_name".',
    )
    parser.add_argument(
        "--tissue",
        default=None,
        help=(
            "Filter rows to biosamples containing one or more text terms "
            "(case-insensitive). Comma-separated terms are supported; "
            'e.g. --tissue brain or --tissue blood,brain,liver.'
        ),
    )
    return parser

# Main function
def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S",)

    parser = _build_parser()
    args = parser.parse_args()

    if not args.api_key:
        raise ValueError("Missing API key. Pass --api-key or set ALPHAGENOME_API_KEY.")

    LOGGER.info("Parsing variant input")
    chrom, pos, ref_base, alt_base = _parse_variant_input(args.variant)

    masked_api_key = (
        f"{args.api_key[:6]}...{args.api_key[-4:]}"
        if len(args.api_key) > 10
        else "[provided]"
    )
    LOGGER.info("Input summary")
    LOGGER.info("- variant: %s:%s:%s>%s", chrom, pos, ref_base, alt_base)
    LOGGER.info("- organism: %s", args.organism)
    LOGGER.info("- sequence_length: %s", args.sequence_length)
    LOGGER.info("- genes_only filter: %s", args.genes)
    LOGGER.info("- tissue filter: %s", args.tissue if args.tissue else "None")
    LOGGER.info("- output_root: %s", args.output_root)
    LOGGER.info("- api_key: %s", masked_api_key)

    organism_map = {"human": dna_client.Organism.HOMO_SAPIENS}
    organism = organism_map[args.organism]

    sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[f"SEQUENCE_LENGTH_{args.sequence_length}"]

    LOGGER.info("Initializing AlphaGenome model client")
    dna_model = dna_client.create(args.api_key)

    variant = genome.Variant(chromosome=chrom, position=pos, reference_bases=ref_base, alternate_bases=alt_base, name=f"{chrom}:{pos}:{ref_base}>{alt_base}")

    interval = variant.reference_interval.resize(sequence_length)
    LOGGER.info("Scoring variant on interval: %s", interval)

    variant_scores = dna_model.score_variant(interval=interval, variant=variant, organism=organism, variant_scorers=list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values()))
    
    LOGGER.info("Converting model outputs to tidy dataframe")
    merged_scores = variant_scorers.tidy_scores(variant_scores)
    n_rows_before_filter = int(len(merged_scores))
    LOGGER.info("Rows before filters: %d", n_rows_before_filter)
    n_rows_after_genes_filter = None
    n_rows_after_tissue_filter = None

    if args.genes:
        if "gene_name" not in merged_scores.columns:
            raise ValueError('"gene_name" column not found in scores output.')
        merged_scores = merged_scores.dropna(subset=["gene_name"]).copy()
        n_rows_after_genes_filter = int(len(merged_scores))
        LOGGER.info("Applied genes filter. Rows after filter: %d", n_rows_after_genes_filter)

    if args.tissue:
        if "biosample_name" not in merged_scores.columns:
            raise ValueError('"biosample_name" column not found in scores output.')
        tissue_query = str(args.tissue)
        tissue_terms = [t.strip() for t in tissue_query.split(",") if t.strip()]
        if not tissue_terms:
            raise ValueError("--tissue cannot be empty.")

        tissue_mask = pd.Series(False, index=merged_scores.index)
        biosample_series = merged_scores["biosample_name"].astype(str)
        for term in tissue_terms:
            tissue_mask = tissue_mask | biosample_series.str.contains(term, case=False, na=False, regex=False)

        merged_scores = merged_scores[tissue_mask].copy()
        n_rows_after_tissue_filter = int(len(merged_scores))
        LOGGER.info("Applied tissue filter (%s). Rows after filter: %d", ",".join(tissue_terms), n_rows_after_tissue_filter)

    merged_scores["query_locus"] = f"{chrom}:{pos}"
    merged_scores["ref_base"] = ref_base
    merged_scores["alt_base"] = alt_base
    merged_scores["variant"] = str(variant)

    if args.stdout:
        LOGGER.info("Writing pickled dataframe to stdout")
        pickle.dump(merged_scores, sys.stdout.buffer, protocol=pickle.HIGHEST_PROTOCOL)
        return

    ts = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(args.output_root) / f"{chrom}_{pos}_{ref_base}_{alt_base}_{ts}"
    out_dir.mkdir(parents=True, exist_ok=True)
    LOGGER.info("Created output directory: %s", out_dir)

    csv_path = out_dir / "df_scores.csv"
    parquet_path = out_dir / "df_scores.parquet"
    LOGGER.info("Writing CSV output: %s", csv_path)
    merged_scores.to_csv(csv_path, index=False)

    parquet_written = False
    try:
        LOGGER.info("Writing Parquet output: %s", parquet_path)
        _to_parquet_safe(merged_scores).to_parquet(parquet_path, index=False)
        parquet_written = True
    except Exception as e:
        LOGGER.warning("Could not write parquet output: %s", e)

    summary_path = out_dir / "run_summary.json"
    LOGGER.info("Writing run summary: %s", summary_path)
    summary_path.write_text(
        json.dumps(
            {
                "variant": f"{chrom}:{pos}:{ref_base}>{alt_base}",
                "sequence_length": args.sequence_length,
                "organism": args.organism,
                "genes_only": bool(args.genes),
                "tissue_filter": args.tissue,
                "n_rows_before_filter": n_rows_before_filter,
                "n_rows_after_genes_filter": n_rows_after_genes_filter,
                "n_rows_after_tissue_filter": n_rows_after_tissue_filter,
                "n_rows": int(len(merged_scores)),
                "outputs": {
                    "csv": str(csv_path),
                    "parquet": str(parquet_path) if parquet_written else None,
                },
            },
            indent=2,
        )
        + "\n"
    )

    LOGGER.info("Done. Saved scores to: %s", out_dir)
    LOGGER.info("- %s", csv_path)
    if parquet_written:
        LOGGER.info("- %s", parquet_path)
    LOGGER.info("- %s", summary_path)

if __name__ == "__main__":
    main()
