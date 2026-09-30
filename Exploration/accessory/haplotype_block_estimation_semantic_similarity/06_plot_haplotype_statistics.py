"""Plot haplotype summaries and estimate genomic coverage."""

# Libraries
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import numpy as np
import pandas as pd

# Define ancestry
ancestry = 'African'

# Define input paths depending on ancestry
if ancestry == 'European':
    INPUT_PATH = Path("European_haplotypes/20251122_ld_clusters_low_recombination_intervals.tsv.gz")
    OUTPUT_DIR = Path("European_haplotypes")
elif ancestry == 'African':
    INPUT_PATH = Path("African_haplotypes/20260914_thr0.2_haplotypes.tsv.gz")
    OUTPUT_DIR = Path("African_haplotypes")

# Define plot preferences
PLOT_PREFIX = "20260917"
GENOME_SIZE_BP = 3_200_000_000

# Function to load haplotypes from the input file
def load_haplotypes(input_path=INPUT_PATH):
    """Load the compressed, tab-separated haplotype table."""
    df = pd.read_csv(
        input_path,
        sep="\t",
        compression="infer",
        dtype={"CHR": str},
    )
    required_columns = {"CHR", "BP1", "BP2", "KB", "NSNPS"}
    missing_columns = required_columns.difference(df.columns)
    if missing_columns:
        raise ValueError(
            f"Input table is missing required columns: {sorted(missing_columns)}"
        )
    return df

# Functions for chromosome sorting and ordering
def _chromosome_sort_key(chromosome):
    """Sort autosomes numerically, followed by X, Y, and mitochondrial DNA."""
    chromosome = str(chromosome)
    label = chromosome[3:] if chromosome.lower().startswith("chr") else chromosome
    label = label.upper()
    if label.isdigit():
        return 0, int(label)
    special_order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return 1, special_order.get(label, 10_000), label

# Functions for preparing output paths and plotting
def _chromosome_order(df):
    return sorted(df["CHR"].dropna().astype(str).unique(), key=_chromosome_sort_key)

# Function to prepare output paths for plots
def _prepare_output_path(output_dir, plot_name):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir / f"{PLOT_PREFIX}_{plot_name}.pdf"

# Function to plot the number of haplotypes per chromosome
def plot_haplotypes_per_chromosome(df, output_dir=OUTPUT_DIR):
    """Plot the number of haplotypes (rows) on each chromosome."""
    chromosome_order = _chromosome_order(df)
    counts = (
        df.assign(CHR=df["CHR"].astype(str))
        .groupby("CHR")
        .size()
        .reindex(chromosome_order)
    )
    positions = np.arange(len(chromosome_order))
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar(positions, counts.values, color="#4C78A8")
    ax.set(
        xlabel="Chromosome",
        ylabel="Number of haplotypes",
        title="Number of haplotypes per chromosome",
        xticks=positions,
        xticklabels=chromosome_order,
    )
    ax.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    output_path = _prepare_output_path(output_dir, "n_haplotypes_per_chr")
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig, ax

# Function to plot violin plots by chromosome
def _plot_violin_by_chromosome(df, value_column, ylabel, title, output_dir, plot_name):
    chromosome_order = _chromosome_order(df)
    plot_df = df.dropna(subset=["CHR", value_column]).copy()
    plot_df["CHR"] = plot_df["CHR"].astype(str)
    distributions = [plot_df.loc[plot_df["CHR"] == chromosome, value_column].to_numpy() for chromosome in chromosome_order]
    positions = np.arange(len(chromosome_order))
    fig, ax = plt.subplots(figsize=(12, 5))
    violins = ax.violinplot(
        distributions,
        positions=positions,
        widths=0.8,
        showmeans=False,
        showmedians=False,
        showextrema=False,
    )
    for body in violins["bodies"]:
        body.set_facecolor("#72B7B2")
        body.set_edgecolor("#2F4F4F")
        body.set_alpha(0.8)
    ax.boxplot(
        distributions,
        positions=positions,
        widths=0.16,
        patch_artist=True,
        showfliers=False,
        boxprops={"facecolor": "white", "edgecolor": "#333333"},
        medianprops={"color": "#E45756", "linewidth": 1.5},
        whiskerprops={"color": "#333333"},
        capprops={"color": "#333333"},
    )
    ax.set(
        xlabel="Chromosome",
        ylabel=ylabel,
        title=title,
        xticks=positions,
        xticklabels=chromosome_order,
    )
    ax.set_yscale("log")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    output_path = _prepare_output_path(output_dir, plot_name)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig, ax

# Function to plot the haplotype size distribution across chromosomes
def plot_haplotype_size_distribution(df, output_dir=OUTPUT_DIR):
    """Plot the haplotype-size (KB) distribution for each chromosome."""
    return _plot_violin_by_chromosome(
        df=df,
        value_column="KB",
        ylabel="Haplotype size (KB)",
        title="Haplotype size distribution across chromosomes",
        output_dir=output_dir,
        plot_name="haplotype_size_distribution",
    )

# Function to plot the number of SNPs per haplotype across chromosomes
def plot_snps_per_haplotype(df, output_dir=OUTPUT_DIR):
    """Plot the NSNPS distribution for each chromosome."""
    return _plot_violin_by_chromosome(
        df=df,
        value_column="NSNPS",
        ylabel="Number of SNPs per haplotype",
        title="SNP count per haplotype across chromosomes",
        output_dir=output_dir,
        plot_name="n_snps_per_haplotype",
    )

# Function to calculate the genomic coverage of haplotypes
def calculate_haplotype_coverage(df, genome_size_bp=GENOME_SIZE_BP):
    """
    Merge inclusive BP1-BP2 intervals and estimate their genomic coverage.

    Directly adjacent intervals are merged as part of the same covered region.
    """
    intervals = df[["CHR", "BP1", "BP2"]].dropna().copy()
    intervals["CHR"] = intervals["CHR"].astype(str)
    intervals["BP1"] = pd.to_numeric(intervals["BP1"], errors="raise").astype(np.int64)
    intervals["BP2"] = pd.to_numeric(intervals["BP2"], errors="raise").astype(np.int64)
    if (intervals["BP1"] > intervals["BP2"]).any():
        raise ValueError("At least one haplotype has BP1 greater than BP2")
    if genome_size_bp <= 0:
        raise ValueError("genome_size_bp must be positive")
    covered_bp = 0
    merged_interval_count = 0
    for _, chromosome_intervals in intervals.groupby("CHR", sort=False):
        chromosome_intervals = chromosome_intervals.sort_values(["BP1", "BP2"])
        current_start = None
        current_end = None
        for start, end in chromosome_intervals[["BP1", "BP2"]].itertuples(
            index=False, name=None
        ):
            if current_start is None:
                current_start, current_end = start, end
            elif start <= current_end + 1:
                current_end = max(current_end, end)
            else:
                covered_bp += current_end - current_start + 1
                merged_interval_count += 1
                current_start, current_end = start, end
        if current_start is not None:
            covered_bp += current_end - current_start + 1
            merged_interval_count += 1
    coverage_fraction = covered_bp / genome_size_bp
    return {
        "covered_bp": covered_bp,
        "merged_interval_count": merged_interval_count,
        "genome_size_bp": genome_size_bp,
        "coverage_fraction": coverage_fraction,
        "coverage_percent": coverage_fraction * 100,
    }

# Function to run all analyses and generate plots
def run_all(input_path=INPUT_PATH, output_dir=OUTPUT_DIR):
    """Load the data, create all plots, and print the coverage estimate."""
    df = load_haplotypes(input_path)
    df = df.drop_duplicates()
    figures = [
        plot_haplotypes_per_chromosome(df, output_dir),
        plot_haplotype_size_distribution(df, output_dir),
        plot_snps_per_haplotype(df, output_dir),
    ]
    for fig, _ in figures:
        plt.close(fig)
    coverage = calculate_haplotype_coverage(df)
    print(f"Haplotype rows: {len(df):,}")
    print(f"Merged intervals: {coverage['merged_interval_count']:,}")
    print(f"Covered bases: {coverage['covered_bp']:,}")
    print(
        "Approximate genome coverage: "
        f"{coverage['coverage_fraction']:.6f} "
        f"({coverage['coverage_percent']:.3f}%)"
    )
    return df, coverage

# Main execution block
if __name__ == "__main__":
    haplotypes, coverage = run_all()
