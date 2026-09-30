# Libraries
import re
import pandas as pd
import numpy as np
import os
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster
import gzip
import subprocess
from multiprocessing import Pool
from tqdm.auto import tqdm

# Functions
# Function to split recombination rates into intervals
def recomb_intervals(df, rate_col="Rate(cM/Mb)", chrom_col="Chromosome", pos_col="Position(bp)", threshold=25.0):
    """
    Continuous segmentation of the genome using recombination rate breakpoints.
    For each chromosome:
      - Start at the first Position(bp)
      - Whenever Rate(cM/Mb) > threshold at position P, close the current
        interval at P-1 and start a new one at P
      - Final interval ends at the last Position(bp)

    Output covers the full range [min(pos), max(pos)] for each chromosome
    with no gaps and no overlaps.
    """
    # Ensure sorted
    df = df.sort_values([chrom_col, pos_col]).reset_index(drop=True)
    intervals = []
    # Process each chromosome independently
    for chrom, sub in df.groupby(chrom_col, sort=False):
        sub = sub.sort_values(pos_col)
        positions = sub[pos_col].values
        rates = sub[rate_col].values
        if len(positions) == 0:
            continue
        start = int(positions[0])
        last_pos = int(positions[-1])
        # Loop over all positions on this chromosome
        for pos, rate in zip(positions, rates):
            pos = int(pos)
            if rate > threshold:
                # Close current interval at pos - 1 (if non-empty)
                end = pos - 1
                if end >= start:
                    intervals.append((chrom, start, end))
                # New interval starts at this breakpoint
                start = pos
        # Close last interval to the last observed position
        if start <= last_pos:
            intervals.append((chrom, start, last_pos))
    return pd.DataFrame(intervals, columns=["chrom", "start", "end"])

# Function to calculate LD given a region of interest
def calculateLD(chrom, start, end, genotype_data_path, min_snps=2):
    # define chromosome id
    chrom_id = chrom.replace('chr', '')
    # define input file of interest
    input_file = os.path.join(genotype_data_path, f'{chrom}_african')
    # plink2 command to calculate LD
    plink_command = f'plink2 --pfile {input_file} --chr {chrom_id} --from-bp {start} --to-bp {end} --r2-unphased square --maf 0.005 --out ld_{chrom}_{start}_{end}'
    # execute command
    subprocess.run(plink_command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    # read ld matrix in
    ld_mt = pd.read_csv(f'ld_{chrom}_{start}_{end}.unphased.vcor2', delim_whitespace=True, header=None)
    ld_mt_names = pd.read_csv(f'ld_{chrom}_{start}_{end}.unphased.vcor2.vars', delim_whitespace=True, header=None)
    # clean up temporary files
    os.remove(f'ld_{chrom}_{start}_{end}.unphased.vcor2.vars')
    os.remove(f'ld_{chrom}_{start}_{end}.unphased.vcor2')
    os.remove(f'ld_{chrom}_{start}_{end}.log')
    # set row and column names
    ld_mt.index = ld_mt_names[0]
    ld_mt.columns = ld_mt_names[0]
    # do clustering
    ld_labels = ld_clusters_complete(ld_mt.values, thresh=thresh)
    # get variant ids
    variant_ids = ld_mt.index.tolist()
    # summarize clusters
    cluster_df = cluster_summary(ld_labels, variant_ids)
    # restrict to clusters with at least min_snps
    cluster_df = cluster_df[cluster_df['NSNPS'] >= min_snps].reset_index(drop=True)
    return cluster_df

# Function to do basic clustering on the LD matrix
def ld_clusters_complete(LD, thresh):
    D = 1.0 - LD
    np.fill_diagonal(D, 0.0)
    Z = linkage(squareform(D, checks=False), method='complete')
    labels = fcluster(Z, t=1.0 - thresh, criterion='distance')
    return labels

# variant_id format expected: "CHR:POS:REF:ALT" (e.g., "17:130306:G:A")
VAR_RE = re.compile(r"^([^:]+):(\d+):[^:]+:[^:]+$")

# Parse variant ID into (chrom, pos)
def parse_variant(v):
    m = VAR_RE.match(v)
    if not m:
        raise ValueError(f"Bad variant id: {v}")
    chrom, pos = m.group(1), int(m.group(2))
    return chrom, pos

# Summarize clusters into a DataFrame
def cluster_summary(labels: np.ndarray, variant_ids: list) -> pd.DataFrame:
    if len(labels) != len(variant_ids):
        raise ValueError("labels and variant_ids must have same length")
    # parse chrom/pos once
    chroms, poss = zip(*[parse_variant(v) for v in variant_ids])
    df = pd.DataFrame({
        "cluster": labels,
        "CHR": chroms,
        "POS": poss,
        "SNP": variant_ids
    })
    # (optional) sanity: ensure a single chromosome per LD matrix
    if df["CHR"].nunique() > 1:
        # if you truly run multi-chrom LD (unusual), this still works per cluster
        pass
    # group by cluster → summarize
    g = df.groupby("cluster", sort=False)
    out = g.apply(lambda sub: pd.Series({
        "CHR": sub["CHR"].iloc[0],
        "BP1": int(sub["POS"].min()),
        "BP2": int(sub["POS"].max()),
        "KB":  round((sub["POS"].max() - sub["POS"].min()) / 1000.0, 2),
        "NSNPS": int(len(sub)),
        "SNPS": "|".join(sub.sort_values("POS")["SNP"].tolist())
    })).reset_index(drop=True)
    # sort by CHR then BP1 (or just BP1 if single-chr)
    out = out.sort_values(["CHR", "BP1"]).reset_index(drop=True)
    return out

# Function to process interval (for parallel processing)
def process_interval(interval):
    chrom = interval['chrom']
    start = interval['start']
    end = interval['end']
    try:
        return calculateLD(chrom, start, end, genotype_data_path, min_snps=5)
    except:
        return None

# Data paths
recomb_rates_path = '/project/holstegelab/Share/nicco/reference_files/recombination_rates_hg38.txt.gz'
genotype_data_path = '/project/holstegelab/Share/gwas_array/1000Genome_common_snps/additional_samples/african_samples/'
thresh = 0.2

# Load recombination rate data
recomb_rates = pd.read_csv(recomb_rates_path, sep='\t')
# Split recombination rates into low recombination intervals
recomb_intervals_df = recomb_intervals(recomb_rates, threshold=25)
recomb_intervals_df['size'] = recomb_intervals_df['end'] - recomb_intervals_df['start'] + 1

# Convert dataframe to a list of dicts for processing
intervals = recomb_intervals_df.to_dict(orient='records')

# Run in parallel
all_ld_results = []
n_procs = 8
chunk_size = 8
with Pool(processes=n_procs) as pool:
    for res in tqdm(pool.imap_unordered(process_interval, intervals, chunksize=chunk_size), total=len(intervals), desc="Intervals"):
        if res is not None:
            all_ld_results.append(res)

# Combine results
final_df = pd.concat(all_ld_results, ignore_index=True) if all_ld_results else pd.DataFrame()

# Sort by CHR and BP1
final_df = final_df.sort_values(['CHR', 'BP1']).reset_index(drop=True)

# Write to the sensitivity analysis file
oupath = '/project/holstegelab/Share/nicco/haplotypeBlocks/R2cluster/African_haplotypes/20260914_thr0.2_haplotypes.tsv'
final_df.to_csv(oupath, sep='\t', index=False)

# Write final results to file
#final_df.to_csv('20251122_ld_clusters_low_recombination_intervals.tsv', sep='\t', index=False)
# BGzip the output file
#os.system('bgzip 20251122_ld_clusters_low_recombination_intervals.tsv')
# Tabix
#os.system('tabix -S 1 -s 1 -b 2 -e 2 20251122_ld_clusters_low_recombination_intervals.tsv.gz')