# Script to count the number of hits of independent 

# Libraries
import sys
import os
import pandas as pd
import numpy as np
from pyliftover import LiftOver
import pickle

# Define ancestry as global parameter
ancestry = 'African'

# Data paths
if ancestry == 'European':
    clusters = 'European_haplotypes/20251122_ld_clusters_low_recombination_intervals.tsv.gz'
    gwas = 'OpenGWAS_sumstats/European'
    ai_based = '20260714_trait_clusters_summary_AI_thresholds.csv'
elif ancestry == 'African':
    clusters = 'African_haplotypes/20260914_thr0.2_haplotypes.tsv'
    gwas = 'OpenGWAS_sumstats/African'
    ai_based = '20260916_trait_clusters_summary_AI_thresholds_African.csv'
threshold = 0.5

# Read clusters
clusters_df = pd.read_csv(clusters, sep='\t')

# Read independent GWAS loci
ai_based_df = pd.read_csv(ai_based, sep=',')
ai_based_df_filtered = ai_based_df[ai_based_df['Threshold'] == threshold]
indep_loci = ai_based_df_filtered['Final Representative file'].tolist()
len(indep_loci)

# Then read the relative GWASs
for gw in indep_loci:
    print(f'Processing GWAS: {indep_loci.index(gw)+1} out of {len(indep_loci)}       ', end='\r')
    # get gwas of interest
    gw_file = os.path.join(gwas, f'{gw}')
    # define tools and output file based on ancestry
    toolcat = 'zcat' if ancestry == 'European' else 'cat'
    toolgrep = 'zgrep' if ancestry == 'European' else 'grep'
    output_file = '20260714_All_indep_gwas_sumstats_AI.txt' if ancestry == 'European' else '20260917_All_indep_gwas_sumstats_AI_African.txt'
    # check if it's the first
    if indep_loci.index(gw) == 0:
        cmd = f'{toolcat} {gw_file} > {output_file}'
        os.system(cmd)
    else:
        cmd = f'{toolgrep} -v "position" {gw_file} >> {output_file}'
        os.system(cmd)

# Read combined GWAS summary statistics
combined_gwas_df = pd.read_csv(f'{output_file}', sep='\t')

# Liftover combined gwas to hg38
lo = LiftOver('hg19', 'hg38')
# Function to do liftover
def liftover_position(row):
    chrom = 'chr' + str(row['chr'])
    pos = int(row['position'])
    lifted = lo.convert_coordinate(chrom, pos)
    if lifted:
        return lifted[0][1]
    else:
        return np.nan
# Apply liftover
combined_gwas_df['position_hg38'] = combined_gwas_df.apply(liftover_position, axis=1)
# Convert position to integer, and remove NAs
combined_gwas_df['position_hg38'] = combined_gwas_df['position_hg38'].astype('Int64')
combined_gwas_df = combined_gwas_df.dropna(subset=['position_hg38'])

# We need to index 2 files -- 1 by trait and 1 by chr/position
# Sort liftovered gwas by trait id and position
combined_gwas_df = combined_gwas_df.sort_values(['id', 'position_hg38'], kind='mergesort').reset_index(drop=True)
output_bytrait = '20260714_All_indep_gwas_sumstats_AI_hg38_byTrait.txt' if ancestry == 'European' else '20260917_All_indep_gwas_sumstats_AI_hg38_byTrait_African.txt'
# Save liftovered gwas
combined_gwas_df.to_csv(f'{output_bytrait}', sep='\t', index=False)
# Compress with bgzip
os.system(f'~/.conda/envs/newtreat/bin/bgzip {output_bytrait}')
# Tabix
os.system(f'~/.conda/envs/newtreat/bin/tabix -S 1 -s 1 -b 14 -e 14 {output_bytrait}.gz')

# Then sort liftovered gwas by chr and position
combined_gwas_df = combined_gwas_df.sort_values(['chr', 'position_hg38'], kind='mergesort').reset_index(drop=True)  # stable sort
output_bychr = '20260714_All_indep_gwas_sumstats_AI_hg38.txt' if ancestry == 'European' else '20260917_All_indep_gwas_sumstats_AI_hg38_African.txt'
# Save liftovered gwas
combined_gwas_df.to_csv(f'{output_bychr}', sep='\t', index=False)
# Compress with bgzip
os.system(f'~/.conda/envs/newtreat/bin/bgzip {output_bychr}')
# Tabix
os.system(f'~/.conda/envs/newtreat/bin/tabix -S 1 -s 3 -b 14 -e 14 {output_bychr}.gz')
# Remove previous uncompressed file
#os.system('rm 20260714_All_indep_gwas_sumstats_AI.txt')

# Read combined GWAS summary statistics
combined_gwas_df = pd.read_csv(f'{output_bychr}.gz', sep='\t')

# Normalize types
clusters = clusters_df.copy()
gwas = combined_gwas_df.copy()
clusters['CHR'] = clusters['CHR'].astype(str)
gwas['chr'] = gwas['chr'].astype(str)
clusters[['BP1','BP2']] = clusters[['BP1','BP2']].astype(np.int64)

# Prepare outputs
clusters['assoc_count'] = 0
clusters['rsids'] = [[] for _ in range(len(clusters))]
clusters['traits'] = [[] for _ in range(len(clusters))]

# Pre-sorted per-chr arrays
by_chr = {}
for c, g in gwas.groupby('chr', sort=False):
    g2 = g.sort_values('position_hg38', kind='mergesort')  # stable
    by_chr[c] = {'pos': g2['position_hg38'].to_numpy(copy=False), 'rsid': g2['rsid'].astype(str).to_numpy(copy=False), 'trait': g2['trait'].astype(str).to_numpy(copy=False)}

# Alternative to match snps by exact positions in SNPS column -- likely do again from here
chrs = clusters['CHR'].unique().tolist()
for i, c in enumerate(chrs, 1):
    if c not in by_chr:
        continue
    # positions/rsids/traits from dictionary
    pos = by_chr[c]['pos']
    rs  = by_chr[c]['rsid']
    tr  = by_chr[c]['trait']
    # rows belonging to this chromosome
    idx = clusters.index[clusters['CHR'] == c]
    # ---- NEW: match SNPs by exact positions in SNPS column ----
    for row_i in idx:
        print(f"\rProcessing row {row_i}/{len(idx)} on chr {c}      ", end='\r')
        # extract SNP position list from SNPS column
        snp_list = clusters.at[row_i, 'SNPS'].split("|")
        snp_positions = [int(s.split(":")[1]) for s in snp_list]
        # find which dictionary positions match
        mask = np.isin(pos, snp_positions)
        # number of matching SNPs
        clusters.at[row_i, "assoc_count"] = mask.sum()
        if mask.any():
            clusters.at[row_i, "rsids"] = list(rs[mask])
            clusters.at[row_i, "traits"] = list(tr[mask])
        else:
            clusters.at[row_i, "rsids"] = []
            clusters.at[row_i, "traits"] = []

# Example: most associated intervals with their variants
top = clusters.sort_values('assoc_count', ascending=False).head(20)
tail = clusters.sort_values('assoc_count', ascending=False).tail(10)

# column rsids and traits to semicolon-separated strings
bck = clusters.copy()
clusters['rsids'] = clusters['rsids'].apply(lambda x: ';'.join(x) if x else '')
clusters['traits'] = clusters['traits'].apply(lambda x: '___'.join(x) if x else '')

# Save the backup as a pickle file
output_pkl = '20260714_ld_clusters_with_gwas_backup_latest.pkl' if ancestry == 'European' else '20260917_ld_clusters_with_gwas_backup_latest_African.pkl'
with open(output_pkl, 'wb') as f:
    pickle.dump(bck, f)
# Sort clusters by CHR and BP1
clusters = clusters.sort_values(['CHR', 'BP1'], kind='mergesort')
# Save final output
output_tsv = '20260714_ld_clusters_with_gwas_AI_latest.tsv' if ancestry == 'European' else '20260917_ld_clusters_with_gwas_AI_latest_African.tsv'
clusters.to_csv(output_tsv, sep='\t', index=False)
# Compress with bgzip
os.system(f'~/.conda/envs/newtreat/bin/bgzip {output_tsv}')
os.system(f'~/.conda/envs/newtreat/bin/tabix -S 1 -s 1 -b 2 -e 2 {output_tsv}.gz')
