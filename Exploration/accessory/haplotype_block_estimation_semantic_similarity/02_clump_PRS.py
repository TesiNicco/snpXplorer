# Script to do clumping of All GWAS - This is used for making the PRS-ready variant sets
# We should now make this a little different and ancestry-aware
# We can use GWASum.py which was designed to do exactly this

# Libraries
import pandas as pd
import sqlite3
import os
import sys
from pathlib import Path
from typing import Optional
from functools import lru_cache
from liftover import get_lifter

# Executables
gwasum = 'python /project/holstegelab/Share/nicco/GWASum/GWASum.py'
threads = 6

# Input directory
ancestry = "European"
if ancestry == "African":
    sumstats_dir = Path(f"/project/holstegelab/Share/nicco/haplotypeBlocks/R2cluster/OpenGWAS_sumstats/African")
    output_dir = Path(f"/project/holstegelab/Share/nicco/haplotypeBlocks/R2cluster/OpenGWAS_sumstats/African/GWASum")
    genotype_data = Path(f"/project/holstegelab/Share/gwas_array/1000Genome_common_snps/additional_samples/african_samples")
elif ancestry == "European":
    sumstats_dir = Path(f"/project/holstegelab/Share/nicco/haplotypeBlocks/R2cluster/OpenGWAS_sumstats/European")
    output_dir = Path(f"/project/holstegelab/Share/nicco/haplotypeBlocks/R2cluster/OpenGWAS_sumstats/European/GWASum")
    genotype_data = Path(f"/project/holstegelab/Share/gwas_array/TOPMED/data/2025_allBatches_combined/snp_qc_imputation/Imputed_unscRambled_genotypes")

# Create output directory if it doesn't exist
output_dir.mkdir(parents=True, exist_ok=True)

# List summary stats
if ancestry == 'European':
    file_list = [f for f in os.listdir(sumstats_dir) if f.endswith('.txt.gz')]
    file_list = [x for x in file_list if not x.startswith('All')]
    file_list = [x for x in file_list if not x.startswith('eqtl')]
else:
    file_list = [f for f in os.listdir(sumstats_dir) if f.endswith('.txt')]

# Iterate over each summary stats
counter = 0
for file_name in file_list:
    counter += 1
    print(f"Processing file {counter}/{len(file_list)}: {file_name}")
    # Define full path to the summary stats file
    full_path = sumstats_dir / file_name
    # Define output file path
    if ancestry == 'European':
        output_file = output_dir / file_name.replace('.txt.gz', '')
    else:
        output_file = output_dir / file_name.replace('.txt', '')        
    if output_file.exists():
        continue
    # Compose the GWASum.py command
    if ancestry == "African":
        gwasum_command = f"{gwasum} --input {full_path} --chrom chr --pos position --effect ea --other nea --beta beta --se se --p p --geno  {genotype_data} --geno-pattern chr{{chr}}_african --output {output_file} --frq eaf --n n --clumping-p 5e-5 --clumping-r2 0.05 --clumping-kb 500 --cojo-p 5e-5 --cojo-wind 1000 --threads {threads} --liftover"
    elif ancestry == "European":
        gwasum_command = f"{gwasum} --input {full_path} --chrom chr --pos position --effect ea --other nea --beta beta --se se --p p --geno  {genotype_data} --output {output_file} --frq eaf --n n --clumping-p 5e-5 --clumping-r2 0.05 --clumping-kb 500 --cojo-p 5e-5 --cojo-wind 1000 --threads {threads} --liftover"
    os.system(gwasum_command + " > /dev/null 2>&1")
 
# Then get the list of files again
file_list_command = f"find {output_dir} -name 'combined_conditional_results.txt'"
file_list = [x for x in os.popen(file_list_command).read().split('\n') if x != ""]
len(file_list)

# Define output directory for reformatted files
reformat_dir = output_dir / "reformatted"
reformat_dir.mkdir(parents=True, exist_ok=True)

# Iterate over each clumped file and reformat
counter = 0
for file in file_list:
    counter += 1
    print(f"[{counter}/{len(file_list)}] Reformatting files")
    # Get trait id, which is the parent-parent-parent directory name
    trait_id = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(file))))
    # Read clumped file
    tmp = pd.read_csv(file, sep='\t')
    # We miss the alternative allele, so we need to get that info
    conditional_files = os.path.dirname(file).replace('conditional_results', 'conditional_files')
    get_all_files = f"cat {conditional_files}/*"
    all_files_content = [x for x in os.popen(get_all_files).read().split('\n') if x != ""]
    # Convert to dataframe
    all_files_df = pd.DataFrame([x.split('\t') for x in all_files_content[1:]], columns=all_files_content[0].split('\t'))
    all_files_df = all_files_df[all_files_df['SNP'] != "SNP"]
    # Merge with all_files_df -- add only A1 and A2 alleles
    tmp = tmp.merge(all_files_df[['SNP', 'A1', 'A2']], left_on='SNP', right_on='SNP', how='left')
    # Assign the OTHER_ALLELE as the allele among A1, A2 that is not the refA
    tmp['OTHER_ALLELE'] = tmp.apply(lambda row: row['A2'] if row['A1'] == row['refA'] else row['A1'], axis=1)
    tmp = tmp.drop(columns=['A1', 'A2'])
    # Select relevant columns
    tmp = tmp[['SNP', 'Chr', 'bp', 'refA', 'OTHER_ALLELE', 'b', 'se', 'p' , 'freq_geno']]
    # Rename columns
    tmp = tmp.rename(columns={'bp': 'POS', 'Chr': 'CHROM', 'refA': 'EFFECT_ALLELE', 'b': 'BETA', 'se': 'SE', 'p': 'P', 'freq_geno': 'EFFECT_ALLELE_FREQUENCY'})
    # Define output file
    output_file = reformat_dir / f"{trait_id}_clumped_reformatted.txt"
    # Write reformatted file
    tmp.to_csv(output_file, sep='\t', index=False)
