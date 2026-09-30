# Script to generate a database per ancestry of LD

# Libraries
import os
import sqlite3
import pandas as pd

# Define ancestry
ancestry = 'African'

# Define output path for the database
output_folder = f"./LD_database/{ancestry}"
os.makedirs(output_folder, exist_ok=True)

# Define the path to the genotyping data
if ancestry == 'African':
    geno_data_path = f"/project/holstegelab/Share/gwas_array/1000Genome_common_snps/additional_samples/african_samples"
elif ancestry == 'European':
    geno_data_path = f"/project/holstegelab/Share/gwas_array/TOPMED/data/2025_allBatches_combined/snp_qc_imputation/Imputed_unscRambled_genotypes"

# If ancestry is european, need to make a phenotype file
if ancestry == 'European':
    pheno_data_path = f"/project/holstegelab/Share/gwas_array/TOPMED/phenotypes/Phenotypes_array_20250530.txt"
    pheno_data = pd.read_csv(pheno_data_path, sep="\t")
    pheno_data_eur = pheno_data[pheno_data['PC1'].notna()]
    pheno_data_eur_healthy = pheno_data_eur[pheno_data_eur['Diagnosis'].isin(['Control_LASA', 'Centenarian', 'SCD', 'Control_100plus', 'Control_other_twin', 'Control_path'])]
    pheno_data_eur_healthy['ID_GWAS'].to_csv(f"./LD_database/{ancestry}/pheno_data_eur_healthy.txt", sep="\t", index=False, header=False)

# Functions
# Create 2-Mb intervals every 1 Mb, clipped at the chromosome boundary.
def create_intervals(boundary, step=1000000):
    intervals = []
    for start in range(1, boundary + 1, step):
        end = min(start + 2 * step - 1, boundary)
        intervals.append((start, end))
    return intervals

# Save the already-filtered LD rows using the TOPMed database schema.
def save_ld_to_database(tmp_ld, db_path, batch_size=50000):
    connection = sqlite3.connect(db_path)
    try:
        connection.executescript("""
            CREATE TABLE IF NOT EXISTS variants (
                id INTEGER PRIMARY KEY,
                uniq TEXT NOT NULL UNIQUE,
                pos INTEGER NOT NULL
            );
            CREATE TABLE IF NOT EXISTS ld (
                v1 INTEGER NOT NULL,
                v2 INTEGER NOT NULL,
                r2_milli INTEGER NOT NULL,
                dist_bp INTEGER NOT NULL,
                PRIMARY KEY (v1, v2)
            ) WITHOUT ROWID;
            CREATE INDEX IF NOT EXISTS idx_variants_pos ON variants(pos);
            CREATE INDEX IF NOT EXISTS idx_ld_v1 ON ld(v1);
            CREATE INDEX IF NOT EXISTS idx_ld_v2 ON ld(v2);
        """)
        for offset in range(0, len(tmp_ld), batch_size):
            batch = tmp_ld.iloc[offset:offset + batch_size]
            variants = {}
            pairs = []
            for pa, a, pb, b, r2 in batch[
                ['POS_A', 'ID_A', 'POS_B', 'ID_B', 'UNPHASED_R2']
            ].itertuples(index=False, name=None):
                a = str(a).split(':', 1)[1]
                b = str(b).split(':', 1)[1]
                variants[a] = int(pa)
                variants[b] = int(pb)
                pairs.append((a, b, int(round(float(r2) * 1000)), abs(int(pb) - int(pa))))
            if not pairs:
                continue
            with connection:
                connection.executemany(
                    'INSERT OR IGNORE INTO variants (uniq, pos) VALUES (?, ?)',
                    variants.items())
                ids = {}
                keys = list(variants)
                for first in range(0, len(keys), 500):
                    subset = keys[first:first + 500]
                    placeholders = ','.join('?' for _ in subset)
                    ids.update(connection.execute(
                        f'SELECT uniq, id FROM variants WHERE uniq IN ({placeholders})',
                        subset).fetchall())
                connection.executemany(
                    'INSERT OR IGNORE INTO ld (v1, v2, r2_milli, dist_bp) VALUES (?, ?, ?, ?)',
                    ((min(ids[a], ids[b]), max(ids[a], ids[b]), r2, distance)
                     for a, b, r2, distance in pairs))
    finally:
        connection.close()

# Iterate over chromosomes 1 to 22
for chr in range(20, 23):
    chr = str(chr)
    window = 1000000
    db_path = os.path.join(output_folder, f"ld_chr{chr}.sqlite")
    if ancestry == 'African':
        cmd_boundary = f"tail -n 1 {geno_data_path}/chr{chr}_african.pvar"
    elif ancestry == 'European':
        cmd_boundary = f"tail -n 1 {geno_data_path}/chr{chr}.dose.unscrambled.pvar"
    boundary = int(os.popen(cmd_boundary).read().strip().split()[1])
    intervals = create_intervals(boundary)
    for start, end in intervals:
        print(f"[{intervals.index((start, end))+1}/{len(intervals)}] Processing chr{chr} interval: {start}-{end}")
        if ancestry == 'African':
            cmd_plink2 = f"plink2 --pfile {geno_data_path}/chr{chr}_african --chr {chr} --from-bp {start} --to-bp {end} --maf 0.005 --r2-unphased --ld-window-kb 1000 --ld-window 99999999 --ld-window-r2 0.2 --out {output_folder}/chr{chr}_{start}_{end} > /dev/null 2>&1"
        elif ancestry == 'European':
            cmd_plink2 = f"plink2 --pfile {geno_data_path}/chr{chr}.dose.unscrambled --keep ./LD_database/{ancestry}/pheno_data_eur_healthy.txt --chr {chr} --from-bp {start} --to-bp {end} --maf 0.005 --r2-unphased --ld-window-kb 1000 --ld-window 99999999 --ld-window-r2 0.2 --out {output_folder}/chr{chr}_{start}_{end} > /dev/null 2>&1"
        os.system(cmd_plink2)
        if not os.path.exists(f"{output_folder}/chr{chr}_{start}_{end}.vcor"):
            os.remove(f"{output_folder}/chr{chr}_{start}_{end}.log")
            continue
        tmp_ld = pd.read_csv(f"{output_folder}/chr{chr}_{start}_{end}.vcor", delim_whitespace=True)
        core_end = min(start + window - 1, boundary)
        tmp_ld = tmp_ld[tmp_ld['POS_A'].between(start, core_end)]
        save_ld_to_database(tmp_ld, db_path)
        os.remove(f"{output_folder}/chr{chr}_{start}_{end}.vcor")
        os.remove(f"{output_folder}/chr{chr}_{start}_{end}.log")
    if os.path.exists(db_path):
        connection = sqlite3.connect(db_path)
        try:
            connection.execute("ANALYZE")
            connection.commit()
        finally:
            connection.close()
