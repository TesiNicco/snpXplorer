# Script to parse eqtl and sqtl data

# libraries
import os
import pandas as pd
import gzip
import pickle

# functions
# function to read and store all sqtls
def readSqtls(all_tissues):
    all_sqtls = {}              # output dictionary where the keys are the chromosomes
    for f in all_tissues:       # loop on all tissues
        f = f.rstrip()
        tissue = f.split('/')[-1].split('.')[0]                            # get tissue
        print("## Working on tissue --> %s" %(tissue))
        with gzip.open(f, 'r') as fopen:    # read data
            for line in fopen:
                if line.startswith(b'variant_id'):
                    pass
                else:
                    line = line.rstrip().split(b'\t')
                    snp_id, info, effect, p = line[0], line[1], line[7], float(line[-1])
                    if p <= 0.05:                                               # check if it's a significant sqtl
                        allele1, allele2 = snp_id.split(b'_')[2:4]          # also extract alleles
                        if len(allele1) + len(allele2) == 2:
                            gene = info.split(b':')[-1]
                            chrom = snp_id.split(b'_')[0]
                            if chrom in all_sqtls.keys():
                                all_sqtls[chrom].append([snp_id, gene, effect, p, tissue])
                            else:
                                all_sqtls[chrom] = [[snp_id, gene, effect, p, tissue]]
    return(all_sqtls)

# function to save obj to file
def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

# function to load obj from file
def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

# main
# main argument is the folder where the single-tissue data is stored
path_sqtl = '/home/nfs/ntesi/bulkALL/niccolo/GTEx_Analysis_v8_sQTL'

# 1. read all files
all_tissues = list(os.popen("ls %s/*sqtl*" %(path_sqtl)))

# 2. main loop through files to store them -- use function
all_sqtls = readSqtls(all_tissues)
save_obj(all_sqtls, "%s/summary_sqtls" %(path_sqtl))
all_sqtls = load_obj("%s/summary_sqtls" %(path_sqtl))

# 3. finally write data
for chrom in all_sqtls.keys():
    print("## Working on chromosome --> %s" %(chrom.decode("utf-8")))
    tmp = all_sqtls[chrom]          # extract data for each chromosome
    df = pd.DataFrame.from_dict(tmp)    # convert to dataframe
    df = df.rename(columns={0: "snp_id", 1: "gene", 2: "effect", 3: "pvalue", 4: "tissue"})
    df = df.sort_values(by=['snp_id'])
    str_df = df.applymap(lambda x: x.decode() if isinstance(x, bytes) else x)
    str_df.to_csv('%s/%s_summary_sqtls.txt.gz' %(path_sqtl, chrom.decode("utf-8")), sep = '\t', header = False, index = False)


    

