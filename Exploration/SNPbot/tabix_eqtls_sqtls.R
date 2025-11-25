# Script to create sqlite db for eqtls

# load necessary libraries
library(data.table)
library(stringr)

# set working directory -- eqtls
setwd('/Users/nicco/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/GitHub/snpXplorer/Data/databases/eQTLs')
tp = 'eqtls'
# set working directory -- sqtls
setwd('/Users/nicco/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/GitHub/snpXplorer/Data/databases/sQTLs')
tp = 'sqtls'

# list files in directory
files <- list.files()

# remove Ensemble_to_GeneName.txt from files
files <- files[files != 'Ensemble_to_GeneName.txt']

# read Ensemble_to_GeneName.txt
gene_map <- fread('Ensemble_to_GeneName.txt', header=TRUE, sep="\t")

# make sure the numbers are not in scientific notation
options(scipen=999)

# iterate over files
for (file in files) {
    print(paste0("Processing file: ", file))
    # read file
    dt <- fread(file, header=TRUE, sep="\t")
    # split variant_id into chrom, pos, ref, alt
    dt$chrom = str_replace_all(str_split_fixed(dt$variant_id, "_", 5)[,1], 'chr', '')
    dt$pos = as.numeric(str_split_fixed(dt$variant_id, "_", 5)[,2])
    dt$ref = str_split_fixed(dt$variant_id, "_", 5)[,3]
    dt$alt = str_split_fixed(dt$variant_id, "_", 5)[,4]
    # split gene_id by '.' and keep first part
    dt$gene_id2 = str_split_fixed(dt$gene_id, "\\.", 2)[,1]
    # merge with gene_map to get gene names
    dt = merge(dt, gene_map, by.x = "gene_id2", by.y = "Gene stable ID", all.x = TRUE)
    # select relevant columns
    sb = dt[, c('chrom', 'pos', 'ref', 'alt', 'gene_id', 'Gene name', 'tissue', 'tss_distance', 'maf', 'pval_nominal', 'slope')]
    # sort by position
    sb = sb[order(sb$pos), ]
    # write table
    if (tp == 'eqtls'){
        fwrite(sb, file = paste0('chr', sb$chrom[1], '_eQTL.txt'), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
        # bgzip
        system(paste0('bgzip chr', sb$chrom[1], '_eQTL.txt'))
        # tabix
        system(paste0('tabix -S 1 -s 1 -b 2 -e 2 chr', sb$chrom[1], '_eQTL.txt.gz'))
    } else if (tp == 'sqtls'){
        fwrite(sb, file = paste0('chr', sb$chrom[1], '_sQTL.txt'), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
        # bgzip
        system(paste0('bgzip chr', sb$chrom[1], '_sQTL.txt'))
        # tabix
        system(paste0('tabix -S 1 -s 1 -b 2 -e 2 chr', sb$chrom[1], '_sQTL.txt.gz'))
    }
}