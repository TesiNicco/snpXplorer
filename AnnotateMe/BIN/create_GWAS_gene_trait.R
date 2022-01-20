# Script to parse the GWAS catalog file upon new updates to create the GWAS-gene-trait dataset #
################################################################################################

# Libraries
library(data.table)
library(stringr)
library(parallel)

# Function
# function to extract GWAS-gene-trait associations for a given chromosome
extractInfo = function(chrom, data){
    # subset to chromosome of interest
    sb = data[which(data$CHR_ID == chrom),]
    # create and empty dataframe with the same column, but empty
    new = data.frame(gene = as.character(), trait = as.character(), study = as.character(), jounal = as.character())
    # main loop on lines
    for (i in 1:nrow(sb)){
        # extract all genes for each row -- some row can have multiple genes
        all_genes = str_split(data$MAPPED_GENE[i], " - ")[[1]]
        # put results in a temporary df
        tmp_df = data.frame(gene = all_genes, trait = rep(data$MAPPED_TRAIT[i], length(all_genes)), study = rep(data$"FIRST AUTHOR"[i], length(all_genes)), journal = rep(data$JOURNAL[i], length(all_genes)))
        # attache to main dataframe
        new = rbind(new, tmp_df)
    }
    return(new)
}

# Main
# 1. When a new update of the GWAS catalog is given and copied into snpXplorer/AnnotateMe/INPUTS_OTHER/, we need to briefly clean it to remove duplicates
    # read original file for comparison
    original = fread("GWAS_catalog_20200310.txt.gz", h=T, stringsAsFactors=F, quote="")
    # then read the new one
    new = fread("GWAS_catalog_20211210.txt.gz", h=T, stringsAsFactors=F, quote="")
    # compare dimensions -- make sure the column are the same
    dim(original)
    dim(new)
    colnames(original) == colnames(new)
    # exclude duplicated rows from new file
    new = new[!duplicated(new),]
    # finally write the file
    write.table(new, "GWAS_catalog_20211210.txt", row.names=F, sep="\t", quote=F)
    # and gzip it
    system("gzip GWAS_catalog_20211210.txt")

# 2. After modifying the main file, we should look at the GWAS-gene-trait dataset and update that as well
    # read original file for comparison
    original = fread("Gwas_catalog_Gene_Traits.txt", h=T, stringsAsFactors=F, quote="")
    # read the updated gwas catalog associations
    updated_gwas = fread("GWAS_catalog_20211210.txt.gz", h=T, stringsAsFactors=F, quote="")
    # in multiprocessing, run the function to create the GWAS-gene-trait dataset
    results = rbindlist(mclapply(X = seq(1, 22), FUN = extractInfo, data = updated_gwas, mc.cores = 2))
    # remove duplicates if they are there
    results = results[!duplicated(results),]
    # write output
    write.table(results, "20211210_Gwas_catalog_Gene_Traits.txt", quote=F, row.names=F, sep = "\t")
    # Done!
