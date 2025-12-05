# libraries
library(stringr)

# basic paths
MAIN = "/Annotation/"
MAIN_SNP = "/Annotation/RUNS/"
args = commandArgs(trailingOnly=TRUE)

## match the sqtls for each chromosome -- for multiprocessing purpose
matchSQTL <- function(i, mapping, ensembl, interesting_tissues, random_num){
    cat(paste0("## Doing sQTL annotation of chr", i, "\n"))
    tmp = mapping[which(mapping$chr == i), ]
    # write temporary output with eqtl ids to match
    write.table(tmp$code_gtex, paste0(MAIN, "INPUTS_OTHER/summary_sqtls/tmp_chr", i, "_", random_num, ".txt"), quote=F, row.names=F, col.names=F)
    # command to grep
    cmd = paste0("zgrep -F -f ", MAIN, "INPUTS_OTHER/summary_sqtls/tmp_chr", i, "_", random_num, ".txt ", MAIN, "INPUTS_OTHER/summary_sqtls/chr", i, "_summary_sqtls.txt.gz")
    tmp_res <- tryCatch({ system(cmd, intern = T) }, error=function(cond){ return(NA) }, warning=function(cond){ return(NA) })    
    # remove temporary data
    system(paste0("rm ", MAIN, "INPUTS_OTHER/summary_sqtls/tmp_chr", i, "_", random_num, ".txt"))
    # create variables for assignment
    tmp$sqtl = NA
    tmp$sqtl_tissue = NA
    # main loop to assign eqtls to snps
    if (!is.na(tmp_res[1])){
        # modify structure of data
        tmp_res = data.frame(stringr::str_split_fixed(tmp_res, "\t", 5), stringsAsFactors = F)
        # restrict to tissues of interest
        if (!('All_tissues' %in% interesting_tissues)){ tmp_res = tmp_res[which(tmp_res$X5 %in% str_split(interesting_tissues, ',')[[1]]), ] }
        for (snp in unique(tmp_res$X1)){
            snp_id = paste0(stringr::str_replace_all(paste(stringr::str_split_fixed(snp, "_", 5)[, 1:2], collapse = "_"), "chr", ""), "_")
            tissues = tmp_res$X5[which(tmp_res$X1 == snp)]
            genes = str_split_fixed(tmp_res$X2[which(tmp_res$X1 == snp)], '\\.', 2)[, 1]
            tissues_genes_df = data.frame(tissue = tissues, gene = genes)
            tissues_genes_df = merge(tissues_genes_df, ensembl, by.x = "gene", by.y = "Gene stable ID", all.x = T)
            # take union of all tissues and genes
            all_tissues = paste0(tissues_genes_df$tissue, collapse = ",")
            all_genes = paste0(tissues_genes_df$"Gene name", collapse = ",")
            # finally assign to relative snp
            tmp$sqtl[which(tmp$code_gtex == snp_id)] = all_genes
            tmp$sqtl_tissue[which(tmp$code_gtex == snp_id)] = all_tissues
        }
    }
    return(tmp)
}

# read arguments
snps_info_path = args[1]
interesting_tissues = args[2]
random_num = args[3]
# run function to read snps info
load(snps_info_path)
# read mapping between ensembl and gene name
ensembl = data.table::fread(paste0(MAIN, "INPUTS_OTHER/Ensemble_to_GeneName.txt"), h=T, sep = "\t")
# we can use the codes generated for the eqtl analysis to grep here -- no need to liftover again
all_chroms = unique(out.gtex$chr)
all_res_sqtl = data.table::rbindlist(parallel::mclapply(X = all_chroms, FUN = matchSQTL, mapping = out.gtex, ensembl = ensembl, interesting_tissues = interesting_tissues, random_num = random_num, mc.cores = 1))

# replace empty cells with NA -- this has to be different if the number of rows is 1 only
if (nrow(all_res_sqtl) == 1){
  all_res_sqtl = as.data.frame(t(apply(all_res_sqtl, 2, function(x) gsub("^$|^ $", NA, x))))
} else {
  all_res_sqtl = as.data.frame(apply(all_res_sqtl, 2, function(x) gsub("^$|^ $", NA, x)))
}
# we can also use LDlink -- this will give results for the snps of interest + all snps in LD with that
#out.annot.sb = out.annot[grep("rs", out.annot$ID),]
#miss = out.annot[which(!(out.annot$locus %in% out.annot.sb$locus)),]
#gtex_info = LDexpress(out.annot.sb$ID, pop = "ALL", tissue = interesting_tissues, r2d = "r2", p_threshold = 0.05, token = "b2735a858324")

# outputs
out.gtex = all_res_sqtl
save(out.gtex, file = snps_info_path)

