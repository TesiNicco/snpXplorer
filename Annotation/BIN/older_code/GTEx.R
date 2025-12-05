# library
library(LDlinkR)
library(stringr)

# basic paths
MAIN = "/Annotation/"
MAIN_SNP = "/Annotation/RUNS/"
args = commandArgs(trailingOnly=TRUE)

## function to match eqtl for each chromosome
matchEQTL <- function(i, mapping, ensembl, interesting_tissues, random_num){
  cat(paste0("## Doing eQTL annotation of chr", i, "\n"))
  tmp = mapping[which(mapping$chr == i), ]
  # write temporary output with eqtl ids to match
  write.table(tmp$code_gtex, paste0(MAIN, "INPUTS_OTHER/summary_eqtls/tmp_chr", i, "_", random_num, ".txt"), quote=F, row.names=F, col.names=F)
  # command to grep
  cmd = paste0("zgrep -F -f ", MAIN, "INPUTS_OTHER/summary_eqtls/tmp_chr", i, "_", random_num, ".txt ", MAIN, "INPUTS_OTHER/summary_eqtls/chr", i, "_summary_eqtls.txt.gz")
  #tmp_res = system(cmd, intern = T)
  tmp_res <- tryCatch({ system(cmd, intern = T) }, error=function(cond){ return(NA) }, warning=function(cond){ return(NA) })    
  # remove temporary data
  system(paste0("rm ", MAIN, "INPUTS_OTHER/summary_eqtls/tmp_chr", i, "_", random_num, ".txt"))
  # create variables for assignment
  tmp$eqtl = NA
  tmp$eqtl_tissue = NA
  # main loop to assign eqtls to snps
  if (!is.na(tmp_res[1])){
    # modify structure of data
    tmp_res = data.frame(stringr::str_split_fixed(tmp_res, "\t", 2), stringsAsFactors = F)
    for (i in 1:nrow(tmp_res)){
      # take snp id
      snp_id = paste0(stringr::str_replace_all(paste(stringr::str_split_fixed(tmp_res$X1[i], "_", 5)[, 1:2], collapse = "_"), "chr", ""), "_")
      # extract tissues and genes
      tissues_genes = stringr::str_split(unlist(strsplit(tmp_res$X2[i], ";")), "_")
      tissues = c()
      genes = c()
      for (x in tissues_genes){ tissues = c(tissues, paste(x[1:(length(x)-3)], collapse = "_")); genes = c(genes, stringr::str_split_fixed(x[length(x)-2], "\\.", 2)[1]) }
      tissues_genes_df = data.frame(tissue = tissues, gene = genes)
      tissues_genes_df = merge(tissues_genes_df, ensembl, by.x = "gene", by.y = "Gene stable ID", all.x = T)
      # filter according to tissues of interest
      if (!("all_tissues" %in% interesting_tissues)){ tissues_genes_df = tissues_genes_df[which(tissues_genes_df$tissue %in% interesting_tissues),] }
      # take union of all tissues and genes
      all_tissues = paste0(tissues_genes_df$tissue, collapse = ",")
      all_genes = paste0(tissues_genes_df$"Gene name", collapse = ",")
      # finally assign to relative snp
      tmp$eqtl[which(tmp$code_gtex == snp_id)] = all_genes
      tmp$eqtl_tissue[which(tmp$code_gtex == snp_id)] = all_tissues
    }
  }
  return(tmp)
}

## grep all eqtl in the new implementation -- september 2021
GTEx_me_generic_allTissues_newImplementation <- function(mapping, interesting_tissues){
  # first, need to change reference as gtex v8 is hg38
  df <- data.frame(chr=paste("chr", mapping$chr, sep=""), start=mapping$pos, end=as.numeric(mapping$pos)+1)
  # change to GR class object
  gr <- GenomicRanges::makeGRangesFromDataFrame(df)
  # set chain file
  chain <- rtracklayer::import.chain(paste(MAIN, "BIN/hg19ToHg38.over.chain", sep=""))
  # change coordinates
  gr_hg38 <- rtracklayer::liftOver(gr, chain)
  # back to dataframe and clean it
  df_hg38 <- as.data.frame(gr_hg38)
  # before doing the match below, manage in case there is some mismatch
  df_in_hg38 = df[df_hg38$group, ]
  df_in_hg38$chr = stringr::str_split_fixed(df_in_hg38$chr, "chr", 2)[, 2]
  sb_df <- data.frame(chr_hg38 = as.character(df_hg38$seqnames), pos_hg38 = as.character(df_hg38$start), locus = paste(df_in_hg38$chr, df_in_hg38$start, sep=":"), code_gtex = paste(as.character(df_in_hg38$chr), "_", as.character(df_hg38$start), "_", sep=""))
  # merge with mapping using locus column
  mapping <- merge(mapping, sb_df, by="locus", all.x = T)
  # read conversion file from ensembl to gene name
  ensembl = data.table::fread(paste0(MAIN, "INPUTS_OTHER/Ensemble_to_GeneName.txt"), h=T, sep = "\t")
  # run function for each chromosome using mclapply
  chr_list = unique(mapping$chr)
  all_res_eqtl = data.table::rbindlist(parallel::mclapply(X = chr_list, FUN = matchEQTL, mapping = mapping, ensembl = ensembl, interesting_tissues = interesting_tissues, random_num = random_num, mc.cores = 1))
  return(all_res_eqtl)
}

# read arguments
snps_info_path = args[1]
interesting_tissues = str_split(args[2], ',')[[1]]
random_num = args[3]
# run function to read snps info
load(snps_info_path)
# run function
out.gtex <- GTEx_me_generic_allTissues_newImplementation(mapping = out.annot, interesting_tissues = interesting_tissues)
# replace empty cells with NA -- this has to be different if the number of rows is 1 only
if (nrow(out.gtex) == 1){
  out.gtex = as.data.frame(t(apply(out.gtex, 2, function(x) gsub("^$|^ $", NA, x))))
} else {
  out.gtex = as.data.frame(apply(out.gtex, 2, function(x) gsub("^$|^ $", NA, x)))
}
# we can also use LDlink -- this will give results for the snps of interest + all snps in LD with that
#out.annot.sb = out.annot[grep("rs", out.annot$ID),]
#miss = out.annot[which(!(out.annot$locus %in% out.annot.sb$locus)),]
#gtex_info = LDexpress(out.annot.sb$ID, pop = "ALL", tissue = interesting_tissues, r2d = "r2", p_threshold = 0.05, token = "b2735a858324")

# outputs
save(out.gtex, file = snps_info_path)

