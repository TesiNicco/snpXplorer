# basic paths
MAIN = "/root/snpXplorer/AnnotateMe/"
MAIN_SNP = "/root/snpXplorer/snpXplorer_v3/"
args = commandArgs(trailingOnly=TRUE)

# i is the index of the snp_list vector of snps -- genome_version is either hg19 or hg38 -- input_type is locus for chr:pos -- tab for chr pos -- adjusted for faster computations (library-wise)
function_caddAnnot_v6 <- function(i, snp_list, genome_version, input_type, random_num){
  # get data depending on input type
  if (input_type == "locus"){
    input = snp_list[i]
  }

  # get genome version
  if (genome_version %in% c("hg19", "GRCh37")){
    genV = "GRCh37-v1.6_inclAnno"
  } else if (genome_version %in% c("hg38", "GRCh38")){
    genV = "GRCh38-v1.6_inclAnno"
  }

  # run cadd
  system(paste0("curl -i https://cadd.gs.washington.edu/api/v1.0/", genV, "/", input, " | tail -n +8 > tmp_", random_num, ".json"))

  # read json file
  json_data <- rjson::fromJSON(file = paste0("tmp_", random_num, ".json"))

  # parse file -- there are multiple objects, let's see what is different
  tmp_df = as.data.frame(matrix(data = NA, nrow = 1, ncol = 18))
  colnames(tmp_df) <- c("Alt", "AnnoType", "CCDS", "CDSpos", "Chrom", "ConsDetail", "ConsScore", "Consequence", "CpG", "Dist2Mutation",
                        "EnsembleRegulatoryFeature", "GeneID", "GeneName", "PHRED", "Pos", "Ref", "Type", "Locus")

  # loop on file
  for (i in 1:length(json_data)){
    tmp_dt = json_data[[i]]
    for (j in 1:length(colnames(tmp_df))){
      xx = tmp_dt[colnames(tmp_df)[j]][[1]]
      if (is.null(xx)){ tmp_df[i, j] <- NA } else {tmp_df[i, j] <- xx}
    }
  }

  # also get annotation information
  info <- paste0(unique(tmp_df$Consequence[!is.na(tmp_df$Consequence)]), collapse=",")
  g <- paste0(unique(tmp_df$GeneName[!is.na(tmp_df$GeneName)]), collapse=",")

  # define results df
  res = data.frame(Chrom = NA, Pos = NA, Ref = NA, Alt = NA, Type = NA, Conseq = NA, Genes = NA, Phred = NA)

  # check if culprit gene was available
  info_splt <- unlist(strsplit(info, ","))
  culprit <- unique(info_splt %in% c("NON_SYNONYMOUS", "SYNONYMOUS", "STOP_GAINED", "MISSENSE"))    #if the variant is not missense, synonymous, stop_gained
  if (TRUE %in% culprit){
    grp <- grep(TRUE, culprit)
    info <- paste0(unique(tmp_df$Consequence[grp], collapse=","))
    g <- paste0(unique(tmp_df$GeneName[grp], collapse=","))
    res$Conseq <- info
    res$Genes <- g
  } else {
    if (info != "") { res$Conseq <- info }
    if (g != "") { res$Genes <- g }
  }

  # finally assign other informations
  res$Chrom <- paste0(unique(tmp_df$Chrom), collapse = ",")
  res$Pos <- paste0(unique(tmp_df$Pos), collapse = ",")
  res$Ref <- paste0(unique(tmp_df$Ref), collapse = ",")
  res$Alt <- paste0(unique(tmp_df$Alt), collapse = ",")
  res$Type <- paste0(unique(tmp_df$Type), collapse = ",")
  res$Phred <- max(as.numeric(tmp_df$PHRED))
  res$Locus <- paste(res$Chrom, res$Pos, sep=":")

  return(res)
}

## function to read and store cadd files -- API cadd v2 -- adjusted for faster computations (library-wise)
readcadd_v2 <- function(data, MAIN, ld_info, ref_genome, random_num){
  # run function in mp using cadd api depending on ref genome
  # annotate the snp_list
  snp_list = paste(data$chr, data$pos, sep=":")
  all_res = data.table::rbindlist(parallel::mclapply(X = 1:length(snp_list), FUN = function_caddAnnot_v6, snp_list = snp_list, genome_version = "hg19", input_type = "locus", random_num = random_num, mc.cores = 1))
  all_res = all_res[, c("Locus", "Conseq", "Genes", "Phred")]
  colnames(all_res) = c("locus", "snp_conseq", "snp_conseq_gene", "phred")
  return(all_res)
}

## function to put information about snp consequences with the permutation sets
mergeInfo_generic <- function(info.sb){
  #assign flag for coding snp --> culprit gene known
  info.sb$coding_snp <- NA
  for (j in 1:nrow(info.sb)){
    if (!is.na(info.sb$snp_conseq[j])){
      info <- strsplit(info.sb$snp_conseq[j], ",")[[1]]
      culprit <- unique(info %in% c("NON_SYNONYMOUS", "SYNONYMOUS", "STOP_GAINED", "MISSENSE"))    #if the variant is not missense, synonymous, stop_gained
      if (TRUE %in% culprit){
        info.sb$coding_snp[j] <- "yes"
      }
    }
  }
  return(info.sb)
}

# read arguments
random_num = args[1]
ref_genome = args[2]
ld_path = args[3]
snps_info_path = args[4]
# run function to read snps info
load(snps_info_path)
load(ld_path)
# run cadd annotation
cadd <- readcadd_v2(data, MAIN, ld_info, ref_genome, random_num)
# merge cadd with data
data = merge(data, cadd, by = "locus")
cadd$coding_snp = NA
# then run the mergeInfo_generic to add coding_snps [yes/no]
out.annot <- mergeInfo_generic(info.sb=data)
out.annot <- out.annot[!duplicated(out.annot$locus),]

# outputs
save(out.annot, file = snps_info_path)
