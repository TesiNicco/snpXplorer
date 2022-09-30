# libraries
library(parallel)
library(data.table)
library(stringr)
library(bedr)

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

# function to do the lookup in the entire cadd dataset (offline)
lookup_cadd_updated = function(i, data, random_num){
    # restrict to chromosome of interest
    tmp = data[which(data$chr == i),]
    # we are going to grep information using chromosome;position; --> write that file
    tmp$towrite = paste0(tmp$chr, ';', tmp$pos, ';')
    write.table(tmp$towrite, paste0('RESULTS_', random_num, '/tmp_chr', i, '_togrep.txt'), quote=F, row.names=F, col.names=F)
    # command to grep information
    cmd = paste0('zgrep -f RESULTS_', random_num, '/tmp_chr', i, '_togrep.txt ', MAIN, '/INPUTS_OTHER/CADD_updated/cadd_annotation_chr', i, '.txt.gz'); caddannot = system(cmd, intern=T)
    # then we should loop across matches to extract relevant information
    caddannot = data.frame(str_split_fixed(caddannot, ';', 13)); colnames(caddannot) = c('chrom', 'pos', 'ref', 'alt', 'type', 'annotype', 'consequence', 'genename', 'conseq_details', 'conseq_score', 'dist_mutation', 'geneid', 'phred')
    # loop on the lines and prioritize results based on highest phred score
    prioritized = list()
    for (x in 1:nrow(caddannot)){
        tmp = caddannot[x, ]; tmp_phred = as.numeric(unlist(strsplit(tmp$phred, '\\|'))); highest_phred_index = which(tmp_phred == max(tmp_phred))[1]
        tmp_df = tmp[, c('chrom', 'pos', 'ref', 'alt', 'type')]; tmp_df$consequence = unlist(strsplit(tmp$consequence, '\\|'))[highest_phred_index]
        tmp_df$genename = unlist(strsplit(tmp$genename, '\\|'))[highest_phred_index]; tmp_df$phred = max(tmp_phred)
        prioritized[[(length(prioritized) + 1)]] = tmp_df
    }
    prioritized = rbindlist(prioritized)
    res = data.frame(locus = paste0(prioritized$chrom, ':', prioritized$pos), snp_conseq = prioritized$consequence, snp_conseq_gene = prioritized$genename, phred = prioritized$phred)
    return(res)
}

# function to put information about snp consequences with the permutation sets
mergeInfo_generic_updated <- function(info.sb){
  #assign flag for coding snp --> culprit gene known
  info.sb$coding_snp <- NA
  for (j in 1:nrow(info.sb)){
    if (!is.na(info.sb$snp_conseq[j])){
      info <- strsplit(info.sb$snp_conseq[j], ",|\\|")[[1]]
      culprit <- unique(info %in% c("NON_SYNONYMOUS", "SYNONYMOUS", "STOP_GAINED", "MISSENSE", 'missense', 'synonymous', 'non_synonymous', 'stop_gained', 'coding_sequence', 'start_lost'))    #if the variant is not missense, synonymous, stop_gained
      if (TRUE %in% culprit){
        info.sb$coding_snp[j] <- "yes"
      }
    }
  }
  return(info.sb)
}

# function to use tabix
lookup_cadd_tabix = function(i, data, random_num){
  # restrict to chromosome of interest
  tmp = data[which(data$chr == i),]
  tmp$roi = paste0(tmp$chr, ':', as.numeric(tmp$pos) - 1, '-', tmp$pos)
  tmp = tmp[!is.na(tmp$pos),]
  # tabix command
  caddinfo = tabix(tmp$roi, paste0(MAIN, '/INPUTS_OTHER/CADD_updated/cadd_annotation_chr', i, '_tabix.txt.gz'), check.chr = F, verbose = FALSE)
  colnames(caddinfo) = c('chrom', 'pos', 'ref', 'alt', 'type', 'annotype', 'consequence', 'genename', 'conseq_details', 'conseq_score', 'dist_mutation', 'geneid', 'phred')
  # loop on the lines and prioritize results based on highest phred score
  prioritized = list()
  for (x in 1:nrow(caddinfo)){
    tmp = caddinfo[x, ]; tmp_phred = as.numeric(unlist(strsplit(tmp$phred, '\\|'))); highest_phred_index = which(tmp_phred == max(tmp_phred))[1]
    tmp_df = tmp[, c('chrom', 'pos', 'ref', 'alt', 'type')]; tmp_df$consequence = unlist(strsplit(tmp$consequence, '\\|'))[highest_phred_index]
    tmp_df$genename = unlist(strsplit(tmp$genename, '\\|'))[highest_phred_index]; tmp_df$phred = max(tmp_phred)
    prioritized[[(length(prioritized) + 1)]] = tmp_df
  }
  prioritized = rbindlist(prioritized)
  res = data.frame(locus = paste0(prioritized$chrom, ':', prioritized$pos), snp_conseq = prioritized$consequence, snp_conseq_gene = prioritized$genename, phred = prioritized$phred)
  return(res)
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
#cadd <- readcadd_v2(data, MAIN, ld_info, ref_genome, random_num)       # previous command based on the online grep of CADD
# run function in multiprocessing to annotate snps with cadd using tabix
#cadd_annotation_results = rbindlist(mclapply(unique(data$chr), lookup_cadd_updated, data = data, random_num = random_num, mc.cores = 2))
cadd_annotation_results = rbindlist(mclapply(unique(data$chr), lookup_cadd_tabix, data = data, random_num = random_num, mc.cores = 2))
# then merge with data and annotate whether the snp is coding or not
data = merge(data, cadd_annotation_results, by = 'locus', all.x = T)
data_codingInfo = mergeInfo_generic_updated(data)
# finally remove duplicates and save
out.annot = data_codingInfo[!duplicated(data_codingInfo$locus),]

# outputs
save(out.annot, file = snps_info_path)
######################################################################
