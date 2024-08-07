# basic paths
MAIN = "/Annotation/"
MAIN_SNP = "/Annotation/RUNS/"
args = commandArgs(trailingOnly=TRUE)
library(stringr)
library(parallel)
library(data.table)
library(bedr)

## use tabix instead of grepping
findTabix <- function(i, d, MAIN){
  d$chrom = stringr::str_split_fixed(d$V1, ':', 2)[, 1]
  d$pos = stringr::str_split_fixed(d$V1, ':', 2)[, 2]
  # restrict to chromosome of interest
  tmp = d[which(d$chrom == i),]
  # sort by position
  tmp$pos = as.numeric(tmp$pos)
  tmp = tmp[order(tmp$pos),]
  tmp$roi = paste0(tmp$chrom, ':', as.numeric(tmp$pos) - 1, '-', tmp$pos)
  tmp = tmp[!is.na(tmp$pos),]
  # tabix command
  info = tabix(tmp$roi, paste0(MAIN, '/INPUTS_OTHER/1000G_frequencies/chr', i, '_locus.afreq.txt.gz'), check.chr = F, verbose = FALSE)
  return(info)
}

## function to read input set of snp given the name and the input type -- adjusted for faster computations (library-wise)
readSNPs <- function(fname, ftype, MAIN, ref_version, analysis_type){
  ## read input file
  d <- data.table::fread(fname, h=F, fill=T, stringsAsFactors=F)

  ## remove empty lines
  d <- d[!apply(d == "", 1, all),]

  ## check number of rows -- if larger than 1000 SNPs, stop
  if (nrow(d) > 6000 & analysis_type == "enrichment"){
    d = NA
    ftype = 99
  } else if (nrow(d) > 100000 & analysis_type != "enrichment"){
    d = NA
    ftype = 99
  }

  ## check type: type 1 --> chr:pos; type 2 --> chr   pos; type 3 --> rsid; type 4 --> rsid_A1
  if (ftype != 99){
    if (ftype == 1){
      # if requested reference was hg38, lift to hg19
      if (ref_version == "GRCh38"){
        tmp <- stringr::str_split_fixed(d$V1, ":", 2)
        d$chr <- as.numeric(as.character(tmp[, 1]))
        d$pos <- as.numeric(as.character(tmp[, 2]))
        colnames(d) <- c("locus", "chr", "pos")
        df <- data.frame(chr=paste("chr", d$chr, sep=""), start=d$pos, end=d$pos+1)
        # change to GR class object
        gr <- GenomicRanges::makeGRangesFromDataFrame(df)
        # set chain file
        chain <- rtracklayer::import.chain(paste(MAIN, "BIN/hg38ToHg19.over.chain", sep=""))
        # change coordinates
        gr_hg38 <- rtracklayer::liftOver(gr, chain)
        # back to dataframe and clean it
        df_hg38 <- as.data.frame(gr_hg38)
        # before doing the match below, manage in case there is some mismatch
        df$group = seq(1, nrow(df))
        final = merge(df, df_hg38, by = "group")
        final$chr_n = stringr::str_split_fixed(final$chr, "chr", 2)[, 2]
        rm(d)
        system(paste0('mv ', fname, ' ', str_replace_all(fname, '_input_', '_originalInput_')))
        d = data.frame(locus = paste(final$chr_n, final$start.y, sep = ":"), chr = final$chr_n, pos = as.numeric(final$start.y))
        d$chr = NULL
        d$pos = NULL
        colnames(d) = 'V1'
        write.table(d$V1, fname, quote=F, row.names=F, col.names = F)
      }
      chrom_list = unique(stringr::str_split_fixed(d$V1, ':', 2)[, 1])
      info = rbindlist(mclapply(chrom_list, findTabix, d = d, MAIN = MAIN, mc.cores = 2))
      colnames(info) <- c("chr", "pos", "ID", "ref", "alt", "ALT_FREQS", "n", "locus")
      miss <- d$V1[which(!(d$V1 %in% info$locus))]
      info = info[, c("chr", "pos", "ID", "ALT_FREQS")]
      info$locus <- paste(info$chr, info$pos, sep=":")
      miss_df <- data.frame(chr = rep("NA", length(miss)), pos = rep(NA, length(miss)), ID = miss, ALT_FREQS = rep(NA, length(miss)), locus = rep(NA, length(miss)))
      d <- rbind(info, miss_df)
      d <- d[!duplicated(d$ID),]
    } else if (ftype == 2){
      d$locus <- paste(d$V1, d$V2, sep=":")
      colnames(d) <- c("chr", "pos", "locus")
      d = d[, 'locus']
      colnames(d) = 'V1'
      # if requested reference was hg38, lift to hg19
      if (ref_version == "GRCh38"){
        df <- data.frame(chr=paste("chr", d$chr, sep=""), start=d$pos, end=d$pos+1)
        # change to GR class object
        gr <- GenomicRanges::makeGRangesFromDataFrame(df)
        # set chain file
        chain <- rtracklayer::import.chain(paste(MAIN, "BIN/hg38ToHg19.over.chain", sep=""))
        # change coordinates
        gr_hg38 <- rtracklayer::liftOver(gr, chain)
        # back to dataframe and clean it
        df_hg38 <- as.data.frame(gr_hg38)
        # before doing the match below, manage in case there is some mismatch
        df$group = seq(1, nrow(df))
        final = merge(df, df_hg38, by = "group")
        final$chr_n = stringr::str_split_fixed(final$chr, "chr", 2)[, 2]
        sb = data.frame(locus = paste(final$chr_n, final$start.y, sep = ":"), chr = final$chr_n, pos = as.numeric(final$start.y))
        d = sb
        d$chr = NULL
        d$pos = NULL
        colnames(d) = 'V1'
        write.table(d$V1, fname, quote=F, row.names=F, col.names = F)
      } else {
        write.table(d$V1, fname, quote=F, row.names=F, col.names=F)
      }
      chrom_list = unique(stringr::str_split_fixed(d$V1, ':', 2)[, 1])
      info = rbindlist(mclapply(chrom_list, findTabix, d = d, MAIN = MAIN, mc.cores = 2))
      colnames(info) <- c("chr", "pos", "ID", "ref", "alt", "ALT_FREQS", "n", "locus")
      miss <- d$V1[which(!(d$V1 %in% info$locus))]
      info = info[, c("chr", "pos", "ID", "ALT_FREQS")]
      info$locus <- paste(info$chr, info$pos, sep=":")
      miss_df <- data.frame(chr = rep("NA", length(miss)), pos = rep(NA, length(miss)), ID = miss, ALT_FREQS = rep(NA, length(miss)), locus = rep(NA, length(miss)))
      d <- rbind(info, miss_df)
      d <- d[!duplicated(d$ID),]
    } else if (ftype %in% c(3, 4)){
      if (ftype == 4){ d <- as.data.frame(stringr::str_split_fixed(d$V1, "_", 2)); d$V2 <- NULL }
      #ref <- fread(paste(MAIN, "INPUTS_OTHER/1000G_frequencies/chrAll.afreq.gz", sep=""), h=T, showProgress=FALSE)
      #info <- ref[which(ref$ID %in% d$V1), c("#CHROM", "POS", "ID", "ALT_FREQS")]
      # before grepping rsid, we should remove empty spaces otherwise there will be issues
      colnames(d) = "V1"
      d$V1 = str_replace_all(d$V1, " ", "")
      write.table(d, fname, quote=F, row.names=F, col.names=F)
      # grep rsid
      cmd <- paste0("zgrep -w -F -f ", fname, " ", MAIN, "INPUTS_OTHER/1000G_frequencies/chrAll.afreq.gz")
      info = system(cmd, intern=T)
      info = as.data.frame(stringr::str_split_fixed(info, "\t", 7))
      colnames(info) <- c("chr", "pos", "ID", "ref", "alt", "ALT_FREQS", "n")
      miss <- d$V1[which(!(d$V1 %in% info$ID))]
      info = info[, c("chr", "pos", "ID", "ALT_FREQS")]
      info$locus <- paste(info$chr, info$pos, sep=":")
      miss_df <- data.frame(chr = rep("NA", length(miss)), pos = rep(NA, length(miss)), ID = miss, ALT_FREQS = rep(NA, length(miss)), locus = rep(NA, length(miss)))
      d <- rbind(info, miss_df)
      d <- d[!duplicated(d$ID),]
    }
    d = d[which(d$chr != "#CHROM"),]
    #d = d[!is.na(d$pos),]
  }
  return(d)
}

## function to liftover
liftOver_fun <- function(df){
  df$chr <- as.numeric(as.character(df$chr))
  df$pos <- as.numeric(as.character(df$pos))
  tmp_df <- data.frame(chr=paste("chr", df$chr, sep=""), start=df$pos, end=df$pos+1)
  # change to GR class object
  gr <- GenomicRanges::makeGRangesFromDataFrame(tmp_df)
  # set chain file
  chain <- rtracklayer::import.chain(paste(MAIN, "BIN/hg38ToHg19.over.chain", sep=""))
  # change coordinates
  gr_hg38 <- rtracklayer::liftOver(gr, chain)
  # back to dataframe and clean it
  df_hg38 <- as.data.frame(gr_hg38)
  # before doing the match below, manage in case there is some mismatch
  df$group = seq(1, nrow(df))
  final = merge(df, df_hg38, by = "group")
  final = final[, c('chr', 'pos', 'ref', 'alt', 'locus', 'rsid', 'start')]
  colnames(final) = c('chr', 'pos_hg38', 'ref', 'alt', 'locus_hg38', 'rsid', 'pos_hg19')
  return(final)
}

## function to use topmed information to extract snp information
readSNPs_alternative <- function(fname, ftype, MAIN, ref_version, analysis_type){
  ## read input file
  d <- data.table::fread(fname, h=F, fill=T)

  ## remove empty lines
  d <- d[!apply(d == "", 1, all),]
  colnames(d) = "V1"
  d$V1 = str_replace_all(d$V1, " ", "")
  write.table(d, fname, quote=F, row.names=F, col.names=F)
  # grep rsid
  cmd <- paste0("zgrep -w -F -f ", fname, " ", MAIN, "INPUTS_OTHER/1000G_frequencies/MARKED_RSID_TO_POSITION_dbSNPv151_and_dbSNPv153VEP_20200626.txt.gz")
  info = system(cmd, intern=T)
  info = as.data.frame(stringr::str_split_fixed(info, "\t", 3))
  if (is.data.frame(info)){
    df = data.frame(stringr::str_split_fixed(info$V1, ":", 4))
    colnames(df) = c("chr", "pos", "ref", "alt")
    df$chr = str_replace_all(df$chr, 'chr', '')
    df$locus = info$V1
    df$rsid = info$V2
    df_lifted = liftOver_fun(df)
    df_lifted = df_lifted[, c('chr', 'pos_hg19', 'rsid', 'ref', 'alt')]
    df_lifted$ALT_FREQS = NA
    df_lifted$n = NA
    colnames(df_lifted) <- c("chr", "pos", "ID", "ref", "alt", "ALT_FREQS", "n")
    miss <- d$V1[which(!(d$V1 %in% df_lifted$ID))]
    info = df_lifted[, c("chr", "pos", "ID", "ALT_FREQS")]
    info$locus <- paste(info$chr, info$pos, sep=":")
    miss_df <- data.frame(chr = rep("NA", length(miss)), pos = rep(NA, length(miss)), ID = miss, ALT_FREQS = rep(NA, length(miss)), locus = rep(NA, length(miss)))
    d <- rbind(info, miss_df)
    d <- d[!duplicated(d$ID),]
    data = d
  } else {
    data = NA
  }
  return(data)
}

## function to use topmed information in case of hg38 and position prompted
readSNPs_alternative_position <- function(data_na, MAIN_SNP, random_num){
  # grep information if available
  write.table(paste0('chr', data_na$ID, ':'), paste0(MAIN_SNP, 'RESULTS_', random_num, '/tmp_snppos.txt'), quote=F, row.names=F, col.names=F)
  cmd <- paste0("zgrep -F -f ", MAIN_SNP, 'RESULTS_', random_num, '/tmp_snppos.txt ', MAIN, "INPUTS_OTHER/1000G_frequencies/MARKED_RSID_TO_POSITION_dbSNPv151_and_dbSNPv153VEP_20200626.txt.gz")
  info = system(cmd, intern=T)
  info_df = data.frame(stringr::str_split_fixed(info, '\t', 3))
  colnames(info_df) = c('id', 'rsid', 'rsid2')
  info_df
  info_df$chr = str_replace_all(stringr::str_split_fixed(info_df$id, ':', 4)[, 1], 'chr', '')
  info_df$pos = str_split_fixed(info_df$id, ':', 4)[, 2]
  info_df$chr <- as.numeric(as.character(info_df$chr))
  info_df$pos <- as.numeric(as.character(info_df$pos))
  tmp_df <- data.frame(chr=paste("chr", info_df$chr, sep=""), start=info_df$pos, end=info_df$pos+1)
  # change to GR class object
  gr <- GenomicRanges::makeGRangesFromDataFrame(tmp_df)
  # set chain file
  chain <- rtracklayer::import.chain(paste(MAIN, "BIN/hg38ToHg19.over.chain", sep=""))
  # change coordinates
  gr_hg37 <- rtracklayer::liftOver(gr, chain)
  # back to dataframe and clean it
  df_hg37 <- as.data.frame(gr_hg37)
  df_hg37 = df_hg37[, c('group', 'seqnames', 'start')]; colnames(df_hg37) = c('group', 'chr_hg37', 'pos_hg37')
  # before doing the match below, manage in case there is some mismatch
  info_df$group = seq(1, nrow(info_df))
  final = merge(info_df, df_hg37, by = "group")
  final$ALT_FREQS = NA
  final$locus_hg38 = paste(final$chr, final$pos, sep=":")
  remaining_missing = data_na[which(!(data_na$ID %in% final$locus_hg38)),]
  final = final[, c('chr', 'pos_hg37', 'rsid', 'ALT_FREQS')]
  colnames(final) = c('chr', 'pos', 'ID', 'ALT_FREQS')
  final$locus = paste(final$chr, final$pos, sep=":")
  final = rbind(final, remaining_missing)
  return(final)
}

# read arguments
fname = args[1]
ftype = args[2]
ref_version = args[3]
analysis_type = args[4]
random_num = args[5]
outpath = paste0("/Annotation/RUNS/RESULTS_", random_num, "/tmp_snpsInfo.RData")
# run function to read
data <- try(readSNPs(fname, ftype, MAIN, ref_version, analysis_type), silent = T)

# check if there were matches, in case there were no matches, try with the other dataset
if (is.data.frame(data) && nrow(data) == 0 && ftype == 3){ data = try(readSNPs_alternative(fname, ftype, MAIN, ref_version, analysis_type), silent = T) }
# also check with the other dataset for the missing annotations in case was requested ftype 3
if (ftype == 3 && is.data.frame(data)){
  missings = data[is.na(data$pos),]
  if (nrow(missings) >0){
    newname = str_replace_all(fname, ".txt", "_miss.txt")
    write.table(missings$ID, newname, quote=F, row.names=F, col.names=F)
    miss_data = tryCatch({ x = readSNPs_alternative(newname, ftype, MAIN, ref_version, analysis_type) }, error=function(cond) { return(NA) })
    if (is.data.frame(miss_data)){
        data = data[!is.na(data$pos),]
        data = rbind(data, miss_data)
    }
    system(paste0("rm ", newname))
  }
}
# if the reference genome was hg38, we can check that
#if (ref_version == 'GRCh38' && ftype == 1){
  #data_na = data[is.na(data$pos),]
  #data_clean = data[!is.na(data$pos),]
  #add_data = readSNPs_alternative_position(data_na, MAIN_SNP, random_num)
  #data = rbind(data_clean, add_data)
#}
if (is.data.frame(data) && ftype == 1){
  data = data.frame(data)
  missing = data[is.na(data$pos),]
  nonmissing = data[!is.na(data$pos),]
  if (nrow(missing) >0){
    missing$chr = str_replace_all(str_split_fixed(missing$ID, ':', 2)[, 1], 'chr', '')
    missing$pos = str_split_fixed(missing$ID, ':', 2)[, 2]
    missing$locus = paste(missing$chr, missing$pos, sep=":")
    data = rbind(missing, nonmissing)
  } else {
    data = nonmissing
  }
}
save(data, file = outpath)
cat(outpath)
