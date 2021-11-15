# basic paths
MAIN = "/root/snpXplorer/AnnotateMe/"
MAIN_SNP = "/root/snpXplorer/snpXplorer_v3/"
args = commandArgs(trailingOnly=TRUE)

## function to read input set of snp given the name and the input type -- adjusted for faster computations (library-wise)
readSNPs <- function(fname, ftype, MAIN, ref_version, analysis_type){
  ## read input file
  d <- data.table::fread(fname, h=F)

  ## check number of rows -- if larger than 1000 SNPs, stop
  if (nrow(d) > 1000 & analysis_type == "enrichment"){
    d = NA
    ftype = 99
  } else if (nrow(d) > 10000 & analysis_type != "enrichment"){
    d = NA
    ftype = 99
  }

  ## check type: type 1 --> chr:pos; type 2 --> chr   pos; type 3 --> rsid; type 4 --> rsid_A1
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
      sb = data.frame(locus = paste(final$chr_n, final$start.y, sep = ":"), chr = final$chr_n, pos = as.numeric(final$start.y))
      d = sb
      write.table(d$locus, fname, quote=F, row.names=F, col.names = F)
    }
    info <- system(paste0("grep -w -F -f ", fname, " ", MAIN, "INPUTS_OTHER/1000G_frequencies/chrAll_locus.afreq"), intern = T)
    info = as.data.frame(stringr::str_split_fixed(info, "\t", 8))
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
      write.table(d$locus, fname, quote=F, row.names=F, col.names=F)
    } else {
      write.table(d$locus, fname, quote=F, row.names=F, col.names=F)
    }
    info <- system(paste0("grep -w -F -f ", fname, " ", MAIN, "INPUTS_OTHER/1000G_frequencies/chrAll_locus.afreq"), intern = T)
    info = as.data.frame(stringr::str_split_fixed(info, "\t", 8))
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
    # grep rsid
    cmd <- paste0("grep -w -F -f ", fname, " ", MAIN, "INPUTS_OTHER/1000G_frequencies/chrAll.afreq")
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
  return(d)
}

# read arguments
fname = args[1]
ftype = args[2]
ref_version = args[3]
analysis_type = args[4]
random_num = args[5]
outpath = paste0("/root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/tmp_snpsInfo.RData")
# run function to read
data <- try(readSNPs(fname, ftype, MAIN, ref_version, analysis_type), silent = T)
data = data[!is.na(data$pos),]
save(data, file = outpath)
cat(outpath)
