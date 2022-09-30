# harmonize TOPMED annotation data and make them suitable for tabix
# library
  library(GenomicRanges)
  library(liftOver)
  library(data.table)
  library(stringr)
# functions
# function to liftover data
  liftOver_data <- function(chrom, start, end, p = NULL, from){
    if (is.null(p)){ df <- data.frame(chr=paste("chr", chrom, sep=""), start=start, end=end) } else { df <- data.frame(chr=paste("chr", chrom, sep=""), start=start, end=end, p = p) }
    cat('--- adding tag to SNPs\n')
    df = df[order(df$start),]
    df$group = seq(1, nrow(df))
    # change to GR class object
    cat('--- creating GR object\n')
    gr <- makeGRangesFromDataFrame(df)
    # set chain file
    if (from == 'hg19') { chain <- import.chain("/root/snpXplorer/data/databases/hg19ToHg38.over.chain") } else { chain <- import.chain("/root/snpXplorer/data/databases/hg38ToHg19.over.chain") }
    # change coordinates
    cat(paste0('--- lifting over from ', from, '\n'))
    gr_lifted <- liftOver(gr, chain)
    # back to dataframe and clean it
    df_lifted <- as.data.frame(gr_lifted)
    df_lifted = df_lifted[which(df_lifted$seqnames == paste0('chr', chrom[1])),]
    cat('--- merging lifted object\n')
    df_lifted = df_lifted[, c("group", "start", "end")]
    df_lifted = merge(df_lifted, df, by = "group")
    colnames(df_lifted) = c('group', 'start_hg19', 'end_hg19', 'chr', 'start_hg38', 'end_hg38')
    df_final = df_lifted[, c('start_hg19', 'start_hg38')]
    df_final = df_final[!duplicated(df_final$start_hg38),]
    return(df_final)
  }

# define input and output directories
data_path = '/root/snpXplorer/AnnotateMe/INPUTS_OTHER/ld_topmed'
out_dir = '/root/snpXplorer/data/databases/topmed_ld_and_annotation'

# this is for the annotation
for (chr in seq(12, 19)){
    cat(paste0('## started with chromosome --> ', chr, '\n'))
    # read data
    cat('## reading topmed data\n')
    topmed = fread(paste0(data_path, '/EUR_chr', chr, '_no_filter_0.2_1000000_info_annotation.csv.gz'), sep = ',')
    # sort by position
    cat('## sorting topmed data\n')
    topmed = topmed[order(topmed$Position),]
    # liftover
    cat('## lifting over topmed data\n')
    topmed_lifted = liftOver_data(chrom = rep(chr, nrow(topmed)), start = topmed$Position, end = topmed$Position + 1, from = 'hg38')
    # merge with entire data
    cat('## merging topmed data\n')
    topmed_annotated = merge(topmed, topmed_lifted, by.x = 'Position', by.y = 'start_hg38', all.x = T)
    topmed_annotated$chrom = chr
    # sort by position
    topmed_annotated = topmed_annotated[order(topmed_annotated$Position),]
    topmed_annotated_hg19 = topmed_annotated[!is.na(topmed_annotated$start_hg19),]
    topmed_annotated_hg19 = topmed_annotated_hg19[order(topmed_annotated_hg19$start_hg19),]
    # write tables
    cat('## writing output table\n')
    write.table(topmed_annotated, paste0(out_dir, '/chr', chr, '_topmed_annotation_hg38.txt'), quote=F, row.names = F, sep = "\t")
    write.table(topmed_annotated_hg19, paste0(out_dir, '/chr', chr, '_topmed_annotation_hg37.txt'), quote=F, row.names = F, sep = "\t")
    # compress
    cat('## compressing output table\n')
    system(paste0('bgzip ', out_dir, '/chr', chr, '_topmed_annotation_hg38.txt'))
    system(paste0('bgzip ', out_dir, '/chr', chr, '_topmed_annotation_hg37.txt'))
    # and tabix
    cat('## tabix output table\n')
    chr_index = which(colnames(topmed_annotated) == 'chrom')
    bp_index = which(colnames(topmed_annotated) == 'Position')
    bp19_index = which(colnames(topmed_annotated_hg19) == 'start_hg19')
    system(paste0('tabix -S 1 -s ', chr_index, ' -b ', bp_index, ' -e -', bp_index, ' ', out_dir, '/chr', chr, '_topmed_annotation_hg38.txt.gz'))
    system(paste0('tabix -S 1 -s ', chr_index, ' -b ', bp19_index, ' -e -', bp19_index, ' ', out_dir, '/chr', chr, '_topmed_annotation_hg37.txt.gz'))
    cat(paste0('## done with chromosome --> ', chr, '\n\n'))
    # clean memory
    rm(topmed); rm(topmed_annotated); rm(topmed_annotated_hg19); rm(topmed_lifted); gc()
}
