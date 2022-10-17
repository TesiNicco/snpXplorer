# libraries
library(data.table)
library(stringr)
library(liftOver)

# tabix all files
path = '/root/snpXplorer/data'
all_files = system(paste0('ls -d ', path, '/*/'), intern = T)
for (folder in all_files){
    if (folder != '/root/snpXplorer/data/databases/' && folder != '/root/snpXplorer/data/TABIX/' && folder != '/root/snpXplorer/data/GTEX/' && folder != '/root/snpXplorer/data/GWAS_catalog/'){
        filelist = system(paste0('ls ', folder, 'chr*gz'), intern = T)
        cat(paste0('** working on --> ', folder, '\n'))
        for (f in filelist){
            cat(paste0('**** working on --> ', f, '\n'))
            tmp = fread(f, h=T, stringsAsFactors = F)
            colnames(tmp) = c('CHR', 'POS', 'P')
            tmp = tmp[order(tmp$CHR, tmp$POS),]
            # annotate the snps now
            tmp_annot = fread(paste0('~/snpXplorer/data/databases/topmed_ld_and_annotation/chr', tmp$CHR[1], '_topmed_annotation_hg37.txt.gz'), h = T, stringsAsFactors = F)
            tmp_annot_red = tmp_annot[, c('Position', 'start_hg19', 'rsID', 'MAF', 'REF', 'ALT')]
            tmp_annotated = merge(tmp, tmp_annot_red, by.x = 'POS', by.y = 'start_hg19', all.x = T)
            # sort
            tmp_annotated = tmp_annotated[order(tmp_annotated$POS),]
            colnames(tmp_annotated) = c('POS', 'CHR', 'P', 'POS_HG38', 'RSID', 'MAF', 'REF', 'ALT')
            # write
            write.table(tmp_annotated, str_replace_all(f, '.gz', ''), quote = F, row.names = F, sep = '\t')
            system(paste0('bgzip -f ', str_replace_all(f, '.gz', '')))
            system(paste0('tabix -S 1 -s 2 -b 1 -e -1 ', f))
            # clean
            rm(tmp_annot); rm(tmp_annot_red); gc()
        }
    }
}
#########
path = '/root/snpXplorer/data'
all_files = system(paste0('ls -d ', path, '/*/'), intern = T)
for (folder in all_files){
    if (folder != '/root/snpXplorer/data/databases/' && folder != '/root/snpXplorer/data/TABIX/' && folder != '/root/snpXplorer/data/GTEX/' && folder != '/root/snpXplorer/data/GWAS_catalog/'){
        filelist = system(paste0('ls ', folder, 'chr*gz'), intern = T)
        cat(paste0('** working on --> ', folder, '\n'))
        for (f in filelist){
            cat(paste0('**** working on ', f, '\n'))
            tmp = fread(f, stringsAsFactors = F)
            tmp = tmp[!is.na(tmp$POS_HG38),]
            tmp = tmp[order(tmp$POS_HG38),]
            tmp$POS = NULL
            colnames(tmp)[3] = "POS"
            write.table(tmp, str_replace_all(f, '.txt.gz', '_hg38.txt'), quote=F, row.names = F, sep = "\t")
            system(paste0('bgzip -f ', str_replace_all(f, '.txt.gz', '_hg38.txt')))
            system(paste0('tabix -S 1 -s 1 -b 3 -e -3 ', str_replace_all(f, '.txt.gz', '_hg38.txt.gz')))
        }
    }
}
#########

# tabix eqtls (data is in /Users/nicco/Documents/GitHub/snpXplorer/bin/snpxplorer_v4/data/databases/eqtls_snpxplorer)
  liftOver_data <- function(chrom, start, end, p = NULL, from){
    if (is.null(p)){ df <- data.frame(chr=paste("chr", chrom, sep=""), start=start, end=end) } else { df <- data.frame(chr=paste("chr", chrom, sep=""), start=start, end=end, p = p) }
    cat('--- adding tag to SNPs\n')
    df = df[order(df$start),]
    df$group = seq(1, nrow(df))
    # change to GR class object
    cat('--- creating GR object\n')
    gr <- makeGRangesFromDataFrame(df)
    # set chain file
    if (from == 'hg19') { chain <- import.chain("../hg19ToHg38.over.chain") } else { chain <- import.chain("../hg38ToHg19.over.chain") }
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

eqtls_files = system('ls chr*', intern = T)
for (f in eqtls_files[2:23]){
    cat(paste('** working on ', f, '\n'))
    d = data.table::fread(f, h=T, stringsAsFactors = F)
    colnames(d) = c('pos', 'a1', 'a2', 'tissue', 'ensg', 'effect', 'p')
    d$ensg = stringr::str_split_fixed(d$ensg, '\\.', 2)[, 1]
    # add genename
    d_annot = merge(d, mapping.ensembl, by.x = 'ensg', by.y = 'ensembl', all.x = T)
    # then liftover
    chrom = stringr::str_split_fixed(f, '_', 2)[, 1]; chrom = stringr::str_replace_all(chrom, 'chr', '')
    d_annot_lifted = liftOver_data(chrom = rep(chrom, nrow(d_annot)), start = d_annot$pos, end = d_annot$pos + 1, p = NULL, from = 'hg38')
    d_annot_lifted_info = merge(d_annot, d_annot_lifted, by.x = 'pos', by.y = 'start_hg38')
    d_annot_lifted_info$chr = chrom
    d_annot_lifted_info = d_annot_lifted_info[order(d_annot_lifted_info$pos),]
    # write and manage hg38
    fout = stringr::str_replace_all(f, '.txt.gz', '_hg38.txt')
    write.table(d_annot_lifted_info, fout, quote=F, row.names = F, sep = '\t')
    system(paste0('bgzip ', fout))
    system(paste0('tabix -S 1 -s 10 -b 1 -e -1 ', fout, '.gz'))
    # write and manage hg19
    fout = stringr::str_replace_all(f, '.txt.gz', '_hg37.txt')
    d_annot_lifted_info$pos = NULL
    colnames(d_annot_lifted_info)[8] = 'pos'
    d_annot_lifted_info = d_annot_lifted_info[order(d_annot_lifted_info$pos),]
    write.table(d_annot_lifted_info, fout, quote=F, row.names = F, sep = '\t')
    system(paste0('bgzip ', fout))
    system(paste0('tabix -S 1 -s 9 -b 8 -e -8 ', fout, '.gz'))
}

# tabix sqtls (data is in /Users/nicco/Documents/GitHub/snpXplorer/bin/snpxplorer_v4/data/databases/summary_sqtls)
sqtls_files = system('ls chr*', intern = T)
for (f in sqtls_files[2:22]){
    cat(paste('** working on ', f, '\n'))
    d = data.table::fread(f, h=T, stringsAsFactors = F)
    colnames(d) = c('pos', 'a1', 'a2', 'ensg', 'effect', 'p', 'tissue')
    d$ensg = stringr::str_split_fixed(d$ensg, '\\.', 2)[, 1]
    # add genename
    d_annot = merge(d, mapping.ensembl, by.x = 'ensg', by.y = 'ensembl', all.x = T)
    # then liftover
    chrom = stringr::str_split_fixed(f, '_', 2)[, 1]; chrom = stringr::str_replace_all(chrom, 'chr', '')
    d_annot_lifted = liftOver_data(chrom = rep(chrom, nrow(d_annot)), start = d_annot$pos, end = d_annot$pos + 1, p = NULL, from = 'hg38')
    d_annot_lifted_info = merge(d_annot, d_annot_lifted, by.x = 'pos', by.y = 'start_hg38')
    d_annot_lifted_info$chr = chrom
    d_annot_lifted_info = d_annot_lifted_info[order(d_annot_lifted_info$pos),]
    # write and manage hg38
    fout = stringr::str_replace_all(f, '.txt.gz', '_hg38.txt')
    write.table(d_annot_lifted_info, fout, quote=F, row.names = F, sep = '\t')
    system(paste0('bgzip ', fout))
    system(paste0('tabix -S 1 -s 10 -b 1 -e -1 ', fout, '.gz'))
    # write and manage hg19
    fout = stringr::str_replace_all(f, '.txt.gz', '_hg37.txt')
    d_annot_lifted_info$pos = NULL
    colnames(d_annot_lifted_info)[8] = 'pos'
    d_annot_lifted_info = d_annot_lifted_info[order(d_annot_lifted_info$pos),]
    write.table(d_annot_lifted_info, fout, quote=F, row.names = F, sep = '\t')
    system(paste0('bgzip ', fout))
    system(paste0('tabix -S 1 -s 9 -b 8 -e -8 ', fout, '.gz'))
}

#########
load('data/databases/snps_info/chrAll_snps_info.RData')
for (chr in seq(1, 22)){
    snps_info_all[[chr]]$CHR = chr
}
all_chroms = rbindlist(snps_info_all)
all_chroms = all_chroms[order(all_chroms$CHR, all_chroms$POS), ]
all_chroms_hg38 = all_chroms[order(all_chroms$CHR, all_chroms$POS_HG38), ]
write.table(all_chroms, 'data/databases/snps_info/chrAll_snps_info_hg19.txt', quote=F, row.names = F, sep = "\t")
write.table(all_chroms_hg38, 'data/databases/snps_info/chrAll_snps_info_hg38.txt', quote=F, row.names = F, sep = "\t")