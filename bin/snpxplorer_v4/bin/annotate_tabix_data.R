# libraries
library(data.table)
library(stringr)

# tabix all files
path = '/root/snpXplorer/data'
all_files = system(paste0('ls -d ', path, '/*/'), intern = T)
for (folder in all_files){
    if (folder != '/root/snpXplorer/data/databases/' && folder != '/root/snpXplorer/data/TABIX/'){
        filelist = system(paste0('ls ', folder, 'chr*'), intern = T)
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

load('data/databases/snps_info/chrAll_snps_info.RData')
for (chr in seq(1, 22)){
    snps_info_all[[chr]]$CHR = chr
}
all_chroms = rbindlist(snps_info_all)
all_chroms = all_chroms[order(all_chroms$CHR, all_chroms$POS), ]
all_chroms_hg38 = all_chroms[order(all_chroms$CHR, all_chroms$POS_HG38), ]
write.table(all_chroms, 'data/databases/snps_info/chrAll_snps_info_hg19.txt', quote=F, row.names = F, sep = "\t")
write.table(all_chroms_hg38, 'data/databases/snps_info/chrAll_snps_info_hg38.txt', quote=F, row.names = F, sep = "\t")