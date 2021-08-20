#################################################
## print tool message
cat("\n########################################\n")
cat("########## This is AnnotateMe ##########\n")
cat("# Tool to perform variant-gene mapping #\n")
cat("# and gene-set pathway enrichment.     #\n")
cat("# Hope you find the tool useful.       #\n")
cat("# Please report any issue or bug to:   #\n")
cat("# n.tesi@amsterdamumc.nl               #\n")
cat("########################################\n\n")

# CLEAN ENVIRONMENT
rm(list=ls())

# LIBRARIES
cat("## Loading Packages\n")
suppressPackageStartupMessages({
  library(data.table)
  library(htmlwidgets)
  library(NbClust)
  library(tidytext)
  library(stringr)
  library(webshot)
  library(viridis)
  library(dendextend)
  library(lme4)
  library(wordcloud2)
  library(dynamicTreeCut)
  library(ggsci)
  library(RColorBrewer)
  library(gprofiler2)
  library(rvest)
  library(GOSemSim)
  library(GO.db)
  library(org.Hs.eg.db)
  library(pheatmap)
  library(circlize)
  library(dplyr)
  library(plotrix)
  library(ggplot2)
  library(tibble)
  library(devtools)
  library(treemap)
  library(basicPlotteR)
  library(gwascat)
  library(GenomicRanges)
  library(rtracklayer)
  library(Homo.sapiens)
  library(BiocGenerics)
  library(liftOver)
  library(clusterProfiler)
  library(rjson)
  library(parallel)
  library(rsnps)
  #library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
  library(GenomicRanges)
  library(AllelicImbalance)
})
args = commandArgs(trailingOnly=TRUE)

# LOAD FUNCTIONS
## read main snp file and do the necessary adjustments
MAIN = "/root/snpXplorer/AnnotateMe/"
source(paste0(MAIN, "BIN/20210819_annotateMe_functions.R"))

# MAIN
## Read arguments
fname <- args[1]
#fname <- "/Users/home/Desktop/SNPbrowser/gitHub_version/snpXplorer_v2/RESULTS_65963/annotateMe_input_47224.txt"
ftype <- args[2]
#ftype <- 3
username <- args[3]
# gene-sets for enrichment analysis
analysis_mode = unlist(strsplit(as.character(args[4]), ","))
# tissues of interest
interesting_tissues = unlist(strsplit(as.character(args[5]), ","))
# reference genome (in case input is not rsid)
ref_version = as.character(args[6])

## create folder for results -- add a random number
random_num <- sample(x = seq(1, 100000), size = 1, replace = F)
cmd = paste("mkdir /root/snpXplorer/snpXplorer_v2/RESULTS_", random_num, sep="")
system(cmd)

## send email myself to notify that someone made a request
cmd_mail <- paste("sendEmail -f snpxplorer@gmail.com -t n.tesi@amsterdamumc.nl -u 'AnnotateMe request sent' -m 'Hello, \n a request to AnnotateMe was just sent from ", username, ". \n The following settings were requested: \n input --> ", fname, "\n input_type --> ", ftype, "\n analysis_mode --> ", analysis_mode, "\n interest_tissue --> ", interesting_tissues, "\n ref_version --> ", ref_version, "\n output_folder --> ", random_num, "\n \n AnnotateMe' -S /usr/sbin/sendmail")
#cmd_mail <- paste("sendEmail -f snpXplorer@gmail.com -t n.tesi@amsterdamumc.nl -u 'AnnotateMe request sent' -m 'Hello, \n a request to AnnotateMe was just sent from ", username, ". \n \n AnnotateMe' -s smtp.gmail.com:25 -xu snpXplorer@gmail.com -xp snpXplorer22101991!", sep="")
system(cmd_mail)

## also put input file with the list of snps in the folder
cmd = paste("cp /root/snpXplorer/snpXplorer_v2/", fname, " /root/snpXplorer/snpXplorer_v2/RESULTS_", random_num, "/", sep="")
system(cmd)

## read input data
cat("## Reading SNPs and positions\n")
data <- readSNPs(fname, ftype, MAIN, ref_version)
# manual settings for missing data
# if (ftype == 3){
#   data <- as.data.frame(data)
#   data[which(data$ID == "rs143080277"), c("chr", "pos", "ALT_FREQS", "locus")] <- c(2, 106366056, 0.005, "2:106366056")
#   data[which(data$ID == "rs3822030"), c("chr", "pos", "ALT_FREQS", "locus")] <- c(4, 987343, 0.47, "4:987343")
#   data[which(data$ID == "rs75932628"), c("chr", "pos", "ALT_FREQS", "locus")] <- c(6, 41129252, 0.001, "6:41129252")
#   data[which(data$ID == "rs60755019"), c("chr", "pos", "ALT_FREQS", "locus")] <- c(6, 41149008, 0.01, "6:41149008")
#   data[which(data$ID == "rs3919533"), c("chr", "pos", "ALT_FREQS", "locus")] <- c(12, 7162801, 0.154, "12:7162801")
#   data[which(data$ID == "rs616338"), c("chr", "pos",  "ALT_FREQS", "locus")] <- c(17, 47297297, 0.006, "17:47297297")
#   data[which(data$ID == "rs141749679"), c("chr", "pos",  "ALT_FREQS", "locus")] <- c(1, 109888432, 0.003, "1:109888432")
# }

## separate snps to annotate from the NAs
miss <- data[which(is.na(data$chr)),]
data <- data[!is.na(data$chr),]

## find all variants in LD with the input variants: this will help with the annotation
ld_info = data.frame(chr=NULL, pos=NULL)
# if (ftype == 3){
#   try(ld_info <- lapply(1:nrow(data), findLD, data=data), silent = T)
#   try(ld_info <- rbindlist(ld_info), silent = T)
# } else {
#   ld_info <- data.frame(chr=NULL, pos=NULL)
# }
## for AD variants, one is missing, add it manually
# tmp_df <- data.frame(chr = 17, pos = 44353222, ID = "rs2732703", ALT_FREQS = 0.1365, locus = "17:44353222")
# data <- rbind(data, tmp_df)

## read additional needed files: gtex, ensemble-gene_name mapping and gene positions
cat("## Loading all genes and gtex\n")
load("/root/snpXplorer/AnnotateMe/INPUTS_OTHER/annotationFiles.RData")
genes <- fread(paste(MAIN, "INPUTS_OTHER/NCBI37.3.gene.loc", sep=""), h=F, stringsAsFactors=F)
colnames(genes) <- c("gene_id", "chr", "start_tss", "stop_tss", "strand", "gene_name")
gtex <- fread(paste(MAIN, "INPUTS_OTHER/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz", sep=""))
mapping.ens <- fread(paste(MAIN, "INPUTS_OTHER/Ensemble_to_GeneName.txt", sep=""))
colnames(mapping.ens) <- c("ensemble", "gene")

## manage cadd files
cat("## Loading CADD scores\n")
#cadd <- mclapply(unique(data$chr), readcadd, data=data, MAIN=MAIN, ld_info=ld_info, mc.cores=1)
#cadd <- rbindlist(cadd)
# new cadd v6
cadd <- readcadd_v2(data, MAIN, ld_info, ref_genome)

## in case there were missings, add them back
#cadd$locus <- paste(cadd$"#Chrom", cadd$Pos, sep=":")
#miss_cadd <- data[which(!(data$locus %in% cadd$locus)),]
miss_cadd <- data[which(!(data$locus %in% cadd$locus)),]

## main function for functional annotation
cat("## Start annotation\n")
res <- AnnotateMe(data, genes, gtex, mapping.ens, ftype, MAIN, ld_info, cadd)
annot <- res[[1]]
geneList <- res[[2]]
# also nice to output the SVs close to the snps
all_sv = rbindlist(mclapply(1:nrow(annot), function_findSVs, annot = annot, all_str_hg38 = all_str_hg38, mc.cores=1))
write.table(all_sv, paste0("RESULTS_", random_num, "/SNP_and_SV_overlap.txt"), quote=F, row.names=F, sep="\t")
cat("## Annotation is done. Now making plots and functional enrichment analysis.\n")

# checked until here
## need now to check the genes associated with the variants: are real genes?
## try to do a round of gene-set analysis to see the mismapped genes
if (length(unique(geneList)) >1){
  # this maybe not necessary anymore as the function for gene-set enrichment has changed
  res_polished <- PolishAnnotation(annot, geneList)
  annot <- res_polished[[1]]
  geneList <- res_polished[[2]]
  
  ## save annotations in a folder called RESULTS
  write.table(annot, paste("RESULTS_", random_num, "/snp_annotation.txt", sep=""), quote=F, row.names=F, sep="\t")
  write.table(geneList, paste("RESULTS_", random_num, "/snp_annotation_geneList.txt", sep=""), quote=F, row.names=F, col.names=F)
  
  ## plot mapping characteristics
  annot$source_finalGenes[which(annot$snp_conseq %in% c("NON_SYNONYMOUS", "SYNONYMOUS"))] = "coding"
  annot$ALT_FREQS = as.numeric(annot$ALT_FREQS)
  pdf(paste("RESULTS_", random_num, "/snp_gene_mapping.pdf", sep=""), height=12.35, width=10)
  genesPerSNP <- plotMapping(annot, MAIN)
  invisible(dev.off())
  
  ## annotate using GWAS catalog
  cat("## Lookup in GWAS catalog\n")
  gwas.cat.res <- GWAScat(annot, geneList, MAIN)
  annot.perSNP <- gwas.cat.res[[1]]
  annot.perGene <- gwas.cat.res[[2]]
  
  ## functional annotation based on sampling
  ## each iteration, sample 1 gene from the pool of genes associated with a variant and do functional annotation with it
  ## at the end, do the mean of the pvalues
  cat("## Functional enrichment analysis: sampling N=1000 times\n")
  n.sampl <- 500
  # clusterProfiler is very slow and does not work with the sampling system
  # for clusterProfiler, seems that i have to run once without multiprocessing before running the whole sampling
  #tmp <- enrichGO(gene = geneList, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
  # first, sample gene-sets
  cat("  Deriving gene-set for the sampling\n")
  gene_sampling_dsets <- mclapply(1:n.sampl, samplingGset, mapping=annot, mc.cores=1)
  # finally gene-set analysis for all sampling dsets
  cat("  Gene-set overlap analysis for the sampling\n")
  # first, let's define the databases
  dbs = c()
  for (i in analysis_mode){
    if (i == "default"){
      #     dbs <- c(dbs, "GO_Biological_Process_2018")
      dbs <- c(dbs, "GO:BP")
    } else if (i == "KEGG"){
      #     dbs <- c(dbs, "KEGG_2019_Human")
      dbs <- c(dbs, "KEGG")
    } else if (i == "wiki"){
      #     dbs <- c(dbs, "WikiPathways_2019_Human")
      dbs <- c(dbs, "WP")
    } else if (i == "Reactome"){
      #     dbs <- c(dbs, "Reactome_2016")
      dbs <- c(dbs, "REAC")
    }
  }
  enrich.res <- mclapply(1:n.sampl, overlapAnalysis, gene_list=gene_sampling_dsets, source_gset = dbs, mc.cores=1)
  # save enrichment results for troubleshooting
  #save(enrich.res, file = paste("RESULTS_", random_num, "/enrichment_results.RData", sep=""))
  #enrich.res <- mclapply(1:n.sampl, newGeneSetEnrichment, gene_list=gene_sampling_dsets, mc.cores=1)
  #library(enrichR)
  #enrich.res <- mclapply(1:n.sampl, overlapAnalysis_enrichR, gene_list=gene_sampling_dsets, dbs = dbs, mc.cores=1)
  cat("  Merging and cleaning results\n")
  sampling.res <- mergeSampling(enrich.res)
  #sampling.res <- mergeSampling_enrichR(enrich.res, dbs)
  # in case i use gost (check for genes-pathway relationships) from here i need to change
  # first save all enrichment results
  #for (i in 1:length(sampling.res)){
  #write.table(sampling.res[[2]], file = paste0("RESULTS_", random_num, "/geneSet_enrichment_results.txt"), quote=F, row.names=F, sep="\t")
  #}
  # for other enrichment sources like kegg, reactome etc, need to make a plot
  if (length(dbs) >1){
    plotOtherEnrichments(sampling.res, dbs, random_num, type="gost")
  } else if (dbs[1] != "GO:BP"){
    plotOtherEnrichments(sampling.res, dbs, random_num, type = "gost")
  }
  # need to see whether GO was selected from the user --> for revigo things
  if (length(grep("GO", dbs)) >0){
    # get go data
    #go_data = sampling.res[[grep("GO", dbs)]]
    go_data = sampling.res[[2]]
    # extract term id
    # go_data$term_id = NA
    # for (i in 1:nrow(go_data)){
    #   go_data$term_id[i] <- unlist(strsplit(go_data$term[i], "[()]"))[length(unlist(strsplit(go_data$term[i], "[()]")))]
    # }
    colnames(go_data) <- c("term_name", "term_id", "avg_p", "log10p")
    # write go_data
    write.table(go_data, file = paste0("RESULTS_", random_num, "/geneSet_enrichment_results_and_clusters.txt"), quote=F, row.names = F, sep = "\t")
    # add source information to select only go
    go_data$source_gset = str_split_fixed(go_data$term_id, ":", 2)[, 1]
    go_data <- go_data[which(go_data$source_gset == "GO"),]
    # last thing is to run REVIGO and make the plot
    revigo_res <- revigo(avg_pvalues = go_data, thr = 0.01)
    # maybe nice to also do the alternative to revigo -- manual
    semsim_results = alternative_revigo_results(MAIN, random_num, go_data)
    functional_clusters = semsim_results[[1]]
    lin_matrix <- semsim_results[[2]]
    n_clust = semsim_results[[3]]

    if (length(functional_clusters) == 1){
      write.table(go_data, file = paste0("RESULTS_", random_num, "/geneSet_enrichment_results_and_clusters.txt"), quote=F, row.names=F, sep="\t")
    } else {
      # save the whole gene-set enrichment analysis
      tmp = functional_clusters
      tmp$term_name = NULL
      go_data = merge(go_data, tmp, by.x = "term_id", by.y = "term", all.x = T)
      go_data = go_data[order(go_data$avg_p),]
      go_data$cluster = NULL
      write.table(go_data, file = paste0("RESULTS_", random_num, "/geneSet_enrichment_results_and_clusters.txt"), quote=F, row.names=F, sep="\t")
      
      # let's try to make some wordcloud images
      if (!is.na(functional_clusters)){
        mostFreq_words <- CountFrequency_words(functional_clusters, n_clust)
      }
    }
  }
  cat("\n## Analysis done. Hope results make sense :)\n")
  
  # before sending everytihng, also add the file description in the folder
  system(paste0("cp /root/snpXplorer/snpXplorer_v2/www/snpXplorer_output_description.pdf RESULTS_", random_num, "/"))
  # finally compress result folder and send it
  cmd_compress <- paste("tar -czf AnnotateMe_results_", random_num, ".tar.gz RESULTS_", random_num, "/", sep="")
  system(cmd_compress)
  #cmd_mail <- paste("sendEmail -f snpXplorer@gmail.com -t ", username, " -u 'AnnotateMe results' -m 'Dear user, \n thanks so much for using snpXplorer and AnnotateMe. \n We hope you find the tool useful. \n AnnotateMe team.' -a 'AnnotateMe_results_", random_num, ".tar.gz' -cc n.tesi@amsterdamumc.nl -s smtp.gmail.com:25 -xu snpXplorer@gmail.com -xp snpXplorer22101991!", sep="")
  cmd_mail <- paste0("sendEmail -f snpxplorer@gmail.com -t ", username, " -u 'AnnotateMe results' -m 'Dear user, \n thanks so much for using snpXplorer and AnnotateMe. \n We hope you find the tool useful. \n AnnotateMe team.' -a 'AnnotateMe_results_", random_num, ".tar.gz' -cc n.tesi@amsterdamumc.nl -S /usr/sbin/sendmail")
  system(cmd_mail)
  # finally also remove input annotateME list of snps
  cmd = paste("rm ", fname, sep="")
  system(cmd)
  cmd = paste("rm -rf RESULTS_", random_num, "/", sep="")
  system(cmd)
  
} else {
  # no enough genes to polish and do gene-set analysis
  
  ## save annotations in a folder called RESULTS
  write.table(annot, paste("RESULTS_", random_num, "/snp_annotation.txt", sep=""), quote=F, row.names=F, sep="\t")
  write.table(geneList, paste("RESULTS_", random_num, "/snp_annotation_geneList.txt", sep=""), quote=F, row.names=F, col.names=F)
  
  ## plot mapping characteristics
  pdf(paste("RESULTS_", random_num, "/snp_gene_mapping.pdf", sep=""), height=12.35, width=10)
  genesPerSNP <- plotMapping(annot, MAIN)
  invisible(dev.off())
  
  ## annotate using GWAS catalog
  cat("## Lookup in GWAS catalog\n")
  gwas.cat.res <- GWAScat(annot, geneList, MAIN)
  annot.perSNP <- gwas.cat.res[[1]]
  annot.perGene <- gwas.cat.res[[2]]
  
  # finally compress result folder and send it
  cmd_compress <- paste("tar -czf AnnotateMe_results_", random_num, ".tar.gz RESULTS_", random_num, "/", sep="")
  system(cmd_compress)
  #cmd_mail <- paste("sendEmail -f snpXplorer@gmail.com -t ", username, " -u 'AnnotateMe results' -m 'Dear user, \n thanks so much for using snpXplorer and AnnotateMe. \n Unfortunately, due to a low number of genes found in your SNPs, we could not perform gene-set enrichment analysis. \n We hope you find the tool useful. \n AnnotateMe team.' -a 'AnnotateMe_results_", random_num, ".tar.gz' -cc n.tesi@amsterdamumc.nl -s smtp.gmail.com:25 -xu snpXplorer@gmail.com -xp snpXplorer22101991!", sep="")
  cmd_mail <- paste0("sendEmail -f snpxplorer@gmail.com -t ", username, " -u 'AnnotateMe results' -m 'Dear user, \n thanks so much for using snpXplorer and AnnotateMe. \n Unfortunately, due to a low number of genes found in your SNPs, we could not perform gene-set enrichment analysis. \n We hope you find the tool useful. \n AnnotateMe team.' -a 'AnnotateMe_results_", random_num, ".tar.gz' -cc n.tesi@amsterdamumc.nl -S /usr/sbin/sendmail")
  system(cmd_mail)
  # finally also remove input annotateME list of snps
  cmd = paste("rm ", fname, sep="")
  system(cmd)
  cmd = paste("rm -rf RESULTS_", random_num, "/", sep="")
  system(cmd)
  
}

cat("###########################################")

##################################################
##################################################
# PART THAT IS USED ONE TIME ONLY

## this part parses the GWAS catalog and create a dataframe of gene~traits -- it is saved as additional file
# gwas.sb <- gwas[, c("REPORTED GENE(S)", "MAPPED_GENE", "MAPPED_TRAIT", "FIRST AUTHOR", "JOURNAL")]
# gwas.sb <- as.matrix(gwas.sb)
# intervals <- ceiling(seq(nrow(gwas.sb)/10, nrow(gwas.sb), nrow(gwas.sb)/10))
# tmp_f <- function(i, gwas.sb, intervals){
# update
# if (i %in% intervals){ print(paste("  Performed", grep(i, intervals)*10, "% of the lines.."), sep="") }
# isolate line
# l <- gwas.sb[i, ]
# all.g <- c(unlist(strsplit(l[1], ",")), unlist(strsplit(l[2], ",")))
# all.g <- all.g[!duplicated(all.g)]
# if (TRUE %in% table(all.g %in% geneList)){
#    tmp.mt <- data.frame(gene = all.g, trait = rep(l[3], length(all.g)), study = rep(l[4], length(all.g)), journal = rep(l[5], length(all.g)))
#    return(tmp.mt)
# } else {
#    tmp.mt <- data.frame(gene = NA, trait = NA, study = NA, journal = NA)
#    return(tmp.mt)
#}
#}
# res_genes <- mclapply(1:nrow(gwas.sb), tmp_f, gwas.sb=gwas.sb, intervals=intervals, mc.cores=1)     # annotate every row (~200k rows)
# all.genes <- as.data.frame(rbindlist(res_genes))       # merge results
# all.genes <- all.genes[!is.na(all.genes), ]       # exclude NAs
# write.table(all.genes, "INPUTS_OTHER/Gwas_catalog_Gene_Traits.txt", quote=F, row.names=F, sep="\t")

###################################################


