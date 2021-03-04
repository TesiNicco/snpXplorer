##################################################
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
library(dendextend)
library(parallel)
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
})
args = commandArgs(trailingOnly=TRUE)

# FUNCTIONS
## function to find variant consequences from CADD
SNP_conseq_CHROM_generic <- function(i, df, cadd, ld_info){  
    # make subset in the input variants according to chromosome
    sb <- df[which(df$chr == i),]

    # take subset from bigger file
    d <- cadd[which(cadd$"#Chrom" == i),]
    
    # also take ld informations
    ld <- ld_info[which(ld_info$CHR_A == i),]

    # in case of duplicates for the same variant, isolate them and take the one with non-missing GeneName field
    sb$snp_conseq <- NA
    sb$snp_conseq_gene <- NA
    for (x in 1:nrow(sb)){
        # find snps information as annotated in cadd
        all <- d[which(d$Pos == sb$pos[x]), ]
        
        # also get annotation information
        info <- paste0(unique(all$Consequence[!is.na(all$Consequence)]), collapse=",")
        g <- paste0(unique(all$GeneName[!is.na(all$GeneName)]), collapse=",")
        
        # check if culprit gene was available
        info_splt <- unlist(strsplit(info, ","))
        culprit <- unique(info_splt %in% c("NON_SYNONYMOUS", "SYNONYMOUS", "STOP_GAINED", "MISSENSE"))    #if the variant is not missense, synonymous, stop_gained
        if (TRUE %in% culprit){
          grp <- grep(TRUE, culprit)
          info <- paste0(unique(all$Consequence[grp], collapse=","))
          g <- paste0(unique(all$GeneName[grp], collapse=","))
          sb$snp_conseq[x] <- info
          sb$snp_conseq_gene[x] <- g
        } else {
          if (info != "") { sb$snp_conseq[x] <- info } else { sb$snp_conseq[x] <- NA }
          if (g != "") { sb$snp_conseq_gene[x] <- g } else { sb$snp_conseq_gene[x] <- NA}

          ld_snps <- ld[which(ld$BP_A == unique(all$Pos)), "BP_B"]
          ld_snps <- ld_snps[which(ld_snps$BP_B != unique(all$Pos)),]
          # get infor for all the remaining variants
          cadd_ld <- d[which(d$Pos %in% ld_snps$BP_B),]
          cadd_ld_info <- rbind(cadd_ld[grep("SYNON", cadd_ld), ], cadd_ld[grep("MISS", cadd_ld), ], cadd_ld[grep("STOP_G", cadd_ld), ]) 
          if (nrow(cadd_ld_info) >0){
            sb$snp_conseq[x] <- paste0(unique(cadd_ld_info$Consequence[!is.na(cadd_ld_info$Consequence)]), collapse = ",")
            sb$snp_conseq_gene[x] <- paste0(unique(cadd_ld_info$GeneName[!is.na(cadd_ld_info$GeneName)]), collapse=",")
          }
        }
      }

    return(sb)
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

## grep all eqtl from my permutation dataset
GTEx_me_generic <- function(mapping, gtex, mapping.ens){
    # need to lift over as the v8 from gtex is wrt build 38
    # create data frame with chromosome and positions
    df <- data.frame(chr=paste("chr", mapping$chr, sep=""), start=mapping$pos, end=as.numeric(mapping$pos)+1)
    # change to GR class object
    gr <- makeGRangesFromDataFrame(df)
    # set chain file
    chain <- import.chain(paste(MAIN, "BIN/hg19ToHg38.over.chain", sep=""))
    # change coordinates
    gr_hg38 <- liftOver(gr, chain)
    # back to dataframe and clean it
    df_hg38 <- as.data.frame(gr_hg38)
    sb_df <- data.frame(chr_hg38 = as.character(df_hg38$seqnames), pos_hg38 = as.character(df_hg38$start), locus = paste(mapping$chr, df$start, sep=":"), code_gtex = paste(as.character(df$chr), "_", as.character(df_hg38$start), "_", sep=""))
    # merge with mapping using locus column
    mapping <- merge(mapping, sb_df, by="locus")

    # adjust gtex data for merging
    tmp <- str_split_fixed(gtex$variant_id, "_", 5)
    gtex$newcode <- paste(tmp[, 1], "_", tmp[, 2], "_", sep="")

    mapping$eqtl_blood <- NA
    gtex.sb <- gtex[which(gtex$newcode %in% mapping$code_gtex), ]

    for (j in 1:nrow(mapping)){
        if (is.na(mapping$coding_snp[j])){ 
            #look for eqtl
            eqtl <- gtex.sb[which(gtex.sb$newcode == mapping$code_gtex[j]), ]
            if (nrow(eqtl) >0){ #if there are eqtl, save results and stop here
                #grep gene name
                id1 <- str_split_fixed(eqtl$gene_id, "\\.", 2)
                gene_name <- mapping.ens[which(mapping.ens$ensemble %in% id1[, 1]), "gene"]
                if (nrow(gene_name) >1){ 
                    mapping$eqtl_blood[j] <- paste0(gene_name$gene, collapse=",")
                } else if (nrow(gene_name) ==1) { 
                    mapping$eqtl_blood[j] <- as.character(gene_name[, 1]) 
                } 
            } else {
                mapping$eqtl_blood[j] <- NA
            }
        }
    }
    return(mapping)
}

## function to run positional mapping for the variants with no annotation from consequences and eqtl
positionalAnn_generic <- function(mapping, genes){
    mapping$positional_mapping <- NA        #add column to be filled
    #tmp.exc <- mapping[!(is.na(mapping$coding_snp) & is.na(mapping$eqtl_blood)), ]      #identify snps to not look at
    #mapping <- mapping[is.na(mapping$coding_snp) & is.na(mapping$eqtl_blood), ]     #make clean dataset
    tmp.exc <- mapping[!(is.na(mapping$coding_snp)), ]      #identify snps to not look at
    mapping <- mapping[is.na(mapping$coding_snp), ]     #make clean dataset
    mapping$pos <- as.numeric(mapping$pos)
    for (j in 1:nrow(mapping)){
        w = 500000   #set window threshold -- this value will be doubled (upstream and downstream)
        sb.gene <- genes[which(genes$chr == mapping$chr[j]),]       #restrict to chromosome of use
        sb.gene <- sb.gene[which((abs(mapping$pos[j] - sb.gene$start_tss <= w)) | (abs(mapping$pos[j] - sb.gene$stop_tss <= w))), ]       #restrict within 1MB region
        sb.gene$distance_start <- abs(mapping$pos[j] - sb.gene$start_tss)        #calculate distance from start tss
        sb.gene$distance_stop <- abs(mapping$pos[j] - sb.gene$stop_tss)      #calculate distance from stop tss
        sb.gene$distance_start[which(mapping$pos[j] >= sb.gene$start_tss & mapping$pos[j] <= sb.gene$stop_tss)] <- 0      #also assign 0 for intronic positions
        sb.gene$min_distance <- apply(X=sb.gene[, c(7, 8)], MARGIN=1, FUN=min)      #function to get the minimum distance between t-start-s and t-stop-site 
        sb.gene$code <- apply(X=sb.gene[, c(9)], MARGIN=1, FUN=function(x){return(round(x/(w/10), 0))})       #assign code (intron=0, <50kb=1, <100kb=2, etc)
        sb.gene <- sb.gene[order(sb.gene$code),]
        if (nrow(sb.gene) > 0){  #if there are hits, check overlap with gene indicated from cadd, save results and stop here
            #tmp <- c(sb.gene$gene_name[which(sb.gene$code == min(sb.gene$code))], unlist(strsplit(mapping$snp_conseq_gene[j], ",")))
            #tmp <- tmp[!duplicated(tmp) & !is.na(tmp)]
            #mapping$positional_mapping[j] <- paste0(tmp, collapse=",")
            # maybe including also cadd introduce too much variablility? try without here
            mapping$positional_mapping[j] <- paste0(sb.gene$gene_name[which(sb.gene$code == min(sb.gene$code))], collapse=",")
        } else {
            mapping$positional_mapping[j] <- mapping$snp_conseq_gene[j]
        }
    }
    mapping <- rbind(mapping, tmp.exc)      #re-join dataset
    return(mapping)
}

## function to clean file after mapping procedure and output the list of genes
getGeneList_mod_generic <- function(mapping){
    #put information about which source will be used to get genes: priority is: 1-coding_snp -- 2-eqtl -- 3-positional
    mapping$source_finalGenes <- NA
    mapping$geneList <- NA

    #assign source and gene list for coding variants
    mapping$source_finalGenes[which(!is.na(mapping$coding_snp))] <- "coding"
    mapping$geneList[which(!is.na(mapping$coding_snp))] <- mapping$snp_conseq_gene[which(!is.na(mapping$coding_snp))]
    
    # before asigning, need to exclude pseudogenes (grep a . in the gene name)
    pseudo <- mapping[grep("\\.", mapping$eqtl_blood),]
    if (nrow(pseudo) >1){
      for (k in 1:nrow(pseudo)){ 
        tmp <- unlist(strsplit(pseudo$eqtl_blood[k], ",")) 
        if (length(tmp) >1){
          todel <- grep("\\.", tmp)
          tmp <- tmp[-todel]
          pseudo$eqtl_blood[k] <- paste0(tmp, collapse=",")
        } else {
          pseudo$eqtl_blood[k] <- NA
        }
        mapping[which(mapping$locus == pseudo$locus[k]), ] <- pseudo[k, ]
      }
    }
    
    # assign source for eqtl hits
    mapping$source_finalGenes[which(!is.na(mapping$eqtl_blood))] <- "eqtl+cadd"
    for (l in 1:nrow(mapping)){
      if ((!is.na(mapping$source_finalGenes[l])) & mapping$source_finalGenes[l] == "eqtl+cadd"){
        if (is.na(mapping$snp_conseq_gene[l]) | (mapping$snp_conseq_gene[l] == mapping$eqtl_blood[l])){
          mapping$geneList[l] <- mapping$eqtl_blood[l]
        } else {
          tmp_eqtl_list <- unlist(strsplit(mapping$eqtl_blood[l], ","))
          tmp_cadd_list <- unlist(strsplit(mapping$snp_conseq_gene[l], ","))
          tmp_all_list <- c(tmp_eqtl_list, tmp_cadd_list)
          is_pseudo <- grep("\\.", tmp_all_list)
          if (length(is_pseudo) >0){
              tmp_all_list <- tmp_all_list[-is_pseudo]
          }
          tmp_full_list <- paste0(unique(tmp_all_list), collapse=",")
          mapping$geneList[l] <- tmp_full_list
        }
      }
    }

    #assign source for position hits
    mapping$source_finalGenes[which(is.na(mapping$coding_snp) & is.na(mapping$eqtl_blood))] <- "positional"
    #mapping$geneList[which(is.na(mapping$coding_snp) & is.na(mapping$eqtl_blood))] <- mapping$positional_mapping[which(is.na(mapping$coding_snp) & is.na(mapping$eqtl_blood))]
    # for these annotations, also consider the gene suggested by cadd
    tmp_posit <- mapping[which(mapping$source_finalGenes == "positional"),]
    for (k in 1:nrow(tmp_posit)){
      # identify cadd genes
      if (is.na(tmp_posit$snp_conseq_gene[k])){
        cadd_genes = NA
      } else {
        cadd_genes <- unlist(strsplit(tmp_posit$snp_conseq_gene[k], ","))
      }
      # identify positional genes
      pos_genes <- unlist(strsplit(tmp_posit$positional_mapping[k], ","))
      # check for pseudogenes and in case remove them
      grp.pseu <- grep("\\.", cadd_genes)
      if (length(grp.pseu) >0){ cadd_genes <- cadd_genes[-grp.pseu] }
      if ((length(cadd_genes) ==0) || (is.na(cadd_genes))){ 
        mapping$snp_conseq_gene[which(mapping$locus == tmp_posit$locus[k])] <- NA
        mapping$geneList[which(mapping$locus == tmp_posit$locus[k])] <- paste0(pos_genes, collapse = ",")
      } else {
        tmp_common <- intersect(cadd_genes, pos_genes)
        if (length(tmp_common) >0){
          mapping$geneList[which(mapping$locus == tmp_posit$locus[k])] <- paste0(tmp_common, collapse = ",")
        } else {
          mapping$geneList[which(mapping$locus == tmp_posit$locus[k])] <- paste0(c(pos_genes, cadd_genes), collapse = ",")
        }
      }
    }
    
    
    #assign finally those that did not map anywhere
    mapping$source_finalGenes[which(is.na(mapping$coding_snp) & is.na(mapping$eqtl_blood) & is.na(mapping$positional_mapping))] <- "missing"

    tmp1 <- mapping$geneList[which(!is.na(mapping$geneList))]
    tmp2 <- strsplit(tmp1, ",")
    all.genes <- c()
    for (i in 1:length(tmp2)) { all.genes <- c(all.genes, tmp2[[i]]) }
    all.genes <- all.genes[!is.na(all.genes)]
    all.genes <- all.genes[which(!(all.genes %in% c("NA", "character(0)")))]
    all.genes <- all.genes[!duplicated(all.genes)]
    all.genes <- str_split_fixed(all.genes, "\\.", 2) 
    all.genes <- all.genes[, 1]

    #merge results
    l = list(all.genes, mapping)
    
    return(l)
}

## function to read and store cadd files
readcadd <- function(i, data, MAIN, ld_info){
    # command for reading file
    fname <- paste(MAIN, "INPUTS_OTHER/CADD/chr", i, "_cadd.txt.gz", sep="")
    f <- fread(fname, h=T, showProgress=FALSE)

    # define list of all snps to look at in this chromosome
    snps_to_take <- unique(c(data$pos[which(data$chr == i)], ld_info$BP_A[which(ld_info$CHR_A == i)], ld_info$BP_B[which(ld_info$CHR_A == i)]))
    
    f <- f[which(f$Pos %in% snps_to_take), ]     #keep only same data as in reference data (ukb)
    cat(paste("  Done with chromosome", i, "\n"))
    return(f)
}

## function to add variant frequency from 1000Genome
addFreq <- function(i, data, MAIN){
    # subset on the chromosome
    sb <- data[which(data$chr == i), ]

    # read corresponding file
    ref <- fread(paste(MAIN, "INPUTS_OTHER/1000G_frequencies/chr", i, ".afreq.gz", sep=""), h=T, showProgress=FALSE)
    info <- ref[which(ref$POS %in% sb$pos), c("POS", "ID", "ALT_FREQS")]

    # add to file
    sb <- merge(sb, info, by.x="pos", by.y="POS", all.x=T)

    return(sb)
}

## main function to guide the functional annotation of variants
AnnotateMe <- function(data, genes, gtex, mapping.ens, ftype, MAIN, ld_info, cadd){
    # run functional annotation with cadd
    chroms <- unique(data$chr)
    cat("  Adding CADD annotation\n")
    out <- mclapply(chroms, SNP_conseq_CHROM_generic, df=data, cadd=cadd, ld_info=ld_info, mc.cores=1)
    # merge results together
    snps.conseq <- rbindlist(out)
    # manually change APOE4 cause for some reasons it is missing
    snps.conseq[which(snps.conseq$pos == 45411941), "snp_conseq"] <- "MISSENSE"
    snps.conseq[which(snps.conseq$pos == 45411941), "snp_conseq_gene"] <- "APOE"

    # need to add frequencies
    cat("  Adding variant frequencies\n")
    if (ftype %in% c(1, 2)){
        out <- mclapply(chroms, addFreq, data=snps.conseq, MAIN=MAIN, mc.cores=1)
        snps.conseq <- rbindlist(out)
    }
    
    # need to re-couple with the permutation sets now
    out.annot <- mergeInfo_generic(info.sb=snps.conseq)

    # now run the GTEx and positional annotations
    cat("  GTEx annotation for variants of interest\n")
    out.gtex <- GTEx_me_generic(mapping=out.annot, gtex=gtex, mapping.ens=mapping.ens)

    #for the remaining, need to run the positional mapping
    cat("  Positional annotation for variants of interest\n")
    full.annot <- positionalAnn_generic(mapping=out.gtex, genes=genes)

    #extract gene list 
    res.clean <- getGeneList_mod_generic(mapping=full.annot)
    geneList <- res.clean[[1]]
    final.annot <- res.clean[[2]]
    
    return(list(final.annot, geneList))
}

## function to annotate using GWAS catalog
GWAScat <- function(annot, geneList, MAIN){
    # read gwas catalog
    gwas <- fread(paste(MAIN, "INPUTS_OTHER/GWAS_catalog_20200310.txt.gz", sep=""), h=T, showProgress=FALSE)

    # do 2 ways of merging to maximize results: rsid and chr:pos
    gwas$LOCUS <- paste(gwas$CHR_ID, gwas$CHR_POS, sep=":")
    m1 <- gwas[which(gwas$SNPS %in% annot$ID), ]
    m2 <- gwas[which(gwas$LOCUS %in% annot$locus), ]
    all.m <- rbind(m1, m2)
    all.m <- all.m[, c("SNPS", "MAPPED_TRAIT", "FIRST AUTHOR", "JOURNAL")]
    all.m <- all.m[!duplicated(all.m),]

    # check number of matches
    if (nrow(all.m) >0){
        # separate data to be plotted
        traits <- as.data.frame(table(all.m$"MAPPED_TRAIT"))
        traits <- traits[order(-traits$Freq),]
        traits$id <- seq(1, nrow(traits))

        # run function for lolliplot
        pdf(paste("RESULTS_", random_num, "/gwas_cat_snps_overlap.pdf", sep=""), height=6, width=6)
        lolliPlot_snp(traits, all.m)
        invisible(dev.off())
    } else {
        cat("  !!! No matches for SNP in GWAS catalog\n")
        traits <- NA
    }

    # also check genes directly -- first read genes~trait dataframe
    all.genes <- fread(paste(MAIN, "INPUTS_OTHER/Gwas_catalog_Gene_Traits.txt", sep=""), h=T, sep="\t")
    
    # check overlap as a background
    overlapping.genes <- all.genes[which(all.genes$gene %in% geneList),]        # find the overlapping genes with my gene list
    overlapping.genes <- overlapping.genes[!duplicated(overlapping.genes), ]    # exclude duplicated rows
    print(overlapping.genes)
    # also here check number of matching genes
    if (nrow(overlapping.genes) >0){
        # sampling: each iteration sample 1 gene from the pool of genes associated with each variant and do the overlap
        # at the end, do the mean of the overlaps
        n.samp <- 500
        samp.res <- lapply(1:n.samp, sampleAnnot, overlapping.genes=overlapping.genes, annot=annot)       #sample 1000 times 1 gene per snp and do the overlap 1000 times
        options(warn=-1)
        tb <- Reduce(function(x, y) merge(x, y, by="Var1", all.x=T, all.y=T), samp.res)
        tb$SUM <- rowSums(tb[, 2:ncol(tb)], na.rm=T)
        tb <- tb[, c("Var1", "SUM")]
        tb$MEAN <- tb$SUM/n.samp
        tb <- tb[order(-tb$MEAN),]
        tb <- tb[!is.na(tb$Var1),]

        # finally lolliplot for genes
        pdf(paste("RESULTS_", random_num, "/gwas_cat_genes_overlap.pdf", sep=""), height=6, width=6)
        lolliPlot_gene(tb, overlapping.genes, geneList)
        invisible(dev.off())
        options(warn=0)
    } else {
        cat("  !!! No matches for Genes in GWAS catalog\n")
        tb <- NA
    }

    ls <- list(traits, tb)
    return(ls)
}

## function to make lolli-plot for genes
lolliPlot_gene <- function(tb, overlapping.genes, geneList){
    # will plot top 10
    sb.tr <- head(tb, 10)
    
    # initialize variable
    sb.tr$upd_name <- NA
    sb.tr$study <- NA

    # need to add some \n in the name description
    for (i in 1:nrow(sb.tr)){
        # count the characters
        tmp.string <- as.character(sb.tr$Var1[i])
        tmp.study <- unique(overlapping.genes[which(overlapping.genes$trait == tmp.string), "study"])
        tmp.study <- as.data.frame(str_split_fixed(tmp.study$study, " ", 2))
        sb.tr$study[i] <- nrow(tmp.study)

        # calculate in how many lines to distribute the name -- max is 23 characters (including spaces) per line
        n.lines <- ceiling(nchar(tmp.string)/25)
        
        # divide by space depending on n.lines
        tmp <- unlist(strsplit(as.character(sb.tr$Var1[i]), " "))
        if (n.lines >1){ 
            x <- split(tmp, ceiling(seq_along(tmp)/n.lines))
            for (j in 1:length(x)){ x[[j]] <- paste(x[[j]], collapse=" ") }
            # add name to data frame
            sb.tr$upd_name[i] <- paste(x, collapse="\n")
        } else {
            sb.tr$upd_name[i] <- paste(tmp, collapse=" ")
        }
    }

    # uppercase first letter of sentence
    sb.tr$upd_name <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", sb.tr$upd_name, perl = TRUE)

    # graphical parameters
    par(mar=c(5, 11, 4, 3))
    
    # basic empty plot
    plot(0, 0, pch=16, col="white", xlim=c(0, ceiling(max(sb.tr$MEAN))), ylim=c(1, 10), xlab="Mean Overlapping Genes", ylab="", cex.axis=1.25, cex.lab=1.50, yaxt="none", bty='n', xaxt='none')

    # axis
    axis(side=1, at=seq(0, ceiling(max(sb.tr$MEAN)), ceiling(max(sb.tr$MEAN)/5)), label=seq(0, ceiling(max(sb.tr$MEAN)), ceiling(max(sb.tr$MEAN)/5)), cex.lab=1.25, cex.axis=1.25)    
    
    # add grid
    for (i in seq(0, ceiling(max(sb.tr$MEAN)), ceiling(max(sb.tr$MEAN))/10)){ abline(v=i, lwd=0.4, col="grey80") }
    for (i in seq(1, 10)){ segments(x0=0, y0=i, x1=ceiling(max(sb.tr$MEAN)), y1=i, lwd=0.4, col="grey80") }
    
    # set up colors
    colz <- colorRampPalette(c("grey80", "orange", "deepskyblue3", "darkolivegreen3"))
    colorz <- colz(nrow(sb.tr))
    
    # add lollipops
    c <- 1
    for (i in seq(10, 1)){
        segments(x0=0, y0=i, x1=sb.tr$MEAN[c], y1=i, lwd=3, col=colorz[c])
        points(x=sb.tr$MEAN[c], y=i, pch=16, col=colorz[c], cex=3)
        text(x=-ceiling(max(sb.tr$MEAN))*0.05, y=i, labels=sb.tr$upd_name[c], cex=0.80, adj=1, xpd=T)
        c <- c+1
    }
}

## function to make lolli-plot for snps
lolliPlot_snp <- function(traits, all.m){
    # will plot top 10 -- isolate them here
    traits <- traits[order(-traits$Freq), ]
    sb.tr <- head(traits, 10)

    # initialize variable
    sb.tr$upd_name <- NA
    sb.tr$study <- NA

    # need to add some \n in the name description
    for (i in 1:nrow(sb.tr)){
        # count the characters
        tmp.string <- as.character(sb.tr$Var1[i])
        tmp.study <- unique(all.m[which(all.m$MAPPED_TRAIT == tmp.string), c("FIRST AUTHOR")])
        tmp.study <- as.data.frame(str_split_fixed(tmp.study$'FIRST AUTHOR', " ", 2))
        sb.tr$study[i] <- nrow(tmp.study)

        # calculate in how many lines to distribute the name -- max is 23 characters (including spaces) per line
        n.lines <- ceiling(nchar(tmp.string)/25)
        
        # divide by space depending on n.lines
        tmp <- unlist(strsplit(as.character(sb.tr$Var1[i]), " "))
        if (n.lines >1){ 
            x <- split(tmp, ceiling(seq_along(tmp)/n.lines))
            for (j in 1:length(x)){ x[[j]] <- paste(x[[j]], collapse=" ") }
            # add name to data frame
            sb.tr$upd_name[i] <- paste(x, collapse="\n")
        } else {
            sb.tr$upd_name[i] <- paste(tmp, collapse=" ")
        }
    }

    # uppercase first letter of sentence
    sb.tr$upd_name <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", sb.tr$upd_name, perl = TRUE)

    # graphical parameters
    par(mar=c(5, 11, 4, 3))
    # basic empty plot
    plot(0, 0, pch=16, col="white", xlim=c(0, ceiling(max(sb.tr$Freq))), ylim=c(1, 10), xlab="Overlapping SNPs", ylab="", cex.axis=1.25, cex.lab=1.50, yaxt="none", bty='n')
    # add grid
    for (i in seq(0, ceiling(max(sb.tr$Freq)), ceiling(max(sb.tr$Freq))/10)){ abline(v=i, lwd=0.4, col="grey80") }
    for (i in seq(1, 10)){ abline(h=i, lwd=0.4, col="grey80") }
    # set up colors
    colz <- colorRampPalette(c("red", "orange", "dark green", "blue"))
    colorz <- colz(nrow(sb.tr))
    # add lollipops
    c <- 1
    for (i in seq(10, 1)){
        segments(x0=0, y0=i, x1=sb.tr$Freq[c], y1=i, lwd=3, col=colorz[c])
        points(x=sb.tr$Freq[c], y=i, pch=16, col=colorz[c], cex=3)
        text(x=-ceiling(max(sb.tr$Freq))*0.05, y=i, labels=sb.tr$upd_name[c], cex=0.80, adj=1, xpd=T)
        c <- c+1
    }
}

## function to sample for multiple genes associated with each variant
sampleAnnot <- function(i, overlapping.genes, annot){
    all.genes <- strsplit(annot$geneList, ",")
    tmp_f <- function(x, all.genes){ tmp <- all.genes[[x]]; g <- sample(x=tmp, size=1); return(g)}
    gset <- unlist(lapply(1:length(all.genes), tmp_f, all.genes=all.genes))
    overl <- overlapping.genes[which(overlapping.genes$gene %in% gset),]
    # check how many matches
    if (nrow(overl) >0){
        df <- as.data.frame(table(overl$trait))     # make table of the trait variable
        df <- df[order(-df$Freq),]      # order the data frame according to occurrence
    } else {
        df <- data.frame(Var1=NA, Freq=NA)
    }
    return(df)
}

## function to plot mapping charcteristics
plotMapping <- function(mapping.after, MAIN){
    #set how many plots to make and the margins
    layout(matrix(c(1, 2, 3, 3, 4, 4, 4, 4), nrow = 4, ncol = 2, byrow = TRUE))
    par(mar=c(3, 3, 3, 3))

    #color palette
    colz <- pal_lancet(palette = "lanonc")(5)
    myPalette <- brewer.pal(3, "Set2")

    #plot 1 -- barplot of the annotation sources
    labs <- c()
    if (nrow(mapping.after[which(mapping.after$source_finalGenes == "coding"),]) >0){ lab.cod <- paste("Coding\n(N=", nrow(mapping.after[which(mapping.after$source_finalGenes == "coding"),]), ")", sep=""); labs <- c(labs, lab.cod) }
    if (nrow(mapping.after[which(mapping.after$source_finalGenes == "eqtl+cadd"),]) >0){ lab.eqtl <- paste("eQTL\n(N=", nrow(mapping.after[which(mapping.after$source_finalGenes == "eqtl+cadd"),]), ")", sep=""); labs <- c(labs, lab.eqtl) }
    if (nrow(mapping.after[which(mapping.after$source_finalGenes == "positional"),]) >0){ lab.pos <- paste("Position\n(N=", nrow(mapping.after[which(mapping.after$source_finalGenes == "positional"),]), ")", sep=""); labs <- c(labs, lab.pos) }
    pie(table(mapping.after$source_finalGenes), col=myPalette, labels=labs, lwd=1.5, cex=1.50)
    
    #plot 2 -- histogram of snp-gene number
    red.na <- mapping.after[which(is.na(mapping.after$geneList)), c("locus", "geneList")]
    red <- mapping.after[which(!is.na(mapping.after$geneList)), c("locus", "geneList")]
    #output data frame
    df <- as.data.frame(matrix(data=NA, nrow=nrow(mapping.after), ncol=2))
    colnames(df) <- c("locus", "n_genes")
    for (i in 1:nrow(red)){
        #split geneList column
        n.gens <- strsplit(as.character(red$geneList[i]), ",")
        df[i, ] <- c(as.character(red$locus[i]), length(n.gens[[1]]))
    }
    if (nrow(red.na) >0){ for (j in 1:nrow(red.na)){ df[(i+j),] <- c(as.character(red.na$locus[j]), 0) } }
    #get max
    df$n_genes <- as.numeric(df$n_genes)
    mx <- max(df$n_genes)
    #then the plot
    par(mar=c(5, 8, 5, 10))
    barplot(table(df$n_genes), xlab="Genes per variant", ylab="Frequency", main="", cex.axis=1.25, col=colz, cex.lab=1.5, cex.names=1.25)

    #plot 3 --  histogram of distribution per chromosome
    par(mar=c(4, 8, 3, 4))
    colorz <- colorRampPalette(c("red", "orange", "dark green", "blue"))
    colorz.chr <- colorz(22)
    tmp <- as.data.frame(table(mapping.after$chr))
    tmp$Var1 <- as.numeric(as.character(tmp$Var1))
    if (nrow(tmp) != 22){ for (i in 1:22){ if (!(i %in% tmp$Var1)){ t <- data.frame(Var1=i, Freq=0); tmp <- rbind(tmp, t) } } }
    tmp <- tmp[order(tmp$Var1),]
    barplot(tmp$Freq, col=colorz.chr, names=tmp$Var1, xlab="Chromosome", ylab="Genes per chromosome", cex.lab=1.50, main="", cex.axis=1.4, cex.main=2)

    #plot 4 -- circular plot integrated visualization
    # read chromosome length
    l <- read.table(paste(MAIN, "INPUTS_OTHER/chromosomes_length_hg19.txt", sep=""), h=F)
    tmp.df <- data.frame(V1="AX", V2=50000000)
    l <- rbind(tmp.df, l)
    
    # circle circumference is the sum of all chromosome lengths
    R <- sum(l$V2)/1000000
    
    # then derive radius
    r <- R/(2*pi)
    
    # find corresponding angles for the chromosomes
    ang <- (l$V2/1000000)/r
    
    # set colors
    colz <- colorRampPalette(c("yellow", "light green", "dark green", "light blue", "navy", "orange", "red"))
    colorz <- colz(22)
    col.bk <- rep(c("white", "grey90"), 11)
    
    # set pch
    annot$pch <- 21
    annot$pch[which(annot$source_finalGenes == "coding")] <- 23
    annot$pch[which(annot$source_finalGenes == "eqtl")] <- 24
  
    par(mar=c(4, 18, 4, 18))
    plot(0, 0, pch=16, col="white", xlim=c(-r, r), ylim=c(-r, r), bty='n', xlab="", ylab="", xaxt='none', yaxt="none")
    draw.circle(x = 0, y = 0, radius = r, border=NA)
    cum.angle <- pi/2-(ang[1]/2)
    # first loop to add sectors
    for (i in 1:23){ 
        # get degrees
        dg1 <- cum.angle*180/pi
        dg2 <- (cum.angle+ang[i])*180/pi
        # draw background and borders
        draw.sector(start.degree = dg1, end.degree = dg2, rou1 = r, rou2 = 0, clock.wise = FALSE, col=alpha(col.bk[i], 0.6), border = NA)
        cum.angle <- cum.angle + ang[i]
    }
  
    # add axis for maf
    a = r*0.85
    b = r*0.30
    d = r*0.575
    c = r*0.01
    segments(x0 = 0, y0 = b, x1 = 0, y1 = a, lwd=1.25, xpd=T)
    # add maf circles for reference
    draw.circle(x = 0, y = 0, radius = a-a*0.025, border="grey40", lty=2)
    draw.circle(x = 0, y = 0, radius = b-b*0.025, border="grey40", lty=2, nv = 1000)
    draw.circle(x = 0, y = 0, radius = d-d*0.025, border="grey40", lty=2, nv = 1000)
    # add axis ticks
    for (i in seq(0, 0.5, 0.1)){
        trx <- (b-a)*(i/0.5)+a
        segments(x0 = -c, y0 = trx, x1 = c, y1 = trx, lwd=1.25, xpd=T)
        # add text
        text(x = 0, y = trx, labels = i, pos=2, offset = 0.2, cex=0.75, font=2)
    }
    dg1 <- (pi/2-(ang[1]/2))*180/pi
    dg2 <- ((pi/2+(ang[1]/2)))*180/pi
    draw.sector(start.degree = dg1, end.degree = dg2, rou1 = r, rou2 = r-r*0.1, clock.wise = FALSE, col = NA)
    text(x = 0, y = r*0.95, labels = "MAF", font=2, adj=0.5, cex = 0.70)
  
    # second loop to add points
    cum.angle <- pi/2+(ang[1]/2)
    for (i in 1:22){
        # calculate points for the segment of the lolliplot -- point0 should be maf
        pp <- annot[which(annot$chr == i),]
        pp$ALT_FREQS[which(pp$ALT_FREQS >0.5)] <- 1-pp$ALT_FREQS[which(pp$ALT_FREQS >0.5)]
    
        # for the frequency, need to normalize it in [r*0.20-r*0.80] [b, a]
        pp$maf_norm <- (b-a)*(pp$ALT_FREQS/0.5)+a
    
        # get the angle
        pp$ang <- (pp$pos/1000000)/r + cum.angle
        pp$x <- pp$maf_norm*cos(pp$ang)
        pp$y <- pp$maf_norm*sin(pp$ang)
        segments(x0 = pp$x, y0 = pp$y, x1 = r*cos(pp$ang), y1 = r*sin(pp$ang), lwd=0.85, col="grey40")
        points(x = pp$x, y = pp$y, type = "p", pch=pp$pch, col="black", bg=colorz[i], cex=1.75)
    
        # draw sector finally
        # get degrees
        dg1 <- cum.angle*180/pi
        dg2 <- (cum.angle+ang[i+1])*180/pi
        draw.sector(start.degree = dg1, end.degree = dg2, rou1 = r, rou2 = r-r*0.1, clock.wise = FALSE, col=colorz[i])

        # add label
        rgm <- cum.angle+(ang[i+1]/2)
        text(x = r*cos(rgm)*0.95, y = r*sin(rgm)*0.95, labels = i, font=2, adj=0.5)
    
        # increment at the end
        cum.angle <- cum.angle + ang[i+1]
    }
    # finally the legend
    legend(x = -r/2, y = -r*1.05, xpd=T, legend = c("Coding", "eQTL", "Position"), pch=c(23, 24, 21), pt.bg="grey40", pt.cex = 1.50, cex=1.50, bty='n', ncol=3)

    return(df)
}

## function to perform overlap analysis over the sampling datasets
overlapAnalysis <- function(i, gene_list, source_gset){
    g <- gene_list[[i]]
    
    res <- suppressMessages(gost(g, organism = "hsapiens", ordered_query = FALSE, multi_query = FALSE, significant = FALSE, exclude_iea = TRUE,
        measure_underrepresentation = FALSE, evcodes = FALSE, user_threshold = 1, correction_method = "fdr",
        domain_scope = "annotated", custom_bg = NULL, numeric_ns = "", sources = source_gset))
    return(res$result)
}

## function to merge sampling procedure for the functional annotation -- clusterProfiler
mergeSampling_clusProf <- function(enrich.res){
  # extract terms common to all datasets
  all <- rbindlist(tr)
  counts <- as.data.frame(table(all$ID))
  common <- as.character(counts[which(counts$Freq == length(tr)), "Var1"])
  
  # function to do mean
  tmp_f <- function(i, all){
    # subset to term of interest
    sb <- all[which(all$ID == i),]
    
    # order by number of characters of the gene IDs
    sb <- sb[order(-nchar(sb$geneID)),]
    
    # put in separate dataframe
    df <- data.frame(ID=i, Description=unique(sb$Description), p_adj = mean(sb$p.adjust), q_value = mean(sb$qvalue), geneID=sb$geneID[1])
    
    return(df)
  }
  
  # run in multiprocessing
  terms.avg <- mclapply(common, tmp_f, all=all, mc.cores=1)
  all_avg <- rbindlist(terms.avg)
  
  # order by p
  all_avg <- all_avg[order(all_avg$p_adj),]
  dim(all_avg[which(all_avg$q_value <= 0.01),])
  dim(all_avg[which(all_avg$p_adj <= 0.05),])
  
  # write output for further analyses
  write.table(all_avg[, c("ID", "p_adj")], paste("/Users/nicco/RESULTS_", random_num, "/revigo_inp_clusterProfiler.txt", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
  
  # also write the full list of enrichment results
  write.table(all_avg, paste("/Users/nicco/RESULTS_", random_num, "/enrichent_results_sampling.txt", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
  
  return(all_avg)
}

## function to merge sampling procedure for the functional annotation
mergeSampling <- function(enrich.res){
    n.sampl = length(enrich.res)
    all.enrich <- do.call("rbind", enrich.res)      #put all results together
    sign.only <- all.enrich[which(all.enrich$p_value <= 0.05),]         #isolate only significant findings
    if (nrow(sign.only) >0){
      sign.table <- as.data.frame(table(sign.only$term_name))
      info <- sign.only[, c("term_name", "term_id")]
      sign.table <- sign.table[order(-sign.table$Freq),]
      sign.table <- merge(sign.table, info, by.x="Var1", by.y="term_name")
    } else { sign.table <- data.frame(chr=NULL, position=NULL)}
    
    # also try with averaging pvalues
    if ("intersection" %in% colnames(all.enrich)){
      terms <- all.enrich[, c("term_name", "term_id", "intersection", "p_value")]
      slow <- TRUE
    } else {
      terms <- all.enrich[, c("term_name", "term_id", "p_value")]
      slow = FALSE
    }
    # see how many times each term was tested and resctrict to common (tested all times)
    ttt = as.data.frame(table(terms$term_id))
    ttt <- ttt[which(ttt$Freq == n.sampl),]
    terms <- terms[which(terms$term_id %in% ttt$Var1),]
    #terms <- terms[!duplicated(terms$term_name),]
    avgP <- function(i, terms, term_list, slow){
      if (isTRUE(slow)){
        all_intersections <- paste0(unique(unlist(strsplit(paste0(terms$intersection[which(terms$term_id == term_list[i])], collapse = ","), ","))), collapse = ",")
        df <- data.frame(term_name = unique(terms$term_name[which(terms$term_id == term_list[i])]), term_id = term_list[i], avgP = mean(terms$p_value[which(terms$term_id == term_list[i])]), intersection = all_intersections)
      } else {
        df <- data.frame(term_name = unique(terms$term_name[which(terms$term_id == term_list[i])]), term_id = term_list[i], avgP = mean(terms$p_value[which(terms$term_id == term_list[i])]))
      }
      return(df) 
    }
    term_list <- unique(terms$term_id)
    all.terms.avg <- rbindlist(mclapply(1:length(term_list), avgP, terms=terms, term_list = term_list, slow = slow, mc.cores=2))
    all.terms.avg <- all.terms.avg[order(all.terms.avg$avgP),]
    all.terms.avg$log10P <- -log10(all.terms.avg$avgP)
    print(head(all.terms.avg))
    #write.table(all.terms.avg, paste("RESULTS_", random_num, "/enrichent_results_sampling.txt", sep=""), quote=F, row.names=F, sep="\t")
    
    # if (nrow(all.terms.avg[which(all.terms.avg$avgP <= 0.05),]) >0){
    #   pdf(paste("RESULTS_", random_num, "/treemap_sampling_PRS_genes.pdf", sep=""), height=10, width=10)
    #   treemap(all.terms.avg[which(all.terms.avg$avgP <= 0.05),], index="term_name", vSize="log10P", type="index")
    #   invisible(dev.off())
    # }

    l <- list(sign.table, all.terms.avg)
    return(l)
}

## function to do clustering of GO based on semantic similarity
semSim <- function(avg_pvalues){
    hsGO <- godata('org.Hs.eg.db', ont="BP")        #set database
    goSim("GO:0007399", "GO:0007399", semData=hsGO, measure="Wang") # example for two terms

    #extract only go terms significant
    onlyGO <- avg_pvalues[grep("GO:", avg_pvalues$term_id), ]
    onlyGO_sign <- onlyGO[which(onlyGO$avgP <= 0.05),]
    go1 = as.character(onlyGO_sign$term_id)
    go2 = as.character(onlyGO_sign$term_id)
    mt <- mgoSim(go1, go2, semData=hsGO, measure="Wang", combine=NULL)        # example for two vectors of terms
    pdf(paste("RESULTS_", random_num, "/PRS_ANNOTATION/corrplot_significant_terms_GO.pdf", sep=""), height=12, width=12)
    par(mfrow=c(1,1), mar=c(5, 5, 5, 5))
    pheatmap(mt)
    invisible(dev.off())
}

## function to read input set of snp given the name and the input type
readSNPs <- function(fname, ftype, MAIN){
    ## read input file
    d <- fread(fname, h=F)

    ## check type: type 1 --> chr:pos; type 2 --> chr   pos; type 3 --> rsid; type 4 --> rsid_A1
    if (ftype == 1){
        tmp <- str_split_fixed(d$V1, ":", 2)
        d$chr <- as.numeric(as.character(tmp[, 1]))
        d$pos <- as.numeric(as.character(tmp[, 2]))
        colnames(d) <- c("locus", "chr", "pos")
    } else if (ftype == 2){
        d$locus <- paste(d$V1, d$V2, sep=":")
        colnames(d) <- c("chr", "pos", "locus")
    } else if (ftype %in% c(3, 4)){
        if (ftype == 4){ d <- as.data.frame(str_split_fixed(d$V1, "_", 2)); d$V2 <- NULL }
        ref <- fread(paste(MAIN, "INPUTS_OTHER/1000G_frequencies/chrAll.afreq.gz", sep=""), h=T, showProgress=FALSE)
        info <- ref[which(ref$ID %in% d$V1), c("#CHROM", "POS", "ID", "ALT_FREQS")]
        miss <- d$V1[which(!(d$V1 %in% info$ID))]
        colnames(info) <- c("chr", "pos", "ID", "ALT_FREQS")
        info$locus <- paste(info$chr, info$pos, sep=":")
        miss_df <- data.frame(chr = rep(NA, length(miss)), pos = rep(NA, length(miss)), ID = miss, ALT_FREQS = rep(NA, length(miss)), locus = rep(NA, length(miss)))
        d <- rbind(info, miss_df)
    }
    return(d)
}

## function to create n sampling dsets
samplingGset <- function(i, mapping){
    all.genes <- strsplit(mapping$geneList, ",")
    tmp_f <- function(x, all.genes){ tmp <- all.genes[[x]]; g <- sample(x=tmp, size=1); return(g)}
    gset <- lapply(1:length(all.genes), tmp_f, all.genes=all.genes)
    gset.clean <- unlist(gset)
    return(gset.clean)
}

## function to run REVIGO
revigo <- function(avg_pvalues, thr){
    # take significant only
    sig <- avg_pvalues[which(avg_pvalues$avg_p <= thr),]
    sig <- sig[, c("term_id", "avg_p")]

    # do a first check and in case increase p to 0.05    
    if (nrow(sig) <= 30){
      thr = 0.05
      sig = avg_pvalues[which(avg_pvalues$avg_p <= thr),]
      sig <- sig[, c("term_id", "avg_p")]
    }
    # do a second check and in case increase p to 0.10 -- never more than this    
    if (nrow(sig) <= 30){
      thr = 0.10
      sig = avg_pvalues[which(avg_pvalues$avg_p <= thr),]
      sig <- sig[, c("term_id", "avg_p")]
    }
    
    # should I check the number of significant terms and in case restrict it?
    # let's try to go to fdr 1%
    #sig <- avg_pvalues[which(avg_pvalues$avg_p <= 0.01),]
    
    #check if any significant
    if (nrow(sig) >0){
        # output revigo input
        write.table(sig[, c("term_id", "avg_p")], paste("RESULTS_", random_num, "/revigo_inp.txt", sep=""), quote=F, row.names=F, col.names=F, sep="\t")    

        # set parameters for REVIGO
        # distance.meas can be: Lin -- SIMREL
        # clust.simil can be: 0.7 -- 0.5 -- 0.4
        distance.meas <- "Lin"
        clust.simil <- 0.40

        # run python script to get results
        cat("## Running REVIGO\n")
        cmd = paste("python ", MAIN, "BIN/parseREVIGO.py ", clust.simil, " ", distance.meas, " ", random_num, sep="")
        system(cmd)

        # read output
        d <- fread(paste("RESULTS_", random_num, "/revigo_out.csv", sep=""), h=T, sep=",")

        # remove input as it is not too informative
        #cmd <- "rm RESULTS/revigo_inp.txt"
        #system(cmd)

        pdf(paste("RESULTS_", random_num, "/clustering_GO_terms.pdf", sep=""), height=11, width = 11)
        revigoPlot(d)
        invisible(dev.off())
    } else {
        print(thr)
        cat("  !!! No significant GO terms enriched\n")
        d <- "!!! No significant GO terms enriched in the submitted list of SNPs."
        write.table(d, paste("RESULTS_", random_num, "/revigo_out.csv", sep=""), sep="\t")
    }
    return(d)
}

## function to plot REVIGO results -- to finish
revigoPlot <- function(one.data){
    # take only terms that were not merged
    one.data <- one.data[which(one.data$eliminated == 0),]

    # do some adjustements
    one.data$plot_X <- as.numeric(one.data$plot_X)
    one.data$plot_Y <- as.numeric(one.data$plot_Y)
    
    ## base plot
    par(mar=c(6, 5, 5, 10))
    plot(0, 0, pch=16, col="white", xlim=c(min(one.data$plot_X)-1, max(one.data$plot_X)+1), ylim=c(min(one.data$plot_Y)-1, max(one.data$plot_Y)+1), xlab="Semantic space X", ylab="Semantic space Y", cex.lab=1.50)
    
    ## grid
    for (i in seq(min(one.data$plot_X)-1, max(one.data$plot_X)+1, (max(one.data$plot_X)-min(one.data$plot_X)+2)/10)){abline(v=i, lwd=0.4, col="grey80")}
    for (i in seq(min(one.data$plot_Y)-1, max(one.data$plot_Y)+1, (max(one.data$plot_Y)-min(one.data$plot_Y)+2)/10)){abline(h=i, lwd=0.4, col="grey80")}
    
    ## colors
    colz <- colorRampPalette(c("red", "orange", "yellow"))
    colorz <- colz(nrow(one.data))
    one.data <- one.data[order(one.data$'log10 p-value'),]
    one.data$col <- colorz

    ## points
    points(x = one.data$plot_X, y = one.data$plot_Y, pch=16, cex=one.data$plot_size+5, col=alpha(colorz, 0.65), lwd=2)

    # need to add some \n in the name description
    one.data$upd_name <- NA
    for (i in 1:nrow(one.data)){
        # count the characters
        tmp.string <- as.character(one.data$description[i])
        
        # calculate in how many lines to distribute the name -- max is 23 characters (including spaces) per line
        n.lines <- ceiling(nchar(tmp.string)/25)
        
        # divide by space depending on n.lines
        tmp <- unlist(strsplit(as.character(one.data$description[i]), " "))
        if (n.lines >1){ 
            x <- split(tmp, ceiling(seq_along(tmp)/n.lines))
            for (j in 1:length(x)){ x[[j]] <- paste(x[[j]], collapse=" ") }
            # add name to data frame
            one.data$upd_name[i] <- paste(x, collapse="\n")
        } else {
            one.data$upd_name[i] <- paste(tmp, collapse=" ")
        }
    }

    # uppercase first letter of sentence
    one.data$upd_name <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", one.data$upd_name, perl = TRUE)

    ## set how many data to annotate -- top n
    if (nrow(one.data) >=25){ n <- 15 } else { n <- nrow(one.data) }
    
    ## set position for annotation
    sb <- one.data[1:n, ]
    addTextLabels(sb$plot_X, sb$plot_Y, sb$upd_name, col.label="black")

    ## legend
    mx.x <- max(one.data$plot_X)+2 + max(one.data$plot_X)*0.2
    mx.y <- max(one.data$plot_Y)+1
    step <- 0.25
    text(x = mx.x, y = mx.y, labels = "-Log10(P)", cex=1.20, xpd=T, font=2)
    gradient.rect(xleft = mx.x-(step*2), ybottom = mx.y-(step*9), xright = mx.x+(step*2), ytop = mx.y-(step*2), col=rev(colorz), nslices = length(colorz), gradient = "horizontal")
    text(x = mx.x+(step*3), y = mx.y-(step*2), labels = round(max(abs(one.data$'log10 p-value')), 1), adj = 0.5, xpd=T, font=2, cex=0.80)
    text(x = mx.x+(step*3), y = mx.y-(step*5.5), labels = round(min(abs(one.data$'log10 p-value')), 1), adj = 0.5, xpd=T, font=2, cex=0.80)
    text(x = mx.x+(step*3), y = mx.y-(step*9), labels = round(median(abs(one.data$'log10 p-value')), 1), adj = 0.5, xpd=T, font=2, cex=0.80)
    
    ## legend for circle size
    text(x = mx.x, y = mx.y-(step*15), labels = "Term size", cex=1.20, font=2, xpd=T)
    points(x = mx.x, y = mx.y-(step*18), pch=16, col="grey80", cex=min(one.data$plot_size+5), xpd=T)
    points(x = mx.x, y = mx.y-(step*23), pch=16, col="grey80", cex=median(one.data$plot_size+5), xpd=T)
    points(x = mx.x, y = mx.y-(step*29), pch=16, col="grey80", cex=max(one.data$plot_size+5), xpd=T)
}

## find variants LD
findLD <- function(i, data){
  # first, restrict to chromosome of choice
  sb <- data[i, ]

  # parameters for LD
  w = 250
  r = 0.10
  
  # command
  plink_dir = "/Users/home/Desktop/SNPbrowser/gitHub_version/public_version_BIN/plink"
  kg_dir = "/Users/home/Desktop/SNPbrowser/gitHub_version/data/databases/1000G_eur/chr"
  cmd = paste(plink_dir, " --bfile ", kg_dir, sb$chr, "_eur --r2 --ld-snp ", sb$ID, " --ld-window-kb ", w, " --out ld_", sb$ID, sep="")
  system(cmd, ignore.stdout = T)
  
  # read back
  ld <- fread(paste("ld_", sb$ID, ".ld", sep=""), h=T)
  
  # clean
  system(paste("rm ld_", sb$ID, ".*", sep=""))
  
  return(ld)
}

## function to polish variant-gene mapping by removing genes that are not recognized in gene-set overlap analysis
PolishAnnotation <- function(annot, geneList){
  # run gene set overlap analysis
  res <- suppressMessages(gost(geneList, organism = "hsapiens", ordered_query = FALSE, multi_query = FALSE, significant = FALSE, exclude_iea = FALSE,
                               measure_underrepresentation = FALSE, evcodes = FALSE, user_threshold = 1, correction_method = "bonferroni",
                               domain_scope = "annotated", custom_bg = NULL, numeric_ns = "", sources = c("GO:BP")))
  
  # isolate mapped and mismapped genes
  used <- names(as.data.frame(res$meta$genes_metadata$query$query_1$mapping))
  used <- str_replace_all(string = used, pattern = "\\.", replacement = "-")
  not_used <- geneList[which(!(geneList %in% used))]
  
  # exclude not used genes from geneList
  geneList <- geneList[which(!(geneList %in% not_used))]
  
  # need to update the variant-gene annotation removing the genes that were not used
  for (g in not_used){
    # find to which variant it belonged to
    tmp_ind <- grep(g, annot$geneList)
    
    # extract gene list, exclude gene and re-assign
    g_list <- unlist(strsplit(annot$geneList[tmp_ind], ","))
    g_list <- paste0(g_list[which(g_list != g)], collapse=",")
    annot$geneList[tmp_ind] <- g_list
  }
  annot$geneList[which(annot$geneList == "")] <- NA
  
  return(list(annot, geneList))
}

## function to do the gene-set enrichment analysis using the new method
newGeneSetEnrichment <- function(i, gene_list){
  g <- gene_list[[i]]
  
  ego2 <- enrichGO(gene = g, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
  res <- ego2@result
  return(res)
}

## function to do the gene-set enrichment analysis using enrichR function
overlapAnalysis_enrichR <- function(i, gene_list, dbs){
  tmp = gene_list[[i]]
  enriched = enrichr(tmp, dbs)
  return(enriched)
}

## function to merge the sampling based gene-set enrichment
mergeSampling_enrichR <- function(enrich.res, dbs){
  n_sources = length(dbs)
  n.sampl <- length(enrich.res)
  # put all values (of all permutations) into a single df
  for (i in 1:length(enrich.res)){
    for (j in 1:n_sources){
      if (i == 1){
        enrich.res[[i]][[j]]$iteration <- i
      } else {
        enrich.res[[i]][[j]]$iteration <- i
        enrich.res[[1]][[j]] <- rbind(enrich.res[[1]][[j]], enrich.res[[i]][[j]])
      }
    }
  }
  
  # then do the means
  final_list = list()
  for (x in 1:n_sources){
    counts <- as.data.frame(table(enrich.res[[1]][[x]]$Term))
    counts = counts[which(counts$Freq == n.sampl),]
    df = as.data.frame(matrix(data=NA, nrow=1, ncol=4))
    colnames(df) <- c("source", "term", "avg_p", "genes")
    for (i in 1:nrow(counts)){
      tmp = enrich.res[[1]][[x]][which(enrich.res[[1]][[x]]$Term == counts$Var1[i]),]
      df[i, ] = c(dbs[x], tmp$Term[1], sum(tmp$Adjusted.P.value)/n.sampl, paste0(unique(unlist(strsplit(tmp$Genes, ";"))), collapse = ","))
    }
    final_list[[x]] <- df[order(df$avg_p),]
  }
  return(final_list)
}

## function to perform alternative version of revigo -- do semantic similarity myself
alternative_revigo_results <- function(MAIN, random_num, go_data){
  # check if revigo input is there
  filelist = system(paste0("ls RESULTS_", random_num, "/"), intern = T)
  if ("revigo_inp.txt" %in% filelist){
    
    # run the python script
    cmd = paste0("python3 ", MAIN, "BIN/Alternative_REVIGO.py RESULTS_", random_num, "/revigo_inp.txt RESULTS_", random_num, "/alternative_Lin_distance.txt ")
    system(cmd)
    
    # read output back
    data = fread(paste0("RESULTS_", random_num, "/alternative_Lin_distance.txt"), h=T)
    
    # Plot heatmap and save it
    png(paste0("RESULTS_", random_num, "/pheatmap_lin_distance.png"), height=7, width=7, res=300, units="in")
    pheatmap(data, show_rownames = F, show_colnames = F, main="Semantic similarity matrix of GO terms")
    dev.off()
    
    hr1 <- hclust(as.dist(1-data), method="ward.D2", members = NULL)
    
    # loop until at least 3 clusters are found
    n_clust = 0
    minModuleSize = 15
    while(n_clust < 3){
      dynamicMods = cutreeDynamic(dendro = hr1, distM = data,
                                  minClusterSize = minModuleSize, deepSplit = TRUE, 
                                  verbose = 0, respectSmallClusters = T)
      tb <- as.data.frame(table(dynamicMods))
      n_clust <- max(as.numeric(as.character(tb$dynamicMods)))
      df_cluster_terms = data.frame(term = colnames(data), cluster = dynamicMods)
      df_cluster_terms$cluster_in_dendro = NA
      minModuleSize = minModuleSize - 1
    }
    
    # Hierarchical clustering on the distance matrix
    hr <- as.dendrogram(hclust(as.dist(1-data), method="ward.D2", members = NULL))
    hc <- hclust(as.dist(1-data), method="ward.D2", members = NULL)
    order_dendro <- hc$order
    cluster_number_dendro = 1
    for (i in 1:length(order_dendro)){
      # get element
      tmp_go = colnames(data)[order_dendro[i]]
      # check whether it was already annotated
      if (is.na(df_cluster_terms$cluster_in_dendro[which(df_cluster_terms$term == tmp_go)])){
        # get previous cluster number
        prev_cluster = df_cluster_terms$cluster[which(df_cluster_terms$term == tmp_go)]
        df_cluster_terms$cluster_in_dendro[which(df_cluster_terms$cluster == prev_cluster)] <- cluster_number_dendro
        cluster_number_dendro <- cluster_number_dendro + 1
      }
      # 
    }
    df_cluster_terms$term <- as.character(df_cluster_terms$term)
    df_cluster_terms = df_cluster_terms[match(hr1$labels[hr1$order], df_cluster_terms$term),]
    
    if (nrow(data) >= 120){
      plt_w = 28
    } else {
      plt_w = 17
    }
    # Plot dendrogram and rectangles of the clusters
    df_cluster_terms$col <- NA
    library(viridis)
    palett = viridis(n = max(df_cluster_terms$cluster_in_dendro)+1, option = "plasma")[1:max(df_cluster_terms$cluster_in_dendro)]
    for (i in 1:max(df_cluster_terms$cluster_in_dendro)){
      df_cluster_terms$col[which(df_cluster_terms$cluster_in_dendro == i)] <- palett[i]
    }
    png(paste0("RESULTS_", random_num, "/dendrogram_GOterms.png"), height=7, width=plt_w, res=300, units="in")
    d1=color_branches(dend = hr, clusters = df_cluster_terms$cluster_in_dendro, groupLabels = T, col = palett)
    d1 %>% set("labels_col", df_cluster_terms$col) %>% set("branches_lwd", 2) %>% set("labels_cex", 0.80) %>% plot(main = paste0("Dynamic Cut-tree algorithm ~ Ward.D2 ~ min. size=", minModuleSize))
    dev.off()
    
    # Merge with terms description
    go_data_desc = go_data[, c("term_name", "term_id")]
    clusters <- merge(df_cluster_terms, go_data_desc, by.x="term", by.y="term_id")
    clusters$col <- NULL
  } else {
    clusters = NA
    data = NA
    n_clust = NA
  }  
  # return
  return(list(clusters, data, n_clust))
}

# function to count the number of occurrences of each word given a vector of sentences
CountFrequency_words <- function(functional_clusters, n_clust){
  # Read the description of all GO:BP -- this will be the background to remove redundant words
  all_go_bp <- fread("/Users/home/Desktop/SNPbrowser/gitHub_version/AnnotateMe/BIN/go_terms_BP_all.txt", h=F, stringsAsFactors = F, sep="\t")
  
  colnames(functional_clusters) <- c("term", "cluster", "cluster_in_dendro", "term.y")
  # need to remove the (GO id) from the term name
  # for (i in 1:nrow(functional_clusters)){
  #   functional_clusters$term.y[i] <- unlist(strsplit(functional_clusters$term.y[i], "[()]"))[1:(length(unlist(strsplit(functional_clusters$term.y[i], "[()]"))) -1)]
  # }
  
  # Put all different words datasets in a list
  all_dset <- list()
  for (i in 1:n_clust){
    if (i == 1){
      all_dset[[i]] <- all_go_bp$V1
      all_dset[[i+1]] <- functional_clusters$term.y[which(functional_clusters$cluster_in_dendro == i)]
    } else {
      all_dset[[i+1]] <- functional_clusters$term.y[which(functional_clusters$cluster_in_dendro == i)]
    }
  }

  # Define output
  word_frequency_dset <- list()
  
  # Main loop to calculate word frequency
  for (i in 1:length(all_dset)){
    # Convert to dataframe and make sure there are characters only
    text <- data.frame(text=all_dset[[i]])
    text$text <- as.character(text$text)
    
    # Covert dataframe to tibble
    text_df <- tibble(line = 1:nrow(text), text = text)
    text_df <- mutate(text_df, text = text$text)
    
    # Count word occurrence
    x <- text_df %>% unnest_tokens(word, text) %>%    # split words
      anti_join(stop_words) %>%    # take out "a", "an", "the", etc.
      count(word, sort = TRUE)    # count occurrences
    
    # Convert back to dataframe
    x <- as.data.frame(x)
    colnames(x) <- c("word", "freq")
    
    # Assign to output
    word_frequency_dset[[i]] <- x
  }
  
  # Extract background and calculate % and take top 5%
  bkg <- word_frequency_dset[[1]]
  bkg$perc <- bkg$freq/nrow(all_go_bp)
  top_5_pc <- bkg[which(bkg$perc >= 0.025),]
  
  # Now exclude these frequent words from the clusters names
  for (i in 2:length(word_frequency_dset)){
    sb <- word_frequency_dset[[i]]
    sb <- sb[which(!(sb$word %in% top_5_pc$word)),]
    word_frequency_dset[[i]] <- sb
  }
  
  # Try some plot -- seems like it is not possible to plot multiple plots in one --> make 4 plots and save them from "Export"
  for (i in 2:length(word_frequency_dset)){ 
    myplot <- wordcloud2(data = word_frequency_dset[[i]]) 
    # save it in html
    saveWidget(myplot, "tmp.html", selfcontained = F)
    # conert to png
    webshot("tmp.html", paste0("RESULTS_", random_num, "/cluster_", i-1, "_wordcloud.png"), delay = 20, vwidth = 1000, vheight = 1000)
  }
  
  return(word_frequency_dset)
}

## function to plot enrichment results other than go
plotOtherEnrichments <- function(sampling.res, dbs, random_num, type){
  if (type == "enrichR"){
    library(enrichplot)
    library(viridis)
    # exclude GO
    pdf(paste0("RESULTS_", random_num, "/Enrichment_results_noGO.pdf"), height = 12, width = 12, onefile = T)
    par(mar = c(6, 16, 4, 4))
    for (x in 1:length(dbs)){
      if (dbs[x] != "GO_Biological_Process_2018"){
        tmp = sampling.res[[x]]
        tmp_sig <- tmp[which(tmp$avg_p <= 0.10),]
        tmp_sig$logP <- -log10(as.numeric(tmp_sig$avg_p))
        max_p <- max(tmp_sig$logP, 3)
        pos <- barplot(tmp_sig$logP, horiz = T, col=viridis(n = nrow(tmp_sig), option = "plasma"), xlab = "-log10(p)", xlim=c(0, max_p), main = tmp_sig$source[1], cex.lab=1.60, cex.axis = 1.25)
        abline(v = -log10(0.05), lty=2, col="red", lwd=2)
        abline(v = -log10(0.10), lty=2, col="deepskyblue3", lwd=2)
        for (i in 1:length(pos)){
          # look into spacing
          n.lines <- ceiling(nchar(tmp_sig$term[i])/32)
          
          # divide by space depending on n.lines
          tmp <- unlist(strsplit(as.character(tmp_sig$term[i]), " "))
          if (n.lines >1){ 
            kk <- split(tmp, ceiling(seq_along(tmp)/n.lines))
            for (j in 1:length(kk)){ kk[[j]] <- paste(kk[[j]], collapse=" ") }
            # add name to data frame
            tmp_sig$term[i] <- paste(kk, collapse="\n")
          } else {
            tmp_sig$term[i] <- paste(tmp, collapse=" ")
          }
          text(x = 0, y = pos[i, 1], labels = tmp_sig$term[i], xpd=T, pos=2, offset = 1)
          legend("topright", legend = c("FDR<10%", "FDR<5%"), lty = c(2,2), col=c("deepskyblue3", "red"), cex=1.5, ncol=1, bty="n", lwd=3)
        }
      }
    }
    dev.off()
  } else if (type == "gost"){
    library(viridis)
    # determine source of enrichment
    tmp <- sampling.res[[2]]
    tmp$source_gset <- str_split_fixed(tmp$term_id, ":", 2)[, 1]
    # exclude GO
    tmp <- tmp[which(tmp$source_gset != "GO"),]
    pdf(paste0("RESULTS_", random_num, "/Enrichment_results.pdf"), height = 12, width = 12, onefile = T)
    for (x in unique(tmp$source_gset)){
      dt = tmp[which(tmp$source_gset == x),]
      tmp_sig <- dt[which(dt$avgP <= 0.10),]
      if (nrow(tmp_sig) >=20){
        if (nrow(tmp_sig) > 45){
          tmp_sig <- tmp_sig[order(tmp_sig$avgP),]
          tmp_sig <- head(tmp_sig, 45)
        }
        par(mar = c(6, 17, 4, 4))
        cex_text = 0.50
        splt_words = 100
      } else {
        par(mar = c(6, 17, 4, 4))
        cex_text = 1
        splt_words = 40
      }
      max_p <- max(ceiling(tmp_sig$log10P), 3)
      if (nrow(tmp_sig) >0){
        pos <- barplot(tmp_sig$log10P, horiz = T, col=viridis(n = nrow(tmp_sig), option = "plasma"), xlab = "-log10(p)", xlim=c(0, max_p), main = dt$source_gset[1], cex.lab=1.60, cex.axis = 1.25)
        abline(v = -log10(0.05), lty=2, col="red", lwd=2)
        abline(v = -log10(0.10), lty=2, col="deepskyblue3", lwd=2)
        tmp_sig$term_name <- as.character(tmp_sig$term_name)
        for (i in 1:length(pos)){
          # look into spacing
          n.lines <- ceiling(nchar(tmp_sig$term_name[i])/splt_words)
          
          # divide by space depending on n.lines
          tmp_w <- unlist(strsplit(as.character(tmp_sig$term_name[i]), " "))
          if (n.lines >1){ 
            pieces = floor(length(tmp_w)/n.lines)
            if (n.lines == 2){
              p1 = paste(tmp_w[1:pieces], collapse = " ")
              p2 = paste(tmp_w[(pieces+1):length(tmp_w)], collapse = " ")
              tmp_sig$term_name[i] <- paste(p1, p2, sep="\n")
            } else if (n.lines == 3){
              p1 = paste(tmp_w[1:pieces], collapse = " ")
              p2 = paste(tmp_w[(pieces+1):(pieces*2)], collapse = " ")
              p3 = paste(tmp_w[(pieces*2+1):length(tmp_w)], collapse = " ")
              tmp_sig$term_name[i] <- paste(p1, p2, p3, sep="\n")
            }
            # kk <- split(tmp_w, ceiling(seq_along(tmp_w)/n.lines))
            # for (j in 1:length(kk)){ kk[[j]] <- paste(kk[[j]], collapse=" ") }
            # # add name to data frame
            # tmp_sig$term_name[i] <- paste(kk, collapse="\n")
          } else {
            tmp_sig$term_name[i] <- paste(tmp_w, collapse=" ")
          }
          text(x = 0, y = pos[i, 1], labels = tmp_sig$term_name[i], xpd=T, pos=2, offset = 1, cex = cex_text)
          legend("topright", legend = c("p<0.05", "p<0.10"), lty = c(2, 2), col=c("red", "deepskyblue3"), cex=1.5, ncol=2, bty="n", lwd=3)
        }
      }else {
        plot(0,0, pch=16, col="white", xlim=c(0, 1), ylim=c(0,1), main=dt$source_gset[1], xlab="", xaxt="none", yaxt="none", ylab="", bty="n")
        text(x = 0.5, y = 0.5, labels = "No significant hits")
      }
    }
    dev.off()
  }
}

## function to find sv close to input snps and genes
function_findSVs <- function(i, annot, all_str_hg38){
  # get locus
  tmp_chr = as.character(annot$chr_hg38[i])
  tmp_pos = as.numeric(as.character(annot$pos_hg38[i]))
  w = 10000
  # restrict sv of the same chromosome and pm 1MB
  all_str_hg38$chr = as.character(all_str_hg38$chr)
  all_str_hg38$end_pos <- as.numeric(all_str_hg38$end_pos)
  tmp_sv = all_str_hg38[which(all_str_hg38$chr == tmp_chr),]
  # find svs
  # snp inside sv
  sv_tmp_inside = tmp_sv[which(tmp_sv$start_pos <= tmp_pos & as.numeric(tmp_sv$end_pos) >= tmp_pos),]
  sv_tmp_inside$SV_SNP_relation = "SNP_inside_SV"
  sv_tmp_inside$dst = NA
  tmp_sv <- tmp_sv[which(!(tmp_sv$start_pos %in% sv_tmp_inside$start_pos)),]
  # look at upstream (negative distance from start)
  tmp_sv$dst <- tmp_sv$start_pos - tmp_pos
  upstream = tmp_sv[which(tmp_sv$dst >0 & tmp_sv$dst <w),]
  upstream$SV_SNP_relation = "SNP_upstream_SV"
  # then downstream  
  tmp_sv$dst <- as.numeric(tmp_sv$end_pos) - tmp_pos
  downstream = tmp_sv[which(tmp_sv$dst <0 & tmp_sv$dst >(-w)),]
  downstream$SV_SNP_relation = "SNP_downstream_SV"
  # merge all
  all_sv = rbind(sv_tmp_inside, upstream, downstream)
  # clean up, add variant information and close
  if (nrow(all_sv) >0){
    all_sv$col <- NULL
    colnames(all_sv) <- c("chr_hg38", "start_pos_SV_hg38", "end_pos_SV_hg38", "diff_allele_size_SV", "SV_type", "SV_source", "SV_SNP_relation", "Distance_SNP_SV")
    all_sv$SV_source[which(all_sv$SV_source == "jasper")] <- "Linthrost_et_al_2020"
    all_sv$SV_source[which(all_sv$SV_source == "audano")] <- "Audano_et_al_2019"
    all_sv$SV_source[which(all_sv$SV_source == "chaisson")] <- "Chaisson_et_al_2019"
    all_sv$SNP_locus_hg38 = paste(annot$chr_hg38[i], annot$pos_hg38[i], sep=":")
    all_sv$geneList <- annot$geneList[i]
  } else {
    all_sv = data.frame("chr_hg38" = NA, "start_pos_SV_hg38" = NA, "end_pos_SV_hg38" = NA, "diff_allele_size_SV" = NA, "SV_type" = NA, "SV_source" = NA, "SV_SNP_relation" = NA, "Distance_SNP_SV" = NA, "SNP_locus_hg38" = paste(annot$chr_hg38[i], annot$pos_hg38[i], sep=":"), "geneList" = annot$geneList[i])
  }
  return(all_sv)
}

# MAIN
## read main snp file and do the necessary adjustments
MAIN = "/Users/home/Desktop/SNPbrowser/gitHub_version/AnnotateMe/"
fname <- args[1]
#fname <- "/Users/home/Desktop/SNPbrowser/gitHub_version/snpXplorer_v2/RESULTS_2150/annotateMe_input_13031.txt"
ftype <- args[2]
#ftype <- 1
username <- args[3]
analysis_mode = unlist(strsplit(as.character(args[4]), ","))

## send email myself to notify that someone made a request
cmd_mail <- paste("sendEmail -f snpXplorer@gmail.com -t n.tesi@amsterdamumc.nl -u 'AnnotateMe request sent' -m 'Hello, \n a request to AnnotateMe was just sent from ", username, ". \n \n AnnotateMe' -s smtp.gmail.com:25 -xu snpXplorer@gmail.com -xp snpXplorer22101991!", sep="")
system(cmd_mail)

## create folder for results -- add a random number
random_num <- sample(x = seq(1, 100000), size = 1, replace = F)
cmd = paste("mkdir RESULTS_", random_num, sep="")
system(cmd)
## also put input file with the list of snps in the folder
cmd = paste("cp ", fname, " RESULTS_", random_num, "/", sep="")
system(cmd)

## read input data
data <- readSNPs(fname, ftype, MAIN)
# manual settings for missing data
# data <- as.data.frame(data)
# data[which(data$ID == "rs143080277"), c("chr", "pos", "ALT_FREQS", "locus")] <- c(2, 106366056, 0.005, "2:106366056")
# data[which(data$ID == "rs3822030"), c("chr", "pos", "ALT_FREQS", "locus")] <- c(4, 987343, 0.47, "4:987343")
# data[which(data$ID == "rs75932628"), c("chr", "pos", "ALT_FREQS", "locus")] <- c(6, 41129252, 0.001, "6:41129252")
# data[which(data$ID == "rs60755019"), c("chr", "pos", "ALT_FREQS", "locus")] <- c(6, 41149008, 0.01, "6:41149008")
# data[which(data$ID == "rs3919533"), c("chr", "pos", "ALT_FREQS", "locus")] <- c(12, 7162801, 0.154, "12:7162801")
# data[which(data$ID == "rs616338"), c("chr", "pos",  "ALT_FREQS", "locus")] <- c(17, 47297297, 0.006, "17:47297297")
# ftype <- 1

## separate snps to annotate from the NAs
miss <- data[which(is.na(data$chr)),]
data <- data[!is.na(data$chr),]

## find all variants in LD with the input variants: this will help with the annotation
if (ftype == 3){
  try(ld_info <- lapply(1:nrow(data), findLD, data=data), silent = T)
  try(ld_info <- rbindlist(ld_info), silent = T)
} else {
  ld_info <- data.frame(chr=NULL, pos=NULL)
}
## for AD variants, one is missing, add it manually
# tmp_df <- data.frame(chr = 17, pos = 44353222, ID = "rs2732703", ALT_FREQS = 0.1365, locus = "17:44353222")
# data <- rbind(data, tmp_df)

## read additional needed files: gtex, ensemble-gene_name mapping and gene positions
cat("## Loading all genes and gtex\n")
load("/Users/home/Desktop/SNPbrowser/gitHub_version/data/databases/annotationFiles.RData")
genes <- fread(paste(MAIN, "INPUTS_OTHER/NCBI37.3.gene.loc", sep=""), h=F, stringsAsFactors=F)
colnames(genes) <- c("gene_id", "chr", "start_tss", "stop_tss", "strand", "gene_name")
gtex <- fread(paste(MAIN, "INPUTS_OTHER/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz", sep=""))
mapping.ens <- fread(paste(MAIN, "INPUTS_OTHER/Ensemble_to_GeneName.txt", sep=""))
colnames(mapping.ens) <- c("ensemble", "gene")

## manage cadd files
cat("## Loading CADD scores\n")
cadd <- mclapply(unique(data$chr), readcadd, data=data, MAIN=MAIN, ld_info=ld_info, mc.cores=1)
cadd <- rbindlist(cadd)

## in case there were missings, add them back
cadd$locus <- paste(cadd$"#Chrom", cadd$Pos, sep=":")
miss_cadd <- data[which(!(data$locus %in% cadd$locus)),]

## main function for functional annotation
cat("## Start annotation\n")
res <- AnnotateMe(data, genes, gtex, mapping.ens, ftype, MAIN, ld_info, cadd)
annot <- res[[1]]
geneList <- res[[2]]
# also nice to output the SVs close to the snps
all_sv = rbindlist(mclapply(1:nrow(annot), function_findSVs, annot = annot, all_str_hg38 = all_str_hg38, mc.cores=2))
write.table(all_sv, paste0("RESULTS_", random_num, "/SNP_and_SV_overlap.txt"), quote=F, row.names=F, sep="\t")
cat("## Annotation is done. Now making plots and functional enrichment analysis.\n")

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
  enrich.res <- mclapply(1:n.sampl, overlapAnalysis, gene_list=gene_sampling_dsets, source_gset = dbs, mc.cores=2)
  #enrich.res <- mclapply(1:n.sampl, newGeneSetEnrichment, gene_list=gene_sampling_dsets, mc.cores=2)
  #library(enrichR)
  #enrich.res <- mclapply(1:n.sampl, overlapAnalysis_enrichR, gene_list=gene_sampling_dsets, dbs = dbs, mc.cores=2)
  cat("  Merging and cleaning results\n")
  sampling.res <- mergeSampling(enrich.res)
  #sampling.res <- mergeSampling_enrichR(enrich.res, dbs)
  # in case i use gost (check for genes-pathway relationships) from here i need to change
  # first save all enrichment results
  #for (i in 1:length(sampling.res)){
  write.table(sampling.res[[2]], file = paste0("RESULTS_", random_num, "/geneSet_enrichment_results.txt"), quote=F, row.names=F, sep="\t")
  #}
  # for other enrichment sources like kegg, reactome etc, need to make a plot
  if (length(dbs) >1){
    plotOtherEnrichments(sampling.res, dbs, random_num, type="gost")
  } else if (dbs[1] != "GO_Biological_Process_2018"){
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
    # add source information to select only go
    go_data$source_gset = str_split_fixed(go_data$term_id, ":", 2)[, 1]
    go_data <- go_data[which(go_data$source_gset == "GO"),]
    # last thing is to run REVIGO and make the plot
    try(revigo_res <- revigo(avg_pvalues = go_data, thr = 0.01), silent = T)
    # maybe nice to also do the alternative to revigo -- manual
    semsim_results = alternative_revigo_results(MAIN, random_num, go_data)
    functional_clusters = semsim_results[[1]]
    lin_matrix <- semsim_results[[2]]
    n_clust = semsim_results[[3]]
    
    # let's try to make some wordcloud images
    if (!is.na(functional_clusters)){
      mostFreq_words <- CountFrequency_words(functional_clusters, n_clust)
    }
  }
  cat("\n## Analysis done. Hope results make sense :)\n")
  
  # finally compress result folder and send it
  cmd_compress <- paste("tar -czf AnnotateMe_results_", random_num, ".tar.gz RESULTS_", random_num, "/", sep="")
  system(cmd_compress)
  cmd_mail <- paste("sendEmail -f snpXplorer@gmail.com -t ", username, " -u 'AnnotateMe results' -m 'Dear user, \n thanks so much for using snpXplorer and AnnotateMe. \n We hope you find the tool useful. \n AnnotateMe team.' -a 'AnnotateMe_results_", random_num, ".tar.gz' -cc n.tesi@amsterdamumc.nl -s smtp.gmail.com:25 -xu snpXplorer@gmail.com -xp snpXplorer22101991!", sep="")
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
  cmd_mail <- paste("sendEmail -f snpXplorer@gmail.com -t ", username, " -u 'AnnotateMe results' -m 'Dear user, \n thanks so much for using snpXplorer and AnnotateMe. \n Unfortunately, due to a low number of genes found in your SNPs, we could not perform gene-set enrichment analysis. \n We hope you find the tool useful. \n AnnotateMe team.' -a 'AnnotateMe_results_", random_num, ".tar.gz' -cc n.tesi@amsterdamumc.nl -s smtp.gmail.com:25 -xu snpXplorer@gmail.com -xp snpXplorer22101991!", sep="")
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
# res_genes <- mclapply(1:nrow(gwas.sb), tmp_f, gwas.sb=gwas.sb, intervals=intervals, mc.cores=3)     # annotate every row (~200k rows)
# all.genes <- as.data.frame(rbindlist(res_genes))       # merge results
# all.genes <- all.genes[!is.na(all.genes), ]       # exclude NAs
# write.table(all.genes, "INPUTS_OTHER/Gwas_catalog_Gene_Traits.txt", quote=F, row.names=F, sep="\t")

###################################################


