# basic paths
MAIN = "/root/snpXplorer/AnnotateMe/"
MAIN_SNP = "/root/snpXplorer/snpXplorer_v3/"
args = commandArgs(trailingOnly=TRUE)

## function to annotate using GWAS catalog -- adjusted for faster computations (library-wise)
GWAScat <- function(annot, geneList, MAIN, random_num){
  # read gwas catalog
  gwas <- data.table::fread(paste(MAIN, "INPUTS_OTHER/GWAS_catalog_20211210.txt.gz", sep=""), h=T, showProgress=FALSE, quote="", stringsAsFactors=F)

  # do 2 ways of merging to maximize results: rsid and chr:pos
  gwas$LOCUS <- paste(gwas$CHR_ID, gwas$CHR_POS, sep=":")
  m1 <- gwas[which(gwas$SNPS %in% annot$ID), ]
  m2 <- gwas[which(gwas$LOCUS %in% annot$locus), ]
  all.m <- rbind(m1, m2)
  tmp = all.m
  all.m <- all.m[, c("SNPS", "MAPPED_TRAIT", "MAPPED_TRAIT_URI", "FIRST AUTHOR", "PUBMEDID", "JOURNAL")]
  all.m <- all.m[!duplicated(all.m),]

  # check number of matches
  if (nrow(all.m) >0){
    # separate data to be plotted
    traits <- as.data.frame(table(all.m$"MAPPED_TRAIT"))
    traits <- traits[order(-traits$Freq),]
    traits$id <- seq(1, nrow(traits))

    # run function for lolliplot
    pdf(paste("RESULTS_", random_num, "/gwas_cat_snps_overlap.pdf", sep=""), height=7, width=7)
    lolliPlot_snp(traits, all.m)
    invisible(dev.off())

    # finally write the table
    write.table(x = tmp, file = paste0("RESULTS_", random_num, "/gwas_cat_snps_overlap.txt"), quote=F, row.names=F, sep = "\t")
  } else {
    cat("  !!! No matches for SNP in GWAS catalog\n")
    traits <- NA

    tmp = "!!! No matches for SNP in GWAS catalog\n"
    write.table(x = tmp, file = paste0("RESULTS_", random_num, "/gwas_cat_snps_overlap.txt"), quote=F, row.names=F, sep = "\t")
  }

  # also check genes directly -- first read genes~trait dataframe
  all.genes <- data.table::fread(paste(MAIN, "INPUTS_OTHER/20211210_Gwas_catalog_Gene_Traits.txt", sep=""), h=T, sep="\t", quote="", stringsAsFactors=F)

  # check overlap as a background
  overlapping.genes <- all.genes[which(all.genes$gene %in% geneList),]        # find the overlapping genes with my gene list
  overlapping.genes <- overlapping.genes[!duplicated(overlapping.genes), ]    # exclude duplicated rows
  #print(overlapping.genes)
  # also here check number of matching genes
  if (nrow(overlapping.genes) >0){
    # sampling: each iteration sample 1 gene from the pool of genes associated with each variant and do the overlap
    # at the end, do the mean of the overlaps
    n.samp <- 500
    samp.res <- lapply(1:n.samp, sampleAnnot, overlapping.genes=overlapping.genes, annot=annot)       #sample 1000 times 1 gene per snp and do the overlap 1000 times
    options(warn=-1)
    tb <- BiocGenerics::Reduce(function(x, y) merge(x, y, by="Var1", all.x=T, all.y=T), samp.res)
    tb$SUM <- rowSums(tb[, 2:ncol(tb)], na.rm=T)
    tb <- tb[, c("Var1", "SUM")]
    tb$MEAN <- tb$SUM/n.samp
    tb$MEAN = ceiling(tb$MEAN)
    tb <- tb[order(-tb$MEAN),]
    tb <- tb[!is.na(tb$Var1),]

    # finally lolliplot for genes
    pdf(paste("RESULTS_", random_num, "/gwas_cat_genes_overlap.pdf", sep=""), height=7, width=7)
    lolliPlot_gene(tb, overlapping.genes, geneList)
    invisible(dev.off())
    options(warn=0)

    # also sve the table -- but we need to add gwas catalog information
    tmp_info = gwas[, c("MAPPED_TRAIT", "MAPPED_TRAIT_URI", "GENOTYPING TECHNOLOGY", "STUDY ACCESSION")]
    tmp_info = tmp_info[!duplicated(tmp_info$MAPPED_TRAIT),]
    tb = merge(tb, tmp_info, by.x="Var1", by.y="MAPPED_TRAIT")

    # write table as output
    write.table(x = tb, file = paste0("RESULTS_", random_num, "/gwas_cat_genes_overlap.txt"), quote=F, row.names=F, sep = "\t")
  } else {
    cat("  !!! No matches for Genes in GWAS catalog\n")
    tb <- NA
    tmp = "!!! No matches for SNP in GWAS catalog\n"
    write.table(x = tmp, file = paste0("RESULTS_", random_num, "/gwas_cat_genes_overlap.txt"), quote=F, row.names=F, sep = "\t")
  }

  ls <- list(traits, tb)
  return(ls)
}

## function to sample for multiple genes associated with each variant -- adjusted for faster computations (library-wise)
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

## function to make lolli-plot for genes -- adjusted for faster computations (library-wise)
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
    tmp.study <- as.data.frame(stringr::str_split_fixed(tmp.study$study, " ", 2))
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
  par(mar=c(5, 15, 4, 3))

  # basic empty plot
  plot(0, 0, pch=16, col="white", xlim=c(0, ceiling(max(sb.tr$MEAN))), ylim=c(1, 10), xlab="Fraction of Genes associated with traits", ylab="", cex.axis=1.25, cex.lab=1.50, yaxt="none", bty='n', xaxt='none')

  # axis
  axis(side=1, at=seq(0, ceiling(max(sb.tr$MEAN)), ceiling(max(sb.tr$MEAN)/5)), label=seq(0, ceiling(max(sb.tr$MEAN)), ceiling(max(sb.tr$MEAN)/5)), cex.lab=1.25, cex.axis=1.25)

  # y-axis label
  text(x = -ceiling(max(sb.tr$MEAN))*0.75, y = 5, labels = "GWAS-catalog traits", srt = 90, cex = 2, xpd=T)

  # add grid
  for (i in seq(0, ceiling(max(sb.tr$MEAN)), ceiling(max(sb.tr$MEAN))/10)){ abline(v=i, lwd=0.4, col="grey80") }
  for (i in seq(1, 10)){ segments(x0=0, y0=i, x1=ceiling(max(sb.tr$MEAN)), y1=i, lwd=0.4, col="grey80") }

  # set up colors
  colz <- grDevices::colorRampPalette(c("grey80", "orange", "deepskyblue3", "darkolivegreen3"))
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

## function to make lolli-plot for snps -- adjusted for faster computations (library-wise)
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
    tmp.study <- as.data.frame(stringr::str_split_fixed(tmp.study$'FIRST AUTHOR', " ", 2))
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
  par(mar=c(5, 15, 4, 3))
  # basic empty plot
  plot(0, 0, pch=16, col="white", xlim=c(0, ceiling(max(sb.tr$Freq))), ylim=c(1, 10), xlab="Fraction of SNPs associated with traits", ylab="", cex.axis=1.25, cex.lab=1.50, yaxt="none", bty='n')
  # add grid
  for (i in seq(0, ceiling(max(sb.tr$Freq)), ceiling(max(sb.tr$Freq))/10)){ abline(v=i, lwd=0.4, col="grey80") }
  for (i in seq(1, 10)){ abline(h=i, lwd=0.4, col="grey80") }
  # y-axis label
  text(x = -ceiling(max(sb.tr$Freq))*0.75, y = 5, labels = "GWAS-catalog traits", srt = 90, cex = 2, xpd=T)
  # set up colors
  colz <- grDevices::colorRampPalette(c("red", "orange", "dark green", "blue"))
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

# read arguments
snps_info_path = args[1]
random_num = args[2]
load(snps_info_path)
annot <- final_res[[1]]
geneList <- final_res[[2]]
gwascat_res = GWAScat(annot, geneList, MAIN, random_num)
