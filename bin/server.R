# SERVER SNP BROWSER

# LIBRARIES
#####
suppressPackageStartupMessages({
library(shiny)
library(data.table)
library(shinyBS)
library(stringr)
library(shinyWidgets)
library(viridis)
library(ggplot2)
library(liftOver)
library(colourpicker)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(rvest)
library(shinyjs)
library(stringr)
library(plotrix)
})
par(family  = "Arial")

# FUNCTIONS
#####
## Function to plot
function.plot <- function(snp.info, y.lim, type, plt.type, windows.number, smooth.par,
                          int.locus, gwas, colorPoint, ld, input, genV, inpStrVar, pop_interest_ld, recomb_yn, dotSize_yn){
  # First, let's prepare everything that is needed
  ## Prepare title
  title <- function.title(gwas, int.locus, type, snp.info)
  ## Get recombination tracks if these are requested
  recomb <- findRecomb(snp.info, genV)
  ## Assign dot size
  snp.info <- function.pointSize(dat = snp.info, range = seq(2, 7, 0.5))

  ## Check whether LD needs to be computed and in case calculate it
  ld_yesNo <- "no"
  if (ld == "Most significant in window"){
    ld.info <- functionLD(snp.info, snp=NULL, pop_interest = pop_interest_ld)
    #merge with association data
    ld.snp <- merge(ld.info, snp.info, by.x="BP_B", by.y="pos")
    ld_yesNo <- "yes"
  } else if (ld == "Input variant"){
    ld.info <- functionLD(snp.info, snp=input$pos, pop_interest = pop_interest_ld)
    #merge with association data
    ld.snp <- merge(ld.info, snp.info, by.x="BP_B", by.y="pos")
    ld_yesNo <- "yes"
  }
  ## Also, in case the genome is hg38, now need to do the liftover on the snp.info
  # if genome version is hg38, need to liftover
  if (genV == "GRCh38 (hg38)"){
    df <- data.frame(chr=paste("chr", snp.info$chr, sep=""), start=snp.info$pos, end=snp.info$pos+1)
    # change to GR class object
    gr <- makeGRangesFromDataFrame(df)
    # set chain file
    chain <- import.chain( "../data/databases/hg19ToHg38.over.chain")
    # change coordinates
    gr_hg38 <- liftOver(gr, chain)
    # back to dataframe and clean it
    df_hg38 <- as.data.frame(gr_hg38)
    sb <- data.frame(chr=rep(snp.info$chr[1], nrow(df_hg38)), pos=df_hg38$start, p=snp.info[df_hg38$group, "p"], log10p=snp.info[df_hg38$group, "-log10(P-value)"], size=snp.info[df_hg38$group, "size"])
    colnames(sb) <- c("chr", "pos", "p", "-log10(P-value)", "size")
    snp.info <- sb
    if (ld_yesNo == "yes"){
      df <- data.frame(chr=paste("chr", ld.snp$CHR_A[1], sep=""), start=ld.snp$BP_B, end=ld.snp$BP_B+1)
      gr <- makeGRangesFromDataFrame(df)
      chain <- import.chain( "../data/databases/hg19ToHg38.over.chain")
      gr_hg38 <- liftOver(gr, chain)
      df_hg38 <- as.data.frame(gr_hg38)
      sb <- data.frame(BP_B=df_hg38$start, col=ld.snp[df_hg38$group, "col"], size=ld.snp[df_hg38$group, "size"], p=ld.snp[df_hg38$group, "p"])
      sb$"-log10(P-value)" <- -log10(sb$p)
      ld.snp <- sb
    }
  }
  ## Select genes that are in the window to adjust y-axis
  genes <- function.dynamicGene(snp.info = snp.info, genV)
  if (nrow(genes) > 0){ min.y <- min(genes$y) } else { min.y <- 0 }

  # Then start with the plot
  ## layout is divided in 4 plots: sizes = plt1: 2, plt2: 1, plt3: 2, plt4: 2
  #layout(matrix(c(1,1, 2, 3,3, 4,4,4), nrow = 8, ncol = 1, byrow = TRUE))
  layout(matrix(c(1,1,1,1,1,1,
                  1,1,1,1,1,1,
                  2,2,2,2,2,2,
                  3,3,3,3,3,3,
                  3,3,3,3,3,3,
                  4,6,6,6,6,6,
                  5,7,7,7,7,7,
                  5,7,7,7,7,7,
                  5,7,7,7,7,7), nrow = 9, ncol = 6, byrow = T))
  ## set margins
  par(mar=c(0, 5, 4, 5))

  ## global parameter for the proportions
  pr <- 12
  ## empty plot as background
  plot(x = 0, y = 0, xlab='', cex.lab=2, xaxt='none', ylab="", ylim=c(0, y.lim), cex.axis = 1.5,
       pch=16, col="white", cex=2, type = "p", xaxs="i", yaxt='none', xlim=c(min(snp.info$pos), max(snp.info$pos)), main=title,
       cex.main=2.50, bty='n')

  ## add grid: 10 lines hor. and vert.
  for (x in seq(0, y.lim, (y.lim-min.y*y.lim/10)/10)){abline(h=x, lwd=0.6, col="grey80")}
  for (x in seq(min(snp.info$pos), max(snp.info$pos), (max(snp.info$pos)-min(snp.info$pos))/10)){ segments(x0 = x, y0 = 0, x1 = x, y1 = y.lim, col = "grey80", lwd=0.6) }

  # add recombination rates: for this, need to normalize between 0 and 100 the recombination rate
  recomb$norm.rate <- (y.lim - 1) * (recomb$"Rate(cM/Mb)" / 100)
  y.axis.recomb <- seq(0, 100, 25)
  y.axis.norm <- (y.lim - 1) * ((y.axis.recomb - min(y.axis.recomb))/(max(y.axis.recomb) - min(y.axis.recomb)))
  if (recomb_yn == "Yes"){
    points(recomb$"Position(bp)", recomb$norm.rate, type="l", lwd=1.5, col="darkolivegreen3")
  }

  #add significance lines and corresponding legend -- for now two at 0.05 and genome-wide 5e-8
  abline(h=-log10(0.05), lty=2, col=alpha("darkgreen", 1))
  abline(h=-log10(5e-8), lty=2, col=alpha("purple", 1))
  legend("topleft", bty='n', title = "Annotation lines", legend = c("p=0.05", "p=5e-8"), lty=c(2,2), lwd=c(2,2), col=c("darkgreen", "purple"), cex=1.5, ncol = 2)

  #axis for recombination rates on the right
  if (recomb_yn == "Yes"){
    axis(side = 4, at = y.axis.norm, labels=seq(0, 100, 25), col='darkolivegreen3', col.axis = "darkolivegreen3", cex.axis=1.5, xpd=T)
    text(x = max(snp.info$pos), y = y.lim/5*4, "Recombination rate (cM/Mb)",srt = -90,
       col='darkolivegreen3', xpd=T, pos = 4, offset = 4, cex=2, font=2)
  }
  text(x = min(snp.info$pos), y = y.lim/3*2, "-Log10(P-value)", srt = 90,
       xpd=T, pos = 2, offset = 4, cex=2, font=2)

  # y-axis for pvalue
  axes.x <- ceiling(seq(0, y.lim, y.lim/5))
  axis(side = 2, at = axes.x, labels = axes.x, cex.axis=2)

  ## manage the point-type plot
  if (plt.type == "Points"){
    #then points
    snp.info$"-log10(P-value)" <- -log10(as.numeric(snp.info$p))
    if (dotSize_yn == "No"){ snp.info$size = 2 }
    points(x = snp.info$pos, y = snp.info$"-log10(P-value)", cex.lab=1.5, xaxt='none',
           pch=16, col=alpha(colorPoint, 0.6), cex=snp.info$size, type = "p", xaxs="i", yaxt='none', xlim=c(min(snp.info$pos), max(snp.info$pos)))
    #in case there is ld to be plotted, plot it here
    if (ld != "No LD"){
      if (ld %in% c("Most significant in window", "Input variant")){
        if (dotSize_yn == "No"){ ld.snp$size = 2 }
        points(x = ld.snp$BP_B, y = ld.snp$"-log10(P-value)", cex.lab=1.5, xaxt='none',
               pch=16, col=alpha(ld.snp$col, 1), cex=ld.snp$size, type = "p", xaxs="i", yaxt='none', xlim=c(min(snp.info$pos), max(snp.info$pos)))
      }
      #finally legend
      legend("top", bty="n", legend = c(">0.8", ">0.6", ">0.4", ">0.2"), col=c("red", "orange", "yellow", "deepskyblue3"),
             pt.cex = 2.5, cex=1.50, ncol = 4, pch=16, title = "R2 linkage")
    }

    #if input type was a single snp (either position or rs id), then color the searched variant differently
    if (type %in% c("Position", "Rs ID")){
      #restrict to snp of interest, then plot it
      if (int.locus != "Type position..."){
        snp.interest <- snp.info[which(snp.info$locus == int.locus),]
        points(x = snp.interest$pos, y = snp.interest$"-log10(P-value)", pch=23, lwd=2, col=alpha(colorPoint, 0.8), cex=snp.interest$size)
      }
    }

    ###### HERE STARTS THE DENSITY PLOT
  } else {
    #Sliding window approach -- then loess on sliding window values
    out <- function.DensityLinePvalue(snp.info = snp.info, wind.n = windows.number)
    #print(dim(out))
    lo <- loess(out$pvalue ~ out$window, span=smooth.par)
    xl <- seq(min(out$window), max(out$window), (max(out$window) - min(out$window))/100)
    pred <- predict(lo, xl)
    pred[which(pred < 0)] <- 0

    #add limits -- left and right for polygon function
    xl <- c(xl[1], xl, xl[length(xl)])
    pred <- c(0, pred, 0)
    #lines(x = xl, y = pred, col='navy', lwd=4)

    # add density
    polygon(x = xl, y = pred, col = alpha(colorPoint, 0.6), lwd=3, xaxs="i", border = colorPoint)

    #if input type was a single snp (either position or rs id), then color the searched variant differently -- here is a bar
    if (type %in% c("Position", "Rs ID")){
      #restrict to snp of interest, then plot it
      if (int.locus != "Type position..."){
        snp.interest <- snp.info[which(snp.info$locus == int.locus),]

        #need to grep the height in order to plot the variant -- the height is derived from the prediction
        chr.pos <- str_split_fixed(int.locus, ":", 2)
        pos.only <- as.numeric(chr.pos[, 2])
        pred.pos <- predict(lo, pos.only)

        #add to plot
        points(x = snp.interest$pos, y = pred.pos, type="h", col=alpha("yellow", 0.8), lwd=4, xaxs="i")
      }
    }
  }

  #print("plot1 ok")

  ## global parameter for the proportions
  pr <- 12
  ## empty plot as background
  par(mar=c(5, 5, 0, 5))
  ymax <- 3
  sb <- snp.info
  if (nrow(genes) >ymax){ genes$y <- abs(genes$y); ymax <- max(genes$y)+0.5 }
  plot(x = 0, y = 0, xlab='Chromosomal position (Mb)', cex.lab=2, xaxt='none', ylab="", ylim=c(0, ymax), cex.axis = 1.5,
       pch=16, col="white", cex=2, type = "p", xaxs="i", yaxt='none', xlim=c(min(sb$pos), max(sb$pos)),
       cex.main=2.50, bty='n', yaxt='none')

  text(x = min(sb$pos), y = ymax/5*4, "Gene track", srt = 90,
       xpd=T, pos = 2, offset = 4, cex=2, font=2)

  #manage axes
  axes <- seq(min(sb$pos), max(sb$pos), (max(sb$pos) - min(sb$pos))/7)
  axes.labels <- round(axes/1000000, 3)
  axis(side = 1, at=axes, cex.axis=2, labels=axes.labels)

  # genes
  genes$y <- genes$y-0.5
  if (nrow(genes) >35){
    genes$cex <- 0.40; genes_lwd=0.6; tmp_pr = 0.02
  } else if (nrow(genes) >20){
    genes$cex <- 0.75; genes_lwd=0.6; tmp_pr = 0.03
  } else if (nrow(genes) >15){
    genes$cex <- 0.90; genes_lwd = 1; tmp_pr = 0.04
  } else if (nrow(genes) >10){
    genes$cex <- 1.2; genes_lwd = 1.5; tmp_pr = 0.05
  } else if (nrow(genes) <=5){
    genes$cex <- 2; genes_lwd = 2.5; tmp_pr = 0.085
  } else {
    genes$cex <- 1.5; genes_lwd = 2; tmp_pr = 0.07
  }

  if (min.y != 0){
    genes$y <- abs(genes$y)
    #manage gene names
    for (g in 1:nrow(genes)){
      #main gene line -- full transcription sequence
      segments(x0 = genes$txStart[g], y0 = genes$y[g], x1 = genes$txEnd[g], y1 = genes$y[g], lwd=genes_lwd)
      #need to divide exones from introns
      start <- as.data.frame(t(str_split_fixed(genes$exonStarts[g], ",", genes$exonCount[g]+1)))
      end <- as.data.frame(t(str_split_fixed(genes$exonEnds[g], ",", genes$exonCount[g]+1)))
      exons <- cbind(start, end)
      colnames(exons) <- c("start", "end")
      exons$start <- as.numeric(as.character(exons$start))
      exons$end <- as.numeric(as.character(exons$end))
      #main loop over exons
      for (j in 1:nrow(exons)){
        rect(xleft=exons$start[j], ybottom=genes$y[g]-(ymax/4*(tmp_pr)), xright=exons$end[j],
             ytop = genes$y[g]+(ymax/4*(tmp_pr)), col='grey80', lwd=0.80)
      }
      text(x = genes$txStart[g] + (genes$txEnd[g] - genes$txStart[g])/2, y = genes$y[g]+ymax*tmp_pr,
           labels=genes$"#geneName"[g], font=3, cex=genes$cex[g])
      #        if (genes$strand[g] == "+"){
      #          arrows(x0 = genes$txEnd[g], y0 = genes$y[g], x1 = genes$txEnd[g] + (10*2/100), y1 = genes$y[g],
      #                 length=0.1, lwd=2, col='coral')
      #        } else if (genes$strand[g] == "-"){
      #          arrows(x0 = genes$txStart[g], y0 = genes$y[g], x1 = genes$txStart[g] - (10*2/100), y1 = genes$y[g],
      #                 length=0.1, lwd=2, col='coral')
      #        }
    }
  }

  #print("plot2 ok")
  ################################
  # HERE IS THE THIRD PLOT -- STRUCTURAL VARIATIONS
  par(mar=c(5, 5, 7, 5))
  plot(x = 0, y = 0, xlab='Chromosomal position (Mb)', cex.lab=2, xaxt='none',
         ylab="", ylim=c(0, 10), cex.axis = 2, pch=16, col="white", cex=2, type = "p",
       xaxs="i", yaxt='none', xlim=c(min(snp.info$pos), max(snp.info$pos)),
       cex.main=2.50, bty='n')

  # y-label on the side
  text(x = min(snp.info$pos), y = 8, "Structural variations", srt = 90,
       xpd=T, pos = 2, offset = 4, cex=2, font=2)

    # add grid
  for (x in seq(0, 10, 2)) {abline(h=x, lwd=0.6, col="grey80")}
  for (x in seq(min(snp.info$pos), max(snp.info$pos), (max(snp.info$pos)-min(snp.info$pos))/10)){
    segments(x0 = x, y0 = 0, x1 = x, y1 = 10, col = "grey80", lwd=0.6)
  }

  # manage axes
  axes <- seq(min(snp.info$pos), max(snp.info$pos), (max(snp.info$pos) - min(snp.info$pos))/7)
  axes.labels <- round(axes/1000000, 3)
  axis(side = 1, at=axes, cex.axis=2, labels=axes.labels)
  axes.x <- ceiling(seq(0, y.lim, y.lim/5))

  # take structural variants in the window
  strVar <- findStr_variats(snp.info, genV, inpStrVar)
  # plot as segments (only if there are variants to plot)
  if (nrow(strVar) >15){ strVar$cex <- 1.40 } else { strVar$cex <- 1.80 }
  if (nrow(strVar) > 0){
    for (i in 1:nrow(strVar)){
      segments(x0 = strVar$start_pos[i], y0 = strVar$y_plot[i], x1 = as.numeric(strVar$end_pos[i]),
               y1 = strVar$y_plot[i], lwd = 3, col=strVar$col[i], lty = strVar$lty[i])
      segments(x0 = strVar$start_pos[i], y0 = strVar$y_plot[i]-0.15, x1 = strVar$start_pos[i],
               y1 = strVar$y_plot[i]+0.15, lwd = 3, col=strVar$col[i], lty = strVar$lty[i])
      segments(x0 = as.numeric(strVar$end_pos[i]), y0 = strVar$y_plot[i]-0.15, x1 = as.numeric(strVar$end_pos[i]),
               y1 = strVar$y_plot[i]+0.15, lwd = 3, col=strVar$col[i], lty = strVar$lty[i])
      if (nrow(strVar) <30){
        text(x = strVar$middle[i], y = strVar$y_plot[i]+1, labels = strVar$diff_alleles[i],
           cex = strVar$cex[i], font=3, xpd=T, col=strVar$col[i])
      }
      tmp_lg <- strVar[!duplicated(strVar$type),]
    }
    legend(x = min(snp.info$pos), y = 13, legend = tmp_lg$type, lty = 1, col=tmp_lg$col, xpd=T, cex=1.50, bty='n', lwd=4, ncol=4)
  } else {
     text(x = min(snp.info$pos) + (max(snp.info$pos) - min(snp.info$pos))/2, y = 5,
          labels = "There are no Structural variants to be plotted here", cex=2, adj=0.5)
  }

  #print("plot3 ok")
  # THIRD PLOT -- GENE EXPRESSION FROM GTEX
  # empty plot first
  par(mar=c(4, 2, 10, 0))
  plot(0, 0, xlim=c(0, 55), ylim=c(0, 25), xaxt='none', yaxt='none', xlab="", ylab="",
       main="", cex.main=2.50, col="white", bty="n")
  # run function to take genes (if there are genes)
  if (nrow(genes) > 0){
    gtex.info <- findExpr_gtex(genes, gtex.db)
    heat <- as.data.frame(gtex.info[[1]])
    g <- gtex.info[[2]]
    if (is.na(heat) & is.na(g)){
      par(mar=c(4, 8, 10, 6))
      plot(0, 0, xlim=c(0, 55), ylim=c(0, 25), bty="n", xaxt='none', yaxt='none', xlab="", ylab="",
           main="", cex.main=2.50, col="white")
      plot(0, 0, xlim=c(0, 55), ylim=c(0, 25), bty="n", xaxt='none', yaxt='none', xlab="", ylab="",
           main="RNA expression from GTEx", cex.main=2.50, col="white")
      text(x = 27.5, y = 25, labels = "No GTEx information for the genes in this region.",
           cex = 2, xpd=T)
    } else {

      # then plot as heatmap
      colors <- viridis(n = 101, option = "plasma")

      # scale heat between 0 and 1
      raw <- heat
      heat <- (heat - min(heat))/(max(heat) - min(heat))

      # dendrogram for the genes
      mx <- 22
      if (nrow(heat) >1){
        hr = as.dendrogram(hclust(d = dist(heat, method = "euclidean"), method = "ward.D2"))
        hc = hclust(d = dist(heat, method = "euclidean"), method = "ward.D2")
        ordered_labels = g[hc$order]
        par(mar=c(4, 2, 0, 0))
        nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), cex = 1.25, col = "navy")
        plot(hr, horiz = T, xaxt="none", yaxt ="none", nodePar = nodePar, leaflab = "none", ylim=c(-mx+nrow(heat)+1.5, nrow(heat)+1.5))
      } else {
        plot(0, 0, xaxt="none", yaxt ="none", ylim=c(0, 1), xlim =c(0,1), main="", xlab="", ylab="", bty="n", pch=16, col="white")
        text(x = 0.5, y = 1, labels = "Too few genes\nfor clustering", xpd=T)
        ordered_labels <- g
      }

      # dendrogram for the tissues
      tr = t(heat)
      hr = as.dendrogram(hclust(dist(tr, method = "euclidean"), method = "ward.D2"))
      hc = hclust(dist(tr, method = "euclidean"), method = "ward.D2")
      ordered_labels_tissues = rownames(tr)[hc$order]
      par(mar=c(0, 0, 10, 6))
      nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), cex = 1.25, col = "red")
      plot(hr, horiz = F, xaxt="none", yaxt ="none", nodePar = nodePar, leaflab = "none", xlim=c(0.5, 55.5), main="RNA expression from GTEx", cex.main=2.50)
      # legend of the heatmap goes here
      mx_height = max(hc$height)
      gradient.rect(xleft = 55.5, ybottom = mx_height, xright = 60, ytop = mx_height+mx_height*0.15, col = colors)
      leg.lab <- c(-1, 0, 1)
      pos <- c(55.5, 57.75, 60)
      for (i in 1:length(pos)){
         text(x = pos[i], y = (mx_height+mx_height*0.25), labels = leg.lab[i], cex=1, xpd=T, font=2)
      }

      # background
      par(mar=c(4, 0, 0, 6))
      # reorder the genes to match dendrogram order
      rownames(heat) <- g
      heat = heat[match(ordered_labels, rownames(heat)),]
      tmp = as.data.frame(t(heat))
      tmp = tmp[match(ordered_labels_tissues, rownames(tmp)),]
      heat <- as.data.frame(t(tmp))

      plot(0, 0, xlim=c(0, 55), ylim=c(0, mx), bty="n", xaxt='none', yaxt='none', xlab="", ylab="",
           main="", cex.main=2.50, col="white")
      heat <- as.matrix(heat)
      # loop
      for (i in 1:nrow(heat)){
        for (j in 1:ncol(heat)){
          # value
          v <- round(as.numeric(heat[i, j])*100)
          rect(xleft = j-1, ybottom = mx-i-1, xright = j, ytop = mx-i, col = colors[v+1],
               border = "white", lwd = 1.5)
        }
        text(x = 54, y = mx-i-0.5, labels = g[i], cex=1.50, font=4, xpd=T, adj=0)
      }
      for (j in 1:ncol(heat)){
        text(x = j-0.5, y = mx-i-1.5, labels = ordered_labels_tissues[j], cex=1, font=2, adj=1, xpd=T, srt=90)
      }
    }
  } else {
    plot(0, 0, xlim=c(0, 55), ylim=c(0, 25), bty="n", xaxt='none', yaxt='none', xlab="", ylab="",
         main="", cex.main=2.50, col="white")
    plot(0, 0, xlim=c(0, 55), ylim=c(0, 25), bty="n", xaxt='none', yaxt='none', xlab="", ylab="",
         main="RNA expression from GTEx", cex.main=2.50, col="white")
    text(x = 27.5, y = 25, labels = "No GTEx information for the genes in this region.",
         cex = 2, xpd=T)
  }
  #print("hello")

}

## function to extract expression data from GTEX
findExpr_gtex <- function(genes, gtex.db){
  #subset to take those in the window
  sb <- gtex.db[which(gtex.db$Description %in% genes$`#geneName`), ]

  #polish
  if (nrow(sb) >=1){
    sb <- sb[which(rowSums(sb[, 3:ncol(sb)]) >0),]
    gene.name <- sb$Description
    sb <- sb[, 3:ncol(sb)]

  # try to scale by column
    if (nrow(sb) >1){ sb <- apply(X = sb, MARGIN = 2, FUN = scale) }

    return(list(sb, gene.name))
  } else {
    return(list(NA, NA))
  }
}

## function to calculate LD within the window
functionLD <- function(snp.info, snp, pop_interest){
  chrom <- snp.info$chr[1]
  min.p <- min(snp.info$pos)
  max.p <- max(snp.info$pos)

  #define window length
  dist <- max.p - min.p

  # if LD was requested, we need to check which populations
  pop_file = fread("../data/databases/1000G/people.txt", h=T, stringsAsFactors = F, sep="\t")
  # depending on the selected populations, create a file with suitable individuals
  # first check for "all" individuals -- in case no group is selected, by default uses european
  if (length(pop_interest) == 0){ pop_interest = c("ALL_eur") }
  if ("ALL_eur" %in% pop_interest){  pop_interest = c(pop_interest, "FIN", "GBR", "IBS", "TSI", "FIN", "CEU") }
  if ("ALL_amr" %in% pop_interest){ pop_interest = c(pop_interest, "CLM", "MXL", "PEL", "PUR") }
  if ("ALL_afr" %in% pop_interest){ pop_interest = c(pop_interest, "ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI") }
  if ("ALL_eas" %in% pop_interest){ pop_interest = c(pop_interest, "CDX", "CHB", "CHS", "JPT", "KHV") }
  if ("ALL_sas" %in% pop_interest){ pop_interest = c(pop_interest, "BEB", "GIH", "ITU", "PJL", "STU") }
  # remove duplicates
  pop_interest = pop_interest[!duplicated(pop_interest)]
  tmp_pop = pop_file[which(pop_file$`Population code` %in% pop_interest),]
  tmp_df = data.frame("FID" = tmp_pop$`Sample name`, "IID" = tmp_pop$`Sample name`)
  write.table(tmp_df, "tmp_populations.txt", quote=F, row.names=F, col.names=T, sep="\t")

  # command for plink: calculate LD wrt most sign variant within plot window
  # to account for snps that are not in 1000G, we iterate until at least a match is found
  match = FALSE
  topSNP_index = 1
  while(match == FALSE){
    #identify most significant variant
    if (is.null(snp)){
      snp.info <- snp.info[order(as.numeric(snp.info$p)),]
      mostSign <- snp.info[topSNP_index, ]
      #print(mostSign)
    } else {
      snp.pos <- str_split_fixed(snp, ":", 2)
      position <- as.numeric(snp.pos[, 2])
      mostSign <- snp.info[which(snp.info$pos == position),]
      #print(mostSign)
    }

    cmd_grep <- paste("grep -w ", mostSign$pos, " ../data/databases/1000G/chr", chrom, ".bim", sep="")
    info <- system(cmd_grep, intern = TRUE)
    if (length(info) >0){
      match = TRUE
    } else {
      topSNP_index <- topSNP_index + 1
    }
  }

  # extract variant information and define the command to execute plink
  x <- as.data.frame(str_split_fixed(info, "\t", 6))
  colnames(x) <- c("chr", "rsid", "p", "pos", "a1", "a2")
  cmd <- paste("./plink --bfile ../data/databases/1000G/chr", chrom, " --keep tmp_populations.txt --r --ld-snp ", x$rsid, " --ld-window-kb ", dist, " --out tmp", sep="")
  #print(cmd)
  system(cmd)

  #read file back
  ld <- fread("tmp.ld", h=T)
  ld$col <- NA
  ld$R <- abs(ld$R)
  ld$col[which(abs(ld$R) >= 0.2)] <- "deepskyblue3"
  ld$col[which(abs(ld$R) >= 0.4)] <- "yellow"
  ld$col[which(abs(ld$R) >= 0.6)] <- "orange"
  ld$col[which(abs(ld$R) >= 0.8)] <- "red"

  #eliminate file
  system("rm tmp.*")

  return(ld)
}

## function to extract structural variants from dedicated file
findStr_variats <- function(snp.info, genV, inpStrVar){
  chrom <- snp.info$chr[1]
  min.p <- min(snp.info$pos)
  max.p <- max(snp.info$pos)

  #extract interval of interest
  if (genV != "GRCh38 (hg38)"){
    sb <- all_str[which(all_str$chr == paste("chr", chrom, sep="") & all_str$start_pos >= min.p & all_str$end_pos <= max.p),]
    sb <- sb[order(sb$start_pos),]
  } else {
    all_str_hg38$end_pos <- as.numeric(all_str_hg38$end_pos)
    sb <- all_str_hg38[which(all_str_hg38$chr == paste("chr", chrom, sep="") & all_str_hg38$start_pos >= min.p & all_str_hg38$end_pos <= max.p),]
    sb <- sb[order(sb$start_pos),]
  }
  # then choose which one to plot
  if (inpStrVar == "jasper"){
    sb <- sb[which(sb$source == "jasper"),]
  } else if (inpStrVar == "chaisson"){
    sb <- sb[which(sb$source == "chaisson"),]
  } else if (inpStrVar == "audano"){
    sb <- sb[which(sb$source == "audano")]
  }
  # also if need to plot all three, assign different lty
  if (inpStrVar == "all"){
    sb$lty = 1
    sb$lty[which(sb$source == "chaisson")] <- 2
    sb$lty[which(sb$source == "audano")] <- 3
  }
  #depending on how many there are, assign a y-axis value (0-5)
  if (nrow(sb) > 0){
    sb$y_plot <- NA
    sb$middle <- NA
    c <- 1
    for (i in 1:nrow(sb)){
      sb$y_plot[i] <- c
      #also add middle point
      sb$middle[i] <- as.numeric(as.character(sb$start_pos[i])) + (as.numeric(as.character(sb$end_pos[i])) - as.numeric(as.character(sb$start_pos[i])))/2
      c <- c + 2
      if (c > 9){
        c <- 1
      }
    }
    sb$chr <- as.character(sb$chr)
    sb$end_pos <- as.numeric(as.character(sb$end_pos))
  }
  return(sb)
}

## read recombination map and give coordinates of interests as output
findRecomb <- function(snp.info, genV){
  chrom <- snp.info$chr[1]
  min.p <- min(snp.info$pos)
  max.p <- max(snp.info$pos)

  # take chromosome file
  if (genV == "GRCh38 (hg38)"){
    rf <- genetic.map.hg38[[chrom]]
  } else {
    rf <- genetic.map[[chrom]]
  }

  #extract interval of interest
  recomb <- rf[which(rf$"Position(bp)" >= min.p & rf$"Position(bp)" <= max.p),]

  return(recomb)
}

## function to assign dot size -- define range
function.pointSize <- function(dat, range){
  dat$size <- 2
  for (x in range){
    dat$size[which(-log10(as.numeric(dat$p)) >= x)] <- x
  }
  dat$size[which(dat$size > 6)] <- 6
  return(dat)
}

## function to manage titles
function.title <- function(gwas, int.locus, type, snp.info){
  t = ""
  if (length(gwas) == 1){
    if (gwas == "example"){t = "Example ~ IGAP"} else {t = gwas}

    if (type == "Locus, Gene, RsID"){
      if (int.locus[[1]] %in% c("locus", "rsid")){
        title = paste(t, " ~ chr", snp.info$chr[1], ":", floor(min(snp.info$pos) + (max(snp.info$pos) - min(snp.info$pos))/2), sep = "")
      } else if (int.locus[[1]] == "gene_name"){
        if (int.locus[[1]] == "example"){
          title = paste(t, " ~ chr", snp.info$chr[1], ":", floor(min(snp.info$pos) + (max(snp.info$pos) - min(snp.info$pos))/2), sep = "")
        } else {
          title = paste(t, " ~ ", toupper(as.character(int.locus[[2]])), " ~ chr", snp.info$chr[1], sep="")
        }
      } else {
        title <- paste0(t, " ~ chr", snp.info$chr[1], " ~ Random position")
      }
    } else if (type == "Manual scroll"){
      title <- paste(t, " ~ chr", snp.info$chr[1], ":", floor(min(snp.info$pos) + (max(snp.info$pos) - min(snp.info$pos))/2), sep = "")
    }
  } else {
    # else=multiplot
    if (type == "Locus, Gene, RsID"){
      if (int.locus[[1]] %in% c("locus", "rsid")){
        title = paste("Multiplot ~", snp.info$chr[1], ":", floor(min(snp.info$pos) + (max(snp.info$pos) - min(snp.info$pos))/2), sep = "")
      } else if (int.locus[[1]] == "gene_name"){
        title = paste("Multiplot ~ ", toupper(as.character(int.locus[[2]])), " ~ chr", snp.info$chr[1], sep="")
      } else if (int.locus[[1]] == "example"){
          title = paste("Multiplot ~ chr", snp.info$chr[1], ":", floor(min(snp.info$pos) + (max(snp.info$pos) - min(snp.info$pos))/2), sep = "")
      }
    } else if (type == "Manual scroll"){
      title <- paste("Multiplot ~ chr", snp.info$chr[1], sep="")
    }
  }
  return(title)
}

## plot in case of multiple files
function.multiPlot <- function(snp.info, list_for_loop, y.lim, type, plt.type, windows.number, smooth.par, int.locus,
                               col_list, genV, inpStrVar, recomb_yn, dotSize_yn){
  #print("1. preparing for plot!")
  # The first thing is to check if there are no association data for a gwas in the region of interest
  todel <- c()
  for (x in length(snp.info):1){
    if ((is.null(snp.info[[x]])) || (nrow(snp.info[[x]]) == 0)){ snp.info[[x]] <- NULL; todel <- c(todel, x) }
  }
  if (length(list_for_loop) >4){ ncol_leg = 4 } else { ncol_leg = length(list_for_loop) }
  pch_v = c()
  col_v = c()
  col_p = c()
  for (k in 1:length(list_for_loop)){
    if (k %in% todel){
      list_for_loop[[k]] <- paste0(list_for_loop[[k]], "=NA")
      pch_v <- c(pch_v, NA)
      col_v <- c(col_v, "white")
    } else {
      pch_v <- c(pch_v, 16)
      col_v <- c(col_v, col_list[[k]])
      col_p <- c(col_p, col_list[[k]])
    }
  }
  # First, let's prepare everything that is needed
  ## Prepare title
  title <- function.title(gwas = list_for_loop, int.locus, type, snp.info[[1]])
  ## Get recombination tracks
  sb <- snp.info[[1]]
  recomb <- findRecomb(sb, genV)
  ## Assign dot size
  for (x in 1:length(snp.info)){
    tmp <- function.pointSize(dat = snp.info[[x]], range = seq(2, 7, 0.5))
    snp.info[[x]] <- tmp
  }
  ## Also, in case the genome is hg38, now need to do the liftover on the snp.info
  if (genV == "GRCh38 (hg38)"){
    for (x in 1:length(snp.info)){
      tmp <- snp.info[[x]]
      df <- data.frame(chr=paste("chr", tmp$chr, sep=""), start=tmp$pos, end=tmp$pos+1)
      gr <- makeGRangesFromDataFrame(df)
      chain <- import.chain( "../data/databases/hg19ToHg38.over.chain")
      gr_hg38 <- liftOver(gr, chain)
      df_hg38 <- as.data.frame(gr_hg38)
      sb <- data.frame(chr=rep(tmp$chr[1], nrow(df_hg38)), pos=df_hg38$start, p=tmp[df_hg38$group, "p"], log10p=tmp[df_hg38$group, "-log10(P-value)"], size=tmp[df_hg38$group, "size"])
      colnames(sb) <- c("chr", "pos", "p", "-log10(P-value)", "size")
      snp.info[[x]] <- sb
    }
  }

  ## Select genes that are in the window to adjust y-axis
  genes <- function.dynamicGene(snp.info = snp.info[[1]], genV)
  if (nrow(genes) > 0){ min.y <- min(genes$y) } else { min.y <- 0 }

  ## layout is divided in 4 plots: sizes = plt1: 2, plt2: 1, plt3: 2, plt4: 2
  #layout(matrix(c(1,1, 2, 3,3, 4,4,4), nrow = 8, ncol = 1, byrow = TRUE))
  layout(matrix(c(1,1,1,1,1,1,
                  1,1,1,1,1,1,
                  2,2,2,2,2,2,
                  3,3,3,3,3,3,
                  3,3,3,3,3,3,
                  4,6,6,6,6,6,
                  5,7,7,7,7,7,
                  5,7,7,7,7,7,
                  5,7,7,7,7,7), nrow = 9, ncol = 6, byrow = T))

  ## set margins
  par(mar=c(0, 5, 4, 5))

  ## global parameter for the proportions
  pr <- 12
  #print("2. plot #1")

  ## empty plot as background
  plot(x = 0, y = 0, xlab='Chromosomal position (Mb)', cex.lab=2, xaxt='none', ylab="", ylim=c(0, y.lim), cex.axis = 1.5,
       pch=16, col="white", cex=2, type = "p", xaxs="i", yaxt='none', xlim=c(min(sb$pos), max(sb$pos)), main=title,
       cex.main=2.50, bty='n')

  ## add grid: 10 lines hor. and vert.
  for (x in seq(0, y.lim, (y.lim-min.y*y.lim/10)/10)){abline(h=x, lwd=0.6, col="grey80")}
  for (x in seq(min(sb$pos), max(sb$pos), (max(sb$pos)-min(sb$pos))/10)){ segments(x0 = x, y0 = 0, x1 = x, y1 = y.lim, col = "grey80", lwd=0.6) }

  # add recombination rates: for this, need to normalize between 0 and 100 the recombination rate
  recomb$norm.rate <- (y.lim - 1) * (recomb$"Rate(cM/Mb)" / 100)
  y.axis.recomb <- seq(0, 100, 25)
  y.axis.norm <- (y.lim - 1) * ((y.axis.recomb - min(y.axis.recomb))/(max(y.axis.recomb) - min(y.axis.recomb)))
  if (recomb_yn == "Yes"){
    points(recomb$"Position(bp)", recomb$norm.rate, type="l", lwd=1.5, col="darkolivegreen3")
  }

  # add significance lines and corresponding legend -- for now two at 0.05 and genome-wide 5e-8
  abline(h=-log10(0.05), lty=2, col=alpha("darkgreen", 1))
  abline(h=-log10(5e-8), lty=2, col=alpha("purple", 1))
  legend("topleft", bty='n', title = "Annotation lines", legend = c("p=0.05", "p=5e-8"), lty=c(2,2), lwd=c(2,2), col=c("darkgreen", "purple"), cex=1.25, ncol = 1)

  #manage axes
  axes <- seq(min(sb$pos), max(sb$pos), (max(sb$pos) - min(sb$pos))/7)
  axes.labels <- round(axes/1000000, 3)
  axes.x <- ceiling(seq(0, y.lim, y.lim/5))
  axis(side = 2, at = axes.x, labels = axes.x, cex.axis=2)

  #axis for recombination rates on the right
  if (recomb_yn == "Yes"){
    axis(side = 4, at = y.axis.norm, labels=seq(0, 100, 25), col='darkolivegreen3', col.axis = "darkolivegreen3", cex.axis=1.5, xpd=T)
    text(x = max(sb$pos), y = y.lim/5*4, "Recombination rate (cM/Mb)",srt = -90,
       col='darkolivegreen3', xpd=T, pos = 4, offset = 4, cex=2, font=2)
  }
  text(x = min(sb$pos), y = y.lim/3*2, "-Log10(P-value)", srt = 90,
       xpd=T, pos = 2, offset = 4, cex=2, font=2)

  #print("start density plot..")

  if (plt.type == "Points"){
    # then points
    for (x in 1:length(snp.info)){
      tmp <- snp.info[[x]]
      if (dotSize_yn == "No"){ tmp$size = 2 }
      tmp$"-log10(P-value)" <- -log10(tmp$p)
      points(x = tmp$pos, y = tmp$"-log10(P-value)",
             pch=16, col=alpha(col_p[x], 0.6), cex=tmp$size, type = "p", xaxs="i")
    }

    #if input type was a single snp (either position or rs id), then color the searched variant differently
    # snp.interest.flag <- 0
    # if (type %in% c("Position", "Rs ID")){
    #   #restrict to snp of interest, then plot it
    #   if (int.locus != "Type position..."){
    #     snp.interest.flag <- 1
    #     snp.interest <- snp.info[which(snp.info$locus == int.locus),]
    #     snp.interest.addF <- snp.info.f2[which(snp.info.f2$locus == int.locus),]
    #
    #     points(x = snp.interest$pos, y = -log10(snp.interest$p), pch=23, col="black", lwd=2, bg=alpha(colorPoint, 1), cex=snp.interest$size)
    #     points(x = snp.interest.addF$pos, y = -log10(as.numeric(snp.interest.addF$p)), pch=23, col="black", lwd=2, bg=alpha(colorPoint2, 1), cex=snp.interest.addF$size)
    #   }
    # }

    #add legend now
    # if (snp.interest.flag == 0){
    #   legend("topright", legend = c(lab1, lab2), col = c(colorPoint, colorPoint2), pch=16,
    #          cex=1.50, ncol=2, xpd=T, bty='n')
    # } else {
    #   legend("topright", legend = c(lab1, paste(lab1, " - Input", sep=""), lab2, paste(lab2, " - Input", sep="")),
    #          col = c(colorPoint, colorPoint, colorPoint2, colorPoint2), pch=c(16, 23, 16, 23), pt.lwd = 1.5,
    #          cex=1.50, ncol=2, xpd=T, bty='n')
    # }
    # add legend now
    legend("topright", legend=list_for_loop, pch=pch_v, col=col_v, cex=1.25, bty='n', ncol=ncol_leg)

    ###### HERE STARTS THE DENSITY PLOT
  } else {
    # main loop across all data to include, but first define output lists
    xl_list = list()
    pred_list = list()
    options(warn = -1)
    for (i in 1:length(snp.info)){
      #print(paste0("working on ", list_for_loop[i]))
      # get the data
      tmp <- snp.info[[i]]
      #print(head(tmp))
      # Sliding window approach -- then loess on sliding window values
      #print("ok before function")
      out <- function.DensityLinePvalue(snp.info = tmp, wind.n = windows.number)
      #print("ok after function")
      lo <- loess(out$pvalue ~ out$window, span=smooth.par)
      xl <- seq(min(out$window), max(out$window), (max(out$window) - min(out$window))/100)
      pred <- predict(lo, xl)
      pred[which(pred < 0)] <- 0
      # add limits -- left and right for polygon function
      xl <- c(xl[1], xl, xl[length(xl)])
      pred <- c(0, pred, 0)
      # finally put results into lists
      xl_list[[i]] <- xl
      pred_list[[i]] <- pred
      #print("ok until the end")
    }

    # then put the pvalue densities
    for (i in 1:length(xl_list)){ polygon(x = xl_list[[i]], y = pred_list[[i]], col = alpha(col_list[[i]], 0.4), border = col_list[[i]], lwd=3, xaxs="i", xpd=T) }
    options(warn = 0)

    #if input type was a single snp (either position or rs id), then color the searched variant differently -- here is a bar
    snp.interest.flag <- 0
    # if (type %in% c("Position", "Rs ID")){
    #   #restrict to snp of interest, then plot it
    #   if (int.locus != "Type position..."){
    #     snp.interest.flag <- 1
    #     snp.interest <- snp.info[which(snp.info$locus == int.locus),]
    #     snp.interest.add <- snp.info.f2[which(snp.info.f2$locus == int.locus),]
    #
    #     #need to grep the height in order to plot the variant -- the height is derived from the prediction
    #     chr.pos <- str_split_fixed(int.locus, ":", 2)
    #     pos.only <- as.numeric(chr.pos[, 2])
    #     pred.pos <- predict(lo, pos.only)
    #     pred.pos.add <- predict(lo.add, pos.only)
    #
    #     #add to plot
    #     points(x = snp.interest.add$pos, y = pred.pos.add, type="h", col=alpha("pink", 0.8), lwd=4, xaxs="i")
    #     points(x = snp.interest$pos, y = pred.pos, type="h", col=alpha("yellow", 0.8), lwd=4, xaxs="i")
    #
    #   }
    # }

    #add legend now
    # if (snp.interest.flag == 0){
    #   legend("topright", legend = c(lab1, lab2), col = c(colorPoint, colorPoint2), pch=16,
    #          cex=1.50, ncol=2, xpd=T, bty='n')
    # } else {
    #   legend("topright", legend = c(lab1, paste(lab1, " - Input", sep=""), lab2, paste(lab2, " - Input", sep="")),
    #          col = c(colorPoint, "lightblue", colorPoint2, "yellow"), pch=c(16, 23, 16, 23), pt.lwd = 1.5,
    #          cex=1.50, ncol=2, xpd=T, bty='n')
    # }
    legend("topright", legend=list_for_loop, pch=pch_v, col=col_v, cex=1.25, bty='n', ncol=ncol_leg)
  }
  #print("3. plot #2")
  #############################
  # HERE IS FOR THE SECOND PLOT -- GENE ANNOTATION
  ## global parameter for the proportions
  pr <- 12
  ## empty plot as background
  par(mar=c(5, 5, 0, 5))
  ymax <- 3
  if (nrow(genes) >ymax){ genes$y <- abs(genes$y); ymax <- max(genes$y)+0.5 }
  plot(x = 0, y = 0, xlab='Chromosomal position (Mb)', cex.lab=2, xaxt='none', ylab="", ylim=c(0, ymax), cex.axis = 1.5,
       pch=16, col="white", cex=2, type = "p", xaxs="i", yaxt='none', xlim=c(min(sb$pos), max(sb$pos)),
       cex.main=2.50, bty='n', yaxt='none')

  text(x = min(sb$pos), y = ymax/5*4, "Gene track", srt = 90,
       xpd=T, pos = 2, offset = 4, cex=2, font=2)

  #manage axes
  axes <- seq(min(sb$pos), max(sb$pos), (max(sb$pos) - min(sb$pos))/7)
  axes.labels <- round(axes/1000000, 3)
  axis(side = 1, at=axes, cex.axis=2, labels=axes.labels)

  # genes
  genes$y <- genes$y-0.5
  if (nrow(genes) >35){
    genes$cex <- 0.40; genes_lwd=0.6; tmp_pr = 0.02
  } else if (nrow(genes) >20){
    genes$cex <- 0.75; genes_lwd=0.6; tmp_pr = 0.03
  } else if (nrow(genes) >15){
    genes$cex <- 0.90; genes_lwd = 1; tmp_pr = 0.04
  } else if (nrow(genes) >10){
    genes$cex <- 1.2; genes_lwd = 1.5; tmp_pr = 0.05
  } else if (nrow(genes) <=5){
    genes$cex <- 2; genes_lwd = 2.5; tmp_pr = 0.085
  } else {
    genes$cex <- 1.5; genes_lwd = 2; tmp_pr = 0.07
  }

  if (min.y != 0){
    genes$y <- abs(genes$y)
    #manage gene names
    for (g in 1:nrow(genes)){
      #main gene line -- full transcription sequence
      segments(x0 = genes$txStart[g], y0 = genes$y[g], x1 = genes$txEnd[g], y1 = genes$y[g], lwd=genes_lwd)
      #need to divide exones from introns
      start <- as.data.frame(t(str_split_fixed(genes$exonStarts[g], ",", genes$exonCount[g]+1)))
      end <- as.data.frame(t(str_split_fixed(genes$exonEnds[g], ",", genes$exonCount[g]+1)))
      exons <- cbind(start, end)
      colnames(exons) <- c("start", "end")
      exons$start <- as.numeric(as.character(exons$start))
      exons$end <- as.numeric(as.character(exons$end))
      #main loop over exons
      for (j in 1:nrow(exons)){
        rect(xleft=exons$start[j], ybottom=genes$y[g]-(ymax/4*(tmp_pr)), xright=exons$end[j],
             ytop = genes$y[g]+(ymax/4*(tmp_pr)), col='grey80', lwd=0.80)
      }
      text(x = genes$txStart[g] + (genes$txEnd[g] - genes$txStart[g])/2, y = genes$y[g]+ymax*tmp_pr,
           labels=genes$"#geneName"[g], font=3, cex=genes$cex[g])
      #        if (genes$strand[g] == "+"){
      #          arrows(x0 = genes$txEnd[g], y0 = genes$y[g], x1 = genes$txEnd[g] + (10*2/100), y1 = genes$y[g],
      #                 length=0.1, lwd=2, col='coral')
      #        } else if (genes$strand[g] == "-"){
      #          arrows(x0 = genes$txStart[g], y0 = genes$y[g], x1 = genes$txStart[g] - (10*2/100), y1 = genes$y[g],
      #                 length=0.1, lwd=2, col='coral')
      #        }
    }
  }
  #print("4. plot #3")
  ################################
    # HERE IS FOR THE THIRD PLOT -- STRUCTURAL VARIATIONS
  par(mar=c(5, 5, 7, 5))
  plot(x = 0, y = 0, xlab='Chromosomal position (Mb)', cex.lab=2, xaxt='none',
         ylab="", ylim=c(0, 10), cex.axis = 2, pch=16, col="white", cex=2, type = "p", xaxs="i",
         yaxt='none', xlim=c(min(sb$pos), max(sb$pos)), cex.main=2.50, bty='n')

    #add grid
    for (x in seq(0, 10, 2)) {abline(h=x, lwd=0.6, col="grey80")}
    for (x in seq(min(sb$pos), max(sb$pos), (max(sb$pos)-min(sb$pos))/10)){
      segments(x0 = x, y0 = 0, x1 = x, y1 = 10, col = "grey80", lwd=0.6)
    }

    # add label on y-axis
    text(x = min(sb$pos), y = 8, "Structural variations", srt = 90,
         xpd=T, pos = 2, offset = 4, cex=2, font=2)

    #manage axes
    axes <- seq(min(sb$pos), max(sb$pos), (max(sb$pos) - min(sb$pos))/7)
    axes.labels <- round(axes/1000000, 3)
    axis(side = 1, at=axes, cex.axis=2, labels=axes.labels)
    axes.x <- ceiling(seq(0, y.lim, y.lim/5))

    #take structural variants in the window
    strVar <- findStr_variats(sb, genV, inpStrVar)

    #plot as segments (only if there are variants to plot)
    # plot as segments (only if there are variants to plot)
    if (nrow(strVar) >15){ strVar$cex <- 1.40 } else { strVar$cex <- 1.80 }
    if (nrow(strVar) > 0){
      for (i in 1:nrow(strVar)){
        segments(x0 = strVar$start_pos[i], y0 = strVar$y_plot[i], x1 = strVar$end_pos[i],
                 y1 = strVar$y_plot[i], lwd = 3, col=strVar$col[i])
        segments(x0 = strVar$start_pos[i], y0 = strVar$y_plot[i]-0.15, x1 = strVar$start_pos[i],
                 y1 = strVar$y_plot[i]+0.15, lwd = 3, col=strVar$col[i])
        segments(x0 = strVar$end_pos[i], y0 = strVar$y_plot[i]-0.15, x1 = strVar$end_pos[i],
                 y1 = strVar$y_plot[i]+0.15, lwd = 3, col=strVar$col[i])
        if (nrow(strVar) <30){
          text(x = strVar$middle[i], y = strVar$y_plot[i]+1, labels = strVar$diff_alleles[i],
               cex = strVar$cex[i], font=3, xpd=T, col=strVar$col[i])
        }
        tmp_lg <- strVar[!duplicated(strVar$type),]
      }
      legend(x = min(sb$pos), y = 13, legend = tmp_lg$type, lty = 1, col=tmp_lg$col, xpd=T, cex=1.50, bty='n', lwd=4, ncol=4)
    } else {
      text(x = min(sb$pos) + (max(sb$pos) - min(sb$pos))/2, y = 5,
           labels = "There are no Structural variants to be plotted here", cex=2, adj=0.5)
    }
    #print("5. plot #4")

    # THIRD PLOT -- GENE EXPRESSION FROM GTEX
    # empty plot first
    par(mar=c(4, 2, 10, 0))
    plot(0, 0, xlim=c(0, 55), ylim=c(0, 25), xaxt='none', yaxt='none', xlab="", ylab="",
         main="", cex.main=2.50, col="white", bty="n")
    # run function to take genes (if there are genes)
    if (nrow(genes) > 0){
      gtex.info <- findExpr_gtex(genes, gtex.db)
      heat <- as.data.frame(gtex.info[[1]])
      g <- gtex.info[[2]]
      if (is.na(heat) & is.na(g)){
        par(mar=c(4, 8, 10, 6))
        plot(0, 0, xlim=c(0, 55), ylim=c(0, 25), bty="n", xaxt='none', yaxt='none', xlab="", ylab="",
             main="PLOT 4", cex.main=2.50, col="white")
        plot(0, 0, xlim=c(0, 55), ylim=c(0, 25), bty="n", xaxt='none', yaxt='none', xlab="", ylab="",
             main="RNA expression from GTEx", cex.main=2.50, col="white")
        text(x = 27.5, y = 25, labels = "No GTEx information for the genes in this region.",
             cex = 2, xpd=T)
      } else {

        # then plot as heatmap
        colors <- viridis(n = 101, option = "plasma")

        # scale heat between 0 and 1
        raw <- heat
        heat <- (heat - min(heat))/(max(heat) - min(heat))

        # dendrogram for the genes
        mx <- 22
        if (nrow(heat) >1){
          hr = as.dendrogram(hclust(d = dist(heat, method = "euclidean"), method = "ward.D2"))
          hc = hclust(d = dist(heat, method = "euclidean"), method = "ward.D2")
          ordered_labels = g[hc$order]
          par(mar=c(4, 2, 0, 0))
          nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), cex = 1.25, col = "navy")
          plot(hr, horiz = T, xaxt="none", yaxt ="none", nodePar = nodePar, leaflab = "none", ylim=c(-mx+nrow(heat)+1.5, nrow(heat)+1.5))
        } else {
          plot(0, 0, xaxt="none", yaxt ="none", ylim=c(0, 1), xlim =c(0,1), main="", xlab="", ylab="", bty="n", pch=16, col="white")
          text(x = 0.5, y = 1, labels = "Too few genes\nfor clustering", xpd=T)
          ordered_labels <- g
        }

        # dendrogram for the tissues
        tr = t(heat)
        hr = as.dendrogram(hclust(dist(tr, method = "euclidean"), method = "ward.D2"))
        hc = hclust(dist(tr, method = "euclidean"), method = "ward.D2")
        ordered_labels_tissues = rownames(tr)[hc$order]
        par(mar=c(0, 0, 10, 6))
        nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), cex = 1.25, col = "red")
        plot(hr, horiz = F, xaxt="none", yaxt ="none", nodePar = nodePar, leaflab = "none", xlim=c(0.5, 55.5), main="RNA expression from GTEx", cex.main=2.50)
        # legend of the heatmap goes here
        mx_height = max(hc$height)
        gradient.rect(xleft = 55.5, ybottom = mx_height, xright = 60, ytop = mx_height+mx_height*0.15, col = colors)
        leg.lab <- c(-1, 0, 1)
        pos <- c(55.5, 57.75, 60)
        for (i in 1:length(pos)){
          text(x = pos[i], y = (mx_height+mx_height*0.25), labels = leg.lab[i], cex=1, xpd=T, font=2)
        }

        # background
        par(mar=c(4, 0, 0, 6))
        # reorder the genes to match dendrogram order
        rownames(heat) <- g
        heat = heat[match(ordered_labels, rownames(heat)),]
        tmp = as.data.frame(t(heat))
        tmp = tmp[match(ordered_labels_tissues, rownames(tmp)),]
        heat <- as.data.frame(t(tmp))

        plot(0, 0, xlim=c(0, 55), ylim=c(0, mx), bty="n", xaxt='none', yaxt='none', xlab="", ylab="",
             main="", cex.main=2.50, col="white")
        heat <- as.matrix(heat)
        # loop
        for (i in 1:nrow(heat)){
          for (j in 1:ncol(heat)){
            # value
            v <- round(as.numeric(heat[i, j])*100)
            rect(xleft = j-1, ybottom = mx-i-1, xright = j, ytop = mx-i, col = colors[v+1],
                 border = "white", lwd = 1.5)
          }
          text(x = 54, y = mx-i-0.5, labels = g[i], cex=1.50, font=4, xpd=T, adj=0)
        }
        for (j in 1:ncol(heat)){
          text(x = j-0.5, y = mx-i-1.5, labels = ordered_labels_tissues[j], cex=1, font=2, adj=1, xpd=T, srt=90)
        }
      }
    } else {
      plot(0, 0, xlim=c(0, 55), ylim=c(0, 25), bty="n", xaxt='none', yaxt='none', xlab="", ylab="",
           main="RNA expression from GTEx", cex.main=2.50, col="white")
      text(x = 27.5, y = 25, labels = "No GTEx information for the genes in this region.",
           cex = 2, xpd=T)
    }
}

## function to check if data is from the required chromosome
function.catchChromosome <- function(chr, gwas, res_example){
  if (gwas == "example"){
    #take path
    chr_toMatch <- as.numeric(chr)-15

    #read data
    dat <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 3))
    try(dat <- res_example[[chr_toMatch]], silent = T)
    colnames(dat) <- c("chr", "pos", "p")
    dat <- dat[!which(is.na(dat$p)),]
    chrom = dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$"-log10(P-value)" <- -log10(dat$p)

  } else {
    #in these cases, the GWAS is usually the name of the folder
    #take path
    fname = paste("../data/", gwas, "/chr", as.character(chr), "_", gwas, ".txt.gz", sep="")

    #read data
    dat <- fread(fname, h=T, stringsAsFactors = F)
    colnames(dat) <- c('chr', 'pos', 'p')
    dat$p <- as.numeric(dat$p)
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$'-log10(P-value)' <- -log10(dat$p)
  }
  return(dat)
}

## function to manage position as input type
function.InputPos <- function(dat, window, snp_locus, gwas){
  # split locus to find close positions and define intervals where to find snps
  chr <- snp_locus$chr
  pos.max <- snp_locus$pos + window
  pos.min <- snp_locus$pos - window

  # check if data is relative to that chromosome, and in case change it
  if (!(gwas %in% c("loaded", "toLoad"))){ dat <- function.catchChromosome(chr, gwas, res_example) }

  # take data of interest
  snp.info <- dat[which(dat$pos >= pos.min & dat$pos <= pos.max), ]

  return(snp.info)
}

## function to manage manual scroll as input type
function.InputManualScroll <- function(dat, window, input.scroll, input.chrom, gwas){
  #print(gwas)
  #initial position is all variable
  pos.init <- as.numeric(input.scroll)

  #define limits
  lower <- pos.init - window
  upper <- pos.init + window

  #check if data is relative to that chromosome, and in case change it
  if (!(gwas %in% c("loaded", "toLoad"))){ dat <- function.catchChromosome(input.chrom, gwas, res_example) }

  #take data of interest
  snp.info <- dat[which(dat$pos >= lower & dat$pos <= upper), ]
  return(snp.info)
}

## Function to manage genes as input
function.InputGenes <- function(gene){
  if (gene != "Type gene symbol..."){
    # make sure the gene input is uppercase
    gene <- toupper(gene)
    gene <- paste("^", gene, "$", sep="")

    # identify gene coordinates
    gene.info <- gene.db[grep(gene, gene.db$"#geneName"),]
    gene.info <- gene.info[order(-gene.info$txEnd),]
    gene.info <- gene.info[!duplicated(gene.info$"#geneName"),]

    # check if the gene actually exists
    if (nrow(gene.info) == 0){ gene.info <- "gene_not_in_list" }

  } else {
    #if there is no input gene, then take a random integer (now in chr19) and plot a random gene
    gene.db <- gene.db[which(gene.db$chrom == "chr21"),]
    gene.info <- gene.db[ceiling(nrow(gene.db)/4), ]
  }

  return(gene.info)
}

## Function to find position in case input is rsid
function.rsIDasInput <- function(target){
  # set background snps
  snps_hg19 <- SNPlocs.Hsapiens.dbSNP144.GRCh37
  # find snp information in hg19
  snp_info <- NULL
  try(snp_info <- as.data.frame(snpsById(snps_hg19, target[[2]])), silent = T)
  if (!is.null(snp_info)){
    snp_info <- data.frame(chr=snp_info$seqnames, pos=snp_info$pos)
  } else {
    snp_info <- "snp_not_in_list"
  }
  return(snp_info)
}

## function to find gwas hits from inputted gene
function.GWASfromGene <- function(dat, gene.info, window, gwas){
  # gene information
  chr <- gene.info$chrom
  chr.n <- str_split_fixed(chr, "chr", 2)[, 2]
  start <- gene.info$txStart - window
  end <- gene.info$txEnd + window

  # check if data is relative to that chromosome, and in case change it
  if (!(gwas %in% c("loaded", "toLoad"))){ dat <- function.catchChromosome(chr.n, gwas, res_example) }

  # define snps of interest
  snp.info <- dat[which(dat$pos >= start & dat$pos <= end),]

  return(snp.info)
}

## function to find density line of p-value across chromosomal position
function.DensityLinePvalue <- function(snp.info, wind.n){
  #define min and max
  min.wind <- min(snp.info$pos)
  max.wind <- max(snp.info$pos)

  #defining sliding window size
  interval <- ceiling((max.wind - min.wind)/wind.n)

  #define output
  out = matrix(data = NA, nrow = wind.n, ncol = 3)
  colnames(out) <- c("window", "pvalue", "error")

  #counter for output assignement
  counter <- 1

  #print("starting loop")
  #loop using sliding window
  #print(seq(min.wind, max.wind, interval))
  for (i in seq(min.wind, max.wind, interval)){
    if (counter <= wind.n){
      #define internal maximum -- f
      f <- i + interval
      #take info in the window of interest
      data.sbs <- snp.info[which((snp.info$pos >= i) & (snp.info$pos <= f)),]
      data.sbs$p <- as.numeric(data.sbs$p)
      data.sbs <- data.sbs[order(data.sbs$p),]
      #check number of rows
      if (nrow(data.sbs) > 0){
        #take maximum pvalue point
        out[counter, ] <- c(((f - i)/2) + i, -log10(data.sbs$p[1]), sd(na.omit(data.sbs$p)))
      } else {
        out[counter, ] <- c(((f - i)/2) + i, 0, 0)
      }
      counter <- counter + 1
    }
  }
  #remove NA and substitute them with 0
  out[which(is.na(out[, "error"]) == TRUE), "error"] <- 0
  out <- as.data.frame(out)

  return(out)
}

## function to parse gene location file and grep only region of interest
function.dynamicGene <- function(snp.info, genV){
  #define searching space for genes
  min.x <- min(snp.info$pos)
  max.x <- max(snp.info$pos)

  if (genV == "GRCh38 (hg38)"){
    #find genes in interval (min.x -- max.x) -- then clean up a bit
    gene.loc.res <- subset(genes.hg38, genes.hg38$chrom == paste("chr", snp.info$chr[1], sep=""))
    gene.loc.res$in.int <- 0
    gene.loc.res$in.int[which((gene.loc.res$txStart >= min.x) & (gene.loc.res$txEnd <= max.x))] <- 1
    gene.loc.res$in.int[which((gene.loc.res$txStart <= min.x) & (gene.loc.res$txEnd >= min.x))] <- 1
    gene.loc.res$in.int[which((gene.loc.res$txStart <= max.x) & (gene.loc.res$txEnd >= max.x))] <- 1
    genes <- gene.loc.res[which(gene.loc.res$in.int == 1),]
    genes <- genes[order(-genes$txEnd),]
    genes <- genes[!duplicated(genes$"#geneName"),]
  } else {
    #find genes in interval (min.x -- max.x) -- then clean up a bit
    gene.loc.res <- subset(gene.db, gene.db$chrom == paste("chr", snp.info$chr[1], sep=""))
    gene.loc.res$in.int <- 0
    gene.loc.res$in.int[which((gene.loc.res$txStart >= min.x) & (gene.loc.res$txEnd <= max.x))] <- 1
    gene.loc.res$in.int[which((gene.loc.res$txStart <= min.x) & (gene.loc.res$txEnd >= min.x))] <- 1
    gene.loc.res$in.int[which((gene.loc.res$txStart <= max.x) & (gene.loc.res$txEnd >= max.x))] <- 1
    genes <- gene.loc.res[which(gene.loc.res$in.int == 1),]
    genes <- genes[order(-genes$txEnd),]
    genes <- genes[!duplicated(genes$"#geneName"),]
  }

  #define position in the plot
  n = ceiling(nrow(genes)/2)
  if (n != 0){
    v <- -1
    genes$y <- NA
    for (i in 1:nrow(genes)){
      genes$y[i] <- v
      v <- v-1
      if (v == -n-1){ v = -1 }
    }

  }
  return(genes)
}

## Function to manage which input data to plot
function.manageInput <- function(inp, supp_f, res_example){
  # Example as input
  if (length(inp) == 0){
    dat <- res_example[[1]]
    colnames(dat) <- c("chr", "pos", "p")
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom = dat$chr[1]
    dat$"-log10(P-value)" <- -log10(as.numeric(dat$p))
    gwas = "example"

    # User file as input
  } else if ("toLoad" %in% inp){
    # read file
    dat <- fread(supp_f, h=T, stringsAsFactors = F)

    # analyse header and slim data
    dat = identiHeader(dat)
    dat$p <- as.numeric(dat$p)
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$"-log10(P-value)" <- -log10(as.numeric(dat$p))
    gwas = "loaded"
    #print(head(dat))

  } else if (inp == 'IGAP') {
    dat <- fread('../data/IGAP/chr19_IGAP.txt.gz', h=T, stringsAsFactors=F)
    colnames(dat) <- c('chr', 'pos', 'p')
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$'-log10(P-value)' <- -log10(dat$p)
    gwas <- 'IGAP'
  } else if (inp == 'CAD') {
    dat <- fread('../data/CAD/chr19_CAD.txt.gz', h=F, stringsAsFactors = F)
    colnames(dat) <- c('chr', 'pos', 'p')
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$'-log10(P-value)' <- -log10(dat$p)
    gwas <- 'CAD'
  } else if (inp == 'CAD_Diabetics') {
    dat <- fread('../data/CAD_Diabetics/chr19_CAD_Diabetics.txt.gz', h=T, stringsAsFactors=F)
    colnames(dat) <- c('chr', 'pos', 'p')
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$'-log10(P-value)' <- -log10(dat$p)
    gwas <- 'CAD_Diabetics'
  } else if (inp == '100plus') {
    dat <- fread('../data/100plus/chr19_100plus.txt.gz', h=T, stringsAsFactors=F)
    colnames(dat) <- c('chr', 'pos', 'p')
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$'-log10(P-value)' <- -log10(dat$p)
    gwas <- '100plus'
  } else if (inp == 'GR@ACE') {
    dat <- fread('../data/GR@ACE/chr19_GR@ACE.txt.gz', h=T, stringsAsFactors=F)
    colnames(dat) <- c('chr', 'pos', 'p')
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$'-log10(P-value)' <- -log10(dat$p)
    gwas <- 'GR@ACE'
  } else if (inp == 'UKBaging') {
    dat <- fread('../data/UKBaging/chr19_UKBaging.txt.gz', h=T, stringsAsFactors=F)
    colnames(dat) <- c('chr', 'pos', 'p')
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$'-log10(P-value)' <- -log10(dat$p)
    gwas <- 'UKBaging'
  } else if (inp == 'exomeADES') {
    dat <- fread('../data/exomeADES/chr19_exomeADES.txt.gz', h=T, stringsAsFactors=F)
    colnames(dat) <- c('chr', 'pos', 'p')
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$'-log10(P-value)' <- -log10(dat$p)
    gwas <- 'exomeADES'
  } else if (inp == 'LVV') {
    dat <- fread('../data/LVV/chr19_LVV.txt.gz', h=T, stringsAsFactors=F)
    colnames(dat) <- c('chr', 'pos', 'p')
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$'-log10(P-value)' <- -log10(dat$p)
    gwas <- 'LVV'
  } else if (inp == 'EADB_AD') {
    dat <- fread('../data/EADB_AD/chr19_EADB_AD.txt.gz', h=T, stringsAsFactors=F)
    dat <- dat[, c("chr", "pos", "p")]
    colnames(dat) <- c('chr', 'pos', 'p')
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$'-log10(P-value)' <- -log10(dat$p)
    gwas <- 'EADB_AD'
  } else if (inp == 'COVID') {
    dat <- fread('../data/COVID/chr19_COVID.txt.gz', h=T, stringsAsFactors=F)
    colnames(dat) <- c('chr', 'pos', 'p')
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$'-log10(P-value)' <- -log10(dat$p)
    gwas <- 'COVID'
  } else if (inp == "Height") {
    dat <- fread('../data/Height/chr19_Height.txt.gz', h=T, stringsAsFactors=F)
    colnames(dat) <- c('chr', 'pos', 'p')
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$'-log10(P-value)' <- -log10(dat$p)
    gwas <- 'Height'
  } else if (inp == "SBP") {
    dat <- fread('../data/SBP/chr19_SBP.txt.gz', h=T, stringsAsFactors=F)
    colnames(dat) <- c('chr', 'pos', 'p')
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$'-log10(P-value)' <- -log10(dat$p)
    gwas <- 'SBP'
  } else if (inp == "BMI") {
    dat <- fread('../data/BMI/chr19_BMI.txt.gz', h=T, stringsAsFactors=F)
    colnames(dat) <- c('chr', 'pos', 'p')
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$'-log10(P-value)' <- -log10(dat$p)
    gwas <- 'BMI'
  } else if (inp == "Breast_cancer") {
    dat <- fread('../data/Breast_cancer/chr19_Breast_cancer.txt.gz', h=T, stringsAsFactors=F)
    colnames(dat) <- c('chr', 'pos', 'p')
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$'-log10(P-value)' <- -log10(dat$p)
    gwas <- 'Breast_cancer'
  } else if (inp == "proxy_AD") {
    dat <- fread('../data/proxy_AD/chr19_proxy_AD.txt.gz', h=T, stringsAsFactors=F)
    colnames(dat) <- c('chr', 'pos', 'p')
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$'-log10(P-value)' <- -log10(dat$p)
    gwas <- 'proxy_AD'
  } else if (inp == "Lupus") {
    dat = fread("../data/Lupus/chr19_Lupus.txt.gz", h=F, stringsAsFactors = F)
    colnames(dat) <- c("chr", "pos", "p")
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$"-log10(P-value)" <- -log10(dat$p)
    gwas = "Lupus"
  } else if (inp == "Education") {
    dat = fread("../data/Education/chr19_Education.txt.gz", h=F, stringsAsFactors = F)
    colnames(dat) <- c("chr", "pos", "p")
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$"-log10(P-value)" <- -log10(dat$p)
    gwas = "Education"
  } else if (inp == "Inflammation"){
    dat = fread("../data/Inflammation/chr19_Inflammation.txt.gz", h=F, stringsAsFactors = F)
    colnames(dat) <- c("chr", "pos", "p")
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$"-log10(P-value)" <- -log10(dat$p)
    gwas = "Inflammation"
  } else if (inp == "Myeloproliferative"){
    dat = fread("../data/Myeloproliferative/chr19_Myeloproliferative.txt.gz", h=F, stringsAsFactors = F)
    colnames(dat) <- c("chr", "pos", "p")
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$"-log10(P-value)" <- -log10(dat$p)
    gwas = "Myeloproliferative"
  } else if (inp == "Bone_density"){
    dat = fread("../data/Bone_density/chr19_Bone_density.txt.gz", h=F, stringsAsFactors = F)
    colnames(dat) <- c("chr", "pos", "p")
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$"-log10(P-value)" <- -log10(dat$p)
    gwas = "Bone_density"
  } else if (inp == "Prostate"){
    dat = fread("../data/Prostate/chr19_Prostate.txt.gz", h=F, stringsAsFactors = F)
    colnames(dat) <- c("chr", "pos", "p")
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$"-log10(P-value)" <- -log10(dat$p)
    gwas = "Prostate"
  } else if (inp == "Lung"){
    dat = fread("../data/Lung/chr19_Lung.txt.gz", h=F, stringsAsFactors = F)
    colnames(dat) <- c("chr", "pos", "p")
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$"-log10(P-value)" <- -log10(dat$p)
    gwas = "Lung"
  } else if (inp == "Leukemia"){
    dat = fread("../data/Leukemia/chr19_Leukemia.txt.gz", h=F, stringsAsFactors = F)
    colnames(dat) <- c("chr", "pos", "p")
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$"-log10(P-value)" <- -log10(dat$p)
    gwas = "Leukemia"
  } else if (inp == "Autism"){
    dat = fread("../data/Autism/chr19_Autism.txt.gz", h=F, stringsAsFactors = F)
    colnames(dat) <- c("chr", "pos", "p")
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$"-log10(P-value)" <- -log10(dat$p)
    gwas = "Autism"
  } else if (inp == "Asthma"){
    dat = fread("../data/Asthma/chr19_Asthma.txt.gz", h=F, stringsAsFactors = F)
    colnames(dat) <- c("chr", "pos", "p")
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$"-log10(P-value)" <- -log10(dat$p)
    gwas = "Asthma"
  } else if (inp == "Depression"){
    dat = fread("../data/Depression/chr19_Depression.txt.gz", h=F, stringsAsFactors = F)
    colnames(dat) <- c("chr", "pos", "p")
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$"-log10(P-value)" <- -log10(dat$p)
    gwas = "Depression"
  } else if (inp == "Vitamin_D"){
    dat = fread("../data/Vitamin_D/chr19_Vitamin_D.txt.gz", h=F, stringsAsFactors = F)
    colnames(dat) <- c("chr", "pos", "p")
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$"-log10(P-value)" <- -log10(dat$p)
    gwas = "Vitamin_D"
  } else if (inp == "Diabetes"){
    dat = fread("../data/Diabetes/chr19_Diabetes.txt.gz", h=F, stringsAsFactors = F)
    colnames(dat) <- c("chr", "pos", "p")
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$"-log10(P-value)" <- -log10(dat$p)
    gwas = "Diabetes"
  } else if (inp == "pTAU"){
    dat = fread("../data/pTAU/chr19_pTAU.txt.gz", h=F, stringsAsFactors = F)
    colnames(dat) <- c("chr", "pos", "p")
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$"-log10(P-value)" <- -log10(dat$p)
    gwas = "pTAU"
  } else if (inp == "Multivariate_Longevity"){
    dat = fread("../data/Multivariate_Longevity/chr19_Multivariate_Longevity.txt.gz", h=F, stringsAsFactors = F)
    colnames(dat) <- c("chr", "pos", "p")
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$"-log10(P-value)" <- -log10(dat$p)
    gwas = "Multivariate_Longevity"
  } else if (inp == "Alzheimer_million"){
    dat = fread("../data/Alzheimer_million/chr19_Alzheimer_million.txt.gz", h=F, stringsAsFactors = F)
    colnames(dat) <- c("chr", "pos", "p")
    dat <- dat[!which(is.na(dat$p)),]
    dat$p[which(dat$p == 0)] <- 0.00000001
    chrom <- dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$"-log10(P-value)" <- -log10(dat$p)
    gwas = "Alzheimer_million"
  }#else here to add repositories

  return(list(dat, chrom, gwas))
}

## function to identify header in uploaded file
identiHeader <- function(dat){
  head <- toupper(names(dat))
  #find chromosome, position and pvalue
  chrom.idx <- NULL
  pos.idx <- NULL
  p.idx <- NULL
  p.kw <- c("P", "PVAL", "P_VALUE", "PVALUE", "P_VAL")
  chrom.idx = grep("chr", head, ignore.case = T)
  pos.idx = grep("pos", head, ignore.case = T)
  if (length(pos.idx) == 0){pos.idx = grep("bp", head, ignore.case = T)}
  for (x in 1:length(head)){if (head[x] %in% p.kw){ p.idx <- x } }

  #check for correctness
  run = FALSE
  if (length(chrom.idx) + length(pos.idx) + length(p.idx) == 3){
    run <- TRUE
    #print("## header correctly read")

    #now rename columns
    head[chrom.idx] <- "chr"
    head[pos.idx] <- "pos"
    head[p.idx] <- "p"
    colnames(dat) <- head
    dat <- dat[, c("chr", "pos", "p")]
  } else {
    #print("!!!There was a problem with the input!!!!")
  }
  return(dat)
}

## Function to identify target input from the user
identiTargetRegion <- function(region){
  # Check if it is a locus (chr:pos)
  if (region == "Type locus, rsID or gene name"){
    target <- list("example")
  } else if (grepl(":", region) == TRUE){
    tmp <- str_split_fixed(region, ":", 2)
    target <- list("locus", as.numeric(tmp[, 1]), as.numeric(tmp[, 2]))
  } else if (grepl("rs", region) == TRUE){
    target <- list("rsid", region)
  } else {
    target <- list("gene_name", region)
  }
  return(target)
}

## Function to output the table of top 10 associations
plotTable <- function(input){
  ############################
  #create flag for uploaded input file
  inFile <- input$inp_f
  path_f <- "None"
  if (!is.null(inFile)){
    path_f = inFile$datapath
    #print("file loaded")
  } else {
    path_f <- "None"
  }
  ##########################

  #############################
  all_gwas <- c(input$gwas_neuro, input$gwas_cardio, input$gwas_immune, input$gwas_cancer, input$gwas_physio, input$gwas_from_file)
  res = function.manageInput(as.character(all_gwas[1]), path_f, res_example)
  dat = as.data.frame(res[[1]])
  chrom = as.numeric(res[[2]])
  gwas = as.character(res[[3]])
  #############################
  if (input$sel == "Locus, Gene, RsID"){
    # If input is a target region, need to identify which was the input (locus, gene or rsid)
    target <- identiTargetRegion(input$target)
    target_type <- target[[1]]

    # check number of GWAS to be plotted
    all_gwas <- c(input$gwas_neuro, input$gwas_cardio, input$gwas_immune, input$gwas_cancer, input$gwas_physio, input$gwas_from_file)
    if (length(all_gwas) < 2){
      # check input type and extract snp information and gwas info accordingly
      # example case
      if (target_type == "example"){
        snp_locus = data.frame(chr=16, pos=dat[ceiling(nrow(dat)/4), "pos"])
        res <- function.InputPos(dat = dat, window = input$x, snp_locus, gwas)

        # rsid case
      } else if (target_type == "rsid"){
        snp_locus <- function.rsIDasInput(target)
        res <- function.InputPos(dat = dat, window = input$x, snp_locus, gwas)

        # chromosome:position case
      } else if (target_type == "locus") {
        snp_locus <- data.frame(chr=target[[2]], pos=target[[3]])
        res <- function.InputPos(dat = dat, window = input$x, snp_locus, gwas)

        # gene_name case
      } else if (target_type == "gene_name"){
        gene.info <- function.InputGenes(gene = target[[2]])
        res <- function.GWASfromGene(dat = dat, gene.info = gene.info, window = input$x, gwas)
      }

      # assign to toplot object
      toplot <- res

      # The following in case multiple gwas need to be plotted
    } else {
      #check whether uploaded file is among those to show
      if (path_f == "None"){
        # loop over the gwas to plot
        # define some outputs -- the first is the names of the gwases
        list_for_loop <- all_gwas
        # the second is the color for each gwas
        col_list <- list()
        for (gw in 1:length(list_for_loop)){
          if (gw == 1){
            # check input type and extract snp information and gwas info accordingly
            # example case
            if (target_type == "example"){
              snp_locus = data.frame(chr=16, pos=dat[ceiling(nrow(dat)/4), "pos"])
              res <- function.InputPos(dat = dat, window = input$x, snp_locus, gwas)

              # rsid case
            } else if (target_type == "rsid"){
              snp_locus <- function.rsIDasInput(target)
              res <- function.InputPos(dat = dat, window = input$x, snp_locus, gwas)

              # chromosome:position case
            } else if (target_type == "locus") {
              snp_locus <- data.frame(chr=target[[2]], pos=target[[3]])
              res <- function.InputPos(dat = dat, window = input$x, snp_locus, gwas)

              # gene_name case
            } else if (target_type == "gene_name"){
              gene.info <- function.InputGenes(gene = target[[2]])
              res <- function.GWASfromGene(dat = dat, gene.info = gene.info, window = input$x, gwas)
            }
            res$Study = gwas
            list_for_results <- res
            col_list[[gw]] <- input$col
          } else {
            # other datasets
            tmp <- as.character(all_gwas[gw])
            res = function.manageInput(tmp, path_f, res_example)
            dat_tmp = as.data.frame(res[[1]])
            chrom_tmp = as.numeric(res[[2]])
            # check input type and extract snp information and gwas info accordingly
            # example case
            if (target_type == "example"){
              snp_locus = data.frame(chr=16, pos=dat_tmp[ceiling(nrow(dat_tmp)/4), "pos"])
              res <- function.InputPos(dat = dat_tmp, window = input$x, snp_locus, tmp)

              # rsid case
            } else if (target_type == "rsid"){
              snp_locus <- function.rsIDasInput(target)
              res <- function.InputPos(dat = dat_tmp, window = input$x, snp_locus, tmp)

              # chromosome:position case
            } else if (target_type == "locus") {
              snp_locus <- data.frame(chr=target[[2]], pos=target[[3]])
              res <- function.InputPos(dat = dat_tmp, window = input$x, snp_locus, tmp)

              # gene_name case
            } else if (target_type == "gene_name"){
              gene.info <- function.InputGenes(gene = target[[2]])
              res <- function.GWASfromGene(dat = dat_tmp, gene.info = gene.info, window = input$x, tmp)
            }
            res$Study = tmp
            list_for_results <- rbind(list_for_results, res)
            if (gw == 2){ col_list[[gw]] <- input$col2 } else if (gw == 3){ col_list[[gw]] <- input$col3 } else if (gw == 4){ col_list[[gw]] <- input$col4 } else if (gw == 5){ col_list[[gw]] <- input$col5 }
          }
        }
        toplot <- list_for_results
        # This section is the one with respect to the loaded file
      } else {
        # loop over the gwas to plot
        # define some outputs -- the first is the names of the gwases
        list_for_loop <- all_gwas
        # the second is the actual data to plot
        list_for_results <- list()
        # the third is the color for each gwas
        col_list <- list()
        for (gw in 1:length(list_for_loop)){
          if (gw == 1){
            # check input type and extract snp information and gwas info accordingly
            # example case
            if (target_type == "example"){
              snp_locus = data.frame(chr=16, pos=dat[ceiling(nrow(dat)/4), "pos"])
              res <- function.InputPos(dat = dat, window = input$x, snp_locus, gwas)

              # rsid case
            } else if (target_type == "rsid"){
              snp_locus <- function.rsIDasInput(target)
              res <- function.InputPos(dat = dat, window = input$x, snp_locus, gwas)

              # chromosome:position case
            } else if (target_type == "locus") {
              snp_locus <- data.frame(chr=target[[2]], pos=target[[3]])
              res <- function.InputPos(dat = dat, window = input$x, snp_locus, gwas)

              # gene_name case
            } else if (target_type == "gene_name"){
              gene.info <- function.InputGenes(gene = target[[2]])
              res <- function.GWASfromGene(dat = dat, gene.info = gene.info, window = input$x, gwas)
            }
            list_for_results[[gw]] <- res
            col_list[[gw]] <- input$col
          } else {
            # other datasets
            # check if the second dataset is the loaded file
            if (list_for_loop[[gw]] != "toLoad"){
              tmp <- as.character(all_gwas[gw])
              res = function.manageInput(tmp, path_f, res_example)
              dat_tmp = as.data.frame(res[[1]])
              chrom_tmp = as.numeric(res[[2]])
              # check input type and extract snp information and gwas info accordingly
              # example case
              if (target_type == "example"){
                snp_locus = data.frame(chr=16, pos=dat_tmp[ceiling(nrow(dat_tmp)/4), "pos"])
                res <- function.InputPos(dat = dat_tmp, window = input$x, snp_locus, tmp)

                # rsid case
              } else if (target_type == "rsid"){
                snp_locus <- function.rsIDasInput(target)
                res <- function.InputPos(dat = dat_tmp, window = input$x, snp_locus, tmp)

                # chromosome:position case
              } else if (target_type == "locus") {
                snp_locus <- data.frame(chr=target[[2]], pos=target[[3]])
                res <- function.InputPos(dat = dat_tmp, window = input$x, snp_locus, tmp)

                # gene_name case
              } else if (target_type == "gene_name"){
                gene.info <- function.InputGenes(gene = target[[2]])
                res <- function.GWASfromGene(dat = dat_tmp, gene.info = gene.info, window = input$x, tmp)
              }
              list_for_results[[gw]] <- res
              if (gw == 2){ col_list[[gw]] <- input$col2 } else if (gw == 3){ col_list[[gw]] <- input$col3 } else if (gw == 4){ col_list[[gw]] <- input$col4 } else if (gw == 5){ col_list[[gw]] <- input$col5 }

            # otherwise if we need to plot the loaded file
            } else {
              tmp <- all_gwas[gw]
              res = function.manageInput(tmp, path_f, res_example)
              dat_tmp = as.data.frame(res[[1]])
              chrom_tmp = as.numeric(res[[2]])
              if (target_type == "example"){
                snp_locus = data.frame(chr=16, pos=dat_tmp[ceiling(nrow(dat_tmp)/4), "pos"])
                res <- function.InputPos(dat = dat_tmp, window = input$x, snp_locus, tmp)

              # rsid case
              } else if (target_type == "rsid"){
                snp_locus <- function.rsIDasInput(target)
                res <- function.InputPos(dat = dat_tmp, window = input$x, snp_locus, tmp)

              # chromosome:position case
              } else if (target_type == "locus"){
                snp_locus <- data.frame(chr=target[[2]], pos=target[[3]])
                res <- function.InputPos(dat = dat_tmp, window = input$x, snp_locus, tmp)

                # gene_name case
              } else if (target_type == "gene_name"){
                gene.info <- function.InputGenes(gene = target[[2]])
                res <- function.GWASfromGene(dat = dat_tmp, gene.info = gene.info, window = input$x, tmp)
              }
            list_for_results[[gw]] <- res
            if (gw == 2){ col_list[[gw]] <- input$col2 } else if (gw == 3){ col_list[[gw]] <- input$col3 } else if (gw == 4){ col_list[[gw]] <- input$col4 } else if (gw == 5){ col_list[[gw]] <- input$col5 }
            }
          }
        }

        # merge results together
        for (i in 1:length(list_for_results)){
          tmp <- list_for_results[[i]]
          if (nrow(tmp) >0){
            tmp$Study <- list_for_loop[i]
            list_for_results[[i]] <- tmp
          } else {
            list_for_results[[i]] <- NULL
          }
        }
        toplot <- rbindlist(list_for_results)

        #################################
      }
    }

    # This in case you want manual scroll
  } else if (input$sel == "Manual scroll"){
    if (length(all_gwas) < 2){
      # manage different manual scroll input
      snp.info <- function.InputManualScroll(dat = dat, window = input$x, input.scroll = input$manual.pos, input.chrom = input$manual.chrom, gwas)

      # assign to toplot object
      toplot <- snp.info

      # The following in case multiple gwas need to be plotted
    } else {
      # check whether uploaded file is among those to show
      if (path_f == "None"){
        # loop over the gwas to plot
        # define some outputs -- the first is the names of the gwases to be plotted
        list_for_loop <- all_gwas
        # the second i the actual data to plot
        list_for_results <- list()
        # the third is the color for each gwas
        col_list <- list()
        for (gw in 1:length(list_for_loop)){
          if (gw == 1){
            tmp <- as.character(all_gwas[gw])
            dat_tmp <- function.InputManualScroll(dat = dat, window = input$x, input.scroll = input$manual.pos, input.chrom = input$manual.chrom, tmp)
            list_for_results[[gw]] <- dat_tmp
            col_list[[gw]] <- input$col
          } else {
            tmp <- as.character(all_gwas[gw])
            res <- function.manageInput(tmp, path_f, res_example)
            dat_tmp <- as.data.frame(res[[1]])
            chrom_tmp <- as.data.frame(res[[2]])
            dat_tmp <- function.InputManualScroll(dat = dat_tmp, window = input$x, input.scroll = input$manual.pos, input.chrom = input$manual.chrom, tmp)
            list_for_results[[gw]] <- dat_tmp
            if (gw == 2){ col_list[[gw]] <- input$col2 } else if (gw == 3){ col_list[[gw]] <- input$col3 } else if (gw == 4){ col_list[[gw]] <- input$col4 } else if (gw == 5){ col_list[[gw]] <- input$col5 }
          }
        }

        # merge results together
        for (i in 1:length(list_for_results)){
          tmp <- list_for_results[[i]]
          if (nrow(tmp) >0){
            tmp$Study <- all_gwas[i]
            list_for_results[[i]] <- tmp
          } else {
            list_for_results[[i]] <- NULL
          }
        }
        toplot <- rbindlist(list_for_results)
        # In the following case, the file is loaded -- to be completed
      } else {
        list_for_loop <- all_gwas
        # the second is the actual data to plot
        list_for_results <- list()
        # the third is the color for each gwas
        col_list <- list()
        for (gw in 1:length(list_for_loop)){
          if (gw == 1){
            tmp <- as.character(all_gwas[gw])
            dat_tmp <- function.InputManualScroll(dat = dat, window = input$x, input.scroll = input$manual.pos, input.chrom = input$manual.chrom, tmp)
            list_for_results[[gw]] <- dat_tmp
            col_list[[gw]] <- input$col
          } else {
            tmp <- as.character(all_gwas[gw])
            res <- function.manageInput(tmp, path_f, res_example)
            dat_tmp <- as.data.frame(res[[1]])
            chrom_tmp <- as.data.frame(res[[2]])
            dat_tmp <- function.InputManualScroll(dat = dat_tmp, window = input$x, input.scroll = input$manual.pos, input.chrom = input$manual.chrom, tmp)
            list_for_results[[gw]] <- dat_tmp
            if (gw == 2){ col_list[[gw]] <- input$col2 } else if (gw == 3){ col_list[[gw]] <- input$col3 } else if (gw == 4){ col_list[[gw]] <- input$col4 } else if (gw == 5){ col_list[[gw]] <- input$col5 }
          }
        }

        # merge results together
        for (i in 1:length(list_for_results)){
          tmp <- list_for_results[[i]]
          if (nrow(tmp) >0){
            tmp$Study <- list_for_loop[i]
            list_for_results[[i]] <- tmp
          } else {
            list_for_results[[i]] <- NULL
          }
        }
        toplot <- rbindlist(list_for_results)

        #######################################
      }
    }
  }

  # for the table, add the study name as a column
  if (length(all_gwas) <2){
    if (nrow(toplot) >0){
      toplot$Study <- gwas
    }
  }
  # finally reorder
  toplot <- toplot[order(toplot$p),]
  toplot$p <- NULL
  #print(head(toplot))
  return(toplot)
}

## Function to output the table of SVs in the region
plotTable_SVs <- function(input){
  ############################
  #create flag for uploaded input file
  inFile <- input$inp_f
  path_f <- "None"
  if (!is.null(inFile)){
    path_f = inFile$datapath
    #print("file loaded")
  } else {
    path_f <- "None"
  }
  ##########################

  #############################
  all_gwas <- c(input$gwas_neuro, input$gwas_cardio, input$gwas_immune, input$gwas_cancer, input$gwas_physio, input$gwas_from_file)
  res = function.manageInput(as.character(all_gwas[1]), path_f, res_example)
  dat = as.data.frame(res[[1]])
  chrom = as.numeric(res[[2]])
  gwas = as.character(res[[3]])
  #############################
  if (input$sel == "Locus, Gene, RsID"){
    # If input is a target region, need to identify which was the input (locus, gene or rsid)
    target <- identiTargetRegion(input$target)
    target_type <- target[[1]]

    # check input type and extract snp information and gwas info accordingly
    # example case
    if (target_type == "example"){
      snp_locus = data.frame(chr=16, pos=dat[ceiling(nrow(dat)/4), "pos"])
    # rsid case
    } else if (target_type == "rsid"){
      snp_locus <- function.rsIDasInput(target)
    # chromosome:position case
    } else if (target_type == "locus") {
      snp_locus <- data.frame(chr=target[[2]], pos=target[[3]])
    # gene_name case
    } else if (target_type == "gene_name"){
      gene.info <- function.InputGenes(gene = target[[2]])
    }

    # assign to toplot object
    if (target_type == "gene_name"){
      tmp <- gene.info[, c("chrom", "txStart", "txEnd")]
      toplot <- data.frame(chrom = c(tmp$chrom, tmp$chrom), pos = c(as.numeric(tmp$txStart) - as.numeric(input$x), as.numeric(tmp$txEnd) + as.numeric(input$x)))
      colnames(toplot) <- c("chr", "pos")
    } else {
      tmp <- snp_locus
      toplot <- data.frame(chrom = c(tmp$chr, tmp$chr), pos = c(as.numeric(tmp$pos) - as.numeric(input$x), as.numeric(tmp$pos) + as.numeric(input$x)))
      colnames(toplot) <- c("chr", "pos")
    }

    # This in case you want manual scroll
  } else if (input$sel == "Manual scroll"){
    tmp <- function.InputManualScroll(dat = dat, window = input$x, input.scroll = input$manual.pos, input.chrom = input$manual.chrom, gwas[1])
    toplot <- data.frame(chr = c(unique(tmp$chr), unique(tmp$chr)), pos = c(min(as.numeric(tmp$pos)), max(as.numeric(tmp$pos))))
  }

  # remove the chr in front of the chromosome number
  if (length(grep("chr", toplot$chr)) >0){
    toplot$chr <- str_split_fixed(toplot$chr, "chr", 2)[, 2]
  }
  # once we have where to look, let's find SVs
  tmp_svs <- findStr_variats(snp.info = toplot, genV = input$genV, inpStrVar = input$strVar_inp)
  tmp_svs$col <- NULL
  tmp_svs$y_plot <- NULL
  tmp_svs$middle <- NULL
  tmp_svs$source[which(tmp_svs$source == "jasper")] <- "Linthorst_et_al_2020"
  tmp_svs$source[which(tmp_svs$source == "audano")] <- "Audano_et_al_2019"
  tmp_svs$source[which(tmp_svs$source == "chaisson")] <- "Chaisson_et_al_2019"
  return(tmp_svs)
}

## Function to output the table of LD patterns
plotTable_LD <- function(input){
  ############################
  #create flag for uploaded input file
  inFile <- input$inp_f
  path_f <- "None"
  if (!is.null(inFile)){
    path_f = inFile$datapath
    #print("file loaded")
  } else {
    path_f <- "None"
  }
  ##########################

  #############################
  all_gwas <- c(input$gwas_neuro, input$gwas_cardio, input$gwas_immune, input$gwas_cancer, input$gwas_physio, input$gwas_from_file)
  res = function.manageInput(as.character(all_gwas[1]), path_f, res_example)
  dat = as.data.frame(res[[1]])
  chrom = as.numeric(res[[2]])
  gwas = as.character(res[[3]])
  #############################
  if (input$sel == "Locus, Gene, RsID"){
    # If input is a target region, need to identify which was the input (locus, gene or rsid)
    target <- identiTargetRegion(input$target)
    target_type <- target[[1]]

    # check number of GWAS to be plotted
    all_gwas <- c(input$gwas_neuro, input$gwas_cardio, input$gwas_immune, input$gwas_cancer, input$gwas_physio, input$gwas_from_file)
    if (length(all_gwas) < 2){
      # check input type and extract snp information and gwas info accordingly
      # example case
      if (target_type == "example"){
        snp_locus = data.frame(chr=16, pos=dat[ceiling(nrow(dat)/4), "pos"])
        res <- function.InputPos(dat = dat, window = input$x, snp_locus, gwas)

        # rsid case
      } else if (target_type == "rsid"){
        snp_locus <- function.rsIDasInput(target)
        res <- function.InputPos(dat = dat, window = input$x, snp_locus, gwas)

        # chromosome:position case
      } else if (target_type == "locus") {
        snp_locus <- data.frame(chr=target[[2]], pos=target[[3]])
        res <- function.InputPos(dat = dat, window = input$x, snp_locus, gwas)

        # gene_name case
      } else if (target_type == "gene_name"){
        gene.info <- function.InputGenes(gene = target[[2]])
        res <- function.GWASfromGene(dat = dat, gene.info = gene.info, window = input$x, gwas)
      }

      # assign to toplot object
      toplot <- res
    }
    # This in case you want manual scroll
  } else if (input$sel == "Manual scroll"){
    if (length(input$gwas) < 2){
      # manage different manual scroll input
      snp.info <- function.InputManualScroll(dat = dat, window = input$x, input.scroll = input$manual.pos, input.chrom = input$manual.chrom, gwas)

      # assign to toplot object
      toplot <- snp.info
    }
  }

  # reorder
  snp.info <- toplot[order(toplot$p),]

  # the look into LD
  # first check and store the populations selected in case LD calculation is requested
  pop_interest_ld = c(input$pop_afr, input$pop_amr, input$pop_eas, input$pop_eur, input$pop_sas)

  ## Check whether LD needs to be computed and in case calculate it
  ld_yesNo <- "no"
  if (input$Linkage == "Most significant in window"){
    ld.info <- functionLD(snp.info, snp=NULL, pop_interest = pop_interest_ld)
    #merge with association data
    ld.snp <- merge(ld.info, snp.info, by.x="BP_B", by.y="pos")
    ld_yesNo <- "yes"
  } else if (input$Linkage == "Input variant"){
    ld.info <- functionLD(snp.info, snp=input$pos, pop_interest = pop_interest_ld)
    #merge with association data
    ld.snp <- merge(ld.info, snp.info, by.x="BP_B", by.y="pos")
    ld_yesNo <- "yes"
  }

  ## Also, in case the genome is hg38, now need to do the liftover on the snp.info
  # if genome version is hg38, need to liftover
  if (input$genV == "GRCh38 (hg38)"){
    df <- data.frame(chr=paste("chr", snp.info$chr, sep=""), start=snp.info$pos, end=snp.info$pos+1)
    # change to GR class object
    gr <- makeGRangesFromDataFrame(df)
    # set chain file
    chain <- import.chain( "../data/databases/hg19ToHg38.over.chain")
    # change coordinates
    gr_hg38 <- liftOver(gr, chain)
    # back to dataframe and clean it
    df_hg38 <- as.data.frame(gr_hg38)
    sb <- data.frame(chr=rep(snp.info$chr[1], nrow(df_hg38)), pos=df_hg38$start, p=snp.info[df_hg38$group, "p"], log10p=snp.info[df_hg38$group, "-log10(P-value)"], size=snp.info[df_hg38$group, "size"])
    colnames(sb) <- c("chr", "pos", "p", "-log10(P-value)", "size")
    snp.info <- sb
    if (ld_yesNo == "yes"){
      df <- data.frame(chr=paste("chr", ld.snp$CHR_A[1], sep=""), start=ld.snp$BP_B, end=ld.snp$BP_B+1)
      gr <- makeGRangesFromDataFrame(df)
      chain <- import.chain( "../data/databases/hg19ToHg38.over.chain")
      gr_hg38 <- liftOver(gr, chain)
      df_hg38 <- as.data.frame(gr_hg38)
      sb <- data.frame(BP_B=df_hg38$start, col=ld.snp[df_hg38$group, "col"], size=ld.snp[df_hg38$group, "size"], p=ld.snp[df_hg38$group, "p"])
      sb$"-log10(P-value)" <- -log10(sb$p)
      ld.snp <- sb
    }
  }

  # Clean and we're done
  ld.snp$pch <- NULL
  ld.snp$col <- NULL
  return(ld.snp)
}

## Function to output the table of the top eQTL associations
plotEqtl_table_blood <- function(snps_in_interval, eQTL, mapping.ensembl, GENOME){
  # first exclude duplicates based on position
  snps_in_interval <- snps_in_interval[!duplicated(snps_in_interval$Position),]
  # do not check all snps -- only top 50 let's say
  snps_in_interval = head(snps_in_interval, 100)
  # then take eqtl in interval
  if (GENOME == "GRCh37 (hg19)"){
    eqtl_in_interval <- eQTL[which((eQTL$chr %in% paste0("chr", snps_in_interval$Chr)) & (eQTL$pos_hg37 %in% snps_in_interval$Position)),]
    eqtl_in_interval$Locus_hg37 <- paste(eqtl_in_interval$chr, eqtl_in_interval$pos_hg37, sep=":")
    tmp = snps_in_interval[, c("Position", "ID")]
    eqtl_in_interval = merge(eqtl_in_interval, tmp, by.x="pos_hg37", by.y = "Position")
  } else {
    eqtl_in_interval <- eQTL[which((eQTL$chr %in% paste0("chr", snps_in_interval$Chr)) & (eQTL$pos_hg38 %in% snps_in_interval$Position)),]
    eqtl_in_interval$Locus_hg38 <- paste(eqtl_in_interval$chr, eqtl_in_interval$pos_hg38, sep=":")
    tmp = snps_in_interval[, c("Position", "ID")]
    eqtl_in_interval = merge(eqtl_in_interval, tmp, by.x="pos_hg38", by.y = "Position")
  }
  # if there are entries, nedd to assign gene name
  if (nrow(eqtl_in_interval) >0){
    # take gene without the .XX and the merge with eqtl file
    eqtl_in_interval$ensemble <- str_split_fixed(eqtl_in_interval$gene_id, "\\.", 2)[, 1]
    eqtl_in_interval <- merge(eqtl_in_interval, mapping.ensembl, by.x="ensemble", by.y="Gene stable ID")
    eqtl_in_interval$log_adjP <- -log10(as.numeric(eqtl_in_interval$pval_beta))
    # clean and we are done
    if (GENOME == "GRCh38 (hg38)"){
      eqtl_in_interval <- eqtl_in_interval[, c("Locus_hg38", "ID", "a1", "a2", "slope", "log_adjP", "Gene name")]
    } else {
      eqtl_in_interval <- eqtl_in_interval[, c("Locus_hg37", "ID", "a1", "a2", "slope", "log_adjP", "Gene name")]
    }
    # order by p and return top 15 entries
    colnames(eqtl_in_interval) <- c("Locus", "ID", "A1", "A2", "Effect", "P", "Gene")
    eqtl_in_interval$"A1/A2" <- paste(eqtl_in_interval$A1, eqtl_in_interval$A2, sep="/")
    eqtl_in_interval$A1 <- NULL
    eqtl_in_interval$A2 <- NULL
    eqtl_in_interval <- eqtl_in_interval[order(-eqtl_in_interval$"P"),]

  } else {
    eqtl_in_interval <- data.frame(Locus="No", "A1/A2"="eQTLs", Effect="in displayed", "P"="region", Gene="-")
  }
  return(eqtl_in_interval)
}

## Function to find rsids from a CADD file we have
assignRSID <- function(data, genV){
  if (genV != "GRCh38 (hg38)"){
    kg = fread(paste0("../AnnotateMe/INPUTS_OTHER/1000G_frequencies/chr", data$Chr[1], ".afreq.gz"), h=T, stringsAsFactors = F)
    data = merge(data, kg, by.x="Position", by.y = "POS", all.x = T)
  } else {
    kg = fread(paste0("../AnnotateMe/INPUTS_OTHER/1000G_frequencies/chr", data$Chr[1], ".afreq.gz"), h=T, stringsAsFactors = F)
    data = merge(data, kg, by.x="Position", by.y = "POS", all.x = T)
    data$group = seq(1, nrow(data))
    df <- data.frame(chr=paste("chr", data$Chr, sep=""), start=data$Position, end=data$Position+1)
    # change to GR class object
    gr <- makeGRangesFromDataFrame(df)
    # set chain file
    chain <- import.chain( "../data/databases/hg19ToHg38.over.chain")
    # change coordinates
    gr_hg38 <- liftOver(gr, chain)
    # back to dataframe and clean it
    df_hg38 <- as.data.frame(gr_hg38)
    # merge with rsid
    df_hg38 = merge(df_hg38, data, by = "group")
    data = df_hg38[, c("Chr", "start", "ID", "REF", "ALT", "-log10(P)", "Study")]
    colnames(data) = c("Chr", "Position", "ID", "REF", "ALT", "-log10(P)", "Study")
  }
  data = data[order(-data$"-log10(P)"),]
  data = data[, c("Chr", "Position", "ID", "REF", "ALT", "-log10(P)", "Study")]
  data$Alleles = paste(data$REF, data$ALT, sep="/")
  data$REF = NULL
  data$ALT = NULL
  return(data)
}

## Function to output the table of the top eQTL associations -- this should look into other tissues than blood
plotEqtl_table_all <- function(snps_in_interval, eQTL, mapping.ensembl, GENOME, tissues){
  # first exclude duplicates based on position
  snps_in_interval <- snps_in_interval[!duplicated(snps_in_interval$Position),]
  # consider only top 100 to make it faster
  snps_in_interval = head(snps_in_interval, 100)
  # then take eqtl in interval
  if (GENOME == "GRCh37 (hg19)"){
    df <- data.frame(chr=paste("chr", snps_in_interval$Chr, sep=""), start=snps_in_interval$Position, end=snps_in_interval$Position+1)
    df$group = seq(1, nrow(df))
    # change to GR class object
    gr <- makeGRangesFromDataFrame(df)
    # set chain file
    chain <- import.chain( "../data/databases/hg19ToHg38.over.chain")
    # change coordinates
    gr_hg38 <- liftOver(gr, chain)
    # back to dataframe and clean it
    df_hg38 <- as.data.frame(gr_hg38)
    df_hg38 = merge(df_hg38, df, by = "group")
    df_hg38 = merge(df_hg38, snps_in_interval, by.x = "start.y", by.y = "Position")
    sb = df_hg38[, c("chr", "start.x", "start.y", "ID", "Alleles", "Study", "-log10(P)")]
    colnames(sb) <- c("chr", "pos_hg38", "pos_hg19", "ID", "Alleles", "Study", "-log10(P)")
  } else {
    snps_in_interval$pos_hg19 = NA
    sb = snps_in_interval[, c("Chr", "Position", "pos_hg19", "ID", "Alleles", "Study", "-log10(P)")]
    colnames(sb) <- c("chr", "pos_hg38", "pos_hg19", "ID", "Alleles", "Study", "-log10(P)")
    sb$chr = paste0("chr", sb$chr)
  }

  # for the grep, need to create some ids
  sb$IDs = paste0(sb$chr, "_", sb$pos_hg38, "_")
  # read gtex file for the correct chromosome
  tmp_gtx = fread(paste0("/root/snpXplorer/AnnotateMe/INPUTS_OTHER/summary_eqtls/", unique(sb$chr), "_summary_eqtls.txt.gz"), h=F, sep="\t", stringsAsFactors = F)
  tmp_gtx$pos = as.numeric(str_split_fixed(tmp_gtx$V1, "_", 5)[, 2])
  tmp_gtx = merge(tmp_gtx, sb, by.x = "pos", by.y="pos_hg38")
  # then grep on the tissue
  if ("All_tissues" %in% tissues){
    pattern_tis = paste0(c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary",
                         "Artery_Tibial", "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia",
                         "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus",
                         "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1",
                         "Brain_Substantia_nigra", "Breast_Mammary_Tissue", "Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes",
                         "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis",
                         "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Kidney_Cortex", "Liver", "Lung", "Minor_Salivary_Gland", "Muscle_Skeletal",
                         "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg",
                         "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood"), collapse = "|")
  } else {
    pattern_tis = paste0(tissues, collapse = "|")
  }
  matches_tis = tmp_gtx[grep(pattern = pattern_tis, x = tmp_gtx$V2),]
  matches_tis$V1 = as.character(matches_tis$V1)
  matches_tis$V2 = as.character(matches_tis$V2)
  # temporary function that does the grep
  tmp_f = function(i, matches_tis, pattern_tis, mapping.ensembl){
    tmp_df = as.data.frame(matrix(data=NA, nrow=0, ncol=7))
    colnames(tmp_df) = c("Locus", "ID", "Effect", "P", "Gene", "A1/A2", "Tissue")
    tmp = unlist(strsplit(x = matches_tis$V2[i], split = ";"))[grep(pattern = pattern_tis, x = unlist(strsplit(x = matches_tis$V2[i], split = ";")))]
    for (x in tmp){
      gene_name = as.character(mapping.ensembl[grep(str_split_fixed(unlist(strsplit(x, "_"))[length(unlist(strsplit(x, "_")))-2], "\\.", 2)[, 1], mapping.ensembl$`Gene stable ID`), "Gene name"])
      if (length(gene_name) > 0){
        loc = as.character(paste0(str_split_fixed(as.character(matches_tis$V1[i]), "_", 5)[1:2], collapse = ":"))
        es = as.character(unlist(strsplit(x, "_"))[length(unlist(strsplit(x, "_")))-1])
        pval = as.numeric(unlist(strsplit(x, "_"))[length(unlist(strsplit(x, "_")))])
        alleles = paste0(str_split_fixed(matches_tis$V1[i], "_", 5)[3:4], collapse = "/")
        tiss = paste0(unlist(strsplit(x, "_"))[1:(grep("ENS", unlist(strsplit(x, "_")))-1)], collapse = "_")
        rsid = as.character(matches_tis$ID[i])
        tmp_df2 = data.frame(Locus = loc, ID = rsid, Effect = es, P = pval, Gene = gene_name, "A1/A2" = alleles, Tissue = tiss)
        tmp_df = rbind(tmp_df, tmp_df2)
      }
    }
    return(tmp_df)
  }

  if (nrow(matches_tis) >0){
    all_eqtls = lapply(X = 1:nrow(matches_tis), FUN = tmp_f, matches_tis = matches_tis, pattern_tis = pattern_tis, mapping.ensembl = mapping.ensembl)
    eqtl_in_interval <- rbindlist(all_eqtls)
    eqtl_in_interval$P = -log10(as.numeric(eqtl_in_interval$P))
    eqtl_in_interval = eqtl_in_interval[order(-eqtl_in_interval$P),]
  } else {
    eqtl_in_interval <- data.frame(Locus="No", "A1/A2"="eQTLs", Effect="in displayed", "P"="region", Gene="-")
  }

  return(eqtl_in_interval)
}

##########
##########
# MAIN APP
# Let's load the annotation data for speeding-up the plots
# save.image("../data/databases/annotationFiles.RData")
load("../data/databases/annotationFiles.RData")

shinyServer(
  function(input, output, session) {
    # Main function to manage all data and conditions -- the output is the plot object
    plt <- function(input, output, session){

      ############################
      # Create flag for uploaded input file
      inFile <- input$inp_f
      if (!is.null(inFile)){ path_f = inFile$datapath } else { path_f <- "None" }
      ##########################

      ##########################
      # Define which input files to plot based on user choice
      if (length(as.character(input$gwas_from_file[1])) == 0){
        all_gwases = c(input$gwas_neuro, input$gwas_cardio, input$gwas_immune, input$gwas_cancer, input$gwas_physio)
        res = function.manageInput(as.character(all_gwases[1]), path_f, res_example)
      } else {
        res = function.manageInput(as.character(input$gwas_from_file[1]), path_f, res_example)
      }
      dat = as.data.frame(res[[1]])
      chrom = as.numeric(res[[2]])
      gwas = as.character(res[[3]])
      n.inputs = 1
      ####################

      # Create a variable to control the error and the plot
      plot_error = FALSE
      plot_error_sex = FALSE

      # Store the populations selected in case LD calculation is requested
      pop_interest_ld = c(input$pop_afr, input$pop_amr, input$pop_eas, input$pop_eur, input$pop_sas)

      # Then we need to look at the browsing options: now you can choose between "locus, gene, rsid" and "manual scroll"
      if (input$sel == "Locus, Gene, RsID"){
        # If input is a target region, need to identify which was the input (locus, gene or rsid)
        target <- identiTargetRegion(input$target)
        target_type <- target[[1]]

        # check number of GWAS to be plotted
        all_gwas <- c(input$gwas_neuro, input$gwas_cardio, input$gwas_immune, input$gwas_cancer, input$gwas_physio, input$gwas_from_file)
        if (length(all_gwas) < 2){
          # check input type and extract snp information and gwas info accordingly
          # example case
          if (target_type == "example"){
            snp_locus = data.frame(chr=16, pos=dat[ceiling(nrow(dat)/4), "pos"])
            res <- function.InputPos(dat = dat, window = input$x, snp_locus, gwas)
          # rsid case
          } else if (target_type == "rsid"){
            snp_locus <- function.rsIDasInput(target)
            print(snp_locus)
            if (snp_locus != "snp_not_in_list"){
              if (snp_locus$chr %in% c("X", "x", "Y", "y")){
                plot_error = TRUE
                plot_error_sex = TRUE
              } else {
                res <- function.InputPos(dat = dat, window = input$x, snp_locus, gwas)
                if ((nrow(res) == 0) || (snp_locus$chr != res$chr)){
                  plot_error = TRUE
                  plot_error_sex = FALSE
                }
              }
            } else {
              plot_error = TRUE
              plot_error_sex = FALSE
            }

          # chromosome:position case
          } else if (target_type == "locus") {
            snp_locus <- data.frame(chr=target[[2]], pos=target[[3]])
            if (is.na(snp_locus$chr)){
              plot_error = TRUE
              plot_error_sex = TRUE
            } else {
              res <- function.InputPos(dat = dat, window = input$x, snp_locus, gwas)
              if ((nrow(res) == 0) || (snp_locus$chr != res$chr)){
                plot_error = TRUE
                plot_error_sex = FALSE
              }
            }

          # gene_name case
          } else if (target_type == "gene_name"){
            gene.info <- function.InputGenes(gene = target[[2]])
            if (gene.info != "gene_not_in_list"){
              if (gene.info$chrom %in% c("chrx", "chrX", "chrY", "chry")){
                plot_error = TRUE
                plot_error_sex = TRUE
              } else {
                res <- function.GWASfromGene(dat = dat, gene.info = gene.info, window = input$x, gwas)
                if ((nrow(res) == 0) || (paste0("chr", res$chr) != gene.info$chrom)){
                  plot_error = TRUE
                  plot_error_sex = FALSE
                }
              }
            } else {
              plot_error = TRUE
              plot_error_sex = FALSE
            }
          }
          # plot
          if (plot_error == TRUE){
            par(mar=c(4, 8, 10, 6))
            plot(0, 0, xlim=c(0, 55), ylim=c(0, 25), bty="n", xaxt='none', yaxt='none', xlab="", ylab="",
                 main=paste0(input$target), cex.main=2.50, col="white")
            if (is.null(all_gwas)){
              text(x = 27.5, y = 25, labels = "Ops, looks like you forgot to check at least an input dataset.",
                   cex = 2, xpd=T)
            } else if (plot_error_sex == TRUE){
              text(x = 27.5, y = 25, labels = "Ops, looks like you are interested in sexual chromosomes.\nSorry, as of today, we do not store sexual chromosome data.",
                   cex = 2, xpd=T)
          } else {
              text(x = 27.5, y = 25, labels = "No Gene/SNP found in RefSeq/1000Genomes.\nIf you are visualizing your own data, make sure this region is covered.",
                 cex = 2, xpd=T) }
          } else {
            function.plot(snp.info = res, y.lim = input$y, type = input$sel, plt.type = input$ploType,
                        windows.number = input$sliding.window, smooth.par = input$smooth, int.locus = target,
                        gwas = gwas, col = input$col, ld = input$Linkage, input = input, genV = input$genV, inpStrVar = input$strVar_inp, pop_interest_ld = pop_interest_ld, recomb_yn = input$recomb_yn, dotSize_yn = input$dotSize_yn)
          }

          # The following in case multiple gwas need to be plotted
        } else {
          #check whether uploaded file is among those to show
          if (path_f == "None"){
            # loop over the gwas to plot
            # define some outputs -- the first is the names of the gwases
            list_for_loop <- all_gwas
            # the second is the actual data to plot
            list_for_results <- list()
            # the third is the color for each gwas
            col_list <- list()
            for (gw in 1:length(list_for_loop)){
              if (gw == 1){
                # check input type and extract snp information and gwas info accordingly
                # example case
                if (target_type == "example"){
                  snp_locus = data.frame(chr=16, pos=dat[ceiling(nrow(dat)/4), "pos"])
                  res <- function.InputPos(dat = dat, window = input$x, snp_locus, gwas)

                  # rsid case
                } else if (target_type == "rsid"){
                  snp_locus <- function.rsIDasInput(target)
                  if (snp_locus != "snp_not_in_list"){
                    if (snp_locus$chr %in% c("X", "x", "Y", "y")){
                      plot_error = TRUE
                      plot_error_sex = TRUE
                    } else {
                      res <- function.InputPos(dat = dat, window = input$x, snp_locus, gwas)
                      list_for_results[[gw]] <- res
                      col_list[[gw]] <- input$col
                    }
                  } else {
                    plot_error = TRUE
                  }

                  # chromosome:position case
                } else if (target_type == "locus") {
                  snp_locus <- data.frame(chr=target[[2]], pos=target[[3]])
                  if (is.na(snp_locus$chr)){
                    plot_error = TRUE
                    plot_error_sex = TRUE
                  } else {
                    res <- function.InputPos(dat = dat, window = input$x, snp_locus, gwas)
                    list_for_results[[gw]] <- res
                    col_list[[gw]] <- input$col
                  }

                  # gene_name case
                } else if (target_type == "gene_name"){
                  gene.info <- function.InputGenes(gene = target[[2]])
                  if (gene.info != "gene_not_in_list"){
                    if (gene.info$chrom %in% c("chrx", "chrX", "chry", "chrY")){
                      plot_error = TRUE
                      plot_error_sex = TRUE
                    } else {
                      res <- function.GWASfromGene(dat = dat, gene.info = gene.info, window = input$x, gwas)
                      list_for_results[[gw]] <- res
                      col_list[[gw]] <- input$col
                    }
                  } else {
                    plot_error = TRUE
                  }
                }
              } else if (plot_error == FALSE) {
                # other datasets
                tmp <- as.character(all_gwas[gw])
                res = function.manageInput(tmp, path_f, res_example)
                dat_tmp = as.data.frame(res[[1]])
                chrom_tmp = as.numeric(res[[2]])
                # check input type and extract snp information and gwas info accordingly
                # example case
                if (target_type == "example"){
                  snp_locus = data.frame(chr=16, pos=dat_tmp[ceiling(nrow(dat_tmp)/4), "pos"])
                  res <- function.InputPos(dat = dat_tmp, window = input$x, snp_locus, tmp)

                  # rsid case
                } else if (target_type == "rsid"){
                  snp_locus <- function.rsIDasInput(target)
                  res <- function.InputPos(dat = dat_tmp, window = input$x, snp_locus, tmp)

                  # chromosome:position case
                } else if (target_type == "locus") {
                  snp_locus <- data.frame(chr=target[[2]], pos=target[[3]])
                  res <- function.InputPos(dat = dat_tmp, window = input$x, snp_locus, tmp)

                  # gene_name case
                } else if (target_type == "gene_name"){
                  gene.info <- function.InputGenes(gene = target[[2]])
                  res <- function.GWASfromGene(dat = dat_tmp, gene.info = gene.info, window = input$x, tmp)
                }
                list_for_results[[gw]] <- res
                if (gw == 2){ col_list[[gw]] <- input$col2 } else if (gw == 3){ col_list[[gw]] <- input$col3 } else if (gw == 4){ col_list[[gw]] <- input$col4 } else if (gw == 5){ col_list[[gw]] <- input$col5 }
              }
            }
            # plot
            if (plot_error == TRUE){
              par(mar=c(4, 8, 10, 6))
              plot(0, 0, xlim=c(0, 55), ylim=c(0, 25), bty="n", xaxt='none', yaxt='none', xlab="", ylab="",
                   main=paste0(input$target), cex.main=2.50, col="white")
              if (plot_error_sex == TRUE){
                text(x = 27.5, y = 25, labels = "Ops, looks like you are interested in sexual chromosomes.\nSorry, as of today, we do not store sexual chromosome data.",
                     cex = 2, xpd=T)
              } else {
                text(x = 27.5, y = 25, labels = "No Gene/SNP found in RefSeq/1000Genomes.",
                     cex = 2, xpd=T)
              }
            } else {
              function.multiPlot(snp.info = list_for_results, list_for_loop, y.lim = input$y, type = input$sel,
                                 plt.type = input$ploType, windows.number = input$sliding.window, smooth.par = input$smooth,
                                 int.locus = target, col_list, input$genV, inpStrVar = input$strVar_inp, recomb_yn = input$recomb_yn, dotSize_yn = input$dotSize_yn)
            }

            # This section is the one with respect to the loaded file
          } else {
            # loop over the gwas to plot
            # define some outputs -- the first is the names of the gwases
            list_for_loop <- all_gwas
            list_for_results <- list()
            col_list <- list()
            for (gw in 1:length(list_for_loop)){
              if (gw == 1){
                # check input type and extract snp information and gwas info accordingly
                # example case -- in this case the loaded file is the one leading, then start with that
                if (target_type == "example"){
                  snp_locus = data.frame(chr=16, pos=dat[ceiling(nrow(dat)/4), "pos"])
                  # in this case, do the inputPos based on the loaded file
                  loaded_gwas <- which(all_gwas == "toLoad")
                  res <- function.InputPos(dat = dat, window = input$x, snp_locus, list_for_loop[loaded_gwas])
                  chr_tmp <- unique(res$chr[1])
                  pos_tmp <- res[which(res$chr == chr_tmp),]
                  snp_locus = data.frame(chr=chr_tmp, pos=pos_tmp[ceiling(nrow(pos_tmp)/4), "pos"])
                  res <- function.InputPos(dat = dat, window = input$x, snp_locus, list_for_loop[gw])
                  list_for_results[[gw]] <- res
                  col_list[[gw]] <- input$col

                  # rsid case
                } else if (target_type == "rsid"){
                  snp_locus <- function.rsIDasInput(target)
                  if (snp_locus != "snp_not_in_list"){
                    res <- function.InputPos(dat = dat, window = input$x, snp_locus, list_for_loop[gw])
                    list_for_results[[gw]] <- res
                    col_list[[gw]] <- input$col
                  } else {
                    plot_error = TRUE
                  }

                  # chromosome:position case
                } else if (target_type == "locus") {
                  snp_locus <- data.frame(chr=target[[2]], pos=target[[3]])
                  res <- function.InputPos(dat = dat, window = input$x, snp_locus, list_for_loop[gw])
                  list_for_results[[gw]] <- res
                  col_list[[gw]] <- input$col

                  # gene_name case
                } else if (target_type == "gene_name"){
                  gene.info <- function.InputGenes(gene = target[[2]])
                  if (gene.info != "gene_not_in_list"){
                    res <- function.GWASfromGene(dat = dat, gene.info = gene.info, window = input$x, list_for_loop[gw])
                    list_for_results[[gw]] <- res
                    col_list[[gw]] <- input$col
                  } else {
                    plot_error = TRUE
                  }
                }
              } else if (plot_error == FALSE) {
                # other datasets
                # check if the second dataset is the loaded file
                if (list_for_loop[gw] != "toLoad"){
                  tmp <- as.character(all_gwas[gw])
                  res = function.manageInput(tmp, path_f, res_example)
                  dat_tmp = as.data.frame(res[[1]])
                  chrom_tmp = as.numeric(res[[2]])
                  # check input type and extract snp information and gwas info accordingly
                  # example case
                  if (target_type == "example"){
                    snp_locus = data.frame(chr=16, pos=dat[ceiling(nrow(dat)/4), "pos"])
                    # in this case, do the inputPos based on the loaded file
                    loaded_gwas <- which(all_gwas == "toLoad")
                    res <- function.InputPos(dat = dat, window = input$x, snp_locus, list_for_loop[loaded_gwas])
                    chr_tmp <- unique(res$chr[1])
                    pos_tmp <- res[which(res$chr == chr_tmp),]
                    snp_locus = data.frame(chr=chr_tmp, pos=pos_tmp[ceiling(nrow(pos_tmp)/4), "pos"])
                    res <- function.InputPos(dat = dat, window = input$x, snp_locus, list_for_loop[gw])
                    list_for_results[[gw]] <- res
                    col_list[[gw]] <- input$col

                    # rsid case
                  } else if (target_type == "rsid"){
                    snp_locus <- function.rsIDasInput(target)
                    res <- function.InputPos(dat = dat_tmp, window = input$x, snp_locus, tmp)

                    # chromosome:position case
                  } else if (target_type == "locus") {
                    snp_locus <- data.frame(chr=target[[2]], pos=target[[3]])
                    res <- function.InputPos(dat = dat_tmp, window = input$x, snp_locus, tmp)

                    # gene_name case
                  } else if (target_type == "gene_name"){
                    gene.info <- function.InputGenes(gene = target[[2]])
                    res <- function.GWASfromGene(dat = dat_tmp, gene.info = gene.info, window = input$x, tmp)
                  }
                  list_for_results[[gw]] <- res
                  if (gw == 2){ col_list[[gw]] <- input$col2 } else if (gw == 3){ col_list[[gw]] <- input$col3 } else if (gw == 4){ col_list[[gw]] <- input$col4 } else if (gw == 5){ col_list[[gw]] <- input$col5 }

                  # otherwise if we need to plot the loaded file
                } else {
                  tmp <- all_gwas[gw]
                  res = function.manageInput(tmp, path_f, res_example)
                  dat_tmp = as.data.frame(res[[1]])
                  chrom_tmp = as.numeric(res[[2]])
                  if (target_type == "example"){
                    snp_locus = data.frame(chr=16, pos=dat_tmp[ceiling(nrow(dat_tmp)/4), "pos"])
                    res <- function.InputPos(dat = dat, window = input$x, snp_locus, list_for_loop[loaded_gwas])
                    chr_tmp <- unique(res$chr[1])
                    pos_tmp <- res[which(res$chr == chr_tmp),]
                    snp_locus = data.frame(chr=chr_tmp, pos=pos_tmp[ceiling(nrow(pos_tmp)/4), "pos"])
                    res <- function.InputPos(dat = dat_tmp, window = input$x, snp_locus, tmp)

                    # rsid case
                  } else if (target_type == "rsid"){
                    snp_locus <- function.rsIDasInput(target)
                    res <- function.InputPos(dat = dat_tmp, window = input$x, snp_locus, tmp)
                    if ((nrow(res) == 0) || (snp_locus$chr != res$chr)){
                      list_for_results[[gw]] <- NULL
                    } else {
                      list_for_results[[gw]] <- res
                    }
                    # chromosome:position case
                  } else if (target_type == "locus"){
                    snp_locus <- data.frame(chr=target[[2]], pos=target[[3]])
                    res <- function.InputPos(dat = dat_tmp, window = input$x, snp_locus, tmp)
                    if ((nrow(res) == 0) || (snp_locus$chr != res$chr)){
                      list_for_results[[gw]] <- NULL
                    } else {
                      list_for_results[[gw]] <- res
                    }

                    # gene_name case
                  } else if (target_type == "gene_name"){
                    gene.info <- function.InputGenes(gene = target[[2]])
                    res <- function.GWASfromGene(dat = dat_tmp, gene.info = gene.info, window = input$x, tmp)
                    if ((nrow(res) == 0) || (paste0("chr", res$chr) != gene.info$chrom)){
                      list_for_results[[gw]] <- NULL
                    } else {
                      list_for_results[[gw]] <- res
                    }
                  }
                  if (gw == 2){ col_list[[gw]] <- input$col2 } else if (gw == 3){ col_list[[gw]] <- input$col3 } else if (gw == 4){ col_list[[gw]] <- input$col4 } else if (gw == 5){ col_list[[gw]] <- input$col5 }
                }
              }
            }
            # plot
            if (plot_error == TRUE){
              par(mar=c(4, 8, 10, 6))
              plot(0, 0, xlim=c(0, 55), ylim=c(0, 25), bty="n", xaxt='none', yaxt='none', xlab="", ylab="",
                   main=paste0(input$target), cex.main=2.50, col="white")
              text(x = 27.5, y = 25, labels = "No Gene/SNP found in RefSeq/1000Genomes.",
                   cex = 2, xpd=T)
            } else {
              function.multiPlot(snp.info = list_for_results, list_for_loop, y.lim = input$y, type = input$sel,
                                 plt.type = input$ploType, windows.number = input$sliding.window, smooth.par = input$smooth,
                                 int.locus = target, col_list, input$genV, inpStrVar = input$strVar_inp, recomb_yn = input$recomb_yn, dotSize_yn = input$dotSize_yn)
            }

           #################################
          }
        }

        # This in case you want manual scroll
      } else if (input$sel == "Manual scroll"){
        all_gwas <- c(input$gwas_neuro, input$gwas_cardio, input$gwas_immune, input$gwas_cancer, input$gwas_physio, input$gwas_from_file)
        #print(paste0("Hello, manual scroll selected --> all_gwas:", all_gwas))
        if (length(all_gwas) < 2){
          # manage different manual scroll input
          snp.info <- function.InputManualScroll(dat = dat, window = input$x, input.scroll = input$manual.pos, input.chrom = input$manual.chrom, gwas)
          if ((nrow(snp.info) == 0) || (snp.info$chr != as.character(input$manual.chrom))){
            plot_error = TRUE
          }
          # plot
          if (plot_error == FALSE){
            function.plot(snp.info = snp.info, y.lim = input$y, type = input$sel, plt.type = input$ploType,
                        windows.number = input$sliding.window, smooth.par = input$smooth, int.locus = input$manual.pos,
                        gwas = gwas, col = input$col, ld=input$Linkage, genV = input$genV, inpStrVar = input$strVar_inp, pop_interest_ld = pop_interest_ld, recomb_yn = input$recomb_yn, dotSize_yn = input$dotSize_yn)
          } else {
            par(mar=c(4, 8, 10, 6))
            plot(0, 0, xlim=c(0, 55), ylim=c(0, 25), bty="n", xaxt='none', yaxt='none', xlab="", ylab="",
                 main=paste0(input$manual.chrom, ":", input$manual.pos), cex.main=2.50, col="white")
            text(x = 27.5, y = 25, labels = "No SNPs in this region.",
                 cex = 2, xpd=T)
          }

          # The following in case multiple gwas need to be plotted
        } else {
          # check whether uploaded file is among those to show
          if (path_f == "None"){
            # loop over the gwas to plot
            # define some outputs -- the first is the names of the gwases to be plotted
            list_for_loop <- all_gwas
            # the second i the actual data to plot
            list_for_results <- list()
            # the third is the color for each gwas
            col_list <- list()
            for (gw in 1:length(list_for_loop)){
              if (gw == 1){
                tmp <- as.character(all_gwas[1])
                dat_tmp <- function.InputManualScroll(dat = dat, window = input$x, input.scroll = input$manual.pos, input.chrom = input$manual.chrom, tmp)
                list_for_results[[gw]] <- dat_tmp
                col_list[[gw]] <- input$col
              } else {
                tmp <- as.character(all_gwas[gw])
                res <- function.manageInput(tmp, path_f, res_example)
                dat_tmp <- as.data.frame(res[[1]])
                chrom_tmp <- as.data.frame(res[[2]])
                dat_tmp <- function.InputManualScroll(dat = dat_tmp, window = input$x, input.scroll = input$manual.pos, input.chrom = input$manual.chrom, tmp)
                list_for_results[[gw]] <- dat_tmp
                if (gw == 2){ col_list[[gw]] <- input$col2 } else if (gw == 3){ col_list[[gw]] <- input$col3 } else if (gw == 4){ col_list[[gw]] <- input$col4 } else if (gw == 5){ col_list[[gw]] <- input$col5 }
              }
            }

            # check if there are missing among the gwas to be plotted
            miss = c()
            for (x in list_for_results){
              if (nrow(x) == 0){ miss <- c(miss, 1)}
            }
            if (length(miss) == length(list_for_results)){
              par(mar=c(4, 8, 10, 6))
              plot(0, 0, xlim=c(0, 55), ylim=c(0, 25), bty="n", xaxt='none', yaxt='none', xlab="", ylab="",
                   main=paste0(input$manual.chrom, ":", input$manual.pos), cex.main=2.50, col="white")
              text(x = 27.5, y = 25, labels = "No SNPs in this region.",
                   cex = 2, xpd=T)

            } else {
            # plot
              function.multiPlot(snp.info = list_for_results, list_for_loop, y.lim = input$y, type = input$sel,
                               plt.type = input$ploType, windows.number = input$sliding.window, smooth.par = input$smooth,
                               int.locus = target, col_list, input$genV, inpStrVar = input$strVar_inp, recomb_yn = input$recomb_yn, dotSize_yn = input$dotSize_yn)
            }
          # In the following case, the file is loaded -- to be completed
          } else {
            # loop over the gwas to plot
            # define some outputs -- the first is the names of the gwases to be plotted
            list_for_loop <- all_gwas
            # the second i the actual data to plot
            list_for_results <- list()
            # the third is the color for each gwas
            col_list <- list()
            for (gw in 1:length(list_for_loop)){
              if (gw == 1){
                tmp <- as.character(all_gwas[1])
                dat_tmp <- function.InputManualScroll(dat = dat, window = input$x, input.scroll = input$manual.pos, input.chrom = input$manual.chrom, tmp)
                list_for_results[[gw]] <- dat_tmp
                col_list[[gw]] <- input$col
              } else {
                tmp <- as.character(all_gwas[gw])
                res <- function.manageInput(tmp, path_f, res_example)
                dat_tmp <- as.data.frame(res[[1]])
                chrom_tmp <- as.data.frame(res[[2]])
                dat_tmp <- function.InputManualScroll(dat = dat_tmp, window = input$x, input.scroll = input$manual.pos, input.chrom = input$manual.chrom, tmp)
                if ((nrow(dat_tmp) == 0) || (dat_tmp$chr != as.character(input$manual.chrom))){
                  list_for_results[[gw]] <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 1))
                } else {
                  list_for_results[[gw]] <- dat_tmp
                }
                if (gw == 2){ col_list[[gw]] <- input$col2 } else if (gw == 3){ col_list[[gw]] <- input$col3 } else if (gw == 4){ col_list[[gw]] <- input$col4 } else if (gw == 5){ col_list[[gw]] <- input$col5 }
              }
            }

            # check if there are missing among the gwas to be plotted
            miss = c()
            for (x in list_for_results){
              if (nrow(x) == 0){ miss <- c(miss, 1)}
            }
            if (length(miss) == length(list_for_results)){
              par(mar=c(4, 8, 10, 6))
              plot(0, 0, xlim=c(0, 55), ylim=c(0, 25), bty="n", xaxt='none', yaxt='none', xlab="", ylab="",
                   main=paste0(input$manual.chrom, ":", input$manual.pos), cex.main=2.50, col="white")
              text(x = 27.5, y = 25, labels = "No SNPs in this region.",
                   cex = 2, xpd=T)

            } else {
              # plot
              function.multiPlot(snp.info = list_for_results, list_for_loop, y.lim = input$y, type = input$sel,
                                 plt.type = input$ploType, windows.number = input$sliding.window, smooth.par = input$smooth,
                                 int.locus = target, col_list, input$genV, inpStrVar = input$strVar_inp, recomb_yn = input$recomb_yn, dotSize_yn = input$dotSize_yn)
              #######################################
            }
          }
        }
      }
    }

    # function to try to avoid grey outting
    autoInvalidate <- reactiveTimer(10000)
    observe({
      autoInvalidate()
      cat(".")
      gc()
    })

    # function to try to avoid grey outting
    # main function to plot -- will plot the above function
    output$hist <- renderPlot(plt(input, output, session), res=95)

    # main function to plot the bug-report section
    output$table_report_last <- renderTable({
      all_gwases = c(input$gwas_neuro, input$gwas_cardio, input$gwas_immune, input$gwas_cancer, input$gwas_physio, input$gwas_from_file)
      if (length(all_gwases) == 0){ all_gwases = NA } else { all_gwases = paste0(all_gwases, collapse = ",") }
      pop_interest_ld = c(input$pop_afr, input$pop_amr, input$pop_eas, input$pop_eur, input$pop_sas)
      if (length(pop_interest_ld) == 0){ pop_interest_ld = NA } else { pop_interest_ld = paste0(pop_interest_ld, collapse = ",") }
      # take results from function
      toreport = data.frame("Input datasets" = paste0(all_gwases, collapse = ","), "Reference" = input$genV, "Browsing option" = input$sel, "Target" = input$target, "Window-scale" = input$x, "P-scale" = input$y,
                            "Plot" = input$ploType, "LD" = input$Linkage, "LD populations" = pop_interest_ld)
      colnames(toreport) = c("Input datasets", "Reference", "Browsing option", "Target", "Window-scale", "P-scale", "Plot", "LD", "LD populations")
      return(toreport)
    })

    # here is the rective object in the bug-report page -- this controls a button that send a email with the bug informations
    output$bugss <- renderPlot({
      all_gwases = c(input$gwas_neuro, input$gwas_cardio, input$gwas_immune, input$gwas_cancer, input$gwas_physio, input$gwas_from_file)
      if (length(all_gwases) == 0){ all_gwases = NA } else { all_gwases = paste0(all_gwases, collapse = ",") }
      pop_interest_ld = c(input$pop_afr, input$pop_amr, input$pop_eas, input$pop_eur, input$pop_sas)
      if (length(pop_interest_ld) == 0){ pop_interest_ld = NA } else { pop_interest_ld = paste0(pop_interest_ld, collapse = ",") }
      # take results from function
      toreport = data.frame("Input datasets" = paste0(all_gwases, collapse = ","), "Reference" = input$genV, "Browsing option" = input$sel, "Target" = input$target, "Window-scale" = input$x, "P-scale" = input$y,
                            "Plot" = input$ploType, "LD" = input$Linkage, "LD populations" = pop_interest_ld)
      colnames(toreport) = c("Input datasets", "Reference", "Browsing option", "Target", "Window-scale", "P-scale", "Plot", "LD", "LD populations")
      toreport$Comments = as.character(input$bug_comments)
      if (length(input$confirm_bug) == 1){
        if (input$confirm_bug == "Yes"){
          # save dataset in bug-report folder -- filename will depend on date and time
          fname = paste0("../bug_report/Bug_report_", str_replace_all(Sys.time(), " ", "_"), ".txt")
          write.table(toreport, fname, quote=F, row.names=F, sep = "\t")
          # finally send an email to me
          cmd_mail <- paste("sendEmail -f snpXplorer@gmail.com -t n.tesi@amsterdamumc.nl -u 'AnnotateMe request sent' -m 'Hello, \n a new bug report has been added.' -s smtp.gmail.com:25 -xu snpXplorer@gmail.com -xp snpXplorer22101991!", sep="")
          system(cmd_mail, wait = F)
          plot(0, pch=16, col="white", bty="n", xlab="", ylab="", xaxt="none", yaxt="none", xlim=c(0, 1), ylim=c(0, 1))
          text(x=0.5, y=0.75, labels="Bug reported! Thanks for helping!", font=2, cex=2.25, xpd=T)
          text(x=0.5, y=0.5, labels="You should now refresh the page", font=4, cex=1.50, xpd=T)
        }
      } else {
        plot(0, pch=16, col="white", bty="n", xlab="", ylab="", xaxt="none", yaxt="none", xlim=c(0, 1), ylim=c(0, 1))
      }
  })

    # main function for the annotateMe
    output$annot <- renderPlot({
      ####################
      # ANNOTATE ME PART
      print(input$email)
      if (input$email != "Type your email address..."){
             # Sample 1 number for randomization
             random_num <- sample(x = seq(1, 100000), size = 1, replace = F)
             # When email address is inputed, take snps and save them
             annotateMe.snplist <- unlist(strsplit(input$snp_list, "\n"))
             write.table(annotateMe.snplist, paste("annotateMe_input_", random_num, ".txt", sep=""), quote=F, row.names=F, col.names = F)
             log_filename = paste0("annotateMe_run_", random_num, ".log")
             # Take also input type
             ftype <- as.character(input$snp_list_type)
             if (ftype == "chr:pos (1:12345678)"){ ftype = 1 } else if (ftype == "chr pos (1 12345678)"){ ftype = 2 } else if (ftype == "rsid (rs12345)"){ ftype = 3 }
             # Take the reference genome version
             ref_version = as.character(input$snp_list_reference)
             ref_version = str_split_fixed(ref_version, " ", 2)[, 1]
             # Save user email
             username <- as.character(input$email)
             # Save analysis setting
             analysis_type = as.character(input$analysis_type)
             if (analysis_type == "Gene-set enrichment analysis"){ analysis_type = "enrichment" } else { analysis_type = "mapping" }
             analysis_mode = paste0(as.character(input$analysis_mode), collapse = ",")
             gtex_tissues = paste0(input$gtex_type, collapse = ",")
             # Then run annotate me externally in background -- this depends on the analysis_type requested
             annotateMe.cmd <- paste0("Rscript /root/snpXplorer/AnnotateMe/BIN/MAIN.R annotateMe_input_", random_num, ".txt ", ftype, " ", username, " ", analysis_type, " ", analysis_mode, " ", gtex_tissues, " ", ref_version, " > ", log_filename)
             print(annotateMe.cmd)
             system(annotateMe.cmd, ignore.stdout = F, wait = F)
             # Finally update the email address so that AnnotateMe is not executed every time user load a new page
             #updateTextInput(session, 'email', label="E-mail", value = "Type your email address...")
             plot(0, pch=16, col="white", bty="n", xlab="", ylab="", xaxt="none", yaxt="none", xlim=c(0, 1), ylim=c(0, 1))
             text(x=0.5, y=0.75, labels="Submission completed!", font=2, cex=1.50, xpd=T)
             text(x=0.5, y=0.5, labels="You will receive results by email.", font=4, cex=1, xpd=T)
           } else {
             plot(0, pch=16, col="white", bty="n", xlab="", ylab="", xaxt="none", yaxt="none", xlim=c(0, 1), ylim=c(0, 1))
             text(x=0.5, y=1, labels="How to run AnnotateMe", font=2, cex=1.50, xpd=T)
             text(x=0.5, y=0.5, labels="1) Paste the list of SNPs of interest\n2) Select input type\n3) Insert email address\n4) Submit!", font=4, cex=1)
           }
       })

    # this should be the snp info on the side
    output$table <- renderTable({
      # take results from function
      toplot_table <- plotTable(input)
      # assign colnames and rsid if there are hits
      if (nrow(toplot_table) >0){
        colnames(toplot_table) <- c("Chr", "Position", "-log10(P)", "Study")
        toplot_table = assignRSID(data = toplot_table, genV = input$genV)
        # plot top 10
        top <- head(toplot_table, 10)
      } else {
        top = data.frame(Chr = "No SNPs", Position = "in the", "-log10(P)" = "region of", Study = "interest")
      }
      return(top)
    })

    # this should be the function to report eQTL associations of variants in the region of interest
    output$eQTL <- renderTable({
      # take results from function
      toplot_table <- plotTable(input)
      colnames(toplot_table) <- c("Chr", "Position", "-log10(P)", "Study")
      # assign rsid
      snps_interest = assignRSID(data = toplot_table, genV = input$genV)
      # look for eqtl with function
      if ((input$gtex_type_tb == "Whole_Blood") & !("All_tissues" %in% input$gtex_type_tb)){
        eqtl_in_interval <- plotEqtl_table_blood(snps_interest, eQTL, mapping.ensembl, input$genV)
        eqtl_in_interval$Tissue = "Whole_Blood"
      } else {
        eqtl_in_interval <- plotEqtl_table_all(snps_interest, eQTL, mapping.ensembl, input$genV, input$gtex_type_tb)
      }
      eqtl_in_interval <- head(eqtl_in_interval, 15)
      return(eqtl_in_interval)
    })

    # this should be the cross-reference to gene-cards of the gene of interest or the closest
    output$genecards_link <- renderUI({
      if (input$sel == "Locus, Gene, RsID"){
        # If input is a target region, need to identify which was the input (locus, gene or rsid)
        target <- identiTargetRegion(input$target)
        target_type <- target[[1]]
        # if variant --> shows dbsnp
        if (target_type == "rsid"){
          link <- paste0("https://www.ncbi.nlm.nih.gov/snp/?term=", as.character(target[[2]]))
          url <- a("Click here!", href = link)
          tagList("dbSNP link:", url)
        # if gene --> show gene cards
        } else if (target_type == "gene_name"){
          link <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", as.character(target[[2]]))
          url <- a("Click here!", href = link)
          tagList("GeneCards link:", url)
        }
      }
    })

    # this should be the cross-reference for GWAS catalog specifically
    output$gwascat_link <- renderUI({
      if (input$sel == "Locus, Gene, RsID"){
        # If input is a target region, need to identify which was the input (locus, gene or rsid)
        target <- identiTargetRegion(input$target)
        target_type <- target[[1]]
        # if variant --> shows dbsnp
        if (target_type %in% c("rsid", "gene_name")){
          link <- paste0("https://www.ebi.ac.uk/gwas/search?query=", as.character(target[[2]]))
          url <- a("Click here!", href = link)
          tagList("GWAS-catalog link:", url)
          # if position, use a range of 100bp up/down
        } else if (target_type == "locus"){
          snp_locus <- data.frame(chr=target[[2]], pos=target[[3]])
          link <- paste0("https://www.ebi.ac.uk/gwas/search?query=", paste0(as.character(target[[2]]), ":", as.character(target[[3]]), "-", as.character(target[[3]]+1)))
          url <- a("Click here!", href = link)
          tagList("GWAS-catalog link:", url)
        }
      }
    })

    # this should be the cross-reference to LD-hub (always nice to put it)
    output$ld_hub_link <- renderUI({
      link <- paste0("http://ldsc.broadinstitute.org")
      url <- a("Click here!", href = link)
      tagList("LD Hub link:", url)
    })

    # this should be table of the studies included
    output$table_info <- renderTable({
      studies_table <- as.data.frame(matrix(data=NA, nrow=18, ncol = 5))
      colnames(studies_table) <- c("Class", "Name", "Trait", "Authors", "Reference")
      studies_table[1, ] <- c("Neurological", "Alzheimer_million", "Alzheimer's disease", "Wightman et al., 2021", "https://doi.org/10.1038/s41588-021-00921-z")
      studies_table[2, ] <- c("Neurological", "GR@ACE", "Alzheimer's disease", "De Rojas et al., 2021", "https://doi.org/10.1038/s41467-021-22491-8")
      studies_table[3, ] <- c("Neurological", "IGAP", "Alzheimer's disease", "Kunkle et al., 2019", "https://doi.org/10.1038/s41588-019-0358-2")
      studies_table[4, ] <- c("Neurological", "proxy_AD", "by-proxy Alzheimer's", "Jansen et al., 2019", "https://doi.org/10.1038/s41588-018-0311-9")
      studies_table[5, ] <- c("Neurological", "Autism", "Autism", "Matoba et al., 2020", "https://doi.org/10.1038/s41398-020-00953-9")
      studies_table[6, ] <- c("Neurological", "Depression", "Depression", "Cai et al., 2020", "https://doi.org/10.1038/s41588-020-0594-5")
      studies_table[7, ] <- c("Cardiovascular", "CAD", "Coronary artery disease", "van der Harst et al., 2017", "https://doi.org/10.1161/CIRCRESAHA.117.312086")
      studies_table[8, ] <- c("Cardiovascular", "CAD_Diabetics", "Coronary artery disease in diabetes", "Fall et al., 2018", "https://doi.org/10.1007/s00125-018-4686-z")
      studies_table[9, ] <- c("Cardiovascular", "Ventricular volumne", "Ventricular volume", "Vojinovic et al., 2018", "https://doi.org/10.1038/s41467-018-06234-w")
      studies_table[10, ] <- c("Cardiovascular", "SBP", "Sistolic blood pressure", "Evangelou et al., 2018", "https://doi.org/10.1038/s41588-018-0205-x")
      studies_table[11, ] <- c("Cardiovascular", "BMI", "Body-mass index", "Yengo et al., 2018", "https://doi.org/10.1093/hmg/ddy271")
      studies_table[12, ] <- c("Cardiovascular", "Diabetes", "Type I Diabetes Mellitus", "Forgetta et al., 2020", "https://doi.org/10.2337/db19-0831")
      studies_table[13, ] <- c("Immunological", "COVID", "COVID-19 critical illness", "Erola Pairo-Castineira et al. 2020", "https://doi.org/10.1038/s41586-020-03065-y")
      studies_table[14, ] <- c("Immunological", "Lupus", "Systemic Lupus", "Yong-Fei Wang et al. 2021", "https://doi.org/10.1038/s41467-021-21049-y")
      studies_table[15, ] <- c("Immunological", "Inflammation", "Inflammatory Biomarkers", "Sanni E. Ruotsalainen et al., 2020", "https://doi.org/10.1038/s41431-020-00730-8")
      studies_table[16, ] <- c("Immunological", "Asthma", "Asthma", "Han et al., 2020", "https://doi.org/10.1038/s41467-020-15649-3")
      studies_table[17, ] <- c("Cancer", "Breast_cancer", "Breast cancer", "Zhang et al., 2020", "https://doi.org/10.1038/s41588-020-0609-2")
      studies_table[18, ] <- c("Cancer", "Myeloproliferative", "Myeloproliferative neoplasm", "Erik L. Bao et al., 2020", "https://doi.org/10.1038/s41586-020-2786-7")
      studies_table[19, ] <- c("Cancer", "Prostate", "Prostate cancer", "Peter N. Fiorica et al., 2020", "https://doi.org/10.1371/journal.pone.0236209")
      studies_table[20, ] <- c("Cancer", "Lung", "Lung cancer", "Sara R. Rashkin et al., 2020", "https://doi.org/10.1038/s41467-020-18246-6")
      studies_table[21, ] <- c("Cancer", "Leukemia", "Lymphocytic Leukemia", "Sara R. Rashkin et al., 2020", "https://doi.org/10.1038/s41467-020-18246-6")
      studies_table[22, ] <- c("Physiological", "Multivariate_Longevity", "Longevity/Lifespan/Healthspan", "Timmers et al., 2020", "https://doi.org/10.1038/s41467-020-17312-3")
      studies_table[23, ] <- c("Physiological", "UKBaging", "Parental longevity", "Timmers et al., 2019", "https://doi.org/10.7554/eLife.39856")
      studies_table[24, ] <- c("Physiological", "Height", "Height", "Yengo et al., 2018", "https://doi.org/10.1093/hmg/ddy271")
      studies_table[25, ] <- c("Physiological", "Education", "Education", "Perline A. Demange et al., 2021", "https://doi.org/10.1038/s41588-020-00754-2")
      studies_table[26, ] <- c("Physiological", "Bone density", "Bone density and fracture risk", "Ida Surakka et al., 2020", "https://doi.org/10.1038/s41467-020-17315-0")
      studies_table[27, ] <- c("Physiological", "Vitamin D", "Vitamin D", "Manousaki et al., 2020", "https://doi.org/10.1016/j.ajhg.2020.01.017")
      return(studies_table)
    }, width = "100%")

    # this should be table of the additional information (1000Genome data)
    output$table_info_1000G <- renderTable({
      studies_table <- as.data.frame(matrix(data=NA, nrow=26, ncol = 6))
      colnames(studies_table) <- c("Population code", "Population Name", "Superpopulation code", "Superpopulation name", "Number of individuals", "Percentage of females")
      studies_table[1, ] <- c("FIN", "Finnish", "EUR", "European Ancestry", "103", "62.13%")
      studies_table[2, ] <- c("GBR", "British", "EUR", "European Ancestry", "104", "52.88%")
      studies_table[3, ] <- c("CEU", "CEPH", "EUR", "European Ancestry", "183", "51.36%")
      studies_table[4, ] <- c("TSI", "Toscani", "EUR", "European Ancestry", "112", "49.11%")
      studies_table[5, ] <- c("IBS", "Iberian", "EUR", "European Ancestry", "157", "48.40%")
      studies_table[6, ] <- c("ACB", "African Caribbean", "AFR", "African Ancestry", "121", "49.58%")
      studies_table[7, ] <- c("ASW", "African Ancestry SW", "AFR", "African Ancestry", "111", "55.85%")
      studies_table[8, ] <- c("ESN", "Esan", "AFR", "African Ancestry", "100", "47.00%")
      studies_table[9, ] <- c("GWD", "Gambian Mandinka", "AFR", "African Ancestry", "113", "51.32%")
      studies_table[10, ] <- c("LWK", "Luhya", "AFR", "African Ancestry", "113", "53.45%")
      studies_table[11, ] <- c("MSL", "Mende", "AFR", "African Ancestry", "89", "51.68%")
      studies_table[12, ] <- c("YRI", "Yoruba", "AFR", "African Ancestry", "186", "45.70%")
      studies_table[13, ] <- c("CLM", "Colombian", "AMR", "American Ancestry", "132", "56.06%")
      studies_table[14, ] <- c("MXL", "Mexican", "AMR", "American Ancestry", "104", "54.81%")
      studies_table[15, ] <- c("PEL", "Peruvian", "AMR", "American Ancestry", "122", "55.37%")
      studies_table[16, ] <- c("PUR", "Puerto Rican", "AMR", "American Ancestry", "139", "49.64%")
      studies_table[17, ] <- c("CDX", "Dai Chinese", "EAS", "East Asian Ancestry", "102", "49.02%")
      studies_table[18, ] <- c("CHB", "Han Chinese", "EAS", "East Asian Ancestry", "108", "55.55%")
      studies_table[19, ] <- c("CHS", "Southern Han Chinese", "EAS", "East Asian Ancestry", "165", "46.67%")
      studies_table[20, ] <- c("JPT", "Japanese", "EAS", "East Asian Ancestry", "105", "45.14%")
      studies_table[21, ] <- c("KHV", "KHV Kinh Vietnamese", "EAS", "East Asian Ancestry", "121", "50.41%")
      studies_table[22, ] <- c("BEB", "Bengali", "SAS", "South Asian Ancestry", "88", "50.00%")
      studies_table[23, ] <- c("GIH", "Gujarati", "SAS", "South Asian Ancestry", "113", "45.13%")
      studies_table[24, ] <- c("ITU", "Telugu", "SAS", "South Asian Ancestry", "103", "41.75%")
      studies_table[25, ] <- c("PJL", "Punjabi", "SAS", "South Asian Ancestry", "113", "47.79%")
      studies_table[26, ] <- c("STU", "Tamil", "SAS", "South Asian Ancestry", "105", "46.67%")
      return(studies_table)
    }, width = "100%")

    # this should be table of the additional information (structural variant datasets)
    output$table_info_SVs <- renderTable({
      studies_table <- as.data.frame(matrix(data=NA, nrow=3, ncol = 3))
      colnames(studies_table) <- c("Study", "Experimental technology", "Reference")
      studies_table[1, ] <- c("Linthorst et al., 2020", "PacBio CLR", "https://doi.org/10.1038/s41398-020-01060-5")
      studies_table[2, ] <- c("Chaisson et al., 2019", "PacBio CLR + Illumina WGS + 10X Genomics + Hi-C", "https://doi.org/10.1038/s41467-018-08148-z")
      studies_table[3, ] <- c("Audano et al., 2019", "PacBio CLR", "https://doi.org/10.1016/j.cell.2018.12.019")
      return(studies_table)
    }, width = "100%")

    # this should be the link to bioRxiv paper
    output$biorxiv_link <- renderUI({
      link <- paste0("https://academic.oup.com/nar/article/49/W1/W603/6287842")
      url <- a("Click here!", href = link)
      tagList("Link to paper:", url)
    })

    # this should be the link to longevity paper
    output$longevity_ms <- renderUI({
      link <- paste0("https://academic.oup.com/biomedgerontology/article/76/5/750/5996044")
      url <- a("Click here!", href = link)
      tagList("Link to the paper:", url)
    })

    # this should be the link to rotation paper
    output$rotation_ms <- renderUI({
      link <- paste0("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8724252/")
      url <- a("Click here!", href = link)
      tagList("Link to the paper:", url)
    })

    # this is the pdf viewer of the documentation
    output$pdf_doc_view <- renderUI({ tags$iframe(style="height:1000px; width:95%", src="snpXplorer_documentation.pdf") })

    # this is the pdf viewer of the biorxiv paper
    output$pdf_biorxiv_view <- renderUI({ tags$iframe(style="height:1000px; width:95%", src="gkab410-3.pdf") })

    # this is to download all SVs in the region
    output$download_SVs <- downloadHandler(
      filename = function() {
        paste0("All_SVs_region.txt")
      },
      content = function(file) {
        toplot_table <- plotTable_SVs(input)
        #colnames(toplot_table) <- c("Chr", "Position", "-log10(P)", "Study")
        write.table(toplot_table, file, row.names = FALSE, quote=F, sep=" ")
      }
    )

    # this is to download all eQTLs in the region
    output$download_eQTL <- downloadHandler(
      filename = function() {
        paste0("All_eQTLs_region.txt")
      },
      content = function(file) {
        # take results from function
        snps_in_interval <- plotTable(input)
        # look for eqtl with function
        if ((input$gtex_type_tb == "Whole_Blood") & !("All_tissues" %in% input$gtex_type_tb)){
          eqtl_in_interval <- plotEqtl_table_blood(snps_in_interval, eQTL, mapping.ensembl, input$genV)
          eqtl_in_interval$Tissue = "Whole_Blood"
        } else {
          eqtl_in_interval <- plotEqtl_table_all(snps_in_interval, eQTL, mapping.ensembl, input$genV, input$gtex_type_tb)
        }
        eqtl_in_interval <- head(eqtl_in_interval, 15)
        return(eqtl_in_interval)
      }
    )

    # this is to download the LD table
    output$download_LDtable <- downloadHandler(
      filename = function() {
        paste0("LD_table.txt")
      },
      content = function(file) {
        toplot_table <- plotTable_LD(input)
        write.table(toplot_table, file, row.names = FALSE, quote=F, sep=" ")
      }
    )

    # this is to download all snp associations
    output$download_SNPs <- downloadHandler(
      filename = function() {
        paste0("All_associations_region.txt")
      },
      content = function(file) {
        toplot_table <- plotTable(input)
        colnames(toplot_table) <- c("Chr", "Position", "-log10(P)", "Study")
        # assign rsid
        toplot_table = assignRSID(data = toplot_table, genV = input$genV)
        write.table(toplot_table, file, row.names = FALSE, quote=F, sep=" ")
      }
    )

    # this is to download sample data -- exploration file 1
    help_explo_1 <- fread("www/trial_dataset_website.txt.gz")
    output$download_help_exploration <- downloadHandler(
      filename = function() {
        paste0("trial_dataset_snpXplorer_exploration.txt")
      },
      content = function(file) {
        write.table(help_explo_1, file, row.names = FALSE, quote=F, sep=" ")
      }
    )

    # this is to download sample data -- exploration file 2
    help_explo_2 <- fread("www/trial_plink_snpxplorer.glm.logistic.gz")
    output$download_help_exploration_plink <- downloadHandler(
      filename = function() {
        paste0("trial_dataset_snpXplorer_exploration_PLINK.txt")
      },
      content = function(file) {
        write.table(help_explo_2, file, row.names = FALSE, quote=F, sep="\t")
      }
    )

    # this is to download sample data -- annotation file 1
    help_annot_1 <- fread("www/trial_annotation_f1.txt", h=F)
    output$download_help_annotation1 <- downloadHandler(
      filename = function() {
        paste0("trial_dataset_snpXplorer_annotation_type1.txt")
      },
      content = function(file) {
        write.table(help_annot_1, file, row.names = FALSE, col.names=F, quote=F)
      }
    )

    # this is to download sample data -- annotation file 1
    help_annot_2 <- fread("www/trial_annotation_f2.txt", h=F)
    output$download_help_annotation2 <- downloadHandler(
      filename = function() {
        paste0("trial_dataset_snpXplorer_annotation_type2.txt")
      },
      content = function(file) {
        write.table(help_annot_2, file, row.names = FALSE, col.names=F, quote=F)
      }
    )

    # this is to download sample data -- annotation file 1
    help_annot_3 <- fread("www/trial_annotation_f3.txt", h=F)
    output$download_help_annotation3 <- downloadHandler(
      filename = function() {
        paste0("trial_dataset_snpXplorer_annotation_type3.txt")
      },
      content = function(file) {
        write.table(help_annot_3, file, row.names = FALSE, col.names=F, quote=F)
      }
    )

    # this is the link to the github page
    url <- a("Visit our Github page", href="https://github.com/TesiNicco/SNPbrowser")
    output$github <- renderUI({ tagList(url) })

    # below all the links to youtube videos tutorials
    output$video1 <- renderUI({
      HTML(paste0('<iframe width="560" height="315" src="https://www.youtube.com/embed/z8kqzTacaS4" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>'))
    })
    output$video2 <- renderUI({
      HTML(paste0('<iframe width="560" height="315" src="https://www.youtube.com/embed/Ai3F-JBQL3U" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>'))
    })
    output$video3 <- renderUI({
      HTML(paste0('<iframe width="560" height="315" src="https://www.youtube.com/embed/eOJTf7tk2Rg" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>'))
    })
    output$video4 <- renderUI({
      HTML(paste0('<iframe width="560" height="315" src="https://www.youtube.com/embed/QP-5XjIEYpI" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>'))
    })

    # this should be to download image
    output$downloadPlot <- downloadHandler(
      file = function() { paste(input$target, "_", paste0(c(input$gwas, input$gwas_from_file), collapse="_"), ".png", sep="") },
      content = function(file) {
        #ggsave(p(), filename = file)
        png(file = file, height = 20, width = 13, res=600, units="in")
        plt(input, output, session)
        dev.off()
      })
})
