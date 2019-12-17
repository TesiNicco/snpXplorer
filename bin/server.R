#SERVER SNP BROWSER

#LIBRARIES
library(shiny)
library(data.table)
library(stringr)
library(ggplot2)
library(colourpicker)
library(rvest)
library(stringr)

#FUNCTIONS
#basic function to plot stuff 
function.plot <- function(snp.info, y.lim, type, plt.type, windows.number, smooth.par, int.locus, gwas, colorPoint){
  par(mar=c(5, 5, 4, 5))
  
  #these are standard values for minimun and maximum recombination rates genome-wide for hg19
  chromatin.lower <- 0
  chromatin.upper <- 100
  
  #prepare title
  title <- function.title(gwas, int.locus, type, snp.info)
  
  #assign dot sizes (function for this)
  snp.info <- function.pointSize(dat = snp.info, range = seq(2, 7, 0.5))

  #independently from plot type, I need the genes that are in the window to adjust y-axis -- here it is
  genes <- function.dynamicGene(snp.info = snp.info)
  if (nrow(genes) > 0){
    min.y <- min(genes$y)
  } else {
    min.y <- 0
  }
  
  #get recombination tracks, function for this
  recomb <- findRecomb(snp.info)
  
  if (plt.type == "Points"){
    #empty plot
    plot(x = 0, y = 0, xlab='Chromosomal position (Mb)', cex.lab=1.75, xaxt='none',
         ylab="", ylim=c(min.y*y.lim/12, y.lim), cex.axis = 1.5, 
         pch=16, col="white", cex=2, type = "p", xaxs="i", yaxt='none', xlim=c(min(snp.info$pos), max(snp.info$pos)), main=title, cex.main=2.50, bty='n')
    
    #add grid
    for (x in seq(0, y.lim, (y.lim-min.y*y.lim/10)/10)){abline(h=x, lwd=0.4, col="grey80")}
    for (x in seq(min(snp.info$pos), max(snp.info$pos), (max(snp.info$pos)-min(snp.info$pos))/10)){
      segments(x0 = x, y0 = 0, x1 = x, y1 = y.lim, col = "grey80", lwd=0.4)
    }
    
    #add recombination rates: for this, need to normalize between 0 and max.y the recombination rate
    recomb$norm.rate <- (y.lim - 1) * ((recomb$"Rate(cM/Mb)" - chromatin.lower) / (chromatin.upper - chromatin.lower))
    y.axis.recomb <- seq(0, chromatin.upper, chromatin.upper/4)
    y.axis.norm <- (y.lim - 1) * ((y.axis.recomb - min(y.axis.recomb))/(max(y.axis.recomb) - min(y.axis.recomb)))
    points(recomb$"Position(bp)", recomb$norm.rate, type="l", lwd=1.5, col="darkolivegreen3")
    
    #add significance lines and corresponding legend -- for now two at 0.05 and genome-wide 5e-8
    abline(h=-log10(0.05), lty=2, col=alpha("darkgreen", 1))
    abline(h=-log10(5e-8), lty=2, col=alpha("purple", 1))
    legend("topleft", bty='n', legend = c("p=0.05", "p=5e-8"), lty=c(2,2), lwd=c(2,2), col=c("darkgreen", "purple"), cex=1.5, ncol = 2)
    
    #then points
    snp.info$"-log10(P-value)" <- -log10(as.numeric(snp.info$p))
    points(x = snp.info$pos, y = snp.info$"-log10(P-value)", cex.lab=1.5, xaxt='none',
           pch=16, col=alpha(colorPoint, 0.6), cex=snp.info$size, type = "p", xaxs="i", yaxt='none', xlim=c(min(snp.info$pos), max(snp.info$pos)))
    
    #if input type was a single snp (either position or rs id), then color the searched variant differently
    if (type %in% c("Position", "Rs ID")){
      #restrict to snp of interest, then plot it
      if (int.locus != "Type position..."){
        snp.interest <- snp.info[which(snp.info$locus == int.locus),]
        points(x = snp.interest$pos, y = snp.interest$"-log10(P-value)", pch=23, lwd=1.5, col=alpha(colorPoint, 0.8), cex=snp.interest$size)
      }
    }
    
    #manage axes
    axes <- seq(min(snp.info$pos), max(snp.info$pos), (max(snp.info$pos) - min(snp.info$pos))/7)
    axes.labels <- round(axes/1000000, 3)
    axis(side = 1, at=axes, cex.axis=1.5, labels=axes.labels)
    axes.x <- ceiling(seq(0, y.lim, y.lim/5))
    axis(side = 2, at = axes.x, labels = axes.x, cex.axis=1.5)
    
    #axis for recombination rates on the right
    axis(side = 4, at = y.axis.norm, labels=seq(0, 100, 25), col='darkolivegreen3', col.axis = "darkolivegreen3", cex.axis=1.5, xpd=T)
    text(x = max(snp.info$pos), y = y.lim/5*4, "Recombination rate (cM/Mb)",srt = -90, 
         col='darkolivegreen3', xpd=T, pos = 4, offset = 4, cex=1.5, font=2)
    text(x = min(snp.info$pos), y = y.lim/3*2, "-Log10(P-value)", srt = 90, 
         xpd=T, pos = 2, offset = 4, cex=1.5, font=2)
    
    if (min.y != 0){
      #manage gene names
      for (g in 1:nrow(genes)){
        #main gene line -- full transcription sequence
        segments(x0 = genes$txStart[g], y0 = genes$y[g]*y.lim/12, x1 = genes$txEnd[g], y1 = genes$y[g]*y.lim/12, lwd=3)
        #need to divide exones from introns
        start <- as.data.frame(t(str_split_fixed(genes$exonStarts[g], ",", genes$exonCount[g]+1)))
        end <- as.data.frame(t(str_split_fixed(genes$exonEnds[g], ",", genes$exonCount[g]+1)))
        exons <- cbind(start, end)
        colnames(exons) <- c("start", "end")
        exons$start <- as.numeric(as.character(exons$start))
        exons$end <- as.numeric(as.character(exons$end))
        #main loop over exons
        for (j in 1:nrow(exons)){
          rect(xleft=exons$start[j], ybottom=genes$y[g]*y.lim/12-(y.lim/4/4*0.15), xright=exons$end[j], 
               ytop = genes$y[g]*y.lim/12+(y.lim/4/4*0.15), col='grey80')
        }
        text(x = genes$txStart[g] + (genes$txEnd[g] - genes$txStart[g])/2, y = genes$y[g]*y.lim/12 + (y.lim/4/4*0.5), 
             labels=genes$"#geneName"[g], font=3, cex=1)
        #        if (genes$strand[g] == "+"){
        #          arrows(x0 = genes$txEnd[g], y0 = genes$y[g], x1 = genes$txEnd[g] + (10*2/100), y1 = genes$y[g], 
        #                 length=0.1, lwd=2, col='coral')
        #        } else if (genes$strand[g] == "-"){
        #          arrows(x0 = genes$txStart[g], y0 = genes$y[g], x1 = genes$txStart[g] - (10*2/100), y1 = genes$y[g], 
        #                 length=0.1, lwd=2, col='coral')
        #        }
      }
    }
    
  } else {
    #Sliding window approach -- then loess on sliding window values
    out <- function.DensityLinePvalue(snp.info = snp.info, wind.n = windows.number)
    lo <- loess(out$pvalue ~ out$window, span=smooth.par)
    xl <- seq(min(out$window), max(out$window), (max(out$window) - min(out$window))/100)
    pred <- predict(lo, xl)
    pred[which(pred < 0)] <- 0
    #add limits -- left and right for polygon function
    xl <- c(xl[1], xl, xl[length(xl)])
    pred <- c(0, pred, 0)
    #lines(x = xl, y = pred, col='navy', lwd=4)
    
    #main plot
    plot(x = 0, y = 0, xlab='Chromosomal position (Mb)', cex.lab=1.5, xaxt='none',
         ylab="-log10(P-value)", ylim=c(min.y*y.lim/12, y.lim), cex.axis = 1.25, xaxs="i", yaxt='none',
         pch=1, col="white", cex=1.75, type = "h", lwd=2, xlim=c(min(xl), max(xl)))
    
    #add grid
    for (x in seq(0, y.lim, (y.lim-min.y*y.lim/10)/10)){abline(h=x, lwd=0.4, col="grey80")}
    for (x in seq(min(snp.info$pos), max(snp.info$pos), (max(snp.info$pos)-min(snp.info$pos))/10)){
      segments(x0 = x, y0 = 0, x1 = x, y1 = y.lim, col = "grey80", lwd=0.4)
    }
    
    polygon(x = xl, y = pred, col = alpha(colorPoint, 0.6), lwd=2, xaxs="i")
    
    #manage axes
    axes <- seq(min(snp.info$pos), max(snp.info$pos), (max(snp.info$pos) - min(snp.info$pos))/12)
    axes.labels <- round(axes/1000000, 2)
    axis(side = 1, at=axes, cex.axis=1.5, labels=axes.labels)
    axes.x <- ceiling(seq(0, y.lim, y.lim/5))
    axis(side = 2, at = axes.x, labels = axes.x)
    
    
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
    
    #manage gene names
    for (g in 1:nrow(genes)){
      #main gene line -- full transcription sequence
      segments(x0 = genes$txStart[g], y0 = genes$y[g]*y.lim/12, x1 = genes$txEnd[g], y1 = genes$y[g]*y.lim/12, lwd=3)
      #need to divide exones from introns
      start <- as.data.frame(t(str_split_fixed(genes$exonStarts[g], ",", genes$exonCount[g]+1)))
      end <- as.data.frame(t(str_split_fixed(genes$exonEnds[g], ",", genes$exonCount[g]+1)))
      exons <- cbind(start, end)
      colnames(exons) <- c("start", "end")
      exons$start <- as.numeric(as.character(exons$start))
      exons$end <- as.numeric(as.character(exons$end))
      #main loop over exons
      for (j in 1:nrow(exons)){
        rect(xleft=exons$start[j], ybottom=genes$y[g]*y.lim/12-(y.lim/4/4*0.15), xright=exons$end[j], 
             ytop = genes$y[g]*y.lim/12+(y.lim/4/4*0.15), col='grey80')
      }
      text(x = genes$txStart[g] + (genes$txEnd[g] - genes$txStart[g])/12, y = genes$y[g]*y.lim/12 + (y.lim/4/4*0.5), 
           labels=genes$"#geneName"[g], font=3, cex=1)
      #     if (genes$strand[g] == "+"){
      #        arrows(x0 = genes$txEnd[g], y0 = genes$y[g], x1 = genes$txEnd[g] + (10*2/100), y1 = genes$y[g], 
      #               length=0.1, lwd=2, col='coral')
      #      } else if (genes$strand[g] == "-"){
      #        arrows(x0 = genes$txStart[g], y0 = genes$y[g], x1 = genes$txStart[g] - (10*2/100), y1 = genes$y[g], 
      #               length=0.1, lwd=2, col='coral')
      #      }
    }
  }
}

#read recombination map and give coordinates of interests as output
findRecomb <- function(snp.info){
  chrom <- snp.info$chr[1]
  min.p <- min(snp.info$pos)
  max.p <- max(snp.info$pos)
  
  #read chromosome file
  inpf <- paste("../data/databases/genetic_map_HapMapII_GRCh37/genetic_map_GRCh37_chr", chrom, ".txt.gz", sep="")
  rf <- fread(inpf, h=T)
  
  #extract interval of interest
  recomb <- rf[which(rf$"Position(bp)" >= min.p & rf$"Position(bp)" <= max.p),]
  
  return(recomb)
  
}

#function to assign dot size -- define range
function.pointSize <- function(dat, range){
  dat$size <- 2
  
  for (x in range){
    dat$size[which(-log10(as.numeric(dat$p)) >= x)] <- x
  }
  
  return(dat)
}

#function to manage titles
function.title <- function(gwas, int.locus, type, snp.info){
  t = ""
  if (length(gwas) == 1){
    if (gwas == "example"){t = "Example ~ IGAP"} else {t = gwas}
  
    if (type %in% c("Position", "RsID")){
      if (int.locus %in% c("Type position...", "Type variant identifier...")){
        title = paste(t, " ~ chr", snp.info$chr[1], ":", floor(min(snp.info$pos) + (max(snp.info$pos) - min(snp.info$pos))/2), sep = "")
      } else {
        title = paste(t, " ~ chr", int.locus, sep = "")
      }
    } else if (type == "Gene"){
      if (int.locus == "Type gene symbol..."){
        title = paste(t, " ~ chr", snp.info$chr[1], ":", floor(min(snp.info$pos) + (max(snp.info$pos) - min(snp.info$pos))/2), sep = "")
      } else {
        title = paste(t, " ~ ", toupper(as.character(int.locus)), " ~ chr", snp.info$chr[1], sep="")
      }
    } else if (type == "Manual scroll"){
      title <- paste(t, " ~ chr", snp.info$chr[1], ":", floor(min(snp.info$pos) + (max(snp.info$pos) - min(snp.info$pos))/2), sep = "")
    }
  } else {
    #else=multiplot
    if (type %in% c("Position", "RsID")){
      if (int.locus == "Type position..."){
        title = paste(gwas[1], " vs. ", gwas[2], ": chr", snp.info$chr[1], ":", floor(min(snp.info$pos) + (max(snp.info$pos) - min(snp.info$pos))/2), sep = "")
      } else {
        title = paste(gwas[1], " vs. ", gwas[2], ": chr", int.locus, sep = "")
      }
    } else if (type == "Gene"){
      print(int.locus)
      if (int.locus == "Type gene symbol..."){
        title = paste(gwas[1], " vs. ", gwas[2], ": chr", snp.info$chr[1], ":", floor(min(snp.info$pos) + (max(snp.info$pos) - min(snp.info$pos))/2), sep = "")
      } else {
        title = paste(gwas[1], " vs. ", gwas[2], " ~ chr", snp.info$chr[1], sep = "")
      }
    } else if (type == "Manual scroll"){
      title <- paste(gwas[1], " vs. ", gwas[2], " ~ chr", snp.info$chr[1], sep="")
    }
  }
  
  return(title)
}

#plot in case of multiple files
function.multiPlot <- function(snp.info, snp.info.f2, y.lim, type, plt.type, windows.number, smooth.par, int.locus, lab1, lab2, 
                               colorPoint, colorPoint2){
  par(mar=c(5, 5, 4, 1))
  
  #prepare title
  title <- function.title(gwas = c(lab1, lab2), int.locus, type, snp.info)

  #assign dot sizes (function for this)
  snp.info <- function.pointSize(dat = snp.info, range = seq(2, 7, 0.5))
  snp.info.f2 <- function.pointSize(dat = snp.info.f2, range = seq(1.5, 7, 0.5))

  #independently from plot type, I need the genes that are in the window to adjust y-axis -- here it is
  genes <- function.dynamicGene(snp.info = snp.info)
  if (nrow(genes) > 0){
    min.y <- min(genes$y)
  } else {
    min.y <- 0
  }

  if (plt.type == "Points"){
    plot(x = 0, y = 0, xlab='Chromosomal position (Mb)', cex.lab=1.75, xaxt='none',
         ylab="-log10(P-value)", ylim=c(min(genes$y)*y.lim/12, y.lim), cex.axis = 1.5, bty='n',
         pch=16, col="white", cex=2, type = "p", xaxs="i", yaxt='none', xlim=c(min(snp.info$pos), max(snp.info$pos)), cex.main=2.50, bty='n', main=title)
    
    #add grid
    for (x in seq(0, y.lim, (y.lim-min.y*y.lim/10)/10)){abline(h=x, lwd=0.4, col="grey80")}
    for (x in seq(min(snp.info$pos), max(snp.info$pos), (max(snp.info$pos)-min(snp.info$pos))/10)){
      segments(x0 = x, y0 = 0, x1 = x, y1 = y.lim, col = "grey80", lwd=0.4)
    }

    #add significance lines and corresponding legend -- for now two at 0.05 and genome-wide 5e-8
    abline(h=-log10(0.05), lty=2, col=alpha("darkgreen", 1))
    abline(h=-log10(5e-8), lty=2, col=alpha("purple", 1))
    legend("topleft", bty='n', legend = c("p=0.05", "p=5e-8"), lty=c(2,2), lwd=c(2,2), col=c("darkgreen", "purple"), cex=1.5, ncol = 2)
    
    #then points - gwas 1
    snp.info$"-log10(P-value)" <- -log10(snp.info$p)
    points(x = snp.info$pos, y = snp.info$"-log10(P-value)",
           pch=16, col=alpha(colorPoint, 0.6), cex=snp.info$size, type = "p", xaxs="i")
    
    #add additional file data points -- gwas 2
    snp.info.f2$"-log10(P-value)" <- -log10(snp.info.f2$p)
    points(x = snp.info.f2$pos, y = snp.info.f2$"-log10(P-value)", pch=16, col=alpha(colorPoint2, 0.6), 
           cex=snp.info.f2$size, xaxs="i", type="p")
    
    #manage axes
    axes <- seq(min(snp.info$pos), max(snp.info$pos), (max(snp.info$pos) - min(snp.info$pos))/7)
    axes.labels <- round(axes/1000000, 3)
    axis(side = 1, at=axes, cex.axis=1.5, labels=axes.labels)
    axes.x <- ceiling(seq(0, y.lim, y.lim/5))
    axis(side = 2, at = axes.x, labels = axes.x, cex.axis=1.5)
    
    #if input type was a single snp (either position or rs id), then color the searched variant differently
    snp.interest.flag <- 0
    if (type %in% c("Position", "Rs ID")){
      #restrict to snp of interest, then plot it
      if (int.locus != "Type position..."){
        snp.interest.flag <- 1
        snp.interest <- snp.info[which(snp.info$locus == int.locus),]
        snp.interest.addF <- snp.info.f2[which(snp.info.f2$locus == int.locus),]
        
        points(x = snp.interest$pos, y = -log10(snp.interest$p), pch=23, col="black", lwd=2, bg=alpha(colorPoint, 1), cex=snp.interest$size)
        points(x = snp.interest.addF$pos, y = -log10(as.numeric(snp.interest.addF$p)), pch=23, col="black", lwd=2, bg=alpha(colorPoint2, 1), cex=snp.interest.addF$size)
      }
    }
    
    #add legend now
    if (snp.interest.flag == 0){
      legend("topright", legend = c(lab1, lab2), col = c(colorPoint, colorPoint2), pch=16, 
             cex=1.50, ncol=2, xpd=T, bty='n')
    } else {
      legend("topright", legend = c(lab1, paste(lab1, " - Input", sep=""), lab2, paste(lab2, " - Input", sep="")), 
             col = c(colorPoint, colorPoint, colorPoint2, colorPoint2), pch=c(16, 23, 16, 23), pt.lwd = 1.5, 
             cex=1.50, ncol=2, xpd=T, bty='n')
    }
    
    if (min.y != 0){
      #manage gene names
      for (g in 1:nrow(genes)){
        #main gene line -- full transcription sequence
        segments(x0 = genes$txStart[g], y0 = genes$y[g]*y.lim/12, x1 = genes$txEnd[g], y1 = genes$y[g]*y.lim/12, lwd=3)
        #need to divide exones from introns
        start <- as.data.frame(t(str_split_fixed(genes$exonStarts[g], ",", genes$exonCount[g]+1)))
        end <- as.data.frame(t(str_split_fixed(genes$exonEnds[g], ",", genes$exonCount[g]+1)))
        exons <- cbind(start, end)
        colnames(exons) <- c("start", "end")
        exons$start <- as.numeric(as.character(exons$start))
        exons$end <- as.numeric(as.character(exons$end))
        #main loop over exons
        for (j in 1:nrow(exons)){
          rect(xleft=exons$start[j], ybottom=genes$y[g]*y.lim/12-(y.lim/4/4*0.15), xright=exons$end[j], 
               ytop = genes$y[g]*y.lim/12+(y.lim/4/4*0.15), col='grey80')
        }
        text(x = genes$txStart[g] + (genes$txEnd[g] - genes$txStart[g])/2, y = genes$y[g]*y.lim/12 + (y.lim/4/4*0.5), 
             labels=genes$"#geneName"[g], font=3, cex=1)
        #        if (genes$strand[g] == "+"){
        #          arrows(x0 = genes$txEnd[g], y0 = genes$y[g]*y.lim/10, x1 = genes$txEnd[g] + (10*2/100), y1 = genes$y[g]*y.lim/10, 
        #                 length=0.1, lwd=2, col='coral')
        #        } else if (genes$strand[g] == "-"){
        #          arrows(x0 = genes$txStart[g], y0 = genes$y[g]*y.lim/10, x1 = genes$txStart[g] - (10*2/100), y1 = genes$y[g]*y.lim/10, 
        #                 length=0.1, lwd=2, col='coral')
        #        }
      }
    }    
  } else {
    #Sliding window approach -- then loess on sliding window values -- this is for data 1
    out <- function.DensityLinePvalue(snp.info = snp.info, wind.n = windows.number)
    lo <- loess(out$pvalue ~ out$window, span=smooth.par)
    xl <- seq(min(out$window), max(out$window), (max(out$window) - min(out$window))/100)
    pred <- predict(lo, xl)
    pred[which(pred < 0)] <- 0
    #add limits -- left and right for polygon function
    xl <- c(xl[1], xl, xl[length(xl)])
    pred <- c(0, pred, 0)
    #lines(x = xl, y = pred, col='navy', lwd=4)
    
    #Sliding window approach -- then loess on sliding window values -- this is for data 2
    out.add <- function.DensityLinePvalue(snp.info = snp.info.f2, wind.n = windows.number)
    lo.add <- loess(out.add$pvalue ~ out.add$window, span=smooth.par)
    xl.add <- seq(min(out.add$window), max(out.add$window), (max(out.add$window) - min(out.add$window))/100)
    pred.add <- predict(lo.add, xl.add)
    pred.add[which(pred.add < 0)] <- 0
    #add limits -- left and right for polygon function
    xl.add <- c(xl.add[1], xl.add, xl.add[length(xl.add)])
    pred.add <- c(0, pred.add, 0)
    #lines(x = xl, y = pred, col='navy', lwd=4)
    
    
    #main plot
    plot(x = 0, y = 0, xlab='Chromosomal position (Mb)', cex.lab=1.5, xaxt='none', bty='n',
         ylab="-log10(P-value)", ylim=c(min.y*y.lim/12, y.lim), cex.axis = 1.25, xaxs="i", yaxt='none',
         pch=1, col="white", cex=1.75, type = "h", lwd=2, xlim=c(min(xl), max(xl)), main=title, cex.main=2.5)
    
    #add grid
    for (x in seq(0, y.lim, (y.lim-min.y*y.lim/10)/10)){abline(h=x, lwd=0.4, col="grey80")}
    for (x in seq(min(snp.info$pos), max(snp.info$pos), (max(snp.info$pos)-min(snp.info$pos))/10)){
      segments(x0 = x, y0 = 0, x1 = x, y1 = y.lim, col = "grey80", lwd=0.4)
    }
    
    #add significance lines and corresponding legend -- for now two at 0.05 and genome-wide 5e-8
    abline(h=-log10(0.05), lty=2, col=alpha("darkgreen", 1))
    abline(h=-log10(5e-8), lty=2, col=alpha("purple", 1))
    legend("topleft", bty='n', legend = c("p=0.05", "p=5e-8"), lty=c(2,2), lwd=c(2,2), col=c("darkgreen", "purple"), cex=1.5, ncol = 2)
    
    polygon(x = xl, y = pred, col = alpha(colorPoint, 0.6), lwd=2, xaxs="i")
    polygon(x = xl.add, y = pred.add, col = alpha(colorPoint2, 0.6), lwd=2, xaxs="i")
    
    #manage axes
    axes <- seq(min(snp.info$pos), max(snp.info$pos), (max(snp.info$pos) - min(snp.info$pos))/7)
    axes.labels <- round(axes/1000000, 3)
    axis(side = 1, at=axes, cex.axis=1.5, labels=axes.labels)
    axes.x <- ceiling(seq(0, y.lim, y.lim/5))
    axis(side = 2, at = axes.x, labels = axes.x, cex.axis=1.5)
    
    
    #if input type was a single snp (either position or rs id), then color the searched variant differently -- here is a bar
    snp.interest.flag <- 0
    if (type %in% c("Position", "Rs ID")){
      #restrict to snp of interest, then plot it
      if (int.locus != "Type position..."){
        snp.interest.flag <- 1
        snp.interest <- snp.info[which(snp.info$locus == int.locus),]
        snp.interest.add <- snp.info.f2[which(snp.info.f2$locus == int.locus),]
        
        #need to grep the height in order to plot the variant -- the height is derived from the prediction
        chr.pos <- str_split_fixed(int.locus, ":", 2)
        pos.only <- as.numeric(chr.pos[, 2])
        pred.pos <- predict(lo, pos.only)
        pred.pos.add <- predict(lo.add, pos.only)
        
        #add to plot
        points(x = snp.interest.add$pos, y = pred.pos.add, type="h", col=alpha("pink", 0.8), lwd=4, xaxs="i")
        points(x = snp.interest$pos, y = pred.pos, type="h", col=alpha("yellow", 0.8), lwd=4, xaxs="i")
        
      }
    }
    
    #manage gene names
    for (g in 1:nrow(genes)){
      #main gene line -- full transcription sequence
      segments(x0 = genes$txStart[g], y0 = genes$y[g]*y.lim/12, x1 = genes$txEnd[g], y1 = genes$y[g]*y.lim/12, lwd=3)
      #need to divide exones from introns
      start <- as.data.frame(t(str_split_fixed(genes$exonStarts[g], ",", genes$exonCount[g]+1)))
      end <- as.data.frame(t(str_split_fixed(genes$exonEnds[g], ",", genes$exonCount[g]+1)))
      exons <- cbind(start, end)
      colnames(exons) <- c("start", "end")
      exons$start <- as.numeric(as.character(exons$start))
      exons$end <- as.numeric(as.character(exons$end))
      #main loop over exons
      for (j in 1:nrow(exons)){
        rect(xleft=exons$start[j], ybottom=genes$y[g]*y.lim/12-(y.lim/4/4*0.15), xright=exons$end[j], 
             ytop = genes$y[g]*y.lim/12+(y.lim/4/4*0.15), col='grey80')
      }
      text(x = genes$txStart[g] + (genes$txEnd[g] - genes$txStart[g])/2, y = genes$y[g]*y.lim/12 + (y.lim/4/4*0.5), 
           labels=genes$"#geneName"[g], font=3, cex=1)
      #     if (genes$strand[g] == "+"){
      #        arrows(x0 = genes$txEnd[g], y0 = genes$y[g], x1 = genes$txEnd[g] + (10*2/100), y1 = genes$y[g], 
      #               length=0.1, lwd=2, col='coral')
      #      } else if (genes$strand[g] == "-"){
      #        arrows(x0 = genes$txStart[g], y0 = genes$y[g], x1 = genes$txStart[g] - (10*2/100), y1 = genes$y[g], 
      #               length=0.1, lwd=2, col='coral')
      #      }
    }

    #add legend now
    if (snp.interest.flag == 0){
      legend("topright", legend = c(lab1, lab2), col = c(colorPoint, colorPoint2), pch=16, 
             cex=1.50, ncol=2, xpd=T, bty='n')
    } else {
      legend("topright", legend = c(lab1, paste(lab1, " - Input", sep=""), lab2, paste(lab2, " - Input", sep="")), 
             col = c(colorPoint, "lightblue", colorPoint2, "yellow"), pch=c(16, 23, 16, 23), pt.lwd = 1.5, 
             cex=1.50, ncol=2, xpd=T, bty='n')
    }
  }
}

#function to check if data is from the required chromosome
function.catchChromosome <- function(chr, gwas){
  if (gwas == "example"){
    #take path
    fname = paste("../data/example/chr", as.character(chr), "_IGAP_2k19.txt.gz", sep="")
    
    #read data
    dat <- fread(fname, h=T, stringsAsFactors = F)
    colnames(dat) <- c("chr", "pos", "p")      
    dat <- dat[!which(is.na(dat$p)),]
    chrom = dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$"-log10(P-value)" <- -log10(dat$p)
    
  } else {
    #in these cases, the GWAS is usually the name of the folder
    #take path
    fname = paste("../data/", gwas, "/chr", as.character(chr), "_", gwas, ".txt", sep="")
    
    #read data
    dat <- fread(fname, h=T, stringsAsFactors = F)
    colnames(dat) <- c('chr', 'pos', 'p')
    dat$p <- as.numeric(dat$p)
    dat <- dat[!which(is.na(dat$p)),]
    chrom <- dat$chr[1]
    dat$'-log10(P-value)' <- -log10(dat$p)
  }
  return(dat)
}

#function to manage position as input type
function.InputPos <- function(dat, window, input.pos, gwas){
  if (input.pos == "Type position..."){
    #initial position is at random pos
    pos.init <- dat[ceiling(nrow(dat)/4), "pos"]
    
    #take data of interest
    snp.info <- dat[which((dat$pos >= pos.init - window) & (dat$pos <= pos.init + window)),]

  } else {
    #split locus to find close positions and define intervals where to find snps
    locus <- str_split_fixed(input.pos, ":", 2)
    chr <- as.numeric(locus[, 1])
    pos.max <- as.numeric(locus[, 2]) + window
    pos.min <- as.numeric(locus[, 2]) - window

    #check if data is relative to that chromosome, and in case change it
    if (!(gwas %in% c("loaded", "toLoad"))){ dat <- function.catchChromosome(chr, gwas) }
    
    #take data of interest
    snp.info <- dat[which(dat$pos >= pos.min & dat$pos <= pos.max), ]
  }
  return(list(snp.info, dat))
}

#function to manage manual scroll as input type
function.InputManualScroll <- function(dat, window, input.scroll, input.chrom, gwas){
  print(gwas)
  #initial position is all variable
  pos.init <- as.numeric(input.scroll)
  
  #define limits
  lower <- pos.init - window
  upper <- pos.init + window
  
  #check if data is relative to that chromosome, and in case change it
  if (!(gwas %in% c("loaded", "toLoad"))){ dat <- function.catchChromosome(input.chrom, gwas) }
  
  #take data of interest
  snp.info <- dat[which(dat$pos >= lower & dat$pos <= upper), ]
  return(snp.info)
}

#function to manage genes as input
function.InputGenes <- function(gene){
  #read gene informations -- first command (commented) was to load directly from server -- now I copied the file into app folder
  gene.db <- fread("../data/databases/hg19_geneListandPos.txt.gz", h=T)
  
  if (gene != "Type gene symbol..."){
    #make sure the gene input is uppercase
    gene <- toupper(gene)
    gene <- paste("^", gene, "$", sep="")
    
    #identify gene coordinates
    gene.info <- gene.db[grep(gene, gene.db$"#geneName"),]
    gene.info <- gene.info[order(-gene.info$txEnd),]
    gene.info <- gene.info[!duplicated(gene.info$"#geneName"),]
    
  } else {
    #if there is no input gene, then take a random integer (now in chr19) and plot a random gene
    gene.db <- gene.db[which(gene.db$chrom == "chr21"),]

    gene.info <- gene.db[ceiling(nrow(gene.db)/4), ]
    
  }
  
  return(gene.info)
}

#function to manage rsID as input
function.rsIDasInput <- function(rsID){
  #check whether there is an input snp, otherwise take a random one -- now is APOE e2
  if (rsID == "Type variant identifier..."){ rsID <- "rs17125944" }
    #get snp information -- chromosome and position basically
    rsid <- rsID
    # Extract the html from the file
    cmd <- paste("wget https://www.ncbi.nlm.nih.gov/snp/", rsid, " .", sep="")
    system(cmd)
    html = read_html(rsid)
    
    # Get all the 'p' nodes (you can do the same for 'table')
    p_nodes <- html %>% html_nodes('tr')
    
    # Get the text from each node
    p_nodes_text <- p_nodes %>% html_text()
    
    # Find the nodes that have the term you are looking for
    match_indeces <- str_detect(p_nodes_text, fixed('GRCh37', ignore_case = TRUE))
    
    # Keep only the nodes with matches
    # Notice that I remove the first match because rvest adds a 
    # 'p' node to the whole file, since it is a text file
    match_p_nodes <- p_nodes[match_indeces][-1]
    
    # Extract the nodes
    x <- as.data.frame(as.character(match_p_nodes[[1]]))
    x.1 <- str_split(string = x$`as.character(match_p_nodes[[1]])`, pattern = "\n")
    x.1 <- as.data.frame(x.1)
    
    # Extract chromosome
    x.2 <- as.character(x.1[grep(pattern = "chr", x = x.1[, 1]), ])
    x.3 <- as.data.frame(str_split(string = x.2, pattern = " "))
    tmp.chr <- as.character(x.3[nrow(x.3), 1])
    chrom <- as.character(str_replace(string = tmp.chr, pattern = "</td>", replacement = ""))
    
    # Extract position
    x.4 <- as.character(x.1[grep(pattern = "&gt;", x = x.1[, 1]), ])
    x.5 <- as.data.frame(str_split(string = x.4, pattern = "<|>|g."))
    tmp.pos1 <- as.character(x.5[4, ])
    tmp.pos2 <- as.data.frame(strsplit(tmp.pos1, ""))
    pos <- paste0(tmp.pos2[1:(nrow(tmp.pos2)-2), 1], collapse = "")
    a1 <- as.character(tmp.pos2[nrow(tmp.pos2)-1, 1])
    tmp.a <- as.character(x.5[5, ])
    tmp.a2 <- as.data.frame(strsplit(tmp.a, ""))
    a2 <- as.character(tmp.a2[2, 1])
    
    # Remove temporary file
    cmd = paste("rm ", rsid, sep="")
    system(cmd)
    locus <- paste(chrom, pos, sep=":")

  return(locus)
}

#function to find gwas hits from inputted gene
function.GWASfromGene <- function(dat, gene.info, window, gwas){
  #gene information
  chr <- gene.info$chrom
  chr.n <- str_split_fixed(chr, "chr", 2)[, 2]
  start <- gene.info$txStart - window
  end <- gene.info$txEnd + window
  
  #check if data is relative to that chromosome, and in case change it
  if (!(gwas %in% c("loaded", "toLoad"))){ dat <- function.catchChromosome(chr.n, gwas) }
  
  #define snps of interest
  snp.info <- dat[which(dat$pos >= start & dat$pos <= end),]

  return(snp.info)
}

#function to find density line of p-value across chromosomal position
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
  
  #loop using sliding window
  for (i in seq(min.wind, max.wind, interval)){
    #define maximum -- f
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
  #remove NA and substitute them with 0
  out[which(is.na(out[, "error"]) == TRUE), "error"] <- 0
  out <- as.data.frame(out)
  
  return(out)
}

#function to parse gene location file and grep only region of interest
function.dynamicGene <- function(snp.info){
  #define searching space for genes
  min.x <- min(snp.info$pos)
  max.x <- max(snp.info$pos)
  
  #input for genes
  gene.db <- fread("../data/databases/hg19_geneListandPos.txt.gz", h=T)
  
  #find genes in interval (min.x -- max.x) -- then clean up a bit
  gene.loc.res <- subset(gene.db, gene.db$chrom == paste("chr", snp.info$chr[1], sep=""))
  gene.loc.res$in.int <- 0
  gene.loc.res$in.int[which((gene.loc.res$txStart >= min.x) & (gene.loc.res$txEnd <= max.x))] <- 1
  gene.loc.res$in.int[which((gene.loc.res$txStart <= min.x) & (gene.loc.res$txEnd >= min.x))] <- 1
  gene.loc.res$in.int[which((gene.loc.res$txStart <= max.x) & (gene.loc.res$txEnd >= max.x))] <- 1
  genes <- gene.loc.res[which(gene.loc.res$in.int == 1),]
  genes <- genes[order(-genes$txEnd),]
  genes <- genes[!duplicated(genes$"#geneName"),]
  
  #define position in the plot
  n = ceiling(nrow(genes)/4)
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

#function to manage which input data to plot -- still 1 dataset only for now
function.manageInput <- function(inp, supp_f){
  if (length(inp) == 0){
    dat <- fread("../data/example/chr21_IGAP_2k19.txt.gz", h=T, stringsAsFactors = F)
    colnames(dat) <- c("chr", "pos", "p")      
    dat <- dat[!which(is.na(dat$p)),]
    chrom = dat$chr[1]
    dat$"-log10(P-value)" <- -log10(as.numeric(dat$p))
    gwas = "example"
    
  } else if ("toLoad" %in% inp){
    #read file
    dat <- fread(supp_f, h=T, stringsAsFactors = F)
    
    #analyse header and slim data
    dat = identiHeader(dat)
    dat$p <- as.numeric(dat$p)
    dat <- dat[!which(is.na(dat$p)),]
    chrom <- dat$chr[1]
    dat$"-log10(P-value)" <- -log10(as.numeric(dat$p))
    gwas = "loaded"
    
	} else if (inp == 'IGAP') {
		dat <- fread('../data/IGAP/chr19_IGAP.txt', h=T, stringsAsFactors=F)
		colnames(dat) <- c('chr', 'pos', 'p')
		dat <- dat[!which(is.na(dat$p)),]
		chrom <- dat$chr[1]
		dat$p <- as.numeric(dat$p)
		dat$'-log10(P-value)' <- -log10(dat$p)
		gwas <- 'IGAP'
	} else if (inp == 'CARDIO') {
		dat <- fread('../data/CARDIO/chr19_CARDIO.txt', h=T, stringsAsFactors=F)
		colnames(dat) <- c('chr', 'pos', 'p')
		dat <- dat[!which(is.na(dat$p)),]
		chrom <- dat$chr[1]
		dat$p <- as.numeric(dat$p)
		dat$'-log10(P-value)' <- -log10(dat$p)
		gwas <- 'CARDIO'
	} else if (inp == '100plus') {
		dat <- fread('../data/100plus/chr19_100plus.txt', h=T, stringsAsFactors=F)
		colnames(dat) <- c('chr', 'pos', 'p')
		dat <- dat[!which(is.na(dat$p)),]
		chrom <- dat$chr[1]
		dat$p <- as.numeric(dat$p)
		dat$'-log10(P-value)' <- -log10(dat$p)
		gwas <- '100plus'
	} else if (inp == 'MetaGrace') {
		dat <- fread('../data/MetaGrace/chr19_MetaGrace.txt', h=T, stringsAsFactors=F)
		colnames(dat) <- c('chr', 'pos', 'p')
		dat <- dat[!which(is.na(dat$p)),]
		chrom <- dat$chr[1]
		dat$p <- as.numeric(dat$p)
		dat$'-log10(P-value)' <- -log10(dat$p)
		gwas <- 'MetaGrace'
	} else if (inp == 'UKBaging') {
		dat <- fread('../data/UKBaging/chr19_UKBaging.txt', h=T, stringsAsFactors=F)
		colnames(dat) <- c('chr', 'pos', 'p')
		dat <- dat[!which(is.na(dat$p)),]
		chrom <- dat$chr[1]
		dat$p <- as.numeric(dat$p)
		dat$'-log10(P-value)' <- -log10(dat$p)
		gwas <- 'UKBaging'
	} else if (inp == 'exomeADES') {
		dat <- fread('../data/exomeADES/chr19_exomeADES.txt', h=T, stringsAsFactors=F)
		colnames(dat) <- c('chr', 'pos', 'p')
		dat <- dat[!which(is.na(dat$p)),]
		chrom <- dat$chr[1]
		dat$p <- as.numeric(dat$p)
		dat$'-log10(P-value)' <- -log10(dat$p)
		gwas <- 'exomeADES'
  } #else here to add repositories


  return(list(dat, chrom, gwas))
  
}

#function to identify header in uploaded file
identiHeader <- function(dat){
  head <- names(dat)
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
    print("## header correctly read")
    
    #now rename columns
    head[chrom.idx] <- "chr"
    head[pos.idx] <- "pos"
    head[p.idx] <- "p"
    colnames(dat) <- head
    
  } else {
      print("!!!There was a problem with the input!!!!")
  }
  return(dat)
}

#MAIN APP
shinyServer(
  function(input, output) {
    #what to plot here
    output$hist <- renderPlot({
      
      ############################
      #create flag for uploaded input file
      inFile <- input$inp_f
      path_f <- "None"
      if (!is.null(inFile)){
        path_f = inFile$datapath
        print("file loaded")
      } else {
        path_f <- "None"
      }
      ##########################

      ####################
      #manage which input to plot in two cases: beginning of app and 1 gwas as input
      res = function.manageInput(as.character(input$gwas[1]), path_f)
      dat = as.data.frame(res[[1]])
      chrom = as.numeric(res[[2]])
      gwas = as.character(res[[3]])
      n.inputs = 1
      ####################
      
      ####################
      #Browsing options -- default is Position -- by default take middle of chromosome
      if (input$sel %in% c("Position", "RsID")){
        
        #check number of GWAS to be plotted
        if (length(input$gwas) < 2){
          
          #check if rsid is the input
          if (input$sel == "RsID"){ 
            locus <- function.rsIDasInput(input$snpID)
            res <- function.InputPos(dat = dat, window = input$x, input.pos = locus, gwas)
          } else {
            res <- function.InputPos(dat = dat, window = input$x, input.pos = input$pos, gwas)
            }
          snp.info <- as.data.frame(res[[1]])
          dat <- as.data.frame(res[[2]])

          #plot
          function.plot(snp.info = snp.info, y.lim = input$y, type = input$sel, plt.type = input$ploType, 
                        windows.number = input$sliding.window, smooth.par = input$smooth, int.locus = input$pos, gwas = gwas, col = input$col)
          ###################
      
          #the following in case multiple gwas need to be plotted
        } else {
          #check whether uploaded file is among those to show
          if (path_f == "None"){
            #first dataset
            d1 <- input$gwas[1]
            
            #check if rsid or position
            if (input$sel == "RsID"){
              locus <- function.rsIDasInput(input$snpID)
              res <- function.InputPos(dat = dat, window = input$x, input.pos = locus, d1)
            } else {
              res <- function.InputPos(dat = dat, window = input$x, input.pos = input$pos, d1)
            }
            snp.info <- as.data.frame(res[[1]])
            dat <- as.data.frame(res[[2]])

            #second dataset
            d2 <- input$gwas[2]
            res = function.manageInput(d2, path_f)
            dat.2 = as.data.frame(res[[1]])
            chrom.2 = as.numeric(res[[2]])
            
            #check if rsid or position
            if (input$sel == "RsID"){
              locus <- function.rsIDasInput(input$snpID)
              res.2 <- function.InputPos(dat = dat.2, window = input$x, input.pos = locus, d2)
            } else {
              res.2 <- function.InputPos(dat = dat.2, window = input$x, input.pos = input$pos, gwas = d2)
            }
            snp.info.2 <- res.2[[1]]
            dat.2 <- res.2[[2]]

            #plot
            function.multiPlot(snp.info = snp.info, snp.info.f2 = snp.info.2, y.lim = input$y, type = input$sel, 
                               plt.type = input$ploType, windows.number = input$sliding.window, smooth.par = input$smooth, 
                               int.locus = input$pos, d1, d2, colorPoint = input$col, colorPoint2 = input$col2)
          } else {
            #first dataset -- this is the loaded one
            index_upl <- grep("toLoad", input$gwas)
            d1 <- input$gwas[index_upl]
            res = function.manageInput(d1, path_f)
            dat = as.data.frame(res[[1]])
            chrom = as.numeric(res[[2]])
            gwas = as.character(res[[3]])
            res <- function.InputPos(dat = dat, window = input$x, input.pos = input$pos, gwas)
            snp.info <- as.data.frame(res[[1]])
            dat <- as.data.frame(res[[2]])

            #plot
            function.multiPlot(snp.info = snp.info, snp.info.f2 = snp.info.2, y.lim = input$y, type = input$sel, 
                               plt.type = input$ploType, windows.number = input$sliding.window, smooth.par = input$smooth, 
                               int.locus = input$pos, gwas, gwas.2, colorPoint = input$col, colorPoint2 = input$col2)
            #################################
          } 
        }
        ####################
        
        ####################
        #this in case you want manual scroll -- not implemented yet
      } else if (input$sel == "Manual scroll"){
        if (length(input$gwas) < 2){
          #manage different manual scroll input
          snp.info <- function.InputManualScroll(dat = dat, window = input$x, input.scroll = input$manual.pos, input.chrom = input$manual.chrom, gwas)
          
          #plot
          function.plot(snp.info = snp.info, y.lim = input$y, type = input$sel, plt.type = input$ploType, 
                        windows.number = input$sliding.window, smooth.par = input$smooth, int.locus = input$manual.pos, gwas = gwas, col = input$col)
          #####################
          
          #the followinf in case multiple gwas need to be plotted
        } else {
          #check whether uploaded file is among those to show
          if (path_f == "None"){
            #first dataset
            d1 <- input$gwas[1]
            snp.info <- function.InputManualScroll(dat = dat, window = input$x, input.scroll = input$manual.pos, input.chrom = input$manual.chrom, d1)

            #second dataset
            d2 <- input$gwas[2]
            res = function.manageInput(d2, path_f)
            dat.2 = as.data.frame(res[[1]])
            chrom.2 = as.numeric(res[[2]])
            snp.info.2 <- function.InputManualScroll(dat = dat.2, window = input$x, input.scroll = input$manual.pos, input.chrom = input$manual.chrom, d2)
            
            #plot
            function.multiPlot(snp.info = snp.info, snp.info.f2 = snp.info.2, y.lim = input$y, type = input$sel, 
                               plt.type = input$ploType, windows.number = input$sliding.window, smooth.par = input$smooth, 
                               int.locus = input$pos, d1, d2, colorPoint = input$col, colorPoint2 = input$col2)
          } else {
            #first dataset -- this is the loaded one
            index_upl <- grep("toLoad", input$gwas)
            d1 <- input$gwas[index_upl]
            res = function.manageInput(d1, path_f)
            dat = as.data.frame(res[[1]])
            chrom = as.numeric(res[[2]])
            gwas = as.character(res[[3]])
            snp.info <- function.InputManualScroll(dat = dat, window = input$x, input.scroll = input$manual.pos, input.chrom = input$manual.chrom, d1)
            
            #second dataset
            index_oth <- which(input$gwas != "toLoad")
            d2 <- input$gwas[index_oth]
            res.2 = function.manageInput(d2, "None")
            dat.2 = as.data.frame(res.2[[1]])
            chrom.2 = as.numeric(res.2[[2]])
            gwas.2 = as.character(res.2[[3]])
            snp.info.2 <- function.InputManualScroll(dat = dat.2, window = input$x, input.scroll = input$manual.pos, input.chrom = input$manual.chrom, d2)
            
            #plot
            function.multiPlot(snp.info = snp.info, snp.info.f2 = snp.info.2, y.lim = input$y, type = input$sel, 
                               plt.type = input$ploType, windows.number = input$sliding.window, smooth.par = input$smooth, 
                               int.locus = input$pos, gwas, gwas.2, colorPoint = input$col, colorPoint2 = input$col2)
            #######################################
          }
        }
        ####################
        
        ####################
        #this is for genes as input
      } else if (input$sel == "Gene"){
        #this in case you only plot 1 gwas
        if (length(input$gwas) < 2){
          
          #find position of input gene
          gene.info <- function.InputGenes(gene = as.character(input$gene))
          
          #find gwas info around the gene
          snp.info <- function.GWASfromGene(dat = dat, gene.info = gene.info, window = input$x, gwas)
          
          #plot
          function.plot(snp.info = snp.info, y.lim = input$y, type = input$sel, plt.type = input$ploType, 
                        windows.number = input$sliding.window, smooth.par = input$smooth, int.locus = input$gene, gwas = gwas, col = input$col)
          ##################
          #this in case you want to plot multiple gwas
        } else {
          #check whether uploaded file is among those to show
          if (path_f == "None"){
            #find position of input gene
            gene.info <- function.InputGenes(gene = as.character(input$gene))

            #first dataset
            d1 <- input$gwas[1]
            #find gwas info around the gene for gwas 1
            snp.info <- function.GWASfromGene(dat = dat, gene.info = gene.info, window = input$x, d1)
            
            #second dataset
            d2 <- input$gwas[2]
            res = function.manageInput(d2, path_f)
            dat.2 = as.data.frame(res[[1]])
            chrom.2 = as.numeric(res[[2]])

            #additional dataset
            snp.info.2 <- function.GWASfromGene(dat = dat.2, gene.info = gene.info, window = input$x, d2)

            #plot
            function.multiPlot(snp.info = snp.info, snp.info.f2 = snp.info.2, y.lim = input$y, type = input$sel, 
                               plt.type = input$ploType, windows.number = input$sliding.window, smooth.par = input$smooth, 
                               int.locus = input$pos, d1, d2, colorPoint = input$col, colorPoint2 = input$col2)
          } else {
            #find position of input gene
            gene.info <- function.InputGenes(gene = as.character(input$gene))
            
            #first dataset -- this is the loaded one
            index_upl <- grep("toLoad", input$gwas)
            d1 <- input$gwas[index_upl]
            res = function.manageInput(d1, path_f)
            dat = as.data.frame(res[[1]])
            chrom = as.numeric(res[[2]])
            gwas = as.character(res[[3]])
            snp.info <- function.GWASfromGene(dat = dat, gene.info = gene.info, window = input$x, d1)
            
            #second dataset
            index_oth <- which(input$gwas != "toLoad")
            d2 <- input$gwas[index_oth]
            res.2 = function.manageInput(d2, "None")
            dat.2 = as.data.frame(res.2[[1]])
            chrom.2 = as.numeric(res.2[[2]])
            gwas.2 = as.character(res.2[[3]])
            snp.info.2 <- function.GWASfromGene(dat = dat.2, gene.info = gene.info, window = input$x, d2)
            
            #plot
            function.multiPlot(snp.info = snp.info, snp.info.f2 = snp.info.2, y.lim = input$y, type = input$sel, 
                               plt.type = input$ploType, windows.number = input$sliding.window, smooth.par = input$smooth, 
                               int.locus = input$pos, gwas, gwas.2, colorPoint = input$col, colorPoint2 = input$col2)
            #######################################
          }
        }
        ####################
        
        ####################
      } 
    })
    
    #this should be to download image
    output$plot <- reactivePlot(function() {
      name <- paste0(input$filename, ".png")
      if(input$savePlot) {
        ggsave(name, plotInput(), type="cairo-png")
      }
      else print(plotInput())
    })
    
    #manage click on points
    output$click_info <- renderPrint({
      
      ############################
      #create flag for uploaded input file
      inFile <- input$inp_f
      path_f <- "None"
      if (!is.null(inFile)){
        path_f = inFile$datapath
        print("file loaded")
      } else {
        path_f <- "None"
      }
      ##########################
      
      ####################
      #manage which input to plot in two cases: beginning of app and 1 gwas as input
      res = function.manageInput(as.character(input$gwas), path_f)
      dat = as.data.frame(res[[1]])
      chrom = as.numeric(res[[2]])
      gwas = as.character(res[[3]])
      dat$"-log10(P-value)" <- -log10(as.numeric(dat$p))
      ####################
      
      if (input$sel == "Position"){
        #check if only 1 gwas is selected
        if (length(input$gwas) < 2){
          #manage different option of position as input type
          res <- function.InputPos(dat = dat, window = input$x, input.pos = input$pos, gwas)
          snp.info <- as.data.frame(res[[1]])
          dat <- as.data.frame(res[[2]])
          snp.info$"-log10(P-value)" <- -log10(snp.info$p)
          if (length(input$plot1_click) == 0 || length(input$plot1_brush) > 0){
            print("Click on SNPs to get info")
          } else {
            nearPoints(df = snp.info, input$plot1_click, addDist = FALSE, xvar = "pos", yvar = "-log10(P-value)")
          }    
          ###################
          
          ###################
          #this in case multiple gwas need to be plotted -- to be implemented
        } else {
          #first dataset
          # d1 <- input$gwas[1]
          # res <- function.InputPos(dat = dat, window = input$x, input.pos = input$pos, d1)
          # snp.info <- as.data.frame(res[[1]])
          # dat <- as.data.frame(res[[2]])
          # 
          # #second dataset
          # d2 <- input$gwas[2]
          # res = function.manageInput(d2)
          # dat.2 = as.data.frame(res[[1]])
          # chrom.2 = as.numeric(res[[2]])
          # 
          # #this was for the input file from user
          # #snp.info.f2 <- function.InputPosGeneric(dat = dat.2, window = input$x, input.pos = input$pos)
          # 
          # res.2 <- function.InputPos(dat = dat.2, window = input$x, input.pos = input$pos, gwas = d2)
          # snp.info.2 <- res.2[[1]]
          # dat.2 <- res.2[[2]]
          # 
          # #plot
          # function.multiPlot(snp.info = snp.info, snp.info.f2 = snp.info.2, y.lim = input$y, type = input$sel, 
          #                    plt.type = input$ploType, windows.number = input$sliding.window, smooth.par = input$smooth, int.locus = input$pos,
          #                    d1, d2)
        }
        ####################
        
        ####################
        #this in case you want manual scroll -- not implemented yet
      } else if (input$sel == "Gene"){
        #this in case you only plot 1 gwas
        if (length(input$gwas) < 2){
          #find position of input gene
          gene.info <- function.InputGenes(gene = as.character(input$gene))
          
          #find gwas info around the gene
          snp.info <- function.GWASfromGene(dat = dat, gene.info = gene.info, window = input$x, gwas)
          
          snp.info$"-log10(P-value)" <- -log10(as.numeric(snp.info$p))
          
          if (length(input$plot1_click) == 0 || length(input$plot1_brush) > 0){
            print("Click on SNPs to get info")
          } else {
            nearPoints(df = snp.info, input$plot1_click, addDist = FALSE, xvar = "pos", yvar = "-log10(P-value)")
          }    
          ##################
          
          ################## 
          #this in case you want to plot multiple gwas -- to be implemented
        } else {
          #find position of input gene
          # gene.info <- function.InputGenes(gene = as.character(input$gene))
          # 
          # #define the two gwas
          # d1 <- input$gwas[1]
          # d2 <- input$gwas[2]
          # 
          # #find gwas info around the gene for gwas 1
          # snp.info <- function.GWASfromGene(dat = dat, gene.info = gene.info, window = input$x, d1)
          # 
          # #second dataset
          # res = function.manageInput(d2)
          # dat.2 = as.data.frame(res[[1]])
          # chrom.2 = as.numeric(res[[2]])
          # 
          # #additional dataset
          # snp.info.2 <- function.GWASfromGene(dat = dat.2, gene.info = gene.info, window = input$x, d2)
          # 
          # #plot
          # function.multiPlot(snp.info = snp.info, snp.info.f2 = snp.info.2, y.lim = input$y, type = input$sel, 
          #                    plt.type = input$ploType, windows.number = input$sliding.window, smooth.par = input$smooth, int.locus = input$pos,
          #                    d1, d2)
        }
        ####################
      }
    })
    
    #manage brush on points
    output$brush_info <- renderPrint({
      
      ############################
      #create flag for uploaded input file
      inFile <- input$inp_f
      path_f <- "None"
      if (!is.null(inFile)){
        path_f = inFile$datapath
        print("file loaded")
      } else {
        path_f <- "None"
      }
      ##########################
      
      ####################
      #manage which input to plot in two cases: beginning of app and 1 gwas as input
      res = function.manageInput(as.character(input$gwas), path_f)
      dat = as.data.frame(res[[1]])
      chrom = as.numeric(res[[2]])
      gwas = as.character(res[[3]])
      dat$"-log10(P-value)" <- -log10(as.numeric(dat$p))
      ####################
      
      if (input$sel == "Position"){
        #check if only 1 gwas is selected
        if (length(input$gwas) < 2){
          #manage different option of position as input type
          res <- function.InputPos(dat = dat, window = input$x, input.pos = input$pos, gwas)
          snp.info <- as.data.frame(res[[1]])
          dat <- as.data.frame(res[[2]])
          
          snp.info$"-log10(P-value)" <- -log10(snp.info$p)
          if (length(input$plot1_brush) == 0){
            print("Brush on SNPs to get info")
          } else {
            brushedPoints(snp.info, input$plot1_brush, xvar = "pos", yvar = "-log10(P-value)")
          }    
          ###################
          
          ###################
          #this in case multiple gwas need to be plotted -- to be implemented
        } else {
          #first dataset
          # d1 <- input$gwas[1]
          # res <- function.InputPos(dat = dat, window = input$x, input.pos = input$pos, d1)
          # snp.info <- as.data.frame(res[[1]])
          # dat <- as.data.frame(res[[2]])
          # 
          # #second dataset
          # d2 <- input$gwas[2]
          # res = function.manageInput(d2)
          # dat.2 = as.data.frame(res[[1]])
          # chrom.2 = as.numeric(res[[2]])
          # 
          # #this was for the input file from user
          # #snp.info.f2 <- function.InputPosGeneric(dat = dat.2, window = input$x, input.pos = input$pos)
          # 
          # res.2 <- function.InputPos(dat = dat.2, window = input$x, input.pos = input$pos, gwas = d2)
          # snp.info.2 <- res.2[[1]]
          # dat.2 <- res.2[[2]]
          # 
          # #plot
          # function.multiPlot(snp.info = snp.info, snp.info.f2 = snp.info.2, y.lim = input$y, type = input$sel, 
          #                    plt.type = input$ploType, windows.number = input$sliding.window, smooth.par = input$smooth, int.locus = input$pos,
          #                    d1, d2)
        }
        ####################
        
        ####################
        #this in case you want manual scroll -- not implemented yet
      } else if (input$sel == "Gene"){
        #this in case you only plot 1 gwas
        if (length(input$gwas) < 2){
          #find position of input gene
          gene.info <- function.InputGenes(gene = as.character(input$gene))
          
          #find gwas info around the gene
          snp.info <- function.GWASfromGene(dat = dat, gene.info = gene.info, window = input$x, gwas)
          
          snp.info$"-log10(P-value)" <- -log10(as.numeric(snp.info$p))
          
          if (length(input$plot1_brush) == 0){
            print("Click on SNPs to get info")
          } else {
            brushedPoints(snp.info, input$plot1_brush, xvar = "pos", yvar = "-log10(P-value)")
          }    
          ##################
          
          ################## 
          #this in case you want to plot multiple gwas -- to be implemented
        } else {
          #find position of input gene
          # gene.info <- function.InputGenes(gene = as.character(input$gene))
          # 
          # #define the two gwas
          # d1 <- input$gwas[1]
          # d2 <- input$gwas[2]
          # 
          # #find gwas info around the gene for gwas 1
          # snp.info <- function.GWASfromGene(dat = dat, gene.info = gene.info, window = input$x, d1)
          # 
          # #second dataset
          # res = function.manageInput(d2)
          # dat.2 = as.data.frame(res[[1]])
          # chrom.2 = as.numeric(res[[2]])
          # 
          # #additional dataset
          # snp.info.2 <- function.GWASfromGene(dat = dat.2, gene.info = gene.info, window = input$x, d2)
          # 
          # #plot
          # function.multiPlot(snp.info = snp.info, snp.info.f2 = snp.info.2, y.lim = input$y, type = input$sel, 
          #                    plt.type = input$ploType, windows.number = input$sliding.window, smooth.par = input$smooth, int.locus = input$pos,
          #                    d1, d2)
        }
        ####################
      }
    })
  })





