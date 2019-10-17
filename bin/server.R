#SERVER SNP BROWSER

#LIBRARIES
library(shiny)
library(data.table)
library(stringr)
library(ggplot2)

#FUNCTIONS
#basic function to plot stuff 
function.plot <- function(snp.info, y.lim, type, plt.type, windows.number, smooth.par, int.locus, gwas){
  par(mar=c(5, 5, 4, 1))
  
  #prepare title
  t <- function.title(gwas)
  title = paste(t, " ~ chr", snp.info$chr[1], ": ", min(snp.info$pos), " - ", max(snp.info$pos), sep = "")

  #assign dot sizes (function for this)
  snp.info <- function.pointSize(dat = snp.info, range = seq(2, 7, 0.5))

  #independently from plot type, I need the genes that are in the window to adjust y-axis -- here it is
  genes <- function.dynamicGene(snp.info = snp.info)
  if (nrow(genes) > 0){
    min.y <- min(genes$y)
  } else {
    min.y <- 0
  }
  
  if (plt.type == "Points"){
    #empty plot
    plot(x = 0, y = 0, xlab='Chromosomal position (Mb)', cex.lab=1.75, xaxt='none',
         ylab="-log10(P-value)", ylim=c(min.y*y.lim/10, y.lim), cex.axis = 1.5,
         pch=16, col="white", cex=2, type = "p", xaxs="i", yaxt='none', xlim=c(min(snp.info$pos), max(snp.info$pos)), main=title, cex.main=2.50, bty='n')
    
    #add grid
    for (x in seq(0, y.lim, (y.lim-min.y*y.lim/10)/10)){abline(h=x, lwd=0.4, col="grey80")}
    for (x in seq(min(snp.info$pos), max(snp.info$pos), (max(snp.info$pos)-min(snp.info$pos))/10)){
      segments(x0 = x, y0 = 0, x1 = x, y1 = y.lim, col = "grey80", lwd=0.4)
    }
    
    #add significance lines and corresponding legend -- for now two at 0.05 and genome-wide 5e-8
    abline(h=-log10(0.05), lty=2, col=alpha("darkgreen", 1))
    abline(h=-log10(5e-8), lty=2, col=alpha("purple", 1))
    legend("topleft", bty='n', legend = c("p=0.05", "p=5e-8"), lty=c(2,2), lwd=c(2,2), col=c("darkgreen", "purple"), cex=1.5, ncol = 2)
    
    #then points
    snp.info$"-log10(P-value)" <- -log10(as.numeric(snp.info$p))
    points(x = snp.info$pos, y = snp.info$"-log10(P-value)", xlab='Chromosomal position (Mb)', cex.lab=1.5, xaxt='none',
                      pch=16, col=alpha("navy", 0.6), cex=snp.info$size, type = "p", xaxs="i", yaxt='none', xlim=c(min(snp.info$pos), max(snp.info$pos)))
    
    #if input type was a single snp (either position or rs id), then color the searched variant differently
    if (type %in% c("Position", "Rs ID")){
      #restrict to snp of interest, then plot it
      if (int.locus != "Type position..."){
        snp.interest <- snp.info[which(snp.info$locus == int.locus),]
        points(x = snp.interest$pos, y = snp.interest$"-log10(P-value)", pch=16, col=alpha("orange", 0.8), cex=snp.interest$size)
      }
    }
    
    #manage axes
    axes <- seq(min(snp.info$pos), max(snp.info$pos), (max(snp.info$pos) - min(snp.info$pos))/7)
    axes.labels <- round(axes/1000000, 3)
    axis(side = 1, at=axes, cex.axis=1.5, labels=axes.labels)
    axes.x <- ceiling(seq(0, y.lim, y.lim/5))
    axis(side = 2, at = axes.x, labels = axes.x, cex.axis=1.5)
    
    if (min.y != 0){
      #manage gene names
      for (g in 1:nrow(genes)){
        #main gene line -- full transcription sequence
        segments(x0 = genes$txStart[g], y0 = genes$y[g]*y.lim/10, x1 = genes$txEnd[g], y1 = genes$y[g]*y.lim/10, lwd=3)
        #need to divide exones from introns
        start <- as.data.frame(t(str_split_fixed(genes$exonStarts[g], ",", genes$exonCount[g]+1)))
        end <- as.data.frame(t(str_split_fixed(genes$exonEnds[g], ",", genes$exonCount[g]+1)))
        exons <- cbind(start, end)
        colnames(exons) <- c("start", "end")
        exons$start <- as.numeric(as.character(exons$start))
        exons$end <- as.numeric(as.character(exons$end))
        #main loop over exons
        for (j in 1:nrow(exons)){
          rect(xleft=exons$start[j], ybottom=genes$y[g]*y.lim/10-(y.lim/4/4*0.15), xright=exons$end[j], 
               ytop = genes$y[g]*y.lim/10+(y.lim/4/4*0.15), col='grey80')
        }
        text(x = genes$txStart[g] + (genes$txEnd[g] - genes$txStart[g])/2, y = genes$y[g]*y.lim/10 + (y.lim/4/4*0.5), 
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
         ylab="-log10(P-value)", ylim=c(min.y*y.lim/10, y.lim), cex.axis = 1.25, xaxs="i", yaxt='none',
         pch=1, col="white", cex=1.75, type = "h", lwd=2, xlim=c(min(xl), max(xl)))
    
    polygon(x = xl, y = pred, col = alpha("navy", 0.6), lwd=2, xaxs="i")
    
    #manage axes
    axes <- seq(min(snp.info$pos), max(snp.info$pos), (max(snp.info$pos) - min(snp.info$pos))/7)
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
      segments(x0 = genes$txStart[g], y0 = genes$y[g]*y.lim/10, x1 = genes$txEnd[g], y1 = genes$y[g]*y.lim/10, lwd=3)
      #need to divide exones from introns
      start <- as.data.frame(t(str_split_fixed(genes$exonStarts[g], ",", genes$exonCount[g]+1)))
      end <- as.data.frame(t(str_split_fixed(genes$exonEnds[g], ",", genes$exonCount[g]+1)))
      exons <- cbind(start, end)
      colnames(exons) <- c("start", "end")
      exons$start <- as.numeric(as.character(exons$start))
      exons$end <- as.numeric(as.character(exons$end))
      #main loop over exons
      for (j in 1:nrow(exons)){
        rect(xleft=exons$start[j], ybottom=genes$y[g]*y.lim/10-(y.lim/4/4*0.15), xright=exons$end[j], 
             ytop = genes$y[g]*y.lim/10+(y.lim/4/4*0.15), col='grey80')
      }
      text(x = genes$txStart[g] + (genes$txEnd[g] - genes$txStart[g])/2, y = genes$y[g]*y.lim/10 + (y.lim/4/4*0.5), 
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

#function to assign dot size -- define range
function.pointSize <- function(dat, range){
  dat$size <- 2
  
  for (x in range){
    dat$size[which(-log10(as.numeric(dat$p)) >= x)] <- x
  }
  
  return(dat)
}

#function to manage titles
function.title <- function(gwas){
  t = ""
  if (gwas == "ad"){t = "AD vs. CTR"} else if (gwas == "age"){t = "CHC vs. CTR"} else if (gwas == "igap"){t = "IGAP"}
  
  return(t)
}

#plot in case of multiple files
function.multiPlot <- function(snp.info, snp.info.f2, y.lim, type, plt.type, windows.number, smooth.par, int.locus, lab1, lab2){
  par(mar=c(5, 5, 4, 1))
  
  #prepare title
  t1 <- function.title(lab1)
  t2 <- function.title(lab2)
  title = paste(t1, " AND ", t2, " ~ chr", snp.info$chr[1], ": ", min(snp.info$pos), " - ", max(snp.info$pos), sep = "")
  
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
         ylab="-log10(P-value)", ylim=c(min(genes$y)*y.lim/10, y.lim), cex.axis = 1.5,
         pch=16, col="white", cex=2, type = "p", xaxs="i", yaxt='none', xlim=c(min(snp.info$pos), max(snp.info$pos)), main=title, cex.main=2.50, bty='n')
    
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
           pch=16, col=alpha("navy", 0.6), cex=snp.info$size, type = "p", xaxs="i")
    
    #add additional file data points -- gwas 2
    snp.info.f2$"-log10(P-value)" <- -log10(snp.info.f2$p)
    points(x = snp.info.f2$pos, y = snp.info.f2$"-log10(P-value)", pch=16, col=alpha("red", 0.6), 
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
        
        points(x = snp.interest$pos, y = -log10(snp.interest$p), pch=16, col=alpha("lightblue", 1), cex=snp.interest$size)
        points(x = snp.interest.addF$pos, y = -log10(as.numeric(snp.interest.addF$p)), pch=16, col=alpha("orange", 1), cex=snp.interest.addF$size)
      }
    }
    
    #add legend now
    if (snp.interest.flag == 0){
      legend("topright", legend = c(t1, t2), col = c("navy","red"), pch=16, 
             cex=1.50, ncol=2, xpd=T, bty='n')
      
    } else {
      legend("topright", legend = c(t1, paste(t1, " - Input", sep=""), t2, paste(t2, " - Input", sep="")), col = c("navy", "lightblue", "red", "orange"), pch=16, 
             cex=1.50, ncol=2, xpd=T, bty='n')
    }
    
    if (min.y != 0){
      #manage gene names
      for (g in 1:nrow(genes)){
        #main gene line -- full transcription sequence
        segments(x0 = genes$txStart[g], y0 = genes$y[g]*y.lim/10, x1 = genes$txEnd[g], y1 = genes$y[g]*y.lim/10, lwd=3)
        #need to divide exones from introns
        start <- as.data.frame(t(str_split_fixed(genes$exonStarts[g], ",", genes$exonCount[g]+1)))
        end <- as.data.frame(t(str_split_fixed(genes$exonEnds[g], ",", genes$exonCount[g]+1)))
        exons <- cbind(start, end)
        colnames(exons) <- c("start", "end")
        exons$start <- as.numeric(as.character(exons$start))
        exons$end <- as.numeric(as.character(exons$end))
        #main loop over exons
        for (j in 1:nrow(exons)){
          rect(xleft=exons$start[j], ybottom=genes$y[g]*y.lim/10-(y.lim/4/4*0.15), xright=exons$end[j], 
               ytop = genes$y[g]*y.lim/10+(y.lim/4/4*0.15), col='grey80')
        }
        text(x = genes$txStart[g] + (genes$txEnd[g] - genes$txStart[g])/2, y = genes$y[g]*y.lim/10 + (y.lim/4/4*0.5), 
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
    plot(x = 0, y = 0, xlab='Chromosomal position (Mb)', cex.lab=1.5, xaxt='none',
         ylab="-log10(P-value)", ylim=c(min.y*y.lim/10, y.lim), cex.axis = 1.25, xaxs="i", yaxt='none',
         pch=1, col="white", cex=1.75, type = "h", lwd=2, xlim=c(min(xl), max(xl)), main=title, cex.main=2.5)
    
    #add significance lines and corresponding legend -- for now two at 0.05 and genome-wide 5e-8
    abline(h=-log10(0.05), lty=2, col=alpha("darkgreen", 1))
    abline(h=-log10(5e-8), lty=2, col=alpha("purple", 1))
    legend("topleft", bty='n', legend = c("p=0.05", "p=5e-8"), lty=c(2,2), lwd=c(2,2), col=c("darkgreen", "purple"), cex=1.5, ncol = 2)
    
    polygon(x = xl, y = pred, col = alpha("navy", 0.6), lwd=2, xaxs="i")
    polygon(x = xl.add, y = pred.add, col = alpha("red", 0.6), lwd=2, xaxs="i")
    
    #manage axes
    axes <- seq(min(snp.info$pos), max(snp.info$pos), (max(snp.info$pos) - min(snp.info$pos))/7)
    axes.labels <- round(axes/1000000, 2)
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
      segments(x0 = genes$txStart[g], y0 = genes$y[g]*y.lim/10, x1 = genes$txEnd[g], y1 = genes$y[g]*y.lim/10, lwd=3)
      #need to divide exones from introns
      start <- as.data.frame(t(str_split_fixed(genes$exonStarts[g], ",", genes$exonCount[g]+1)))
      end <- as.data.frame(t(str_split_fixed(genes$exonEnds[g], ",", genes$exonCount[g]+1)))
      exons <- cbind(start, end)
      colnames(exons) <- c("start", "end")
      exons$start <- as.numeric(as.character(exons$start))
      exons$end <- as.numeric(as.character(exons$end))
      #main loop over exons
      for (j in 1:nrow(exons)){
        rect(xleft=exons$start[j], ybottom=genes$y[g]*y.lim/10-(y.lim/4/4*0.15), xright=exons$end[j], 
             ytop = genes$y[g]*y.lim/10+(y.lim/4/4*0.15), col='grey80')
      }
      text(x = genes$txStart[g] + (genes$txEnd[g] - genes$txStart[g])/2, y = genes$y[g]*y.lim/10 + (y.lim/4/4*0.5), 
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
      legend("topright", legend = c(t1, t2), col = c("navy","red"), pch=16, 
             cex=1.50, ncol=2, xpd=T, bty='n')
      
    } else {
      legend("topright", legend = c(t1, paste(t1, " - Input", sep=""), t2, paste(t2, " - Input", sep="")), col = c("navy", "lightblue", "red", "orange"), pch=16, 
             cex=1.50, ncol=2, xpd=T, bty='n')
    }
    
  }
}

#function to check if data is from the required chromosome
function.catchChromosome <- function(chr, gwas){
  if (gwas == "ad"){
    #take path
    fname = paste("/Users/nicco/Desktop/2k19_work/SNPbrowser/data/AD_CTR/chr", as.character(chr), ".PHENO1.glm.logistic", sep="")
    
    #read data
    dat <- fread(fname, h=T)
    colnames(dat) <- c("chr", "pos", "locus", "ref", "alt", "a1", "a1_frq", "a1_case_frq", "a1_ctr_frq", "r2", "test", "n", "beta", "se", "z-stat", "p")      
    dat <- dat[, c("chr", "pos", "locus", "ref", "alt", "a1", "test", "n", "beta",
                   "se", "z-stat", "p")]
    dat <- dat[!which(is.na(dat$p)),]
  } else if (gwas == "age"){
    #take path
    fname = paste("/Users/nicco/Desktop/2k19_work/SNPbrowser/data/CHC_CTR/chr", as.character(chr), ".PHENO1.glm.logistic", sep="")
    
    #read data
    dat <- fread(fname, h=T)
    colnames(dat) <- c("chr", "pos", "locus", "ref", "alt", "a1", "a1_frq", "a1_case_frq", "a1_ctr_frq", "r2", "test", "n", "beta", "se", "z-stat", "p")      
    dat <- dat[, c("chr", "pos", "locus", "ref", "alt", "a1", "test", "n", "beta",
                   "se", "z-stat", "p")]
    dat <- dat[!which(is.na(dat$p)),]
    
  } else if (gwas == "igap"){
    #take path
    fname = paste("/Users/nicco/Desktop/2k19_work/SNPbrowser/data/IGAP_2k19/chr", as.character(chr), "_IGAP_2k19.txt", sep="")
    
    #read data
    dat <- fread(fname, h=T, stringsAsFactors = F)
    colnames(dat) <- c("iid", "locus", "chr", "pos", "a1", "a2", "beta", "se", "p", "freq_a1", "OR", "95ci", "rsid", "stage")      
    dat <- dat[, c("chr", "pos", "locus", "a1", "a2", "beta", "se", "p")]
    dat <- dat[!which(is.na(dat$p)),]
  } else if (gwas == "metaAge"){
    #take path
    fname = paste("/Users/nicco/Desktop/2k19_work/SNPbrowser/data/META_AGING/chr", as.character(chr), "_meta.txt", sep="")
    
    #read data
    dat <- fread(fname, h=T, stringsAsFactors = F)
    colnames(dat) <- c("chr", "pos", "a1", "a2", "id", "z", "p", "Nsum", "Neff", "dir")      
    dat$locus <- paste(dat$chr, dat$pos, sep=":")
    dat <- dat[, c("chr", "pos", "locus", "a1", "a2", "p")]
    chrom = dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$"-log10(P-value)" <- -log10(as.numeric(dat$p))
    
  }
  return(dat)
}

#function to manage position as input type
function.InputPos <- function(dat, window, input.pos, gwas){
  if (input.pos == "Type position..."){
    #initial position is halfway
    #pos.init <- sample(x = seq(1, max(dat$pos)), size = 1)
    pos.init <- dat[ceiling(nrow(dat)/4), "pos"]
    
    #take data of interest
    snp.info <- dat[which((dat$pos >= pos.init - window) & (dat$pos <= pos.init + window)),]
    if (gwas != "metaAge") {snp.info <- snp.info[which(is.na(snp.info$beta) == FALSE),]  }
    
  } else {
    #split locus to find close positions and define intervals where to find snps
    locus <- str_split_fixed(input.pos, ":", 2)
    interval <- seq(as.numeric(locus[, 2]) - window, as.numeric(locus[, 2]) + window)
    chr <- as.numeric(locus[, 1])
    
    #check if data is relative to that chromosome, and in case change it
    dat <- function.catchChromosome(chr, gwas)
    
    #take data of interest
    snp.info <- dat[which((dat$chr == chr) & (dat$pos %in% interval)), ]
    if (gwas != "metaAge") {snp.info <- snp.info[which(is.na(snp.info$beta) == FALSE),]  }
    
  }
  return(list(snp.info, dat))
}

#function to manage position as input type -- this is for generic gwas data -- has to have at least the following columns:
# chr, pos, p -- other columns for now will be discarded
function.InputPosGeneric <- function(dat, window, input.pos){
  if (input.pos == "Type position..."){
    #initial position is halfway
    pos.init <- dat[ceiling(nrow(dat)/2), "pos"]
    
    #create a "locus" variable as chr:pos
    dat$locus <- paste(dat$chr, dat$pos, sep=":")
    
    #take data of interest
    ind.snp <- grep(pos.init, dat$locus)
    snp.info <- dat[seq(ind.snp - window, ind.snp + window),]
    
  } else {
    #split locus to find close positions and define intervals where to find snps
    locus <- str_split_fixed(input.pos, ":", 2)
    interval <- seq(as.numeric(locus[, 2]) - window, as.numeric(locus[, 2]) + window)
    chr <- as.numeric(locus[, 1])
    
    #take data of interest
    snp.info <- dat[which((dat$chr == chr) & (dat$pos %in% interval)), ]
    
  }
  return(snp.info)
}

#function to manage manual scroll as input type
function.InputManualScroll <- function(dat, window, input.scroll){
  #initial position is all variable
  pos.init <- as.numeric(input.scroll)
  
  #find closest variant
  vectorPos <- seq(pos.init - window, pos.init + window)
  closest <- dat[which(dat$pos %in% vectorPos),]
  closest$diff <- pos.init - closest$pos
  closest <- closest[order(abs(closest$diff)),]
  closest <- closest[which(is.na(closest$p) == FALSE),]
  #  pos.init  <- closest$locus[1]
  
  #take data of interest
  #  ind.snp <- grep(pos.init, dat$locus)
  #  snp.info <- dat[seq(ind.snp - window, ind.snp + window),]
  #  snp.info <- snp.info[which(is.na(snp.info$p) == FALSE)]
  
  return(closest)
}

#function to manage manual scroll as input type -- this is for generic gwas data -- has to have at least the following columns:
# che, pos, p -- other columns for now will be discarded
function.InputManualScrollGeneric <- function(dat, window, input.scroll){
  #initial position is all variable
  pos.init <- as.numeric(input.scroll)
  
  #need to have locus column with chr:pos
  dat$locus <- paste(dat$chr, dat$pos, sep=':')
  
  #find closest variant
  vectorPos <- seq(pos.init - window, pos.init + window)
  closest <- dat[which(dat$pos %in% vectorPos),]
  closest$diff <- pos.init - closest$pos
  closest <- closest[order(abs(closest$diff)),]
  closest$p <- as.numeric(closest$p) 
  closest <- closest[which(is.na(closest$p) == FALSE),]
  
  #  pos.init  <- closest$locus[1]
  #take data of interest
  #  ind.snp <- grep(pos.init, dat$locus)
  #  snp.info <- dat[seq(ind.snp - window, ind.snp + window),]
  
  return(closest)
}

#function to manage genes as input
function.InputGenes <- function(gene){
  #read gene informations -- first command (commented) was to load directly from server -- now I copied the file into app folder
  #gene.db <- fread('sshpass -p "Ciccia22101991!" ssh ntesi@linux-bastion.tudelft.nl "cat /tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/databases/geneLocation_ncbi_hg19/hg19_geneListandPos"')
  gene.db <- fread("/Users/nicco/Desktop/2k19_work/SNPbrowser/data/databases/hg19_geneListandPos", h=T)
  
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
    gene.db <- gene.db[which(gene.db$chrom == "chr19"),]
    #random.numer <- sample(x = seq(1, nrow(gene.db)), size = 1)
    #take corresponding gene
    #gene.info <- gene.db[random.numer, ]
    
    gene.info <- gene.db[ceiling(nrow(gene.db)/4), ]
    
  }
  
  return(gene.info)
}

#function to manage rsID as input
function.rsIDasInput <- function(rsID){
  #check whether there is an input snp, otherwise take a random one -- now is APOE e2
  if (rsID == "Type variant identifier..."){
    #get snp information -- chromosome and position basically
    snp.info <- fread('sshpass -p "Ciccia22101991!" ssh ntesi@linux-bastion.tudelft.nl "grep -w rs7412 /tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/databases/variants_DB/1kGenome/plink/chrAll_red.bim"')
    
    #take position out
    pos <- colnames(snp.info)[2]
    
  } else {
    
    #get snp information -- chromosome and position basically
    cmd <- paste('sshpass -p "Ciccia22101991!" ssh ntesi@linux-bastion.tudelft.nl "grep -w ', rsID, ' /tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/databases/variants_DB/1kGenome/plink/chrAll_red.bim"', sep='')
    snp.info <- fread(cmd = cmd)
    
    #take position out
    pos <- colnames(snp.info)[2]
  }
  
  return(pos)
}

#function to find gwas hits from inputted gene
function.GWASfromGene <- function(dat, gene.info, window, gwas){
  #gene information
  chr <- gene.info$chrom
  chr.n <- str_split_fixed(chr, "chr", 2)[, 2]
  start <- gene.info$txStart - window
  end <- gene.info$txEnd + window
  
  #check if data is relative to that chromosome, and in case change it
  dat <- function.catchChromosome(chr.n, gwas)
  
  #define snps of interest
  snp.info <- dat[which((dat$chr == chr.n) & (dat$pos >= start) & (dat$pos <= end)),]
  if (gwas != "metaAge") {snp.info <- snp.info[which(is.na(snp.info$beta) == FALSE),]  }
  
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
  gene.db <- fread("/Users/nicco/Desktop/2k19_work/SNPbrowser/data/databases/hg19_geneListandPos", h=T)
  
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
  if (nrow(genes) > 3){
    genes$y <- seq(-1, -4, -1)
  } else {
    genes$y <- seq(-1, -3, -1)
  }
  
  return(genes)
}

#function to manage which input data to plot -- still 1 dataset only for now
function.manageInput <- function(inp){
  if (length(inp) == 0){
    dat <- fread("/Users/nicco/Desktop/2k19_work/SNPbrowser/data/AD_CTR/chr19.PHENO1.glm.logistic", h=T, stringsAsFactors = F)
    colnames(dat) <- c("chr", "pos", "locus", "ref", "alt", "a1", "a1_frq", "a1_case_frq", "a1_ctr_frq", "r2", "test", "n", "beta", "se", "z-stat", "p")      
    dat <- dat[, c("chr", "pos", "locus", "ref", "alt", "a1", "test", "n", "beta",
                   "se", "z-stat", "p")]
    dat <- dat[!which(is.na(dat$p)),]
    dat$"-log10(P-value)" <- -log10(as.numeric(dat$p))
    chrom = dat$chr[1]
    gwas = "ad"
    
  } else if (inp == "ad"){
    
    #read input for basic plot
    dat <- fread("/Users/nicco/Desktop/2k19_work/SNPbrowser/data/AD_CTR/chr19.PHENO1.glm.logistic", h=T, stringsAsFactors = F)
    colnames(dat) <- c("chr", "pos", "locus", "ref", "alt", "a1", "a1_frq", "a1_case_frq", "a1_ctr_frq", "r2", "test", "n", "beta", "se", "z-stat", "p")      
    dat <- dat[, c("chr", "pos", "locus", "ref", "alt", "a1", "test", "n", "beta",
                   "se", "z-stat", "p")]
    dat <- dat[!which(is.na(dat$p)),]
    dat$"-log10(P-value)" <- -log10(as.numeric(dat$p))
    chrom = dat$chr[1]
    gwas = "ad"
    
  } else if (inp == "age"){
    dat <- fread("/Users/nicco/Desktop/2k19_work/SNPbrowser/data/CHC_CTR/chr19.PHENO1.glm.logistic", h=T, stringsAsFactors = F)
    colnames(dat) <- c("chr", "pos", "locus", "ref", "alt", "a1", "a1_frq", "a1_case_frq", "a1_ctr_frq", "r2", "test", "n", "beta", "se", "z-stat", "p")      
    dat <- dat[, c("chr", "pos", "locus", "ref", "alt", "a1", "test", "n", "beta",
                   "se", "z-stat", "p")]
    dat <- dat[!which(is.na(dat$p)),]
    dat$"-log10(P-value)" <- -log10(as.numeric(dat$p))
    chrom = dat$chr[1]
    gwas = "age"
    
  } else if (inp == "igap"){
    dat <- fread("/Users/nicco/Desktop/2k19_work/SNPbrowser/data/IGAP_2k19/chr19_IGAP_2k19.txt", h=T, stringsAsFactors = F)
    colnames(dat) <- c("iid", "locus", "chr", "pos", "a1", "a2", "beta", "se", "p", "freq_a1", "OR", "95ci", "rsid", "stage")      
    dat <- dat[, c("chr", "pos", "locus", "a1", "a2", "beta", "se", "p")]
    dat <- dat[!which(is.na(dat$p)),]
    chrom = dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$"-log10(P-value)" <- -log10(as.numeric(dat$p))
    gwas = "igap"
    
  } else if (inp == "toLoad"){
    # dat <- fread("/Users/nicco/Desktop/2k19_work/SNPbrowser/data/trial.ass.txt", h=T, check.names = F)
    # colnames(dat) <- c("chr", "pos", "locus", "ref", "alt", "a1", "test", "n", "beta",
    #                    "se", "z-stat", "p")
  } else if (inp == "metaAge"){
    dat <- fread("/Users/nicco/Desktop/2k19_work/SNPbrowser/data/META_AGING/chr19_meta.txt", h=T, stringsAsFactors = F)
    colnames(dat) <- c("chr", "pos", "a1", "a2", "id", "z", "p", "Nsum", "Neff", "dir")      
    dat$locus <- paste(dat$chr, dat$pos, sep=":")
    dat <- dat[, c("chr", "pos", "locus", "a1", "a2", "p")]
    chrom = dat$chr[1]
    dat$p <- as.numeric(dat$p)
    dat$"-log10(P-value)" <- -log10(as.numeric(dat$p))
    gwas = "metaAge"

  }
  
  return(list(dat, chrom, gwas))
  
}

#MAIN APP
shinyServer(
  function(input, output) {
    #what to plot here
    output$hist <- renderPlot({
      
      ####################
      #manage which input to plot in two cases: beginning of app and 1 gwas as input
      res = function.manageInput(as.character(input$gwas))
      dat = as.data.frame(res[[1]])
      chrom = as.numeric(res[[2]])
      gwas = as.character(res[[3]])
      ####################
      
      ####################
      # #additional input file -- for the moment is 1 only -- LEAVE IT FOR NOW
      # inFile <- input$inp_f
      # n.inputs = 1
      # if (is.null(inFile) == F){
      #   n.inputs = 2
      #   add.f <- fread(inFile$datapath)
      # }
      ####################
      
      ####################
      #Browsing options -- default is Position -- by default take middle of chromosome
      if (input$sel == "Position"){
        #check if only 1 gwas is selected
        if (length(input$gwas) < 2){
          #manage different option of position as input type
          res <- function.InputPos(dat = dat, window = input$x, input.pos = input$pos, gwas)
          snp.info <- as.data.frame(res[[1]])
          dat <- as.data.frame(res[[2]])
          
          #plot
          function.plot(snp.info = snp.info, y.lim = input$y, type = input$sel, plt.type = input$ploType, 
                        windows.number = input$sliding.window, smooth.par = input$smooth, int.locus = input$pos, gwas = gwas)
          ###################
          
          ###################
          #this in case multiple gwas need to be plotted
        } else {
          #first dataset
          d1 <- input$gwas[1]
          res <- function.InputPos(dat = dat, window = input$x, input.pos = input$pos, d1)
          snp.info <- as.data.frame(res[[1]])
          dat <- as.data.frame(res[[2]])
          
          #second dataset
          d2 <- input$gwas[2]
          res = function.manageInput(d2)
          dat.2 = as.data.frame(res[[1]])
          chrom.2 = as.numeric(res[[2]])
          
          #this was for the input file from user
          #snp.info.f2 <- function.InputPosGeneric(dat = dat.2, window = input$x, input.pos = input$pos)
          
          res.2 <- function.InputPos(dat = dat.2, window = input$x, input.pos = input$pos, gwas = d2)
          snp.info.2 <- res.2[[1]]
          dat.2 <- res.2[[2]]
          
          #plot
          function.multiPlot(snp.info = snp.info, snp.info.f2 = snp.info.2, y.lim = input$y, type = input$sel, 
                             plt.type = input$ploType, windows.number = input$sliding.window, smooth.par = input$smooth, int.locus = input$pos,
                             d1, d2)
        }
        ####################
        
        ####################
        #this in case you want manual scroll -- not implemented yet
      } else if (input$sel == "Manual scroll"){
        if (n.inputs == 1){
          #manage different manual scroll input
          snp.info <- function.InputManualScroll(dat = dat, window = input$x, input.scroll = input$all)
          
          #plot
          function.plot(snp.info = snp.info, y.lim = input$y, type = input$sel, plt.type = input$ploType, 
                        windows.number = input$sliding.window, smooth.par = input$smooth, int.locus = NA)
        } else {
          #default dataaset
          snp.info <- function.InputManualScrollGeneric(dat = dat, window = input$x, input.scroll = input$all)
          
          #additional dataset
          snp.info.f2 <- function.InputManualScrollGeneric(dat = add.f, window = input$x, input.scroll = input$all)
          
          #plot
          function.multiPlot(snp.info = snp.info, snp.info.f2 = snp.info.f2, y.lim = input$y, type = input$sel, 
                             plt.type = input$ploType, windows.number = input$sliding.window, smooth.par = input$smooth, int.locus = input$pos)
          
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
                        windows.number = input$sliding.window, smooth.par = input$smooth, int.locus = NA, gwas = gwas)
          ##################
          
          ##################
          #this in case you want to plot multiple gwas
        } else {
          #find position of input gene
          gene.info <- function.InputGenes(gene = as.character(input$gene))
          
          #define the two gwas
          d1 <- input$gwas[1]
          d2 <- input$gwas[2]
          
          #find gwas info around the gene for gwas 1
          snp.info <- function.GWASfromGene(dat = dat, gene.info = gene.info, window = input$x, d1)
          
          #second dataset
          res = function.manageInput(d2)
          dat.2 = as.data.frame(res[[1]])
          chrom.2 = as.numeric(res[[2]])
          
          #additional dataset
          snp.info.2 <- function.GWASfromGene(dat = dat.2, gene.info = gene.info, window = input$x, d2)
          
          #plot
          function.multiPlot(snp.info = snp.info, snp.info.f2 = snp.info.2, y.lim = input$y, type = input$sel, 
                             plt.type = input$ploType, windows.number = input$sliding.window, smooth.par = input$smooth, int.locus = input$pos,
                             d1, d2)
        }
        ####################
        
        ####################
      } else if (input$sel == "Rs ID"){
        if (n.inputs == 1){
          #in case n.inputs == 1 and input type is rsID
          #get locus of the variant
          variant.locus <- function.rsIDasInput(rsID = input$snpID)
          
          #manage different option of position as input type
          snp.info <- function.InputPos(dat = dat, window = input$x, input.pos = variant.locus)
          
          #plot
          function.plot(snp.info = snp.info, y.lim = input$y, type = input$sel, plt.type = input$ploType, 
                        windows.number = input$sliding.window, smooth.par = input$smooth, int.locus = variant.locus)
          
        } else {
          #in case n.inputs > 1 and input type is rsID
          #get locus of the variant
          variant.locus <- function.rsIDasInput(rsID = input$snpID)
          
          #manage different option of position as input type
          snp.info <- function.InputPos(dat = dat, window = input$x, input.pos = variant.locus)
          
          #additional dataset
          snp.info.f2 <- function.InputPosGeneric(dat = add.f, window = input$x, input.pos = variant.locus)
          
          #plot
          function.multiPlot(snp.info = snp.info, snp.info.f2 = snp.info.f2, y.lim = input$y, type = input$sel, 
                             plt.type = input$ploType, windows.number = input$sliding.window, smooth.par = input$smooth, int.locus = input$pos)
        }
      }
    })
  
    #manage click on points
    output$click_info <- renderPrint({

      ####################
      #manage which input to plot in two cases: beginning of app and 1 gwas as input
      res = function.manageInput(as.character(input$gwas))
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
      
      ####################
      #manage which input to plot in two cases: beginning of app and 1 gwas as input
      res = function.manageInput(as.character(input$gwas))
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





