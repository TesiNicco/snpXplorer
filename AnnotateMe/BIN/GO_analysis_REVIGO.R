# basic paths
MAIN = "/root/snpXplorer/AnnotateMe/"
MAIN_SNP = "/root/snpXplorer/snpXplorer_v3/"
args = commandArgs(trailingOnly=TRUE)

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
  if (nrow(sig) >1){
    # output revigo input
    write.table(sig[, c("term_id", "avg_p")], paste("RESULTS_", random_num, "/revigo_inp.txt", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
    path_finp = paste0("/root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/revigo_inp.txt")

    # set parameters for REVIGO
    # distance.meas can be: Lin -- SIMREL
    # clust.simil can be: 0.7 -- 0.5 -- 0.4
    distance.meas <- "LIN"
    clust.simil <- 0.40

    # run python script to get results
    cat("## Running REVIGO\n")
    cmd = paste0("Rscript ", MAIN, "BIN/parseREVIGO_new.R ", path_finp, " ", clust.simil, " ", distance.meas, " ", random_num)
    system(cmd, wait = T, ignore.stdout = T, ignore.stderr = T)

    # read output
    d <- data.table::fread(paste("/root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/revigo_out.csv", sep=""), h=T, sep=",")
    tmp = sig[, c("term_id", "avg_p")]
    d = merge(d, tmp, by.x = "TermID", by.y = "term_id")
    # remove input as it is not too informative
    #cmd <- "rm RESULTS/revigo_inp.txt"
    #system(cmd)

    pdf(paste("RESULTS_", random_num, "/clustering_GO_terms_REVIGO.pdf", sep=""), height=11, width = 11)
    revigoPlot(d)
    dev.off()
  } else {
    print(thr)
    cat("  !!! No significant GO terms enriched\n")
    d <- "!!! No significant GO terms enriched in the submitted list of SNPs."
    write.table(d, paste("RESULTS_", random_num, "/revigo_out.csv", sep=""), sep="\t")
  }
  return("Done")
}

## function to plot REVIGO results -- updated
revigoPlot <- function(one.data){
  # take only terms that were not merged
  one.data <- one.data[which(one.data$Eliminated == FALSE),]

  # do some adjustements
  one.data$PlotX <- as.numeric(one.data$PlotX)
  one.data$PlotY <- as.numeric(one.data$PlotY)

  ## base plot
  plot(0, 0, pch=16, col="white", xlim=c(min(one.data$PlotX)-1, max(one.data$PlotX)+1), ylim=c(min(one.data$PlotY)-1, max(one.data$PlotY)+1), xlab="Semantic space X", ylab="Semantic space Y", cex.lab=1.50, xpd=T)

  ## grid
  for (i in seq(min(one.data$PlotX)-1, max(one.data$PlotX)+1, (max(one.data$PlotX)-min(one.data$PlotX)+2)/10)){abline(v=i, lwd=0.4, col="grey80")}
  for (i in seq(min(one.data$PlotY)-1, max(one.data$PlotY)+1, (max(one.data$PlotY)-min(one.data$PlotY)+2)/10)){abline(h=i, lwd=0.4, col="grey80")}

  ## colors
  colz <- viridis::viridis(n = nrow(one.data), option = "plasma")
  #colz <- colorRampPalette(c("red", "orange", "yellow"))
  #colorz <- colz(nrow(one.data))
  #one.data <- one.data[order(one.data$'log10 p-value'),]
  one.data <- one.data[order(one.data$Value),]
  #one.data$col <- colorz
  one.data$col <- colz

  ## adjust sizes
  one.data$size = 3
  one.data$size[which(one.data$LogSize >3)] <- one.data$LogSize[which(one.data$LogSize >3)] + 3
  ## points
  points(x = one.data$PlotX, y = one.data$PlotY, pch=16, cex=one.data$size, col=ggplot2::alpha(one.data$col, 0.65), lwd=2)

  # need to add some \n in the name description
  one.data$upd_name <- NA
  for (i in 1:nrow(one.data)){
    # count the characters
    tmp.string <- as.character(one.data$Name[i])

    # calculate in how many lines to distribute the name -- max is 23 characters (including spaces) per line
    n.lines <- ceiling(nchar(tmp.string)/25)

    # divide by space depending on n.lines
    tmp <- unlist(strsplit(as.character(one.data$Name[i]), " "))
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

  ## set how many data to annotate -- top n (ordered by logsize)
  one.data = one.data[order(one.data$Value),]
  n <- 8

  ## set position for annotation
  sb = head(one.data, n)
  basicPlotteR::addTextLabels(sb$PlotX, sb$PlotY, sb$upd_name, col.label="black")

  ## legend
  mx.x <- max(one.data$PlotX)+2 + max(one.data$PlotX)*0.2
  mx.y <- max(one.data$PlotY)-2
  step <- 0.25
  text(x = mx.x, y = mx.y, labels = "-Log10(P)", cex=0.80, xpd=T, font=2)
  plotrix::gradient.rect(xleft = mx.x-(step*2), ybottom = mx.y-(step*9), xright = mx.x+(step*2), ytop = mx.y-(step*2), col=rev(colz), nslices = length(colz), gradient = "horizontal")
  text(x = mx.x+(step*4), y = mx.y-(step*2), labels = round(max(abs(one.data$Value)), 1), adj = 0.5, xpd=T, font=2, cex=0.80)
  text(x = mx.x+(step*4), y = mx.y-(step*5.5), labels = round(median(abs(one.data$Value)), 1), adj = 0.5, xpd=T, font=2, cex=0.80)
  text(x = mx.x+(step*4), y = mx.y-(step*9), labels = round(min(abs(one.data$Value)), 1), adj = 0.5, xpd=T, font=2, cex=0.80)
  #text(x = mx.x, y = mx.y+(step*5)+2.8, labels = "B", xpd=T, font=2, cex=3)

  ## legend for circle size
  # text(x = mx.x, y = mx.y-(step*15), labels = "Term size", cex=1.20, font=2, xpd=T)
  # points(x = mx.x, y = mx.y-(step*18), pch=16, col="grey80", cex=min(one.data$PlotSize+5), xpd=T)
  # points(x = mx.x, y = mx.y-(step*23), pch=16, col="grey80", cex=median(one.data$PlotSize+5), xpd=T)
  # points(x = mx.x, y = mx.y-(step*29), pch=16, col="grey80", cex=max(one.data$PlotSize+5), xpd=T)
}

# READ ARGUMENTS AND RUN FUNCTION
random_num = args[1]
load(paste0("RESULTS_", random_num, "/tmp_enrichRes.RData"))
go_data = sampling.res[[2]]
colnames(go_data)[1:4] <- c("term_name", "term_id", "avg_p", "log10p")
# write go_data
write.table(go_data, file = paste0("RESULTS_", random_num, "/geneSet_enrichment_results_and_clusters.txt"), quote=F, row.names = F, sep = "\t")
# add source information to select only go
go_data$source_gset = stringr::str_split_fixed(go_data$term_id, ":", 2)[, 1]
go_data <- go_data[which(go_data$source_gset == "GO"),]
# last thing is to run REVIGO and make the plot
revigo_res <- function(go_data, thr = 0.01) {
    out <- tryCatch({revigo(avg_pvalues = go_data, thr = 0.01)},
            error=function(cond) {
                message("Connection to REVIGO does not seem to work..Skipping analysis!")
                # Choose a return value in case of error
                return(NA)
            })
    return(out)
}
final_res = revigo_res(go_data, thr = 0.01)
