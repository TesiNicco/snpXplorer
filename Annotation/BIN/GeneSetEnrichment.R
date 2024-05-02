# libraries
library(stringr)

# basic paths
MAIN = "/Annotation/"
MAIN_SNP = "/Annotation/RUNS/"
args = commandArgs(trailingOnly=TRUE)

## function to create n sampling dsets
samplingGset <- function(i, mapping){
  all.genes <- strsplit(mapping$geneList, ",")
  tmp_f <- function(x, all.genes){ tmp <- all.genes[[x]]; g <- sample(x=tmp, size=1); return(g)}
  gset <- lapply(1:length(all.genes), tmp_f, all.genes=all.genes)
  gset.clean <- unlist(gset)
  return(gset.clean)
}

## function to perform overlap analysis over the sampling datasets
overlapAnalysis <- function(i, gene_list, source_gset){
  g <- gene_list[[i]]
  g = unique(g)
  res <- suppressMessages(gprofiler2::gost(g, organism = "hsapiens", ordered_query = FALSE, multi_query = FALSE, significant = FALSE, exclude_iea = TRUE,
                               measure_underrepresentation = FALSE, evcodes = FALSE, user_threshold = 1, correction_method = "fdr",
                               domain_scope = "annotated", custom_bg = NULL, numeric_ns = "", sources = source_gset))
  return(res$result)
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
    terms <- all.enrich[, c("term_name", "term_id", "intersection_size", "p_value")]
    slow <- TRUE
  } else {
    terms <- all.enrich[, c("term_name", "term_id", "p_value")]
    slow = FALSE
  }
  # see how many times each term was tested and resctrict to common (tested all times)
  ttt = as.data.frame(table(terms$term_id))
  ttt <- ttt[which(ttt$Freq == n.sampl),]

  # if for some reason the number of samplings was not 500 (e.g. small number of genes), need to adjust for this
  if (nrow(ttt) == 0){
    ttt = as.data.frame(table(terms$term_id))
    mx = max(ttt$Freq)
    ttt <- ttt[which(ttt$Freq == mx),]
  }

  terms <- terms[which(terms$term_id %in% ttt$Var1),]
  #terms <- terms[!duplicated(terms$term_name),]
  avgP <- function(i, terms, term_list, slow){
    if (isTRUE(slow)){
      #all_intersections <- paste0(unique(unlist(strsplit(paste0(terms$intersection[which(terms$term_id == term_list[i])], collapse = ","), ","))), collapse = ",")
      #df <- data.frame(term_name = unique(terms$term_name[which(terms$term_id == term_list[i])]), term_id = term_list[i], avgP = mean(terms$p_value[which(terms$term_id == term_list[i])]), intersection = all_intersections)
      df <- data.frame(term_name = unique(terms$term_name[which(terms$term_id == term_list[i])]), term_id = term_list[i], avgP = mean(terms$p_value[which(terms$term_id == term_list[i])]))
    } else {
      df <- data.frame(term_name = unique(terms$term_name[which(terms$term_id == term_list[i])]), term_id = term_list[i], avgP = mean(terms$p_value[which(terms$term_id == term_list[i])]))
    }
    return(df)
  }
  term_list <- unique(terms$term_id)
  all.terms.avg <- data.table::rbindlist(parallel::mclapply(1:length(term_list), avgP, terms=terms, term_list = term_list, slow = slow, mc.cores=1))
  all.terms.avg <- all.terms.avg[order(all.terms.avg$avgP),]
  all.terms.avg$log10P <- -log10(all.terms.avg$avgP)

  # if intersection were requested, add them back here
  if (slow == TRUE){
    all.terms.avg$intersection = NA
    for (i in 1:nrow(all.terms.avg)){
      cat(paste0('** processing line ', i, '                         \r'))
      all.terms.avg$intersection[i] = paste(unique(unlist(strsplit(paste(all.enrich$intersection[which(all.enrich$term_id == all.terms.avg$term_id[i])], collapse = ','), ','))), collapse = ',')
    }
  }

  #print(head(all.terms.avg))
  #write.table(all.terms.avg, paste("RESULTS_", random_num, "/enrichent_results_sampling.txt", sep=""), quote=F, row.names=F, sep="\t")

  # if (nrow(all.terms.avg[which(all.terms.avg$avgP <= 0.05),]) >0){
  #   pdf(paste("RESULTS_", random_num, "/treemap_sampling_PRS_genes.pdf", sep=""), height=10, width=10)
  #   treemap(all.terms.avg[which(all.terms.avg$avgP <= 0.05),], index="term_name", vSize="log10P", type="index")
  #   invisible(dev.off())
  # }

  l <- list(sign.table, all.terms.avg)
  return(l)
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

# READ ARGUMENTS AND RUN FUNCTION
snps_info_path = args[1]
random_num = args[2]
analysis_mode = args[3]
analysis_mode = unlist(strsplit(analysis_mode, ","))
load(snps_info_path)
annot <- final_res[[1]]
geneList <- final_res[[2]]
# number of iterations
n.sampl <- 300
# first, sample gene-sets
gene_sampling_dsets <- parallel::mclapply(1:n.sampl, samplingGset, mapping=annot, mc.cores=3)
# finally gene-set analysis for all sampling dsets
# first, let's define the databases
dbs = c()
for (i in analysis_mode){
    if (i %in% c("default", "Default")){
        dbs <- c(dbs, "GO:BP")
    } else if (i == "KEGG"){
        dbs <- c(dbs, "KEGG")
    } else if (i %in% c("wiki", "Wiki")){
        dbs <- c(dbs, "WP")
    } else if (i == "Reactome"){
        dbs <- c(dbs, "REAC")
    }
}
enrich.res <- parallel::mclapply(1:n.sampl, overlapAnalysis, gene_list=gene_sampling_dsets, source_gset = dbs, mc.cores=1)
sampling.res <- mergeSampling(enrich.res)
# for other enrichment sources like kegg, reactome etc, need to make a plot
if (length(dbs) >1){
    plotOtherEnrichments(sampling.res, dbs, random_num, type="gost")
} else if (dbs[1] != "GO:BP"){
    plotOtherEnrichments(sampling.res, dbs, random_num, type = "gost")
}
# save sampling res
sampling.res[[(length(sampling.res)+1)]] = dbs
save(sampling.res, file = paste0("/Annotation/RUNS/RESULTS_", random_num, "/tmp_enrichRes.RData"))
# also write text output
write.table(sampling.res[[2]], paste0("/Annotation/RUNS/RESULTS_", random_num, "/Enrichment_results.txt"), quote=F, row.names=F, sep = "\t")

