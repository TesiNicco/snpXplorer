# basic paths
MAIN = "/root/snpXplorer/AnnotateMe/"
MAIN_SNP = "/root/snpXplorer/snpXplorer_v3/"
args = commandArgs(trailingOnly=TRUE)

## function to run positional mapping for the variants with no annotation from consequences and eqtl -- adjusted for faster computations (library-wise)
positionalAnn_generic <- function(mapping){
  # read database of genes
  genes <- data.table::fread(paste(MAIN, "INPUTS_OTHER/NCBI37.3.gene.loc", sep=""), h=F, stringsAsFactors=F)
  colnames(genes) <- c("gene_id", "chr", "start_tss", "stop_tss", "strand", "gene_name")
  #add column to be filled
  mapping$positional_mapping <- NA
  #identify snps to look at
  tmp.exc <- mapping[!(is.na(mapping$coding_snp) & is.na(mapping$eqtl)), ]
  tmp.todo <- mapping[is.na(mapping$coding_snp) & is.na(mapping$eqtl), ]
  tmp.todo$pos <- as.numeric(tmp.todo$pos)
  if (nrow(tmp.todo) >0){
    for (j in 1:nrow(tmp.todo)){
        w = 500000   #set window threshold -- this value will be doubled (upstream and downstream)
        sb.gene <- genes[which(genes$chr == tmp.todo$chr[j]),]       #restrict to chromosome of use
        sb.gene <- sb.gene[which((abs(tmp.todo$pos[j] - sb.gene$start_tss <= w)) | (abs(tmp.todo$pos[j] - sb.gene$stop_tss <= w))), ]       #restrict within 1MB region
        sb.gene$distance_start <- abs(tmp.todo$pos[j] - sb.gene$start_tss)        #calculate distance from start tss
        sb.gene$distance_stop <- abs(tmp.todo$pos[j] - sb.gene$stop_tss)      #calculate distance from stop tss
        sb.gene$distance_start[which(tmp.todo$pos[j] >= sb.gene$start_tss & tmp.todo$pos[j] <= sb.gene$stop_tss)] <- 0      #also assign 0 for intronic positions
        sb.gene$min_distance <- apply(X=sb.gene[, c(7, 8)], MARGIN=1, FUN=min)      #function to get the minimum distance between t-start-s and t-stop-site
        sb.gene$code <- apply(X=sb.gene[, c(9)], MARGIN=1, FUN=function(x){return(round(x/(w/10), 0))})       #assign code (intron=0, <50kb=1, <100kb=2, etc)
        sb.gene <- sb.gene[order(sb.gene$code),]
        if (nrow(sb.gene) > 0){  #if there are hits, check overlap with gene indicated from cadd, save results and stop here
        #tmp <- c(sb.gene$gene_name[which(sb.gene$code == min(sb.gene$code))], unlist(strsplit(tmp.todo$snp_conseq_gene[j], ",")))
        #tmp <- tmp[!duplicated(tmp) & !is.na(tmp)]
        #tmp.todo$positional_tmp.todo[j] <- paste0(tmp, collapse=",")
        # maybe including also cadd introduce too much variablility? try without here
        tmp.todo$positional_mapping[j] <- paste0(sb.gene$gene_name[which(sb.gene$code == min(sb.gene$code))], collapse=",")
        } else {
        tmp.todo$positional_mapping[j] <- tmp.todo$snp_conseq_gene[j]
        }
    }
    mapping <- rbind(tmp.exc, tmp.todo)      #re-join dataset
  } else {
      mapping = tmp.exc
  }
  return(mapping)
}

## function to clean file after mapping procedure and output the list of genes -- adjusted for faster computations (library-wise)
getGeneList_mod_generic <- function(mapping){
  #put information about which source will be used to get genes: priority is: 1-coding_snp -- 2-eqtl -- 3-positional
  mapping$source_finalGenes <- NA
  mapping$geneList <- NA

  #assign source and gene list for coding variants
  mapping$source_finalGenes[which(!is.na(mapping$coding_snp))] <- "coding"
  mapping$geneList[which(!is.na(mapping$coding_snp))] <- mapping$snp_conseq_gene[which(!is.na(mapping$coding_snp))]

  # before assigning, need to exclude pseudogenes (grep a . in the gene name)
  pseudo <- mapping[grep("\\.", mapping$eqtl),]
  if (nrow(pseudo) >1){
    for (k in 1:nrow(pseudo)){
      tmp <- unlist(strsplit(pseudo$eqtl[k], ","))
      tmp_tis = unlist(strsplit(pseudo$eqtl_tissue[k], ","))
      if (length(tmp) >1){
        todel <- grep("\\.", tmp)
        tmp <- tmp[-todel]
        tmp_tis = tmp_tis[-todel]
        pseudo$eqtl[k] <- paste0(tmp, collapse=",")
        pseudo$eqtl_tissue[k] <- paste0(tmp_tis, collapse = ",")
      } else {
        pseudo$eqtl[k] <- NA
        pseudo$eqtl_tissue[k] <- NA
      }
      mapping[which(mapping$locus == pseudo$locus[k]), ] <- pseudo[k, ]
    }
  }

  # assign source for eqtl hits
  mapping$source_finalGenes[which(!is.na(mapping$eqtl))] <- "eqtl+cadd"
  for (l in 1:nrow(mapping)){
    if ((!is.na(mapping$source_finalGenes[l])) & mapping$source_finalGenes[l] == "eqtl+cadd"){
      if (is.na(mapping$snp_conseq_gene[l]) | (mapping$snp_conseq_gene[l] == mapping$eqtl[l])){
        mapping$geneList[l] <- mapping$eqtl[l]
      } else {
        tmp_eqtl_list <- unlist(strsplit(mapping$eqtl[l], ","))
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
  mapping$source_finalGenes[which(is.na(mapping$coding_snp) & is.na(mapping$eqtl))] <- "positional"
  #mapping$geneList[which(is.na(mapping$coding_snp) & is.na(mapping$eqtl_blood))] <- mapping$positional_mapping[which(is.na(mapping$coding_snp) & is.na(mapping$eqtl_blood))]
  # for these annotations, also consider the gene suggested by cadd
  tmp_posit <- mapping[which(mapping$source_finalGenes == "positional"),]
  if (nrow(tmp_posit) > 0){
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
  }


  #assign finally those that did not map anywhere
  mapping$source_finalGenes[which(is.na(mapping$coding_snp) & is.na(mapping$eqtl) & is.na(mapping$positional_mapping))] <- "missing"

  tmp1 <- mapping$geneList[which(!is.na(mapping$geneList))]
  tmp2 <- strsplit(tmp1, ",")
  all.genes <- c()
  for (i in 1:length(tmp2)) { all.genes <- c(all.genes, tmp2[[i]]) }
  all.genes <- all.genes[!is.na(all.genes)]
  all.genes <- all.genes[which(!(all.genes %in% c("NA", "character(0)")))]
  all.genes <- all.genes[!duplicated(all.genes)]
  all.genes <- stringr::str_split_fixed(all.genes, "\\.", 2)
  all.genes <- all.genes[, 1]

  #merge results
  l = list(all.genes, mapping)

  return(l)
}

## function to remove NAs from the genelists
removeNAgeneList <- function(annot){
  for (i in 1:nrow(annot)){
    tmp = unlist(strsplit(annot$geneList[i], ","))
    annot$geneList[i] = paste(tmp[which(tmp != "NA")], collapse = ",")
  }
  annot[which(annot$geneList == ""), "geneList"] <- NA
  return(annot)
}

# read arguments
snps_info_path = args[1]
# run function to read snps info
load(snps_info_path)
# run functions
full.annot <- positionalAnn_generic(mapping=out.gtex)
res.clean <- getGeneList_mod_generic(mapping=full.annot)
annot = removeNAgeneList(res.clean[[2]])
final_res = list(annot, res.clean[[1]])

# outputs
save(final_res, file = snps_info_path)
