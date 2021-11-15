# basic paths
MAIN = "/root/snpXplorer/AnnotateMe/"
MAIN_SNP = "/root/snpXplorer/snpXplorer_v3/"
args = commandArgs(trailingOnly=TRUE)

## function to find sv close to input snps and genes -- adjusted for faster computations (library-wise)
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

# READ ARGUMENTS AND RUN FUNCTION
snps_info_path = args[1]
random_num = args[2]
load(snps_info_path)
annot <- final_res[[1]]
geneList <- final_res[[2]]
all_str_hg38 = data.table::fread("/root/snpXplorer/AnnotateMe/INPUTS_OTHER/StrVar_hg38.txt.gz", h=T, sep = "\t", stringsAsFactors=F)
all_sv = data.table::rbindlist(parallel::mclapply(1:nrow(annot), function_findSVs, annot = annot, all_str_hg38 = all_str_hg38, mc.cores=1))
write.table(all_sv, paste0("RESULTS_", random_num, "/SNP_and_SV_overlap.txt"), quote=F, row.names=F, sep="\t")
