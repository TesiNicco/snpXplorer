## FIND POSITION TO PLOT
  # browsing option is 'Locus, Gene, RsID' 
  directInput = function(target, window, gene.db, genes.hg38, ref_version, SNPlocs.Hsapiens.dbSNP144.GRCh37, MAIN_PATH){
    span = 100000
    if (target %in% c('Type locus, rsID or gene name', '')){
      chrom = 16; pos1 = 12000000 - window; pos2 = 12500000 + window; snp_interest = NA
    } else if (length(grep(":", target)) >0){    # input is chr:pos
      chrom = stringr::str_split_fixed(target, ":", 2)[, 1]
      chrom = stringr::str_replace_all(chrom, 'chr', '')
      pos1 = as.numeric(stringr::str_split_fixed(target, ":", 2)[, 2]) - window
      pos2 = as.numeric(stringr::str_split_fixed(target, ":", 2)[, 2]) + window
      snp_interest = as.numeric(stringr::str_split_fixed(target, ":", 2)[, 2])
    } else if (length(grep("rs", target)) >0){
      snp_info <- NULL
      try(snp_info <- as.data.frame(snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37, target)), silent = T)
      if (!is.null(snp_info)){ chrom = snp_info$seqnames; pos1 = as.numeric(snp_info$pos) - window; pos2 = as.numeric(snp_info$pos) + window; snp_interest = snp_info$pos } else { chrom = NA; pos1 = NA; pos2 = NA; snp_interest = NA }
      if (ref_version == 'GRCh38 (hg38)'){ data_lifted = liftOver_data(chrom = as.character(chrom), start = pos1, end = pos2, type = "interval", p = NULL, MAIN_PATH, from = 'hg19'); chrom = data_lifted[[1]]; pos1 = data_lifted[[2]]; pos2 = data_lifted[[3]] }
    } else {
      if (ref_version == "GRCh37 (hg19)"){ tmp = gene.db[which(gene.db$"#geneName" == toupper(target)),] } else { tmp = genes.hg38[which(genes.hg38$"#geneName" == toupper(target)),] }
      if (nrow(tmp) >0){ chrom = tmp$chrom[1]; pos1 = as.numeric(tmp$txStart[1]) - window; pos2 = as.numeric(tmp$txEnd[1]) + window } else { chrom = NA; pos1 = NA; pos2 = NA }
      snp_interest = NA
    }
    df = data.frame(chrom = str_replace_all(chrom, 'chr', ''), start = pos1, end = pos2)
    df$chrom = as.character(df$chrom)
    results = list(df, snp_interest)
    return(results)
  }

## FIND RSID OF SNPS
  findRsID <- function(region_of_interest, reference_genome, span, MAIN_PATH, mode){
    #tmp = fread(paste0(MAIN_PATH, 'data/databases/snps_info/chr', region_of_interest$chrom, '_snps_info.txt.gz'), h=T, stringsAsFactors=F)
    region_to_search = paste0(region_of_interest$chr, ':', region_of_interest$start, '-', region_of_interest$end)
    if (reference_genome == 'GRCh37 (hg19)'){ 
      target_ref = paste0(MAIN_PATH, '/data/databases/snps_info/chrAll_snps_info_hg19.txt.gz')
    } else {
      target_ref = paste0(MAIN_PATH, 'data/databases/snps_info/chrAll_snps_info_hg38.txt.gz')
    }
    tmp = tabix(region_to_search, target_ref, check.chr = F, verbose = FALSE)
    colnames(tmp) = c("POS", "POS_HG38", "ID", "REF", "ALT", "ALT_FREQS", "CHR")
    return(tmp)
  }

## EXTRACT GENES FROM A REGION
  findGenes <- function(region_of_interest, reference_genome, genes_hg19, genes_hg38){
    span = 1000
    if (reference_genome == 'GRCh37 (hg19)'){ tmp = genes_hg19 } else { tmp = genes_hg38 }
    genes_in_region = tmp[which(tmp$chrom == paste0('chr', region_of_interest$chrom)),]
    genes_in_region = genes_in_region[which((genes_in_region$txStart >= (region_of_interest$start - span)) & (genes_in_region$txEnd <= (region_of_interest$end + span))),]
    if (nrow(genes_in_region) >0){
      colnames(genes_in_region)[which(colnames(genes_in_region) == "#geneName")] = 'Gene'
      genes_in_region = genes_in_region[order(genes_in_region$txStart),]
      # add y-position for the plot
      n = ceiling(nrow(genes_in_region)/2.5)
      if (n != 0){v <- -1; genes_in_region$y <- NA; for (i in 1:nrow(genes_in_region)){ genes_in_region$y[i] <- v; v <- v-1; if (v == -n-1){ v = -1 } } }
    } else {
      genes_in_region = data.frame(chrom = as.character(), pos = as.character())
    }
    return(genes_in_region)
  }
## EXTRACT SV TO PLOT
  findSVs <- function(region_of_interest, reference_genome, span_value, sv_source, MAIN_PATH){
    if (reference_genome == 'GRCh37 (hg19)'){ ref = paste0(MAIN_PATH, '/data/databases/Structural_variants/str_set_hg19.txt.gz') } else { ref = paste0(MAIN_PATH, '/data/databases/Structural_variants/str_set_hg38.txt.gz') }
    region_to_search = paste0('chr', region_of_interest$chrom, ':', region_of_interest$start, '-', region_of_interest$end)
    svs = tabix(region_to_search, ref, verbose = FALSE)
    if (reference_genome == 'GRCh37 (hg19)'){
      colnames(svs) = c('chrom', 'start_pos', 'end_pos', 'type', 'col', 'diff_alleles', 'source')
    } else {
      colnames(svs) = c('chrom', 'start_pos', 'end_pos', 'diff_alleles', 'type', 'col', 'source')
    }
    #if (reference_genome == 'GRCh37 (hg19)'){ tmp = all_str } else { tmp = all_str_hg38; tmp$end_pos = as.numeric(tmp$end_pos) }
    #svs = tmp[which(tmp$chr == paste0('chr', region_of_interest$chrom)),]
    #svs = svs[which(svs$start_pos >= (region_of_interest$start - span_value) & svs$end_pos <= (region_of_interest$end + span_value)),]
    if (sv_source != 'all'){ svs = svs[which(svs$source == sv_source),] }
    if (nrow(svs) >0){
      svs$end_pos_plus = as.numeric(svs$end_pos) + as.numeric(svs$diff_alleles)
      svs = svs[order(svs$start_pos),]
      # add y-position for the plot
      n = ceiling(nrow(svs)/2.5)
      if (n != 0){v <- -1; svs$y <- NA; for (i in 1:nrow(svs)){ svs$y[i] <- v; v <- v-1; if (v == -n-1){ v = -1 } } }
      svs$source[which(svs$source == "jasper")] <- "Linthorst et al., 2020"
      svs$source[which(svs$source == "audano")] <- "Audano et al., 2019"
      svs$source[which(svs$source == "chaisson")] <- "Chaisson et al., 2019"
    } else {
      svs = data.frame(chrom = as.character(), pos = as.character())
    }
    return(svs)
  }
## EXTRACT DATA TO PLOT
  # function to extract data to plot
  extractDataForPlot = function(gwas_to_plot, region_of_interest, span, colors, res_example, reference_genome, MAIN_PATH){
    region_to_search = paste0(region_of_interest$chrom, ':', region_of_interest$start, '-', region_of_interest$end)
    data_to_plot = list()                 # initialize the list that will contain data to plot
    main_path = paste0(MAIN_PATH, 'data/')
    if (is.null(gwas_to_plot)){           # this is the case when no input gwas are selected --> example data
      if (region_of_interest$chrom %in% c(16, 17, 18, 19, 20, 21)){
        tmp = res_example
        colnames(tmp) = c("pos", "chrom", "p", "pos_hg38", "rsid", "maf", "ref", "alt")
        # in case, we need to liftover
        if (reference_genome == 'GRCh38 (hg38)'){ 
          #data_lifted = liftOver_data(chrom = tmp$chrom, start = tmp$pos, end = tmp$pos + 1, type = "gwas", p = tmp$p, MAIN_PATH, from = 'hg19'); chrom_lf = data_lifted[[1]]; pos_lf = data_lifted[[2]]; p_lf = data_lifted[[3]]; tmp = data.table(chrom = chrom_lf, pos = pos_lf, p = p_lf) }
          tmp = tmp[which((tmp$pos_hg38 >= region_of_interest$start - span) & (tmp$pos_hg38 <= region_of_interest$end + span)),]
          tmp$pos = NULL; colnames(tmp)[3] = 'pos'; tmp = tmp[!is.na(tmp$pos),]
        } else {
          # then match based on the position
          tmp = tmp[which((tmp$pos >= region_of_interest$start - span) & (tmp$pos <= region_of_interest$end + span)),]       
        }
        tmp$col = colors[1]
        tmp$name = 'Example'
        data_to_plot[[(length(data_to_plot) + 1)]] = tmp
        plotError = FALSE
      } else {
        plotError = TRUE; data_to_plot = NULL
      }
    } else {
      for (i in 1:length(gwas_to_plot)){
        if (length(grep('/var/', gwas_to_plot[i])) >0){
          tmp = fread(gwas_to_plot[i], h=T, stringsAsFactors=F)
          res = identiHeader(tmp); plotError = res[[1]]; tmp = res[[2]]
          tmp$name = 'Uploaded GWAS'
          tmp = tmp[which(tmp$chrom == region_of_interest$chrom),]
        } else {
          if (reference_genome == 'GRCh38 (hg38)'){
            tmp = tabix(region_to_search, paste0(main_path, gwas_to_plot[i], '/chr', region_of_interest$chrom, '_', gwas_to_plot[i], '_hg38.txt.gz'), check.chr = F, verbose = FALSE)
            colnames(tmp) = c('chrom', 'p', "pos", "rsid", "maf", "ref", "alt")
          } else {
            tmp = tabix(region_to_search, paste0(main_path, gwas_to_plot[i], '/chr', region_of_interest$chrom, '_', gwas_to_plot[i], '.txt.gz'), check.chr = F, verbose = FALSE)
            colnames(tmp) = c('pos', 'chrom', 'p', "pos_hg38", "rsid", "maf", "ref", "alt")
          }
          plotError = FALSE
          tmp$name = str_replace_all(gwas_to_plot[i], '_', ' ')
        }
        #tmp = tmp[which((tmp$pos >= region_of_interest$start - span) & (tmp$pos <= region_of_interest$end + span)),]
        tmp$col = colors[i]; tmp$pos = as.numeric(tmp$pos); tmp$p = as.numeric(tmp$p)
        data_to_plot[[i]] = tmp
      }
    }
    results = list(plotError, data_to_plot)
    return(results)
  }

  # function to identify header in uploaded file
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
    plotError = FALSE
    if (length(chrom.idx) + length(pos.idx) + length(p.idx) == 3){
      head[chrom.idx] <- "chr"; head[pos.idx] <- "pos"; head[p.idx] <- "p"; colnames(dat) <- head
      dat <- dat[, c("chr", "pos", "p")]
      colnames(dat) = c("chrom", "pos", "p")
      dat$chrom = str_replace_all(dat$chrom, 'chr', '')
      dat$pos = as.numeric(as.character(dat$pos))
      dat$p = as.numeric(as.character(dat$p))
    } else {
      plotError = TRUE
      dat = data.frame(chrom = as.character(), pos = as.character(), p = as.character())
    }
    res = list(plotError, dat)
    return(res)
  }

## EXTRACT RECOMBINATION RATES
  # function to extract recombination rates to plot
  extractRecombination = function(region_of_interest, reference_genome, span, MAIN_PATH){
    #if (reference_genome == 'GRCh37 (hg19)'){ tmp = genetic_map } else { tmp = genetic_map_hg38 }
    if (reference_genome == 'GRCh37 (hg19)'){ ref = paste0(MAIN_PATH, '/data/databases/Recombination_rates/recombination_rates_hg19.txt.gz') } else { ref = paste0(MAIN_PATH, '/data/databases/Recombination_rates/recombination_rates_hg38.txt.gz') }
    region_to_search = paste0('chr', region_of_interest$chrom, ':', as.character(region_of_interest$start), '-', as.character(region_of_interest$end))
    tmp = tabix(region_to_search, ref, verbose = FALSE)
    #tmp = tmp[[as.numeric(region_of_interest$chrom)]]
    colnames(tmp) = c("chr", "pos", "combined", "cm")
    tmp$pos = as.numeric(tmp$pos); tmp$combined = as.numeric(tmp$combined)
    #tmp = tmp[which((tmp$pos >= region_of_interest$start - span) & (tmp$pos <= region_of_interest$end + span)),]
    return(tmp)
  }

## DENSITY OF PVALUES
  # function to find density line of p-value across chromosomal position
  DensityLinePvalue <- function(snp.info, wind.n, smooth.par){
    # define results
    res = list()
    # define intervals
    intervals = seq(min(snp.info$pos), max(snp.info$pos), length.out = wind.n)
    # loop through values of the intervals
    for (i in 1:wind.n){
      # take data in intervals
      sb = snp.info[which(snp.info$pos >= intervals[i] & snp.info$pos <= intervals[i+1]),]
      sb$p = as.numeric(sb$p); sb = sb[!is.na(sb$p),]
      if (nrow(sb) >0){
        # create results
        tmp_res = data.frame(window = intervals[i] + (intervals[i + 1] - intervals[i])/2, pvalue = -log10(min(sb$p)), error = as.numeric(quantile(sb$p, na.rm=T, probs = 0.025)))
        res[[i]] = tmp_res
      } else {
        tmp_res = data.frame(window = intervals[i] + (intervals[i + 1] - intervals[i])/2, pvalue = 0, error = 0)
      }
    }
    res = rbindlist(res)
    # remove NA and substitute them with 0
    res$error[is.na(res$error)] = 0
    # local regression on maximal values per bin and prediction
    options(warn = -1)
    lo <- loess(res$pvalue ~ res$window, span=smooth.par)
    xl <- seq(min(res$window), max(res$window), (max(res$window) - min(res$window))/100)
    pred <- predict(lo, xl)
    pred[which(pred < 0)] <- 0
    # add limits -- left and right for polygon function
    xl <- c(xl[1], xl, xl[length(xl)])
    pred <- c(0, pred, 0)
    options(warn = 0)
    # make object to plot
    toplot = data.frame(x = xl, y = pred)
    return(toplot)
  }

## PLOTS
  # function to plot points and recombination rates -- based on plotly
  Plot_plotly <- function(reference_genome, region_of_interest, rsid_region, snps_data, snp_interest, recomb_data, significance, pos_start, pos_end, plot_type, recomb, genes_in_region, svs_in_region, showExons, ld){
    # PLOT 1 IS THE MAIN SNP-PLOT
      # initialize the plot
      fig1 <- plot_ly()
      # plot data depending on whether Scatter or Density is requested
      for (snp_group in snps_data){
        if (plot_type == 'Scatter'){
          snp_group$pos = as.numeric(snp_group$pos)
          snp_group$p = as.numeric(snp_group$p)
          if (reference_genome == 'GRCh37 (hg19)'){ snp_group = merge(snp_group, rsid_region, by.x = 'pos', by.y = "POS", all.x = T) } else { snp_group = merge(snp_group, rsid_region, by.x = 'pos', by.y = "POS_HG38", all.x = T) }
          snp_group$labels = paste0("<b>ID:<b> ", snp_group$ID, "<br><b>Position:<b> %{x} <br><b>Log(p):<b> %{y} <br><b>REF:<b> ", snp_group$REF, "<br><b>ALT:<b> ", snp_group$ALT, "<br><b>MAF:<b> ", snp_group$ALT_FREQS)
          fig1 = fig1 %>% add_trace(data=snp_group, name=unique(snp_group$name), x=~pos, y=~-log10(p), type='scatter', mode='markers', yaxis='y1', marker=list(color=~col, size=8), hovertemplate = ~labels)
          # LD here below
          if (!is.null(ld)){
            ld = merge(ld, snp_group, by.x = 'BP_B', by.y = 'pos')
            ld$labels = paste0("<b>ID:<b> ", ld$ID, "<br><b>Position:<b> %{x} <br><b>Log(p):<b> %{y} <br><b>REF:<b> ", ld$REF, "<br><b>ALT:<b> ", ld$ALT, "<br><b>MAF:<b> ", ld$ALT_FREQS, "<br><b>LD with:<b> ", ld$SNP_A, "<br><b>R2:<b> ", ld$R2)
            if (nrow(ld) >0){ fig1 = fig1 %>% add_trace(data=ld, name='Linkage disequilibrium', x=~BP_B, y=~-log10(p), type='scatter', mode='markers', yaxis='y1', marker=list(color=~colo, size=10, symbol = 'diamond'), hovertemplate = ~labels) }
          }
          # Input SNP here below
          if (!is.na(snp_interest)){ tmp = snp_group[which(snp_group$pos == snp_interest),]; fig1 = fig1 %>% add_trace(data=tmp, name='SNP of interest', x=~pos, y=~-log10(p), type='scatter', mode='markers', yaxis='y1', marker=list(color=~col, size=16), hovertemplate = "<b>Position:<b> %{x} <br><b>Log(p):<b> %{y}") }
        } else {
          df = densityLinePvalue(snp_group, 20, 0.1)
          fl = paste0('rgba(', paste(as.vector(col2rgb(unique(snp_group$col))), collapse = ", "), ', 0.5)')
          fig1 = fig1 %>% add_trace(data = df, name=unique(snp_group$name), x = ~x, y = ~y, type = 'scatter', mode = 'lines', fill='tozeroy', fillcolor=fl, hovertemplate = "<b>Position:<b> %{x} <br><b>Log(p):<b> %{y}")
        }
      }
      # add data from second dataframe recomb_data
      if (recomb == "Yes"){ 
        fig1 = fig1 %>% add_trace(x=recomb_data$pos, y=recomb_data$combined, name="Recombination rate", type = 'scatter', mode = 'lines', yaxis = 'y2', line = list(color = 'orange', width = 2), hovertemplate = "<b>Position:<b> %{x} <br><b>Rate:<b> %{y}") 
      }
      # show figure
      plt_title = paste0('chr', region_of_interest$chrom, ':', format(pos_start, scientific = F, digits = 0), '-', format(pos_end, scientific = F, digits = 0))
      fig1 = fig1 %>% layout(plot_bgcolor='#e5ecf6', title = list(text = plt_title, xanchor = 'center'), 
                    yaxis2 = list(zerolinecolor = '#ffff', gridcolor = 'ffff', zerolinewidth = 4, title = "Recombination rates", autorange = FALSE, range = c(0, 100), overlaying = "y2", side = "right"),
                    yaxis = list(zerolinecolor = '#ffff', gridcolor = 'ffff', zerolinewidth = 4, title = "-log10(P-value)", autorange = FALSE, range = c(0, significance)), 
                    xaxis = list(zerolinecolor = '#ffff', gridcolor = 'ffff', zerolinewidth = 4, title = "Genomic position (bp)", autorange = FALSE, range = c(pos_start, pos_end)))
      fig1 = fig1 %>% layout(legend = list(bgcolor = 'rgba(0,0,0,0)'))
    # PLOT 2 IS THE GENE TRACK
      fig2 <- plot_ly()
      y_space = 0.125
      genes_in_region$hover_label = paste0("<b>Gene name:<b> ", genes_in_region$Gene, "<br><b>Transcript:<b> ", genes_in_region$name, "<br><b>Strand:<b> ", genes_in_region$strand, "<br><b>TxStart:<b> ", genes_in_region$txStart, "<br><b>TxEnd:<b> ", genes_in_region$txEnd, "<br><b>Exons:<b> ", genes_in_region$exonCount, "<extra></extra>")
      for (i in 1:nrow(genes_in_region)){
        fig2 = fig2 %>% add_trace(x=c(genes_in_region$txStart[i], genes_in_region$txEnd[i]), y=rep(genes_in_region$y[i], 2), type = 'scatter', mode = 'lines', showlegend = F, line = list(width = 4), hovertemplate = genes_in_region$hover_label[i])
        if (nrow(genes_in_region) <= 20){ txt_size = 15 } else if (nrow(genes_in_region) <= 30){ txt_size = 12 } else if (nrow(genes_in_region) <= 40){ txt_size = 10 } else if (nrow(genes_in_region) <= 50){ txt_size = 6 }
        if (nrow(genes_in_region) <= 50){ fig2 = fig2 %>% add_text(x=(genes_in_region$txStart[i] + (genes_in_region$txEnd[i] - genes_in_region$txStart[i])/2), y=genes_in_region$y[i]+y_space*1.3, text = genes_in_region$Gene[i], textposition = 'middle', textfont = list(color = '#000000', size = txt_size), showlegend = F, hoverinfo="none") }
        exon_start = as.numeric(str_split(genes_in_region$exonStarts[i], ',')[[1]]); exon_end = as.numeric(str_split(genes_in_region$exonEnds[i], ',')[[1]])
        exon_start = exon_start[!is.na(exon_start)]; exon_end = exon_end[!is.na(exon_end)]
        if (showExons == "Yes"){ for (j in 1:length(exon_start)){ fig2 = fig2 %>% add_trace(x=c(exon_start[j], exon_end[j], exon_end[j], exon_start[j]), y = c(genes_in_region$y[i]+y_space, genes_in_region$y[i]+y_space, genes_in_region$y[i]-y_space, genes_in_region$y[i]-y_space), type = 'scatter', mode = 'lines', fill = 'tozeroy', fillcolor = 'black', line = list(width = 1, color = 'black'), showlegend = F, hoverinfo="none") } }
      }
      fig2 = fig2 %>% layout(plot_bgcolor='#e5ecf6', yaxis = list(zeroline = F, showticklabels=FALSE, gridcolor = 'ffff', title = "Gene track", autorange = FALSE, range = c(min(genes_in_region$y)-1, 0)), xaxis = list(zerolinecolor = '#ffff', gridcolor = 'ffff', zerolinewidth = 4, title = "Genomic position (bp)", autorange = FALSE, range = c(pos_start, pos_end)))
    # PLOT 3 IS THE SV
      fig3 <- plot_ly()
      svs_in_region$labels = paste0("<b>Start position:<b> ", svs_in_region$start_pos, "<br><b>End position:<b> ", svs_in_region$end_pos, "<br><b>Max. allele size:<b> ", svs_in_region$diff_alleles, "<br><b>Type:<b> ", svs_in_region$type, "<br><b>Source:<b> ", svs_in_region$source, "<extra></extra>")
      if (nrow(svs_in_region) >0){ 
        tmp_svs = svs_in_region[!duplicated(svs_in_region$type),]; tmp_svs$fake_x = 0; tmp_svs$fake_y = 0 # get the single types of SVs plotted --> this is to make the legend
        fig3 = fig3 %>% add_trace(data=tmp_svs, x=~fake_x, y=~fake_y, name=~type, type='scatter', mode='markers', marker=list(color=~col, symbol = 'square'), hoverinfo="none")
        for (i in 1:nrow(svs_in_region)){
          if (svs_in_region$diff_alleles[i] >= 3000){
            fig3 = fig3 %>% add_trace(x=c(svs_in_region$start_pos[i], svs_in_region$end_pos_plus[i]), y=rep(svs_in_region$y[i], 2), line = list(color = svs_in_region$col[i], width = 8), type = 'scatter', mode = 'lines', showlegend = F, hovertemplate = svs_in_region$labels[i]) 
          } else {
            fig3 = fig3 %>% add_trace(x=c(svs_in_region$start_pos[i], svs_in_region$end_pos_plus[i]), y=rep(svs_in_region$y[i], 2), type = 'scatter', mode = 'markers', marker=list(color=svs_in_region$col[i], size=8, symbol = 'square'), showlegend = F, hovertemplate = svs_in_region$labels[i])
          }
        }
      }
      fig3 = fig3 %>% layout(plot_bgcolor='#e5ecf6', 
                      yaxis = list(zeroline = F, showticklabels=FALSE, gridcolor = 'ffff', title = "Structural variants", autorange = FALSE, range = c(min(svs_in_region$y)-1, 0)), 
                      xaxis = list(zeroline = F, gridcolor = 'ffff', title = "Genomic position (bp)", autorange = FALSE, range = c(pos_start, pos_end)))
    # COMBINE THE FIGURES
      fig_final <- subplot(fig1, fig2, fig3, nrows = 3, heights = c(0.5, 0.25, 0.25), shareX = TRUE, titleX = TRUE, titleY = TRUE, margin = 0.02)    
      fig_final = fig_final %>% layout(legend = list(x = 0, y = 1, orientation = 'h', font = list(size = 13, color = "#000")), margin = list(l = 75, r = 100, b = 75, t = 75, pad = 4))
      fig_final = fig_final %>% config(toImageButtonOptions = list(format = "png", filename = "snpXplorer_plot", width = 1280, height = 960, scale = 3))
    return(fig_final)
  }

  # function to plot points and recombination rates -- based on ggplot
  Plot <- function(reference_genome, region_of_interest, snps_data, snp_interest, recomb_data, significance, pos_start, pos_end, plot_type, recomb, genes_in_region, svs_in_region, showExons, ld){
    # PLOT 1 IS THE MAIN SNP-PLOT
      # put data together for the plot
      data_combined = rbindlist(snps_data)
      # change position and pvalue to be numbers
      data_combined$pos = as.numeric(data_combined$pos); data_combined$p = as.numeric(data_combined$p)
      # set colors
      group_colors = c(); for (i in unique(data_combined$col)){ group_colors = c(group_colors, i = i) }; names(group_colors) = unique(data_combined$col2)
      # plot based on plot-type
      if (plot_type == 'Scatter'){
        # optionally add rsid for the snps that are plotted
        #if (reference_genome == 'GRCh37 (hg19)'){ snp_group = merge(snp_group, rsid_region, by.x = 'pos', by.y = "POS", all.x = T) } else { snp_group = merge(snp_group, rsid_region, by.x = 'pos', by.y = "POS_HG38", all.x = T) }
        # main plot 1
        plot_snp = ggplot(data = data_combined, aes(x = pos, y = -log10(p), color = col)) + geom_point(size = 3) + ylab('-Log10(P-value)') + xlab('') + ylim(0, significance) + 
          ggtitle(paste0(region_of_interest$chrom, ':', region_of_interest$start, '-', region_of_interest$end)) + xlim(region_of_interest$start, region_of_interest$end) +
          scale_color_manual(values=group_colors, name = 'GWAS', labels = unique(data_combined$name))
          theme(legend.position = "top", axis.text.x = element_blank(), axis.ticks.x = element_blank())
      } else {
        data_combined$bin = NA; data_combined$bin_value = NA; intervals = seq(region_of_interest$start, region_of_interest$end, length.out = 15)
        for (x in 1:length(intervals)){
          for (t in unique(data_combined$name)){
            data_combined$bin[which(data_combined$pos >= intervals[x] & data_combined$pos < intervals[x+1] & data_combined$name == t)] = rep(intervals[x], nrow(data_combined[which(data_combined$pos >= intervals[x] & data_combined$pos < intervals[x+1] & data_combined$name == t),]))
            data_combined$bin_value[which(data_combined$pos >= intervals[x] & data_combined$pos < intervals[x+1] & data_combined$name == t)] = rep(min(data_combined$p[which(data_combined$pos >= intervals[x] & data_combined$pos < intervals[x+1] & data_combined$name == t)]), nrow(data_combined[which(data_combined$pos >= intervals[x] & data_combined$pos < intervals[x+1] & data_combined$name == t),]))
          }
        }
        data_combined = data_combined[!is.na(data_combined$bin),]
        plot_snp = ggplot() + stat_smooth(data = data_combined, aes(x = bin, y = -log10(bin_value), fill = col, alpha = 0.5), geom = "area", method = 'loess', se = T) + 
          scale_fill_manual(values=group_colors, name = 'GWAS', labels = unique(data_combined$name)) + ylab('-Log10(P-value)') + xlab('') + ylim(0, significance) + 
          ggtitle(paste0(region_of_interest$chrom, ':', region_of_interest$start, '-', region_of_interest$end)) + scale_alpha(guide = 'none') + xlim(region_of_interest$start, region_of_interest$end) +
          theme(legend.position = "top", axis.text.x = element_blank(), axis.ticks.x = element_blank())
      }
      # add recombination rates if they are requested -- first need to scale the recombination rates
      if (recomb == 'Yes'){
        recomb_data$combined_scaled = (significance * (recomb_data$combined)) / 100
        plot_snp = plot_snp + geom_line(data = recomb_data, aes(x = pos, y = combined_scaled), color = 'orange') + xlim(region_of_interest$start, region_of_interest$end) +
          scale_y_continuous('-Log10(P-value)', limits = c(0, significance), sec.axis = sec_axis(~.*(100/significance), name="Recombination rates")) + 
          theme(axis.title.y.right = element_text(color = "orange"), axis.text.y.right = element_text(colour = "orange"), axis.ticks.y.right = element_line(color = "orange"))
      }
      plot_snp = plot_snp + theme(plot.margin = margin(1, 1, 0, 1, "pt"), plot.title = element_text(size=22), axis.text=element_text(size=14), axis.title=element_text(size=16, face="bold")) + 
        guides(colour = guide_legend(override.aes = list(size=8)))
    # PLOT 2 IS THE GENE TRACK
      # create new variables for the genes (start position based on strand, middle point)
      genes_in_region$y_space = genes_in_region$y + 0.05*abs(min(genes_in_region$y))
      genes_in_region$middle = genes_in_region$txStart + (genes_in_region$txEnd - genes_in_region$txStart)/2
      genes_in_region$start_strand = ifelse(genes_in_region$strand == '+', genes_in_region$txStart, genes_in_region$txEnd)
      genes_in_region$end_strand = ifelse(genes_in_region$strand == '+', genes_in_region$txEnd, genes_in_region$txStart)
      # plot gene segment without text and exons
      plot_genes = ggplot() + geom_segment(data = genes_in_region, aes(x = start_strand, y = y, xend = end_strand, yend = y, colour = strand, size = 2)) + 
        xlim(region_of_interest$start, region_of_interest$end) + theme(legend.position = 'top', axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
        xlab('') + ylab('Gene track') + ylim(min(genes_in_region$y) - 0.5, max(genes_in_region$y) + 0.5) + scale_color_discrete(name = 'Gene Strand')
      # if there are not too many genes, show the names
      if (nrow(genes_in_region) <= 20){ txt_size = 6 } else if (nrow(genes_in_region) <= 30){ txt_size = 4 } else if (nrow(genes_in_region) <= 40){ txt_size = 2 } else if (nrow(genes_in_region) <= 50){ txt_size = 1 }
      if (nrow(genes_in_region) <= 50){ plot_genes = plot_genes + annotate(geom = 'text', x = genes_in_region$middle, y = genes_in_region$y_space, label = genes_in_region$Gene, fontface = 'italic', size = txt_size) }
      # if exons are requested, plot them
      if (showExons == 'Yes'){
        exons_df = data.frame(start = as.numeric(unlist(strsplit(genes_in_region$exonStarts, ','))), end = as.numeric(unlist(strsplit(genes_in_region$exonEnds, ','))), y_pos = rep(genes_in_region$y, genes_in_region$exonCount), strand = rep(genes_in_region$strand, genes_in_region$exonCount))
        exons_df$col = ifelse(exons_df$strand == '+', '#55BCC2', '#E87D72')
        plot_genes = plot_genes + annotate('segment', x = exons_df$start, y = exons_df$y_pos, xend = exons_df$end, yend = exons_df$y_pos, size = 8, colour = exons_df$col) + 
          scale_size(guide = 'none')
      }
      plot_genes = plot_genes + theme(plot.margin = margin(0, 1, 0, 1, "pt"), axis.text=element_text(size=14), axis.title=element_text(size=16, face="bold")) + scale_size(guide = 'none') + 
        guides(colour = guide_legend(override.aes = list(size=8)))
    # PLOT 3 IS THE SV
      # plot small svs as dots, big svs as segments
      small_svs = svs_in_region[which(as.numeric(svs_in_region$diff_alleles) < 3000),]; big_svs = svs_in_region[which(as.numeric(svs_in_region$diff_alleles) > 3000),]
      if (nrow(small_svs) >0){
        plot_sv = ggplot() + geom_point(data = small_svs, aes(x = as.numeric(start_pos), y = y, colour = type, size = 3), shape = 15) + xlim(region_of_interest$start, region_of_interest$end) + 
          theme(legend.position = "top", axis.text.y = element_blank(), axis.ticks.y = element_blank()) + scale_color_discrete(name = 'SV Type', labels = unique(svs_in_region$type)) + 
          ylab('Structural variation') + xlab('Genomic Position (bp)') + scale_size(guide = 'none')
      } else if (nrow(big_svs) >0){
        plot_sv = ggplot() + geom_segment(data = big_svs, aes(x = as.numeric(start_pos), y = y, xend = as.numeric(end_pos), yend = y, size = 3, color = type)) + xlim(region_of_interest$start, region_of_interest$end) + 
          theme(legend.position = "top", axis.text.y = element_blank(), axis.ticks.y = element_blank()) + scale_color_discrete(name = 'SV Type', labels = unique(svs_in_region$type)) + 
          ylab('Structural variation') + xlab('Genomic Position (bp)') + scale_size(guide = 'none')
      }
      if (nrow(big_svs) >0){
        plot_sv = plot_sv + geom_segment(data = big_svs, aes(x = as.numeric(start_pos), y = y, xend = as.numeric(end_pos), yend = y, size = 3, color = type)) 
        plot_sv = plot_sv + theme(plot.margin = margin(0, 1, 1, 1, "pt")) + scale_size(guide = 'none')
      }
      plot_sv = plot_sv + theme(plot.margin = margin(0, 1, 1, 1, "pt"), axis.text=element_text(size=14), axis.title=element_text(size=16, face="bold")) + scale_size(guide = 'none') +
        guides(colour = guide_legend(override.aes = list(size=8)))
    # COMBINE THE FIGURES
      combined <- (plot_snp / plot_genes / plot_sv) + plot_layout(heights = c(2, 1, 1), guides = "collect") & theme(legend.position = "top", legend.key.size = unit(1, 'cm'), legend.text = element_text(size=12), legend.title = element_text(size=14), legend.box="vertical")
        #fig_final = fig_final %>% config(toImageButtonOptions = list(format = "png", filename = "snpXplorer_plot", width = 1280, height = 960, scale = 3))
    return(combined)
  }

  # function to plot points and recombination rates -- based on ggplot
  Plot_My <- function(reference_genome, region_of_interest, snps_data, snp_interest, recomb_data, significance, pos_start, pos_end, plot_type, recomb, genes_in_region, svs_in_region, showExons, ld){
    # PLOT 1 IS THE MAIN SNP-PLOT
      # put data together for the plot
      data_combined = rbindlist(snps_data)
      # change position and pvalue to be numbers
      data_combined$pos = as.numeric(data_combined$pos); data_combined$p = as.numeric(data_combined$p)
      # set colors
      group_colors = c(); for (i in unique(data_combined$col)){ group_colors = c(group_colors, i = i) }; names(group_colors) = unique(data_combined$col2)
      # set layout for the plot
      par(mar = c(0, 5, 4, 5))
      layout(matrix(c(1,1,2,3), nrow = 4, ncol = 1, byrow = T))
      # basic plot 1
      plot(x = 0, y = 0, cex = 1.50, pch = 16, col = 'white', ylab = '-Log10(P-value)', xlab = '', ylim = c(0, significance), 
        main = paste0(region_of_interest$chrom, ':', region_of_interest$start, '-', region_of_interest$end), 
        xlim = c(region_of_interest$start, region_of_interest$end), xaxt = 'none', cex.lab = 2.25, cex.axis = 1.80, cex.main = 4)
      # rectange for color panel
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey95")
      # grid
      grid(nx = NULL, ny = NULL, lty = 1, col = "white", lwd = 1.5)
      # add recombination rates if they are requested -- first need to scale the recombination rates
      if (recomb == 'Yes'){
        recomb_data$combined_scaled = (significance * (recomb_data$combined)) / 100
        lines(recomb_data$pos, recomb_data$combined_scaled, col = 'orange', lwd = 1.25)
        recomb_axis <- significance * ((seq(0, 100, 25) - min(seq(0, 100, 25)))/(max(seq(0, 100, 25)) - min(seq(0, 100, 25))))
        axis(side = 4, at = recomb_axis, labels = seq(0, 100, 25), col= 'orange', col.axis = "orange", xpd=T, cex.axis = 1.80)
        text(x = max(region_of_interest$end), y = significance/3*2, "Recombination rate (cM/Mb)", srt = -90, col='orange', xpd=T, pos = 4, offset = 6, cex = 2.25)
      }
      # plot points
      if (plot_type == 'Scatter'){
        # points        
        points(x = data_combined$pos, y = -log10(data_combined$p), cex = 3, pch = 16, col = alpha(data_combined$col, 0.8))
        # ld in case it's requested
        if (!is.na(ld)){
          ld_df = merge(ld, data_combined, by.x = 'SNP_B', by.y = 'rsid')
          points(x = ld_df$pos, y = -log10(ld_df$p), cex = 3, pch = 16, col = alpha(ld_df$colo, 0.8))
        }
        # annotate top snp
        top = data_combined[which(data_combined$p == min(data_combined$p)),]
        text(x = top$pos, y = -log10(top$p), label = top$rsid, pos = 3, offset = 1, cex = 1.50)
        # legend
        legend('topright', legend = unique(data_combined$name), col = group_colors, pch = 16, ncol = length(group_colors), bg = 'grey95', pt.cex = 3, cex = 1.80)
      # else plot densities
      } else {
        # densities
        for (i in 1:length(snps_data)){
          toplot <- DensityLinePvalue(snp.info = snps_data[[i]], wind.n = 50, smooth.par = 0.2)
          polygon(x = toplot$x, y = toplot$y, col = alpha(group_colors[i], 0.6), lwd=3, xaxs="i", border = group_colors[i])
        }
        # legend
        legend('topright', legend = unique(data_combined$name), col = group_colors, lty = 1, ncol = length(group_colors), bg = 'grey95', lwd = 4, cex = 1.80)
      }
    # PLOT 2 IS THE GENE TRACK
      # parse gene-data
      genes_in_region$y_space = genes_in_region$y + 0.05*abs(min(genes_in_region$y))
      genes_in_region$middle = genes_in_region$txStart + (genes_in_region$txEnd - genes_in_region$txStart)/2
      genes_in_region$start_strand = ifelse(genes_in_region$strand == '+', genes_in_region$txStart, genes_in_region$txEnd)
      genes_in_region$end_strand = ifelse(genes_in_region$strand == '+', genes_in_region$txEnd, genes_in_region$txStart)
      genes_in_region$color = ifelse(genes_in_region$strand == '+', 'navy', 'coral')
      # basic plot 2
      par(mar = c(0, 5, 1, 5))
      plot(x = 0, y = 0, cex = 1.50, pch = 16, col = 'white', ylab = 'Gene track', xlab = '', ylim = c(min(genes_in_region$y) - 0.5, max(genes_in_region$y) + 0.5),
        xlim = c(region_of_interest$start, region_of_interest$end), xaxt = 'none', cex.lab = 2.25, yaxt = 'none')
      # rectange for color panel
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey95")
      # grid
      grid(nx = NULL, ny = NULL, lty = 1, col = "white", lwd = 1.5)
      # plot gene segment without text and exons
      segments(x0 = genes_in_region$start_strand, y0 = genes_in_region$y_space, x1 = genes_in_region$end_strand, y1 = genes_in_region$y_space, lwd = 10, col = genes_in_region$color)
      # if there are not too many genes, show the names
      if (nrow(genes_in_region) <= 10){ txt_size = 2.25 } else if (nrow(genes_in_region) <= 20){ txt_size = 1.75 } else if (nrow(genes_in_region) <= 30){ txt_size = 1.25 } else if (nrow(genes_in_region) <= 40){ txt_size = 0.75 } else if (nrow(genes_in_region) <= 50){ txt_size = 0.25 }
      if (nrow(genes_in_region) <= 50){ text(x = genes_in_region$middle, y = genes_in_region$y_space, label = genes_in_region$Gene, cex = txt_size, pos = 3, offset = 1.5) }
      # legend for strand
      legend('topright', legend = c('Forward strand', 'Reverse strand'), col = c('navy', 'coral'), lty = 1, ncol = 2, bg = 'grey95', lwd = 4, cex = 1.80)
      # if exons are requested, plot them
      if (showExons == 'Yes'){
        exons_df = data.frame(start = as.numeric(unlist(strsplit(genes_in_region$exonStarts, ','))), end = as.numeric(unlist(strsplit(genes_in_region$exonEnds, ','))), y_pos = rep(genes_in_region$y_space, genes_in_region$exonCount), strand = rep(genes_in_region$strand, genes_in_region$exonCount))
        exons_df$color = ifelse(exons_df$strand == '+', 'navy', 'coral')
        size = ((max(genes_in_region$y) + 0.5) - (min(genes_in_region$y) - 0.5)) * 0.025
        rect(xleft = exons_df$start, ybottom = exons_df$y_pos - size, xright = exons_df$end, ytop = exons_df$y_pos + size, col = 'white', border = 'black')
      }
    # PLOT 3 IS THE SV
      # basic plot 3
      par(mar = c(0, 5, 1, 5))
      plot(x = 0, y = 0, cex = 1.50, pch = 16, col = 'white', ylab = 'Genomic Position (bp)', xlab = '', ylim = c(min(genes_in_region$y) - 0.5, max(genes_in_region$y) + 0.5),
        xlim = c(region_of_interest$start, region_of_interest$end), cex.lab = 2.25, yaxt = 'none', cex.axis = 1.80)
      # rectange for color panel
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey95")
      # grid
      grid(nx = NULL, ny = NULL, lty = 1, col = "white", lwd = 1.5)
      # plot gene segment without text and exons
      svs_in_region$col = as.character(svs_in_region$col)
      svs_in_region$start_pos = as.numeric(svs_in_region$start_pos); svs_in_region = svs_in_region[!is.na(svs_in_region$start_pos),]; svs_in_region$end_pos = as.numeric(svs_in_region$end_pos);
      segments(x0 = as.numeric(svs_in_region$start_pos), y0 = svs_in_region$y, x1 = as.numeric(svs_in_region$end_pos), y1 = svs_in_region$y, lwd = 10, col = svs_in_region$col)
      # legend
      legend('topright', legend = unique(svs_in_region$type), col = unique(svs_in_region$col), lty = 1, ncol = length(unique(svs_in_region$type)), bg = 'grey95', lwd = 4, cex = 1.80)
  }
 
  # function to plot error
  PlotError <- function(title){
    p <- plotly_empty(type = "scatter", mode = "markers") %>% config(displayModeBar = FALSE) %>% layout( title = list(text = title, yref = "paper", y = 0.5))
    return(p)
  }

  # function to calculate density
  densityLinePvalue <- function(snps_data, wind.n, span_value){
    min.wind <- min(snps_data$pos)          # define minimum position
    max.wind <- max(snps_data$pos)          # define maximum position
    interval <- ceiling((max.wind - min.wind)/wind.n)   # define windows
    out = matrix(data = NA, nrow = wind.n, ncol = 3)    # define output
    colnames(out) <- c("window", "pvalue", "error")
    #counter for output assignement
    counter <- 1
    for (i in seq(min.wind, max.wind, interval)){
      if (counter <= wind.n){
        #define internal maximum -- f
        f <- i + interval
        #take info in the window of interest
        data.sbs <- snps_data[which((snps_data$pos >= i) & (snps_data$pos <= f)),]
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
    out[which(is.na(out[, "error"]) == TRUE), "error"] <- 0
    out <- as.data.frame(out)
    lo <- loess(out$pvalue ~ out$window, span=span_value)
    xl <- seq(min(out$window), max(out$window), (max(out$window) - min(out$window))/100)
    pred <- predict(lo, xl)
    pred[which(pred < 0)] <- 0
    #add limits -- left and right for polygon function
    xl <- c(xl[1], xl, xl[length(xl)])
    pred <- c(0, pred, 0)
    df = data.frame(x = xl, y = pred)
    return(df)
  }

## LIFTOVER FUNCTION
  # function to liftover data
  liftOver_data <- function(chrom, start, end, type, p = NULL, MAIN_PATH, from){
    if (is.null(p)){ df <- data.frame(chr=paste("chr", chrom, sep=""), start=start, end=end) } else { df <- data.frame(chr=paste("chr", chrom, sep=""), start=start, end=end, p = p) }
    df$group = seq(1, nrow(df))
    # change to GR class object
    gr <- makeGRangesFromDataFrame(df)
    # set chain file
    if (from == 'hg19') { chain <- import.chain(paste0(MAIN_PATH, "data/databases/hg19ToHg38.over.chain")) } else { chain <- import.chain(paste0(MAIN_PATH, "data/databases/hg38ToHg19.over.chain")) }
    # change coordinates
    gr_hg38 <- liftOver(gr, chain)
    # back to dataframe and clean it
    df_hg38 <- as.data.frame(gr_hg38)
    if (type == "interval"){
      chrom = as.character(chrom); start = as.numeric(df_hg38$start); end = as.numeric(df_hg38$end)
      results = list(chrom, start, end)
    } else {
      df_hg38 = df_hg38[, c("group", "start", "end")]
      df_hg38 = merge(df_hg38, df, by = "group")
      results = list(rep(chrom[1], nrow(df_hg38)), df_hg38$start.x, df_hg38$p)
    }
    return(results)
  }

## LD FUNCTIONS
  # function to manage whether LD should be done
  findLD <- function(ld_type, pop_interest, data_to_plot, rsid_region, reference_genome, region_of_interest, MAIN_PATH, snp_interest){
    pop_file = fread(paste0(MAIN_PATH, "data/databases/1000G/people.txt"), h=T, stringsAsFactors = F, sep="\t")     # read all populations
    # first check for "all" individuals -- in case no group is selected, by default uses european
    if (length(pop_interest) == 0){ pop_interest = c("ALL_eur") }
    if ("ALL_eur" %in% pop_interest){ pop_interest = c(pop_interest, "FIN", "GBR", "IBS", "TSI", "FIN", "CEU") }
    if ("ALL_amr" %in% pop_interest){ pop_interest = c(pop_interest, "CLM", "MXL", "PEL", "PUR") }
    if ("ALL_afr" %in% pop_interest){ pop_interest = c(pop_interest, "ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI") }
    if ("ALL_eas" %in% pop_interest){ pop_interest = c(pop_interest, "CDX", "CHB", "CHS", "JPT", "KHV") }
    if ("ALL_sas" %in% pop_interest){ pop_interest = c(pop_interest, "BEB", "GIH", "ITU", "PJL", "STU") }
    pop_interest = pop_interest[!duplicated(pop_interest)]; tmp_pop = pop_file[which(pop_file$`Population code` %in% pop_interest),]; tmp_df = data.frame("FID" = tmp_pop$`Sample name`, "IID" = tmp_pop$`Sample name`)
    write.table(tmp_df, "tmp_populations.txt", quote=F, row.names=F, col.names=T, sep="\t")   # write file with individuals of interest based on population
    # now choose which snp to find LD for
    if (ld_type == 'Input variant'){
      all_data = rbindlist(data_to_plot); target_snp = all_data[which(all_data$pos == snp_interest),]
      #if (reference_genome == 'GRCh37 (hg19)'){ target_snp = merge(target_snp, rsid_region, by.x = 'pos', by.y = 'POS') } else { target_snp = merge(target_snp, rsid_region, by.x = 'pos', by.y = 'POS_HG38') }
      target_snp = target_snp$rsid
    } else {
      all_data = rbindlist(data_to_plot); top_snps = head(all_data[order(all_data$p)], 100)     # extract top snps
      target_snp = top_snps[!is.na(top_snps$rsid),]; target_snp = target_snp$rsid[1]
      #if (reference_genome == 'GRCh37 (hg19)'){ top_snps = merge(top_snps, rsid_region, by.x = 'pos', by.y = 'POS') } else { top_snps = merge(top_snps, rsid_region, by.x = 'pos', by.y = 'POS_HG38') }
      #target_snp = top_snps[!is.na(top_snps$ID),]; target_snp = head(target_snp$ID[order(target_snp$p)], 1)
    }
    system(paste0("./plink --bfile ", MAIN_PATH, "data/databases/1000G/chr", region_of_interest$chrom, " --keep tmp_populations.txt --r2 --ld-snp ", target_snp, " --ld-window-kb 250000 --out tmp_ld"))
    ld <- fread("tmp_ld.ld", h=T)        # read ld file back
    # assign colors to the original data depending on LD (r2)
    ld$colo <- NA; ld$colo[which(abs(ld$R) >= 0.2)] <- "deepskyblue3"; ld$colo[which(abs(ld$R) >= 0.4)] <- "yellow"; ld$colo[which(abs(ld$R) >= 0.6)] <- "orange"; ld$colo[which(abs(ld$R) >= 0.8)] <- "red"
    # finally add genome hg38 position in case this is needed
    if (reference_genome == 'GRCh38 (hg38)'){ 
      ld = merge(ld, all_data, by.x = 'SNP_B', by.y = 'rsid'); ld = ld[, c('SNP_A', 'SNP_B', 'pos', 'R2', 'colo')]; colnames(ld) = c('SNP_A', 'SNP_B', 'BP_B', 'R2', 'colo')
    } else { 
      ld = ld[, c('SNP_A', 'SNP_B', 'BP_B', 'R2', 'colo')]
    }
    system(paste0("rm ", MAIN_PATH, "snpXplorer_v4/tmp_*"))
    return(ld)
  }

  findLD_topmed <- function(ld_type, data_to_plot, rsid_region, reference_genome, region_of_interest){}

## PLOT GTEX
  # function to plot gtex expression information
  plotGTEx <- function(gtex, color_palette, tissues_interest){
    if (color_palette == 'Default'){
      col1 = jcolors('pal9')[5]; col2 = jcolors('pal9')[1]; colfunc <- colorRampPalette(c(col1, col2)); colfunc = colfunc(256)
    } else if (color_palette == "Viridis"){
      colfunc = viridis(256)
    } else if (color_palette == 'Plasma'){
      colfunc = viridis(256, option = 'plasma')
    } else if (color_palette == 'Blues'){
      colfunc <- colorRampPalette(c('light blue', 'navy')); colfunc = colfunc(256)
    } else if (color_palette == 'Reds'){
      colfunc <- colorRampPalette(c('pink', 'dark red')); colfunc = colfunc(256)
    }
    colfunc_clusters <- colorRampPalette(c('purple', 'yellow', 'red')); col_clusters = colfunc_clusters(31)
    labels = gtex$Description
    tissues = data.frame(tissue_name = colnames(gtex), tissue_class = NA); tissues$tissue_class = str_split_fixed(tissues$tissue_name, " - ", 2)[, 1]
    if (tissues_interest[1] != 'All tissues'){ gtex = as.data.frame(gtex[, ..tissues_interest]); tissues = tissues[which(tissues$tissue_name %in% tissues_interest),] } else { gtex = gtex[, 3:ncol(gtex)] }
    rownames(gtex) = labels
    tmp_for_plot = as.matrix(t(gtex)); colnames(tmp_for_plot) = labels
    cb_grid <- setup_colorbar_grid(nrows = 1, y_start = 0.60, y_length = 0.9, y_spacing = 0.6)
    toolt_opt <- setup_tooltip_options(row = TRUE, col = TRUE, value = TRUE, prepend_row = "Tissue: ", prepend_col = "Gene: ", prepend_value = "Median TPM: ")
    if (ncol(tmp_for_plot) <=1){
        tmp = matrix(data=NA, nrow=1, ncol = 1)
        #main_heatmap(tmp) %>% add_col_title("Ops, too few genes. Please enlarge your region of interest!", side= "top")
        return(pheatmap(tmp))
    } else {
      return(pheatmap(tmp_for_plot))
      #fig = main_heatmap(tmp_for_plot, colors = colfunc, name = "Median TPM", colorbar_grid = cb_grid, tooltip = toolt_opt) %>% 
      #  add_row_clusters(title = "Tissue", name = 'Tissue', factor(tissues$tissue_class), colors = col_clusters) %>% add_col_clustering() %>% add_row_clustering() %>%
      #  add_row_title("Tissues") %>% add_col_labels()
    }
    #return(fig)
  }

## FIND EQTLS
  # function to extract eqtls
  extractEqtls <- function(region_of_interest, reference_genome, MAIN_PATH, tissues_interest){
    region_to_search = paste0(region_of_interest$chrom, ':', region_of_interest$start, '-', region_of_interest$end)
    print(region_to_search)
    if (reference_genome == 'GRCh37 (hg19)'){
      eqtls = tabix(region_to_search, paste0(MAIN_PATH, 'data/databases/eqtls_snpxplorer/chr', region_of_interest$chrom, '_summary_eqtls_hg37.txt.gz'), check.chr = F, verbose = FALSE)
      if (nrow(eqtls) >0){
        colnames(eqtls) = c('ensg', 'a1', 'a2', 'tissue', 'effect', 'p', 'gene', 'pos', 'chr')
      } else {
        eqtls = data.frame(ensg = as.character(), a1 = as.character(), a2 = as.character(), tissue = as.character(), effect = as.character(), p = as.character(), gene = as.character(), pos = as.character(), chr = as.character())
      }
    } else {
      eqtls = tabix(region_to_search, paste0(MAIN_PATH, 'data/databases/eqtls_snpxplorer/chr', region_of_interest$chrom, '_summary_eqtls_hg38.txt.gz'), check.chr = F, verbose = FALSE)
      if (nrow(eqtls) >0){
        colnames(eqtls) = c('pos', 'ensg', 'a1', 'a2', 'tissue', 'effect', 'p', 'gene', 'start_hg19', 'chr')
      } else {
        eqtls = data.frame(pos = as.character(), ensg = as.character(), a1 = as.character(), a2 = as.character(), tissue = as.character(), effect = as.character(), p = as.character(), gene = as.character(), start_hg19 = as.character(), chr = as.character())
      }
    }
    # restrict by tissue
    if (!('All tissues' %in% tissues_interest)){ 
      tissues_interest = str_replace_all(tissues_interest, '-', ''); tissues_interest = str_replace_all(tissues_interest, '  ', '_'); tissues_interest = str_replace_all(tissues_interest, ' ', '_')
      tissues_interest = str_replace_all(tissues_interest, '\\(', ''); tissues_interest = str_replace_all(tissues_interest, '\\)', '') 
      eqtls = eqtls[which(eqtls$tissue %in% tissues_interest),]
    } 
    # reduce and output
    eqtls = eqtls[, c('chr', 'pos', 'a1', 'a2', 'tissue', 'gene', 'effect', 'p')]
    colnames(eqtls) = c('Chr', 'Pos', 'A1', 'A2', 'Tissue', 'Gene', 'Effect', 'P')
    return(eqtls)
  }

## FIND SQTLS
  # function to extract sqtls
  extractSqtls <- function(region_of_interest, reference_genome, MAIN_PATH, tissues_interest){
    region_to_search = paste0(region_of_interest$chrom, ':', region_of_interest$start, '-', region_of_interest$end)
    if (reference_genome == 'GRCh37 (hg19)'){
      sqtls = tabix(region_to_search, paste0(MAIN_PATH, 'data/databases/summary_sqtls/chr', region_of_interest$chrom, '_summary_sqtls_hg37.txt.gz'), check.chr = F, verbose = FALSE)
      if (nrow(sqtls) >0){
        colnames(sqtls) = c('ensg', 'a1', 'a2', 'effect', 'p', 'tissue', 'gene', 'pos', 'chr')
      } else {
        sqtls = data.frame(ensg = as.character(), a1 = as.character(), a2 = as.character(), tissue = as.character(), effect = as.character(), p = as.character(), gene = as.character(), pos = as.character(), chr = as.character())
      }
    } else {
      sqtls = tabix(region_to_search, paste0(MAIN_PATH, 'data/databases/summary_sqtls/chr', region_of_interest$chrom, '_summary_sqtls_hg38.txt.gz'), check.chr = F, verbose = FALSE)
      if (nrow(sqtls) >0){
        colnames(sqtls) = c('pos', 'ensg', 'a1', 'a2', 'effect', 'p', 'tissue', 'gene', 'start_hg19', 'chr')
      } else {
        sqtls = data.frame(pos = as.character(), ensg = as.character(), a1 = as.character(), a2 = as.character(), tissue = as.character(), effect = as.character(), p = as.character(), gene = as.character(), start_hg19 = as.character(), chr = as.character())
      }
    }
    # restrict by tissue
    if (!('All tissues' %in% tissues_interest)){ 
      tissues_interest = str_replace_all(tissues_interest, '-', ''); tissues_interest = str_replace_all(tissues_interest, '  ', '_'); tissues_interest = str_replace_all(tissues_interest, ' ', '_')
      tissues_interest = str_replace_all(tissues_interest, '\\(', ''); tissues_interest = str_replace_all(tissues_interest, '\\)', '') 
      sqtls = eqtls[which(sqtls$tissue %in% tissues_interest),]
    } 
    # reduce and output
    sqtls = sqtls[, c('chr', 'pos', 'a1', 'a2', 'tissue', 'gene', 'effect', 'p')]
    colnames(sqtls) = c('Chr', 'Pos', 'A1', 'A2', 'Tissue', 'Gene', 'Effect', 'P')
    return(sqtls)
  }

## FIND GWASCATALOG INFO
  # function to extract Gwas cat info
  findGWAScat <- function(region_of_interest, reference_genome, MAIN_PATH, gwascat_type){
    region_to_search = paste0(region_of_interest$chrom, ':', region_of_interest$start, '-', region_of_interest$end)
    if (reference_genome == 'GRCh37 (hg19)'){
      gwascat = tabix(region_to_search, paste0(MAIN_PATH, 'data/databases/GWAS_catalog/Gwas_catalog_hg19.txt.gz'), check.chr = F, verbose = FALSE)
      colnames(gwascat) = c('chr', 'pos', 'pos_hg38', 'rsid', 'gene', 'trait', 'pubmed_id', 'p')
    } else {
      gwascat = tabix(region_to_search, paste0(MAIN_PATH, 'data/databases/GWAS_catalog/Gwas_catalog_hg38.txt.gz'), check.chr = F, verbose = FALSE)
      colnames(gwascat) = c('chr', 'pos_hg19', 'pos', 'rsid', 'gene', 'trait', 'pubmed_id', 'p')
    }
    if (nrow(gwascat) >0){
      tmp = gwascat[, c('chr', 'pos', 'rsid', 'p', 'gene', 'trait', 'pubmed_id')]
      colnames(tmp) = c('Chr', 'Pos', 'RsID', 'P', 'Gene', 'Trait', 'Pubmed')
      if (gwascat_type == 'SNPs'){
        tmp = tmp[order(as.numeric(tmp$P)),]
        tmp = tmp[!duplicated(tmp$Pos),]
      } else {
        tmp = tmp[order(as.numeric(tmp$P)),]
        tmp = tmp[!duplicated(tmp$Gene),]
      }
    } else {
      tmp = data.frame(Chr = as.character(), Pos = as.character(), RsID = as.character(), P = as.character(), Gene = as.character(), Trait = as.character(), Pubmed = as.character())
    }
    return(tmp)
  }

## IDENTIFY TARGET REGION FOR THE CROSS-REFERENCE LINK
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

## CHECK IF STRING IS EMAIL ADDRESS
  isValidEmail <- function(x) {
    grepl("\\<[A-Z0-9._%+-]+@[A-Z0-9.-]+\\.[A-Z]{2,}\\>", as.character(x), ignore.case=TRUE)
  }

## NOT USED
  # function to plot densities and recombination rates
  plot_density <- function(snps_data, recomb_data, significance, pos_start, pos_end){
    # initialize the plot
    fig <- plot_ly()
    # add data from first dataframe Df1
    for (snp_group in snps_data){
    }
    # add data from second dataframe Df2
    fig <- fig %>% add_trace(data=recomb_data, name="Recombination", x = ~pos, y = ~combined, type = 'scatter', mode = 'lines', yaxis = 'y2')
    # show figure
    fig <- fig %>% layout(plot_bgcolor='#e5ecf6', title = "Regional plots", 
                          yaxis2 = list(zerolinecolor = '#ffff', gridcolor = 'ffff', zerolinewidth = 4, title = "Recombination rates", autorange = FALSE, range = c(0, 100), overlaying = "y", side = "right"),
                          yaxis = list(zerolinecolor = '#ffff', gridcolor = 'ffff', zerolinewidth = 4, title = "-log10(P-value)", autorange = FALSE, range = c(0, significance)), 
                          xaxis = list(zerolinecolor = '#ffff', gridcolor = 'ffff', zerolinewidth = 4, title = "Genomic position (bp)", autorange = FALSE, range = c(pos_start, pos_end)))
    return(fig)
  }
