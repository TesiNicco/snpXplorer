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
      if (ref_version == 'GRCh38 (hg38)'){ data_lifted = liftOver_data(chrom = as.character(chrom), start = pos1, end = pos2, type = "interval", MAIN_PATH); chrom = data_lifted[[1]]; pos1 = data_lifted[[2]]; pos2 = data_lifted[[3]] }
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
  findRsID <- function(region_of_interest, reference_genome, span, MAIN_PATH){
    tmp = fread(paste0(MAIN_PATH, 'data/databases/snps_info/chr', region_of_interest$chrom, '_snps_info.txt.gz'), h=T, stringsAsFactors=F)
    if (reference_genome == 'GRCh37 (hg19)'){ 
      tmp = tmp[which((tmp$POS >= region_of_interest$start - span) & (tmp$POS <= region_of_interest$end + span)),] 
    } else {
      tmp = tmp[which((tmp$POS_HG38 >= region_of_interest$start - span) & (tmp$POS_HG38 <= region_of_interest$end + span)),] 
    }
    return(tmp)
  }

## EXTRACT GENES FROM A REGION
  findGenes <- function(region_of_interest, reference_genome, genes_hg19, genes_hg38){
    span = 100000
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
  findSVs <- function(region_of_interest, reference_genome, all_str, all_str_hg38, span_value, sv_source){
    if (reference_genome == 'GRCh37 (hg19)'){ tmp = all_str } else { tmp = all_str_hg38; tmp$end_pos = as.numeric(tmp$end_pos) }
    svs = tmp[which(tmp$chr == paste0('chr', region_of_interest$chrom)),]
    svs = svs[which(svs$start_pos >= (region_of_interest$start - span_value) & svs$end_pos <= (region_of_interest$end + span_value)),]
    if (sv_source != 'all'){ svs = svs[which(svs$source == sv_source),] }
    if (nrow(svs) >0){
      svs$end_pos_plus = svs$end_pos + svs$diff_alleles
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
    data_to_plot = list()                 # initialize the list that will contain data to plot
    main_path = paste0(MAIN_PATH, 'data/')
    if (is.null(gwas_to_plot)){           # this is the case when no input gwas are selected --> example data
      if (region_of_interest$chrom %in% c(16, 17, 18, 19, 20, 21)){
        tmp = res_example[[(as.numeric(region_of_interest$chrom) - 15)]]
        colnames(tmp) = c("chrom", "pos", "p")
        # in case, we need to liftover
        if (reference_genome == 'GRCh38 (hg38)'){ data_lifted = liftOver_data(chrom = tmp$chrom, start = tmp$pos, end = tmp$pos + 1, type = "gwas", p = tmp$p); chrom_lf = data_lifted[[1]]; pos_lf = data_lifted[[2]]; p_lf = data_lifted[[3]]; tmp = data.table(chrom = chrom_lf, pos = pos_lf, p = p_lf) }
        # then match based on the position
        tmp = tmp[which((tmp$pos >= region_of_interest$start - span) & (tmp$pos <= region_of_interest$end + span)),]       
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
          tmp = fread(paste0(main_path, gwas_to_plot[i], '/chr', region_of_interest$chrom, '_', gwas_to_plot[i], '.txt.gz'), h=T, stringsAsFactors=F)
          colnames(tmp) = c('chrom', 'pos', 'p')
          # in case, we need to liftover
          if (reference_genome == 'GRCh38 (hg38)'){ data_lifted = liftOver_data(chrom = tmp$chrom, start = tmp$pos, end = tmp$pos + 1, type = "gwas", p = tmp$p, MAIN_PATH); chrom_lf = data_lifted[[1]]; pos_lf = data_lifted[[2]]; p_lf = data_lifted[[3]]; tmp = data.table(chrom = chrom_lf, pos = pos_lf, p = p_lf) }
          plotError = FALSE
          tmp$name = str_replace_all(gwas_to_plot[i], '_', ' ')
        }
        tmp = tmp[which((tmp$pos >= region_of_interest$start - span) & (tmp$pos <= region_of_interest$end + span)),]
        tmp$col = colors[i]; tmp$pos = as.numeric(tmp$pos); 
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
  extractRecombination = function(region_of_interest, reference_genome, genetic_map, genetic_map_hg38, span){
    if (reference_genome == 'GRCh37 (hg19)'){ tmp = genetic_map } else { tmp = genetic_map_hg38 }
    tmp = tmp[[as.numeric(region_of_interest$chrom)]]
    colnames(tmp) = c("chr", "pos", "combined", "cm")
    tmp$pos = as.numeric(tmp$pos); tmp$combined = as.numeric(tmp$combined)
    tmp = tmp[which((tmp$pos >= region_of_interest$start - span) & (tmp$pos <= region_of_interest$end + span)),]
    return(tmp)
  }

## PLOTS
  # function to plot points and recombination rates
  Plot <- function(reference_genome, region_of_interest, rsid_region, snps_data, snp_interest, recomb_data, significance, pos_start, pos_end, plot_type, recomb, genes_in_region, svs_in_region, showExons, ld){
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
          if (!is.null(ld)){
            ld = merge(ld, snp_group, by.x = 'BP_B', by.y = 'pos')
            ld$labels = paste0("<b>ID:<b> ", ld$ID, "<br><b>Position:<b> %{x} <br><b>Log(p):<b> %{y} <br><b>REF:<b> ", ld$REF, "<br><b>ALT:<b> ", ld$ALT, "<br><b>MAF:<b> ", ld$ALT_FREQS, "<br><b>LD with:<b> ", ld$SNP_A, "<br><b>R2:<b> ", ld$R2)
            if (nrow(ld) >0){ fig1 = fig1 %>% add_trace(data=ld, name='Linkage disequilibrium', x=~BP_B, y=~-log10(p), type='scatter', mode='markers', yaxis='y1', marker=list(color=~colo, size=10, symbol = 'diamond'), hovertemplate = ~labels) }
          }
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
      fig_final <- subplot(fig1, fig2, fig3, nrows = 3, heights = c(0.6, 0.2, 0.2), shareX = TRUE, titleX = TRUE, titleY = TRUE, margin = 0.02)    
      fig_final = fig_final %>% layout(legend = list(x = 0.25, y = 1, orientation = 'h', font = list(size = 13, color = "#000")), margin = list(l = 75, r = 100, b = 75, t = 75, pad = 4))
      fig_final = fig_final %>% config(toImageButtonOptions = list(format = "png", filename = "snpXplorer_plot", width = 1280, height = 960, scale = 3))
    return(fig_final)
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
  liftOver_data <- function(chrom, start, end, type, p = NULL, MAIN_PATH){
    df <- data.frame(chr=paste("chr", chrom, sep=""), start=start, end=end, p = p)
    df$group = seq(1, nrow(df))
    # change to GR class object
    gr <- makeGRangesFromDataFrame(df)
    # set chain file
    chain <- import.chain(paste0(MAIN_PATH, "data/databases/hg19ToHg38.over.chain"))
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
  findLD <- function(ld_type, pop_interest, data_to_plot, rsid_region, reference_genome, region_of_interest, MAIN_PATH){
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
      # to implement
    } else {
      all_data = rbindlist(data_to_plot); top_snps = head(all_data[order(all_data$p)], 100)     # extract top snps
      if (reference_genome == 'GRCh37 (hg19)'){ top_snps = merge(top_snps, rsid_region, by.x = 'pos', by.y = 'POS') } else { top_snps = merge(top_snps, rsid_region, by.x = 'pos', by.y = 'POS_HG38') }
      target_snp = top_snps[!is.na(top_snps$ID),]; target_snp = head(target_snp$ID[order(target_snp$p)], 1)
    }
    system(paste0("./plink --bfile ", MAIN_PATH, "data/databases/1000G/chr", region_of_interest$chrom, " --keep tmp_populations.txt --r2 --ld-snp ", target_snp, " --ld-window-kb 250000 --out tmp_ld"))
    ld <- fread("tmp_ld.ld", h=T)        # read ld file back
    # assign colors to the original data depending on LD (r2)
    ld$colo <- NA; ld$colo[which(abs(ld$R) >= 0.2)] <- "deepskyblue3"; ld$colo[which(abs(ld$R) >= 0.4)] <- "yellow"; ld$colo[which(abs(ld$R) >= 0.6)] <- "orange"; ld$colo[which(abs(ld$R) >= 0.8)] <- "red"
    # finally add genome hg38 position in case this is needed
    if (reference_genome == 'GRCh38 (hg38)'){ 
      ld = merge(ld, rsid_region, by.x = 'SNP_B', by.y = 'ID'); ld = ld[, c('SNP_A', 'SNP_B', 'POS_HG38', 'R2', 'colo')]; colnames(ld) = c('SNP_A', 'SNP_B', 'BP_B', 'R2', 'colo')
    } else { 
      ld = ld[, c('SNP_A', 'SNP_B', 'BP_B', 'R2', 'colo')]
    }
    system("rm tmp_.*")
    return(ld)
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
