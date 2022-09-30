  ## PLOTLY FUNCTION FOR PLOT 2 (GTEX)
    output$plot2 <- renderIheatmap({
      # get genes in the region and relative GTEx information
      genes_interest = genes_region()
      color_palette = input$heat_cols; tissues_interest = input$tissues
      if (nrow(genes_interest) >0 && ('All tissues' %in% tissues_interest || length(tissues_interest) >=2)){
        # define tissues and genes to show
        tmp_gtex = gtex.db[which(gtex.db$Description %in% genes_interest$Gene), ]
        fig = plotGTEx(tmp_gtex, color_palette, tissues_interest)
      } else {
        tmp = matrix(data=NA, nrow=1, ncol = 1)
        main_heatmap(tmp) %>% add_col_title("Ops, too few genes of too few tissues selected!", side= "top")
      }
    })

  ## PLOT TABLE OF SNPS
    output$table <- renderTable({
      plotError_for_table = plotError_for_table()
      if (plotError_for_table == TRUE){ 
        df = data.frame("Rs ID" = as.character(), "Pos" = as.character(), "Ref" = as.character(), "Alt" = as.character(), "MAF" = as.character(), "Log(p)" = as.character(), "Study" = as.character()) 
      } else { 
        rsid_region = rsid_region_for_table()
        reference_gen = reference_genome()
        df = rbindlist(data_for_table()); df$p = as.numeric(df$p); df$p = -log10(df$p); top = head(df, 80); 
        if (reference_gen == 'GRCh37 (hg19)'){ top = merge(top, rsid_region, by.x = 'pos', by.y = 'POS', all.x = T) } else { top = merge(top, rsid_region, by.x = 'pos', by.y = 'POS_HG38', all.x = T) }
      }
      if (nrow(df) >0){
          top = top[, c('ID', 'pos', 'REF', 'ALT', 'ALT_FREQS', 'p', 'name')]
          top = top[order(-top$p),]
          colnames(top) = c('Rs ID', 'Pos', 'Ref', 'Alt', 'MAF', 'Log(p)', 'Study')
        } else {
          top = data.frame(RsID = "No", POS = 'SNPs', Ref = "in", Alt = "the", MAF = "region", p = "of", Study = "interest")
          colnames(top) = c('Rs ID', 'Pos', 'Ref', 'Alt', 'MAF', 'Log(p)', 'Study')
        }
        return(top)
    }, bordered = TRUE, striped = TRUE)
  
  ## PLOT TABLE OF GWAS-CATALOG
    output$gwascat_table <- renderTable({
      plotError_for_table = plotError_for_table()
      if (plotError_for_table == TRUE){
        tmp_gwascat = data.frame(Chr = as.character(), Pos = as.character(), P = as.character())
      } else {
        reference_gen = reference_genome()
        region_of_interest = plotted_region()
        gwascat_type = gwascat_mode()
        genes_in_region = genes_region()
        df = rbindlist(data_for_table()); df$p = as.numeric(df$p); df = df[order(as.numeric(df$p)),]; df$p = -log10(df$p)
        if (gwascat_type == 'SNPs'){
          tmp_gwascat = gwascat[[as.numeric(region_of_interest$chrom)]]
          if (reference_gen == 'GRCh37 (hg19)') { tmp_gwascat = tmp_gwascat[which(tmp_gwascat$'Pos (hg19)' %in% df$pos),]; tmp_gwascat = tmp_gwascat[order(as.numeric(tmp_gwascat$'P-value')),] } else { tmp_gwascat = tmp_gwascat[which(tmp_gwascat$'Pos (hg38)' %in% df$pos),]; tmp_gwascat = tmp_gwascat[order(as.numeric(tmp_gwascat$'P-value')),] }
          tmp_gwascat = tmp_gwascat[!duplicated(tmp_gwascat[c('SNP', 'Pubmed ID')]),]
        } else {
          tmp_gwascat = gwascat_genes[which(gwascat_genes$gene %in% genes_in_region$Gene),]
          colnames(tmp_gwascat) = c('Gene', 'Trait', 'Study', 'Journal')
        }
      }
      if (nrow(tmp_gwascat) >0){
        top = head(tmp_gwascat, 50)
        if (gwascat_type == 'SNPs'){ top = top[, c('SNP', 'Gene', 'Trait', 'Pubmed ID', 'P-value')] }
      } else {
        if (gwascat_type == 'SNPs'){ top = data.frame('SNP' = 'No', 'Gene' = 'entries', 'Trait' = 'in', 'Pubmed ID' = 'GWAS', 'P-value' = 'Catalog') } else { top = data.frame('Gene' = 'No entries', 'Trait' = 'in the', 'Study' = 'GWAS', 'Journal' = 'Catalog') } 
      }
      return(top)
    }, bordered = TRUE, striped = TRUE)

  ## PLOT TABLE OF STRUCTURAL VARIANTS
    output$sv_table <- renderTable({
      sv_info = svs_for_table()
      region_of_interest = plotted_region()
      if (nrow(sv_info) >0){
        sb = sv_info[which(sv_info$end_pos >= region_of_interest$start & sv_info$start_pos <= region_of_interest$end),]
        tmp = data.frame(Start = sb$start_pos, End = sb$end_pos, Type = sb$type, Size = sb$diff_alleles, Source = sb$source)
      } else {
        tmp = data.frame(Start = "No", End = "SV", Type = "in the", Size = "region of", Source = "interest!")
      }
      return(tmp)     
    }, bordered = TRUE, striped = TRUE)
  
  ## DOWNLOAD SV TABLE
    output$download_SVs <- downloadHandler(
      filename = function(){ reference_gen = reference_genome(); if (reference_gen == 'GRCh37 (hg19)'){ paste0("All_SVs_region_hg19.txt") } else { paste0("All_SVs_region_hg38.txt") } },
      content = function(file) {
        sv_info = svs_for_table()
        region_of_interest = plotted_region()
        if (nrow(sv_info) >0){
          sb = sv_info[which(sv_info$end_pos >= region_of_interest$start & sv_info$start_pos <= region_of_interest$end),]
          tmp = data.frame(Start = sb$start_pos, End = sb$end_pos, Type = sb$type, Size = sb$diff_alleles, Source = sb$source)
        } else {
          tmp = data.frame(Start = "No", End = "SV", Type = "in the", Size = "region of", Source = "interest!")
        }
        write.table(tmp, file, row.names = FALSE, quote=F, sep="\t")
      })
  
  ## DOWNLOAD GWAS CATALOG TABLE
    output$download_GWASCat <- downloadHandler(
      filename = function(){ gwascat_type = gwascat_mode(); reference_gen = reference_genome(); 
        if (reference_gen == 'GRCh37 (hg19)'){ genome = 'hg19' } else { genome = 'hg38' }
        if (gwascat_type == 'SNPs'){ paste0("SNPs_information_from_GWASCatalog_", genome, ".txt") } else { paste0("Genes_information_from_GWASCatalog_", genome, ".txt") } },
      content = function(file) {
        plotError_for_table = plotError_for_table()
        if (plotError_for_table == TRUE){
          tmp_gwascat = data.frame(Chr = as.character(), Pos = as.character(), P = as.character())
        } else {
          reference_gen = reference_genome()
          region_of_interest = plotted_region()
          gwascat_type = gwascat_mode()
          genes_in_region = genes_region()
          df = rbindlist(data_for_table()); df$p = as.numeric(df$p); df = df[order(as.numeric(df$p)),]; df$p = -log10(df$p)
          if (gwascat_type == 'SNPs'){
            tmp_gwascat = gwascat[[as.numeric(region_of_interest$chrom)]]
            if (reference_gen == 'GRCh37 (hg19)') { tmp_gwascat = tmp_gwascat[which(tmp_gwascat$'Pos (hg19)' %in% df$pos),]; tmp_gwascat = tmp_gwascat[order(as.numeric(tmp_gwascat$'P-value')),] } else { tmp_gwascat = tmp_gwascat[which(tmp_gwascat$'Pos (hg38)' %in% df$pos),]; tmp_gwascat = tmp_gwascat[order(as.numeric(tmp_gwascat$'P-value')),] }
            tmp_gwascat = tmp_gwascat[!duplicated(tmp_gwascat[c('SNP', 'Pubmed ID')]),]
          } else {
            tmp_gwascat = gwascat_genes[which(gwascat_genes$gene %in% genes_in_region$Gene),]
            colnames(tmp_gwascat) = c('Gene', 'Trait', 'Study', 'Journal')
          }
        }
        if (nrow(tmp_gwascat) >0){
          if (gwascat_type == 'SNPs'){ tmp_gwascat = tmp_gwascat[, c('SNP', 'Gene', 'Trait', 'Pubmed ID', 'P-value')] }
        } else {
          if (gwascat_type == 'SNPs'){ tmp_gwascat = data.frame('SNP' = 'No', 'Gene' = 'entries', 'Trait' = 'in', 'Pubmed ID' = 'GWAS', 'P-value' = 'Catalog') } else { tmp_gwascat = data.frame('Gene' = 'No entries', 'Trait' = 'in the', 'Study' = 'GWAS', 'Journal' = 'Catalog') } 
        }
        write.table(tmp_gwascat, file, row.names = FALSE, quote=F, sep = "\t")
      })

  ## DOWNLOAD SNPS TABLE
    output$download_SNPsTable <- downloadHandler(
      filename = function(){ reference_gen = reference_genome(); if (reference_gen == 'GRCh37 (hg19)'){ paste0("Top_SNPs_region_hg19.txt") } else { paste0("Top_SNPs_region_hg38.txt") } },
      content = function(file){
        plotError_for_table = plotError_for_table()
        if (plotError_for_table == TRUE){ 
          df = data.frame("Rs ID" = as.character(), "Pos" = as.character(), "Ref" = as.character(), "Alt" = as.character(), "MAF" = as.character(), "Log(p)" = as.character(), "Study" = as.character()) 
        } else { 
          rsid_region = rsid_region_for_table()
          reference_gen = reference_genome()
          df = rbindlist(data_for_table()); df$p = as.numeric(df$p); df = df[order(as.numeric(df$p)),]; df$p = -log10(df$p) 
          if (reference_gen == 'GRCh37 (hg19)'){ df = merge(df, rsid_region, by.x = 'pos', by.y = 'POS', all.x = T) } else { df = merge(df, rsid_region, by.x = 'pos', by.y = 'POS_HG38', all.x = T) }
        }
        if (nrow(df) >0){
            df = df[, c('ID', 'pos', 'REF', 'ALT', 'ALT_FREQS', 'p', 'name')]
            colnames(df) = c('Rs ID', 'Pos', 'Ref', 'Alt', 'MAF', 'Log(p)', 'Study')
          } else {
            df = data.frame(RsID = "No", POS = 'SNPs', Ref = "in", Alt = "the", MAF = "region", p = "of", Study = "interest")
            colnames(df) = c('Rs ID', 'Pos', 'Ref', 'Alt', 'MAF', 'Log(p)', 'Study')
          }
          write.table(df, file, row.names = F, quote = F, sep = "\t")
      })

  ## DOWNLOAD IMAGE OF THE PLOT
    output$download_plot <- downloadHandler(
      filename = function(){ format = figureFormat(); if (format == 'PDF'){ paste0('snpXplorer_plot.pdf') } else { paste('snpXplorer_plot.tar.gz') } },
      content = function(file){
        fig = mainFigure(); format = figureFormat()
        system('rm snpXplorer_plot.html'); system('rm -rf snpXplorer_plot_files/'); system('rm -rf snpXplorer_plot'); system('rm snpXplorer_plot.tar.gz')      # clean files before saving plot
        htmlwidgets::saveWidget(fig, 'snpXplorer_plot.html', selfcontained = F)      # first save plot in html
        if (format == 'PDF'){
          path = getwd(); system(paste0('cutycapt --min-width=1080 --min-height=720 --url=file://', path, '/snpXplorer_plot.html --out=', file))   # use tool to convert to pdf
        } else {
          system('mkdir snpXplorer_plot'); system('mv snpXplorer_plot_files snpXplorer_plot/'); system('mv snpXplorer_plot.html snpXplorer_plot/')
          system('tar -cvf snpXplorer_plot.tar.gz snpXplorer_plot'); system(paste0('mv snpXplorer_plot.tar.gz ', file))
        }
      })
  ## EQTLS TABLE
    output$eqtls_table <- renderTable({
      plotError_for_table = plotError_for_table()
      if (plotError_for_table == FALSE){
        eqtl_info = eqtls_for_table()
        rsid_region = rsid_region_for_table()
        reference_gen = reference_genome()
        if (reference_gen == 'GRCh37 (hg19)'){ top = merge(eqtl_info, rsid_region, by.x = 'pos_hg19', by.y = 'POS', all.x = T); top = top[, c('ID', 'V2', 'V3', 'V4', 'V6', 'V7', 'gene')] } else { top = merge(eqtl_info, rsid_region, by.x = 'V1', by.y = 'POS_HG38', all.x = T); top = top[, c('ID', 'V2', 'V3', 'V4', 'V6', 'V7', 'gene')] }
      } else {
        top = data.frame(ID = 'No', A1 = 'eQTLs', A2 = 'in', Tissue = 'the', Effect = 'region', P = 'of', Gene = 'interest')
      }
      if (nrow(top) >0){
        colnames(top) = c('ID', 'A1', 'A2', 'Tissue', 'Effect', 'P', 'Gene'); top = top[order(top$P),]
      }
      return(top)
    })
  ## DOWNLOAD EQTLS TABLE
    output$download_eQTLsTable <- downloadHandler(
      filename = function(){ reference_gen = reference_genome(); if (reference_gen == 'GRCh37 (hg19)'){ paste0("All_eQTL_region_hg19.txt") } else { paste0("All_eQTL_region_hg38.txt") } },
      content = function(file) {
        plotError_for_table = plotError_for_table()
        if (plotError_for_table == FALSE){
          eqtl_info = eqtls_for_table()
          rsid_region = rsid_region_for_table()
          reference_gen = reference_genome()
          if (reference_gen == 'GRCh37 (hg19)'){ top = merge(eqtl_info, rsid_region, by.x = 'pos_hg19', by.y = 'POS', all.x = T); top = top[, c('ID', 'pos_hg19', 'V2', 'V3', 'V4', 'V6', 'V7', 'gene')] } else { top = merge(eqtl_info, rsid_region, by.x = 'V1', by.y = 'POS_HG38', all.x = T); top = top[, c('ID', 'V1', 'V2', 'V3', 'V4', 'V6', 'V7', 'gene')] }
        } else {
          top = data.frame(ID = 'No', Position = '', A1 = 'eQTLs', A2 = 'in', Tissue = 'the', Effect = 'region', P = 'of', Gene = 'interest')
        }
        if (nrow(top) >0){ colnames(top) = c('ID', 'Position', 'A1', 'A2', 'Tissue', 'Effect', 'P', 'Gene'); top = top[order(top$P),] }
        write.table(top, file, row.names = FALSE, quote=F, sep="\t")
      })

  ## SQTLS TABLE
    output$sqtls_table <- renderTable({
      plotError_for_table = plotError_for_table()
      if (plotError_for_table == FALSE){
        sqtl_info = sqtls_for_table()
        rsid_region = rsid_region_for_table()
        reference_gen = reference_genome()
        if (reference_gen == 'GRCh37 (hg19)'){ top = merge(sqtl_info, rsid_region, by.x = 'pos_hg19', by.y = 'POS', all.x = T); top = top[, c('ID', 'V2', 'V3', 'V7', 'V5', 'V6', 'gene')] } else { top = merge(sqtl_info, rsid_region, by.x = 'V1', by.y = 'POS_HG38', all.x = T); top = top[, c('ID', 'V2', 'V3', 'V7', 'V5', 'V6', 'gene')] }
      } else {
        top = data.frame(ID = 'No', A1 = 'sQTLs', A2 = 'in', Tissue = 'the', Effect = 'region', P = 'of', Gene = 'interest')
      }
      if (nrow(top) >0){
        colnames(top) = c('ID', 'A1', 'A2', 'Tissue', 'Effect', 'P', 'Gene'); top = top[order(top$P),]
      }
      return(top)
    })
  ## DOWNLOAD SQTLS TABLE
    output$download_sQTLsTable <- downloadHandler(
      filename = function(){ reference_gen = reference_genome(); if (reference_gen == 'GRCh37 (hg19)'){ paste0("All_sQTL_region_hg19.txt") } else { paste0("All_sQTL_region_hg38.txt") } },
      content = function(file) {
        plotError_for_table = plotError_for_table()
        if (plotError_for_table == FALSE){
          sqtl = sqtls_for_table()
          rsid_region = rsid_region_for_table()
          reference_gen = reference_genome()
          if (reference_gen == 'GRCh37 (hg19)'){ top = merge(sqtl, rsid_region, by.x = 'pos_hg19', by.y = 'POS', all.x = T); top = top[, c('ID', 'pos_hg19', 'V2', 'V3', 'V7', 'V5', 'V6', 'gene')] } else { top = merge(sqtl, rsid_region, by.x = 'V1', by.y = 'POS_HG38', all.x = T); top = top[, c('ID', 'V1', 'V2', 'V3', 'V7', 'V5', 'V6', 'gene')] }
        } else {
          top = data.frame(ID = 'No', Position = '', A1 = 'sQTLs', A2 = 'in', Tissue = 'the', Effect = 'region', P = 'of', Gene = 'interest')
        }
        if (nrow(top) >0){ colnames(top) = c('ID', 'Position', 'A1', 'A2', 'Tissue', 'Effect', 'P', 'Gene'); top = top[order(top$P),] }
        write.table(top, file, row.names = FALSE, quote=F, sep="\t")
      })
  ## BUG REPORT SECTION
    observeEvent(input$send_bug, {
      if (input$send_bug != 0){
        gwas_showed = gwas_to_show(); reference_gen = reference_genome(); plotError = plotError_for_table; user_message = as.character(input$bug_comments)
         ld_pop = ld_populations(); if (is.null(ld_pop)){ ld_pop = "Not requested" }
        # Find information about the current run
          toreport = data.frame("Input datasets" = paste0(gwas_showed, collapse = ","), "Reference" = reference_gen, "Browsing option" = input$browsing_options, "Target" = input$target, "Window-scale" = input$x_axis, "P-scale" = input$y_axis,
                              "Plot" = input$plot_type, "LD" = input$linkage_type, "LD populations" = ld_pop, 'User message' = user_message)
          print(toreport)
          colnames(toreport) = c("Input datasets", "Reference", "Browsing option", "Target", "Window-scale", "P-scale", "Plot", "LD", "LD populations", "User message")
        # Write table
          fname = paste0(MAIN_PATH, "bug_report/Bug_report_", str_replace_all(Sys.time(), " ", "_"), ".txt")
          write.table(toreport, fname, quote=F, row.names=F, sep = "\t")
        # Show alert
          show_alert( title = "Bug submitted!", text = "Thank you for your feedback!", type = "success")        
        # Finally send an email to me
          system(paste0("sendEmail -f n.tesi@amsterdamumc.nl -t snpxplorer@gmail.com -u 'snpXplorer bug report' -m 'Hi Nicco, \nA new bug has just been added. Go and check it out! \n snpXplorer team.' -cc n.tesi@amsterdamumc.nl snpxplorer@gmail.com -S /usr/sbin/sendmail"))
      }
    })
  ## GENE-CARDS LINK     
    output$genecards_link <- renderUI({
      if (input$browsing_options == "Locus, Gene, RsID"){
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
  ## GWAS-CATALOG LINK
    output$gwascat_link <- renderUI({
      if (input$browsing_options == "Locus, Gene, RsID"){
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
        } else {
          link <- paste0("https://www.ebi.ac.uk/gwas")
          url <- a("Click here!", href = link)
          tagList("GWAS-catalog link:", url)
        }
      }
    })
  ## LD-HUB LINK
    output$ld_hub_link <- renderUI({
      link <- paste0("http://ldsc.broadinstitute.org")
      url <- a("Click here!", href = link)
      tagList("LD Hub link:", url)
    })
  ## TABLE OF STUDIES INCLUDED (DOCUMENTATION TAB)
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
  ## TABLE OF THE 1000GENOME PROJECT (DOCUMENTATION TAB)
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
  ## TABLE OF STRUCTURAL VARIANTS DATASETS
    output$table_info_SVs <- renderTable({
      studies_table <- as.data.frame(matrix(data=NA, nrow=3, ncol = 3))
      colnames(studies_table) <- c("Study", "Experimental technology", "Reference")
      studies_table[1, ] <- c("Linthorst et al., 2020", "PacBio CLR", "https://doi.org/10.1038/s41398-020-01060-5")
      studies_table[2, ] <- c("Chaisson et al., 2019", "PacBio CLR + Illumina WGS + 10X Genomics + Hi-C", "https://doi.org/10.1038/s41467-018-08148-z")
      studies_table[3, ] <- c("Audano et al., 2019", "PacBio CLR", "https://doi.org/10.1016/j.cell.2018.12.019")
      return(studies_table)
    }, width = "100%")
  ## LINK TO NAR PAPER
    output$biorxiv_link <- renderUI({
      link <- paste0("https://academic.oup.com/nar/article/49/W1/W603/6287842")
      url <- a("Click here!", href = link)
      tagList("Link to paper:", url)
    })
  ## LINK TO LONGEVITY PAPER
    output$longevity_ms <- renderUI({
      link <- paste0("https://academic.oup.com/biomedgerontology/article/76/5/750/5996044")
      url <- a("Click here!", href = link)
      tagList("Link to the paper:", url)
    })
  ## LINK TO ROTATION PAPER
    output$rotation_ms <- renderUI({
      link <- paste0("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8724252/")
      url <- a("Click here!", href = link)
      tagList("Link to the paper:", url)
    })
  ## PDF VIEWER OF DOCUMENTATION
    output$pdf_doc_view <- renderUI({ tags$iframe(style="height:1000px; width:95%", src="snpXplorer_documentation.pdf") })
  ## PDF VIEWER OF NAR PAPER
    output$pdf_biorxiv_view <- renderUI({ tags$iframe(style="height:1000px; width:95%", src="gkab410-3.pdf") })
  ## DOWNLOAD SAMPLE DATA 1
    help_explo_1 <- fread("www/trial_dataset_website.txt.gz")
    output$download_help_exploration <- downloadHandler( filename = function() { paste0("trial_dataset_snpXplorer_exploration.txt") },
      content = function(file) { write.table(help_explo_1, file, row.names = FALSE, quote=F, sep=" ") })
  ## DOWNLOAD SAMPLE DATA 2
    help_explo_2 <- fread("www/trial_plink_snpxplorer.glm.logistic.gz")
    output$download_help_exploration_plink <- downloadHandler( filename = function() { paste0("trial_dataset_snpXplorer_exploration_PLINK.txt") },
      content = function(file) { write.table(help_explo_2, file, row.names = FALSE, quote=F, sep="\t") })
  ## DOWNLOAD SAMPLE DATA 3
    help_annot_1 <- fread("www/trial_annotation_f1.txt", h=F)
    output$download_help_annotation1 <- downloadHandler( filename = function() { paste0("trial_dataset_snpXplorer_annotation_type1.txt") },
      content = function(file) { write.table(help_annot_1, file, row.names = FALSE, col.names=F, quote=F) })
  ## DOWNLOAD SAMPLE DATA 4
    help_annot_2 <- fread("www/trial_annotation_f2.txt", h=F)
    output$download_help_annotation2 <- downloadHandler( filename = function() { paste0("trial_dataset_snpXplorer_annotation_type2.txt") },
      content = function(file) { write.table(help_annot_2, file, row.names = FALSE, col.names=F, quote=F) })
  ## DOWNLOAD SAMPLE DATA 5
    help_annot_3 <- fread("www/trial_annotation_f3.txt", h=F)
    output$download_help_annotation3 <- downloadHandler( filename = function() { paste0("trial_dataset_snpXplorer_annotation_type3.txt") },
      content = function(file) { write.table(help_annot_3, file, row.names = FALSE, col.names=F, quote=F) })
  ## LINK TO GITHUB PAGE
    url <- a("Visit our Github page", href="https://github.com/TesiNicco/SNPbrowser"); output$github <- renderUI({ tagList(url) })
  ## LINK TO YOUTUBE VIDEOS
    output$video1 <- renderUI({ HTML(paste0('<iframe width="560" height="315" src="https://www.youtube.com/embed/z8kqzTacaS4" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>')) })
    output$video2 <- renderUI({ HTML(paste0('<iframe width="560" height="315" src="https://www.youtube.com/embed/Ai3F-JBQL3U" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>')) })
    output$video3 <- renderUI({ HTML(paste0('<iframe width="560" height="315" src="https://www.youtube.com/embed/eOJTf7tk2Rg" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>')) })
    output$video4 <- renderUI({ HTML(paste0('<iframe width="560" height="315" src="https://www.youtube.com/embed/QP-5XjIEYpI" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>')) })
