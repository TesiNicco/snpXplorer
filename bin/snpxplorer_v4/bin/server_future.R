## LIBRARIES
  suppressPackageStartupMessages({
    library(shiny)
    library(bedr)
    library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
    library(pheatmap)
    library(patchwork)
    library(plotly)
    library(jcolors)
    library(viridis)
    library(iheatmapr)
    library(data.table)
    library(stringr)
    library(future)
  })

## ANNOTATIONS AND FUNCTIONS
  #MAIN_PATH = '/Users/nicco/Documents/GitHub/snpXplorer/bin/snpxplorer_v4/'
  #source(paste0(MAIN_PATH, "bin/functions.R"))
  MAIN_PATH = '/root/snpXplorer/'
  source(paste0(MAIN_PATH, "snpXplorer_v4/functions.R"))
  load(paste0(MAIN_PATH, "data/databases/20220803_annotationFiles.RData"))
  #load(paste0(MAIN_PATH, "data/databases/snps_info/chrAll_snps_info.RData"))

## DEFINE FUTURE CLASS
  plan(multicore)

## MAIN SERVER FUNCTION
server <- function(input, output) {  
  ## REACTIVE OBJECTS
    data_for_table = reactiveVal()
    plotError_for_table = reactiveVal()
    rsid_region_for_table = reactiveVal()
    reference_genome = reactiveVal()
    plotted_region = reactiveVal()
    gwascat_mode = reactiveVal()
    genes_region = reactiveVal()
    svs_for_table = reactiveVal()
    mainFigure = reactiveVal()
    figureFormat = reactiveVal()
    eqtls_for_table = reactiveVal()
    sqtls_for_table = reactiveVal()
    gwas_to_show = reactiveVal()
    gwascat_res = reactiveVal()
    ld_populations = reactiveVal()

  ## REACTIVE OBJECT THAT OBSERVES THE ANNOTATION SECTION
    observeEvent(input$run_annotation, {
      if (input$run_annotation != 0){
        message('** Annotation run was requested!\n')
        # Sample a number for randomization
          random_num <- sample(x = seq(1, 100000), size = 1, replace = F)
        # Take snps and save them
          annotateMe.snplist <- unlist(strsplit(input$snp_list, "\n"))
        # Define log file
          log_filename = paste0("annotateMe_run_", random_num, ".log")
        # Take other arguments
          ftype <- as.character(input$snp_list_type)
          if (ftype == "chr:pos (1:12345678)"){ ftype = 1 } else if (ftype == "chr pos (1 12345678)"){ ftype = 2 } else if (ftype == "rsid (rs12345)"){ ftype = 3 }
          ref_version = as.character(input$snp_list_reference)
          ref_version = str_split_fixed(ref_version, " ", 2)[, 1]
          username <- as.character(input$email)
          analysis_type = as.character(input$analysis_type)
          if (analysis_type == "gsea"){ analysis_type = "enrichment" } else { analysis_type = "mapping" }
          if (analysis_type == "gsea"){ analysis_mode = paste0(as.character(input$analysis_mode), collapse = ",") } else { analysis_mode = "default" }
          gtex_tissues = paste0(input$gtex_type, collapse = ",")
        # Then run annotate me externally in background -- this depends on the analysis_type requested
          annotateMe.cmd <- paste0("Rscript /root/snpXplorer/AnnotateMe/BIN/MAIN.R annotateMe_input_", random_num, ".txt ", ftype, " ", username, " ", analysis_type, " ", analysis_mode, " ", gtex_tissues, " ", ref_version, " > ", log_filename)
          message(paste0('** ', annotateMe.cmd))
      }
      if (username == 'Type your email address...'){
        show_alert(title = "Please insert your e-mail address!", text = "E-mail address is a mandatory field!", type = "error")
      } else {
        write.table(annotateMe.snplist, paste("/root/snpXplorer/snpXplorer_v3/annotateMe_input_", random_num, ".txt", sep=""), quote=F, row.names=F, col.names = F)
        system(annotateMe.cmd, ignore.stdout = F, wait = F)
        show_alert( title = "Analysis submitted!", text = "You should get a confirmation e-mail soon!", type = "success")        
      }
    })

  ## REACTIVE OBJECT THAT AVOIDS GREY OUTTING
    autoInvalidate <- reactiveTimer(10000)
    observe({ autoInvalidate(); cat("."); gc() })

  ## PLOTLY FUNCTION FOR PLOT 1 (SNPS, GENES AND SVS)
    output$plot <- renderPlot({
      ## DEFINE THE REACTIVE VARIABLES
        browsing_options = input$browsing_options; target = input$target; x_axis = input$x_axis; significance = input$y_axis
        reference_genome = input$reference_gen; sv_source = input$strVar_inp; all_colors = c(input$col1, input$col2, input$col3, input$col4, input$col5)
        ld_type = input$linkage_type; ld_pop = c(input$pop_eur, input$pop_afr, input$pop_amr, input$pop_eas, input$pop_sas)
        tissues_interest = input$tissues; ld_populations(ld_pop); gwascat_type = input$gwascat_type;

      ## FIND POSITIONS TO PLOT
        if (browsing_options == 'Locus, Gene, RsID'){ results = directInput(target, x_axis, gene.db, genes.hg38, reference_genome, SNPlocs.Hsapiens.dbSNP144.GRCh37, MAIN_PATH); region_of_interest = results[[1]]; snp_interest = results[[2]] } else { region_of_interest = data.frame(chrom = input$manual_chrom, start = input$manual_pos - input$x_axis, end = input$manual_pos + input$x_axis); snp_interest = NA }
        if (is.na(region_of_interest$start)){ plotError = TRUE } else { plotError = FALSE }
        message('** positions of interest FOUND!')

      ## FIND EQTL
        if (plotError == FALSE) { future_eqtl = future({ extractEqtls(region_of_interest, reference_genome, MAIN_PATH, tissues_interest) }) }
      
      ## FIND SQTL
        if (plotError == FALSE) { future_sqtl = future({ extractSqtls(region_of_interest, reference_genome, MAIN_PATH, tissues_interest) }) }
        
      ## FIND DATA TO PLOT IN CASE THE INPUT BROWSING OPTION WAS VALID
        if ((plotError == FALSE) & (is.null(input$upload_file) == FALSE) & (length(input$owned_gwas) >0)) { uploaded_file = input$upload_file$datapath } else { uploaded_file = NULL }
        gwas_to_plot = c(input$neurol_gwas, input$cardio_gwas, input$immune_gwas, input$cancer_gwas, input$physio_gwas, uploaded_file)
        gwas_to_show(gwas_to_plot)

      ## EXTRACT DATA FROM DATA TO PLOT IN CASE THE INPUT BROWSING OPTION WAS VALID
        if (plotError == FALSE) { future_res = future({ extractDataForPlot(gwas_to_plot, region_of_interest, span = 1000, colors = all_colors, res_example, reference_genome, MAIN_PATH) }) } else { data_to_plot = NA }

      ## FIND RECOMBINATION RATES
        if (plotError == FALSE) { future_rec = future({ extractRecombination(region_of_interest, reference_genome, span = 1000, MAIN_PATH) }) } 

      ## FIND GENES IN THE REGION
        future_genes = future({ findGenes(region_of_interest, reference_genome, genes_hg19 = gene.db, genes_hg38 = genes.hg38) })

      ## FIND SVs IN THE REGION
        future_svs = future({ findSVs(region_of_interest, reference_genome, span_value = 1000, sv_source, MAIN_PATH) })

      ## FIND GWAS CATALOG ENTRIES
        future_gwascat = future({ findGWAScat(region_of_interest, reference_genome, MAIN_PATH, gwascat_type) })
      ## MANAGE LD IF THIS WAS REQUESTED
        if (plotError == FALSE) { 
          results = value(future_res); 
          plotError = results[[1]]; data_to_plot = results[[2]];
          #rsid_region = value(future_rsid)
          if (plotError == FALSE && ld_type != 'No LD'){ future_ld = future({ findLD(ld_type, ld_pop, data_to_plot, rsid_region, reference_genome, region_of_interest, MAIN_PATH, snp_interest) }) }
        }

      ## GATHER ALL DATA BEFORE THE PLOT
        if (plotError == FALSE){
          recomb_data = value(future_rec)
          print('ok1')
          genes_in_region = value(future_genes)
          print('ok2')
          svs_in_region = value(future_svs)
          print('ok3')
          if (ld_type != 'No LD'){ ld = value(future_ld) } else { ld = NULL }
          eqtls_info = value(future_eqtl)
          print('ok4')
          sqtls_info = value(future_sqtl)
          print('ok5')
          gwascat_info = value(future_gwascat)
          print('ok6')
        }
      
      ## UPDATE REACTIVE OBJECT CONTAINING DATA TO PLOT
        if (plotError == FALSE){
          data_for_table(data_to_plot)
          plotError_for_table(plotError)
          #rsid_region_for_table(rsid_region)
          reference_genome(input$reference_gen)
          plotted_region(region_of_interest)
          genes_region(genes_in_region)
          svs_for_table(svs_in_region)
          eqtls_for_table(eqtls_info)
          sqtls_for_table(sqtls_info)
          gwascat_res(gwascat_info)
        }

      ## PLOT
        if (plotError == TRUE){
          main_fig = PlotError('Ops, no data in this region!')
        } else {
          message(paste0('** Chrom --> ', region_of_interest$chrom, '\n** Pos --> ', region_of_interest$start, '-', region_of_interest$end, '\n** Data --> ', paste(gwas_to_plot, collapse = ","), '\n\n'))
          main_fig = Plot(reference_genome, region_of_interest = region_of_interest, snps_data = data_to_plot, snp_interest = snp_interest, recomb_data = recomb_data, significance, 
            pos_start = region_of_interest$start, pos_end = region_of_interest$end, plot_type = input$plot_type, recomb = input$recomb_yn, genes_in_region, svs_in_region, showExons = input$exons, ld)
        }
        #mainFigure(main_fig); figureFormat(input$plot_format)
        main_fig
    })

  ## PLOT FUNCTION FOR PLOT 2 (GTEX)
    output$plot2 <- renderPlot({
      # get genes in the region and relative GTEx information
      genes_interest = genes_region()
      color_palette = input$heat_cols; tissues_interest = input$tissues
      if (nrow(genes_interest) >1 && ('All tissues' %in% tissues_interest || length(tissues_interest) >=2)){
        # define tissues and genes to show
        tmp_gtex = gtex.db[which(gtex.db$Description %in% genes_interest$Gene), ]
        #fig = plotGTEx(tmp_gtex, color_palette, tissues_interest)
        plotGTEx(tmp_gtex, color_palette, tissues_interest)
      } else {
        plot(0, 0, pch = 16, col = 'white', xlim = c(0, 1), ylim = c(0, 1), xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = '', bty = 'none')
        text(x = 0.5, y = 0.5, labels = 'Ops, too few genes for clustering\nTip: enlarge window ;)', cex = 3)
      }
    })

  ## PLOT TABLE OF SNPS
    output$table <- renderTable({
      plotError_for_table = plotError_for_table()
      if (plotError_for_table == TRUE){ 
        df = data.frame("RsID" = "No", "Chr" = "SNPs", "Pos" = "in", "P" = "the", "MAF" = "region", "Ref" = "of", "Alt" = "interest", "Study" = "!") 
      } else { 
        tmp_data = data_for_table()
        df = rbindlist(data_for_table())
        df$p = as.numeric(df$p); df = df[order(df$p),]; top = head(df, 80)
      }
      if (nrow(df) >0){
        df = df[, c("rsid", 'chrom', "pos", 'p', "maf", "ref", "alt", "name")]        
        colnames(df) = c('RsID', 'Chr', 'Pos', 'P', 'MAF', 'Ref', 'Alt', 'Study')
        # reformat
        df$MAF = round(as.numeric(df$MAF), 2); df$P = formatC(as.numeric(df$P), format = "e", digits = 1)
      } else {
        df = data.frame("RsID" = "No", "Chr" = "SNPs", "Pos" = "in", "P" = "the", "MAF" = "region", "Ref" = "of", "Alt" = "interest", "Study" = "!") 
      }
      return(df)
    }, bordered = TRUE, striped = TRUE)
  
  ## PLOT TABLE OF GWAS-CATALOG
    output$gwascat_table <- renderTable({
      plotError_for_table = plotError_for_table()
      if (plotError_for_table == TRUE){
        tmp_gwascat = data.frame(Chr = as.character(), Pos = as.character(), RsID = as.character(), P = as.character(), Gene = as.character(), Trait = as.character(), Pubmed = as.character())
      } else {
        tmp_gwascat = gwascat_res()
        genes_in_region = genes_region()
        df = rbindlist(data_for_table()); df$p = as.numeric(df$p); df = df[order(as.numeric(df$p)),]; df$p = -log10(df$p)
      }
      return(tmp_gwascat)
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
  
  ## PLOT TABLE OF EQTLS
    output$eqtls_table <- renderTable({
      # gather eqtl data and snp data showed and plotError
      eqtls_info = eqtls_for_table()
      snp_info = data_for_table()[[1]]
      plotError_for_table = plotError_for_table()
      if (plotError_for_table == TRUE | nrow(eqtls_info) == 0){
        tmp = data.frame(Chr = 'No', Pos = 'eQTLs', A1 = 'in', A2 = 'the', Tissue = 'region', Gene = 'of', Effect = 'interest', P = '!')
      } else {
        tmp_snp = snp_info[, c('rsid', 'pos', 'maf')]; colnames(tmp_snp) = c('RsID', 'Pos', 'MAF')
        tmp_snp$Pos = as.character(tmp_snp$Pos); eqtls_info$Pos = as.character(eqtls_info$Pos)
        tmp = merge(tmp_snp, eqtls_info, by = 'Pos')
        # adjust formats
        tmp$MAF = round(as.numeric(tmp$MAF), 2); tmp$Effect = round(as.numeric(tmp$Effect), 2)
        tmp$P = formatC(as.numeric(tmp$P), format = "e", digits = 1)
      }
      return(tmp)
    }, bordered = TRUE, striped = TRUE)

  ## PLOT TABLE OF SQTLS
    output$sqtls_table <- renderTable({
      # gather eqtl data and snp data showed and plotError
      sqtls_info = sqtls_for_table()
      snp_info = data_for_table()[[1]]
      plotError_for_table = plotError_for_table()
      if (plotError_for_table == TRUE | nrow(sqtls_info) == 0){
        tmp = data.frame(Chr = 'No', Pos = 'sQTLs', A1 = 'in', A2 = 'the', Tissue = 'region', Gene = 'of', Effect = 'interest', P = '!')
      } else {
        tmp_snp = snp_info[, c('rsid', 'pos', 'maf')]; colnames(tmp_snp) = c('RsID', 'Pos', 'MAF')
        tmp_snp$Pos = as.character(tmp_snp$Pos); sqtls_info$Pos = as.character(sqtls_info$Pos)
        tmp = merge(tmp_snp, sqtls_info, by = 'Pos')
        # adjust formats
        tmp$MAF = round(as.numeric(tmp$MAF), 2); tmp$Effect = round(as.numeric(tmp$Effect), 2)
        tmp$P = formatC(as.numeric(tmp$P), format = "e", digits = 1)
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
    #output$download_GWASCat <- downloadHandler(
    #  filename = function(){ gwascat_type = gwascat_mode(); reference_gen = reference_genome(); 
    #    if (reference_gen == 'GRCh37 (hg19)'){ genome = 'hg19' } else { genome = 'hg38' }
    #    if (gwascat_type == 'SNPs'){ paste0("SNPs_information_from_GWASCatalog_", genome, ".txt") } else { paste0("Genes_information_from_GWASCatalog_", genome, ".txt") } },
    #  content = function(file) {
    #    plotError_for_table = plotError_for_table()
    #    if (plotError_for_table == TRUE){
    #      tmp_gwascat = data.frame(Chr = as.character(), Pos = as.character(), P = as.character())
    #    } else {
    #      reference_gen = reference_genome()
    #      region_of_interest = plotted_region()
    #      gwascat_type = gwascat_mode()
    #      genes_in_region = genes_region()
    #      df = rbindlist(data_for_table()); df$p = as.numeric(df$p); df = df[order(as.numeric(df$p)),]; df$p = -log10(df$p)
    #      if (gwascat_type == 'SNPs'){
    #        tmp_gwascat = gwascat[[as.numeric(region_of_interest$chrom)]]
    #        if (reference_gen == 'GRCh37 (hg19)') { tmp_gwascat = tmp_gwascat[which(tmp_gwascat$'Pos (hg19)' %in% df$pos),]; tmp_gwascat = tmp_gwascat[order(as.numeric(tmp_gwascat$'P-value')),] } else { tmp_gwascat = tmp_gwascat[which(tmp_gwascat$'Pos (hg38)' %in% df$pos),]; tmp_gwascat = tmp_gwascat[order(as.numeric(tmp_gwascat$'P-value')),] }
    #        tmp_gwascat = tmp_gwascat[!duplicated(tmp_gwascat[c('SNP', 'Pubmed ID')]),]
    #      } else {
    #        tmp_gwascat = gwascat_genes[which(gwascat_genes$gene %in% genes_in_region$Gene),]
    #        colnames(tmp_gwascat) = c('Gene', 'Trait', 'Study', 'Journal')
    #      }
    #    }
    #    if (nrow(tmp_gwascat) >0){
    #      if (gwascat_type == 'SNPs'){ tmp_gwascat = tmp_gwascat[, c('SNP', 'Gene', 'Trait', 'Pubmed ID', 'P-value')] }
    #    } else {
    #      if (gwascat_type == 'SNPs'){ tmp_gwascat = data.frame('SNP' = 'No', 'Gene' = 'entries', 'Trait' = 'in', 'Pubmed ID' = 'GWAS', 'P-value' = 'Catalog') } else { tmp_gwascat = data.frame('Gene' = 'No entries', 'Trait' = 'in the', 'Study' = 'GWAS', 'Journal' = 'Catalog') } 
    #    }
    #    write.table(tmp_gwascat, file, row.names = FALSE, quote=F, sep = "\t")
    #  })

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

    output$ld_hub_link <- renderUI({
      link <- paste0("http://ldsc.broadinstitute.org")
      url <- a("Click here!", href = link)
      tagList("LD Hub link:", url)
    })
}

