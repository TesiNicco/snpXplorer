## LIBRARIES
  suppressPackageStartupMessages({
    library(shiny)
    library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
    library(plotly)
    library(data.table)
    library(stringr)
    library(future)
  })

## ANNOTATIONS AND FUNCTIONS
  MAIN_PATH = '/Users/nicco/Documents/GitHub/snpXplorer/bin/snpxplorer_v4/'
  source(paste0(MAIN_PATH, "bin/functions.R"))
  load(paste0(MAIN_PATH, "data/databases/20220310_annotationFiles.RData"))
  #MAIN_PATH = '/root/snpXplorer/'
  #source(paste0(MAIN_PATH, "snpXplorer_v4/functions.R"))
  #load(paste0(MAIN_PATH, "data/databases/20220310_annotationFiles.RData"))

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

  ## PLOTLY FUNCTION FOR PLOT 1 (SNPS, GENES AND SVS)
    output$plot <- renderPlotly({
      ## DEFINE FUTURE CLASS
        plan(multicore)

      ## DEFINE THE REACTIVE VARIABLES
        browsing_options = input$browsing_options; target = input$target; x_axis = input$x_axis; significance = input$y_axis
        reference_genome = input$reference_gen; sv_source = input$strVar_inp; all_colors = c(input$col1, input$col2, input$col3, input$col4, input$col5)
        ld_type = input$linkage_type; ld_pop = c(input$pop_eur, input$pop_afr, input$pop_amr, input$pop_eas, input$pop_sas)

      ## FIND POSITIONS TO PLOT
        if (browsing_options == 'Locus, Gene, RsID'){ results = directInput(target, x_axis, gene.db, genes.hg38, reference_genome, SNPlocs.Hsapiens.dbSNP144.GRCh37, MAIN_PATH); region_of_interest = results[[1]]; snp_interest = results[[2]] } else { region_of_interest = data.frame(chrom = input$manual_chrom, start = input$manual_pos - input$x_axis, end = input$manual_pos + input$x_axis); snp_interest = NA }
        if (is.na(region_of_interest$start)){ plotError = TRUE } else { plotError = FALSE }
        message('** positions of interest FOUND!')
      
      ## FIND RSID FOR ANNOTATION
        if ((plotError == FALSE)){ future_rsid = future({ findRsID(region_of_interest = region_of_interest, reference_genome = reference_genome, span = 100000, MAIN_PATH) }) }
        message('** rsid in the region FOUND!')
            
      ## FIND DATA TO PLOT IN CASE THE INPUT BROWSING OPTION WAS VALID
        if ((plotError == FALSE) & (is.null(input$upload_file) == FALSE) & (length(input$owned_gwas) >0)) { uploaded_file = input$upload_file$datapath } else { uploaded_file = NULL }
        gwas_to_plot = c(input$neurol_gwas, input$cardio_gwas, input$immune_gwas, input$cancer_gwas, input$physio_gwas, uploaded_file)

      ## EXTRACT DATA FROM DATA TO PLOT IN CASE THE INPUT BROWSING OPTION WAS VALID
        if (plotError == FALSE) { future_res = future({ extractDataForPlot(gwas_to_plot, region_of_interest, span = 100000, colors = all_colors, res_example, reference_genome, MAIN_PATH) }) } else { data_to_plot = NA }
        message('** data to plot FOUND!')

      ## FIND RECOMBINATION RATES
        if (plotError == FALSE) { future_rec = future({ extractRecombination(region_of_interest, reference_genome, genetic_map = genetic.map, genetic_map_hg38 = genetic.map.hg38, span = 100000) }) } 
        message('** recombination rates FOUND!')

      ## FIND GENES IN THE REGION
        future_genes = future({ findGenes(region_of_interest, reference_genome, genes_hg19 = gene.db, genes_hg38 = genes.hg38) })
        message('** genes in the region FOUND!')

      ## FIND SVs IN THE REGION
        future_svs = future({ findSVs(region_of_interest, reference_genome, all_str, all_str_hg38, span_value = 100000, sv_source) })
        message('** svs in the region FOUND!\n')
      
      ## MANAGE LD IF THIS WAS REQUESTED
        if (plotError == FALSE) { 
          results = value(future_res); plotError = results[[1]]; data_to_plot = results[[2]];
          rsid_region = value(future_rsid)
          if (plotError == FALSE && ld_type != 'No LD'){ future_ld = future({ findLD(ld_type, ld_pop, data_to_plot, rsid_region, reference_genome, region_of_interest, MAIN_PATH) }) }
        }

      ## GATHER ALL DATA BEFORE THE PLOT
        if (plotError == FALSE){
          recomb_data = value(future_rec)
          genes_in_region = value(future_genes)
          svs_in_region = value(future_svs)
          if (ld_type != 'No LD'){ ld = value(future_ld) } else { ld = NULL }
        }
      
      ## UPDATE REACTIVE OBJECT CONTAINING DATA TO PLOT
        if (plotError == FALSE){
          data_for_table(data_to_plot)
          plotError_for_table(plotError)
          rsid_region_for_table(rsid_region)
          reference_genome(input$reference_gen)
          plotted_region(region_of_interest)
          gwascat_mode(input$gwascat_type)
          genes_region(genes_in_region)
          svs_for_table(svs_in_region)
        }

      ## PLOT
        if (plotError == TRUE){
          main_fig = PlotError('Ops, no data in this region!')
        } else {
          main_fig = Plot(reference_genome, region_of_interest = region_of_interest, rsid_region = rsid_region, snps_data = data_to_plot, snp_interest = snp_interest, recomb_data = recomb_data, significance, 
            pos_start = region_of_interest$start, pos_end = region_of_interest$end, plot_type = input$plot_type, recomb = input$recomb_yn, genes_in_region, svs_in_region, showExons = input$exons, ld)
        }
        mainFigure(main_fig); figureFormat(input$plot_format)
        main_fig
    })

  ## PLOTLY FUNCTION FOR PLOT 2 (SQTLS AND EQTLS)
    output$plot2 <- renderPlotly({
      USPersonalExpenditure <- data.frame("Categorie"=rownames(USPersonalExpenditure), USPersonalExpenditure)
      data <- USPersonalExpenditure[,c('Categorie', 'X1960')]

      fig <- plot_ly(data, labels = ~Categorie, values = ~X1960, type = 'pie')
      fig <- fig %>% layout(title = 'United States Personal Expenditures by Categories in 1960',
          xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
          yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
    })

  ## PLOT TABLE OF SNPS
    output$table <- renderTable({
      plotError_for_table = plotError_for_table()
      if (plotError_for_table == TRUE){ 
        df = data.frame("Rs ID" = as.character(), "Pos" = as.character(), "Ref" = as.character(), "Alt" = as.character(), "MAF" = as.character(), "Log(p)" = as.character(), "Study" = as.character()) 
      } else { 
        rsid_region = rsid_region_for_table()
        reference_gen = reference_genome()
        df = rbindlist(data_for_table()); df$p = as.numeric(df$p); df = df[order(as.numeric(df$p)),]; df$p = -log10(df$p); top = head(df, 80); 
        if (reference_gen == 'GRCh37 (hg19)'){ top = merge(top, rsid_region, by.x = 'pos', by.y = 'POS', all.x = T) } else { top = merge(top, rsid_region, by.x = 'pos', by.y = 'POS_HG38', all.x = T) }
      }
      if (nrow(df) >0){
          top = top[, c('ID', 'pos', 'REF', 'ALT', 'ALT_FREQS', 'p', 'name')]
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
}

