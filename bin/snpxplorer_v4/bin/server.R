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
  source("/Users/nicco/Documents/GitHub/snpXplorer/bin/snpxplorer_v4/bin/functions.R")
  load("/Users/nicco/Documents/GitHub/snpXplorer/bin/snpxplorer_v4/data/databases/annotationFiles.RData")

## MAIN SERVER FUNCTION
server <- function(input, output) {  
  ## REACTIVE OBJECTS
    data_for_table = reactiveVal()
    plotError_for_table = reactiveVal()
    rsid_region_for_table = reactiveVal()
    reference_genome = reactiveVal()
    plotted_region = reactiveVal()

  ## PLOTLY FUNCTION FOR PLOT 1 (SNPS, GENES AND SVS)
    output$plot <- renderPlotly({
      ## FIND POSITIONS TO PLOT
        if (input$browsing_options == 'Locus, Gene, RsID'){ results = directInput(target = input$target, input$x_axis, gene.db, genes.hg38, input$reference_gen, SNPlocs.Hsapiens.dbSNP144.GRCh37); region_of_interest = results[[1]]; snp_interest = results[[2]] } else { region_of_interest = data.frame(chrom = input$manual_chrom, start = input$manual_pos - input$x_axis, end = input$manual_pos + input$x_axis); snp_interest = NA }
        if (is.na(region_of_interest$start)){ plotError = TRUE } else { plotError = FALSE }
        message('** positions of interest FOUND!')
      
      ## FIND RSID FOR ANNOTATION
        if ((plotError == FALSE)){ future_rsid = future({ results = findRsID(region_of_interest = region_of_interest, reference_genome = input$reference_gen, span = 500000) }) }
        message('** rsid in the region FOUND!')
      
      ## FIND DATA TO PLOT IN CASE THE INPUT BROWSING OPTION WAS VALID
        if ((plotError == FALSE) & (is.null(input$upload_file) == FALSE) & (length(input$owned_gwas) >0)) { uploaded_file = input$upload_file$datapath } else { uploaded_file = NULL }
        gwas_to_plot = c(input$neurol_gwas, input$cardio_gwas, input$immune_gwas, input$cancer_gwas, input$physio_gwas, uploaded_file)
      ## EXTRACT DATA FROM DATA TO PLOT IN CASE THE INPUT BROWSING OPTION WAS VALID
        all_colors = c(input$col1, input$col2, input$col3, input$col4, input$col5)
        if (plotError == FALSE) { future_res = future({ results = extractDataForPlot(gwas_to_plot, region_of_interest, span = 500000, colors = all_colors, res_example, input$reference_gen) }) } else { data_to_plot = NA }
        message('** data to plot FOUND!')

      ## FIND RECOMBINATION RATES
        if (plotError == FALSE) { future_rec = future({ results = extractRecombination(region_of_interest, reference_genome = input$reference_gen, genetic_map = genetic.map, genetic_map_hg38 = genetic.map.hg38, span = 500000) }) } 
        message('** recombination rates FOUND!')

      ## FIND GENES IN THE REGION
        future_genes = future({ results = findGenes(region_of_interest, reference_genome = input$reference_gen, genes_hg19 = gene.db, genes_hg38 = genes.hg38) })
        message('** genes in the region FOUND!\n')

      ## FIND SVs IN THE REGION
        future_svs = future({ results = findSVs(region_of_interest, reference_genome = input$reference_gen, all_str, all_str_hg38, span_value = 500000, sv_source = input$strVar_inp) })

      ## GATHER ALL DATA BEFORE THE PLOT
        rsid_region = future_rsid$result$value
        results = future_res$result$value; plotError = results[[1]]; data_to_plot = results[[2]]
        recomb_data = future_rec$result$value
        genes_in_region = future_genes$result$value
        svs_in_region = future_svs$result$value

      ## UPDATE REACTIVE OBJECT CONTAINING DATA TO PLOT
        data_for_table(data_to_plot)
        plotError_for_table(plotError)
        rsid_region_for_table(rsid_region)
        reference_genome(input$reference_gen)
        plotted_region(region_of_interest)

      ## PLOT
        if (plotError == TRUE){
          PlotError('Ops, no data in this region!')
        } else {
          main_fig = Plot(reference_genome = input$reference_gen, region_of_interest = region_of_interest, rsid_region = rsid_region, snps_data = data_to_plot, snp_interest = snp_interest, recomb_data = recomb_data, significance = input$y_axis, 
            pos_start = region_of_interest$start, pos_end = region_of_interest$end, plot_type = input$plot_type, 
            recomb = input$recomb_yn, genes_in_region, svs_in_region, showExons = input$exons)
        }
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
      rsid_region = rsid_region_for_table()
      reference_gen = reference_genome()
      if (plotError_for_table == TRUE){ 
        df = data.frame("Rs ID" = as.character(), "Pos" = as.character(), "Ref" = as.character(), "Alt" = as.character(), "MAF" = as.character(), "Log(p)" = as.character(), "Study" = as.character()) 
      } else { 
        df = rbindlist(data_for_table()); df$p = as.numeric(df$p); df = df[order(as.numeric(df$p)),]; df$p = -log10(df$p); top = head(df, 8); 
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
    })
  
  ## PLOT TABLE OF GWAS-CATALOG
    output$gwascat_table <- renderTable({
      plotError_for_table = plotError_for_table()
      reference_gen = reference_genome()
      region_of_interest = plotted_region()
      if (plotError_for_table == TRUE){
        tmp_gwascat = data.frame(Chr = as.character(), Pos = as.character(), P = as.character())
      } else {
        df = rbindlist(data_for_table()); df$p = as.numeric(df$p); df = df[order(as.numeric(df$p)),]; df$p = -log10(df$p)
        tmp_gwascat = gwascat[[region_of_interest$chrom]]
        if (reference_gen == 'GRCh37 (hg19)') { 
          tmp_gwascat = tmp_gwascat[which(tmp_gwascat$'Pos (hg19)' %in% df$pos),]
          tmp_gwascat = tmp_gwascat[order(as.numeric(tmp_gwascat$'P-value')),]
        } else {
          tmp_gwascat = tmp_gwascat[which(tmp_gwascat$'Pos (hg38)' %in% df$pos),]
          tmp_gwascat = tmp_gwascat[order(as.numeric(tmp_gwascat$'P-value')),]
        }
      }
      if (nrow(tmp_gwascat) >0){
        top = head(tmp_gwascat, 10); top = top[, c('SNP', 'Gene', 'Trait', 'Pubmed ID', 'P-value')]
        } else {
        top = data.frame('SNP' = 'No', 'Gene' = 'entries', 'Trait' = 'in', 'Pubmed ID' = 'GWAS', 'P-value' = 'Catalog')
      }
      return(top)
    })
}

