## LIBRARIES
  suppressPackageStartupMessages({
    library(shiny)
    library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
    library(plotly)
    library(data.table)
    library(stringr)
  })

## ANNOTATIONS AND FUNCTIONS
  source("/Users/nicco/Documents/GitHub/snpXplorer/bin/snpxplorer_v4/bin/functions.R")
  load("/Users/nicco/Documents/GitHub/snpXplorer/bin/snpxplorer_v4/data/databases/annotationFiles.RData")

## MAIN SERVER FUNCTION
server <- function(input, output) {  
  # REACTIVE OBJECTS
  data_for_table = reactiveVal()
  plotError_for_table = reactiveVal()

  # PLOTLY FUNCTION
  output$plot <- renderPlotly({
    # FIND POSITIONS TO PLOT
    if (input$browsing_options == 'Locus, Gene, RsID'){ results = directInput(target = input$target, input$x_axis, gene.db, genes.hg38, input$reference_gen, SNPlocs.Hsapiens.dbSNP144.GRCh37); region_of_interest = results[[1]]; snp_interest = results[[2]] } else { region_of_interest = data.frame(chrom = input$manual_chrom, start = input$manual_pos - input$x_axis, end = input$manual_pos + input$x_axis); snp_interest = NA }
    if (is.na(region_of_interest$start)){ plotError = TRUE } else { plotError = FALSE }
    message('** positions of interest FOUND!')

    # FIND DATA TO PLOT IN CASE THE INPUT BROWSING OPTION WAS VALID
    if ((plotError == FALSE) & (is.null(input$upload_file) == FALSE) & (length(input$owned_gwas) >0)) { uploaded_file = input$upload_file$datapath } else { uploaded_file = NULL }
    gwas_to_plot = c(input$neurol_gwas, input$cardio_gwas, input$immune_gwas, input$cancer_gwas, input$physio_gwas, uploaded_file)
    
    # EXTRACT DATA FROM DATA TO PLOT IN CASE THE INPUT BROWSING OPTION WAS VALID
    all_colors = c(input$col1, input$col2, input$col3, input$col4, input$col5)
    if (plotError == FALSE) { results = extractDataForPlot(gwas_to_plot, region_of_interest, span = 500000, colors = all_colors, res_example); plotError = results[[1]]; data_to_plot = results[[2]] } else { data_to_plot = NA}
    message('** data to plot FOUND!')

    # FIND RECOMBINATION RATES
    if (plotError == FALSE) { recomb_data = extractRecombination(region_of_interest, reference_genome = input$reference_gen, genetic_map = genetic.map, genetic_map_hg38 = genetic.map.hg38, span = 500000) }
    message('** recombination rates FOUND!')

    # FIND GENES IN THE REGION
    genes_in_region = findGenes(region_of_interest, reference_genome = input$reference_gen, genes_hg19 = gene.db, genes_hg38 = genes.hg38)
    message('** genes in the region FOUND!\n')

    # FIND SVs IN THE REGION
    svs_in_region = findSVs(region_of_interest, reference_genome = input$reference_gen, all_str, all_str_hg38, span_value = 500000, sv_source = input$strVar_inp)

    # UPDATE REACTIVE OBJECT CONTAINING DATA TO PLOT
    data_for_table(data_to_plot)
    plotError_for_table(plotError)

    # PLOT
    if (plotError == TRUE){
      PlotError('Ops, no data in this region!')
    } else {
      main_fig = Plot(snps_data = data_to_plot, snp_interest = snp_interest, recomb_data = recomb_data, significance = input$y_axis, 
        pos_start = region_of_interest$start, pos_end = region_of_interest$end, plot_type = input$plot_type, 
        recomb = input$recomb_yn, genes_in_region, svs_in_region, showExons = input$exons)
    }
    #main_fig = Plot(snps_data = data_to_plot, snp_interest = snp_interest, recomb_data = recomb_data, significance = 8, pos_start = region_of_interest$start, pos_end = region_of_interest$end, plot_type = 'Scatter', recomb = 'Yes', genes_in_region = genes_in_region)    
  })

  # PLOT TABLE OF SNPS
  output$table <- renderTable({
    plotError_for_table = plotError_for_table()
    if (plotError_for_table == TRUE){ df = data.frame(Chr = as.character(), Position = as.character(), "-log10(P)" = as.character(), Color = as.character(), Study = as.character()) } else { df = rbindlist(data_for_table()); df$p = as.numeric(df$p); df = df[order(as.numeric(df$p)),]; df$p = -log10(df$p) }
    if (nrow(df) >0){
        colnames(df) <- c("Chr", "Position", "-log10(P)", "Color", "Study")
        df = df[, c("Chr", "Position", "-log10(P)", "Study")]
        top <- head(df, 10)
      } else {
        top = data.frame(Chr = "No SNPs", Position = "in the", "-log10(P)" = "region of", Study = "interest")
        colnames(top) <- c("Chr", "Position", "-log10(P)", "Study")
      }
      return(top)
  })
  
}

