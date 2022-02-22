library(shiny)
library(plotly)

# Define server logic required to draw a histogram
server <- function(input, output) {
  # read association data
  d = data.table::fread("/Users/nicco/Downloads/trial_dataset_snpXplorer_exploration.txt", h=T, stringsAsFactors = F)
  # recombination rates
  recomb = data.table::fread("/Users/nicco/Desktop/genetic_map_chr16_combined_b37.txt", h=T, stringsAsFactors = F)
  colnames(recomb) = c("pos", "combined", "cm")
  
  # main function for plot
  output$plot <- renderPlotly({
    # define positions
    pos_start = input$position - input$window
    pos_end = input$position + input$window
    span = 1000000
    
    # prepare for plot
    snps_data = d[which((d$pos >= pos_start - span) & (d$pos <= pos_end + span)),]
    recomb_data = recomb[which((recomb$pos >= pos_start - span) & (recomb$pos <= pos_end + span)),]
    
    # plot
    if (input$plot_type == 'Points'){ 
      main_fig = plot_points(snps_data, recomb_data, input$significance, pos_start, pos_end)
    } else {
      main_fig = plot_density(snps_data, recomb_data, input$significance, pos_start, pos_end, input$sliding.window, input$smooth)
    }
    
  })
}

