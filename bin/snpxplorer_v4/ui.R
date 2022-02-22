library(shiny)
library(plotly)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("snpXplorer v4"),
  
  # Sidebar with a slider input for position and width of the region
  sidebarLayout(
    sidebarPanel(
      sliderInput("position", "Number of bins:", min = 1, max = 250000000, value = 15000000),
      sliderInput("window", "Width of the window:", min = 100, max = 500000, value = 100000),
      sliderInput("significance", "-log10(P-value)", min = 0, max = 250, value = 8),
      radioButtons(inputId = "plot_type", label = "Select plot type", choices = c("Points", "Lines")),
      conditionalPanel(condition = "input.plot_type == 'Lines'", 
                       sliderInput(inputId = 'sliding.window', label = "Number of bins", value = 20, min=20, max=40),
                       sliderInput(inputId = "smooth", label = "Smoothing parameter", value = 0.10, min=0.10, max=0.5))),
    
    # Show a plot of the generated distribution
    mainPanel( plotlyOutput("plot", height = "500px") )
  )
)
