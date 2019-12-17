###SNP BROWSER
###

#LIBRARIES
library(shiny)
library(data.table)
library(rvest)
library(stringr)

library("colourpicker")

#increase maximum size of uploaded file -- now is 60MB
options(shiny.maxRequestSize=60*1024^2)

#USER INTERFACE
shinyServer(
  fluidPage(
    #title of the app
    titlePanel(title = "SNP browser v1.0", windowTitle = "SNPbrowser session"),
    
    hr(),

    fluidRow(
        
      column(5,
          wellPanel(
            h3("Input selection"),
            #select which gwas as input -- default is ad
            checkboxGroupInput("gwas", "GWAS to show:", c("IGAP" = "IGAP", "CARDIO" = "CARDIO", "Your input" = "toLoad"), inline = T),

            #select input file and/or additional files
            fileInput(inputId = "inp_f", label = "Add file", multiple = F, placeholder = "Loading file..."),
            
            #help message
            h6("If you upload your file, make sure data has at least chromosome, position and pvalue columns. 
                      Works with most headers. See doc for more info.")
            )
        ),

        column(2,
          wellPanel(
            
            h3("Genome options"),
            #select genome version -- to be implemented
            radioButtons(inputId = "genV", label = "Select genome version", choices = c("GRCh37 (hg19)", "GRCh38 (hg38)"))
          )
        ),
      
      column(5,
        wellPanel(
          
            h3("Browsing options"),
            #select input type
            radioButtons(inputId = "sel", label = "Select browser type", choices = c("Position", "Gene", "RsID", "Manual scroll")),
            
            br(),
        
            #conditional panel for position
            conditionalPanel(condition = "input.sel == 'Position'", 
                         textInput(inputId = "pos", label = "Position (GRCh37) [Chr:Pos]", value = "Type position...")),
        
            #conditional panel for gene
            conditionalPanel(condition = "input.sel == 'Gene'",
                         textInput(inputId = "gene", label = "Gene", value = "Type gene symbol...")),
    
            #conditional panel for rsID
            conditionalPanel(condition = "input.sel == 'RsID'", 
                         textInput(inputId = "snpID", label = "SNP Rs ID", value = "Type variant identifier...")),
    
            #select chromosome in case of manual scroll -- this has to be fixed and implemented
            conditionalPanel(condition = "input.sel == 'Manual scroll'",
                         sliderInput(inputId = 'manual.chrom', label = "Chromosome", 
                                     value = 11, min=1, max=22, width = "60%"),
                         sliderInput(inputId = 'manual.pos', label = "Position", 
                                     15000000, min = 1,
                                     max = 200000000, width = "100%"))
            )
          )
        ),
    
    hr(),
        
      fluidRow(  
        column(3,
               wellPanel(
                 
                 h3("Visualization options"),
                 
                 fluidRow(
                   column(6,
                          #choose color
                          colourInput("col", "Select colour", "black"),
                          
                          #conditional panel for multiple points
                          conditionalPanel(condition = "input.gwas != undefined && input.gwas.length > 1",
                                           colourInput("col2", "Select color 2", "coral"))
                          )
                 ),
                 
                 #sliding window -- x.axis
                 sliderInput(inputId = "x", label = "Window (bp)",
                             value = 300000, min = 1000, max = 2000000),
                 
                 #sliding widow -- y.axis
                 sliderInput(inputId = "y", label = "-log10 (P-value)", value = 9, min = 0, max = 100),
                 
                 #select type of plot -- lines or points
                 radioButtons(inputId = "ploType", label = "Select plot type", choices = c("Points", "Lines")),
                        
                 #conditional panel for line type -- adjusting smoothing parameters for density line
                 conditionalPanel(condition = "input.ploType == 'Lines'",
                                         sliderInput(inputId = 'sliding.window', label = "N of sliding windows", 
                                                     value = 30, min=10, max=50),
                                         sliderInput(inputId = "smooth", label = "Smoothing parameter",
                                                     value = 0.10, min=0, max=1)),
                        
                 #select type of LD -- no LD, most significant in window, Input variant or Other (other makes appear a text panel)
                 radioButtons(inputId = "Linkage", label = "LD options", choices = c("No LD", "Input variant", "Most significant in window", "Other variant")),
                 
                 #Download button
                 h5("Download your plot"),
                 downloadButton('downloadPlot', 'Download Plot'),
                 
                 #here is how to make the LD "Other" panel appear
                 conditionalPanel(condition = "input.Linkage == 'Other variant'",
                                          textInput(inputId = "ld_var", label = "Input variant [Chr:Pos]", value = "Type position...")),
               
                 
                 #conditional panel for click and brush of the points -- only when points is plot-type
                 conditionalPanel(condition = "input.ploType == 'Points'",
                                         #text Panel for points near click
                                         h5("Points near click"),
                                         verbatimTextOutput("click_info"),

                                         h5("Brushed points"),
                                         verbatimTextOutput("brush_info")
                    )
                  )
               ),
               
        column(9,
          #output
          plotOutput(outputId = "hist", height = "700px", click = "plot1_click", brush = brushOpts(id = "plot1_brush"))
        )
    ),
    hr()
  )
)
