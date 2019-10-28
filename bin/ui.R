###SNP BROWSER
###

#LIBRARIES
library(shiny)
library(data.table)

#need to do this otherwise error -- but weird, don't get why!
dat <- fread("../data/example/chr21_IGAP_2k19.txt.gz", h=T, check.names = F)
colnames(dat) <- c("chr", "pos", "p")      
###

#increase maximum size of uploaded file -- now is 30MB
options(shiny.maxRequestSize=30*1024^2)

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
            checkboxGroupInput("gwas", "GWAS to show:", c("Your input" = "toLoad"), inline = T),
        
            #select input file and/or additional files
            fileInput(inputId = "inp_f", label = "Add file", multiple = F, placeholder = "Loading file...")
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
            radioButtons(inputId = "sel", label = "Select browser type", choices = c("Position", "Gene", "Rs ID", "Manual scroll")),
            
            br(),
        
            #conditional panel for position
            conditionalPanel(condition = "input.sel == 'Position'", 
                         textInput(inputId = "pos", label = "Position (GRCh37) [Chr:Pos]", value = "Type position...")),
        
            #conditional panel for gene
            conditionalPanel(condition = "input.sel == 'Gene'",
                         textInput(inputId = "gene", label = "Gene", value = "Type gene symbol...")),
    
            #conditional panel for rsID
            conditionalPanel(condition = "input.sel == 'Rs ID'", 
                         textInput(inputId = "snpID", label = "SNP Rs ID", value = "Type variant identifier...")),
    
            #select chromosome in case of manual scroll -- this has to be fixed and implemented
            conditionalPanel(condition = "input.sel == 'Manual scroll'",
                         sliderInput(inputId = 'manual.chrom', label = "Chromosome", 
                                     value = 11, min=1, max=22, width = "20%"),
                         sliderInput(inputId = 'manual.pos', label = "Position", 
                                     dat[nrow(dat)/2, "pos"], min = min(dat$pos),
                                     max = max(dat$pos), width = "40%"))
            )
          )
        ),
    
    hr(),
        
      fluidRow(  
        column(3,
               wellPanel(
                 
                 h3("Visualization options"),
                 #sliding window -- x.axis
                 sliderInput(inputId = "x", label = "Window (bp)",
                             value = 300000, min = 1000, max = 1000000),
                 
                 #sliding widow -- y.axis
                 sliderInput(inputId = "y", label = "-log10 (P-value)", value = 9, min = 0, max = 200)
               )
        ),
        
        column(6,
          #output
          plotOutput(outputId = "hist", height = "600px", click = "plot1_click", brush = brushOpts(id = "plot1_brush"))
        ),
    
      column(3,
        wellPanel(
            #select type of plot -- lines or points
            radioButtons(inputId = "ploType", label = "Select plot type", choices = c("Points", "Lines")),
               
            #conditional panel for line type -- adjusting smoothing parameters for density line
            conditionalPanel(condition = "input.ploType == 'Lines'",
                            sliderInput(inputId = 'sliding.window', label = "N of sliding windows", 
                                            value = 30, min=10, max=50),
                            sliderInput(inputId = "smooth", label = "Smoothing parameter",
                                            value = 0.10, min=0, max=1)),
            
            #conditional panel for click and brush of the points -- only when points is plot-type
            conditionalPanel(condition = "input.ploType == 'Points'",
                             #text Panel for points near click
                             h4("Points near click"),
                             verbatimTextOutput("click_info"),
                             br(),
                             
                             h4("Brushed points"),
                             verbatimTextOutput("brush_info"),
                             br()
          )
        )
      )
    ),
    hr()
  )
)
