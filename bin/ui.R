#################################
# libraries
suppressPackageStartupMessages({
  library(shiny)
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(liftOver)
  library(colourpicker)
  library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
  library(rvest)
  library(stringr)
  library(plotrix)
})
# increase maximum size of uploaded file -- now is 300MB
options(shiny.maxRequestSize=600*1024^2)

#################################

#################################
# main
navbarPage("snpXplorer", inverse = T,
           # Exploration TAB
           tabPanel("Exploration",
                    # DNA image
                    fluidRow(column(12, img(src='dna.png', align="center", width="100%"))),
                    
                    # Tutorial image
                    fluidRow(column(12, img(src='tutorial.png', align="center", width="100%"))),
                    
                    hr(),
                    
                    wellPanel(
                    fluidRow(
                      # first put the input selection
                      column(4,
                         # title for input selection
                         h3("Input selection"),
                         h6("Please check up to 5 boxes at the same time."),
                         # select which gwas as input -- try to do based on disease type
                         checkboxGroupInput("gwas_neuro", "Neurological traits", c("IGAP" = "IGAP", "proxy_AD" = "proxy_AD", "Autism" = "Autism", "Depression" = "Depression"), inline = T),
                         checkboxGroupInput("gwas_cardio", "Cardiovascular traits", c("CAD" = "CAD", "Ventricular volume" = "LVV", "SBP" = "SBP", "BMI" = "BMI", "Diabetes" = "Diabetes"), inline = T),
                         checkboxGroupInput("gwas_immune", "Immune-related traits", c("COVID" = "COVID", "Lupus" = "Lupus", "Inflammation" = "Inflammation", "Asthma" = "Asthma"), inline = T),
                         checkboxGroupInput("gwas_cancer", "Cancer-related traits", c("Breast_cancer" = "Breast_cancer", "Myeloproliferative" = "Myeloproliferative", "Prostate" = "Prostate", "Lung" = "Lung", "Leukemia" = "Leukemia"), inline = T),
                         checkboxGroupInput("gwas_physio", "Physiological traits", c("UKBaging" = "UKBaging", "Height" = "Height", "Education" = "Education", "Bone_density" = "Bone_density", "Vitamin_D" = "Vitamin_D"), inline = T)),

                      # next put the genome options and the personal file upload
                      column(4,
                         br(),
                         # select which gwas as input from the user only
                         checkboxGroupInput("gwas_from_file", "Your input", c("Your input" = "toLoad"), inline = T),
                         # select user input file
                         fileInput(inputId = "inp_f", label = "Add file", multiple = F, placeholder = "Loading file", width="300px"),
                         # help message
                         div(style = "font-size: 10px; padding: 0px 0px; margin-top:-2em", h6("If you upload your file, make sure data has at least chromosome, position and pvalue columns. 
                         See doc for more info.")),
                         div(style = "border-top: 2px solid #999999"),
                         h3("Reference Genome"),
                         # select genome options
                         radioButtons(inputId = "genV", label = "Select genome version", choices = c("GRCh37 (hg19)", "GRCh38 (hg38)")),
                         div(style = "border-top: 2px solid #999999"),
                      ),
                         
                      # next column is for the browsing option
                      column(4,
                         h3("Browsing options"),
                         # radio buttons for browsing option
                         radioButtons(inputId = "sel", label = "Select browser type", choices = c("Locus, Gene, RsID", "Manual scroll")),
                               
                         # conditional panel for position
                         conditionalPanel(condition = "input.sel == 'Locus, Gene, RsID'", 
                              textInput(inputId = "target", label = "Insert your target region", value = "Type locus, rsID or gene name", width="400px"),
                              h5("Locus: genomic location of interest. The format is chromosome:position (e.g 1:12345678)"),
                              h5("RsID: variant of interest. The format is rsXXXX (e.g rs7412)"),
                              h5("Gene name: gene of interest. The format is gene symbol (e.g APOE)")),
                         # select chromosome in case of manual scroll -- this has to be fixed and implemented
                         conditionalPanel(condition = "input.sel == 'Manual scroll'",
                              sliderInput(inputId = 'manual.chrom', label = "Chromosome", 
                                          value = 11, min=1, max=22, width = "60%"),
                              sliderInput(inputId = 'manual.pos', label = "Position", 
                                          15000000, min = 1, max = 200000000, width = "100%")),
                               
                         # maybe also good to put a submit button -- to avoid refreshes that are annoying and slow down application
                         div(submitButton(text = "Show locus!", icon("refresh"), width="200px"), style = "font-size:120%")))),
                    
                    # Space before main plot window
                    hr(),
                    
                    # Main plot window and side bar
                    fluidRow(  
                      column(4,
                             wellPanel(
                               h3("Visualization options"),
                               
                               # need a fluid row so that in case there are two GWAS to plot, a line is added for the colors
                               fluidRow(
                                 column(6,
                                        # choose color
                                        colourInput("col", "Select colour", "black"),
                                        
                                        # conditional panel for multiple points
                                        conditionalPanel(condition = "(input.gwas_neuro.length + input.gwas_cardio.length + input.gwas_immune.length + input.gwas_cancer.length + input.gwas_physio.length + input.gwas_from_file.length) > 1",
                                                         colourInput("col2", "Select color 2", "coral")),
                                        conditionalPanel(condition = "(input.gwas_neuro.length + input.gwas_cardio.length + input.gwas_immune.length + input.gwas_cancer.length + input.gwas_physio.length + input.gwas_from_file.length) > 2",
                                                         colourInput("col3", "Select color 3", "navy")),
                                        conditionalPanel(condition = "(input.gwas_neuro.length + input.gwas_cardio.length + input.gwas_immune.length + input.gwas_cancer.length + input.gwas_physio.length + input.gwas_from_file.length) > 3",
                                                         colourInput("col4", "Select color 4", "green")),
                                        conditionalPanel(condition = "(input.gwas_neuro.length + input.gwas_cardio.length + input.gwas_immune.length + input.gwas_cancer.length + input.gwas_physio.length + input.gwas_from_file.length) > 4",
                                                         colourInput("col5", "Select color 5", "yellow")))),
                               
                               # sliding window -- x.axis
                               sliderInput(inputId = "x", label = "Window (bp)", value = 30000, min = 1000, max = 1000000),
                               
                               # sliding widow -- y.axis
                               sliderInput(inputId = "y", label = "-log10 (P-value)", value = 9, min = 0, max = 100),
                               
                               # select type of plot -- lines or points
                               radioButtons(inputId = "ploType", label = "Select plot type", choices = c("Points", "Lines")),
                               
                               # conditional panel for line type -- adjusting smoothing parameters for density line
                               conditionalPanel(condition = "input.ploType == 'Lines'", sliderInput(inputId = 'sliding.window', label = "Number of bins", 
                                                            value = 30, min=10, max=50),
                                                sliderInput(inputId = "smooth", label = "Smoothing parameter",
                                                            value = 0.10, min=0, max=1)),
                               # Download button
                               h5("Download your plot"),
                               downloadButton('downloadPlot', 'Download Plot'),
                               div(style = "margin-top: 20px"),
                               div(style = "border-top: 2px solid #999999"),
                               
                               # for LD i also need a fluid row to allow the user to choose a separate variant to calculate LD from
                               h3("Linkage Disequilibrium (LD)"),
                               fluidRow(
                                 column(12,
                                        # select type of LD -- no LD, most significant in window, Input variant or Other (other makes appear a text panel)
                                        radioButtons(inputId = "Linkage", label = "LD options", choices = c("No LD", "Input variant", "Most significant in window",
                                                                                                            "Other variant")),
                                        # conditional panel for the choice of the population to use
                                        h5("Select population(s) to calculate LD"),
                                        div(checkboxGroupInput("pop_eur", "European", c("CEU" = "CEU", "FIN" = "FIN", "GBR" = "GBR", "IBS" = "IBS", "TSI" = "TSI", "ALL Sub-Pop" = "ALL_eur"), inline = T),  style = "font-size:80%"),
                                        div(checkboxGroupInput("pop_afr", "African", c("ACB" = "ACB", "ASW" = "ASW", "ESN" = "ESN", "GWD" = "GWD", "LWK" = "LWK", "MSL" = "MSL", "YRI" = "YRI", "ALL Sub-Pop" = "ALL_afr"), inline = T),  style = "font-size:80%"),
                                        div(checkboxGroupInput("pop_amr", "American", c("CLM" = "CLM", "MXL" = "MXL", "PEL" = "PEL", "PUR" = "PUR", "ALL Sub-Pop" = "ALL_amr"), inline = T),  style = "font-size:80%"),
                                        div(checkboxGroupInput("pop_eas", "East-Asian", c("CDX" = "CDX", "CHB" = "CHB", "CHS" = "CHS", "JPT" = "JPT", "KHV" = "KHV", "ALL Sub-Pop" = "ALL_eas"), inline = T),  style = "font-size:80%"),
                                        div(checkboxGroupInput("pop_sas", "South-Asian", c("BEB" = "BEB", "GIH" = "GIH", "ITU" = "ITU", "PJL" = "PJL", "STU" = "STU", "ALL Sub-Pop" = "ALL_sas"), inline = T),  style = "font-size:80%"),
                                        h6("** If you don't select any population, by default European samples of the 1000Genome project will be used. Please refer to the Help section for additional information about the populations reported above."),                                        
                                        div(style = "margin-top: 10px"),
                                        h6("** In case the top SNP is not in 1000Genome data, we take the second-most significant, then the third-most significant, etc."),
                                        # conditional panel
                                        conditionalPanel(condition = "input.Linkage == 'Other variant'",
                                                         textInput(inputId = "ld_other_snp", label = "Input variant [Chr:Pos]", value = "Type position...")))),
                               
                               # empty line with line draw before the table
                               h5("Download LD info"),
                               downloadButton("download_LDtable", "Download LD table"),
                               div(style = "margin-top: 20px"),
                               div(style = "border-top: 2px solid #999999"),
                               
                               # then put the choice of the structural variation dataset
                               h3("Structural variations (SV)"),
                               radioButtons(inputId = "strVar_inp", label = "Structural variants data", choices = c("Linthrost et al. 2019" = "jasper", "Chaisson et al. 2019" = "chaisson", "Audano et al 2019" = "audano", "All" = "all")),
                               h5("Download the SVs in the region"),
                               downloadButton("download_SVs", "Download All SVs"),
                               div(style = "margin-top: 20px"),
                               div(style = "border-top: 2px solid #999999"),
                               
                               # this should be a window where to put table of snps info
                               h3("Top SNPs info"),
                               div(tableOutput('table'),  style = "font-size:100%"),
                               h6("Chr: chromosome; Position: position according to genome version chosen; P-value: association p-value (in log10 scale); Study: input study name of the association"),
                               h5("Download the top associations"),
                               downloadButton("download_SNPs", "Download All associations"),
                               div(style = "margin-top: 20px"),
                               div(style = "border-top: 2px solid #999999"),
                               
                               # then add the section about eqtl
                               h3("eQTL associations"),
                               div(tableOutput('eQTL'),  style = "font-size:100%"),
                               h6("Locus: genomic position according to genome version chosen; A1/A2: the two alleles; MAF: minor allele frequency; Effect: normalized effect-size; P-value: association significance (in log10 scale); Gene: gene in eQTL with variant"),
                               h5("Download the top associations"),
                               downloadButton("download_eQTL", "Download Top eQTLs"),
                               div(style = "margin-top: 20px"),
                               div(style = "border-top: 2px solid #999999"),
                               
                               # then add section about cross-references
                               h4("Cross references"),
                               h6("If you specify a Gene or a SNP id, you will find here links to GeneCards and dbSNP."),
                               uiOutput(outputId = "genecards_link"),
                               uiOutput(outputId = "gwascat_link"),
                               uiOutput(outputId = "ld_hub_link"))), 
                      
                      # finally, this is the plotting window on the right
                      column(8, plotOutput(outputId = "hist", height = "2000px"))),
           
                    # final row with images of tudelft and aumc                    
                    fluidRow(
                      column(8,
                             h5("snpXplorer and AnnotateMe have been developed in collaboration with Amsterdam UMC and TU Delft."),
                             h5("Any suggestion or bug report is highly appreciated. Please email n.tesi@amsterdamumc.nl.")),
                      column(2, img(src='tudelft1.png', width="90%")),
                      column(2, img(src='amstUMC.jpg', width="90%")))),
           
           tabPanel("AnnotateMe",
                    # DNA image
                    fluidRow(column(12, img(src='dna.png', align="center", width="100%"))),
                    
                    # Overview of the tool
                    fluidRow(column(12, img(src='AnnotateMe.png', width="100%", align="center"))),

                    # Tutorial image
                    fluidRow(column(12, img(src='tutorial_annot.png', align="center", width="100%"))),
                    
                    hr(),
                    
                    # This within a new fluidRow
                    fluidRow(
                      # text area for the input variants
                      column(3,
                             textAreaInput(inputId = "snp_list", label = "List of variants", value = "Paste variants of interest (chr:pos -- chr  pos -- rsid)", 
                                           height = "350px")),
                      
                      # radio buttons for input type
                      column(2, 
                              radioButtons(inputId = "snp_list_type", label = "Input type", choices = c("chr:pos (1:12345678)", "chr pos (1 12345678)", 
                                                                                              "rsid (rs12345)")),
                              #radioButtons(inputId = "analysis_mode", label = "Analysis type", choices = c("Default (GO:BP)", "KEGG", "Reactome", "Wiki Pathways", "All"))),
                              checkboxGroupInput("analysis_mode", "Source Enrichment", c("Default (GO:BP)" = "default", "KEGG" = "KEGG", "Reactome" = "Reactome", "Wiki Pathways" = "wiki"), inline = F, selected = "default")),
                      # text area for email
                      column(3, 
                             textInput(inputId = "email", label = "E-mail", value = "Type your email address...", width = "400px"),
                             h6("Once analysis is done, you'll receive results by email."),
                             # maybe also good to put a submit button -- to avoid refreshes that are annoying and slow down application
                             submitButton(text = "Submit!", icon("refresh"), width="200px")),
                      
                      # dummy plot
                      column(4, plotOutput(outputId = "annot", height = "350px"))),

                    # final row with images of tudelft and aumc                    
                    fluidRow(
                      column(8,
                             h5("snpXplorer and AnnotateMe have been developed in collaboration with Amsterdam UMC and TU Delft."),
                             h5("Any suggestion or bug report is highly appreciated. Please email n.tesi@amsterdamumc.nl.")),
                      column(2, img(src='tudelft1.png', width="90%")),
                      column(2, img(src='amstUMC.jpg', width="90%")))),
           
           navbarMenu("More",
                      tabPanel("Studies and References",
                               wellPanel(h2("Studies and References", align="center"),
                                         h5("We are including more and more summary statistics.", align="center"),
                                         h5("For particular requests, please contact snpxplorer@gmail.com", align="center")),
                               
                               # references
                               h3("Available summary statistic within snpXplorer"),
                               tableOutput(outputId = 'table_info'),
                               div(style = "margin-top: 50px"),
                               h3("1000 Genome Project population codes"),
                               tableOutput(outputId = 'table_info_1000G'),
                               div(style = "margin-top: 50px"),
                               h3("Structural variants datasets"),
                               tableOutput(outputId = 'table_info_SVs'),
                               
                               # final row with images of tudelft and aumc                    
                               br(), 
                               br(),
                               br(),
                               fluidRow(
                                 column(8,
                                        h5("snpXplorer and AnnotateMe have been developed in collaboration with Amsterdam UMC and TU Delft."),
                                        h5("Any suggestion or bug report is highly appreciated. Please email n.tesi@amsterdamumc.nl.")),
                                 column(2, img(src='tudelft1.png', width="90%")),
                                 column(2, img(src='amstUMC.jpg', width="90%")))),
                      
                      tabPanel("Documentation and code",
                               wellPanel(h2("Documentation and code", align="center"),
                                         h5("Documentation can be viewed in the blox below.", align="center"),
                                         h5("Source code is available at Github", align="center")),
                               column(4, wellPanel(h5("Source code is freely available through our Github page."),
                                                   uiOutput("github"),
                                                   br(),
                                                   h5("snpXplorer can also be installed locally in your machine, but you will need to download summary statistics yourself."))),
                               column(8, uiOutput("pdf_doc_view")),
           
                               # final row with images of tudelft and aumc                    
                               br(), 
                               br(),
                               br(),
                               fluidRow(
                                 column(8,
                                        h5("snpXplorer and AnnotateMe have been developed in collaboration with Amsterdam UMC and TU Delft."),
                                        h5("Any suggestion or bug report is highly appreciated. Please email n.tesi@amsterdamumc.nl.")),
                                 column(2, img(src='tudelft1.png', width="90%")),
                                 column(2, img(src='amstUMC.jpg', width="90%")))),

                      tabPanel("Cite us",
                               wellPanel(h2("Cite us", align="center")),
                               column(4, wellPanel(h4("Please do not forget to cite us if you find snpXplorer useful for your research."),
                                                   br(),
                                                   h4("!! The preprint of snpXplorer is out"),
                                                   uiOutput(outputId = "biorxiv_link"),
                                                   br(),
                                                   h4("snpXplorer functional annotation was used in the following papers:"),
                                                   h5("Niccolo’ Tesi et al., Polygenic risk score of longevity predicts longer survival across an age-continuum. The Journals of Gerontology: Series A, , glaa289, https://doi.org/10.1093/gerona/glaa289"),
                                                   uiOutput(outputId = "longevity_ms"),
                                                   br(),
                                                   h5("Niccolo’ Tesi et al., The effect of Alzheimer's disease-associated genetic variants on longevity. merRxiv, https://doi.org/10.1101/2021.02.02.21250991"),
                                                   uiOutput(outputId = "rotation_ms"))),
                               
                               column(8, uiOutput("pdf_biorxiv_view")),
                               
                      
                               # final row with images of tudelft and aumc                    
                               br(), 
                               br(),
                               br(),
                               fluidRow(
                                 column(8,
                                        h5("snpXplorer and AnnotateMe have been developed in collaboration with Amsterdam UMC and TU Delft."),
                                        h5("Any suggestion or bug report is highly appreciated. Please email n.tesi@amsterdamumc.nl.")),
                                 column(2, img(src='tudelft1.png', width="90%")),
                                 column(2, img(src='amstUMC.jpg', width="90%")))),
                      
                      tabPanel("Help", 
                               wellPanel(
                                 h2("Help page", align="center"),
                                 h4("This page provides sample files to try out both the Exploration and the Annotation part, and a quick-start video.", align = "center")),
                               
                               # add wellpanel for the sample file of the exploration section
                               column(3, wellPanel(
                                 h4("Exploration section"),
                                 hr(),
                                 h5("For the exploration section, please provide a space- or tab-separated file with header."), 
                                 h5("Make sure there are at least chromosome number, position and p-value. snpXplorer will try to understand the columns from the file header."), 
                                 h5("Alternatively, you can use PLINK formatted files. Look at the examples below or download a ready-to-use dataset to try out."),
                                 br()
                               )),
                               
                                # add wellpanel for example file
                                column(2,
                                   h5("Example file", align="center"),
                                   hr(),
                                   h6("chr pos p"),
                                   h6("16 89659 0.19"),
                                   h6("16 90318 0.48"),
                                   h6("16 439873 0.28"),
                                   h6("16 904328 0.01"),
                                   h6(".. .. .."),
                                   hr(),
                                   # Button to download trial and plink-based files
                                   downloadButton("download_help_exploration", "Download sample file")
                                  ),
                               
                                 # add wellpanel for example file
                                 column(7,
                                   h5("Example PLINK file", align="center"),
                                   hr(),
                                   h6("#CHROM	POS	ID	REF	ALT	A1	A1_FREQ	A1_CASE_FREQ	A1_CTRL_FREQ	MACH_R2	TEST	OBS_CT	BETA	SE	Z_STAT	P"),
                                   h6("17	15802438	rs12449443	G	A	A	0.07	0.06	0.07	0.99	ADD	4191	0.01	0.16	0.05	0.96"),
                                   h6("17	29378199	rs57278847	T	C	C	0.01	0.01	0.02	0.31	ADD	4191	16.40	18.15	0.90	0.36"),
                                   h6("17	43286432	rs57847	T	C	C	0.02	0.11	0.42	0.31	ADD	4191	1.40	1.15	0.70	0.16"),
                                   h6("17	22318299	rs5724327	T	C	C	0.31	0.61	0.42	0.11	ADD	4191	6.40	58.12	0.60	0.26"),
                                   h6(".. .. .."),
                                   hr(),
                                   # Button to download trial and plink-based files
                                   downloadButton("download_help_exploration_plink", "Download sample PLINK file")
                                 ),
                               
                               fluidRow(),
                               hr(),
                               br(),
                               
                               # add wellpanel for the sample file of the annotation section
                               column(3, wellPanel(
                                 h4("Annotation section"),
                                 hr(),
                                 h5("The annotation section accepts any set of SNPs (possibly more than 1)."),
                                 h5("Multiple format are accepted: you can choose between rs-identifier or genomic position."),
                                 h5("Please see examples on the right or download a trial dataset."),
                                 br()
                               )),
                               
                               # add wellpanel for example file
                               column(3,
                                      h5("Example file #1 ~ chr:position", align="center"),
                                      hr(),
                                      h6("16:89659"),
                                      h6("16:90318"),
                                      h6("16:439873"),
                                      h6("16:904328"),
                                      h6(".. .. .."),
                                      hr(),
                                      # Button to download trial and plink-based files
                                      downloadButton("download_help_annotation1", "Download sample file")
                               ),
                               
                                       # add wellpanel for example file
                                column(3,
                                       h5("Example file #1 ~ chr position", align="center"),
                                       hr(),
                                       h6("16 89659"),
                                       h6("16 90318"),
                                       h6("16 439873"),
                                       h6("16 904328"),
                                       h6(".. .. .."),
                                       hr(),
                                       # Button to download trial and plink-based files
                                       downloadButton("download_help_annotation2", "Download sample file")
                              ),
                      
                                column(3,
                                        h5("Example file #3 ~ rsID", align="center"),
                                        hr(),
                                        h6("rs7412"),
                                        h6("rs1234"),
                                        h6("rs439873"),
                                        h6("rs904328"),
                                        h6(".. .. .."),
                                        hr(),
                                        # Button to download trial and plink-based files
                                        downloadButton("download_help_annotation3", "Download sample file")
                              ),
                      
                              fluidRow(),
                              hr(),
                            
                              # add wellpanel for the youtube videos
                              column(12, wellPanel(
                                h3("Video tutorials", align = "center"),
                                h4("Have a look at our quick start videos of snpXplorer.", align = "center"),
                                h4("There are tutorials about exploration section and annotation sections.", align = "center"),
                              )),
                              
                              column(6,
                                     h4("Quick Start tutorial", align="center"),
                                     hr(),
                                     uiOutput("video1"),
                                     hr(),
                              ),
                              
                              column(6,
                                     h4("Load your own file tutorial", align="center"),
                                     hr(),
                                     uiOutput("video2"),
                                     hr(),
                              ),
                              
                              column(6,
                                     h4("Functional annotation of your SNP list", align="center"),
                                     hr(),
                                     uiOutput("video3"),
                                     hr(),
                              ),

                              column(6,
                                     h4("Discover LD patterns", align="center"),
                                     hr(),
                                     uiOutput("video4"),
                                     hr(),
                              ),
                              
                               # final row with images of tudelft and aumc                    
                               br(),
                               fluidRow(
                                 column(8,
                                        h5("snpXplorer and AnnotateMe have been developed in collaboration with Amsterdam UMC and TU Delft."),
                                        h5("Any suggestion or bug report is highly appreciated. Please email n.tesi@amsterdamumc.nl.")),
                                 column(2, img(src='tudelft1.png', width="90%")),
                                 column(2, img(src='amstUMC.jpg', width="90%")))))
)
