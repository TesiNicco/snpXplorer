#################################
# libraries
suppressPackageStartupMessages({
  library(shiny)
  library(data.table)
  library(stringr)
  library(shinyWidgets)
  library(ggplot2)
  library(shinyBS)
  library(liftOver)
  library(colourpicker)
  library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
  library(rvest)
  library(stringr)
  library(plotrix)
})
# increase maximum size of uploaded file -- now is 600MB
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
                         h6("By default (i.e. without selecting any input dataset, snpXplorer will load IGAP study (chr16)."),
                         # select which gwas as input -- try to do based on disease type
                         checkboxGroupInput("gwas_neuro", "Neurological traits", c("Alzheimer_million" = "Alzheimer_million", "GR@ACE" = "GR@ACE", "IGAP" = "IGAP", "proxy_AD" = "proxy_AD", "Autism" = "Autism", "Depression" = "Depression", "Ventricular volume" = "LVV"), inline = T),
                         checkboxGroupInput("gwas_cardio", "Cardiovascular traits", c("CAD" = "CAD", "SBP" = "SBP", "BMI" = "BMI", "Diabetes" = "Diabetes"), inline = T),
                         checkboxGroupInput("gwas_immune", "Immune-related traits", c("COVID" = "COVID", "Lupus" = "Lupus", "Inflammation" = "Inflammation", "Asthma" = "Asthma"), inline = T),
                         checkboxGroupInput("gwas_cancer", "Cancer-related traits", c("Breast_cancer" = "Breast_cancer", "Myeloproliferative" = "Myeloproliferative", "Prostate" = "Prostate", "Lung" = "Lung", "Leukemia" = "Leukemia"), inline = T),
                         checkboxGroupInput("gwas_physio", "Physiological traits", c("Multivariate_Longevity" = "Multivariate_Longevity", "UKBaging" = "UKBaging", "Height" = "Height", "Education" = "Education", "Bone_density" = "Bone_density", "Vitamin_D" = "Vitamin_D"), inline = T)),
                         #bsTooltip(id = "gwas_neuro", title = "IGAP: International Genomics of Alzheimer's Project\nproxy_AD: familial Alzheimer's disease", placement = "right", trigger = "hover"),
                         bsPopover("gwas_neuro", title='Neurological traits', content='Alzheimer_million: GWAS of Alzheimer including 1 million individuals -- GR@ACE: Genomic Research @ Fundacio ACE (GWAS of Alzheimers disease) -- IGAP: International Genomics of Alzheimers disease (AD, case-control GWAS) -- Proxy_AD: family history of AD', placement="right",trigger="hover", options = list(container = "body")),
                         bsPopover("gwas_cardio", title='Cardiovascular traits', content='CAD: Coronary Artery Disease -- SBP: Systolic Blood Pressure -- BMI: Body-mass index', placement="right",trigger="hover", options = list(container = "body")),
                         bsPopover("gwas_immune", title='Immune-related traits', content='COVID: severe CODIV-19 reactions -- Lupus: Lupus Erithematosous', placement="right",trigger="hover", options = list(container = "body")),
                         bsPopover("gwas_physio", title='Physiological traits', content='Multivariate_Longevity: Multivariate analysis of Lifespan, Longevity and Healthspan -- UKBaging: parental longevity', placement="right",trigger="hover", options = list(container = "body")),

                      # next put the genome options and the personal file upload
                      column(4,
                         br(),
                         # select which gwas as input from the user only
                         checkboxGroupInput("gwas_from_file", "Your input", c("Your input" = "toLoad"), inline = T),
                         # select user input file
                         fileInput(inputId = "inp_f", label = "Add file", multiple = F, placeholder = "Loading file", width="300px"),
                         # help message
                         bsPopover("gwas_from_file", title='Your file', content='Please make sure data has at least chromosome, position and pvalue columns. See More/Help to download sample files and More/Documentation for additional information.', placement="right",trigger="hover", options = list(container = "body")),
                         #div(style = "font-size: 10px; padding: 0px 0px; margin-top:-2em", h6("If you upload your file, make sure data has at least chromosome, position and pvalue columns.
                         #See More/Help to download sample files for Exploration section. See More/Documentation for additional information.")),
                         div(style = "border-top: 2px solid #999999"),
                         h3("Reference Genome"),
                         # select genome options
                         radioButtons(inputId = "genV", label = "Select genome version", choices = c("GRCh37 (hg19)", "GRCh38 (hg38)")),
                         bsPopover("genV", title='Reference Genome', content='By default, GRCh37 (hg19) is used. snpXplorer uses liftOver tool to change genomic coordinates between reference genome versions.', placement="right",trigger="hover", options = list(container = "body")),
                         #h6("By default, GRCh37 (hg19) is used. snpXplorer uses liftOver tool to change genomic coordinates between reference genome versions."),
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
                         br(),
                         div(submitButton(text = "Show Results!", icon("dna"), width="250px"), style = "font-size:150%"),
                         h5("Please remeber to press the 'Show Results!' button to update the plot.")))),

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
                                        colourpicker::colourInput("col", "Select colour", "black"),
                                        bsPopover("col", title='Colors', content='You can use the colour picker to change the color of the plotted studies. Press Show Results! to update.', placement="right",trigger="hover", options = list(container = "body")),

                                        # conditional panel for multiple points
                                        conditionalPanel(condition = "(input.gwas_neuro.length + input.gwas_cardio.length + input.gwas_immune.length + input.gwas_cancer.length + input.gwas_physio.length + input.gwas_from_file.length) > 1",
                                                         colourpicker::colourInput("col2", "Select color 2", "coral")),
                                        conditionalPanel(condition = "(input.gwas_neuro.length + input.gwas_cardio.length + input.gwas_immune.length + input.gwas_cancer.length + input.gwas_physio.length + input.gwas_from_file.length) > 2",
                                                         colourpicker::colourInput("col3", "Select color 3", "navy")),
                                        conditionalPanel(condition = "(input.gwas_neuro.length + input.gwas_cardio.length + input.gwas_immune.length + input.gwas_cancer.length + input.gwas_physio.length + input.gwas_from_file.length) > 3",
                                                         colourpicker::colourInput("col4", "Select color 4", "green")),
                                        conditionalPanel(condition = "(input.gwas_neuro.length + input.gwas_cardio.length + input.gwas_immune.length + input.gwas_cancer.length + input.gwas_physio.length + input.gwas_from_file.length) > 4",
                                                         colourpicker::colourInput("col5", "Select color 5", "yellow")))),

                               # sliding window -- x.axis
                               sliderInput(inputId = "x", label = "Window (bp)", value = 30000, min = 1000, max = 1000000),
                               bsPopover("x", title='Window x-size', content='You can use the slider to adjust x-axis of the main plot. Press Show Results! to update.', placement="right",trigger="hover", options = list(container = "body")),

                               # sliding widow -- y.axis
                               sliderInput(inputId = "y", label = "-log10 (P-value)", value = 9, min = 0, max = 100),
                               bsPopover("y", title='Window y-size', content='You can use the slider to adjust y-axis of the main plot. Press Show Results! to update.', placement="right",trigger="hover", options = list(container = "body")),

                               # select type of plot -- lines or points
                               radioButtons(inputId = "ploType", label = "Select plot type", choices = c("Points", "Lines")),
                               bsPopover("ploType", title='Plot type', content='Points: each SNP is a dot -- Lines: draw p-value density lines. Press Show Results! to update.', placement="right",trigger="hover", options = list(container = "body")),

                               # show/hide recombination rates
                               radioButtons(inputId = "recomb_yn", label = "Add recombination rates", choices = c("Yes", "No")),
                               bsPopover("recomb_yn", title='Show recombination rates', content='Show (Yes) or hide (No) recombination rates from the main plot. Press Show Results! to update.', placement="right",trigger="hover", options = list(container = "body")),

                               # show/hide increasing dot-sizes
                               radioButtons(inputId = "dotSize_yn", label = "Increase dot-size as significance increases", choices = c("Yes", "No")),
                               bsPopover("dotSize_yn", title='Dot-size as function of significance', content='Increase size of dots depending on significance (Yes). Press Show Results! to update.', placement="right",trigger="hover", options = list(container = "body")),

                               # conditional panel for line type -- adjusting smoothing parameters for density line
                               conditionalPanel(condition = "input.ploType == 'Lines'", sliderInput(inputId = 'sliding.window', label = "Number of bins",
                                                            value = 30, min=10, max=50),
                                                sliderInput(inputId = "smooth", label = "Smoothing parameter",
                                                            value = 0.10, min=0, max=1),
                                                bsPopover("sliding.window", title='Number of bins', content='Will affect p-value density lines. Smaller values generate more precise p-value density lines. Press Show Results! to update.', placement="right",trigger="hover", options = list(container = "body")),
                                                bsPopover("smooth", title='Smoothing parameter', content='Will affect p-value density lines. Smaller values generate more precise p-value density lines. Press Show Results! to update.', placement="right",trigger="hover", options = list(container = "body")),
                               ),

                               # repeat the main button to update the plot
                               div(style = "margin-top: 20px"),
                               div(submitButton(text = "Show Results!", icon("dna"), width="250px"), style = "font-size:150%"),
                               div(style = "margin-top: 20px"),
                               div(style = "border-top: 2px solid #999999"),

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
                                        bsPopover("Linkage", title='LD options', content='Define if LD should be shown (only available when 1 study is selected). If LD should be shown, please indicate in which populations (from 1000Genome Project) LD should be calculated. ', placement="right",trigger="hover", options = list(container = "body")),

                                        # conditional panel for the choice of the population to use
                                        h5("Select population(s) to calculate LD"),
                                        div(checkboxGroupInput("pop_eur", "European", c("CEU" = "CEU", "FIN" = "FIN", "GBR" = "GBR", "IBS" = "IBS", "TSI" = "TSI", "ALL Sub-Pop" = "ALL_eur"), inline = T),  style = "font-size:80%"),
                                        bsPopover("pop_eur", title='European populations', content='CEU: Utah Residents with Northern and Western European Ancestry -- FIN: Finnish -- GBR: British in England and Scotland -- IBS: Iberian Population in Spain -- TSI: Toscani in Italy', placement="right",trigger="hover", options = list(container = "body")),
                                        div(checkboxGroupInput("pop_afr", "African", c("ACB" = "ACB", "ASW" = "ASW", "ESN" = "ESN", "GWD" = "GWD", "LWK" = "LWK", "MSL" = "MSL", "YRI" = "YRI", "ALL Sub-Pop" = "ALL_afr"), inline = T),  style = "font-size:80%"),
                                        bsPopover("pop_afr", title='African populations', content='ACB: African Caribbeans in Barbados -- ASW: Americans of African Ancestry in SW USA -- ESN: Esan in Nigeria -- GWD: Gambian in Western Divisions in the Gambia -- LWK: Luhya in Webuye, Kenya -- MSL: Mende in Sierra Leone -- YRI: Yoruba in Ibadan, Nigeria', placement="right",trigger="hover", options = list(container = "body")),
                                        div(checkboxGroupInput("pop_amr", "American", c("CLM" = "CLM", "MXL" = "MXL", "PEL" = "PEL", "PUR" = "PUR", "ALL Sub-Pop" = "ALL_amr"), inline = T),  style = "font-size:80%"),
                                        bsPopover("pop_amr", title='American populations', content='CLM: Colombians from Medellin, Colombia -- MXL: Mexican Ancestry from Los Angeles USA -- PEL: Peruvians from Lima, Peru -- PUR: Puerto Ricans from Puerto Rico', placement="right",trigger="hover", options = list(container = "body")),
                                        div(checkboxGroupInput("pop_eas", "East-Asian", c("CDX" = "CDX", "CHB" = "CHB", "CHS" = "CHS", "JPT" = "JPT", "KHV" = "KHV", "ALL Sub-Pop" = "ALL_eas"), inline = T),  style = "font-size:80%"),
                                        bsPopover("pop_eas", title='East-Asian populations', content='CDX: Chinese Dai in Xishuangbanna, China -- CHB: Han Chinese in Beijing, China -- CHS: Southern Han Chinese -- JPT: Japanese in Tokyo, Japan -- KHV: Kinh in Ho Chi Minh City, Vietnam', placement="right",trigger="hover", options = list(container = "body")),
                                        div(checkboxGroupInput("pop_sas", "South-Asian", c("BEB" = "BEB", "GIH" = "GIH", "ITU" = "ITU", "PJL" = "PJL", "STU" = "STU", "ALL Sub-Pop" = "ALL_sas"), inline = T),  style = "font-size:80%"),
                                        bsPopover("pop_sas", title='South-Asian populations', content='BEB: Bengali from Bangladesh -- GIH: Gujarati Indian from Houston, Texas -- ITU: Indian Telugu from the UK -- PJL: Punjabi from Lahore, Pakistan -- STU: Sri Lankan Tamil from the UK', placement="right",trigger="hover", options = list(container = "body")),
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
                               bsPopover("strVar_inp", title='Structural variants', content='Third-generation sequencing (long-read sequencing) studies estimating the major Structural Variants (SV) alleles.', placement="right",trigger="hover", options = list(container = "body")),
                               h5("Download the SVs in the region"),
                               downloadButton("download_SVs", "Download All SVs"),
                               div(style = "margin-top: 20px"),
                               div(style = "border-top: 2px solid #999999"),

                               # this should be a window where to put table of snps info
                               h3("Top SNPs info"),
                               div(tableOutput('table'),  style = "font-size:75%"),
                               bsPopover("table", title='Top SNPs table', content='Most significant SNP-associations in the main plot interface. Chr: chromosome; Position: genomic position depending on reference genome choosen; ID: variant identifier (missing if variant is not in 1000Genome or too rare); -log10(P): association p-value (in -log10 scale); Study: name of the study; Alleles: reference and alternative alleles, respectively (from 1000Genome data).', placement="right",trigger="hover", options = list(container = "body")),
                               #h6("Chr: chromosome; Position: position according to genome version chosen; P-value: association p-value (in log10 scale); Study: input study name of the association"),
                               h5("Download the top associations"),
                               downloadButton("download_SNPs", "Download All associations"),
                               div(style = "margin-top: 20px"),
                               div(style = "border-top: 2px solid #999999"),

                               # then add the section about eqtl
                               h3("eQTL associations"),
                               h6("You can add additional tissues by selecting them on the box below. Selecting All_tissues will display results of all GTEx tissues. It may take some time in case you select many tissues and the window-size is big, so please be patient."),
                               multiInput(inputId = "gtex_type_tb", label = "GTEx tissues", choices = c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary",
                                                                                                     "Artery_Tibial", "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia",
                                                                                                     "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus",
                                                                                                     "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1",
                                                                                                     "Brain_Substantia_nigra", "Breast_Mammary_Tissue", "Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes",
                                                                                                     "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis",
                                                                                                     "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Kidney_Cortex", "Liver", "Lung", "Minor_Salivary_Gland", "Muscle_Skeletal",
                                                                                                     "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg",
                                                                                                     "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood", "All_tissues"), selected = "Whole_Blood",
                                          options = list(enable_search = TRUE, non_selected_header = "Choose between:", selected_header = "You have selected:")),
                               div(tableOutput('eQTL'),  style = "font-size:75%"),
                               bsPopover("eQTL", title='Top eQTLs table', content='Locus: genomic position according GRCh38 (hg38); A1/A2: SNP alleles; MAF: minor allele frequency; Effect: normalized effect-size as taken from GTEx; P-value: association significance (in -log10 scale); Gene: gene in eQTL with SNP.', placement="right",trigger="hover", options = list(container = "body")),
                               #h6("Locus: genomic position according to genome version chosen; A1/A2: the two alleles; MAF: minor allele frequency; Effect: normalized effect-size; P-value: association significance (in log10 scale); Gene: gene in eQTL with variant"),
                               div(style = "margin-top: 20px"),
                               div(submitButton(text = "Update eQTL table", icon("dna"), width="250px"), style = "font-size:150%"),
                               div(style = "margin-top: 20px"),
                               h5("Download the top associations"),
                               downloadButton("download_eQTL", "Download All eQTLs"),
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

           tabPanel("Annotation",
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
                      column(2,
                             textAreaInput(inputId = "snp_list", label = "List of variants", value = "Paste variants of interest (chr:pos -- chr  pos -- rsid)",
                                           height = "700px"),
                             bsPopover("snp_list", title='Annotation input', content='Please insert your SNPs of interest here: make sure the format is correct and select the adequate format in the Input type box.', placement="right",trigger="hover", options = list(container = "body"))),


                      # radio buttons for input type
                      column(4,
                              radioButtons(inputId = "snp_list_type", label = "Input type", choices = c("chr:pos (1:12345678)", "chr pos (1 12345678)",
                                                                                              "rsid (rs12345)"), selected = "rsid (rs12345)"),
                              radioButtons("analysis_type", "Analysis type", c("SNP-gene annotation", "Gene-set enrichment analysis")),
                              bsPopover("snp_list_type", title='Input type', content='For compatibility, rsid is the preferred input. If your input is not rsid, please specify the Reference version.', placement="right",trigger="hover", options = list(container = "body")),
                              radioButtons("snp_list_reference", "Reference Genome", choices = c("GRCh37 (hg19)", "GRCh38 (hg38)")),
                              bsPopover("snp_list_reference", title='Reference Genome', content='If your input is different from rsid, then select the correct Reference genome version.', placement="right",trigger="hover", options = list(container = "body")),
                              multiInput(inputId = "gtex_type", label = "GTEx tissues", choices = c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary",
                                                                                                   "Artery_Tibial", "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia",
                                                                                                   "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus",
                                                                                                   "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1",
                                                                                                   "Brain_Substantia_nigra", "Breast_Mammary_Tissue", "Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes",
                                                                                                   "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis",
                                                                                                   "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Kidney_Cortex", "Liver", "Lung", "Minor_Salivary_Gland", "Muscle_Skeletal",
                                                                                                   "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg",
                                                                                                   "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood", "All_tissues"), selected = "Whole_Blood",
                                        options = list(enable_search = TRUE, non_selected_header = "Choose between:", selected_header = "You have selected:"), width = "500px"),
                              h6("By default, only Blood is considered. Adding multiple tissues may impact your analysis as a larger number of genes may be associated with each SNP."),
                              bsPopover("analysis_type", title='Analysis type', content='SNP-gene annotation: only performs SNP-to-gene mapping (max 10.000 SNPs allowed). -- Gene-set enrichment analysis: performs SNP-to-gene mapping and gene-set enrichment analysis (max 1.000 SNPs allowed)', placement = "right", trigger = "hover", options = list(container = "body")),
                              checkboxGroupInput("analysis_mode", "Source Enrichment", c("Default (GO:BP)" = "default", "KEGG" = "KEGG", "Reactome" = "Reactome", "Wiki Pathways" = "wiki"), inline = F, selected = "default"),
                              bsPopover("analysis_mode", title='Source Enrichment', content='Please select here the gene-set databases for the gene-set enrichment analysis.', placement="right",trigger="hover", options = list(container = "body"))),
                      # text area for email
                      column(3,
                             textInput(inputId = "email", label = "E-mail", value = "Type your email address...", width = "400px"),
                             bsPopover("email", title='Your e-mail address', content='We do not do anything with it...apart from sending your annotation results. :)', placement="right",trigger="hover", options = list(container = "body")),
                             h6("Once analysis is done, you'll receive results by email."),
                             h6("Please check your spam folder to get your snpXplorer results."),
                             # maybe also good to put a submit button -- to avoid refreshes that are annoying and slow down application
                             submitButton(text = "Submit!", icon("paper-plane"), width="250px")),

                      # dummy plot
                      column(3, plotOutput(outputId = "annot", height = "350px"))),

                    # final row with images of tudelft and aumc
                    fluidRow(
                      column(8,
                             h5("snpXplorer and AnnotateMe have been developed in collaboration with Amsterdam UMC and TU Delft."),
                             h5("Any suggestion or bug report is highly appreciated. Please email n.tesi@amsterdamumc.nl.")),
                      column(2, img(src='tudelft1.png', width="90%")),
                      column(2, img(src='amstUMC.jpg', width="90%")))),

           navbarMenu("More",
                      tabPanel("Studies and References",
                               # DNA image
                               fluidRow(column(12, img(src='dna.png', align="center", width="100%"))),

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
                               # DNA image
                               fluidRow(column(12, img(src='dna.png', align="center", width="100%"))),

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
                               # DNA image
                               fluidRow(column(12, img(src='dna.png', align="center", width="100%"))),

                               wellPanel(h2("Cite us", align="center")),
                               column(4, wellPanel(h4("Please do not forget to cite us if you find snpXplorer useful for your research."),
                                                   br(),
                                                   h4("!! snpXplorer was published in Nucleic Acid Research"),
                                                   uiOutput(outputId = "biorxiv_link"),
                                                   br(),
                                                   h4("snpXplorer functional annotation was used in the following papers:"),
                                                   h5("Niccolo’ Tesi et al., Polygenic risk score of longevity predicts longer survival across an age-continuum. The Journals of Gerontology: Series A, , glaa289, https://doi.org/10.1093/gerona/glaa289"),
                                                   uiOutput(outputId = "longevity_ms"),
                                                   br(),
                                                   h5("Niccolo’ Tesi et al., The effect of Alzheimer's disease-associated genetic variants on longevity. medRxiv, https://doi.org/10.1101/2021.02.02.21250991"),
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
                               # DNA image
                               fluidRow(column(12, img(src='dna.png', align="center", width="100%"))),

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
                                 column(2, img(src='amstUMC.jpg', width="90%"))))),

           navbarMenu("Bug Report",
                      tabPanel("Report a bug",
                               # DNA image
                               fluidRow(column(12, img(src='dna.png', align="center", width="100%"))),
                               wellPanel(
                                 h2("Bug Report page", align="center")),
                                 column(12,
                                    h4("snpXplorer remembers the settings of your last run. In case you found a problem in the last run, you can now directly report the bug to snpXplorer."),
                                    br(),
                                    tableOutput(outputId = 'table_report_last'),
                                    textInput(inputId = "bug_comments", label = "What did go wrong?", value = "Please add any relevant comment here...", width="350px"),
                                    div(style = "margin-top: 20px"),
                                    div(style = "border-top: 2px solid #999999"),
                                    div(style = "margin-top: 20px"),
                                    radioButtons(inputId = "confirm_bug", label = "Do you want to send this bug report?", choices = c("Yes", "No"), selected = character(0)),
                                    div(submitButton(text = "Report bug!", icon("bug"), width="200px"), style = "font-size:120%"),
                                    column(4, plotOutput(outputId = "bugss", height = "350px"))))),

           tabPanel("What's new",
                    # DNA image
                    fluidRow(column(12, img(src='dna.png', align="center", width="100%"))),
                    wellPanel(
                      h2("What's new", align="center")),
                      column(12,
                             h4("24-12-21: GWAS catalog release was updated to the latest available."),
                             h4("11-11-21: New GWAS added: Largest GWAS of Alzheimer's disease was added."),
                             h4("11-11-21: New GWAS added: Multivariate analysis of Longevity, Healthspan and Lifespan was added."),
                             h4("10-11-21: New Annotation analysis added: possibility to do only SNP-gene mapping without gene-set enrichment analysis. This allows the user to annotate up to 10,000 SNPs."),
                             h4("25-08-21: CADD v1.6 (the most updated) is now used for SNP annotation."),
                             h4("25-08-21: To cope with rare SNPs, we now use data from all individuals of the 1000Genome Project (before, only European individuals and common SNPs were recognized)."),
                             h4("25-07-21: New GWAS added: GR@ACE GWAS of Alzheimer's disease was added.")))
)
