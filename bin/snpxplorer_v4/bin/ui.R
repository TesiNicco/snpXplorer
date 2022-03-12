## Libraries
  suppressPackageStartupMessages({
    library(shiny)
    library(plotly)
    library(shinyWidgets)
    library(shinyBS)
    library(colourpicker)
  })

## User Interface
  navbarPage("snpXplorer", inverse = T,
    tabPanel("Exploration",             # Exploration panel
      fluidRow(column(12, img(src='dna.png', align="center", width="100%"))),             # place DNA image on top of everything -- this is in a single row 
      fluidRow(column(12, img(src='tutorial.png', align="center", width="100%"))),        # then place in another row with the tutorial information
      hr(),
      fluidRow(                                                                         # initialize a single row -- in this space there will be input selection, reference genome, browsing options
        column(3,                                                                       # first column will be about the available input data to visualize (datasets or own dataset)
          wellPanel(
            div(tags$h2("Available GWAS Menu"), align = "center"),                                               # name of the dropdown menu
            div(dropdown(                                                                     # dropdown menu here
              style = "unite", circle = F, size = "lg", icon = icon("database"), status = "default", width = "500px", animate = animateOptions(enter = animations$fading_entrances$fadeInLeftBig, exit = animations$fading_exits$fadeOutLeftBig),
              tooltip = tooltipOptions(title = "Click to see available options!"),
              h6("Please check up to 5 boxes at the same time."),                         # additional text here and next line
              h6("By default (i.e. without selecting any input dataset, snpXplorer will load IGAP study (chr16)."),
              prettyCheckboxGroup(inputId = "neurol_gwas", label = "Neurological GWAS:", choices = c("Alzheimer Million" = "Alzheimer_million", "Alzheimer GR@ACE" = "GR@ACE", "Alzheimer IGAP" = "IGAP", "Alzheimer Proxy" = "proxy_AD", "Autism", "Depression", "Ventricular volume" = "LVV"), inline = TRUE, status = "default", fill = TRUE),
              prettyCheckboxGroup(inputId = "cardio_gwas", label = "Cardiovascular GWAS:", choices = c("CAD", "SBP", "BMI", "Diabetes"), inline = TRUE, status = "danger", fill = TRUE),
              prettyCheckboxGroup(inputId = "immune_gwas", label = "Immunological GWAS:", choices = c("COVID", "Lupus", "Inflammation", "Asthma"), inline = TRUE, status = "success", fill = TRUE),
              prettyCheckboxGroup(inputId = "cancer_gwas", label = "Oncological GWAS:", choices = c("Breast cancer" = "Breast_cancer", "Myeloproliferative", "Prostate", "Lung cancer" = "Lung", "Leukemia"), inline = TRUE, status = "warning", fill = TRUE),
              prettyCheckboxGroup(inputId = "physio_gwas", label = "Physiological GWAS:", choices = c("Multivariate Longevity" = "Multivariate_Longevity", "Parental Longevity" = "UKBaging", "Height", "Education", "Bone density" = "Bone_density", "Vitamin D" = "Vitamin_D"), inline = TRUE, status = "info", fill = TRUE),
              bsPopover("neurol_gwas", title='Neurological GWAS', content='Alzheimer Million: GWAS of Alzheimer including 1 million individuals -- Alzheimer GR@ACE: Genomic Research @ Fundacio ACE (GWAS of Alzheimers disease) -- Alzheimer IGAP: International Genomics of Alzheimers disease (AD, case-control GWAS) -- Alzheimer Proxy: family history of AD', placement="right",trigger="hover", options = list(container = "body")),
              bsPopover("cardio_gwas", title='Cardiovascular GWAS', content='CAD: Coronary Artery Disease -- SBP: Systolic Blood Pressure -- BMI: Body-mass index', placement="right",trigger="hover", options = list(container = "body")),
              bsPopover("immune_gwas", title='Immunological GWAS', content='COVID: severe CODIV-19 reactions -- Lupus: Lupus Erithematosous', placement="right",trigger="hover", options = list(container = "body")),
              bsPopover("physio_gwas", title='Physiological GWAS', content='Multivariate Longevity: Multivariate analysis of Lifespan, Longevity and Healthspan', placement="right",trigger="hover", options = list(container = "body")),
            ), align = "center"),
          ),
        ),
        column(3,                                                                       # second column is for the input file from the user
          wellPanel(
            div(tags$h2("Add your GWAS"), align = "center"),
            div(dropdown(
              style = "unite", size = "lg", icon = icon("upload"), status = "danger", width = "500px", animate = animateOptions(enter = animations$fading_entrances$fadeInLeftBig, exit = animations$fading_exits$fadeOutLeftBig),
              tooltip = tooltipOptions(title = "Click to see available options!"),
              h3("Visualize your own GWAS data"),
              prettyCheckboxGroup(inputId = "owned_gwas", label = "Visualize your own GWAS:", choices = c("Show your data"), inline = TRUE, status = "primary", fill = TRUE),
              fileInput(inputId = "upload_file", label = "Upload your data:", multiple = F, placeholder = "Loading file", width="300px"),
              h5('Please make sure data has at least chromosome, position and p-value columns. A header is required. See More/Help to download sample files and More/Documentation for additional information.'),
            ), align = "center"),
          ),
        ),
        column(3,                                                                       # third column is for the reference genome
          wellPanel(
            div(tags$h2("Reference Genome"), align = "center"),                                              # name of the dropdown menu
            div(dropdown(
              style = "unite", size = "lg", icon = icon("asterisk"), status = "warning", width = "500px", animate = animateOptions(enter = animations$fading_entrances$fadeInRightBig, exit = animations$fading_exits$fadeOutRightBig),
              tooltip = tooltipOptions(title = "Click to see available options!"),
              h3("Choose your reference genome"),
              prettyRadioButtons(inputId = "reference_gen", label = "Available versions:", choices = c("GRCh37 (hg19)", "GRCh38 (hg38)"), inline = TRUE, status = "info", fill = TRUE),
              h5('By default, GRCh37 (hg19) is used. snpXplorer uses liftOver tool to change genomic coordinates between reference genome versions.'),
            ), align = "center"),
          ),
        ),
        column(3,                                                                       # fourth column is for the browsing options
          wellPanel(
            div(tags$h2("Browsing options"), align = "center"),                                              # name of the dropdown menu
            div(dropdown(
              style = "unite", size = "lg", icon = icon("search"), status = "success", width = "400px", animate = animateOptions(enter = animations$fading_entrances$fadeInRightBig, exit = animations$fading_exits$fadeOutRightBig),
              tooltip = tooltipOptions(title = "Click to see available options!"),
              h3("Choose how to browse the genome"),
              prettyRadioButtons(inputId = "browsing_options", label = "Available options:", choices = c("Locus, Gene, RsID", "Manual Scroll"), icon = icon("check"), bigger = TRUE, status = "warning", animation = "jelly"),
              conditionalPanel(condition = "input.browsing_options == 'Locus, Gene, RsID'",            # conditional panel for the 'Locus, Gene, RsID mode'
                  textInput(inputId = "target", label = "Insert your target region", value = "Type locus, rsID or gene name", width="400px"),
                  h5("Locus: genomic location of interest. The format is chromosome:position (e.g 1:12345678)"),
                  h5("RsID: variant of interest. The format is rsXXXX (e.g rs7412)"),
                  h5("Gene name: gene of interest. The format is gene symbol (e.g APOE)")),
              conditionalPanel(condition = "input.browsing_options == 'Manual Scroll'",                # conditional panel for the 'Manual scroll'
                  sliderInput(inputId = 'manual_chrom', label = "Chromosome", value = 11, min=1, max=22, width = "100%"),
                  sliderInput(inputId = 'manual_pos', label = "Position", 15000000, min = 1, max = 250000000, width = "100%")),
            ), align = "center"),
          ),
        ),
      ),
      hr(),
      fluidRow(
        column(3,
          wellPanel(style = "background: #f2f2f2",
            column(6,
              div(tags$h3("Graphical options"), align = "center"),
              div(dropdown(
                style = "unite", size = "lg", icon = icon("gear"), status = "success", width = "500px", animate = animateOptions(enter = animations$fading_entrances$fadeInLeftBig, exit = animations$fading_exits$fadeOutLeftBig),
                tooltip = tooltipOptions(title = "Click to see available options!"),
                h3("Customize your plot"),
                colourpicker::colourInput("col1", "Select colour", "blue"),             # color selection for the different input gwases
                conditionalPanel(condition = "(input.neurol_gwas.length + input.cardio_gwas.length + input.immune_gwas.length + input.cancer_gwas.length + input.physio_gwas.length + input.owned_gwas.length) > 1",
                  colourpicker::colourInput("col2", "Select color 2", "coral")),
                conditionalPanel(condition = "(input.neurol_gwas.length + input.cardio_gwas.length + input.immune_gwas.length + input.cancer_gwas.length + input.physio_gwas.length + input.owned_gwas.length) > 2",
                  colourpicker::colourInput("col3", "Select color 3", "black")),
                conditionalPanel(condition = "(input.neurol_gwas.length + input.cardio_gwas.length + input.immune_gwas.length + input.cancer_gwas.length + input.physio_gwas.length + input.owned_gwas.length) > 3",
                  colourpicker::colourInput("col4", "Select color 4", "green")),
                conditionalPanel(condition = "(input.neurol_gwas.length + input.cardio_gwas.length + input.immune_gwas.length + input.cancer_gwas.length + input.physio_gwas.length + input.owned_gwas.length) > 4",
                  colourpicker::colourInput("col5", "Select color 5", "yellow")),
                bsPopover("col", title='Colors', content='You can use the colour picker to change the color of the plotted studies.', placement="right",trigger="hover", options = list(container = "body")),
                sliderInput(inputId = "x_axis", label = "Window (bp)", value = 50000, min = 1000, max = 1000000, width = "100%"),        # adjust the width of the x-axis
                bsPopover("x_axis", title='Window size', content='You can use the slider to adjust x-axis of the main plot.', placement="right",trigger="hover", options = list(container = "body")),
                sliderInput(inputId = "y_axis", label = "-log10 (P-value)", value = 9, min = 5, max = 100, width = "100%"),              # adjust the width of the y-axis
                bsPopover("y_axis", title='Window y-size', content='You can use the slider to adjust y-axis of the main plot.', placement="right",trigger="hover", options = list(container = "body")),
                radioGroupButtons(inputId = "plot_type", label = "Plot type:", choices = c("Scatter", "Density"), status = "primary", checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon"))),
                bsPopover("ploType", title='Plot type', content='Points: each SNP is a dot -- Lines: draw p-value density lines.', placement="right",trigger="hover", options = list(container = "body")),
                radioGroupButtons(inputId = "recomb_yn", label = "Add recombination rates:", choices = c("Yes", "No"), status = "warning", checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon"))),
                bsPopover("recomb_yn", title='Show recombination rates', content='Show (Yes) or hide (No) recombination rates from the main plot.', placement="right",trigger="hover", options = list(container = "body")),
                radioGroupButtons(inputId = "exons", label = "Add exons:", choices = c("No", "Yes"), status = "danger", checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon"))),
                bsPopover("exons", title='Add gene exons', content='Show (Yes) or hide (No) exons in each gene. Note that the waiting time will be slightly higher.', placement="right",trigger="hover", options = list(container = "body")),
              ), align = "center"),
            ),
            column(6,
              div(tags$h3("Linkage Disequilibrium"), align = "center"),
              div(dropdown(
                style = "unite", size = "lg", icon = icon("link"), status = "royal", width = "500px", animate = animateOptions(enter = animations$fading_entrances$fadeInLeftBig, exit = animations$fading_exits$fadeOutLeftBig),
                tooltip = tooltipOptions(title = "Click to see available options!"),
                h3("Add LD"),
                radioGroupButtons(inputId = "linkage_type", label = "LD type:", choices = c("No LD", "Input variant", "Most significant in window"), status = "danger", checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon"))),      # select LD type
                bsPopover("linkage_type", title='LD type', content='Define if LD should be shown (only available when 1 study is selected). If LD should be shown, please indicate in which populations (from 1000Genome Project) LD should be calculated. ', placement="right",trigger="hover", options = list(container = "body")),
                h5("Select population(s) to calculate LD"),                                   # panels for the choice of the population of interest
                prettyCheckboxGroup(inputId ="pop_eur", label = "European", choices = c("CEU" = "CEU", "FIN" = "FIN", "GBR" = "GBR", "IBS" = "IBS", "TSI" = "TSI", "ALL Sub-Pop" = "ALL_eur"), inline = TRUE, status = "info", fill = TRUE),
                bsPopover("pop_eur", title='European populations', content='CEU: Utah Residents with Northern and Western European Ancestry -- FIN: Finnish -- GBR: British in England and Scotland -- IBS: Iberian Population in Spain -- TSI: Toscani in Italy', placement="right",trigger="hover", options = list(container = "body")),
                prettyCheckboxGroup(inputId = "pop_afr", label = "African", choices = c("ACB" = "ACB", "ASW" = "ASW", "ESN" = "ESN", "GWD" = "GWD", "LWK" = "LWK", "MSL" = "MSL", "YRI" = "YRI", "ALL Sub-Pop" = "ALL_afr"), inline = TRUE, status = "warning", fill = TRUE),
                bsPopover("pop_afr", title='African populations', content='ACB: African Caribbeans in Barbados -- ASW: Americans of African Ancestry in SW USA -- ESN: Esan in Nigeria -- GWD: Gambian in Western Divisions in the Gambia -- LWK: Luhya in Webuye, Kenya -- MSL: Mende in Sierra Leone -- YRI: Yoruba in Ibadan, Nigeria', placement="right",trigger="hover", options = list(container = "body")),
                prettyCheckboxGroup(inputId = "pop_amr", label = "American", c("CLM" = "CLM", "MXL" = "MXL", "PEL" = "PEL", "PUR" = "PUR", "ALL Sub-Pop" = "ALL_amr"), inline = TRUE, status = "danger", fill = TRUE),
                bsPopover("pop_amr", title='American populations', content='CLM: Colombians from Medellin, Colombia -- MXL: Mexican Ancestry from Los Angeles USA -- PEL: Peruvians from Lima, Peru -- PUR: Puerto Ricans from Puerto Rico', placement="right",trigger="hover", options = list(container = "body")),
                prettyCheckboxGroup(inputId = "pop_eas", label = "East-Asian", choices = c("CDX" = "CDX", "CHB" = "CHB", "CHS" = "CHS", "JPT" = "JPT", "KHV" = "KHV", "ALL Sub-Pop" = "ALL_eas"), inline = TRUE, status = "success", fill = TRUE),
                bsPopover("pop_eas", title='East-Asian populations', content='CDX: Chinese Dai in Xishuangbanna, China -- CHB: Han Chinese in Beijing, China -- CHS: Southern Han Chinese -- JPT: Japanese in Tokyo, Japan -- KHV: Kinh in Ho Chi Minh City, Vietnam', placement="right",trigger="hover", options = list(container = "body")),
                prettyCheckboxGroup(inputId = "pop_sas", label = "South-Asian", choices = c("BEB" = "BEB", "GIH" = "GIH", "ITU" = "ITU", "PJL" = "PJL", "STU" = "STU", "ALL Sub-Pop" = "ALL_sas"), inline = TRUE, status = "primary", fill = TRUE),
                bsPopover("pop_sas", title='South-Asian populations', content='BEB: Bengali from Bangladesh -- GIH: Gujarati Indian from Houston, Texas -- ITU: Indian Telugu from the UK -- PJL: Punjabi from Lahore, Pakistan -- STU: Sri Lankan Tamil from the UK', placement="right",trigger="hover", options = list(container = "body")),
                h6("** If you don't select any population, by default European samples of the 1000Genome project will be used. Please refer to the Help section for additional information about the populations reported above."),
                div(style = "margin-top: 10px"),
                h6("** In case the top SNP is not in 1000Genome data, we take the second-most significant, then the third-most significant, etc."),
              ), align = "center"),
            ),
            column(6,
              div(tags$h3("Structural Variants"), align = "center"),
              div(dropdown(
                style = "unite", size = "lg", icon = icon("eraser"), status = "warning", width = "500px", animate = animateOptions(enter = animations$fading_entrances$fadeInLeftBig, exit = animations$fading_exits$fadeOutLeftBig),
                tooltip = tooltipOptions(title = "Click to see available options!"),
                h3("Choose dataset"),
                radioGroupButtons(inputId = "strVar_inp", label = "Structural variants data", choices = c("Linthrost et al. 2019" = "jasper", "Chaisson et al. 2019" = "chaisson", "Audano et al 2019" = "audano", "All" = "all"), status = "primary", checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon"))),
                bsPopover("strVar_inp", title='Structural variants', content='Third-generation sequencing (long-read sequencing) studies estimating the major Structural Variants (SV) alleles.', placement="right",trigger="hover", options = list(container = "body")),
              ), align = "center"),
            ),
            column(6,
              div(tags$h3("GWAS Catalog settings"), align = "center"),
              div(dropdown(
                style = "unite", size = "lg", icon = icon("table"), status = "danger", width = "500px", animate = animateOptions(enter = animations$fading_entrances$fadeInLeftBig, exit = animations$fading_exits$fadeOutLeftBig),
                tooltip = tooltipOptions(title = "Click to see available options!"),
                h3("Choose dataset"),
                radioGroupButtons(inputId = "gwascat_type", label = "GWAS Catalog data", choices = c("SNPs", "Genes"), status = "primary", checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon"))),
                bsPopover("gwascat_type", title='GWAS Catalog', content='Match results in GWAS Catalog by SNPs (default) or by Gene.', placement="right",trigger="hover", options = list(container = "body")),
              ), align = "center"),
            ),
            column(12,
              div(tags$h3("Download center"), align = "center"),
              div(dropdown(
                style = "unite", size = "lg", icon = icon("download"), status = "success", width = "300px", animate = animateOptions(enter = animations$fading_entrances$fadeInLeftBig, exit = animations$fading_exits$fadeOutLeftBig),
                tooltip = tooltipOptions(title = "Click to see available options!"),
                div(h5("SNPs table"), align = "center"),
                div(downloadButton("download_SNPsTable", "Download SNPs"), align = "center"),
                hr(),
                div(h5("GWAS Catalog results"), align = "center"),
                div(downloadButton("download_GWASCat", "Download All Results"), align = "center"),
                hr(),
                div(h5("SVs in the region"), align = "center"),
                div(downloadButton("download_SVs", "Download All SVs"), align = "center"),
                hr(),
                div(h5("LD info"), align = "center"),
                div(downloadButton("download_LDtable", "Download LD table"), align = "center"),
                hr(),
                div(h5("Plot"), align = "center"),
                div(h6("You can get a .png image of the plot by hovering on the plot. Click on the camera icon!"), align = "center"),
                div(h6("Alternatively, here you can get a .pdf plot or the original .html page (preserves hovering information)."), align = "center"),
                radioGroupButtons(inputId = "plot_format", label = "", choices = c("PDF", "HTML"), status = "warning", checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon"))),
                div(downloadButton("download_plot", "Download Plot"), align = "center"),
                div(h6('Note that it will take about 10 seconds to produce your plot!'), align = "center"),
              ), align = "center"),
            ),
            div(h6("ss"), style = "color: #f2f2f2", align = "center"),
            hr(),
            div(h3("Top SNPs info"), align = "center"),
            div(tableOutput('table'), align = "center", style = 'font-size:90%; overflow-x: scroll; overflow-y: scroll; height:250px'),
            bsPopover("table", title='Top SNPs table', content='Most significant SNP-associations in the main plot interface. Chr: chromosome; Position: genomic position depending on reference genome choosen; ID: variant identifier (missing if variant is not in 1000Genome or too rare); -log10(P): association p-value (in -log10 scale); Study: name of the study; Alleles: reference and alternative alleles, respectively (from 1000Genome data).', placement="right",trigger="hover", options = list(container = "body")),
            hr(),
            div(h3("GWAS Catalog SNPs"), align = "center"),
            div(tableOutput('gwascat_table'), align = "center", style = 'font-size:90%; overflow-x: scroll; overflow-y: scroll; height:250px'),
            bsPopover("gwascat_table", title='GWAS Catalog table', content='Most significant associations of SNPs and genes in the plotted region as found in the GWAS Catalog (updated on 2021-12-21).', placement="right",trigger="hover", options = list(container = "body")),
            hr(),
            div(h3("Structural Variants"), align = "center"),
            div(tableOutput('sv_table'), align = "center", style = 'font-size:90%; overflow-x: scroll; overflow-y: scroll; height:250px'),
            bsPopover("sv_table", title='Structual variants table', content='Most relevant structural variations in the region of interest and from the chosen source.', placement="right",trigger="hover", options = list(container = "body")),
          )
        ),
        column(9,
          plotlyOutput(outputId = "plot", height = "1300px"),
        )
      ),
      hr(),
      fluidRow(
        column(3,
          wellPanel(
            div(tags$h2("eQTL and sQTL"), align = "center"),            
          )
        ),
        column(9,
          plotlyOutput(outputId = "plot2", height = "500px"),
        )
      )
    )
  )
