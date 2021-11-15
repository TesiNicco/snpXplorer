#######################################################################################
# For the annotation pipeline of snpXplorer, it is probably best to divide the scripts
# This will in theory help a lot with the memory and resources consumption
# And hopefully there will be no failed jobs
#######################################################################################

###########################################################################################################
# This script is the main script that leads the job, runs external scripts, send emails and reports errors
###########################################################################################################

# LIBRARIES
    suppressPackageStartupMessages({
        library(data.table)
        library(stringr)
        library(GenomicRanges)
        library(liftOver)
    })

# FUNCTIONS
    ## function to wait some time -- to manage memory consumption
    waitForMe <- function(x){
        p1 <- proc.time()
        Sys.sleep(x)
        proc.time() - p1 # The cpu usage should be negligible
    }



# READ ARGUMENTS AND MAIN PATH
    ## SET MAIN PATH
    MAIN = "/root/snpXplorer/AnnotateMe/"
    args = commandArgs(trailingOnly=TRUE)

    # ARGUMENTS
    fname <- args[1]
    #fname <- "annotateMe_input_86734.txt"
    ftype <- args[2]
    #ftype <- 3
    username <- args[3]
    # analysis type -- gene-set enrichment or mapping only
    analysis_type = args[4]
    # gene-sets for enrichment analysis
    analysis_mode = unlist(strsplit(as.character(args[5]), ","))
    # tissues of interest
    interesting_tissues = unlist(strsplit(as.character(args[6]), ","))
    # reference genome (in case input is not rsid)
    ref_version = as.character(args[7])

    # CREATE FOLDER FOR RESULTS -- ADD RANDOM NUMBER AND COPY INPUT FILE IN THERE
    random_num <- sample(x = seq(1, 100000), size = 1, replace = F)
    system(paste("mkdir /root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, sep=""))
    system(paste("mv /root/snpXplorer/snpXplorer_v3/", fname, " /root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/", sep=""))

    # SEND EMAIL TO MYSELF AS A NOTIFICATION SOMEONE REQUESTED A JOB
    cmd_mail <- paste("sendEmail -f snpxplorer@gmail.com -t n.tesi@amsterdamumc.nl -u 'AnnotateMe request sent' -cc snpxplorer@gmail.com -m 'Hello, \n a request to AnnotateMe was just sent from ", username, ". \n The following settings were requested: \n input --> ", fname, "\n input_type --> ", ftype, "\n analysis_mode --> ", analysis_mode, "\n interest_tissue --> ", interesting_tissues, "\n ref_version --> ", ref_version, "\n output_folder --> ", random_num, "\n \n AnnotateMe' -S /usr/sbin/sendmail")
    system(cmd_mail)

    ## START OF THE PIPELINE
    ## Before the start, we need to check whether enough memory is available otherwise it will crash at some point and i have to look at it manually
    ## So we check the available memory, and if that is >10Gb (50% of all available memory) then the process starts
    # first initialize the START variable
    START = FALSE
    while (START == FALSE){
        # calculate available memory
        memfree <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern=T))
        # check if this is >8Gb
        if (memfree > 12000000){
            START = TRUE
            print("Start the job!")
        } else {
            # if there's not enough memory, wait for a minute and try again
            waitForMe(60)
            print("Wait to start the job!")
            START = FALSE
        }
    }

    # IN CASE ENOUGH MEMORY IS AVAILABLE, LET'S START
    if (START == TRUE){
        # FIRST STEP IS TO READ SNPS OF INTEREST AND REPORT ERRORS IN CASE OF WRONG INPUT AND/OR TOO LONG INPUT
        cat("## Reading SNPs and positions\n")
        inpf = paste0("/root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/", fname)
        snps_info_path = system(paste0("Rscript ", MAIN, "BIN_v2/readSNPs.R ", inpf, " ", ftype, " ", ref_version, " ", analysis_type, " ", random_num), intern = T)
        load(snps_info_path)

      # CHECK WHETHER INPUT LIST OF SNPS WAS CORRECT
    if (length(data) == 1){
        system("sendEmail -f snpxplorer@gmail.com -t ", username, " -u 'snpXplorer input error' -m 'Dear user, \n thanks so much for using snpXplorer and its annotation pipeline. \n Unfortunately, an error occurred while reading the input SNPs you provided. Possible reasons include a number of SNPs larger than 1000 or a wrong input type. \n Please correct the input and try again. In the More/Help section of the website you can find example datasets. \n Please do not hesitate to contact us in case of any question. \n snpXplorer team.' -a '", fname, "' -cc n.tesi@amsterdamumc.nl snpxplorer@gmail.com -S /usr/sbin/sendmail")
        # don't remove data, but zip them
        system(paste("tar -czf AnnotateMe_results_", random_num, ".tar.gz RESULTS_", random_num, "/", sep=""))
    } else {
        # FIND ALL VARIANTS IN LD WITH THE INPUT VARIANTS -- this can be done with LDlink functions easily -- question is what to do with that?
        ld_info = data.frame(chr=NULL, pos=NULL)

        # RUN CADD
        cadd_info_path = system(paste0("Rscript ", MAIN, "BIN_v2/CADD.R ", random_num, " ", ref_genome, " ", ld_info, " ", snps_info_path), intern = T)
