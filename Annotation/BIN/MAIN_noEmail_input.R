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
        library(LDlinkR)
    })

# FUNCTIONS
    ## function to wait some time -- to manage memory consumption
    waitForMe <- function(x){
        p1 <- proc.time()
        Sys.sleep(x)
        proc.time() - p1 # The cpu usage should be negligible
    }

    ## function to plot mapping charcteristics -- adjusted for faster computations (library-wise)
    plotMapping <- function(mapping.after, MAIN){
        #set how many plots to make and the margins
        layout(matrix(c(1, 2, 3, 3, 4, 4, 4, 4), nrow = 4, ncol = 2, byrow = TRUE))
        par(mar=c(3, 3, 3, 3))

        #color palette
        colz <- ggsci::pal_lancet(palette = "lanonc")(5)
        myPalette <- RColorBrewer::brewer.pal(3, "Set2")

        # make sure frequency is numeric
        mapping.after$ALT_FREQS = as.numeric(mapping.after$ALT_FREQS)

        #plot 1 -- barplot of the annotation sources
        labs <- c()
        if (nrow(mapping.after[which(mapping.after$source_finalGenes == "coding"),]) >0){ lab.cod <- paste("Coding\n(N=", nrow(mapping.after[which(mapping.after$source_finalGenes == "coding"),]), ")", sep=""); labs <- c(labs, lab.cod) }
        if (nrow(mapping.after[which(mapping.after$source_finalGenes == "eqtl+cadd"),]) >0){ lab.eqtl <- paste("eQTL\n(N=", nrow(mapping.after[which(mapping.after$source_finalGenes == "eqtl+cadd"),]), ")", sep=""); labs <- c(labs, lab.eqtl) }
        if (nrow(mapping.after[which(mapping.after$source_finalGenes == "positional"),]) >0){ lab.pos <- paste("Position\n(N=", nrow(mapping.after[which(mapping.after$source_finalGenes == "positional"),]), ")", sep=""); labs <- c(labs, lab.pos) }
        pie(table(mapping.after$source_finalGenes), col=myPalette, labels=labs, lwd=1.5, cex=1.50)

        #plot 2 -- histogram of snp-gene number
        red.na <- mapping.after[which(is.na(mapping.after$geneList)), c("locus", "geneList")]
        red <- mapping.after[which(!is.na(mapping.after$geneList)), c("locus", "geneList")]
        #output data frame
        df <- as.data.frame(matrix(data=NA, nrow=nrow(mapping.after), ncol=2))
        colnames(df) <- c("locus", "n_genes")
        for (i in 1:nrow(red)){
            #split geneList column
            n.gens <- strsplit(as.character(red$geneList[i]), ",")
            df[i, ] <- c(as.character(red$locus[i]), length(n.gens[[1]]))
        }
        if (nrow(red.na) >0){ for (j in 1:nrow(red.na)){ df[(i+j),] <- c(as.character(red.na$locus[j]), 0) } }
        #get max
        df$n_genes <- as.numeric(df$n_genes)
        mx <- max(df$n_genes)
        #then the plot
        par(mar=c(5, 8, 5, 10))
        barplot(table(df$n_genes), xlab="Genes per variant", ylab="Frequency", main="", cex.axis=1.25, col=colz, cex.lab=1.5, cex.names=1.25)

        #plot 3 --  histogram of distribution per chromosome
        par(mar=c(4, 8, 3, 4))
        colorz <- grDevices::colorRampPalette(c("red", "orange", "dark green", "blue"))
        colorz.chr <- colorz(22)
        tmp <- as.data.frame(table(mapping.after$chr))
        tmp$Var1 <- as.numeric(as.character(tmp$Var1))
        if (nrow(tmp) != 22){ for (i in 1:22){ if (!(i %in% tmp$Var1)){ t <- data.frame(Var1=i, Freq=0); tmp <- rbind(tmp, t) } } }
        tmp <- tmp[order(tmp$Var1),]
        barplot(tmp$Freq, col=colorz.chr, names=tmp$Var1, xlab="Chromosome", ylab="Genes per chromosome", cex.lab=1.50, main="", cex.axis=1.4, cex.main=2)

        #plot 4 -- circular plot integrated visualization
        # read chromosome length
        l <- read.table(paste(MAIN, "INPUTS_OTHER/chromosomes_length_hg19.txt", sep=""), h=F)
        tmp.df <- data.frame(V1="AX", V2=50000000)
        l <- rbind(tmp.df, l)

        # circle circumference is the sum of all chromosome lengths
        R <- sum(l$V2)/1000000

        # then derive radius
        r <- R/(2*pi)

        # find corresponding angles for the chromosomes
        ang <- (l$V2/1000000)/r

        # set colors
        colz <- grDevices::colorRampPalette(c("yellow", "light green", "dark green", "light blue", "navy", "orange", "red"))
        colorz <- colz(22)
        col.bk <- rep(c("white", "grey90"), 11)

        # set pch
        mapping.after$pch <- 21
        mapping.after$pch[which(mapping.after$source_finalGenes == "coding")] <- 23
        mapping.after$pch[which(mapping.after$source_finalGenes == "eqtl")] <- 24

        par(mar=c(4, 18, 4, 18))
        plot(0, 0, pch=16, col="white", xlim=c(-r, r), ylim=c(-r, r), bty='n', xlab="", ylab="", xaxt='none', yaxt="none")
        plotrix::draw.circle(x = 0, y = 0, radius = r, border=NA)
        cum.angle <- pi/2-(ang[1]/2)
        # first loop to add sectors
        for (i in 1:23){
            # get degrees
            dg1 <- cum.angle*180/pi
            dg2 <- (cum.angle+ang[i])*180/pi
            # draw background and borders
            circlize::draw.sector(start.degree = dg1, end.degree = dg2, rou1 = r, rou2 = 0, clock.wise = FALSE, col=ggplot2::alpha(col.bk[i], 0.6), border = NA)
            cum.angle <- cum.angle + ang[i]
        }

        # add axis for maf
        a = r*0.85
        b = r*0.30
        d = r*0.575
        c = r*0.01
        segments(x0 = 0, y0 = b, x1 = 0, y1 = a, lwd=1.25, xpd=T)
        # add maf circles for reference
        plotrix::draw.circle(x = 0, y = 0, radius = a-a*0.025, border="grey40", lty=2)
        plotrix::draw.circle(x = 0, y = 0, radius = b-b*0.025, border="grey40", lty=2, nv = 1000)
        plotrix::draw.circle(x = 0, y = 0, radius = d-d*0.025, border="grey40", lty=2, nv = 1000)
        # add axis ticks
        for (i in seq(0, 0.5, 0.1)){
            trx <- (b-a)*(i/0.5)+a
            segments(x0 = -c, y0 = trx, x1 = c, y1 = trx, lwd=1.25, xpd=T)
            # add text
            text(x = 0, y = trx, labels = i, pos=2, offset = 0.2, cex=0.75, font=2)
        }
        dg1 <- (pi/2-(ang[1]/2))*180/pi
        dg2 <- ((pi/2+(ang[1]/2)))*180/pi
        circlize::draw.sector(start.degree = dg1, end.degree = dg2, rou1 = r, rou2 = r-r*0.1, clock.wise = FALSE, col = NA)
        text(x = 0, y = r*0.95, labels = "MAF", font=2, adj=0.5, cex = 0.70)

        # second loop to add points
        cum.angle <- pi/2+(ang[1]/2)
        for (i in 1:22){
            # calculate points for the segment of the lolliplot -- point0 should be maf
            pp <- mapping.after[which(mapping.after$chr == i),]
            pp$ALT_FREQS[which(pp$ALT_FREQS >0.5)] <- 1-pp$ALT_FREQS[which(pp$ALT_FREQS >0.5)]

            # for the frequency, need to normalize it in [r*0.20-r*0.80] [b, a]
            pp$maf_norm <- (b-a)*(pp$ALT_FREQS/0.5)+a

            # get the angle
            pp$ang <- (as.numeric(pp$pos)/1000000)/r + cum.angle
            pp$x <- pp$maf_norm*cos(pp$ang)
            pp$y <- pp$maf_norm*sin(pp$ang)
            segments(x0 = pp$x, y0 = pp$y, x1 = r*cos(pp$ang), y1 = r*sin(pp$ang), lwd=0.85, col="grey40")
            points(x = pp$x, y = pp$y, type = "p", pch=pp$pch, col="black", bg=colorz[i], cex=1.75)

            # draw sector finally
            # get degrees
            dg1 <- cum.angle*180/pi
            dg2 <- (cum.angle+ang[i+1])*180/pi
            circlize::draw.sector(start.degree = dg1, end.degree = dg2, rou1 = r, rou2 = r-r*0.1, clock.wise = FALSE, col=colorz[i])

            # add label
            rgm <- cum.angle+(ang[i+1]/2)
            text(x = r*cos(rgm)*0.95, y = r*sin(rgm)*0.95, labels = i, font=2, adj=0.5)

            # increment at the end
            cum.angle <- cum.angle + ang[i+1]
        }
        # finally the legend
        legend(x = -r/2, y = -r*1.05, xpd=T, legend = c("Coding", "eQTL", "Position"), pch=c(23, 24, 21), pt.bg="grey40", pt.cex = 1.50, cex=1.50, bty='n', ncol=3)

        return(df)
        }

##############################
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
    analysis_mode_all = as.character(args[5])
    # tissues of interest
    interesting_tissues = unlist(strsplit(as.character(args[6]), ","))
    interesting_tissues_all = as.character(args[6])
    # reference genome (in case input is not rsid)
    ref_version = as.character(args[7])
    # random number should also be inputed
    random_num = as.numeric(args[8])

    # CREATE FOLDER FOR RESULTS -- ADD RANDOM NUMBER AND COPY INPUT FILE IN THERE
    #random_num <- sample(x = seq(1, 100000), size = 1, replace = F)        # no need to run this
    #system(paste("mkdir /root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, sep=""))         # no need to run this -- folder should be already there
    #system(paste("mv /root/snpXplorer/snpXplorer_v3/", fname, " /root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/", sep=""))         # no need to run this -- file should be already there

    # SEND EMAIL TO MYSELF AS A NOTIFICATION SOMEONE REQUESTED A JOB
    cmd_mail <- paste("sendEmail -f n.tesi@amsterdamumc.nl -t ", username, " -cc n.tesi@amsterdamumc.nl snpxplorer@gmail.com -u 'AnnotateMe request sent' -m 'Dear user, \n snpXplorer received an annotation request from you. \n You receive this email to confirm that your request is under processing. A typical job takes about 30 minutes to complete, however, due to the high number of requests, jobs may be delayed. \n\nThe following settings were requested: \n input --> ", fname, "\n input_type --> ", ftype, "\n analysis_type --> ", analysis_type, "\n analysis_mode --> ", analysis_mode_all, "\n interest_tissue --> ", interesting_tissues_all, "\n ref_version --> ", ref_version, "\n output_folder --> ", random_num, "\n \n snpXplorer Team' -S /usr/sbin/sendmail")
    #system(cmd_mail)           # email should not be sent

    ## START OF THE PIPELINE
    ## Before the start, we need to check whether enough memory is available otherwise it will crash at some point and i have to look at it manually
    ## So we check the available memory, and if that is >10Gb (50% of all available memory) then the process starts
    # first initialize the START variable
    START = FALSE
    while (START == FALSE){
        # calculate available memory
        meminfo = system("free -h | sed 's/  */ /g'", intern = T); values = list()
        used = str_split_fixed(meminfo[2], ' ', 7)[,3]; free = str_split_fixed(meminfo[2], ' ', 7)[,4]; cache = str_split_fixed(meminfo[2], ' ', 7)[,6]
        for (v in c(used, free, cache)){
            if (length(grep('G', v)) > 0){
                v = str_replace_all(v, 'G', ''); v = str_replace_all(v, ',', '\\.'); v = as.numeric(v)
            } else if (length(grep('M', v)) > 0){
                v = str_replace_all(v, 'M', ''); v = str_replace_all(v, ',', '\\.'); v = paste0('0.', v); v = as.numeric(v)
            }
            values[[(length(values) + 1)]] = v
        }
        used = values[[1]]; free = values[[2]]; cache = values[[3]]
        #memfree <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern=T))
        memfree = tryCatch({ free + cache }, error = function(e){ 0 }) 
        # check if this is >8Gb
        if (memfree > 8 && analysis_type != "enrichment"){
            START = TRUE
            print("Start the job!")
        } else if (memfree > 8 && analysis_type == "enrichment"){
            START = TRUE
            print("Start the job!")
        } else {
            # if there's not enough memory, wait for a minute and try again
            print("Wait to start the job!")
            waitForMe(300)
            START = FALSE
        }
    }
    # IN CASE ENOUGH MEMORY IS AVAILABLE, LET'S START
    if (START == TRUE){
        # FIRST STEP IS TO READ SNPS OF INTEREST AND REPORT ERRORS IN CASE OF WRONG INPUT AND/OR TOO LONG INPUT
        cat("## Reading SNPs and positions\n")
        inpf = paste0("/root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/", fname)
        # before starting, make a quick parse of the input file as people sometimes use wrong input like comma
        cmd = paste0("sed -i 's/,/\\n/g' ", inpf, " | sed '/^$/d'"); system(cmd)
        cmd = paste0("sed -i 's/chr//g' ", inpf); system(cmd)
        snps_info_path = system(paste0("Rscript ", MAIN, "BIN/readSNPs.R ", inpf, " ", ftype, " ", ref_version, " ", analysis_type, " ", random_num), intern = T)
        load(snps_info_path)

        # CHECK WHETHER INPUT LIST OF SNPS WAS CORRECT
        if (!is.na(data)){ if (length(unique(data$chr)) == 1 && unique(data$chr) == "NA"){ data = NA } }
        if ((is.na(data)) || (length(data) == 1) || (nrow(data) == 0)){
            system(paste0("sendEmail -f n.tesi@amsterdamumc.nl -t ", username, " -u 'snpXplorer input error' -m 'Dear user, \n\n thanks so much for using snpXplorer and its annotation pipeline. \n\n Unfortunately, an error occurred while reading the input SNPs you provided. Possible reasons include: \n- the number of SNP(s) is >1000 for enrichment analysis; \n- the number of SNP(s) is >10000 for mapping analysis; \n- the input type is wrong; \n\n Please correct the input and try again. In the More/Help section of the website you can find example datasets. \n\n Please do not hesitate to contact us in case of any question. \n snpXplorer team.' -a '", inpf, "' -cc n.tesi@amsterdamumc.nl snpxplorer@gmail.com -S /usr/sbin/sendmail"))
            # zip data
            system(paste("tar -czf AnnotateMe_results_", random_num, ".tar.gz RESULTS_", random_num, "/", sep=""))
            # remove data
            system(paste("rm -rf RESULTS_", random_num, "/", sep=""))
        } else {
            # SAVE DATA WITH MISSING ANNOTATION AS THEY WILL BE SKIPPED
            missing_data = data[is.na(data$pos),]
            data = data[!is.na(data$pos),]
            save(data, file = snps_info_path)

            # FIND ALL VARIANTS IN LD WITH THE INPUT VARIANTS -- this can be done with LDlink functions easily -- for now just report SNPs in LD, in the future those results will be implemented in the main results
            ld_info = data.frame(chr=NULL, pos=NULL)
            ld_outpath = paste0("/root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/tmp_ldInfo.RData")
            save(ld_info, file = ld_outpath)
            #LDproxy_batch(data$ID[grep("rs", data$ID)], pop = "ALL", r2d = "r2", token = "b2735a858324", append = TRUE)

            # RUN CADD
            cat("## Running CADD annotation\n")
            system(paste0("Rscript ", MAIN, "BIN/CADD.R ", random_num, " ", ref_version, " ", ld_outpath, " ", snps_info_path), ignore.stdout=T, ignore.stderr=F)
            #load(snps_info_path)

            # GTEX ANNOTATION -- eQTLs
            cat("## Running GTEx annotation ~ eQTLs\n")
            #gtex_info = LDexpress(out.annot$ID, pop = "ALL", tissue = interesting_tissues, r2d = "r2", p_threshold = 0.05, token = "b2735a858324")
            system(paste0("Rscript ", MAIN, "BIN/GTEx.R ", snps_info_path, " ", interesting_tissues_all, " ", random_num), ignore.stdout=T, ignore.stderr=F)
            #load(snps_info_path)

            # GTEX ANNOTATION -- sQTLs
            cat("## Running GTEx annotation ~ sQTLs\n")
            system(paste0("Rscript ", MAIN, "BIN/sQTL.R ", snps_info_path, " ", interesting_tissues_all, " ", random_num), ignore.stdout=T, ignore.stderr=F)
            #load(snps_info_path)

            # POSITIONAL MAPPING
            cat("## Running positional annotation\n")
            system(paste0("Rscript ", MAIN, "BIN/PosMapping.R ", snps_info_path), ignore.stdout=T, ignore.stderr=F)
            load(snps_info_path)
            annot <- final_res[[1]]
            geneList <- final_res[[2]]

            # LOOK-UP GWAS-catalog
            cat("## Looking the GWAS catalog\n")
            #all_traits = LDtrait(annot$ID[grep("rs", annot$ID)], pop = "ALL", r2d = "r2", r2d_threshold = 0.4, token = "b2735a858324")
            #all_traits_snps_only = all_traits[which(all_traits$RS_Number %in% annot$ID & all_traits$P_value <= 5e-8),]
            #write.table(all_traits_snps_only, paste("/root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/lookup_gwas_catalog.txt", sep=""), quote=F, row.names=F, col.names=F)
            system(paste0("Rscript ", MAIN, "BIN/GWAScatalog.R ", snps_info_path, " ", random_num), ignore.stdout=T, ignore.stderr=F)

            # STRUCTUAL VARIATIONS IN THE VICINITY
            cat("## Running SV analysis\n")
            system(paste0("Rscript ", MAIN, "BIN/StrVar.R ", snps_info_path, " ", random_num), ignore.stdout=T, ignore.stderr=F)

            # DONE -- SAVE TABLES
            cat("## Saving tables, plotting and finishing\n")
            # add also missing data when saving results so that user has all data that has inputed
            annot_with_miss = plyr::rbind.fill(annot, missing_data)
            write.table(annot_with_miss, paste("/root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/snp_annotation.txt", sep=""), quote=F, row.names=F, sep="\t")
            write.table(geneList, paste("/root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/snp_annotation_geneList.txt", sep=""), quote=F, row.names=F, col.names=F)

            # PLOT MAPPING CHARACTERISTICS
            pdf(paste("/root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/snp_gene_mapping.pdf", sep=""), height=12.35, width=10)
            genesPerSNP <- plotMapping(annot, MAIN)
            invisible(dev.off())

            # HERE DECIDE WHETHER TO DO GENE-SET ENRICHMENT ANALYSIS OR NOT DEPENDING ON ANALYSIS_TYPE
            if (analysis_type == "enrichment"){
                cat("## Gene-set enrichment analysis was requested. Working on it now\n")

                # DO ENRICHMENT ANALYSIS ONLY IF THERE ARE ENOUGH GENES
                if (length(unique(geneList)) >1){
                    # RUN ENRICHMENT ANALYSIS
                    system(paste0("Rscript ", MAIN, "BIN/GeneSetEnrichment.R ", snps_info_path, " ", random_num, " ", analysis_mode_all))
                    load(paste0("RESULTS_", random_num, "/tmp_enrichRes.RData"))
                    dbs = sampling.res[[3]]

                    # NOW CHECK IF WE NEED TO DO THE GO ANALYSIS
                    if (length(grep("GO", dbs)) >0){
                        cat("## REVIGO analysis\n")
                        # DO REVIGO ANALYSIS
                        system(paste0("Rscript ", MAIN, "BIN/GO_analysis_REVIGO.R ", random_num))

                        # DO ALTERNATIVE TO REVIGO
                        cat("## Semantic Similarity analysis\n")
                        system(paste0("Rscript ", MAIN, "BIN/GO_analysis_alternative.R ", random_num), ignore.stdout=T, ignore.stderr=F)

                        # CLEAN TEMPORARY DATA
                        system(paste0("rm /root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/tmp_*"))
                        #system(paste0("mv ", MAIN, "BIN/combined_query_snp_list.txt /root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/"))
                        #system(paste0("mv /root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/combined_query_snp_list.txt /root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/LD_information.txt"))
                        system(paste0("rm /root/snpXplorer/snpXplorer_v3/tmp_", random_num, ".json"))

                        # FINISH -- SEND DATA BACK TO THE OWNER
                        system(paste0("cp /root/snpXplorer/snpXplorer_v3/www/snpXplorer_output_description.pdf RESULTS_", random_num, "/"))
                        system(paste("tar -czf AnnotateMe_results_", random_num, ".tar.gz RESULTS_", random_num, "/", sep=""))
                        link = system(paste0("curl --upload-file AnnotateMe_results_", random_num, ".tar.gz https://transfer.sh/AnnotateMe_results_", random_num, ".tar.gz"), intern = T)
                        # finally the email
                        system(paste0("sendEmail -f snpxplorer@gmail.com -t ", username, " -u 'snpXplorer results' -m 'Dear user, \n\n thanks so much for using snpXplorer and its annotation pipeline. \n We hope you find the tool useful. \n\n We now provide a WeTransfer link to download your results. Please find the link below: \n", link, " \n\n Best wishes, \n snpXplorer team.' -cc n.tesi@amsterdamumc.nl snpxplorer@gmail.com -S /usr/sbin/sendmail"))
                        system(paste("rm -rf RESULTS_", random_num, "/", sep=""))
                    } else {
                        # CLEAN TEMPORARY DATA
                        system(paste0("rm /root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/tmp_*"))
                        #system(paste0("mv ", MAIN, "BIN/combined_query_snp_list.txt /root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/"))
                        #system(paste0("mv /root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/combined_query_snp_list.txt /root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/LD_information.txt"))
                        system(paste0("rm /root/snpXplorer/snpXplorer_v3/tmp_", random_num, ".json"))

                        # FINISH -- SEND DATA BACK TO THE OWNER
                        system(paste0("cp /root/snpXplorer/snpXplorer_v3/www/snpXplorer_output_description.pdf RESULTS_", random_num, "/"))
                        system(paste("tar -czf AnnotateMe_results_", random_num, ".tar.gz RESULTS_", random_num, "/", sep=""))
                        link = system(paste0("curl --upload-file AnnotateMe_results_", random_num, ".tar.gz https://transfer.sh/AnnotateMe_results_", random_num, ".tar.gz"), intern = T)
                        # finally the email
                        system(paste0("sendEmail -f snpxplorer@gmail.com -t ", username, " -u 'snpXplorer results' -m 'Dear user, \n\n thanks so much for using snpXplorer and its annotation pipeline. \n We hope you find the tool useful. \n\n We now provide a WeTransfer link to download your results. Please find the link below: \n", link, " \n\n Best wishes, \n snpXplorer team.' -cc n.tesi@amsterdamumc.nl snpxplorer@gmail.com -S /usr/sbin/sendmail"))
                        system(paste("rm -rf RESULTS_", random_num, "/", sep=""))
                    }
                } else {
                    # CLEAN TEMPORARY DATA
                    system(paste0("rm /root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/tmp_*"))
                    #system(paste0("mv ", MAIN, "BIN/combined_query_snp_list.txt /root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/"))
                    #system(paste0("mv /root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/combined_query_snp_list.txt /root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/LD_information.txt"))
                    system(paste0("rm /root/snpXplorer/snpXplorer_v3/tmp_", random_num, ".json"))

                    # FINISH -- SEND DATA BACK TO THE OWNER
                    system(paste0("cp /root/snpXplorer/snpXplorer_v3/www/snpXplorer_output_description.pdf RESULTS_", random_num, "/"))
                    system(paste("tar -czf AnnotateMe_results_", random_num, ".tar.gz RESULTS_", random_num, "/", sep=""))
                    link = system(paste0("curl --upload-file AnnotateMe_results_", random_num, ".tar.gz https://transfer.sh/AnnotateMe_results_", random_num, ".tar.gz"), intern = T)
                    # finally the email
                    system(paste0("sendEmail -f snpxplorer@gmail.com -t ", username, " -u 'snpXplorer results' -m 'Dear user, \n\n thanks so much for using snpXplorer and its annotation pipeline. \n Unfortunately, due to a low number of genes found in your SNPs, we could not perform gene-set enrichment analysis. \n\n We now provide a WeTransfer link to download your results. Please find the link below: \n", link, " \n\n Best wishes, \n snpXplorer team.' -cc n.tesi@amsterdamumc.nl snpxplorer@gmail.com -S /usr/sbin/sendmail"))
                    system(paste("rm -rf RESULTS_", random_num, "/", sep=""))
                }
            } else {
                cat("## Gene-set enrichment was NOT requested. Finishing analysis now.\n")
                # CLEAN TEMPORARY DATA
                system(paste0("rm /root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/tmp_*"))
                #system(paste0("mv ", MAIN, "BIN/combined_query_snp_list.txt /root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/"))
                #system(paste0("mv /root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/combined_query_snp_list.txt /root/snpXplorer/snpXplorer_v3/RESULTS_", random_num, "/LD_information.txt"))
                #system(paste0("rm /root/snpXplorer/snpXplorer_v3/tmp_", random_num, ".json"))

                # FINISH -- SEND DATA BACK TO THE OWNER
                system(paste0("cp /root/snpXplorer/snpXplorer_v3/www/snpXplorer_output_description.pdf RESULTS_", random_num, "/"))
                system(paste("tar -czf AnnotateMe_results_", random_num, ".tar.gz RESULTS_", random_num, "/", sep=""))
                # need to upload the file on wetransfer
                link = system(paste0("curl --upload-file AnnotateMe_results_", random_num, ".tar.gz https://transfer.sh/AnnotateMe_results_", random_num, ".tar.gz"), intern = T)
                # then make the email
                system(paste0("sendEmail -f n.tesi@amsterdamumc.nl -t ", username, " -u 'snpXplorer results' -m 'Dear user, \n\n thanks so much for using snpXplorer and its annotation pipeline. \n We hope you find the tool useful. \n\n We now provide a WeTransfer link to download your results. Please find the link below: \n", link, " \n\n Best wishes, \n snpXplorer team.' -cc n.tesi@amsterdamumc.nl snpxplorer@gmail.com -S /usr/sbin/sendmail"))
                system(paste("rm -rf RESULTS_", random_num, "/", sep=""))
            }
        }
    }
