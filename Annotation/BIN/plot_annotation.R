# Script for plots in R

# Load necessary libraries without startup messages
suppressPackageStartupMessages({library(data.table)})
suppressPackageStartupMessages({library(argparse)})
suppressPackageStartupMessages({library(RColorBrewer)})
suppressPackageStartupMessages({library(circlize)})
suppressPackageStartupMessages({library(viridis)})
suppressPackageStartupMessages({library(plotrix)})

# Define command line arguments
parser <- ArgumentParser()
parser$add_argument("-i", required=TRUE, help="Variant annotation as performed with standalone_annotation.py")
parser$add_argument("-g", required=FALSE, default=FALSE, help="GWAS annotation as performed with standalone_annotation.py")
parser$add_argument("-o", required=TRUE, help="Output folder for plots")
args <- parser$parse_args()

# Functions
# Function to plot the distribution of variant consequences
plotMapping <- function(variant_annotation){
    # Set how many plots to make and the margins
    layout(matrix(c(1, 2, 3, 3, 4, 4, 4, 4), nrow = 4, ncol = 2, byrow = TRUE))
    par(mar=c(3, 3, 3, 3))

    # Color palette
    myPalette <- RColorBrewer::brewer.pal(5, "Set2")

    # Make sure frequency is numeric
    variant_annotation$"Alternative Allele Frequency" = as.numeric(variant_annotation$"Alternative Allele Frequency")

    # Plot 1 -- pie of the annotation sources
    labs <- c()
    if (nrow(variant_annotation[which(variant_annotation$Source == "coding"),]) >0){ lab.cod <- paste("Coding\n(N=", nrow(variant_annotation[which(variant_annotation$Source == "coding"),]), ")", sep=""); labs <- c(labs, lab.cod) }
    if (nrow(variant_annotation[which(variant_annotation$Source == "qtl_query"),]) >0){ lab.eqtl <- paste("QTL\n(N=", nrow(variant_annotation[which(variant_annotation$Source == "qtl_query"),]), ")", sep=""); labs <- c(labs, lab.eqtl) }
    if (nrow(variant_annotation[which(variant_annotation$Source == "coding_ld"),]) >0){ lab.pos <- paste("Coding Haplotype\n(N=", nrow(variant_annotation[which(variant_annotation$Source == "coding_ld"),]), ")", sep=""); labs <- c(labs, lab.pos) }
    if (nrow(variant_annotation[which(variant_annotation$Source == "qtl_ld"),]) >0){ lab.pos <- paste("QTL Haplotype\n(N=", nrow(variant_annotation[which(variant_annotation$Source == "qtl_ld"),]), ")", sep=""); labs <- c(labs, lab.pos) }
    if (nrow(variant_annotation[which(variant_annotation$Source == "closest_gene"),]) >0){ lab.pos <- paste("Closest Gene\n(N=", nrow(variant_annotation[which(variant_annotation$Source == "closest_gene"),]), ")", sep=""); labs <- c(labs, lab.pos) }
    pie(table(variant_annotation$Source), col=myPalette, labels=labs, lwd=1.5, cex=1.50)

    # Plot 2 -- histogram of snp-gene number
    red.na <- variant_annotation[which(is.na(variant_annotation$"Most Likely Genes")), c("Locus ID", "Most Likely Genes")]
    red <- variant_annotation[which(!is.na(variant_annotation$"Most Likely Genes")), c("Locus ID", "Most Likely Genes")]
    # Output data frame
    df <- as.data.frame(matrix(data=NA, nrow=nrow(variant_annotation), ncol=2))
    colnames(df) <- c("locus", "n_genes")
    for (i in 1:nrow(red)){
        n.gens <- strsplit(as.character(red$"Most Likely Genes"[i]), ",")
        df[i, ] <- c(as.character(red$"Locus ID"[i]), length(n.gens[[1]]))
    }
    if (nrow(red.na) >0){ for (j in 1:nrow(red.na)){ df[(i+j),] <- c(as.character(red.na$"Locus ID"[j]), 0) } }
    # Get max
    df$n_genes <- as.numeric(df$n_genes)
    mx <- max(df$n_genes)
    # Then the plot
    par(mar=c(5, 8, 5, 10))
    barplot(table(df$n_genes), xlab="Genes per variant", ylab="Frequency", main="", cex.axis=1.25, cex.lab=1.5, cex.names=1.25)

    # Plot 3 --  histogram of distribution per chromosome
    par(mar=c(4, 8, 3, 4))
    # Get 22 colors from a palette
    colorz.chr = viridis(22, option = "plasma")
    tmp <- as.data.frame(table(variant_annotation$Chromosome))
    tmp$Var1 <- as.numeric(as.character(tmp$Var1))
    if (nrow(tmp) != 22){ for (i in 1:22){ if (!(i %in% tmp$Var1)){ t <- data.frame(Var1=i, Freq=0); tmp <- rbind(tmp, t) } } }
    tmp <- tmp[order(tmp$Var1),]
    barplot(tmp$Freq, col=colorz.chr, names=tmp$Var1, xlab="Chromosome", ylab="Genes per chromosome", cex.lab=1.50, main="", cex.axis=1.4, cex.main=2)

    # Plot 4 -- circular plot integrated visualization
    # read chromosome length
    l <- data.frame(chrom = paste0("chr", c(1:23)), length = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 50000000))

    # circle circumference is the sum of all chromosome lengths
    R <- sum(l$length)/1000000

    # then derive radius
    r <- R/(2*pi)

    # find corresponding angles for the chromosomes
    ang <- (l$length/1000000)/r

    # set colors
    colorz <- colorz.chr
    col.bk <- rep(c("white", "grey90"), 11)

    # set pch
    variant_annotation$pch <- 21
    variant_annotation$pch[which(variant_annotation$Source == "coding")] <- 24
    variant_annotation$pch[which(variant_annotation$Source == "qtl_query")] <- 23
    variant_annotation$pch[which(variant_annotation$Source == "coding_ld")] <- 25
    variant_annotation$pch[which(variant_annotation$Source == "qtl_ld")] <- 22

    par(mar=c(4, 18, 4, 18))
    plot(0, 0, pch=16, col="white", xlim=c(-r, r), ylim=c(-r, r), bty='n', xlab="", ylab="", xaxt='none', yaxt="none")
    plotrix::draw.circle(x = 0, y = 0, radius = r, border=NA)
    cum.angle <- pi/2-(ang[1]/2)
    # First loop to add sectors
    for (i in 1:23){
        # Get degrees
        dg1 <- cum.angle*180/pi
        dg2 <- (cum.angle+ang[i])*180/pi
        # Draw background and borders
        circlize::draw.sector(start.degree = dg1, end.degree = dg2, rou1 = r, rou2 = 0, clock.wise = FALSE, col=ggplot2::alpha(col.bk[i], 0.6), border = NA)
        cum.angle <- cum.angle + ang[i]
    }

    # Add axis for maf
    a = r*0.85
    b = r*0.30
    d = r*0.575
    c = r*0.01
    segments(x0 = 0, y0 = b, x1 = 0, y1 = a, lwd=1.25, xpd=T)
    # Add maf circles for reference
    plotrix::draw.circle(x = 0, y = 0, radius = a-a*0.025, border="grey40", lty=2)
    plotrix::draw.circle(x = 0, y = 0, radius = b-b*0.025, border="grey40", lty=2, nv = 1000)
    plotrix::draw.circle(x = 0, y = 0, radius = d-d*0.025, border="grey40", lty=2, nv = 1000)
    # Add axis ticks
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

    # Second loop to add points
    cum.angle <- pi/2+(ang[1]/2)
    for (i in 1:22){
        # Calculate points for the segment of the lolliplot -- point0 should be maf
        pp <- variant_annotation[which(variant_annotation$Chromosome == i),]
        pp$"Alternative Allele Frequency"[which(pp$"Alternative Allele Frequency" >0.5)] <- 1-pp$"Alternative Allele Frequency"[which(pp$"Alternative Allele Frequency" >0.5)]

        # For the frequency, need to normalize it in [r*0.20-r*0.80] [b, a]
        pp$maf_norm <- (b-a)*(pp$"Alternative Allele Frequency"/0.5)+a

        # Get the angle
        pp$ang <- (as.numeric(pp$"Position (hg38)")/1000000)/r + cum.angle
        pp$x <- pp$maf_norm*cos(pp$ang)
        pp$y <- pp$maf_norm*sin(pp$ang)
        segments(x0 = pp$x, y0 = pp$y, x1 = r*cos(pp$ang), y1 = r*sin(pp$ang), lwd=0.85, col="grey40")
        points(x = pp$x, y = pp$y, type = "p", pch=pp$pch, col="black", bg=colorz[i], cex=1.75)

        # Draw sector finally
        # Get degrees
        dg1 <- cum.angle*180/pi
        dg2 <- (cum.angle+ang[i+1])*180/pi
        circlize::draw.sector(start.degree = dg1, end.degree = dg2, rou1 = r, rou2 = r-r*0.1, clock.wise = FALSE, col=colorz[i])

        # Add label
        rgm <- cum.angle+(ang[i+1]/2)
        text(x = r*cos(rgm)*0.95, y = r*sin(rgm)*0.95, labels = i, font=2, adj=0.5)

        # Increment at the end
        cum.angle <- cum.angle + ang[i+1]
    }
    # Finally the legend
    legend(x = -r/2, y = -r*1.05, xpd=T, legend = c("Coding", "QTL", "Coding LD", "QTL LD", "Position"), pch=c(23, 24, 25, 22, 21), pt.bg="grey40", pt.cex = 1.50, cex=1.50, bty='n', ncol=3)

    return(df)
    }

# Function to plot GWAS annotation
plotGWAS_annotation = function(gwas_annot, output_folder){
    # Read clustering of GWAS
    gwas_clustering = fread('../../Data/databases/haplotypes/trait_clusters_summary_AI_thresholds.csv', h=TRUE, sep=",", stringsAsFactors=FALSE)
    #Subset to threshold 0.5
    gwas_clustering = gwas_clustering[which(gwas_clustering$Threshold == 0.5),]
    # Add a column to gwas_annot with the "Final Representative" from gwas_clustering based on grepping the Trait in gwas_clustering$Members
    if (nrow(gwas_annot) ==0){
        return(NULL)
    }
    gwas_annot$Final_Representative = NA
    for (i in 1:nrow(gwas_annot)){
        tmp_trait = gwas_annot$"GWAS Trait"[i]
        match_row <- gwas_clustering[grepl(tmp_trait, gwas_clustering$Members, fixed = TRUE),]
        if (nrow(match_row) >0){
            gwas_annot$Final_Representative[i] = match_row$"Final Representative"[1]
        }
    }
    # Make dataframe of frequencies
    freq_df = as.data.frame(table(gwas_annot$Final_Representative))
    colnames(freq_df) = c("Trait", "Frequency")
    # sort by frequency
    freq_df = freq_df[order(-freq_df$Frequency),]
    # Take top 20
    top_n = 20
    if (nrow(freq_df) < top_n){ top_n = nrow(freq_df) }
    top_freq_df = freq_df[1:top_n,]
    # Exclude parts in parenthesis in the trait names for better visualization
    top_freq_df$Trait_Short = gsub(" \\(.*\\)", "", top_freq_df$Trait)
    # If more than 30 characters, truncate
    top_freq_df$Trait_Short = ifelse(nchar(top_freq_df$Trait_Short) > 30, paste0(substr(top_freq_df$Trait_Short, 1, 27), "..."), top_freq_df$Trait_Short)
    # Reverse order for better visualization
    top_freq_df = top_freq_df[order(-top_freq_df$Frequency),]
    # Set palette
    palett = viridis(n = nrow(top_freq_df), option = 'turbo')
    # Set margins
    # Plot barplot
    png(paste0(output_folder, "/gwas_annotation_top_traits.png"), height=7, width=10, res=300, units="in")
    par(mar=c(12, 5, 4, 2))
    barplot(top_freq_df$Frequency, names.arg = top_freq_df$Trait_Short, las=2, col=palett, main="Top GWAS traits in variant annotation", ylab="Frequency", cex.names=0.85)
    dev.off()

}

# Place arguments in variables
input_file <- args$i
output_folder <- args$o
gwas_annotation <- args$g
# input_file = '/Users/nicco/Downloads/snpXplorer_annotation_123456/variant_annotation_combined.tsv'
# output_folder = '/Users/nicco/Downloads/snpXplorer_annotation_123456/'
# gwas_annotation = '/Users/nicco/Downloads/snpXplorer_annotation_123456/target_annotations/gwas_annotation_target.tsv'

# Create a plot folder
tryCatch({
    plot_folder <- paste0(output_folder)
    dir.create(plot_folder, showWarnings = FALSE, recursive = TRUE)
}, error = function(e) {
    message("Error in creating output folder: ", e)
})

# Plot variant annotation
tryCatch({
    # Load variant annotation
    variant_annotation <- fread(input_file, h=TRUE, sep="\t", stringsAsFactors=FALSE)
    # Plot variant annotation summary
    pdf(paste(output_folder, "/variant_annotation_summary.pdf", sep=""), height=12.35, width=10)
    df_mapping <- plotMapping(variant_annotation)
    invisible(dev.off())
}, error = function(e) {
    message("Error in plotting variant annotation summary: ", e)
})

# Load GWAS annotation
tryCatch({
    gwas_annot <- fread(gwas_annotation, h=TRUE, sep="\t", stringsAsFactors=FALSE)
    # Plot GWAS annotation
    plotGWAS_annotation(gwas_annot, output_folder)
}, error = function(e) {
    message("No GWAS annotation file provided or file could not be read.")
})
