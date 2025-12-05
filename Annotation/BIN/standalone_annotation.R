# Rscript to do clustering and dynamic cut of GO terms based on semantic similarity

# Libraries
# Disable starting messages
suppressPackageStartupMessages(library(dynamicTreeCut))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(pheatmap))

# Argument parser
parser <- ArgumentParser()
parser$add_argument("-d", "--dist_matrix", help="Path to the semantic similarity distance matrix file (TSV format)", required=TRUE)
parser$add_argument("-o", "--output-folder", help="Output folder to save dendrogram plots", required=TRUE)
# Parse arguments
args <- parser$parse_args()

# Functions
# Clustering function
clustering_function = function(data, max_n_clust, hc, hr){
    # Set initial parameters
    n_clust = 0
    minModuleSize = 15
    # Hierarchical clustering on the distance matrix
    while (n_clust < 2 | n_clust > max_n_clust){
        dynamicMods = dynamicTreeCut::cutreeDynamic(dendro = hc, distM = data, minClusterSize = minModuleSize, deepSplit = 0, verbose = 0, respectSmallClusters = T)
        tb <- as.data.frame(table(dynamicMods))
        n_clust <- max(as.numeric(as.character(tb$dynamicMods)))
        df_cluster_terms = data.frame(term = colnames(data), cluster = dynamicMods)
        df_cluster_terms$cluster_in_dendro = NA
        if (n_clust > max_n_clust){
            minModuleSize = minModuleSize + ceiling(minModuleSize*0.10)
        } else {
            minModuleSize = minModuleSize - 1
        }
    }
    # Get cluster in dendrogram order
    order_dendro <- hc$order
    cluster_number_dendro = 1
    # Iterate over ordered terms to assign cluster in dendrogram
    for (i in 1:length(order_dendro)){
        # Get element
        tmp_go = colnames(data)[order_dendro[i]]
        # Check whether it was already annotated
        if (is.na(df_cluster_terms$cluster_in_dendro[which(df_cluster_terms$term == tmp_go)])){
            # Get previous cluster number
            prev_cluster = df_cluster_terms$cluster[which(df_cluster_terms$term == tmp_go)]
            # Assign new cluster number to all terms in the previous cluster
            df_cluster_terms$cluster_in_dendro[which(df_cluster_terms$cluster == prev_cluster)] <- cluster_number_dendro
            cluster_number_dendro <- cluster_number_dendro + 1
        }
    }
    df_cluster_terms$term <- as.character(df_cluster_terms$term)
    df_cluster_terms = df_cluster_terms[match(hc$labels[hc$order], df_cluster_terms$term),]
    return(df_cluster_terms)
}

# Function to plot dendrogram
plot_dendrogram <- function(df_cluster_terms, data, hd, output_folder){    
    # Determine plot width
    if (nrow(data) >= 120){
        plt_w = 28
    } else {
        plt_w = 17
    }
    # Define color palette
    palett = viridis(n = max(df_cluster_terms$cluster_in_dendro)+1, option = 'turbo')[1:max(df_cluster_terms$cluster_in_dendro)]
    # Assign colors to clusters
    df_cluster_terms$col <- NA
    for (i in 1:max(df_cluster_terms$cluster_in_dendro)){
        df_cluster_terms$col[which(df_cluster_terms$cluster_in_dendro == i)] <- palett[i]
    }
    # Derive number of clusters for title
    n_clust = max(df_cluster_terms$cluster_in_dendro)
    # Plot dendrogram and rectangles of the clusters
    outname = paste0(output_folder, "_dendrogram_clusters.png")
    png(outname, height=7, width=plt_w, res=300, units="in")
    d1 = color_branches(dend = hd, clusters = df_cluster_terms$cluster_in_dendro, groupLabels = T, col = palett)
    d1 %>% dendextend::set("labels_col", df_cluster_terms$col) %>% dendextend::set("branches_lwd", 2) %>% dendextend::set("labels_cex", 0.80) %>% plot(main = paste0("Dynamic Cut-tree algorithm ~ Ward.D2 ~ Clusters=", n_clust))
    dev.off()
}

# Arguments to variables
dist_matrix_file <- args$dist_matrix
output_folder <- args$output_folder

# Read distance matrix
data <- as.matrix(read.table(dist_matrix_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE))

# Plot heatmap and save it
pdf(paste0(output_folder, "/pheatmap_lin_distance.pdf"), height=7, width=7)
pheatmap::pheatmap(data, show_rownames = F, show_colnames = F, main="Semantic similarity matrix of GO terms")
dev.off()

# Hierarchical clustering on the distance matrix
hc <- hclust(as.dist(1-data), method="ward.D2", members = NULL)
hd <- as.dendrogram(hclust(as.dist(1-data), method="ward.D2", members = NULL))

# Dynamic clustering loop until clusters are found
# Set different thresholds for how aggressive the clustering should be
aggressive_thresholds <- c(5, 8, 10, 12)
res_list = list()
for (max_n_clust in aggressive_thresholds){
    # Define name for the outputs
    output_prefix = paste0(output_folder, "/clustering_max_", max_n_clust)
    print(paste0("Clustering with max ", max_n_clust, " clusters"))
    df_cluster_terms <- clustering_function(data, max_n_clust, hc, hd)
    # Plot dendrogram
    plot_dendrogram(df_cluster_terms, data, hd, output_prefix)
    # Save clustering results
    write.table(df_cluster_terms, file=paste0(output_prefix, "_clusters.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
}
