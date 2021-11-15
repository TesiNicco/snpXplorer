# basic paths
MAIN = "/root/snpXplorer/AnnotateMe/"
MAIN_SNP = "/root/snpXplorer/snpXplorer_v3/"
args = commandArgs(trailingOnly=TRUE)

## function to perform alternative version of revigo -- do semantic similarity myself
alternative_revigo_results <- function(MAIN, random_num, go_data){
  # check if revigo input is there
  filelist = system(paste0("ls RESULTS_", random_num, "/"), intern = T)
  if ("revigo_inp.txt" %in% filelist){

    # run the python script
    cmd = paste0("python3 ", MAIN, "BIN/Alternative_REVIGO.py RESULTS_", random_num, "/revigo_inp.txt RESULTS_", random_num, "/alternative_Lin_distance.txt ")
    system(cmd)

    # read output back
    data = data.table::fread(paste0("RESULTS_", random_num, "/alternative_Lin_distance.txt"), h=T)

    # Plot heatmap and save it
    pdf(paste0("RESULTS_", random_num, "/pheatmap_lin_distance.pdf"), height=7, width=7)
    pheatmap::pheatmap(data, show_rownames = F, show_colnames = F, main="Semantic similarity matrix of GO terms")
    dev.off()

    hr1 <- hclust(as.dist(1-data), method="ward.D2", members = NULL)

    # loop until at least 3 clusters are found
    n_clust = 0
    minModuleSize = 15
    while(n_clust < 2){
      dynamicMods = dynamicTreeCut::cutreeDynamic(dendro = hr1, distM = data,
                                  minClusterSize = minModuleSize, deepSplit = TRUE,
                                  verbose = 0, respectSmallClusters = T)
      tb <- as.data.frame(table(dynamicMods))
      n_clust <- max(as.numeric(as.character(tb$dynamicMods)))
      df_cluster_terms = data.frame(term = colnames(data), cluster = dynamicMods)
      df_cluster_terms$cluster_in_dendro = NA
      minModuleSize = minModuleSize - 1
    }

    # Hierarchical clustering on the distance matrix
    hr <- as.dendrogram(hclust(as.dist(1-data), method="ward.D2", members = NULL))
    hc <- hclust(as.dist(1-data), method="ward.D2", members = NULL)
    order_dendro <- hc$order
    cluster_number_dendro = 1
    for (i in 1:length(order_dendro)){
      # get element
      tmp_go = colnames(data)[order_dendro[i]]
      # check whether it was already annotated
      if (is.na(df_cluster_terms$cluster_in_dendro[which(df_cluster_terms$term == tmp_go)])){
        # get previous cluster number
        prev_cluster = df_cluster_terms$cluster[which(df_cluster_terms$term == tmp_go)]
        df_cluster_terms$cluster_in_dendro[which(df_cluster_terms$cluster == prev_cluster)] <- cluster_number_dendro
        cluster_number_dendro <- cluster_number_dendro + 1
      }
    }
    df_cluster_terms$term <- as.character(df_cluster_terms$term)
    df_cluster_terms = df_cluster_terms[match(hr1$labels[hr1$order], df_cluster_terms$term),]

    if (nrow(data) >= 120){
      plt_w = 28
    } else {
      plt_w = 17
    }
    # Plot dendrogram and rectangles of the clusters
    df_cluster_terms$col <- NA
    library(viridis)
    library(RColorBrewer)
    palett = brewer.pal(n = max(df_cluster_terms$cluster_in_dendro)+1, name = "Set1")[1:max(df_cluster_terms$cluster_in_dendro)]
    for (i in 1:max(df_cluster_terms$cluster_in_dendro)){
      df_cluster_terms$col[which(df_cluster_terms$cluster_in_dendro == i)] <- palett[i]
    }
    library(dendextend)
    png(paste0("RESULTS_", random_num, "/dendrogram_GOterms.png"), height=7, width=plt_w, res=300, units="in")
    d1=color_branches(dend = hr, clusters = df_cluster_terms$cluster_in_dendro, groupLabels = T, col = palett)
    d1 %>% dendextend::set("labels_col", df_cluster_terms$col) %>% dendextend::set("branches_lwd", 2) %>% dendextend::set("labels_cex", 0.80) %>% plot(main = paste0("Dynamic Cut-tree algorithm ~ Ward.D2 ~ min. size=", minModuleSize))
    dev.off()

    # Merge with terms description
    go_data_desc = go_data[, c("term_name", "term_id")]
    clusters <- merge(df_cluster_terms, go_data_desc, by.x="term", by.y="term_id")
    clusters$col <- NULL
  } else {
    clusters = NA
    data = NA
    n_clust = NA
  }
  # return
  return(list(clusters, data, n_clust))
}

# function to count the number of occurrences of each word given a vector of sentences
CountFrequency_words <- function(functional_clusters, n_clust){
  # Read the description of all GO:BP -- this will be the background to remove redundant words
  all_go_bp <- data.table::fread("/root/snpXplorer/AnnotateMe/BIN/go_terms_BP_all.txt", h=F, stringsAsFactors = F, sep="\t")

  colnames(functional_clusters) <- c("term", "cluster", "cluster_in_dendro", "term.y")
  # need to remove the (GO id) from the term name
  # for (i in 1:nrow(functional_clusters)){
  #   functional_clusters$term.y[i] <- unlist(strsplit(functional_clusters$term.y[i], "[()]"))[1:(length(unlist(strsplit(functional_clusters$term.y[i], "[()]"))) -1)]
  # }

  # Put all different words datasets in a list
  all_dset <- list()
  for (i in 1:n_clust){
    if (i == 1){
      all_dset[[i]] <- all_go_bp$V1
      all_dset[[i+1]] <- functional_clusters$term.y[which(functional_clusters$cluster_in_dendro == i)]
    } else {
      all_dset[[i+1]] <- functional_clusters$term.y[which(functional_clusters$cluster_in_dendro == i)]
    }
  }

  # Define output
  word_frequency_dset <- list()

  # Main loop to calculate word frequency
  for (i in 1:length(all_dset)){
    # Convert to dataframe and make sure there are characters only
    text <- data.frame(text=all_dset[[i]])
    text$text <- as.character(text$text)

    # Covert dataframe to tibble
    text_df <- tibble(line = 1:nrow(text), text = text)
    text_df <- mutate(text_df, text = text$text)

    # Count word occurrence -- for some reasons this does not work anymore (??)
    #x <- text_df %>% unnest_tokens(word, text) %>%    # split words
    #  anti_join(stop_words) %>%    # take out "a", "an", "the", etc.
    #  count(word, sort = TRUE)    # count occurrences
    # Convert back to dataframe
    # x <- as.data.frame(x)
    # colnames(x) <- c("word", "freq")
    # Here's the alternative
    x <- text_df %>% unnest_tokens(word, text) %>%
      anti_join(stop_words)
    frq = as.data.frame(table(x$word))
    frq = frq[order(-frq$Freq),]

    # Assign to output
    word_frequency_dset[[i]] <- frq
  }

  # Extract background and calculate % and take top 5%
  bkg <- word_frequency_dset[[1]]
  bkg$perc <- bkg$Freq/nrow(all_go_bp)
  top_5_pc <- bkg[which(bkg$perc >= 0.025),]

  # Now exclude these frequent words from the clusters names
  for (i in 2:length(word_frequency_dset)){
    sb <- word_frequency_dset[[i]]
    sb <- sb[which(!(sb$Var1 %in% top_5_pc$Var1)),]
    word_frequency_dset[[i]] <- sb
  }

  # Try some plot -- seems like it is not possible to plot multiple plots in one --> make 4 plots and save them from "Export"
  for (i in 2:length(word_frequency_dset)){
    myplot <- wordcloud2(data = word_frequency_dset[[i]])
    # save it in html
    saveWidget(myplot, "tmp.html", selfcontained = F)
    # conert to png
    webshot("tmp.html", paste0("RESULTS_", random_num, "/cluster_", i-1, "_wordcloud.png"), delay = 20, vwidth = 1000, vheight = 1000)
  }

  return(word_frequency_dset)
}

# READ ARGUMENTS AND RUN FUNCTION
suppressPackageStartupMessages({
    library(htmlwidgets)
    library(wordcloud2)
    library(tidytext)
    library(dendextend)
    library(dplyr)
    library(webshot)
    library(data.table)
    library(tibble)
})
random_num = args[1]
load(paste0("RESULTS_", random_num, "/tmp_enrichRes.RData"))
go_data = sampling.res[[2]]
alt_revigo_res <- function(MAIN, random_num, go_data) {
    out <- tryCatch({ semsim_results = alternative_revigo_results(MAIN, random_num, go_data) },
            error=function(cond) {
                message("GO terms analysis encountered an error..Skipping analysis!")
                # Choose a return value in case of error
                return(NA)
            })
    return(out)
}
final_res = alt_revigo_res(MAIN, random_num, go_data)
# FINALLY IN CASE WE HAVE RESULTS DO THE WORDCLOUDS
if (length(final_res) == 1){
    functional_clusters = NULL
    lin_matrix = NULL
    n_clust = 0
} else {
    functional_clusters = final_res[[1]]
    lin_matrix <- final_res[[2]]
    n_clust = final_res[[3]]
}
if (length(functional_clusters) == 1){
    write.table(go_data, file = paste0("RESULTS_", random_num, "/geneSet_enrichment_results_and_clusters.txt"), quote=F, row.names=F, sep="\t")
} else if (length(functional_clusters) > 1) {
    # save the whole gene-set enrichment analysis
    tmp = functional_clusters
    tmp$term_name = NULL
    go_data = go_data[!duplicated(go_data$term_id),]
    go_data = merge(go_data, tmp, by.x = "term_id", by.y = "term", all.x = T)
    go_data = go_data[order(go_data$avg_p),]
    go_data$cluster = NULL
    write.table(go_data, file = paste0("RESULTS_", random_num, "/geneSet_enrichment_results_and_clusters.txt"), quote=F, row.names=F, sep="\t")

    # let's try to make some wordcloud images
    if (!is.na(functional_clusters)){
        mostFreq_words <- CountFrequency_words(functional_clusters, n_clust)
    }
}
