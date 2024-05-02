# Libraries
    library(stringr)
    library(gprofiler2)
    library(pheatmap)
    library(dynamicTreeCut)
    library(viridis)
    library(RColorBrewer)
    library(dendextend)
    library(tibble)
    library(htmlwidgets)
    library(dplyr)
    library(wordcloud2)
    library(tidytext)
    library(webshot)
    
# Functions
    # Function to create n sampling datasets
    samplingGset <- function(i, mapping){
        all.genes <- strsplit(mapping$geneList, ",")
        tmp_f <- function(x, all.genes){ tmp <- all.genes[[x]]; g <- sample(x=tmp, size=1); return(g)}
        gset <- lapply(1:length(all.genes), tmp_f, all.genes=all.genes)
        gset.clean <- unlist(gset)
        return(gset.clean)
    }

    # Function for gene-set enrichment analysis
    overlapAnalysis <- function(i, gene_list, source_gset){
        print(i)
        g <- unique(gene_list[[i]])
        res <- suppressMessages(gprofiler2::gost(g, organism = "hsapiens", ordered_query = FALSE, multi_query = FALSE, significant = FALSE, exclude_iea = TRUE,                             measure_underrepresentation = FALSE, evcodes = FALSE, user_threshold = 1, correction_method = "fdr", domain_scope = "annotated", custom_bg = NULL, numeric_ns = "", sources = source_gset))
        return(res$result)
    }

    # Function for merging gene-set enrichment analysis
    mergeSampling <- function(enrich.res){
        n.sampl = length(enrich.res)
        all.enrich <- do.call("rbind", enrich.res)      #put all results together
        sign.only <- all.enrich[which(all.enrich$p_value <= 0.05),]         #isolate only significant findings
        if (nrow(sign.only) >0){
            sign.table <- as.data.frame(table(sign.only$term_name))
            info <- sign.only[, c("term_name", "term_id")]
            sign.table <- sign.table[order(-sign.table$Freq),]
            sign.table <- merge(sign.table, info, by.x="Var1", by.y="term_name")
        } else { sign.table <- data.frame(chr=NULL, position=NULL)}

        # also try with averaging pvalues
        if ("intersection" %in% colnames(all.enrich)){
            terms <- all.enrich[, c("term_name", "term_id", "intersection_size", "p_value")]
            slow <- TRUE
        } else {
            terms <- all.enrich[, c("term_name", "term_id", "p_value")]
            slow = FALSE
        }
        # see how many times each term was tested and resctrict to common (tested all times)
        ttt = as.data.frame(table(terms$term_id))
        ttt <- ttt[which(ttt$Freq == n.sampl),]

        # if for some reason the number of samplings was not 500 (e.g. small number of genes), need to adjust for this
        if (nrow(ttt) == 0){
            ttt = as.data.frame(table(terms$term_id))
            mx = max(ttt$Freq)
            ttt <- ttt[which(ttt$Freq == mx),]
        }

        terms <- terms[which(terms$term_id %in% ttt$Var1),]
        #terms <- terms[!duplicated(terms$term_name),]
        avgP <- function(i, terms, term_list, slow){
            if (isTRUE(slow)){
            #all_intersections <- paste0(unique(unlist(strsplit(paste0(terms$intersection[which(terms$term_id == term_list[i])], collapse = ","), ","))), collapse = ",")
            #df <- data.frame(term_name = unique(terms$term_name[which(terms$term_id == term_list[i])]), term_id = term_list[i], avgP = mean(terms$p_value[which(terms$term_id == term_list[i])]), intersection = all_intersections)
            df <- data.frame(term_name = unique(terms$term_name[which(terms$term_id == term_list[i])]), term_id = term_list[i], avgP = mean(terms$p_value[which(terms$term_id == term_list[i])]))
            } else {
            df <- data.frame(term_name = unique(terms$term_name[which(terms$term_id == term_list[i])]), term_id = term_list[i], avgP = mean(terms$p_value[which(terms$term_id == term_list[i])]))
            }
            return(df)
        }
        term_list <- unique(terms$term_id)
        all.terms.avg <- data.table::rbindlist(parallel::mclapply(1:length(term_list), avgP, terms=terms, term_list = term_list, slow = slow, mc.cores=1))
        all.terms.avg <- all.terms.avg[order(all.terms.avg$avgP),]
        all.terms.avg$log10P <- -log10(all.terms.avg$avgP)

        # if intersection were requested, add them back here
        if (slow == TRUE){
            all.terms.avg$intersection = NA
            for (i in 1:nrow(all.terms.avg)){
            cat(paste0('** processing line ', i, '                         \r'))
            all.terms.avg$intersection[i] = paste(unique(unlist(strsplit(paste(all.enrich$intersection[which(all.enrich$term_id == all.terms.avg$term_id[i])], collapse = ','), ','))), collapse = ',')
            }
        }

        #print(head(all.terms.avg))
        #write.table(all.terms.avg, paste("RESULTS_", random_num, "/enrichent_results_sampling.txt", sep=""), quote=F, row.names=F, sep="\t")

        # if (nrow(all.terms.avg[which(all.terms.avg$avgP <= 0.05),]) >0){
        #   pdf(paste("RESULTS_", random_num, "/treemap_sampling_PRS_genes.pdf", sep=""), height=10, width=10)
        #   treemap(all.terms.avg[which(all.terms.avg$avgP <= 0.05),], index="term_name", vSize="log10P", type="index")
        #   invisible(dev.off())
        # }

        l <- list(sign.table, all.terms.avg)
        return(l)
    }

    # Function to plot enrichment results other than go
    plotOtherEnrichments <- function(sampling.res, dbs, random_num, type){
        if (type == "enrichR"){
            library(enrichplot)
            library(viridis)
            # exclude GO
            pdf(paste0("RESULTS_", random_num, "/Enrichment_results_noGO.pdf"), height = 12, width = 12, onefile = T)
            par(mar = c(6, 16, 4, 4))
            for (x in 1:length(dbs)){
            if (dbs[x] != "GO_Biological_Process_2018"){
                tmp = sampling.res[[x]]
                tmp_sig <- tmp[which(tmp$avg_p <= 0.10),]
                tmp_sig$logP <- -log10(as.numeric(tmp_sig$avg_p))
                max_p <- max(tmp_sig$logP, 3)
                pos <- barplot(tmp_sig$logP, horiz = T, col=viridis(n = nrow(tmp_sig), option = "plasma"), xlab = "-log10(p)", xlim=c(0, max_p), main = tmp_sig$source[1], cex.lab=1.60, cex.axis = 1.25)
                abline(v = -log10(0.05), lty=2, col="red", lwd=2)
                abline(v = -log10(0.10), lty=2, col="deepskyblue3", lwd=2)
                for (i in 1:length(pos)){
                # look into spacing
                n.lines <- ceiling(nchar(tmp_sig$term[i])/32)

                # divide by space depending on n.lines
                tmp <- unlist(strsplit(as.character(tmp_sig$term[i]), " "))
                if (n.lines >1){
                    kk <- split(tmp, ceiling(seq_along(tmp)/n.lines))
                    for (j in 1:length(kk)){ kk[[j]] <- paste(kk[[j]], collapse=" ") }
                    # add name to data frame
                    tmp_sig$term[i] <- paste(kk, collapse="\n")
                } else {
                    tmp_sig$term[i] <- paste(tmp, collapse=" ")
                }
                text(x = 0, y = pos[i, 1], labels = tmp_sig$term[i], xpd=T, pos=2, offset = 1)
                legend("topright", legend = c("FDR<10%", "FDR<5%"), lty = c(2,2), col=c("deepskyblue3", "red"), cex=1.5, ncol=1, bty="n", lwd=3)
                }
            }
            }
            dev.off()
        } else if (type == "gost"){
            library(viridis)
            # determine source of enrichment
            tmp <- sampling.res[[2]]
            tmp$source_gset <- str_split_fixed(tmp$term_id, ":", 2)[, 1]
            # exclude GO
            tmp <- tmp[which(tmp$source_gset != "GO"),]
            pdf(paste0("RESULTS_", random_num, "/Enrichment_results.pdf"), height = 12, width = 12, onefile = T)
            for (x in unique(tmp$source_gset)){
            dt = tmp[which(tmp$source_gset == x),]
            tmp_sig <- dt[which(dt$avgP <= 0.10),]
            if (nrow(tmp_sig) >=20){
                if (nrow(tmp_sig) > 45){
                tmp_sig <- tmp_sig[order(tmp_sig$avgP),]
                tmp_sig <- head(tmp_sig, 45)
                }
                par(mar = c(6, 17, 4, 4))
                cex_text = 0.50
                splt_words = 100
            } else {
                par(mar = c(6, 17, 4, 4))
                cex_text = 1
                splt_words = 40
            }
            max_p <- max(ceiling(tmp_sig$log10P), 3)
            if (nrow(tmp_sig) >0){
                pos <- barplot(tmp_sig$log10P, horiz = T, col=viridis(n = nrow(tmp_sig), option = "plasma"), xlab = "-log10(p)", xlim=c(0, max_p), main = dt$source_gset[1], cex.lab=1.60, cex.axis = 1.25)
                abline(v = -log10(0.05), lty=2, col="red", lwd=2)
                abline(v = -log10(0.10), lty=2, col="deepskyblue3", lwd=2)
                tmp_sig$term_name <- as.character(tmp_sig$term_name)
                for (i in 1:length(pos)){
                # look into spacing
                n.lines <- ceiling(nchar(tmp_sig$term_name[i])/splt_words)

                # divide by space depending on n.lines
                tmp_w <- unlist(strsplit(as.character(tmp_sig$term_name[i]), " "))
                if (n.lines >1){
                    pieces = floor(length(tmp_w)/n.lines)
                    if (n.lines == 2){
                    p1 = paste(tmp_w[1:pieces], collapse = " ")
                    p2 = paste(tmp_w[(pieces+1):length(tmp_w)], collapse = " ")
                    tmp_sig$term_name[i] <- paste(p1, p2, sep="\n")
                    } else if (n.lines == 3){
                    p1 = paste(tmp_w[1:pieces], collapse = " ")
                    p2 = paste(tmp_w[(pieces+1):(pieces*2)], collapse = " ")
                    p3 = paste(tmp_w[(pieces*2+1):length(tmp_w)], collapse = " ")
                    tmp_sig$term_name[i] <- paste(p1, p2, p3, sep="\n")
                    }
                    # kk <- split(tmp_w, ceiling(seq_along(tmp_w)/n.lines))
                    # for (j in 1:length(kk)){ kk[[j]] <- paste(kk[[j]], collapse=" ") }
                    # # add name to data frame
                    # tmp_sig$term_name[i] <- paste(kk, collapse="\n")
                } else {
                    tmp_sig$term_name[i] <- paste(tmp_w, collapse=" ")
                }
                text(x = 0, y = pos[i, 1], labels = tmp_sig$term_name[i], xpd=T, pos=2, offset = 1, cex = cex_text)
                legend("topright", legend = c("p<0.05", "p<0.10"), lty = c(2, 2), col=c("red", "deepskyblue3"), cex=1.5, ncol=2, bty="n", lwd=3)
                }
            }else {
                plot(0,0, pch=16, col="white", xlim=c(0, 1), ylim=c(0,1), main=dt$source_gset[1], xlab="", xaxt="none", yaxt="none", ylab="", bty="n")
                text(x = 0.5, y = 0.5, labels = "No significant hits")
            }
            }
            dev.off()
        }
    }

    # Function for REVIGO alternative analysis
    semanticSimilarityAnalysis <- function(geneset, max_cluster_n){
        # take significant genesets
        geneset_sig = geneset[which(geneset$avgP <= 0.05),]
        if (nrow(geneset_sig) >0){
            write.table(geneset_sig[, c('term_id', 'avgP')], 'revigo_inp.txt', quote=F, row.names=F, col.names=F, sep="\t")
            # run the python script
            cmd = paste0("python3 /project/holstegelab/Share/nicco/snpxplorer/snpXplorer/Annotation/BIN/Alternative_REVIGO.py revigo_inp.txt alternative_Lin_distance.txt indep")
            system(cmd)
            # read output back
            data = data.table::fread("alternative_Lin_distance.txt", h=T)
            # Plot heatmap and save it
            pdf("pheatmap_lin_distance.pdf", height=7, width=7)
            pheatmap(data, show_rownames = F, show_colnames = F, main="Semantic similarity matrix of GO terms")
            dev.off()
            hr1 <- hclust(as.dist(1-data), method="ward.D2", members = NULL)
            # loop until at least 3 clusters are found
            n_clust = 0
            minModuleSize = 15
            while (n_clust < 2 | n_clust > max_cluster_n){
                dynamicMods = cutreeDynamic(dendro = hr1, distM = data, minClusterSize = minModuleSize, deepSplit = TRUE, verbose = 0, respectSmallClusters = T)
                tb <- as.data.frame(table(dynamicMods))
                n_clust <- max(as.numeric(as.character(tb$dynamicMods)))
                df_cluster_terms = data.frame(term = colnames(data), cluster = dynamicMods)
                df_cluster_terms$cluster_in_dendro = NA
                if (n_clust > max_cluster_n){
                    minModuleSize = minModuleSize + ceiling(minModuleSize*0.10)
                } else {
                minModuleSize = minModuleSize - 1
                }
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
            # graphical parameters for the plot
            if (nrow(data) >= 120){
                plt_w = 28
            } else {
                plt_w = 17
            }
            # Plot dendrogram and rectangles of the clusters
            df_cluster_terms$col <- NA
            palett = viridis(n = max(df_cluster_terms$cluster_in_dendro)+1, option = 'turbo')[1:max(df_cluster_terms$cluster_in_dendro)]
            for (i in 1:max(df_cluster_terms$cluster_in_dendro)){
                df_cluster_terms$col[which(df_cluster_terms$cluster_in_dendro == i)] <- palett[i]
            }
            png("dendrogram_GOterms.png", height=7, width=plt_w, res=300, units="in")
            d1=color_branches(dend = hr, clusters = df_cluster_terms$cluster_in_dendro, groupLabels = T, col = palett)
            d1 %>% dendextend::set("labels_col", df_cluster_terms$col) %>% dendextend::set("branches_lwd", 2) %>% dendextend::set("labels_cex", 0.80) %>% plot(main = paste0("Dynamic Cut-tree algorithm ~ Ward.D2 ~ min. size=", minModuleSize))
            dev.off()
            # Merge with terms description
            go_data_desc = geneset_sig[, c("term_name", "term_id")]
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
        all_go_bp <- data.table::fread("/project/holstegelab/Share/nicco/snpxplorer/snpXplorer/Annotation/BIN/go_terms_BP_all.txt", h=F, stringsAsFactors = F, sep="\t")
        colnames(functional_clusters) <- c("term", "cluster", "cluster_in_dendro", "term.y")
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
        # Plot
        for (i in 2:length(word_frequency_dset)){
            myplot <- wordcloud2(data = word_frequency_dset[[i]])
            # save it in html
            saveWidget(myplot, "tmp.html", selfcontained = F)
            # conert to png
            webshot("tmp.html", paste0("cluster_", i-1, "_wordcloud.png"), delay = 20, vwidth = 1000, vheight = 1000)
        }
        return(word_frequency_dset)
    }

# Arguments
    snps_info_path = 'RESULTS_299/snp_annotation.txt'               # Path to snp annotation data
    analysis_mode = 'Default'                                       # Gene-sets for analysis, can be comma separated
    analysis_mode = unlist(strsplit(analysis_mode, ','))
    n.sampl = 300                                                   # Number of iterations

# Read Annotation data
    annot = read.table(snps_info_path, h=T, stringsAsFactors = F)
    # derive gene list
    geneList = unlist(strsplit(annot$geneList, ','))

# Sample gene-sets
    gene_sampling_dsets <- parallel::mclapply(1:n.sampl, samplingGset, mapping=annot, mc.cores=3)

# Define databases for gene-set enrichment analysis
    dbs = c()
    for (i in analysis_mode){
        if (i %in% c("default", "Default")){
            dbs <- c(dbs, "GO:BP")
        } else if (i == "KEGG"){
            dbs <- c(dbs, "KEGG")
        } else if (i %in% c("wiki", "Wiki")){
            dbs <- c(dbs, "WP")
        } else if (i == "Reactome"){
            dbs <- c(dbs, "REAC")
        }
    }

# Enrichment analysis
    enrich.res <- parallel::mclapply(1:n.sampl, overlapAnalysis, gene_list=gene_sampling_dsets, source_gset = dbs, mc.cores=4)

# Merge sampling sets
    sampling.res <- mergeSampling(enrich.res)

# For other enrichment sources like kegg, reactome etc, need to make a plot
    if (length(dbs) >1){
        plotOtherEnrichments(sampling.res, dbs, random_num, type="gost")
    } else if (dbs[1] != "GO:BP"){
        plotOtherEnrichments(sampling.res, dbs, random_num, type = "gost")
    }

# Check if need to do semantic similarity analysis
    if (length(grep("GO", dbs)) >0){
        cat("## REVIGO analysis\n")
        # DO ALTERNATIVE TO REVIGO
        cat("## Semantic Similarity analysis\n")
        semsim_results = semanticSimilarityAnalysis(sampling.res[[2]], 15)   # adjust second argument to get fewer/more clusters. By default this is 15
        # if there are results, do the wordclouds
        # FINALLY IN CASE WE HAVE RESULTS DO THE WORDCLOUDS
        if (length(semsim_results) == 1){
            functional_clusters = NULL
            lin_matrix = NULL
            n_clust = 0
        } else {
            functional_clusters = semsim_results[[1]]
            lin_matrix <- semsim_results[[2]]
            n_clust = semsim_results[[3]]
        }
        if (length(functional_clusters) == 1){
            write.table(go_data, file = "geneSet_enrichment_results_and_clusters.txt", quote=F, row.names=F, sep="\t")
        } else if (length(functional_clusters) > 1){
            # save the whole gene-set enrichment analysis
            tmp = functional_clusters
            tmp$term_name = NULL
            go_data = sampling.res[[2]]
            go_data = go_data[!duplicated(go_data$term_id),]
            go_data = merge(go_data, tmp, by.x = "term_id", by.y = "term", all.x = T)
            go_data = go_data[order(go_data$avgP),]
            go_data$cluster = NULL
            write.table(go_data, file = "geneSet_enrichment_results_and_clusters.txt", quote=F, row.names=F, sep="\t")
        # let's try to make some wordcloud images
        if (is.data.frame(functional_clusters)){
            mostFreq_words <- CountFrequency_words(functional_clusters, n_clust)
        }
    }
