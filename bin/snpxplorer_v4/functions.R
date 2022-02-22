# function to plot points and recombination rates
plot_points <- function(snps_data, recomb_data, significance, pos_start, pos_end){
  # initialize the plot
  fig <- plot_ly()
  # add data from first dataframe Df1
  fig <- fig %>% add_markers(data=snps_data, name="SNPs", x = ~pos, y = ~-log10(p))
  # add data from second dataframe Df2
  fig <- fig %>% add_trace(data=recomb_data, name="Recombination", x = ~pos, y = ~combined, type = 'scatter', mode = 'lines', yaxis = 'y2')
  # show figure
  fig <- fig %>% layout(plot_bgcolor='#e5ecf6', title = "Regional plots", 
                        yaxis2 = list(zerolinecolor = '#ffff', gridcolor = 'ffff', zerolinewidth = 4, title = "Recombination rates", autorange = FALSE, range = c(0, 100), overlaying = "y", side = "right"),
                        yaxis = list(zerolinecolor = '#ffff', gridcolor = 'ffff', zerolinewidth = 4, title = "-log10(P-value)", autorange = FALSE, range = c(0, significance)), 
                        xaxis = list(zerolinecolor = '#ffff', gridcolor = 'ffff', zerolinewidth = 4, title = "Genomic position (bp)", autorange = FALSE, range = c(pos_start, pos_end)))
  return(fig)
}

# function to calculate density
densityLinePvalue <- function(snps_data, wind.n, span_value){
  min.wind <- min(snps_data$pos)          # define minimum position
  max.wind <- max(snps_data$pos)          # define maximum position
  interval <- ceiling((max.wind - min.wind)/wind.n)   # define windows
  out = matrix(data = NA, nrow = wind.n, ncol = 3)    # define output
  colnames(out) <- c("window", "pvalue", "error")
  #counter for output assignement
  counter <- 1
  for (i in seq(min.wind, max.wind, interval)){
    if (counter <= wind.n){
      #define internal maximum -- f
      f <- i + interval
      #take info in the window of interest
      data.sbs <- snps_data[which((snps_data$pos >= i) & (snps_data$pos <= f)),]
      data.sbs$p <- as.numeric(data.sbs$p)
      data.sbs <- data.sbs[order(data.sbs$p),]
      #check number of rows
      if (nrow(data.sbs) > 0){
        #take maximum pvalue point
        out[counter, ] <- c(((f - i)/2) + i, -log10(data.sbs$p[1]), sd(na.omit(data.sbs$p)))
      } else {
        out[counter, ] <- c(((f - i)/2) + i, 0, 0)
      }
      counter <- counter + 1
    }
  }
  out[which(is.na(out[, "error"]) == TRUE), "error"] <- 0
  out <- as.data.frame(out)
  lo <- loess(out$pvalue ~ out$window, span=span_value)
  xl <- seq(min(out$window), max(out$window), (max(out$window) - min(out$window))/100)
  pred <- predict(lo, xl)
  pred[which(pred < 0)] <- 0
  #add limits -- left and right for polygon function
  xl <- c(xl[1], xl, xl[length(xl)])
  pred <- c(0, pred, 0)
  df = data.frame(x = xl, y = pred)
  return(df)
}

# function to plot densities and recombination rates
plot_density <- function(snps_data, recomb_data, significance, pos_start, pos_end, wind_n, span_value){
  # first get the densities
  df = densityLinePvalue(snps_data, wind_n, span_value)
  # initialize the plot
  fig <- plot_ly()
  # add data from first dataframe Df1
  fig <- fig %>% add_trace(data = df, name="SNPs", x = ~x, y = ~y, type = 'scatter', mode = 'lines', fill = 'tozeroy')
  # add data from second dataframe Df2
  fig <- fig %>% add_trace(data=recomb_data, name="Recombination", x = ~pos, y = ~combined, type = 'scatter', mode = 'lines', yaxis = 'y2')
  # show figure
  fig <- fig %>% layout(plot_bgcolor='#e5ecf6', title = "Regional plots", 
                        yaxis2 = list(zerolinecolor = '#ffff', gridcolor = 'ffff', zerolinewidth = 4, title = "Recombination rates", autorange = FALSE, range = c(0, 100), overlaying = "y", side = "right"),
                        yaxis = list(zerolinecolor = '#ffff', gridcolor = 'ffff', zerolinewidth = 4, title = "-log10(P-value)", autorange = FALSE, range = c(0, significance)), 
                        xaxis = list(zerolinecolor = '#ffff', gridcolor = 'ffff', zerolinewidth = 4, title = "Genomic position (bp)", autorange = FALSE, range = c(pos_start, pos_end)))
  return(fig)
}

# function to plot genes
plot_genes <- function(snps_data){}