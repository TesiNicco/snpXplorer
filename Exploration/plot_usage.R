# Libraries
library(reshape2)
library(ggplot2)

# Create dataframe based on the usage in monitor folder
snpxplorer_data = data.frame(date = c('21/01/26', '04/02/26', '18/02/26', '05/03/26', '27/03/26', '15/04/26', '07/05/26', '27/05/26'),
    Exploration = c(272, 387, 579, 645, 740, 1085, 1344, 1373), 
    Haplotypes = c(34, 53, 69, 78, 97, 106, 133, 149),
    PRS = c(32, 38, 58, 63, 67, 79, 85, 95),
    Single_variant = c(182, 651, 690, 919, 1318, 1434, 1530, 1712),
    Batch_annotation = c(56, 65, 98, 141, 176, 240, 292, 368))

# reshape the data for ggplot
snpxplorer_data_melted = melt(snpxplorer_data, id.vars = 'date', variable.name = 'Category', value.name = 'Count')

# set date as factor
snpxplorer_data_melted$date = factor(snpxplorer_data_melted$date, levels = c('21/01/26', '04/02/26', '18/02/26', '05/03/26', '27/03/26', '15/04/26', '07/05/26', '27/05/26'))

# plot
ggplot(snpxplorer_data_melted, aes(x = Category, y = Count, fill = date)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    labs(title = 'SNPXplorer usage over time and modules', x = 'Modules', y = 'Count (# of jobs)') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 14), axis.title = element_text(size = 16), legend.title = element_blank(), legend.text = element_text(size = 14), legend.position = 'top', title = element_text(size = 18))

