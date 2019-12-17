# FIRST CHECK IF THE REQUIRED PACKAGES ARE INSTALLED AND IN CASE
# THEY ARE NOT, INSTALL THEM
	list.of.packages <- c("data.table", "shiny", "shinyShortcut", "stringr", "ggplot2")
	new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
	if(length(new.packages)) install.packages(new.packages)

# THEN CREATE THE EXECUTABLE FILE
	library(shinyShortcut)

	shinyShortcut(shinyDirectory = "bin/")

# FINALLY RUN THE EXECUTABLE
	cmd = "Rscript bin/.shiny_run/shinyShortcut.r"
	system(cmd)
