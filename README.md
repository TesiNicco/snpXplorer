# SNPbrowser `v1.0`
`SNPbrowser` is a command line tool written in R to display results from Genome-Wide Association Studies (GWAS) with a wide range of possibilities, including overlap of multiple studies. Once executed from the command line, the tool opens and runs on a web application with easy-to-use user interface.

`SNPbrowser` is based on shiny library in R.

## Getting Started
In order to download `SNPbrowser`, you should clone this repository via the commands

```  
git clone https://github.com/TesiNicco/SNPbrowser.git
cd SNPbrowser/
```


`SNPbrowser` requires R to be correctly installed on your machine. To work, few R libraries need to be installed: 
`shiny`
`data.table`
`stringr`
`ggplot2`


If you don't have these packages installed, or you don't know, you can run the following within R:

```
if(!require("package_name"))install.packages("package_name")
```

This will check for the presence of the package and in case it is not present, it will install it.


Once the above has completed, you can exit R and move inside `SNPbrowser` folder. To execute, you can run: 

```
Rscript -e 'library(methods); shiny::runApp("bin/", launch.browser=TRUE)'
```
