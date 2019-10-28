# SNPbrowser `v1.0`
`SNPbrowser` is a command line tool written in R to display results from Genome-Wide Association Studies (GWAS) with a wide range of possibilities, including overlap of multiple studies. Once executed from the command line, the tool opens and runs on a web application with easy-to-use user interface.

`SNPbrowser` is based on shiny library in R.

`SNPbrowser` is undergoing active development and many additional features will be implemented in future versions of the application.

## Getting Started
In order to download `SNPbrowser`, you should clone this repository via the commands

```  
git clone https://github.com/TesiNicco/SNPbrowser.git
cd SNPbrowser/
```


`SNPbrowser` requires R to be correctly installed on your machine. To work, the following R libraries also need to be installed: 
`shiny`, `data.table`, `stringr`, `ggplot2`.


If you don't have these packages installed, or you don't know, you can run the following in R:

```
if(!require(c("shiny", "data.table", "stringr", "ggplot2")))install.packages(c("shiny", "data.table", "stringr", "ggplot2"))
```

This will check for the presence of the package and in case it is not present, it will install it.


Once the above has completed, you can exit R and run the following to open a new window on your browser and start a `SNPbrowser` instance: 

```
Rscript -e 'library(methods); shiny::runApp("bin/", launch.browser=TRUE)'
```

`SNPbrowser` can run on it own at this point. In its most basic form (the one you will download), there are two main folders (bin and data). `bin/` contains the main scripts to execute `SNPbrowser`: `server.R` and `ui.R`. In addition, we provide the `python` script `GetNewGWAS.py` to download additional GWAS summary statistics, split them by chromosomes and make them suitable to be loaded by `SNPbrowser`. The second folder is `data/`, that contains the `databases/` and `example/` folders. `databases/` folder contains annotation files of the genes and recombination rates in human genome build hg19. The `example/` folder contains example data to try out the tool. These data consists of the summary statistics (chromosome 16 to 21) of the publicly available case-control meta-analysis of Alzheimer's disease as published in Kunkle et al., 2019.


It is possible to use your own data and additional data (for example, publicly available datasets), in two ways.
1. The first option to use your data or piblicly available data (most recommended for us) is to use the `GetNewGWAS.py` script. 
2. The second option to use your data or publicly available data is to load them through the `Load` button. In this way, every time you need to visualize data, you need to load files again.

## Adding new GWAS files
To add new GWAS datasets and make them always available within `SNPbrowser`, you can use the `python` script `GetNewGWAS.py`. This script adds as repositories new or existing GWAS dataset so that `SNPbrowser` knows where to find them. For this script to work, make sure you have `python` correctly installed in your system. In addition, make sure you have installed these basic `python` libraries: `os`, `sys`, `gzip`. To understand `GetNewGWAS.py` arguments, enter `SNPbrowser` main folder and type `python bin/GetNewGWAS.py --help`. The followind help message will be displayed:

```
##### ADD NEW GWAS #####
### ARGUMENTS:
### arg1: MODE (--download-split-parse=download+split+add; --split-add=split+add; --add=add; --slim-add=slim+add; --help=this_message)
### Supplementary arguments:
### 	MODE --download-split-add: arg2 --> URL (The url of the GWAS to download)
### 	MODE --download-split-add: arg3 --> TIT (This is the title for the GWAS and folders)
#############################################################################
### 	MODE --split-add: arg2 --> GWAS (This is the name of the GWAS dataset)
### 	MODE --split-add: arg3 --> TIT (This is the name of the folder and GWAS)
#############################################################################
### 	MODE --add: arg2 --> TIT (This is the name of the folder and GWAS)
#############################################################################
### 	MODE --slim-add: arg2 --> GWAS (This is the name of the chromosome 1 of the GWAS dataset)
### 	MODE --slim-add: arg3 --> TIT (This is the name of the folder and GWAS)
#############################################################################
```

There are 4 ways of running the script:
1. `--download-split-add`: select this to download a GWAS dataset, split each chromosome (this greatly speeds up `SNPbrowser`) and add the GWAS as repository. Additional mandatory arguments for this mode are the url of the GWAS to download and the title (a new folder in `data/` will be created with the given name and this name will be used to access the relative GWAS in `SNPbrowser`). For example, typing:
```
python bin/GetNewGWAS.py --download-split-add https://blablablaGWAS.txt BLA
```
will download the corresponding GWAS dataset, create a folder `data/BLA` and will put inside each chromosome of the download GWAS dataset. We provide a list of GWAS URLs that can be used.

2. `--split-add`: select this to split each chromosome (this greatly speeds up `SNPbrowser`) and add the GWAS as repository. Additional mandatory arguments for this mode are the url of the GWAS to split and the title (a new folder in `data/` will be created with the given name and this name will be used to access the relative GWAS in `SNPbrowser`). For example, typing:
```
python bin/GetNewGWAS.py --split-add BLA_GWAS.txt BLA
```
will split the corresponding GWAS dataset, create a folder `data/BLA` and will put inside each chromosome of the download GWAS dataset.

3. `--add`: select this to add the GWAS as repository. We recommend to use this option if you already have you GWAS dataset splitted in single chromosomes and the format is three columns: `chromosome`, `position` and `p-value`. Additional mandatory arguments for this mode is the title (this assumes the existance of a folder in `data/` with the given name (title) and this name will be used to access the relative GWAS in `SNPbrowser`). For example, typing:
```
python bin/GetNewGWAS.py --add BLA
```
will add as repository for `SNPbrowser` the BLA GWAS dataset that underlies summary statistics data in `data/BLA/`.

3. `--slim-add`: select this to clean GWAS data and add them as repository. We recommend to use this option if you already have you GWAS dataset splitted in single chromosomes and the format is NOT three columns: `chromosome`, `position` and `p-value`. To speed up the tool, this mode can be used to take only informative columns. Additional mandatory arguments for this mode is the chromsome 1 of the GWAS dataset (the other chromosome-files will be determined automatically) and this name will be used to access the relative GWAS in `SNPbrowser`). For example, typing:
```
python bin/GetNewGWAS.py --slim-add trial_data/chr1_GWAS.txt BLA
```
will add create a folder `data/BLA/`, put cleaned files in there, and add as repository for `SNPbrowser` the BLA GWAS dataset .
