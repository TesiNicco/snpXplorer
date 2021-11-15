# snpXplorer `v2.0`
`snpXplorer` is an open source web-server to explore associations from Genome-Wide Association Studies (GWAS) available at https://snpxplorer.net/.
`snpXplorer` offers a wide range of visualization options, for example, you can superimpose multiple studies on top of each other, or on top of your study.
Please see the documentation in this repository to get more information about how to best use `snpXplorer` and its options, and how to make it work on your machine.

# Web-server structure
`snpXplorer` offers an `Exploration section` and a `Functional annotation section`. Although we recommend the usage of the web-server, `snpXplorer` can also be installed locally. While this may be more problematic for the Exploration section, as it requires several resources to be downloaded, the Functional annotation section should be easily usable. Please scroll down to section `Standalone version` for more information.

# snpXplorer
Please find the extensive snpXplorer documentation as a pdf in this folder. Just scroll all the way up.

# Stand-alone version
## Getting Started
If you are familar with R and command line, you may want to install snpXplorer Exploration or Functional Annotation sections locally on your machine. In order to download `snpXplorer`, you should clone this repository via the command below. It will take some time because the example datasets are relatively big (~400Mb). The usage of the Functional Annotation section is much easier than that of the Exploration section, as it requires far less additional files to be downloaded, and much less coding required.

```  
git clone https://github.com/TesiNicco/SNPbrowser.git
cd SNPbrowser/
```

`snpXplorer` requires R and Python v3 to be correctly installed on your machine. To work, a number of R packages need to be installed. Please refer to the documentation for more information about the required R and python packages.

## Functional annotation section
The functional annotation section required several inputs to run, most of which are chosen by the user on the web-server while completing the submission of a job. However, the annotation can also be run from the command line. A typical AnnotateMe command has the following syntax:

```
Rscript snpXplorer/AnnotateMe/BIN/MAIN.R annotateMe_input_37611.txt 3 email_address mapping default Kidney_Cortex,Whole_Blood GRCh37 > annotateMe_run_37611.log
```

where:
`snpXplorer/AnnotateMe/BIN/MAIN.R` --> is the executable
`annotateMe_input_37611.txt` --> input file containing SNPs to annotate
`3` --> [1 / 2 / 3]: input type, where 1="19:45411941" -- 2="19 45411941" -- 3="rs429358"
`email_address` --> email address
`mapping` --> [mapping / enrichment]: analysis_type, where mapping only performs SNP-gene mapping, and does not execute gene-set enrichment analysis
`default` --> [default / KEGG / WikiPathways]: analysis_mode, which defines the gene-set for gene-set enrichment analysis (in case analysis_type is enrichment)
`Kidney_Cortex,Whole_Blood` --> GTEx tissues for eQTL annotation
`GRCh37` --> the version of the reference genome to use
`annotateMe_run_37611.log` --> the log file of the run
