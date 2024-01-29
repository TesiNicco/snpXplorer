# snpXplorer
**snpXplorer** was developed to help geneticists and researchers working with genetic data to interpret association signals from Genome-Wide Association Studies (GWAS).  

Practically, **snpXplorer** is an open source web-server running at https://snpxplorer.net: the core of the **snpXplorer** is written in Python as a [Flask application](https://flask.palletsprojects.com/en/3.0.x/).  

## Where do I find more information?
**snpXplorer** was [published in Nucleic Acid Research](https://pubmed.ncbi.nlm.nih.gov/34048563/) in 2021. The paper describes the basic **snpXplorer** functionalities and the computing engine.  
While the original version of **snpXplorer** was written in R as a [Shiny application](https://shiny.posit.co/), **snpXplorer** has been completely rewritten in 2022, with many improvements in terms of speed and reproducibility.  
**snpXplorer** is still under active development.

## What are the main parts of snpXplorer?
**snpXplorer** is characterized by an *Exploration* and an *Annotation* section.
- *Exploration*: This section allows you to browse GWAS association signals across a set of pre-defined GWAS summary statistics. The user can select multiple GWAS traits and these will be superimpose on top of each other. Additional information about Single Nucleotide Polymorphisms (SNP), Structural Variants (SV), Quantitative Trait Loci (QTL) and RNA expression are also available.
- *Annotation*: This section allows to annotate a list of SNPs of interest. The user can use the provided textbox to insert SNP-IDs of interest. SNPs will be annotated to the most likely affected gene integrating public resources of QTL and variant effect predictor. The user can perform both a SNP-gene annotation (*SNP-gene mapping*) or a gene-set enrichment analysis (*Gene-set enrichment analysis*).

## What if I you have a question/comment/feedback?
Comments, questions and feedbacks are very very appreciated.  
To contact us, you can send an email to one of the following address:
- [Contact snpXplorer Team](mailto:snpxplorer@gmail.com)
- [Contact the main author](mailto:n.tesi@amsterdamumc.nl)

## Can I download snpXplorer locally?
**snpXplorer** was primarily developed to being used as a web-server.  
The source code can be downloaded locally, however, for the full capacity of both the *Exploration* and the *Annotation* sections, many additional files are required, and are not provided on GitHub. The reason for not including them on GitHub is storage: GWAS summary statistics and other genome-wide annotation tracks such as expression-QTL, splicing-QTL, and 1000Genome genotypes are very large in size.  

In case the user really wants to have a local copy of snpXplorer on her/his machine, it is still possible, and we can help the user in the setting up.  

Alternatively, if the user has a specific requests that is currently not available in **snpXplorer**, please [contact us](mailto:snpxplorer@gmail.com): we are happy to develop/implement new capabilities.  

## Exploration section
The *Exploration* section of **snpXplorer** allows the user to browse GWAS association signals and compare association statistics across multiple human traits.

### What are the GWAS I can visualize?
At the moment, the choice of the input GWAS to choose is limited to 29 GWAS summary statistics, including:
- Alzheimer's Disease ([Wightman et al.](https://www.nature.com/articles/s41588-021-00921-z), [De Rojas et al.](https://www.nature.com/articles/s41467-021-22491-8), [Kunkle et al.](https://www.nature.com/articles/s41588-019-0358-2), [Jansen et al.](https://www.nature.com/articles/s41588-018-0311-9))
- Autism ([Matoba et al.](https://www.nature.com/articles/s41398-020-00953-9))
- Depression ([Cai et al.](https://www.nature.com/articles/s41588-020-0594-5))
- Coronary Artery Disease ([Van der Harst et al.](https://www.ahajournals.org/doi/10.1161/CIRCRESAHA.117.312086), [Fall et al.](https://link.springer.com/article/10.1007/s00125-018-4686-z))
- Cardiovascular Disease ([Voijnovic et al.](https://www.nature.com/articles/s41467-018-06234-w), [Evangelou et al.](https://www.nature.com/articles/s41588-018-0205-x))
- Body Mass Index ([Yengo et al.](https://academic.oup.com/hmg/article/27/20/3641/5067845?login=false))
- Diabetes ([Forgetta et al.](https://diabetesjournals.org/diabetes/article/69/4/784/40594/Rare-Genetic-Variants-of-Large-Effect-Influence))
- COVID-19 ([Erola Pairo-Castineira et al.](https://www.nature.com/articles/s41586-020-03065-y))
- Inflammation ([Yong-Fei Wang et al.](https://www.nature.com/articles/s41467-021-21049-y), [Ruotsalainen et al.](https://www.nature.com/articles/s41431-020-00730-8), [Han et al.](https://www.nature.com/articles/s41467-020-15649-3))
- Cancer ([Zhang et al.](https://www.nature.com/articles/s41588-020-0609-2), [Bao et al.](https://www.nature.com/articles/s41586-020-2786-7), [Fiorica et al.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0236209), [Rashkin et al.](https://www.nature.com/articles/s41467-020-18246-6))
- Longevity ([Timmers et al.](https://www.nature.com/articles/s41467-020-17312-3), [Timmers et al.](https://elifesciences.org/articles/39856))
- Socio-economic ([Yengo et al.](https://academic.oup.com/hmg/article/27/20/3641/5067845?login=false), [Demange et al.](https://www.nature.com/articles/s41588-020-00754-2), [Surakka et al.](https://www.nature.com/articles/s41467-020-17315-0), [Manousaki et al.](https://linkinghub.elsevier.com/retrieve/pii/S0002929720300173))
- GWAS-Catalog  

### How do I browse the genome?
After deciding on the GWAS trait(s) to visualize, the user can set the desired version of the reference genome. Finally, by interacting with the provided textarea, the user can search for a specific SNP by using the SNP *rsid*, or a specific gene by using the *gene symbol*.  The user can adjust the window around the region of interest, which by default is 25,000 bp up/downtreat the target region.  

### What am I looking at?
After a search is done, **snpXplorer** will show 4 stacked plots:
- Main Association Plot: if you are familiar with regional plots or Manhattan plots, this will be easy for you. The x-axis shows the region of interest (in bp), the y-axis shows the association p-value (in -log10 scale). Each dot is a SNP-association, the higher the dot, the more significant is the association between the SNP and the trait of interest.
- Gene Plot: the plot shows the different genes in the region of interest. Exons can be shown optionally.
- Structural Variants Plot: the plots shows structural variants in the region of interest. The user can select the source dataset from which to take SV from.
- GTEx Expression Plot: tissue-specific gene expression of all genes shown in the region of interest.

### What options can I modify?
Additional visualization options are available, which allow the user to:
- Plot associations as single dots (one dot per SNP), or as densities (smoothing the association signal over the region) [*Plot type*]
- Plot gene exons or not (*Exons*)
- Plot recombination rates (*Recombination*)
- Define the dataset of Structural Variants (SV), from [Linthrost et al.](https://www.nature.com/articles/s41398-020-01060-5), [Audano et al.](https://www.cell.com/cell/fulltext/S0092-8674(18)31633-7?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867418316337%3Fshowall%3Dtrue) or [Chaisson et al.](https://www.nature.com/articles/s41467-018-08148-z).

### What additional outputs do I get?
As output, **snpXplorer** always reports:
- Table including the most significant SNP-associations as shown in the figure (Downloadable)
- Table including the structural variants as shown in the figure (Downloadable)
- GWAS-Catalog table, including all associations of SNPs shown in the figure as reported in the GWAS-Catalog (Downloadable)
- GTEx expression table, including the gene expression of all genes shown in the figure (Downloadable)  
Is it of course possible to download the full plot.  

## Annotation section
The *Annotation section* of **snpXplorer** allows the user to functionally annotate a set of SNPs of interest. To do this:
1. Paste a list of SNPs in the textarea. While multiple formats are accepted for the input SNPs, we recommend to use *rsid*.
2. Choose the input type
3. Define the analysis type: the *Annotation* will only annotate SNPs with the likely affected gene(s), while the *Gene-set Enrichment Analysis* will first annotate SNPs with their likely affected gene(s), and then will perform a gene-set enrichment analysis to investigate the enriched biological pathways. 
4. Select the tissues to be used for QTL-analysis.
5. Add your email
6. Click to submit  
Please note that the number input SNPs is limited to 10,000 in case of a simple *Annotation* analysis, and it is limited to 1,000 SNPs for *Gene-set Enrichment Analysis*. In case you're analysis exceeds these numbers, then the analysis will fail. However, you can still [contact us](mailto:snpxplorer@gmail.com) and we can see the best options for your analysis.

### What is happening behind the hood?
Althought a complete description of what **snpXplorer** does is explained in the [published manuscript](https://pubmed.ncbi.nlm.nih.gov/34048563/), here is a description of the analysis:
1. Input validation: the input SNPs are checked and ensured the input data is correct. If not, the user will receive an email notifying that the input SNPs were not correct.
2. Get SNP information: this procedure uses known databases to retrieve chromosome, position, allele frequecy and rsid of the SNPs of interest
3. Get CADD information: this procedure uses [CADD (Combined Annotation Dependent Depletion)](https://cadd.gs.washington.edu/) to derive the right CADD score for the SNPs of interest.
4. Get QTL information: this procedure query databased of expression-QTL and splicing-QTL for known interactions between the SNPs and genes. It is based on [GTEx](https://gtexportal.org/home/)
5. Get positional information: this procedure checks in the vicinity of the SNPs of interest to find the closest genes.
6. Aggregate information: CADD, QTL and positional information is gathered together. In case of coding SNPs, the affected gene is the one where the mutation happens. In case of non-coding SNPs, then CADD and QTL information is prioritized over the positional mapping. Finally, when no CADD and/or QTL information is available to the SNPs of interest, then the closest gene positionally is used.
7. In case the *Gene-set Enrichment Analysis* is selected, **snpXplorer** will use all genes associated with the input SNPs to find biological pathways enriched. This is done using a sampling-based framework, in order to allow each SNP to associate to multiple genes. The user can define which gene-sets to use between [Gene Ontology](https://geneontology.org/), [Wiki Pathways](https://www.wikipathways.org/), [Reactome](https://reactome.org/) and [KEGG Pathways](https://www.genome.jp/kegg/pathway.html). Gene-set Enrichment Analysis is followed by semantic similarity-based reduction of redundant pathways.

### How do I get my Annotation results?
Once you submit the *Annotation* analysis, you will receive a confirmation email. Please be aware that sometimes this email goes into spam, so you take a look there if you don't receive any email.  

A typical *Annotation* analysis, excluding gene-set enrichment analysis, will take from 5 to 30 minutes depending on the number of SNPs used as input. When *Gene-set Enrichment Analysis* is selected, then it will take longer, up to a few hours depending on the number of SNPs.  

After **snpXplorer** analysis is completed, you will receive a new email. This email contains also information how to download your results. We now implemented a [Download](https://snpxplorer.net/download/) page in **snpXplorer** which allows to download your results. You just need to paste the *run ID* received with the email, click enter, and if the *run id* is correct, a link will appear, allowing you to download the results.

### What additional outputs do I get?
When you run **snpXplorer** *Annotation* section, you will get a number of outputs:
1. SNP-gene mapping table: Table summarizing the SNP to gene mapping. The table will include, for each SNP used as input, the SNP annotations including chromosome, position, frequency, eQTL, sQTL, CADD, and the most likely affected gene(s).
2. SNP-gene mapping figure: Circular figure showing the input SNPs along the chromosomes, their allele frequency and the type of annotation (coding SNP, QTL, or based on position).
3. SV table: Table reporting the structural variants overlapping with the input SNPs.
4. GWAS-Catalog Table: Table reporting any GWAS-Catalog hits regarding the input SNPs.
5. Gene-set Enrichment Analysis (*only when Gene-set Enrichment Analysis is selected*): Table with the significantly enriched pathways.
6. Gene-set Enrichment Figures (*only when Gene-set Enrichment Analysis is selected*): Figures showing (*i*) the semantic similarity-based matrix along with dendrograms and clusters, (*ii*) the enriched terms after semantic similarity reduction, and (*iii*) a wordcloud representation for each cluster with the most common terms.

## What if you have questions?
Please do not hesitate to [contact us](mailto:snpxplorer@gmail.com) for any comment, feedback, or question.