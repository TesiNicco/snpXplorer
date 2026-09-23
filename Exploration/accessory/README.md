# SNPXplorer accessory scripts

This directory contains offline preprocessing, analysis, and reporting scripts associated with SNPXplorer. These files are not imported by the current web application during normal operation. Some of them generate data products that can subsequently be copied into `Data/` and used by the application.

Several legacy scripts use constants and hard-coded paths rather than command-line arguments. Consequently, paths and analysis settings should be reviewed before each run. Python scripts should be run from the repository root unless their paths are updated. `haplotype_region_ld.py` instead has a command-line interface and repository-relative defaults.

## Contents

| Script | Purpose | Runtime component? |
| --- | --- | --- |
| `R2cluster.py` | Identifies proxy haplotypes by clustering variants according to pairwise LD within recombination-defined intervals. | No; offline data generation |
| `haplotype_region_ld.py` | Extracts the SQLite-stored LD between every variant in a supplied set of haplotypes for regional validation and plotting. | No; offline analysis |
| `sumstats_allele_match.py` | Harmonizes the allele orientation of repeated GWAS associations at each genomic locus. | No; offline database preparation |
| `semantic_similarity_LLM.py` | Groups semantically similar GWAS traits and creates trait embeddings and a UMAP projection. | No; offline model/data generation |
| `plot_usage.R` | Draws a bar chart from manually entered cumulative SNPXplorer usage counts. | No; offline reporting |
| `app_bck.py` | Legacy backup of an earlier SNPXplorer Flask application. | No; retained as a historical snapshot |

`.DS_Store` is macOS metadata and is not an analysis script.

## `R2cluster.py`

### Brief description

Defines candidate genomic intervals using recombination-rate breakpoints, calculates an unphased pairwise LD matrix with PLINK 2 for every interval, and applies complete-linkage hierarchical clustering to that matrix. Each retained cluster is used as a proxy for a haplotype: a set of variant alleles expected to be inherited together because of LD.

The script runs the intervals in parallel and retains clusters containing at least five variants in its current top-level workflow.

### Inputs

- Recombination map: `recomb_rates_path`, currently `/project/holstegelab/Share/nicco/reference_files/recombination_rates_hg38.txt.gz`. It must contain `Chromosome`, `Position(bp)`, and `Rate(cM/Mb)` columns.
- Genotype data: `genotype_data_path`, currently a directory of chromosome-specific PLINK 2 datasets. For chromosome `chrN`, the expected prefix is `chrN.dose.unscrambled`, with the associated `.pgen`, `.pvar`, and `.psam` files.
- Sample inclusion list: currently hard-coded inside `calculateLD()` as `/project/holstegelab/Share/nicco/common_SV/results/samples_of_interest.txt`.
- Variant identifiers in the PLINK data must follow `CHR:POS:REF:ALT`, because positions are parsed from this format.
- External executable: `plink2` must be available on `PATH`.

### Outputs

- Main TSV: `oupath`, currently `/project/holstegelab/Share/nicco/haplotypeBlocks/R2cluster/sensitivity_haplotyping/20260910_thr0.5_haplotypes.tsv`.
- Output columns:
  - `CHR`: chromosome.
  - `BP1`, `BP2`: minimum and maximum variant positions in the cluster.
  - `KB`: cluster span in kilobases.
  - `NSNPS`: number of variants in the cluster.
  - `SNPS`: position-sorted variant identifiers separated by `|`.
- Temporary PLINK LD matrix, variant-list, and log files are created in the working directory and removed after each successful interval.

The commented final lines show an optional bgzip/Tabix workflow, but they are not executed.

### Parameters

- `recomb_intervals(..., threshold=25.0)`: starts a new interval wherever the recombination rate is greater than this value, in cM/Mb. The top-level call uses `25`.
- `rate_col`, `chrom_col`, `pos_col`: recombination-map column names; defaults are `Rate(cM/Mb)`, `Chromosome`, and `Position(bp)`.
- `thresh = 0.2`: minimum pairwise LD value implied by complete-linkage clustering. Clustering uses distance `1 - LD` and cuts at `1 - thresh`.
- `min_snps`: minimum number of variants retained per cluster. The function default is `2`, while `process_interval()` explicitly uses `5`.
- PLINK options inside `calculateLD()`:
  - `--r2-unphased square`: calculate a square unphased LD matrix.
  - `--keep`: restrict samples to the hard-coded sample list.
  - `--maf 0.005`: exclude variants with minor-allele frequency below 0.005.
- `n_procs = 6`: number of multiprocessing workers.
- `chunk_size = 6`: intervals supplied to each worker per multiprocessing chunk.

### Important notes

- The output filename contains `thr0.5`, but the effective clustering threshold in the current code is `thresh = 0.2`. The filename should not be used to infer the analysis setting.
- The script has no `if __name__ == "__main__"` guard. Importing it executes the complete analysis.
- `process_interval()` catches all exceptions and returns no result for a failed interval. Therefore, failures are currently silent.

## `haplotype_region_ld.py`

### Brief description

Creates a single, physical-position-ordered LD matrix containing the union of all variants assigned to one or more specified SNPXplorer haplotypes. Haplotype labels such as `A` and `B` are supplied explicitly because those display labels are assigned dynamically by the web application and are not stored in the haplotype table.

LD is read directly from the chromosome-specific SNPXplorer SQLite database. The script retrieves every stored LD edge between the selected variants and does not apply an additional r² threshold. It does not run PLINK or recalculate LD from genotypes. Because an edge absent from SQLite may have been omitted when that database was created, a separate observed-pair mask is always written.

### Inputs

- Haplotype selection, supplied in exactly one of two ways:
  - Repeat `--haplotype LABEL=CHR:START-END`, for example `--haplotype B=7:12194394-12246599`.
  - Use `--manifest` with a tab-separated file containing the columns `label` and `haplotype_id`.
- `--haplotype-table`: bgzip-compressed, Tabix-indexed SNPXplorer haplotype table. The default is `Data/databases/haplotypes/ld_clusters_with_gwas_AI_latest.tsv.gz`. Its first six columns must be `CHR`, `BP1`, `BP2`, `KB`, `NSNPS`, and pipe-separated `SNPS` identifiers in `CHR:POS:REF:ALT` form.
- `--ld-directory`: directory containing chromosome-specific databases named `ld_chr<CHR>.sqlite`. The default is `Data/databases/LD_db`.
- All selected haplotypes must be located on the same chromosome. Overlapping haplotypes are supported, and a variant belonging to multiple haplotypes appears only once in the matrix.

### Outputs

For an `--output-prefix` of `results/tmem106b`, the script writes:

- `results/tmem106b.ld_matrix.tsv.gz`: symmetric r² matrix with variants ordered by genomic position.
- `results/tmem106b.ld_observed_mask.tsv.gz`: symmetric binary matrix identifying pairs actually present in SQLite. The diagonal is one.
- `results/tmem106b.ld_observed_pairs.tsv.gz`: long-form table of stored variant pairs, r² values, and physical distances.
- `results/tmem106b.ld_heatmap.png`: borderless, lower-triangular SNP-by-SNP heatmap in physical-position order, without individual SNP tick labels.
- `results/tmem106b.ld_heatmap_by_haplotype.png`: borderless, lower-triangular SNP-by-SNP heatmap grouped by natural haplotype-label order (`A`, `B`, and so on), with position ordering within each haplotype. Aligned colored annotation strips and enlarged labels identify haplotype blocks to the left and below the matrix.
- `results/tmem106b.ld_heatmap_genomic_position.png`: position-constrained roof-shaped LD plot. The matrix is coordinate-transformed so its diagonal lies on a true genomic-position X-axis. Cells receive a small minimum display width for visibility. Pairs absent from SQLite share the navy “NA and 0–0.2” display category with observed low-LD pairs, matching the conventional heatmap.
- `results/tmem106b.variants.tsv`: physical matrix order, haplotype-grouped plot order, genomic coordinates, complete and primary plot haplotype memberships, SQLite identifiers, and mapping status for every variant.
- `results/tmem106b.haplotypes.tsv`: requested labels and IDs plus declared and retrieved variant counts.
- `results/tmem106b.qc.json`: source paths, mapping counts, overlapping-membership counts, possible/observed pair counts, missing-edge policy, and output paths.

### Parameters

- `--haplotype LABEL=CHR:START-END`: labeled haplotype selection; may be repeated.
- `--manifest FILE`: alternative multi-haplotype input. It is mutually exclusive with `--haplotype`.
- `--haplotype-table FILE`: override the default haplotype-membership table.
- `--ld-directory DIRECTORY`: override the default SQLite directory.
- `--output-prefix PATH`: basename and optional destination directory for all outputs; default `./haplotype_region_ld`.
- `--missing-value {zero,nan}`: represent pairs absent from SQLite as `0` or `NA` in the LD matrix; default `zero`, matching the current web-app convention. The observed-mask file distinguishes absent edges in either mode.
- `--heatmap-dpi INTEGER`: PNG resolution; default `300` and minimum `72`.
- `--verbose`: enable debug-level logging.

Both heatmaps place a compact horizontal LD legend above the matrix and use the same discrete colors as SNPXplorer: navy (≤0.2), light blue (>0.2–0.4), green (>0.4–0.6), orange (>0.6–0.8), and red (>0.8). If a variant belongs to more than one supplied haplotype, it is plotted once and assigned to the earliest natural-sorted label in the grouped heatmap; all memberships remain recorded in `variants.tsv`.

Example for haplotypes A and B:

```bash
python Exploration/accessory/haplotype_region_ld.py \
  --haplotype A=7:START-END \
  --haplotype B=7:12194394-12246599 \
  --output-prefix results/tmem106b_all_haplotypes
```

The literal `START-END` in this example must be replaced by haplotype A's stored coordinates.

## `sumstats_allele_match.py`

### Brief description

Harmonizes repeated association records at each chromosome-position locus in two passes:

1. Count each unordered allele pair at the locus and retain the majority pair. The first encountered ordered orientation of that winning pair becomes the retained `EA`/`NEA` orientation.
2. Write records with the retained pair, leaving matching orientations unchanged and flipping reversed orientations. Records carrying a minority allele pair are removed.

When a record is flipped, `EA` and `NEA` are exchanged, `beta` is negated, and a valid `EAF` is replaced by `1 - EAF`. EAF concordance is logged for auditing but does not decide whether the record is flipped.

### Inputs

- `INPUT_FILE`: `Data/databases/haplotypes/All_indep_gwas_sumstats_AI_hg38.txt.gz`.
- The input must be a headered, gzip-compressed, tab-separated file.
- The current zero-based column mapping is:
  - chromosome: `2`
  - position: `3`
  - effect allele (`EA`): `5`
  - non-effect allele (`NEA`): `6`
  - effect-allele frequency (`EAF`): `7`
  - effect estimate (`beta`): `8`

### Outputs

- `OUTPUT_FILE`: `Data/databases/haplotypes/All_indep_gwas_sumstats_AI_hg38_allele_matches.txt.gz`.
- The output retains the original header and column order.
- Minority-pair records are absent; reversed records contain harmonized alleles, EAF, and beta.
- A temporary `<OUTPUT_FILE>.tmp` is written first and atomically moved to the final path only after a successful pass.
- The logger reports input/output counts, unchanged and flipped records, removed records, affected loci, allele-pair ties, and EAF audit categories.

### Parameters

- `INPUT_FILE`, `OUTPUT_FILE`: input and output paths.
- `PROGRESS_INTERVAL = 1_000_000`: number of records between progress messages.
- `CHROM_INDEX`, `POSITION_INDEX`, `EA_INDEX`, `NEA_INDEX`, `EAF_INDEX`, `BETA_INDEX`: zero-based schema positions. These must be changed if the input layout changes.
- Majority rule: the unordered allele pair with the most records at a locus is retained.
- Tie rule: when pairs have equal maximum counts, the pair encountered first is retained and the tie is logged.
- Orientation rule: the first ordered representation of the winning pair is retained.
- Valid EAF range: inclusive interval `[0, 1]`; missing, nonnumeric, or out-of-range values are treated as unavailable.

Run from the repository root with:

```bash
python Exploration/accessory/sumstats_allele_match.py
```

The generated file is not selected automatically by SNPXplorer merely because it exists. It must replace, or be wired into the application in place of, the active association database.

## `semantic_similarity_LLM.py`

### Brief description

Builds semantic representations of GWAS trait names using two biomedical sentence-transformer models, clusters the traits at several cosine-distance thresholds, chooses representative studies, creates clustering diagnostics, and generates a two-dimensional UMAP projection for display in SNPXplorer.

Trait strings are normalized by lowercasing, removing selected punctuation and GWAS boilerplate phrases, dropping a small set of stopwords, and collapsing whitespace. Each model produces an independently normalized embedding. The embeddings are scaled by the square root of their model weights, concatenated, and normalized again, making cosine similarity equal to the weighted mean of the two models' cosine similarities.

### Inputs

- `data_path`: currently `/project/holstegelab/Share/nicco/haplotypeBlocks/R2cluster/OpenGWAS_sumstats`.
- Every `.gz` file in that directory is considered except files beginning with `eqtl`. The script obtains the trait name from the second field of the first data row using `zcat ... | head -n 2`.
- Association-count table: `20251122_newhaplo/trait_association_counts.txt`, a tab-separated file with at least `Trait Name`, `Association Count`, and `File Name`.
- Transformer models downloaded/loaded by `sentence-transformers`:
  - `FremyCompany/BioLORD-2023`
  - `amrothemich/sapbert-sentence-transformers`
- External commands: `zcat` and `head` must be available.

### Outputs

- `<data_path>/trait_names.txt`: extracted trait names, one per input file.
- `20260731_trait_clusters.pkl`: lists of representatives, cluster maps, labels, raw representative indices, ensemble embeddings, and distance matrices for every tested threshold.
- `20260731_trait_clusters_summary_AI_thresholds.csv`: a tab-separated intermediate cluster summary despite its `.csv` extension.
- `20260731_trait_clusters_summary_AI_thresholds.png`: number and size distribution of clusters across thresholds.
- `20260714_trait_clusters_summary_AI_thresholds.csv`: comma-separated final summary augmented with the representative study having the largest association count.
- `20260802_umap_trait_embeddings.png`: UMAP scatter plot.
- `20260802_trait_embeddings_umap_thr05.pkl`: UMAP coordinates, embeddings, labels, and cluster map intended for the threshold-0.5 display.
- Console diagnostics: five closest distinct traits and five closest clusters to `query_trait`.

Most outputs are written relative to the current working directory. The dated files are not automatically copied or renamed to the filenames read by SNPXplorer under `Data/`.

### Parameters

- `STOPPHRASES`: regular expressions removed from trait names, including phrases such as “GWAS of” and “meta-analysis”.
- `STOPWORDS`: words removed from normalized trait strings.
- `MODEL_NAMES`: sentence-transformer model identifiers.
- `WEIGHTS = [0.7, 0.3]`: relative BioLORD and SapBERT contributions; normalized internally.
- `encode_ensemble(..., batch_size=256, show_progress_bar=True)`: encoding batch size and progress display.
- `cluster_embeddings(..., distance_threshold=0.20)`: agglomerative-clustering cut on cosine distance. Lower values produce more, tighter clusters.
- Tested thresholds: `0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7`.
- Clustering linkage: `average`, operating on a precomputed cosine-distance matrix.
- `reassign_singletons(..., min_similarity=0.40)`: assigns a singleton to its nearest original non-singleton cluster only when cosine similarity is at least 0.40.
- `pick_representatives(..., method='centroid')`: initially chooses the observed trait nearest its cluster centroid. The alternative `shortest` method is implemented but not used by the top-level workflow.
- Final representative rule: among the traits in a semantic cluster, choose the study with the greatest `Association Count`.
- UMAP parameters: `n_neighbors=50`, `min_dist=0.1`, `metric='cosine'`, and `random_state=42`. The embedding array is selected from the threshold-0.5 run.
- `query_trait = "Multiple sclerosis"`: trait used only for printed nearest-neighbor diagnostics.

### Important notes

- The script has no main guard. Importing it executes trait extraction, model loading, all clustering runs, plotting, and serialization.
- The two transformer models may be downloaded on first use and require substantial memory and execution time.
- In the current code, `ens` is explicitly reset to the threshold-0.5 embedding before UMAP, but `labels` and `cluster_map` are not reset after the loop over thresholds. They therefore retain the last loop's threshold-0.7 values when plotting and writing `20260802_trait_embeddings_umap_thr05.pkl`. This should be corrected before regenerating the production UMAP artifact.

## `plot_usage.R`

### Brief description

Creates a grouped bar chart showing cumulative SNPXplorer job counts over time for the Exploration, Haplotypes, PRS, Single-variant, and Batch-annotation modules.

### Inputs

- No external input file is read.
- Dates and cumulative counts are entered directly in the `snpxplorer_data` data frame.
- R packages: `reshape2` and `ggplot2`.

### Outputs

- A `ggplot2` bar chart is printed to the active R graphics device.
- The script does not call `ggsave()` and therefore does not create an image file by itself.

### Parameters

- `date`: measurement dates, stored as `DD/MM/YY` strings.
- `Exploration`, `Haplotypes`, `PRS`, `Single_variant`, `Batch_annotation`: cumulative job counts for each module at each date.
- Factor `levels`: explicit chronological display order for the dates and legend.
- `geom_bar(stat='identity', position='dodge')`: uses the supplied counts and places dates side by side.
- `labs(...)` and `theme(...)`: plot title, axis labels, font sizes, label angle, and legend position.

Run with:

```bash
Rscript Exploration/accessory/plot_usage.R
```

For non-interactive use, add an explicit `ggsave()` call or open a graphics device before constructing the plot.

## `app_bck.py`

### Brief description

A historical backup of an earlier, monolithic SNPXplorer Flask application. It includes user registration and login, locus/GWAS search endpoints, regional GWAS exploration, annotation submission, download routes, LD and haplotype lookup, GTEx expression, GWAS Catalog and structural-variant tracks, and Matplotlib/Plotly plotting helpers.

This file is not imported by the current application and should not be treated as the authoritative implementation. It is documented here so its role is unambiguous and its dependencies can be understood if the snapshot must be inspected or reproduced.

### Inputs

- `data_path`: currently an absolute local path ending in `snpXplorer/Data`; a commented server alternative is `/Data`.
- Local databases beneath `data_path`, including gene tables and `variants_info.sqlite`, GTEx expression, LD, haplotypes, structural variants, recombination maps, and GWAS Catalog resources.
- IEU OpenGWAS metadata and association results through `ieugwaspy`.
- HTTP form/query inputs from the SNPXplorer templates.
- Redis at `redis://localhost:6379` for Flask sessions.
- SQLite user database `sqlite:///db.sqlite`.
- External utilities used by helpers include `tabix`, `zgrep`, and `Rscript`.
- The annotation route expects the external pipeline `/Annotation/BIN/MAIN.R` and writable `/Annotation/RUNS/` paths.

### Outputs

- Rendered Flask pages and JSON responses.
- Server-sent progress events at `/events`.
- Interactive Plotly regional and GTEx plots, plus downloadable PNG versions.
- Downloadable SNP, structural-variant, GWAS Catalog, GTEx, and annotation-result files.
- SQLite user records and Redis-backed session state.
- Temporary annotation inputs and logs under `/Annotation/RUNS/`.

### Parameters

The file has no command-line interface. Its principal configuration and request parameters are:

- `data_path`: root directory for reference databases.
- `SQLALCHEMY_DATABASE_URI`: user database location; currently `sqlite:///db.sqlite`.
- `SESSION_TYPE`, `SESSION_REDIS`, `SESSION_USE_SIGNER`, `SECRET_KEY`: Flask session configuration. The literal backup value `SECRET_KEY='secret'` is unsuitable for deployment.
- `browse`: locus expressed as a gene, rsID, `chromosome:position`, or `chromosome:start-end`.
- `refGenome`: reference build, expected as `GRCh37` or `GRCh38`.
- `window`: number of base pairs added around the selected locus; default `25,000`.
- `gwas_data`: one or more selected OpenGWAS studies.
- `recomb`: include recombination-rate data when set to `Yes`.
- `exons`: include exon rendering when set to `Yes`.
- `plotype`: regional plot type; the default is `Scatter`.
- `ld`: calculate/display LD when set to `Yes`.
- `plotHaplo`: include haplotype blocks when set to `Yes`.
- `sv_source`: selected structural-variant categories. `tr` expands to Tandem Repeat and Satellite; `sine`, `line`, and `ltr` map to their corresponding repeat classes.
- `y_range`: `Auto` or `Manual`; manual mode uses integer `y_min` and `y_max`.
- `QTLtissuesExplo`: selected QTL tissues.
- GWAS/locus search API parameters: `q` is the search string, `page` is the one-based page number, and the fixed page size is 50.
- `extract_ld(..., r2_min)` helpers: minimum LD threshold; partner-query helpers default to `0.0` and accept an optional result `limit`.
- Annotation inputs: variant list, delimiter/input type, reference genome, mapping versus enrichment analysis, enrichment sources, QTL tissues, and email address.
- Flask development server: `app.run()` uses Flask defaults because no host, port, or debug values are supplied.

To start this historical snapshot directly:

```bash
python Exploration/accessory/app_bck.py
```

This is not recommended for production: it initializes databases and OpenGWAS metadata at import time, relies on old templates and infrastructure, contains environment-specific absolute paths, and uses development configuration.
