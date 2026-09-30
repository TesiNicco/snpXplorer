# Haplotype block estimation & semantic trait reduction

Offline pipeline that builds the **haplotype**, **GWAS association**, **LD** and
**PRS** resources used by snpXplorer, for two ancestries (European and African).

## Idea

Public GWAS resources (here: [OpenGWAS](https://gwas.mrcieu.ac.uk/)) contain
thousands of studies, many of which measure the same or almost the same trait
("LDL cholesterol", "LDL-C levels", "Low density lipoprotein cholesterol", ...).
Counting raw associations per genomic region is therefore biased toward traits
that were studied many times. The pipeline addresses this in three parts:

1. **Haplotype blocks from LD.** The genome is split into low-recombination
   intervals using a recombination map. Within each interval, variants are
   grouped by hierarchical clustering on pairwise LD (r²). The result is a set of
   "haplotypes", which are groups of variants in mutual LD.
2. **Semantic de-duplication of traits.** Trait names are embedded with two
   biomedical sentence-embedding models (BioLORD-2023 and SapBERT) and grouped
   by hierarchical clustering on cosine distance. For each group of synonymous
   or near-synonymous traits, only one study is kept: the one with the most
   associations.
3. **Annotation and resources.** Associations from the retained,
   non-redundant studies are lifted to hg38 and intersected with the haplotype
   blocks. This gives, for each haplotype, how many independent traits it is
   associated with. The same GWASs are also clumped/conditioned (GWASum) to
   build PRS-ready variant sets, and ancestry-specific LD databases (SQLite)
   are built for the web app.

```
                 ┌───────────────────────────┐
                 │ 01_gather_associations.py │  OpenGWAS top hits (p < 5e-5)
                 └─────────────┬─────────────┘
          ┌────────────────────┼─────────────────────────┐
          ▼                    ▼                         ▼
 02_clump_PRS.py     04_semantic_reduction_AI.py    03_R2cluster_{EUR,AFR}.py
 (GWASum: clump +    (embed + cluster trait names;  (recombination intervals
  COJO → PRS sets)    pick 1 study per cluster)      + LD clustering → haplotypes)
          │                    │         │                   │
          ▼                    │         ▼                   ├──► 06_plot_haplotype_statistics.py
 09_compress_and_check_PRS.py  │  07_align_semantic_         │
                               │     clustering.py (UMAP)    │
                               ▼                             │
                 05_count_associations_haplotypes.py ◄───────┘
                 (liftover hg19→hg38, count hits per haplotype)
                               │
                               ▼
                 sumstats_allele_match.py (harmonise alleles across studies)

 08_create_LD_db.py      → per-chromosome SQLite LD databases (independent of the above)
 haplotype_region_ld.py  → CLI: LD matrix + heatmaps for chosen haplotypes (uses 05 + 08 outputs)
```

## General notes

- **These scripts are not command-line tools.** Apart from
  `haplotype_region_ld.py`, each script is written to be run top to bottom
  (or interactively, block by block) after you edit the configuration variables
  at the top. The most important one is `ancestry` / `population_of_interest`,
  which is set to `"African"` in most scripts at the moment.
- **Paths.** Most paths are absolute paths on the cluster
  (`/project/holstegelab/Share/...`) or relative to the working directory
  `/project/holstegelab/Share/nicco/haplotypeBlocks/R2cluster/`. Output file
  names include hard-coded date prefixes (for example `20260917_`). When you
  re-run a step, update both the output name and the input names of the
  downstream scripts.
- **Genome builds.** OpenGWAS top hits are **hg19**. The genotype data, the
  recombination map, the haplotypes and the LD databases are **hg38**. Liftover
  happens in `02` (inside GWASum, via `--liftover`) and in `05` (pyliftover).
- **Genotype reference panels.**
  - European: TOPMed-imputed array data (`Imputed_unscRambled_genotypes`,
    plink2 pfiles `chr{N}.dose.unscrambled`), restricted to a sample subset.
  - African: African samples from 1000 Genomes (`chr{N}_african`).
- **Software.** Python with `pandas`, `numpy`, `scipy`, `tqdm`, `ieugwaspy`,
  `pyliftover`, `liftover`, `sentence-transformers`, `scikit-learn`,
  `umap-learn` and `matplotlib`. External tools: `plink2`, `bgzip`, `tabix`,
  and [GWASum](/project/holstegelab/Share/nicco/GWASum/GWASum.py). The NLP
  scripts (`04`, `07`, `09`) were run in the `nlpy` conda environment. `05`
  calls `bgzip`/`tabix` from `~/.conda/envs/newtreat/bin/`.

---

## Per-script documentation

### `01_gather_associations.py` — download GWAS top hits from OpenGWAS

Queries the OpenGWAS metadata for all studies, keeps those of the requested
ancestry that pass a minimum sample size, and downloads every association with
p < 5×10⁻⁵ for each study.

| Parameter | Value | Meaning / rationale |
|---|---|---|
| `population_of_interest` | `"European"` or `"African"` | African expands to `["African American or Afro-Caribbean", "Sub-Saharan African", "African unspecified"]`, which are the three OpenGWAS African labels. |
| Minimum `sample_size` | EUR > 10,000; AFR > 5,000 | Removes underpowered studies. The African threshold is lower because far fewer large African GWASs exist. |
| `pval` in `query.tophits` | `5e-5` | A *suggestive* threshold, looser than genome-wide significance (5e-8). It captures sub-threshold signals that are useful for counting associations per haplotype and for PRS. |
| `clump` | `0` | No LD pruning. All hits are returned, because clumping happens later (in `02`), and `05` needs every variant to count hits per haplotype. |

- **Output:** one file per study, `{output_folder}/{study_id}_hits_p5e-5.txt`
  (tab-separated, OpenGWAS `tophits` columns plus `study_id`).
- **Notes:**
  - Studies are de-duplicated on `id`. Failed queries are printed and skipped.
  - The European branch passes a string to `.isin()`. Pandas rejects this, so
    the European run needs `population_of_interest = ["European"]`. The
    European files used downstream are also `.txt.gz` in
    `OpenGWAS_sumstats/European`, not `.txt` in `OpenGWAS_sumstats`, so they
    were compressed and moved after download.

---

### `02_clump_PRS.py` — ancestry-aware clumping + conditional analysis (PRS sets)

Runs [GWASum](/project/holstegelab/Share/nicco/GWASum/GWASum.py) on each
study's top hits. GWASum lifts the variants to hg38, clumps them against the
ancestry-matched reference panel, and runs a GCTA-COJO conditional analysis.
The script then reformats the conditionally independent hits into a standard
PRS weight file.

| Parameter | Value | Meaning / rationale |
|---|---|---|
| `ancestry` | `"European"` / `"African"` | Selects the input sumstats folder and the LD reference genotypes (African: `--geno-pattern chr{chr}_african`). |
| `threads` | `6` | Threads passed to GWASum. |
| `--clumping-p` | `5e-5` | Index-SNP threshold. Matches the download threshold in `01`, so every downloaded hit can be an index SNP. |
| `--clumping-r2` | `0.05` | Strict LD cut-off, so the retained SNPs are close to independent. This is a common choice for PRS. |
| `--clumping-kb` | `500` | Clumping window of ±500 kb. |
| `--cojo-p` | `5e-5` | Significance threshold for COJO stepwise selection. Kept consistent with the clumping threshold. |
| `--cojo-wind` | `1000` | COJO window (kb). Variants more than 1 Mb apart are assumed independent. |
| `--liftover` | on | Converts the hg19 OpenGWAS hits to hg38 so they match the genotype data. |
| Column mapping | `chr, position, ea, nea, beta, se, p, eaf, n` | OpenGWAS `tophits` column names. |

- **Inputs:** EUR `*.txt.gz`, excluding files that start with `All` or `eqtl`
  (the combined file and eQTL studies). AFR: `*.txt`.
- **Outputs:**
  - `GWASum/{study}/...`: raw GWASum output. Studies whose output directory
    already exists are skipped, so the script can be resumed after an
    interruption.
  - `GWASum/reformatted/{trait_id}_clumped_reformatted.txt` with columns
    `SNP, CHROM, POS, EFFECT_ALLELE, OTHER_ALLELE, BETA, SE, P, EFFECT_ALLELE_FREQUENCY`.
    COJO output lacks the non-effect allele, so `OTHER_ALLELE` is recovered
    from the per-locus `conditional_files` (A1/A2). `EFFECT_ALLELE_FREQUENCY`
    is the frequency in the reference genotypes (`freq_geno`).
- **Notes:** The `trait_id` is taken as the name of the directory three levels
  above `combined_conditional_results.txt`, so this relies on GWASum's output
  layout. GWASum's stdout/stderr is discarded (`> /dev/null 2>&1`); if a study
  produces no output, re-run its command by hand to see the error.

---

### `03_R2cluster_EUR.py` / `03_R2cluster_AFR.py` — haplotype block estimation

Both files contain the same algorithm. They differ only in the genotype panel,
the sample restriction and the parallelism.

**Algorithm**

1. **Recombination segmentation** (`recomb_intervals`). Walk along each
   chromosome of the hg38 recombination map. Every time the rate exceeds
   `threshold` cM/Mb, close the current interval and start a new one. The
   result is a set of contiguous, non-overlapping intervals whose boundaries
   are recombination hotspots, so no LD block crosses a hotspot.
2. **LD per interval** (`calculateLD`). `plink2 --r2-unphased square`
   computes the full r² matrix of all variants in the interval.
3. **Clustering** (`ld_clusters_complete`). Distance is `D = 1 − r²`. The
   script runs **complete-linkage** hierarchical clustering and cuts the tree
   at `1 − thresh`. With complete linkage, **every pair of variants in a
   cluster has r² ≥ `thresh`**, which is a strict definition of a haplotype.
4. **Summary** (`cluster_summary`). Each cluster becomes one row:
   `CHR, BP1, BP2, KB, NSNPS, SNPS` (`SNPS` = `|`-joined `CHR:POS:REF:ALT`
   sorted by position). Clusters with fewer than `min_snps` variants are
   dropped.

| Parameter | EUR | AFR | Meaning / rationale |
|---|---|---|---|
| `threshold` (recomb.) | 25 cM/Mb | 25 cM/Mb | Hotspot threshold. The genome-wide average is about 1 cM/Mb, so 25 cM/Mb only marks strong hotspots. A lower value gives more, smaller intervals and faster LD computation, but can split real blocks. |
| `thresh` (r²) | 0.2 | 0.2 | Minimum pairwise r² inside a haplotype. A higher value gives tighter, smaller blocks. A sensitivity run at 0.5 exists (see the EUR output path). |
| `--maf` | 0.005 | 0.005 | Removes very rare variants, whose r² estimates are unstable. |
| `--keep` | `samples_of_interest.txt` | none | EUR: restricts LD to a defined sample subset. AFR: uses all African 1000G samples. |
| `min_snps` | 5 | 5 | Set in `process_interval` (the function default of 2 is overridden). Blocks with 1–4 variants are discarded as uninformative. |
| `n_procs` / `chunk_size` | 6 / 6 | 8 / 8 | Multiprocessing pool size and `imap_unordered` chunk size. |

- **Inputs:** `recombination_rates_hg38.txt.gz` (columns `Chromosome`,
  `Position(bp)`, `Rate(cM/Mb)`) and plink2 pfiles per chromosome.
- **Outputs:**
  - EUR: currently `sensitivity_haplotyping/20260910_thr0.5_haplotypes.tsv`.
    The production file was `20251122_ld_clusters_low_recombination_intervals.tsv.gz`
    (the bgzip/tabix lines are commented out at the end of the script).
  - AFR: `African_haplotypes/20260914_thr0.2_haplotypes.tsv`.
- **Notes:**
  - `calculateLD` reads `thresh` as a **global variable**, not as an argument.
  - The EUR output file name says `thr0.5`, but `thresh = 0.2` in the script.
    Check that the value and the file name agree before running a sensitivity
    analysis.
  - Any failure in an interval (for example, fewer than 2 variants, or a plink
    error) is silently skipped (`except: return None`).
  - Temporary `ld_*` files are written to the current working directory.
  - The square r² matrix grows quadratically with the number of variants, so
    very long low-recombination intervals are the memory bottleneck.

---

### `04_semantic_reduction_AI.py` — semantic de-duplication of GWAS traits

Groups OpenGWAS trait names by meaning and keeps one study per group.

**Steps**

0. **Collect trait names.** Read the trait name from row 2, column 2 of every
   top-hits file. Write `trait_names.txt` and `trait_names_and_ids.txt` into
   the sumstats folder.
1. **Normalise** (`normalize_trait`). Lower-case the name and strip brackets
   and punctuation. Remove boilerplate phrases ("GWAS of", "genome-wide
   association study", "meta-analysis", "study of", "association with") and
   English stop-words. Acronyms are kept.
2. **Embed** (`encode_ensemble`). Encode each name with both models, L2-normalise
   each embedding, scale it by √weight, concatenate the two, and re-normalise.
   With this construction, the cosine similarity in the combined space equals
   the **weighted mean of the per-model cosine similarities**, while the two
   latent spaces stay separate.
3. **Cluster** (`cluster_embeddings`). Agglomerative clustering on the
   pre-computed cosine distance matrix with **average linkage**, cut at
   `distance_threshold`.
4. **Reassign singletons** (`reassign_singletons`). A trait that ends up alone
   is moved to the nearest *multi-member* cluster (by centroid) if its cosine
   similarity to that cluster is ≥ `singleton_min_similarity`. Singletons are
   never merged with each other.
5. **Pick the representative name** (`pick_representatives`). Choose the
   member closest to the cluster centroid (`'centroid'`) or the shortest
   normalised name (`'shortest'`).
6. **Pick the representative study.** Count the associations in every file
   (`trait_association_counts.txt`). In each cluster, the **study with the
   most associations** becomes the `Final Representative`. Its file name
   (`Final Representative file`) is what `05` uses.
7. **Visualise.** Plot the number of clusters and the cluster-size
   distribution per threshold. Plot a UMAP of the embeddings, and print the
   nearest traits and clusters for a query trait as a diagnostic.

| Parameter | Value | Meaning / rationale |
|---|---|---|
| `MODEL_NAMES` | `FremyCompany/BioLORD-2023`, `amrothemich/sapbert-sentence-transformers` | BioLORD is trained on ontology definitions and captures broader biomedical relationships (for example, a disease and its sub-types). SapBERT is trained on UMLS synonyms and is strong at recognising that two strings name the *same* concept. |
| `WEIGHTS` | `[0.7, 0.3]` | BioLORD dominates the grouping. SapBERT adds a pull toward exact synonyms. |
| `distance_threshold` | scanned over `0.1 … 0.7`; **0.5 used downstream** | Cosine-distance cut-off. A lower value gives more, purer clusters; a higher value gives fewer, broader ones (0.3–0.4 is already fairly aggressive). Results are listed in the table below. At 0.4–0.5, for example, all Alzheimer's disease studies fall into one cluster, but they are split at lower thresholds. |
| `singleton_min_similarity` | `0.40` | A permissive, conservative-in-direction rescue: it only ever *joins* singletons to existing clusters. |
| `pick` | `'centroid'` | The most "typical" name is used as the cluster label. |
| `batch_size` | 256 | Batch size used during encoding. |
| UMAP | `n_neighbors=50, min_dist=0.1, metric='cosine', random_state=42` | Only used for visualisation. |

Number of clusters per threshold (from the script's comments):

| Threshold | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 |
|---|---|---|---|---|---|---|---|
| European (10,604 traits) | 1842 | 1662 | 1343 | 978 | **604** | 319 | 136 |
| African (751 traits) | 236 | 197 | 172 | 159 | **136** | 109 | 66 |

- **Outputs** (written to the working directory, prefix `20260916_`, suffix
  `_{ancestry}`):
  - `trait_clusters_{anc}.pkl`: for every threshold, the reps, cluster maps,
    labels, embeddings and distance matrices.
  - `trait_clusters_summary_AI_thresholds_{anc}.csv`: one row per cluster per
    threshold, with `Cluster Label, Representative, Members (___-joined),
    Number of Members, Threshold, Final Representative, Final Representative
    associations, Final Representative file`. This file is written twice: first
    tab-separated, then overwritten as **comma-separated** with the final
    columns. `05` reads the comma-separated version.
  - `trait_clusters_summary_AI_thresholds_{anc}.png`,
    `umap_trait_embeddings_{anc}.png`, `trait_embeddings_umap_thr05_{anc}.pkl`.
- **Notes:**
  - The embeddings do not depend on the threshold, so `ens_list[4]` is the same
    as every other entry. The `labels` and `cluster_map` used for the UMAP
    colours and saved in `..._thr05_...pkl` are, however, the leftovers from
    the last loop iteration, which is **threshold 0.7**, not 0.5. Set
    `labels = labels_list[4]` and `cluster_map = cluster_map_list[4]` before
    the UMAP section if you want 0.5.
  - In `trait_names_and_ids.txt`, the `id` column is built as
    `f.replace('.txt','') + '.txt.gz'`. This does not reproduce the real file
    name for either ancestry.
  - When no member matches the association counts, `Final Representative` is
    set to `0` rather than an empty string.
  - The diagnostic query is currently `"Multiple sclerosis"` (the comments
    still mention Alzheimer's disease).

---

### `05_count_associations_haplotypes.py` — annotate haplotypes with GWAS hits

Takes the representative study per cluster at the chosen semantic threshold,
merges their top hits, lifts them over to hg38, and counts how many
associations fall on the variants of each haplotype.

| Parameter | Value | Meaning / rationale |
|---|---|---|
| `ancestry` | `"European"` / `"African"` | Selects the input haplotypes, the sumstats and the semantic clustering file. |
| `threshold` | `0.5` | Semantic distance threshold used to select the representative studies (604 EUR / 136 AFR studies). |
| Liftover | `pyliftover` hg19 → hg38 | Only the first mapped coordinate is kept. Unmapped variants are dropped. |

- **Matching rule:** an association is assigned to a haplotype if its hg38
  **position** equals the position of one of the haplotype's variants (from
  the `SNPS` column). Positions in the gaps between a haplotype's variants are
  not counted, and alleles are not compared at this stage (see
  `sumstats_allele_match.py`).
- **Outputs** (EUR prefix `20260714_`, AFR prefix `20260917_` with suffix `_African`):
  - `All_indep_gwas_sumstats_AI.txt`: concatenated top hits (hg19).
  - `All_indep_gwas_sumstats_AI_hg38_byTrait.txt.gz` (+ `.tbi`): sorted by
    study `id`, indexed with `tabix -s 1 -b 14 -e 14`. The study id acts as the
    "sequence" name, so `tabix file <study_id>` returns all hits of one study.
  - `All_indep_gwas_sumstats_AI_hg38.txt.gz` (+ `.tbi`): sorted by chromosome,
    indexed with `tabix -s 3 -b 14 -e 14` (column 3 = `chr`, column 14 =
    `position_hg38`) for region queries.
  - `ld_clusters_with_gwas_AI_latest.tsv.gz` (+ `.tbi`): the haplotype table
    plus `assoc_count`, `rsids` (`;`-joined) and `traits` (`___`-joined).
    This is the haplotype file that snpXplorer reads.
  - `ld_clusters_with_gwas_backup_latest.pkl`: the same table with Python
    lists instead of joined strings.
- **Notes:** The EUR semantic file name (`20260714_...csv`) is older than the
  current `04` output name (`20260916_..._European.csv`), so update one of the
  two before re-running. The tabix column indices assume the OpenGWAS
  `tophits` column order plus `study_id` and `position_hg38`.

---

### `06_plot_haplotype_statistics.py` — haplotype QC plots and genome coverage

Descriptive statistics for a haplotype table.

- **Plots** (PDF, prefix `PLOT_PREFIX = "20260917"`, written to `OUTPUT_DIR`):
  - `n_haplotypes_per_chr`: bar plot of the number of haplotypes per chromosome.
  - `haplotype_size_distribution`: violin and box plot of `KB` per chromosome
    (log scale).
  - `n_snps_per_haplotype`: violin and box plot of `NSNPS` per chromosome
    (log scale).
- **Coverage.** Merges overlapping and adjacent `[BP1, BP2]` intervals per
  chromosome and reports the covered bases as a fraction of `GENOME_SIZE_BP`.

| Parameter | Value | Meaning |
|---|---|---|
| `ancestry` | `"European"` / `"African"` | Selects `INPUT_PATH` / `OUTPUT_DIR`. |
| `GENOME_SIZE_BP` | `3_200_000_000` | Approximate haploid genome size used as the coverage denominator. |
| `PLOT_PREFIX` | `"20260917"` | Output file name prefix. |

- **Note:** The African `INPUT_PATH` expects `20260914_thr0.2_haplotypes.tsv.gz`,
  but `03_R2cluster_AFR.py` writes an uncompressed `.tsv`. Run `bgzip` on it
  first, or change the path. Exact duplicate rows are removed before plotting.

---

### `07_align_semantic_clustering.py` — joint UMAP of European and African traits

Stacks the European and African trait embeddings (from the
`..._trait_embeddings_umap_thr05...pkl` files written by `04`) and fits **one**
UMAP (`n_neighbors=50, min_dist=0.1, metric='cosine', random_state=42`).
This puts both ancestries' traits in the same 2-D space. Both ancestries use
the same models and weights, so their embeddings are directly comparable.

- **Plots:**
  - `20260918_combined_umap_clusters.png`: all traits. Each (ancestry, cluster)
    pair has its own colour from the `turbo` palette. European traits are
    circles, African traits are triangles.
  - `20260918_combined_umap_triglycerides_clusters.png`: only the clusters
    with at least one member containing `query` (case-insensitive substring),
    with the cluster representative in the legend. Change `query` to focus on
    another trait (the comment still says hypertension).
- **Output:** `20260918_combined_umap_embeddings.pkl` (coordinates, both
  datasets, colour map). snpXplorer uses this file.
- **Note:** The cluster labels come from the `..._thr05...pkl` files. As
  explained under `04`, these may actually hold the threshold-0.7 labels.

---

### `08_create_LD_db.py` — per-chromosome SQLite LD database

Pre-computes pairwise LD and stores it in compact SQLite files
(`LD_database/{ancestry}/ld_chr{N}.sqlite`) for fast look-ups in the web app.

**Windowing.** Each chromosome is covered by **2 Mb windows with a 1 Mb step**
(`create_intervals`). plink2 computes LD for every pair up to 1 Mb apart
inside each window. Only pairs whose first variant (`POS_A`) lies in the first
1 Mb of the window (the "core") are stored. Together these rules ensure that
every pair ≤ 1 Mb apart is computed with both variants inside one window and
stored exactly once.

| Parameter | Value | Meaning / rationale |
|---|---|---|
| `ancestry` | `"European"` / `"African"` | Selects the genotype panel. |
| Sample restriction (EUR) | healthy individuals with PCs | `--keep pheno_data_eur_healthy.txt` is built from the phenotype file. It keeps samples with non-missing `PC1` and a control/centenarian diagnosis, so that LD reflects a homogeneous, unaffected European sample. AFR uses all samples. |
| `window` / `step` | 1 Mb core, 2 Mb window | See the windowing paragraph above. |
| `--maf` | `0.005` | Same as in `03`. |
| `--ld-window-kb` | `1000` | Maximum pair distance of 1 Mb. |
| `--ld-window` | `99999999` | Effectively no limit on the number of variants between the two members of a pair. |
| `--ld-window-r2` | `0.2` | Only pairs with r² ≥ 0.2 are stored, which keeps the database small. **Pairs that are absent may have r² < 0.2 or may not have been tested.** |
| `batch_size` | `50000` | Number of rows inserted per transaction. |
| Chromosomes | `range(20, 23)` | **Currently only chr20–22.** Set this to `range(1, 23)` for a full build. |

**Schema**

```sql
variants(id INTEGER PRIMARY KEY, uniq TEXT UNIQUE,  -- "POS:A1:A2" (chrom stripped)
         pos INTEGER)
ld(v1 INTEGER, v2 INTEGER,                          -- v1 < v2, variant ids
   r2_milli INTEGER,                                -- round(r2 * 1000)
   dist_bp INTEGER, PRIMARY KEY (v1, v2)) WITHOUT ROWID
-- indexes on variants(pos), ld(v1), ld(v2); ANALYZE run at the end
```

---

### `09_compress_and_check_PRS.py` — finalise PRS files and the trait table

1. Compares the GWASum output folders with the `reformatted/` files for both
   ancestries and reports the studies that failed reformatting (2 African
   studies at the time of writing; these are skipped).
2. Runs `gzip` on every `*_clumped_reformatted.txt`.
3. Concatenates the European and African `trait_names_and_ids.txt` tables,
   adds an `ancestry` column, and writes
   `OpenGWAS_sumstats/combined_traits_and_ids.txt`, which snpXplorer uses.

No tunable parameters apart from the input folders.

---

### `sumstats_allele_match.py` — harmonise effect alleles across studies

In the merged association file from `05`, the same locus can appear with
different allele orders (EA/NEA swapped) or with different allele pairs. This
script makes the alleles consistent in **two streaming passes**, so the file
never has to fit in memory:

1. **Pass 1.** For each locus (`chr`, `position`), count the unordered allele
   pairs and select the **majority pair**. Its first-seen ordering becomes the
   reference orientation. On a tie, the first pair encountered wins, and a
   warning is logged.
2. **Pass 2.** Write every record:
   - same pair and same orientation → keep unchanged;
   - same pair, opposite orientation → **flip**: swap EA/NEA, set
     `beta → −beta` and `eaf → 1 − eaf`;
   - different pair → **remove** (a minority allele pair, for example a
     multi-allelic site or an indel/SNV clash).

   The log reports, as a sanity check, how many flips have EAFs on opposite
   sides of 0.5 (as expected for a true strand/orientation swap). The two
   passes are cross-checked (the expected and observed numbers of removed
   records must agree), and the output is written atomically through a `.tmp`
   file.

| Parameter | Value | Meaning |
|---|---|---|
| `INPUT_FILE` | `Data/databases/haplotypes/All_indep_gwas_sumstats_AI_hg38.txt.gz` | Relative to the **snpXplorer repository root**. Run the script from there. |
| `OUTPUT_FILE` | `..._allele_matches.txt.gz` | Harmonised output. |
| `CHROM/POSITION/EA/NEA/EAF/BETA_INDEX` | `2, 3, 5, 6, 7, 8` | 0-based column indices in the OpenGWAS `tophits` layout. Update them if the column order changes. |
| `PROGRESS_INTERVAL` | `1_000_000` | Log a progress line every N records. |

---

### `haplotype_region_ld.py` — LD matrix and heatmaps for selected haplotypes (CLI)

A stand-alone command-line tool. It does not modify the web app. Give it one or
more labelled haplotypes on the same chromosome. It looks up their variants in
the haplotype table (with tabix when a `.tbi` index exists, otherwise by
scanning the file), maps the variants to the SQLite LD database from `08`, and
exports every stored LD pair between them.

```bash
# Run from this folder; B is a placeholder, so use real BP1-BP2 values from the table
python haplotype_region_ld.py \
  --haplotype A=7:12194394-12246599 \
  --haplotype B=7:<BP1>-<BP2> \
  --haplotype-table ../../../Data/databases/haplotypes/ld_clusters_with_gwas_AI_latest.tsv.gz \
  --ld-directory   ../../../Data/databases/LD_db \
  --output-prefix  results/tmem106b
```

| Argument | Default | Meaning |
|---|---|---|
| `--haplotype LABEL=CHR:START-END` | – | Repeat once per haplotype. `START-END` must exactly match `BP1`/`BP2` of a row in the table. Commas and a `chr` prefix are accepted. Mutually exclusive with `--manifest`. |
| `--manifest FILE` | – | Tab-separated file with the columns `label` and `haplotype_id`. |
| `--haplotype-table` | `<REPO_ROOT>/Data/databases/haplotypes/ld_clusters_with_gwas_AI_latest.tsv.gz` | Haplotype table (output of `05`). |
| `--ld-directory` | `<REPO_ROOT>/Data/databases/LD_db` | Folder with `ld_chr{N}.sqlite`. |
| `--output-prefix` | `./haplotype_region_ld` | Prefix for all outputs. |
| `--missing-value {zero,nan}` | `zero` | Value written to the matrix for pairs that are not in SQLite. |
| `--heatmap-dpi` | `300` | PNG resolution (minimum 72). |
| `--verbose` | off | Debug logging. |

Constraints: labels and IDs must be unique, all haplotypes must be on the same
chromosome, and each ID must match exactly one row.

**Variant mapping.** Variants are matched to the database by position and by
the allele pair in either order. If several database entries match, the one
with the lowest id is used (`mapping_status = ambiguous_first_id`).

**Outputs** (`<prefix>.*`)

| File | Content |
|---|---|
| `ld_matrix.tsv.gz` | Square r² matrix with variants in physical-position order; diagonal = 1. |
| `ld_observed_mask.tsv.gz` | Same shape; 1 if the pair is stored in SQLite, otherwise 0. |
| `ld_observed_pairs.tsv.gz` | Long format: `variant_1, variant_2, r2, distance_bp`. |
| `ld_heatmap.png` | Lower-triangle heatmap in position order. |
| `ld_heatmap_by_haplotype.png` | The same heatmap with variants grouped by haplotype label (natural sort) and coloured annotation strips. A variant that belongs to several haplotypes is plotted once, under the earliest label. |
| `ld_heatmap_genomic_position.png` | "Roof" (Haploview-style) plot on a true genomic-position axis. Missing pairs are drawn navy. |
| `variants.tsv` | Per-variant metadata, haplotype membership and database mapping status. |
| `haplotypes.tsv` | Per-haplotype span, declared number of SNPs and number of variants retrieved. |
| `qc.json` | Sources, counts (mapped/unmapped/ambiguous variants, observed/absent pairs) and plotting conventions. |

The heatmaps use snpXplorer's five LD bands: 0–0.2 navy, >0.2–0.4 light blue,
>0.4–0.6 green, >0.6–0.8 orange, >0.8–1.0 red.

**Notes**
- The database only stores pairs with r² ≥ 0.2 (see `08`). An absent pair
  therefore means "r² < 0.2 **or** not tested", and cannot be told apart from
  a pair that was never computed. This is why the observed mask is written.
- `REPO_ROOT` is `Path(__file__).parents[2]`. From this folder, that resolves
  to `snpXplorer/Exploration/`, not the repository root, so the **default
  paths do not exist**. Pass `--haplotype-table` and `--ld-directory`
  explicitly, as in the example above, or change the line to `parents[3]`.
