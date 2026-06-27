# GSE156171 RNA-seq analysis with DESeq2

This repository contains a reproducible RNA-seq analysis of **GSE156171**, a cohort of human basal cell carcinoma samples profiled by bulk RNA-seq. The main goal of this workflow is to reproduce and extend the comparison between **Infiltrating** tumors and the remaining histological subtypes grouped as **Other**.

The analysis was implemented in **R** using **DESeq2** for differential expression, together with visualization and enrichment tools for exploratory and functional interpretation.

## Study overview

The dataset comes from GEO accession **GSE156171**, associated with the study *Subtype specific analyses reveal Infiltrative Basal Cell Carcinoma are highly interactive with their environment*. In this cohort, tumors were classified into histological subtypes including **Infiltrating**, **Superficial**, **Solid**, and **Mixed**.

To mirror the contrast provided in the GEO supplementary differential expression file, this repository collapses the subtypes into two groups:

- **Infiltrating**
- **Other** = Superficial + Solid + Mixed

This design allows direct comparison with the published `InfiltVsOther` analysis while keeping the full metadata available for subtype-level exploration.

## Repository structure

The repository is organized as follows:

```text
.
├── README.md
├── scripts/
│   └── Rscript_DESeq2.R
├── metadata/
│   └── GSE156171_metadata.csv
├── results/
│   ├── 00_metadata_used.csv
│   ├── 01_PCA_group.png
│   ├── 02_PCA_subtype.png
│   ├── 03_sample_distance_heatmap.png
│   ├── 04_DESeq2_results_Infiltrating_vs_Other.csv
│   ├── 05_DESeq2_significant_genes.csv
│   ├── 06_GEO_original_DE_results_copy.csv
│   ├── 07_MA_plot.png
│   ├── 08_volcano_plot.png
│   ├── 09_heatmap_top50_DE_genes.png
│   ├── 10_heatmap_genes_of_interest.png
│   ├── 11_UP_GO_BP.csv
│   ├── 11_UP_GO_BP_dotplot.png
│   ├── 11_UP_GO_MF.csv
│   ├── 11_UP_GO_MF_dotplot.png
│   ├── 11_UP_GO_CC.csv
│   ├── 11_UP_GO_CC_dotplot.png
│   ├── 11_UP_KEGG.csv
│   ├── 11_UP_KEGG_dotplot.png
│   ├── 12_DOWN_GO_BP.csv
│   ├── 12_DOWN_GO_BP_dotplot.png
│   ├── 12_DOWN_GO_MF.csv
│   ├── 12_DOWN_GO_MF_dotplot.png
│   ├── 12_DOWN_GO_CC.csv
│   ├── 12_DOWN_GO_CC_dotplot.png
│   ├── 12_DOWN_KEGG.csv
│   ├── 12_DOWN_KEGG_dotplot.png
│   ├── 13_DESeq2_Infiltrating_vs_Superficial.csv
│   ├── 14_summary.txt
│   └── 15_sessionInfo.txt
└── data/
    └── raw/
```

## Input data

The workflow expects three local input files:

```text
GSE156171_raw-counts.csv.gz
GSE156171_DE_InfiltVsOther_results.csv.gz
GSE156171_metadata.csv
```

The count matrix includes gene-level counts plus annotation columns such as chromosome and gene coordinates. The script automatically keeps only the columns matching the sample names listed in the metadata file.

The metadata file must contain at least the following columns:

```csv
sample_id,Sample_geo_accession,subtype
```

Example:

```csv
BCC53,GSM4725965,Superficial
BCC18,GSM4725970,Infiltrating
```

## Analysis workflow

The script performs the following steps:

1. Reads the metadata and count matrix.
2. Detects the gene identifier format.
3. Keeps only sample columns present in the metadata.
4. Collapses histological subtypes into **Infiltrating** versus **Other**.
5. Filters low-count genes.
6. Runs **DESeq2** differential expression analysis.
7. Generates quality-control and clustering plots.
8. Produces a ranked differential expression table and a significance-filtered subset.
9. Performs **GO** and **KEGG** enrichment analyses for upregulated and downregulated genes.
10. Runs an additional subtype-specific contrast: **Infiltrating vs Superficial**.

## Exploratory visualization

A first useful check is the PCA after variance-stabilizing transformation. In this repository, PCA is shown both for the collapsed two-group design and for the original subtype labels.

The grouped PCA helps assess whether infiltrating tumors tend to separate from the other samples globally:

![PCA by group](results/01_PCA_group.png)

Looking at the original labels can also show whether the subtype structure is more nuanced than the binary comparison used in the main DESeq2 model:

![PCA by original subtype](results/02_PCA_subtype.png)

To complement PCA, the sample-to-sample distance heatmap provides a more direct view of clustering and overall similarity between tumors:

![Sample distance heatmap](results/03_sample_distance_heatmap.png)

## Differential expression results

The main differential expression contrast is **Infiltrating vs Other**. Full DESeq2 results are saved in:

- `results/04_DESeq2_results_Infiltrating_vs_Other.csv`

A filtered table of significant genes is available in:

- `results/05_DESeq2_significant_genes.csv`

The MA plot provides an overview of expression change across the dynamic range of the dataset:

![MA plot](results/07_MA_plot.png)

The volcano plot highlights the genes with the strongest combination of effect size and statistical significance in infiltrating tumors relative to the other subtypes:

![Volcano plot](results/08_volcano_plot.png)

## Expression heatmaps

To summarize the most informative genes from the main comparison, the script generates a heatmap of the top differentially expressed genes after transformation and row-wise scaling:

![Top DE genes heatmap](results/09_heatmap_top50_DE_genes.png)

Because the original paper discusses subtype-associated stromal and microenvironment-related genes, the workflow also includes a focused heatmap for selected markers such as **WISP1**, **POSTN**, **COL1A1**, **COL1A2**, **FN1**, **ITGA5**, and **ITGB1** whenever they are present in the results:

![Genes of interest heatmap](results/10_heatmap_genes_of_interest.png)

## Functional enrichment

To move beyond single-gene interpretation, the script performs over-representation analysis separately for genes upregulated and downregulated in infiltrating tumors.

For upregulated genes, the outputs include:

- `results/11_UP_GO_BP.csv`
- `results/11_UP_GO_MF.csv`
- `results/11_UP_GO_CC.csv`
- `results/11_UP_KEGG.csv`

and their corresponding dotplots.

For downregulated genes, the outputs include:

- `results/12_DOWN_GO_BP.csv`
- `results/12_DOWN_GO_MF.csv`
- `results/12_DOWN_GO_CC.csv`
- `results/12_DOWN_KEGG.csv`

and their corresponding dotplots.

An example of the enrichment visualization is shown below:

![GO enrichment example](results/11_UP_GO_BP_dotplot.png)

These analyses help connect the DE results to broader biological programs, especially those related to extracellular matrix remodeling, signaling, and tumor-microenvironment interaction.

## Additional comparison

Besides the main collapsed contrast, the script also computes a direct **Infiltrating vs Superficial** comparison:

- `results/13_DESeq2_Infiltrating_vs_Superficial.csv`

This file can be useful when a more subtype-specific interpretation is preferred over the broader infiltrating-versus-all-others framework.

## How to run

Install the required R packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

install.packages(c(
  "tidyverse",
  "data.table",
  "pheatmap",
  "ggrepel",
  "RColorBrewer"
))

BiocManager::install(c(
  "DESeq2",
  "EnhancedVolcano",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "ashr"
))
```

Then run the pipeline from the repository root:

```bash
Rscript scripts/Rscript_DESeq2.R
```

## Notes on implementation

A few implementation details are worth noting:

- The script detects whether gene identifiers are stored as **SYMBOL** or **ENSEMBL**.
- If Ensembl identifiers include version suffixes, these suffixes are removed before downstream analysis.
- Because the raw count file contains both sample counts and annotation fields, the script selects only the columns matching the metadata sample IDs.
- The main design formula is based on the collapsed group variable rather than the original multi-class subtype annotation.
- No responder/non-responder comparison is included in this repository, because that annotation is not part of the metadata used here.

## Reproducibility files

Two additional files are included for reproducibility and traceability:

- `results/14_summary.txt`
- `results/15_sessionInfo.txt`

The summary file provides a compact description of the analysis output, while `sessionInfo` records the software environment used to run the pipeline.

## Suggested `.gitignore`

A minimal `.gitignore` for this repository is:

```gitignore
.Rhistory
.RData
.Rproj.user/
renv/library/
renv/staging/
*.Rproj
*.gz
/data/raw/
```

If you want to version compressed GEO input files directly in the repository, remove the `*.gz` rule.

## Citation

Villani R, Murigneux V, Alexis J, Sim SL et al. Subtype-Specific Analyses Reveal Infiltrative Basal Cell Carcinomas Are Highly Interactive with their Environment. J Invest Dermatol 2021 Oct;141(10):2380-2390. PMID: 33865912


GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi
ID_GEO: GSE156171
