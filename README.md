# GSE156171 DESeq2 Analysis

This repository contains a reproducible RNA-seq differential expression workflow for **GSE156171**, focused on the comparison **Infiltrating vs Other** basal cell carcinoma subtypes.

## Overview

The analysis uses:

- **DESeq2** for differential expression
- **EnhancedVolcano** for volcano plots
- **pheatmap** for clustering and heatmaps
- **clusterProfiler** for GO and KEGG enrichment
- **org.Hs.eg.db** and **AnnotationDbi** for gene annotation

## Cohort

Source dataset:

- **GEO accession:** GSE156171
- **Title:** *Subtype specific analyses reveal Infiltrative Basal Cell Carcinoma are highly interactive with their environment*

Subtype grouping used in this analysis:

- **Infiltrating** -> case group
- **Superficial, Solid, Mixed** -> collapsed into **Other**

## Input files

Expected local inputs:

```text
GSE156171_raw-counts.csv.gz
GSE156171_DE_InfiltVsOther_results.csv.gz
GSE156171_metadata.csv
```

## Repository structure

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
    └── raw/   # optional, usually excluded from version control
```

## Metadata format

The metadata file must contain at least:

```csv
sample_id,Sample_geo_accession,subtype
```

Example:

```csv
BCC53,GSM4725965,Superficial
BCC18,GSM4725970,Infiltrating
```

## How to run

Install dependencies in R:

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

Run the pipeline:

```bash
Rscript scripts/Rscript_DESeq2.R
```

## Main outputs

- PCA by collapsed group and original subtype
- Sample-to-sample distance heatmap
- DESeq2 results for **Infiltrating vs Other**
- Volcano plot
- Heatmap of top DE genes
- Heatmap of selected genes of interest
- GO enrichment for upregulated and downregulated genes
- KEGG enrichment for upregulated and downregulated genes
- Additional contrast: **Infiltrating vs Superficial**

## Notes

- The script detects whether gene identifiers are **SYMBOL** or **ENSEMBL**.
- If ENSEMBL identifiers contain version suffixes, they are cleaned before analysis.
- The raw counts table may contain annotation columns such as `Chromosome`, `Gene_Start`, `Gene_End`, `Gene_Name`, and `Gene_Type`; the script keeps only columns matching the metadata sample names.

## Suggested .gitignore entries

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

## Citation

If you use this workflow, cite the GEO series and the associated article.
