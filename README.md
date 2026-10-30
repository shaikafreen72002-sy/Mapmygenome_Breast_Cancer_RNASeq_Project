# Breast Cancer RNA-Seq Analysis

This repository contains small Python and R examples of exploratory analysis and differential expression for two breast-cancer-like sample groups: **Luminal A** and **Basal-like**.

> **Important:** The analysis scripts currently generate deterministic, synthetic expression/count matrices. They do **not** read the GEO files in `data/`. The results are illustrative and are not findings from the GEO cohorts or for clinical use.

## Project contents

```text
data/
  GSE62944_series_matrix.txt.gz
  GSE62944_series_matrix.txt/
    GSE62944_series_matrix.txt
  GSE81538_series_matrix.txt.gz
  GSE81538_series_matrix.txt/
    GSE81538_series_matrix.txt
scripts/
  eda_and_plots.py
  deseq2_analysis.R
outputs/
  pca_plot.png
  volcano_plot.png
  significant_degs.csv
  differential_expression_results.csv
```

- `data/` contains GEO Series Matrix files for [GSE62944](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62944) and [GSE81538](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81538). These files are included as data resources but are not consumed by the scripts yet.
- `scripts/eda_and_plots.py` generates a 1,000-gene by 30-sample synthetic expression matrix, with 15 samples in each subtype. It runs PCA, performs per-gene Welch's t-tests, and writes the PCA plot, volcano plot, and significant-gene table.
- `scripts/deseq2_analysis.R` generates a separate synthetic integer count matrix of the same dimensions, then runs DESeq2 and exports the results.
- `outputs/` contains example output files. Running either script regenerates its corresponding outputs.

## Requirements

### Python analysis

Python 3 and the following packages are required:

- `matplotlib`
- `numpy`
- `pandas`
- `scipy`
- `scikit-learn`
- `seaborn`

For example, create an environment and install the dependencies from the repository root:

```powershell
py -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install matplotlib numpy pandas scipy scikit-learn seaborn
```

Run the analysis from the repository root:

```powershell
python scripts/eda_and_plots.py
```

The script writes `outputs/pca_plot.png`, `outputs/volcano_plot.png`, and `outputs/significant_degs.csv`.

### R analysis

Install R and the Bioconductor package `DESeq2`. If needed, install it from an R session:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("DESeq2")
```

Run the analysis from the repository root:

```powershell
Rscript scripts/deseq2_analysis.R
```

The script writes `outputs/differential_expression_results.csv`.

## Analysis notes

- Both scripts use a fixed random seed (`42`), so their synthetic inputs and results are reproducible with compatible software versions.
- The Python script selects genes using an unadjusted `p_value < 0.05` and an absolute `Log2FC > 1` threshold. Its `Log2FC` column is the difference of group means on the generated expression values; no log transformation is performed.
- The R script uses DESeq2, which reports `log2FoldChange`, `pvalue`, and multiple-testing-adjusted `padj` values. With the current subtype ordering, the reported contrast is Basal-like relative to Luminal A.
- The two scripts create different synthetic data and use different statistical methods. Their outputs should not be compared as if they were analyses of the same samples.
- Adapting this project to real cohort analysis requires loading and preprocessing the GEO data, matching sample metadata to expression/count columns, and using an appropriate study design and statistical workflow.
