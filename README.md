---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->





# specmine

<!-- badges: start -->
[![](https://www.r-pkg.org/badges/version/specmine?color=green)](https://cran.r-project.org/package=specmine)
<!-- badges: end -->

The goal of *specmine* is to provide a set of methods for metabolomics
data analysis, including data loading in different formats,
pre-processing, metabolite identification, univariate and multivariate
data analysis, machine learning, feature selection and pathway analysis.
Case studies can be found on the website:
<http://bio.di.uminho.pt/metabolomicspackage/index.html>. This package
suggests 'rcytoscapejs', a package not in mainstream repositories. If
you need to install it, use:
`devtools::install_github('cytoscape/r-cytoscape.js@v0.0.7')`.

## Installation

You can install the released version of *specmine* from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("specmine")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("BioSystemsUM/specmine")
```

## Example

This is a basic example which shows you how to load the namespace of
*specmine* and add it to your search list:

``` r
library(specmine)
```

## WebSpecmine — Running the App

This repository also includes an extended version of the WebSpecmine
Shiny application with new unsupervised learning capabilities (UMAP,
t-SNE, ICA, DBSCAN, HDBSCAN, GMM, Hierarchical Clustering), a
preprocessing module, and an embedding comparison dashboard.

### Prerequisites

**R** (≥ 4.1.0) and optionally **RStudio**.

### 1. Install the required packages

``` r
install.packages(c(
  "shiny",
  "shinydashboard",
  "plotly",
  "DT",
  "readxl",
  "uwot",
  "Rtsne",
  "fastICA",
  "dbscan",
  "mclust"
))
```

Install **specmine** from GitHub (see above) if not already installed.

### 2. Clone this repository

``` bash
git clone https://github.com/PedroFontao/Projeto-Bioinf.git
cd Projeto-Bioinf
```

### 3. Run the app

**Option A — RStudio**

Open `app.R` in RStudio and click **Run App**.

**Option B — R console**

``` r
shiny::runApp("app.R")
```

**Option C — command line**

``` bash
Rscript -e "shiny::runApp('app.R')"
```

The app will open in your browser at `http://127.0.0.1:XXXX`.

### 4. Quick start

1. Go to the **Load Workspace** tab
2. Select a built-in dataset (iris, wine, cachexia, propolis) or upload your own `.csv` / `.xlsx` file
3. Use the **Preprocessing** tab to clean and prepare the data
4. Run embeddings (PCA, UMAP, t-SNE, ICA) in the **Embeddings** tab
5. Apply clustering (DBSCAN, HDBSCAN, GMM, HC) in the **Clustering** tab
6. Compare embedding quality metrics in the **Comparison** tab

> **Note:** The new WebSpecmine modules currently operate as a standalone
> extension and are not yet fully merged into the upstream specmine
> R package. Supported upload formats: `.csv`, `.tsv`, `.txt`, `.xlsx`, `.xls`.
