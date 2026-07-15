
<!-- README.md is generated from README.Rmd. Please edit that file -->

# specmine

<!-- badges: start -->

[![](https://www.r-pkg.org/badges/version/specmine?color=green)](https://cran.r-project.org/package=specmine)
<!-- badges: end -->

*specmine* is an R package for metabolomics and spectral data analysis.
It provides methods for data import, preprocessing, visualization,
univariate and multivariate analysis, machine learning, feature
selection, and pathway analysis.

Project repository: <https://github.com/PedroFontao/specmine>

This package suggests `cyjShiny`, which is available from CRAN. If
needed, the development version can also be installed from GitHub:

``` r
install.packages("cyjShiny")
#> Installing package into '/private/var/folders/y9/p5k4znxj62s8385cyxkjgf6w0000gn/T/RtmpVnMLTT/temp_libpath89df4bf3d01c'
#> (as 'lib' is unspecified)
#> 
#> The downloaded binary packages are in
#>  /var/folders/y9/p5k4znxj62s8385cyxkjgf6w0000gn/T//RtmpZZ3G1y/downloaded_packages
# or
remotes::install_github("cytoscape/cyjShiny", ref = "master", build_vignettes = TRUE)
#> Using GitHub PAT from the git credential store.
#> Downloading GitHub repo cytoscape/cyjShiny@master
#> Skipping 1 packages not available: graph
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#> * checking for file ‘/private/var/folders/y9/p5k4znxj62s8385cyxkjgf6w0000gn/T/RtmpZZ3G1y/remotes1b445a6e18e7/cytoscape-cyjShiny-594c2a7/DESCRIPTION’ ... OK
#> * preparing ‘cyjShiny’:
#> * checking DESCRIPTION meta-information ... OK
#> * installing the package (it is needed to build vignettes)
#> * creating vignettes ... OK
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> Omitted ‘LazyData’ from DESCRIPTION
#> * building ‘cyjShiny_1.0.43.tar.gz’
#> Installing package into '/private/var/folders/y9/p5k4znxj62s8385cyxkjgf6w0000gn/T/RtmpVnMLTT/temp_libpath89df4bf3d01c'
#> (as 'lib' is unspecified)
```

Installation Install the released version from CRAN:

``` r
install.packages("specmine")
#> Installing package into '/private/var/folders/y9/p5k4znxj62s8385cyxkjgf6w0000gn/T/RtmpVnMLTT/temp_libpath89df4bf3d01c'
#> (as 'lib' is unspecified)
#> Warning: package 'specmine' is not available for this version of R
#> 
#> A version of this package for your version of R might be available elsewhere,
#> see the ideas at
#> https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages
```

Install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("PedroFontao/specmine")
#> Warning: `install_github()` was deprecated in devtools 2.5.0.
#> ℹ Please use pak::pak("user/repo") instead.
#> This warning is displayed once per session.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
#> Using GitHub PAT from the git credential store.
#> Downloading GitHub repo PedroFontao/specmine@HEAD
#> Skipping 2 packages not available: impute, genefilter
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#> * checking for file ‘/private/var/folders/y9/p5k4znxj62s8385cyxkjgf6w0000gn/T/RtmpZZ3G1y/remotes1b442478d82d/PedroFontao-specmine-c4aa037/DESCRIPTION’ ... OK
#> * preparing ‘specmine’:
#> * checking DESCRIPTION meta-information ... OK
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> * building ‘specmine_3.1.8.tar.gz’
#> Installing package into '/private/var/folders/y9/p5k4znxj62s8385cyxkjgf6w0000gn/T/RtmpVnMLTT/temp_libpath89df4bf3d01c'
#> (as 'lib' is unspecified)
```

Example Load the package:

``` r
library(specmine)
```
