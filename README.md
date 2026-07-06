
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

This package suggests `rcytoscapejs`, which is not available from
mainstream repositories. If needed, install it with:

``` r
devtools::install_github("cytoscape/r-cytoscape.js@v0.0.7")
#> Warning: `install_github()` was deprecated in devtools 2.5.0.
#> ℹ Please use pak::pak("user/repo") instead.
#> This warning is displayed once per session.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
#> Using GitHub PAT from the git credential store.
#> Skipping install of 'rcytoscapejs' from a github remote, the SHA1 (8649c7c2) has not changed since last install.
#>   Use `force = TRUE` to force installation
```

Installation Install the released version from CRAN:

``` r
install.packages("specmine")
#> Installing package into '/private/var/folders/y9/p5k4znxj62s8385cyxkjgf6w0000gn/T/RtmpVnMLTT/temp_libpath89dfa6e0134'
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
#> Using GitHub PAT from the git credential store.
#> Downloading GitHub repo PedroFontao/specmine@HEAD
#> Skipping 2 packages not available: impute, genefilter
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#> * checking for file ‘/private/var/folders/y9/p5k4znxj62s8385cyxkjgf6w0000gn/T/RtmpsOFoyD/remotes15d345638edec/PedroFontao-specmine-826f9ab/DESCRIPTION’ ... OK
#> * preparing ‘specmine’:
#> * checking DESCRIPTION meta-information ... OK
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> Removed empty directory ‘specmine/system_data’
#> * building ‘specmine_3.1.7.tar.gz’
#> Installing package into '/private/var/folders/y9/p5k4znxj62s8385cyxkjgf6w0000gn/T/RtmpVnMLTT/temp_libpath89dfa6e0134'
#> (as 'lib' is unspecified)
```

Example Load the package:

``` r
library(specmine)
```
