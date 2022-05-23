# ISAnalytics <img src='man/figures/logo.png' align="right"/>

<!-- badges: start -->
[![codecov](https://codecov.io/gh/calabrialab/ISAnalytics/branch/master/graph/badge.svg)](https://codecov.io/gh/calabrialab/ISAnalytics)
[![R build status - bioc](https://github.com/calabrialab/isanalytics/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/calabrialab/isanalytics/actions)
[![BioC status](http://www.bioconductor.org/shields/build/devel/bioc/ISAnalytics.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/ISAnalytics)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

ISAnalytics is an R package developed to analyze gene therapy vector
insertion sites data identified from genomics next generation sequencing
reads for clonal tracking studies.

In gene therapy, stem cells are modified using viral vectors to deliver
the therapeutic transgene and replace functional properties since the
genetic modification is stable and inherited in all cell progeny. The
retrieval and mapping of the sequences flanking the virus-host DNA
junctions allows the identification of insertion sites (IS), essential
for monitoring the evolution of genetically modified cells in vivo. A
comprehensive toolkit for the analysis of IS is required to foster
clonal trackign studies and supporting the assessment of safety and long
term efficacy in vivo. This package is aimed at (1) supporting
automation of IS workflow, (2) performing base and advance analysis for
IS tracking (clonal abundance, clonal expansions and statistics for
insertional mutagenesis, etc.), (3) providing basic biology insights of
transduced stem cells in vivo.

# Installation and options

`ISAnalytics` can be installed quickly in different ways:

-   You can install it via [Bioconductor](http://bioconductor.org)
-   You can install it via GitHub using the package `devtools`

There are always 2 versions of the package active:

-   `RELEASE` is the latest stable version
-   `DEVEL` is the development version, it is the most up-to-date
    version where all new features are introduced

## Installation from bioconductor

RELEASE version:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ISAnalytics")
```

DEVEL version:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
# The following initializes usage of Bioc devel
BiocManager::install(version='devel')
BiocManager::install("ISAnalytics")
```

## Installation from GitHub

RELEASE:

``` r
if (!require(devtools)) {
  install.packages("devtools")
}
devtools::install_github("calabrialab/ISAnalytics",
                         ref = "RELEASE_3_14",
                         dependencies = TRUE,
                         build_vignettes = TRUE)
```

DEVEL:

``` r
if (!require(devtools)) {
  install.packages("devtools")
}
devtools::install_github("calabrialab/ISAnalytics",
                         ref = "master",
                         dependencies = TRUE,
                         build_vignettes = TRUE)
```

## Setting options

`ISAnalytics` has a verbose option that allows some functions to print
additional information to the console while they’re executing. To
disable this feature do:

``` r
# DISABLE
options("ISAnalytics.verbose" = FALSE)
# ENABLE
options("ISAnalytics.verbose" = TRUE)
```

Some functions also produce report in a user-friendly HTML format, to
set this feature:

``` r
# DISABLE HTML REPORTS
options("ISAnalytics.reports" = FALSE)
# ENABLE HTML REPORTS
options("ISAnalytics.reports" = TRUE)
```

# Getting help

For help please contact the maintainer of the package or open an issue
on GitHub.
