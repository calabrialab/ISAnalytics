
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ISAnalytics <a href='https://bioconductor.org/packages/3.12/bioc/html/ISAnalytics.html'><img src='man/figures/isanalytics_logo.png' align="right" height="250" /></a>

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/calabrialab/isanalytics.svg?branch=master)](https://travis-ci.com/calabrialab/isanalytics)
[![codecov](https://codecov.io/gh/calabrialab/ISAnalytics/branch/master/graph/badge.svg)](https://codecov.io/gh/calabrialab/ISAnalytics)
[![R build status -
bioc](https://github.com/calabrialab/isanalytics/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/calabrialab/isanalytics/actions)
[![R build
status](https://github.com/calabrialab/isanalytics/workflows/R-CMD-check/badge.svg)](https://github.com/calabrialab/isanalytics/actions)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
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

# Installation

To install the package from bioconductor:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("ISAnalytics")
```

To install the package from GitHub:

``` r
if (!require(devtools)) {
  install.packages("devtools")
}
devtools::install_github("calabrialab/ISAnalytics", build_vignettes = TRUE)
```

# Current functionality

  - Import integration matrices from files: more info with
    `vignette("How to use import functions", package = "ISAnalytics")`
  - Collision removal: more info with `vignette("Collision removal
    functionality", package = "ISAnalytics")`
  - Aggregation: more info with `vignette("Working with aggregate
    functions", package = "ISAnalytics")`
  - Re-calibration functions: `compute_near_integrations`
  - Analysis functions: `compute_abundance`, `comparison_matrix`, others
    `separate_quant_matrices`
  - Utility functions

# NEWS

# ISAnalytics News

## Changes in version 0.99.11 ()

#### NEW FEATURES

  - Added analysis functions `threshold_filter`, `top_integrations`
  - Added support for multi-quantification matrices in
    `compute_abundance`

#### MINOR FIXES

  - Fixed bug in `comparison_matrix` that ignored custom column names
  - Fixed issues in some documentation pages

## Changes in version 0.99.10 (2020-09-14)

ISanalytics is officially on bioconductor\!

#### NEW FEATURES

  - Added analysis functions `comparison_matrix` and
    `separate_quant_matrices`
  - Added utility function `as_sparse_matrix`
  - Added package logo

#### SIGNIFICANT USER-VISIBLE CHANGES

  - Changed algorithm for `compute_near_integrations`
  - Added support for multi-quantification matrices to
    `remove_collisions`
  - Added usage of lifecycle badges in documentation: users can now see
    if a feature is experimental/maturing/stable etc

#### MINOR FIXES

  - Added fix for `import_single_Vispa2Matrix` to remove non significant
    0 values

## Changes in version 0.99.9 (2020-09-01)

#### NEW FEATURES

  - Added functionality: aggregate functions
  - Added vignette on aggregate functions
  - Added recalibration functions
  - Added first analysis function (compute\_abundance)

#### SIGNIFICANT USER-VISIBLE CHANGES

  - Dropped structure `ISADataFrame`: now the package only uses standard
    tibbles
  - Modified package documentation

## Changes in version 0.99.8 (2020-08-12)

  - Submitted to Bioconductor

# TO DO in future updates

  - [ ] Add vignette for association file usage
  - [ ] Add vignette for re-calibration functionality
  - [ ] Add support for multi-quantification matrices to several
    functions

# Getting help

For help please contact the maintainer of the package or open an issue
on GitHub.
