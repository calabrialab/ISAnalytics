
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ISAnalytics

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/calabrialab/isanalytics.svg?branch=master)](https://travis-ci.com/calabrialab/isanalytics)
[![codecov](https://codecov.io/gh/calabrialab/ISAnalytics/branch/master/graph/badge.svg)](https://codecov.io/gh/calabrialab/ISAnalytics)
[![R build status -
bioc](https://github.com/calabrialab/isanalytics/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/calabrialab/isanalytics/actions)
[![R build
status](https://github.com/calabrialab/isanalytics/workflows/R-CMD-check/badge.svg)](https://github.com/calabrialab/isanalytics/actions)
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

ISAnalytics is currently under development with the goal of being
published on Bioconductor soon.

# Installation

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
  - Analysis functions: `compute_abundance`
  - Utility functions

# NEWS

# ISAnalytics News

## Changes in version 0.99.8 (2020-08-12)

  - Submitted to Bioconductor

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

# TO DO in future updates

  - Add vignette for association file usage

# Getting help

For help please contact the maintainer of the package or open an issue
on GitHub.
