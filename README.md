
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ISAnalytics <a href='https://bioconductor.org/packages/3.12/bioc/html/ISAnalytics.html'><img src='man/figures/logo.png' align="right" height="250" /></a>

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/calabrialab/isanalytics.svg?branch=master)](https://travis-ci.com/calabrialab/isanalytics)
[![codecov](https://codecov.io/gh/calabrialab/ISAnalytics/branch/master/graph/badge.svg)](https://codecov.io/gh/calabrialab/ISAnalytics)
[![R build status -
bioc](https://github.com/calabrialab/isanalytics/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/calabrialab/isanalytics/actions)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![BioC
status](http://www.bioconductor.org/shields/build/release/bioc/ISAnalytics.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/ISAnalytics)
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

# Visit the package website

You can visit the package website to view documentation, vignettes and
more.

-   For the release version: [ISAnalytics Website
    release](https://calabrialab.github.io/ISAnalytics/)
-   For the devel version: [ISAnalytics Website
    dev](https://calabrialab.github.io/ISAnalytics/dev)

# Installation

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
                         ref = "RELEASE_3_12",
                         dependencies = TRUE,
                         build_vignettes = TRUE)

## Safer option for vignette building issue
devtools::install_github("calabrialab/ISAnalytics",
                         ref = "RELEASE_3_12")
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

## Safer option for vignette building issue
devtools::install_github("calabrialab/ISAnalytics",
                         ref = "master")
```

# NEWS

# ISAnalytics 1.1.11 (2021-05-11)

## NEW FUNCTIONALITY

-   `HSC_population_size_estimate` and `HSC_population_plot` allow
    estimates on hematopoietic stem cell population size
-   Importing of Vispa2 stats per pool now has a dedicated function,
    `import_Vispa2_stats`
-   `outlier_filter` and `outliers_by_pool_fragments` offer a mean to
    filter poorly represented samples based on custom outliers tests

## VISIBLE USER CHANGES

-   The argument `import_stats` of `aggregate_metadata` is officially
    deprecated in favor of `import_Vispa2_stats`
-   `aggregate_metadata` is now a lot more flexible on what operations
    can be performed on columns via the new argument
    `aggregating_functions`
-   `import_association_file` allows directly for the import of Vispa2
    stats and converts time points to months and years where not already
    present
-   File system alignment of `import_association_file` now produces 3
    separate columns for paths
-   `separate_quant_matrices` and `comparison_matrix` now do not require
    mandatory columns other than the quantifications - this allows for
    separation or joining also for aggregated matrices

## FIXES

-   Fixed a minor issue in `CIS_volcano_plot` that caused duplication of
    some labels if highlighted genes were provided in input

# ISAnalytics 1.1.10 (2021-04-08)

## FIXES

-   Fixed issue in `compute_near_integrations`: when provided
    recalibration map export path as a folder now the function works
    correctly and produces an automatically generated file name
-   Fixed issue in `aggregate_metadata`: now paths to folder that
    contains Vispa2 stats is looked up correctly. Also, VISPA2 stats
    columns are aggregated if found in the input data frame
    independently from the parameter `import_stats`.

## IMPROVEMENTS

-   `compute_abundance` can now take as input aggregated matrices and
    has additional parameters to offer more flexibility to the user.
    Major updates and improvements also on documentation and
    reproducible examples.
-   Major improvements in function `import_single_Vispa2Matrix`: import
    is now preferentially carried out using `data.table::fread` greatly
    speeding up the process - where not possible `readr::read_delim` is
    used instead
-   Major improvements in function `import_association_file`: greatly
    improved parsing precision (each column has a dedicated type),
    import report now signals parsing problems and their location and
    signals also problems in parsing dates. Report also includes
    potential problems in column names and signals missing data in
    important columns. Added also the possibility to give various file
    formats in input including `*.xls(x)` formats.
-   Function `top_integrations` can now take additional parameters to
    compute top n genes for each specified group
-   Removed faceting parameters in `CIS_volcano_plot` due to poor
    precision (easier to add faceting manually) and added parameters to
    return the data frame that generated the plot as an additional
    result. Also, it is now possible to specify a vector of gene names
    to highlight even if they’re not above the annotation threshold.

## MINOR

-   ISAnalytics website has improved graphic theme and has an additional
    button on the right that leads to the devel (or release) version of
    the website
-   Updated vignettes

## FOR DEVS ONLY

-   Complete rework of test suite to be compliant to testthat v.3

# ISAnalytics 1.1.9 (2021-02-17)

## FIXES

-   Fixed minor issues in internal functions with absolute file paths &
    corrected typos

# ISAnalytics 1.1.8 (2020-02-15)

## FIXES

-   Fixed minor issues in internal functions to optimize file system
    alignment

# ISAnalytics 1.1.7 (2020-02-10)

## FIXES

-   Fixed minor issues in import\_association\_file when checking
    parameters

# ISAnalytics 1.1.6 (2020-02-06)

## UPGRADES

-   It is now possible to save html reports to file from
    import\_parallel\_Vispa2Matrices\_auto and
    import\_parallel\_Vispa2Matrices\_interactive, remove\_collisions
    and compute\_near\_integrations

## FIXES

-   Fixed sample\_statistics: now functions that have data frame output
    do not produce nested tables. Flat tables are ready to be saved to
    file or can be nested.
-   Simplified association file check logic in remove\_collisions: now
    function blocks only if the af doesn’t contain the needed columns

# ISAnalytics 1.1.5 (2020-02-03)

## UPGRADES

-   Upgraded import\_association\_file function: now file alignment is
    not mandatory anymore and it is possible to save the html report to
    file
-   Updated vignettes and documentation

# ISAnalytics 1.1.4 (2020-11-16)

## UPGRADES

-   Greatly improved reports for collision removal function
-   General improvements for all widget reports

# ISAnalytics 1.1.3 (2020-11-10)

## FIXES

-   Further fixes for printing reports when widgets not available
-   Added progress bar to collision processing in `remove_collisions`
-   Updated vignettes

## NEW

-   Added vignette “Using ISAnalytics without RStudio support”

# ISAnalytics 1.1.2 (2020-11-05)

## FIXES

-   Fixed missing restarts for non-blocking widgets

# ISAnalytics 1.1.1 (2020-11-04)

## FIXES

-   Functions that make use of widgets do not interrupt execution
    anymore if errors are thrown while producing or printing the widgets
-   Optimized widget printing for importing functions
-   If widgets can’t be printed and verbose option is active, reports
    are now displayed on console instead (needed for usage in
    environments that do not have access to a browser)
-   Other minor fixes (typos)
-   Bug fixes: fixed a few bugs in importing and recalibration functions
-   Minor fix in import\_association\_file file function: added multiple
    strings to be translated as NA

## IMPORTANT NOTES

-   Vignette building might fail due to the fact that package
    “knitcitations” is temporarily unavailable through CRAN
-   ISAnalytics is finally in release on bioconductor!

# ISAnalytics 0.99.14 (2020-10-21)

-   Minor fixes in tests

# ISAnalytics 0.99.13 (2020-10-19)

## NEW FEATURES

-   Added analysis functions `CIS_grubbs` and `cumulative_count_union`
-   Added plotting functions `CIS_volcano_plot`

# ISAnalytics 0.99.12 (2020-10-04)

## NEW FEATURES

-   Added analysis function `sample_statistics`

## SIGNIFICANT USER-VISIBLE CHANGES

-   `aggregate_values_by_key` has a simplified interface and supports
    multi-quantification matrices

## MINOR CHANGES

-   Updated vignettes
-   `import_parallel_Vispa2Matrices_interactive` and
    `import_parallel_Vispa2Matrices_auto` now have an option to return a
    multi-quantification matrix directly after import instead of a list

# ISAnalytics 0.99.11 (2020-09-21)

## NEW FEATURES

-   Added analysis functions `threshold_filter`, `top_integrations`
-   Added support for multi-quantification matrices in
    `compute_abundance`

## MINOR FIXES

-   Fixed bug in `comparison_matrix` that ignored custom column names
-   Fixed issues in some documentation pages

# ISAnalytics 0.99.10 (2020-09-14)

ISanalytics is officially on bioconductor!

## NEW FEATURES

-   Added analysis functions `comparison_matrix` and
    `separate_quant_matrices`
-   Added utility function `as_sparse_matrix`
-   Added package logo

## SIGNIFICANT USER-VISIBLE CHANGES

-   Changed algorithm for `compute_near_integrations`
-   Added support for multi-quantification matrices to
    `remove_collisions`
-   Added usage of lifecycle badges in documentation: users can now see
    if a feature is experimental/maturing/stable etc

## MINOR FIXES

-   Added fix for `import_single_Vispa2Matrix` to remove non significant
    0 values

# ISAnalytics 0.99.9 (2020-09-01)

## NEW FEATURES

-   Added functionality: aggregate functions
-   Added vignette on aggregate functions
-   Added recalibration functions
-   Added first analysis function (compute\_abundance)

## SIGNIFICANT USER-VISIBLE CHANGES

-   Dropped structure `ISADataFrame`: now the package only uses standard
    tibbles
-   Modified package documentation

# ISAnalytics 0.99.8 (2020-08-12)

-   Submitted to Bioconductor

# TO DO in future updates

-   [ ] Add further plot wrappers
-   [x] New functionality: raw reads flagging/filtering (ISAnalytics
    1.1.11)
-   [x] New functionality: HSCs population size estimate (ISAnalytics
    1.1.11)
-   [ ] New functionality: genomic annotation
-   [ ] Obtain a system to visualize interactive reports based on shiny
    (possibly replacing self-contained HTML reports)

# Getting help

For help please contact the maintainer of the package or open an issue
on GitHub.
