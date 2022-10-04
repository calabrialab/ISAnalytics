
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ISAnalytics <a href='https://bioconductor.org/packages/3.12/bioc/html/ISAnalytics.html'><img src='man/figures/logo.png' align="right" height="250" /></a>

<!-- badges: start -->

[![codecov](https://codecov.io/gh/calabrialab/ISAnalytics/branch/master/graph/badge.svg)](https://codecov.io/gh/calabrialab/ISAnalytics)
[![R build status -
bioc](https://github.com/calabrialab/isanalytics/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/calabrialab/isanalytics/actions)
[![BioC
status](http://www.bioconductor.org/shields/build/devel/bioc/ISAnalytics.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/ISAnalytics)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
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
                         ref = "RELEASE_3_15",
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

# NEWS

<details>
<summary>
Show more
</summary>

# ISAnalytics 1.7.4 (2022-10-04)

## VISIBLE USER CHANGES

-   Progress bars for long processing functions are now implemented via
    the package `progressr`, added a wrapper function for fast enabling
    progress bars, `enable_progress_bars()`
-   Introduced logging for issues in `HSC_population_size_estimate()` -
    signals eventual problems in computing estimates and why

## BUG FIXES AND MINOR CHANGES

-   Fixed minor bugs and typos

# ISAnalytics 1.7.3 (2022-06-17)

## BUG FIXES AND MINOR CHANGES

-   All functions that check for options now have a default value if
    option is not set
-   `CIS_grubbs` function is now faster (removed dependency from
    `psych::describe`)

## NEW

-   New functions `CIS_grubbs_overtime()` and associated plotting
    function `top_cis_overtime_heatmap()` to compute CIS_grubbs test
    over time

# ISAnalytics 1.7.2 (2022-05-23)

## BUG FIXES AND MINOR CHANGES

-   Fixed minor issues in `import_association_file()` - function had
    minor issues when importing \*.xlsx files and missing optional
    columns threw errors
-   Fixed bug in `as_sparse_matrix()` - function failed when trying to
    process an aggregated matrix

## NEW

-   Added 2 new utility functions `export_ISA_settings()` and
    `import_ISA_settings()` that allow a faster workflow setup

# ISAnalytics 1.7.1 (2022-05-04)

## BUG FIXES AND MINOR CHANGES

-   Fixed minor issue in `compute_near_integrations()` - function
    errored when `report_path` argument was set to `NULL`
-   Fixed dplyr warning in `integration_alluvial_plot()` internals
-   Fixed issue with report of VISPA2 stats - report failed due to minor
    error in rmd fragment
-   Internals of `remove_collisions()` use again dplyr internally for
    joining and grouping operations - needed because of performance
    issues with data.table
-   `fisher_scatterplot()` has 2 new arguments that allow the disabling
    of highlighting for some genes even if their p-value is under the
    threshold

# ISAnalytics 1.5.4 (2022-04-20)

## MAJOR CHANGES

-   ISAnalytics has now a new “dynamic vars system” to allow more
    flexibility on user inputs, view the dedicated vignette with
    `vignette("workflow_start", package="ISAnalytics")`
-   All package functions were reviewed to work properly with this
    system

## NEW FEATURES

-   `gene_frequency_fisher()` is a new function of the analysis family
    that allows the computation of Fisher’s exact test p-values on gene
    frequency - `fisher_scatterplot()` is the associated plotting
    function
-   `top_targeted_genes()` is a new function of the analysis family that
    produces the top n targeted genes based on the number of IS
-   `NGSdataExplorer()` is a newly implemented Shiny interface that
    allows the exploration and plotting of data
-   zipped examples were removed from the package to contain size. To
    compensate, the new function `generate_default_folder_structure()`
    generates the standard folder structure with package-included data
    on-demand
-   `transform_columns()` is a new utility function, also used
    internally by other exported functions, that allows arbitrary
    transformations on data frame columns

## MINOR CHANGES

-   `remove_collisions()` now has a dedicated parameter to specify how
    independent samples are identified
-   `compute_near_integration_sites()` now has a parameter called
    `additional_agg_lambda()` to allow aggregation of additional columns
-   `CIS_grubbs()` now signals if there are missing genes in the
    refgenes table and eventually returns them as a df
-   `outlier_filter()` is now able to take multiple tests in input and
    combine them with a given logic. It now also produces an HTML
    report.
-   Several functions now use data.table under the hood
-   Color of the strata containing IS below threshold can now be set in
    `integration_alluvial_plot()`

## BUG FIXES

-   Fixed a minor bug in `import_Vispa2_stats()` - function failed when
    passing `report_path = NULL`
-   Fixed minor issue in `circos_genomic_density()` when trying to use a
    pdf device

## DEPRECATED FUNCTIONS

-   `unzip_file_system()` was made defunct in favor of
    `generate_default_folder_structure()`
-   `cumulative_count_union()` was deprecated and its functionality was
    moved to `cumulative_is()`

# ISAnalytics 1.5.3 (2022-01-13)

## MINOR CHANGES

-   Added arguments `fragmentEstimate_column` and
    `fragmentEstimate_threshold` in `HSC_population_size_estimate()`.
    Slightly revised filtering logic.
-   Updated package logo and website

# ISAnalytics 1.5.2 (2021-12-14)

## NEW (MINOR)

-   Added function to check for annotation problems in IS matrices

## MINOR CHANGES

-   Added argument `max_workers` in function `remove_collisions()`
-   Updated default functions for `aggregate_metadata()`
-   Added annotation issues section in import matrices report

## FIXES

-   Fixed minor issue in internals for file system alignment checks
-   Fixed minor issue in internal call to `import_Vispa2_stats()` from
    `import_association_file()`
-   Added safe computation of sharing in `remove_collisions()`: if
    process fails function doesn’t stop

# ISAnalytics 1.5.1 (2021-10-28)

## FIXES

-   Attempt to fix issues with parallel computation on Windows for some
    plotting functions

# ISAnalytics 1.3.9 (2021-10-25)

## FIXES

-   Fixed issues with function that make use of BiocParallel that
    sometimes failed on Windows platform

# ISAnalytics 1.3.7 (2021-10-20)

## NEW

-   Added new feature `iss_source()`

## FIXES

-   Fixed minor issues in data files `refGenes_mm9` and function
    `compute_near_integrations()`

# ISAnalytics 1.3.6 (2021-10-05)

## NEW

-   Added new feature `purity_filter()`

## FIXES

-   Fixed small issue in printing information in reports

# ISAnalytics 1.3.5 (2021-09-21)

## MAJOR CHANGES

-   Reworked `is_sharing()` function, detailed usage in vignette
    `vignette("sharing_analyses", package = "ISAnalytics")`

## NEW

-   New function `cumulative_is()`
-   New function for plotting sharing as venn/euler diagrams
    `sharing_venn()`

# ISAnalytics 1.3.4 (2021-08 -03)

## FIXES/MINOR UPDATES

-   Fixed issue in tests that lead to broken build
-   Slightly modified included data set for better examples

# ISAnalytics 1.3.3 (2021-07-30)

## MAJOR CHANGES

-   Completely reworked interactive HTML report system, for details take
    a look at the new vignette
    `vignette("report_system", package = "ISAnalytics")`
-   Old `ISAnalytics.widgets` option has been replaced by
    `ISAnalytics.reports`
-   In `remove_collisions()`, removed arguments `seq_count_col`,
    `max_rows_reports` and `save_widget_path`, added arguments
    `quant_cols` and `report_path` (see documentation for details)

## MINOR CHANGES

-   `import_single_Vispa2Matrix()` now allows keeping additional
    non-standard columns
-   `compute_near_integrations()` is now faster on bigger data sets
-   Changed default values for arguments `columns` and `key` in
    `compute_abundance()`
-   `compute_near_integrations()` now produces only re-calibration map
    in \*.tsv format
-   `CIS_grubbs()` now supports calculations for each group specified in
    argument `by`
-   In `sample_statistics()` now there is the option to include the
    calculation of distinct integration sites for each group (if
    mandatory vars are present)

## NEW FUNCTIONALITY

-   Added new plotting function `circos_genomic_density()`

## FIXES

-   Fixed minor issue with NA values in alluvial plots

## DEPRECATIONS

-   `import_parallel_Vispa2Matrices_interactive()` and
    `import_parallel_Vispa2Matrices_auto()` are officially deprecated in
    favor of `import_parallel_Vispa2Matrices()`

## OTHER

-   The package has now a more complete and functional example data set
    for executable examples
-   Reworked documentation

# ISAnalytics 1.3.2 (2021-06-28)

## FIXES

-   Corrected issues in man pages

# ISAnalytics 1.3.1 (2021-06-24)

## NEW FUNCTIONALITY

-   `is_sharing` computes the sharing of IS between groups
-   `sharing_heatmap` allows visualization of sharing data through
    heatmaps
-   `integration_alluvial_plot` allows visualization of integration
    sites distribution in groups over time.
-   `top_abund_tableGrob` can be used in combination with the previous
    function or by itself to obtain a summary of top abundant
    integrations as an R graphic (tableGrob) object that can be combined
    with plots.

## MINOR UPDATES

-   Added more default stats functions to `default_stats`
-   Added optional automatic conversion of time points in months and
    years when importing association file
-   Minor fixes in `generate_Vispa2_launch_AF`

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

-   Fixed minor issues in import_association_file when checking
    parameters

# ISAnalytics 1.1.6 (2020-02-06)

## UPGRADES

-   It is now possible to save html reports to file from
    import_parallel_Vispa2Matrices_auto and
    import_parallel_Vispa2Matrices_interactive, remove_collisions and
    compute_near_integrations

## FIXES

-   Fixed sample_statistics: now functions that have data frame output
    do not produce nested tables. Flat tables are ready to be saved to
    file or can be nested.
-   Simplified association file check logic in remove_collisions: now
    function blocks only if the af doesn’t contain the needed columns

# ISAnalytics 1.1.5 (2020-02-03)

## UPGRADES

-   Upgraded import_association_file function: now file alignment is not
    mandatory anymore and it is possible to save the html report to file
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
-   Minor fix in import_association_file file function: added multiple
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
-   Added first analysis function (compute_abundance)

## SIGNIFICANT USER-VISIBLE CHANGES

-   Dropped structure `ISADataFrame`: now the package only uses standard
    tibbles
-   Modified package documentation

# ISAnalytics 0.99.8 (2020-08-12)

-   Submitted to Bioconductor

</details>

# Getting help

For help please contact the maintainer of the package or open an issue
on GitHub.
