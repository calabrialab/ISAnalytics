# Getting started with ISAnalytics

## Introduction

ISAnalytics is an R package developed to analyze gene therapy vector
insertion sites data identified from genomics next generation sequencing
reads for clonal tracking studies.

## Installation and options

`ISAnalytics` can be installed quickly in different ways:

- You can install it via [Bioconductor](http://bioconductor.org)
- You can install it via GitHub using the package `devtools`

There are always 2 versions of the package active:

- `RELEASE` is the latest stable version
- `DEVEL` is the development version, it is the most up-to-date version
  where all new features are introduced

### Installation from bioconductor

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

### Installation from GitHub

RELEASE:

``` r

if (!require(devtools)) {
  install.packages("devtools")
}
devtools::install_github("calabrialab/ISAnalytics",
                         ref = "RELEASE_3_21",
                         dependencies = TRUE,
                         build_vignettes = TRUE)
```

DEVEL:

``` r

if (!require(devtools)) {
  install.packages("devtools")
}
devtools::install_github("calabrialab/ISAnalytics",
                         ref = "devel",
                         dependencies = TRUE,
                         build_vignettes = TRUE)
```

### Setting options

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

## Setting up the workflow

In the newer version of ISAnalytics, we introduced a “dynamic variables
system”, to allow more flexibility in terms of input formats. Before
starting with the analysis workflow, you can specify how your inputs are
structured so that the package can process them. For more information on
how to do this take a look at
[`vignette("workflow_start", package = "ISAnalytics")`](https://calabrialab.github.io/ISAnalytics/dev/articles/workflow_start.md).

## The first steps

The first steps of the analysis workflow involve the import and parsing
of data and metadata files from disk.

- Import metadata with
  [`import_association_file()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_association_file.md)
  and/or
  [`import_Vispa2_stats()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_Vispa2_stats.md)
- Import data with
  [`import_single_Vispa2Matrix()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_single_Vispa2Matrix.md)
  or
  [`import_parallel_Vispa2Matrices()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_parallel_Vispa2Matrices.md)

Refer to the vignette
[`vignette("workflow_start", package = "ISAnalytics")`](https://calabrialab.github.io/ISAnalytics/dev/articles/workflow_start.md)
for more details.

## Data cleaning and pre-processing

ISAnalytics offers several different functions for cleaning and
pre-processing your data.

- Recalibration: identifies integration events that are near to each
  other and condenses them into a single event whenever appropriate -
  [`compute_near_integrations()`](https://calabrialab.github.io/ISAnalytics/dev/reference/compute_near_integrations.md)
- Outliers identification and removal: identifies samples that are
  considered outliers according to user-defined logic and filters them
  out -
  [`outlier_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/outlier_filter.md)
- Collision removal: identifies collision events between independent
  samples -
  [`remove_collisions()`](https://calabrialab.github.io/ISAnalytics/dev/reference/remove_collisions.md),
  see also the dedicated vignette
  [`vignette("workflow_start", package = "ISAnalytics")`](https://calabrialab.github.io/ISAnalytics/dev/articles/workflow_start.md)
- Filter based on cell lineage purity: identifies and removes
  contamination between different cell types -
  [`purity_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/purity_filter.md)
- Data and metadata aggregation: allows the union of biological samples
  from single pcr replicates or other arbitrary aggregations -
  [`aggregate_values_by_key()`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_values_by_key.md),
  [`aggregate_metadata()`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_metadata.md),
  see also the dedicated vignette
  [`vignette("workflow_start", package = "ISAnalytics")`](https://calabrialab.github.io/ISAnalytics/dev/articles/workflow_start.md)

## Answering biological questions

You can answer very different biological questions by using the provided
functions with appropriate inputs.

- Descriptive statistics:
  [`sample_statistics()`](https://calabrialab.github.io/ISAnalytics/dev/reference/sample_statistics.md)
- IS relative abundance:
  [`compute_abundance()`](https://calabrialab.github.io/ISAnalytics/dev/reference/compute_abundance.md),
  [`integration_alluvial_plot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/integration_alluvial_plot.md)
- Top abundant IS:
  [`top_integrations()`](https://calabrialab.github.io/ISAnalytics/dev/reference/top_integrations.md)
- Top targeted genes:
  [`top_targeted_genes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/top_targeted_genes.md)
- Grubbs test for common insertion sites (CIS):
  [`CIS_grubbs()`](https://calabrialab.github.io/ISAnalytics/dev/reference/CIS_grubbs.md),
  [`CIS_volcano_plot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/CIS_volcano_plot.md)
- Fisher’s exact test for gene frequency and IS distribution on target
  genome:
  [`gene_frequency_fisher()`](https://calabrialab.github.io/ISAnalytics/dev/reference/gene_frequency_fisher.md),
  [`fisher_scatterplot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/fisher_scatterplot.md),
  [`circos_genomic_density()`](https://calabrialab.github.io/ISAnalytics/dev/reference/circos_genomic_density.md)
- Clonal sharing analyses:
  [`is_sharing()`](https://calabrialab.github.io/ISAnalytics/dev/reference/is_sharing.md),
  [`iss_source()`](https://calabrialab.github.io/ISAnalytics/dev/reference/iss_source.md),
  [`sharing_heatmap()`](https://calabrialab.github.io/ISAnalytics/dev/reference/sharing_heatmap.md),
  [`sharing_venn()`](https://calabrialab.github.io/ISAnalytics/dev/reference/sharing_venn.md)
- Estimate HSPCs population size:
  [`HSC_population_size_estimate()`](https://calabrialab.github.io/ISAnalytics/dev/reference/HSC_population_size_estimate.md),
  [`HSC_population_plot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/HSC_population_plot.md)

For more, please refer to the full function reference.

## Working with other kinds of data

ISAnalytics is designed to be flexible concerning input formats, thus it
is suited to process various kinds of data provided the correct dynamic
configuration is set.

We demonstrate this with an example that uses barcodes data. The matrix
is publicly available
[here](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE144340&format=file&file=GSE144340%5FMatrix%5F542%2Etsv%2Egz)
(Ferrari Samuele Jacob Aurelien, 2020), metadata was provided to us by
the authors and it is available in the package additional files.

``` r

library(ISAnalytics)

# Set appropriate data and metadata specs ----
metadata_specs <- tibble::tribble(
    ~names, ~types, ~transform, ~flag, ~tag,
    "ProjectID", "char", NULL, "required", "project_id",
    "SubjectID", "char", NULL, "required", "subject",
    "Tissue", "char", NULL, "required", "tissue",
    "TimePoint", "int", NULL, "required", "tp_days",
    "CellMarker", "char", NULL, "required", "cell_marker",
    "ID", "char", NULL, "required", "pcr_repl_id",
    "SourceFileName", "char", NULL, "optional", NA_character_,
    "Link", "char", NULL, "optional", NA_character_
)
set_af_columns_def(metadata_specs)
#> Warning: Warning: important tags missing
#> ℹ Some tags are required for proper execution of some functions. If these tags are not provided, execution of dependent functions might fail. Review your inputs carefully.
#> ℹ Missing tags: pool_id, fusion_id, tag_seq, vector_id, tag_id, pcr_replicate, vispa_concatenate, proj_folder
#> ℹ To see where these are involved type `inspect_tags(c('pool_id','fusion_id','tag_seq','vector_id','tag_id','pcr_replicate','vispa_concatenate','proj_folder'))`
#> Association file columns specs successfully changed

mandatory_specs <- tibble::tribble(
    ~names, ~types, ~transform, ~flag, ~tag,
    "BarcodeSeq", "char", NULL, "required", NA_character_
)
set_mandatory_IS_vars(mandatory_specs)
#> Warning: Warning: important tags missing
#> ℹ Some tags are required for proper execution of some functions. If these tags are not provided, execution of dependent functions might fail. Review your inputs carefully.
#> ℹ Missing tags: chromosome, locus, is_strand
#> ℹ To see where these are involved type `inspect_tags(c('chromosome','locus','is_strand'))`
#> Mandatory IS vars successfully changed

# Files ----
data_folder <- tempdir()
utils::unzip(zipfile = system.file("testdata", 
                                   "testdata.zip",
                                   package = "ISAnalytics"), 
             exdir = data_folder, overwrite = TRUE)
meta_file <- "barcodes_example_af.tsv.xz"
matrix_file <- "GSE144340_Matrix_542.tsv.xz"

# Data import ----
af <- import_association_file(fs::path(data_folder, meta_file),
    report_path = NULL
)
af
#> # A tibble: 15 × 10
#>    ProjectID    SubjectID Tissue TimePoint CellMarker ID    SourceFileName Link 
#>    <chr>        <chr>     <chr>      <int> <chr>      <chr> <chr>          <chr>
#>  1 PMID32601433 A0        BM            21 Whole      BM_A0 GSE144340_Mat… http…
#>  2 PMID32601433 A0        PB            21 Whole      PB21… GSE144340_Mat… http…
#>  3 PMID32601433 A1        BM            21 Whole      BM_A1 GSE144340_Mat… http…
#>  4 PMID32601433 A1        PB            21 Whole      PB21… GSE144340_Mat… http…
#>  5 PMID32601433 A2        PB            21 Whole      PB21… GSE144340_Mat… http…
#>  6 PMID32601433 A3        PB            21 Whole      PB21… GSE144340_Mat… http…
#>  7 PMID32601433 A4        BM            21 Whole      BM_A4 GSE144340_Mat… http…
#>  8 PMID32601433 A4        PB            21 Whole      PB21… GSE144340_Mat… http…
#>  9 PMID32601433 C0        PB            21 Whole      PB21… GSE144340_Mat… http…
#> 10 PMID32601433 C1        BM            21 Whole      BM_C1 GSE144340_Mat… http…
#> 11 PMID32601433 C1        PB            21 Whole      PB21… GSE144340_Mat… http…
#> 12 PMID32601433 C2        BM            21 Whole      BM_C2 GSE144340_Mat… http…
#> 13 PMID32601433 C2        PB            21 Whole      PB21… GSE144340_Mat… http…
#> 14 PMID32601433 C3        BM            21 Whole      BM_C3 GSE144340_Mat… http…
#> 15 PMID32601433 C3        PB            21 Whole      PB21… GSE144340_Mat… http…
#> # ℹ 2 more variables: TimepointMonths <chr>, TimepointYears <chr>

matrix <- import_single_Vispa2Matrix(fs::path(data_folder, matrix_file),
    sample_names_to = "ID"
)
#> Warning: compression format not supported by fread
#> ℹ File will be read using readr
#> Reading file...
#> ℹ Mode: classic
#> Reshaping...
#> *** File info *** 
#> • --- Annotated: FALSE
#> • --- Dimensions: 31757 x 16
#> • --- Read mode: classic
#> • --- Sample count: 15
matrix
#> # A tibble: 33,776 × 3
#>    BarcodeSeq             ID      Value
#>    <chr>                  <chr>   <dbl>
#>  1 AAAAAAAACACGGAGAACGACG PB21_C3     2
#>  2 AAAAAAAACGCGAACAACTACG PB21_C3     1
#>  3 AAAAAAAACTCAAAAAAGAAAT BM_C3       1
#>  4 AAAAAAAATTTACACAAAGAAA BM_A4       1
#>  5 AAAAAAAATTTTTAAACGTACC BM_A0       1
#>  6 AAAAAACATATCTATAGTTACC BM_A0       1
#>  7 AAAAAAGACGACGATAGGCACG PB21_C1     1
#>  8 AAAAAAGACGTTTATAGGTGTA PB21_A2     1
#>  9 AAAAAAGACTGCGACAAAAGGG BM_A4       1
#> 10 AAAAAAGACTTTGATAACCACG PB21_C3     1
#> # ℹ 33,766 more rows

# Descriptive stats ----
desc_stats <- sample_statistics(matrix, af,
    sample_key = pcr_id_column(),
    value_columns = "Value"
)$metadata |>
    dplyr::rename(distinct_barcodes = "nIS")
desc_stats
#> # A tibble: 15 × 29
#>    ProjectID    SubjectID Tissue TimePoint CellMarker ID    SourceFileName Link 
#>    <chr>        <chr>     <chr>      <int> <chr>      <chr> <chr>          <chr>
#>  1 PMID32601433 A0        BM            21 Whole      BM_A0 GSE144340_Mat… http…
#>  2 PMID32601433 A0        PB            21 Whole      PB21… GSE144340_Mat… http…
#>  3 PMID32601433 A1        BM            21 Whole      BM_A1 GSE144340_Mat… http…
#>  4 PMID32601433 A1        PB            21 Whole      PB21… GSE144340_Mat… http…
#>  5 PMID32601433 A2        PB            21 Whole      PB21… GSE144340_Mat… http…
#>  6 PMID32601433 A3        PB            21 Whole      PB21… GSE144340_Mat… http…
#>  7 PMID32601433 A4        BM            21 Whole      BM_A4 GSE144340_Mat… http…
#>  8 PMID32601433 A4        PB            21 Whole      PB21… GSE144340_Mat… http…
#>  9 PMID32601433 C0        PB            21 Whole      PB21… GSE144340_Mat… http…
#> 10 PMID32601433 C1        BM            21 Whole      BM_C1 GSE144340_Mat… http…
#> 11 PMID32601433 C1        PB            21 Whole      PB21… GSE144340_Mat… http…
#> 12 PMID32601433 C2        BM            21 Whole      BM_C2 GSE144340_Mat… http…
#> 13 PMID32601433 C2        PB            21 Whole      PB21… GSE144340_Mat… http…
#> 14 PMID32601433 C3        BM            21 Whole      BM_C3 GSE144340_Mat… http…
#> 15 PMID32601433 C3        PB            21 Whole      PB21… GSE144340_Mat… http…
#> # ℹ 21 more variables: TimepointMonths <chr>, TimepointYears <chr>,
#> #   Value_sum <dbl>, Value_count <int>, Value_shannon <dbl>,
#> #   Value_simpson <dbl>, Value_invsimpson <dbl>, Value_describe_vars <dbl>,
#> #   Value_describe_n <dbl>, Value_describe_mean <dbl>, Value_describe_sd <dbl>,
#> #   Value_describe_median <dbl>, Value_describe_trimmed <dbl>,
#> #   Value_describe_mad <dbl>, Value_describe_min <dbl>,
#> #   Value_describe_max <dbl>, Value_describe_range <dbl>, …

# Aggregation and new stats ----
agg_key <- c("SubjectID")
agg <- aggregate_values_by_key(matrix, af,
    key = agg_key,
    group = "BarcodeSeq",
    join_af_by = pcr_id_column()
)
agg
#> # A tibble: 33,267 × 3
#>    BarcodeSeq             SubjectID Value_sum
#>    <chr>                  <chr>         <dbl>
#>  1 AAAAAAAACACGGAGAACGACG C3                2
#>  2 AAAAAAAACGCGAACAACTACG C3                1
#>  3 AAAAAAAACTCAAAAAAGAAAT C3                1
#>  4 AAAAAAAATTTACACAAAGAAA A4                1
#>  5 AAAAAAAATTTTTAAACGTACC A0                1
#>  6 AAAAAACATATCTATAGTTACC A0                1
#>  7 AAAAAAGACGACGATAGGCACG C1                1
#>  8 AAAAAAGACGTTTATAGGTGTA A2                1
#>  9 AAAAAAGACTGCGACAAAAGGG A4                1
#> 10 AAAAAAGACTTTGATAACCACG C3                1
#> # ℹ 33,257 more rows

agg_meta_functions <- tibble::tribble(
    ~Column, ~Function, ~Args, ~Output_colname,
    "TimePoint", ~ mean(.x, na.rm = TRUE), NA, "{.col}_avg",
    "CellMarker", ~ length(unique(.x)), NA, "distinct_cell_marker_count",
    "ID", ~ length(unique(.x)), NA, "distinct_id_count"
)
agg_meta <- aggregate_metadata(
    af,
    aggregating_functions = agg_meta_functions,
    grouping_keys = agg_key
)
agg_meta
#> # A tibble: 9 × 4
#>   SubjectID TimePoint_avg distinct_cell_marker_count distinct_id_count
#>   <chr>             <dbl>                      <int>             <int>
#> 1 A0                   21                          1                 2
#> 2 A1                   21                          1                 2
#> 3 A2                   21                          1                 1
#> 4 A3                   21                          1                 1
#> 5 A4                   21                          1                 2
#> 6 C0                   21                          1                 1
#> 7 C1                   21                          1                 2
#> 8 C2                   21                          1                 2
#> 9 C3                   21                          1                 2

agg_stats <- sample_statistics(agg, agg_meta,
    sample_key = agg_key,
    value_columns = "Value_sum"
)$metadata |>
    dplyr::rename(distinct_barcodes = "nIS")
agg_stats
#> # A tibble: 9 × 23
#>   SubjectID TimePoint_avg distinct_cell_marker…¹ distinct_id_count Value_sum_sum
#>   <chr>             <dbl>                  <int>             <int>         <dbl>
#> 1 A0                   21                      1                 2        326467
#> 2 A1                   21                      1                 2        378987
#> 3 A2                   21                      1                 1        124676
#> 4 A3                   21                      1                 1        180497
#> 5 A4                   21                      1                 2        473256
#> 6 C0                   21                      1                 1         59966
#> 7 C1                   21                      1                 2        399316
#> 8 C2                   21                      1                 2        492590
#> 9 C3                   21                      1                 2        341166
#> # ℹ abbreviated name: ¹​distinct_cell_marker_count
#> # ℹ 18 more variables: Value_sum_count <int>, Value_sum_shannon <dbl>,
#> #   Value_sum_simpson <dbl>, Value_sum_invsimpson <dbl>,
#> #   Value_sum_describe_vars <dbl>, Value_sum_describe_n <dbl>,
#> #   Value_sum_describe_mean <dbl>, Value_sum_describe_sd <dbl>,
#> #   Value_sum_describe_median <dbl>, Value_sum_describe_trimmed <dbl>,
#> #   Value_sum_describe_mad <dbl>, Value_sum_describe_min <dbl>, …

# Abundance ----
abundance <- compute_abundance(agg, columns = "Value_sum", key = agg_key)
abundance
#> # A tibble: 33,267 × 5
#>    BarcodeSeq  SubjectID Value_sum Value_sum_RelAbundance Value_sum_PercAbunda…¹
#>    <chr>       <chr>         <dbl>                  <dbl>                  <dbl>
#>  1 AAAAAAAACA… C3                2             0.00000586               0.000586
#>  2 AAAAAAAACG… C3                1             0.00000293               0.000293
#>  3 AAAAAAAACT… C3                1             0.00000293               0.000293
#>  4 AAAAAAAATT… A4                1             0.00000211               0.000211
#>  5 AAAAAAAATT… A0                1             0.00000306               0.000306
#>  6 AAAAAACATA… A0                1             0.00000306               0.000306
#>  7 AAAAAAGACG… C1                1             0.00000250               0.000250
#>  8 AAAAAAGACG… A2                1             0.00000802               0.000802
#>  9 AAAAAAGACT… A4                1             0.00000211               0.000211
#> 10 AAAAAAGACT… C3                1             0.00000293               0.000293
#> # ℹ 33,257 more rows
#> # ℹ abbreviated name: ¹​Value_sum_PercAbundance

reset_dyn_vars_config()
#> Mandatory IS vars reset to default
#> Annotation IS vars reset to default
#> Association file columns specs reset to default
#> ISS stats specs reset to default
#> Matrix suffixes specs reset to default
```

## Using the Shiny interface

The package provides a simple Shiny interface for data exploration and
plotting. To start the interface use:

``` r

NGSdataExplorer()
```

The application main page will show a loading screen for a file. It is
possible to load files also from the R environment, for example, before
opening the app, we can load the included association file:

``` r

data("association_file")
```

Once in the application we can choose `"association_file"` from the R
environment loading option screen and click on “Import data”. Once data
is imported, we can click on the “Explore” tab in the upper navbar: here
we will see 2 tabs, one allows interactive exploration of data in
tabular form, in the other tab we can plot data. It is possible to
customize several different parameters for the plot and finally save it
to file with the dedicated button at the end of the page.

The Shiny interface is still currently under active development and new
features will be added in the near future.

## Ensuring reproducibility of results

Several implemented functions produce static HTML reports that can be
saved on disk, or tabular files. Reports contain the relevant
information on how the function was called, inputs and outputs
statistics, and session info for reproducibility.

## Browse documentation online and keep updated

ISAnalytics has it’s dedicated package website where you can browse the
documentation and vignettes easily, in addition to keeping up to date
with all relevant updates. Visit the website at
<https://calabrialab.github.io/ISAnalytics/>

## Problems?

If you have any issues the documentation can’t solve, get in touch by
opening an issue on GitHub or contacting the maintainers

## Reproducibility

The following session information was included to ensure reproducibility
of results and environment tracking, as recommended by Bioconductor
guidelines.

    #> R version 4.6.1 (2026-06-24)
    #> Platform: x86_64-pc-linux-gnu
    #> Running under: Ubuntu 24.04.4 LTS
    #> 
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    #> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    #> 
    #> locale:
    #>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    #>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    #>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    #> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    #> 
    #> time zone: UTC
    #> tzcode source: system (glibc)
    #> 
    #> attached base packages:
    #> [1] stats     graphics  grDevices utils     datasets  methods   base     
    #> 
    #> other attached packages:
    #> [1] ISAnalytics_1.23.0 RefManageR_1.4.0   BiocStyle_2.40.0  
    #> 
    #> loaded via a namespace (and not attached):
    #>   [1] mnormt_2.1.2        writexl_1.5.4       permute_0.9-10     
    #>   [4] readxl_1.5.0        datamods_1.5.3      rlang_1.3.0        
    #>   [7] magrittr_2.0.5      rio_1.3.0           otel_0.2.0         
    #>  [10] e1071_1.7-17        compiler_4.6.1      mgcv_1.9-4         
    #>  [13] systemfonts_1.3.2   vctrs_0.7.3         stringr_1.6.0      
    #>  [16] pkgconfig_2.0.3     crayon_1.5.3        fastmap_1.2.0      
    #>  [19] backports_1.5.1     utf8_1.2.6          promises_1.5.0     
    #>  [22] phosphoricons_0.2.1 rmarkdown_2.31      tzdb_0.5.0         
    #>  [25] ragg_1.5.2          purrr_1.2.2         bit_4.6.0          
    #>  [28] xfun_0.60           cachem_1.1.0        jsonlite_2.0.0     
    #>  [31] later_1.4.8         BiocParallel_1.46.0 psych_2.6.5        
    #>  [34] cluster_2.1.8.2     parallel_4.6.1      R6_2.6.1           
    #>  [37] bslib_0.11.0        stringi_1.8.7       RColorBrewer_1.1-3 
    #>  [40] shinybusy_0.3.3     parallelly_1.48.0   lubridate_1.9.5    
    #>  [43] jquerylib_0.1.4     cellranger_1.1.0    Rcpp_1.1.2         
    #>  [46] bookdown_0.47       iterators_1.0.14    knitr_1.51         
    #>  [49] future.apply_1.20.2 readr_2.2.0         Matrix_1.7-5       
    #>  [52] splines_4.6.1       httpuv_1.6.17       timechange_0.4.0   
    #>  [55] tidyselect_1.2.1    yaml_2.3.12         vegan_2.7-5        
    #>  [58] codetools_0.2-20    listenv_1.0.0       lattice_0.22-9     
    #>  [61] tibble_3.3.1        plyr_1.8.9          shiny_1.14.0       
    #>  [64] withr_3.0.3         S7_0.2.2            evaluate_1.0.5     
    #>  [67] future_1.70.0       desc_1.4.3          proxy_0.4-29       
    #>  [70] xml2_1.6.0          pillar_1.11.1       BiocManager_1.30.27
    #>  [73] KernSmooth_2.23-26  foreach_1.5.2       generics_0.1.4     
    #>  [76] vroom_1.7.1         hms_1.1.4           ggplot2_4.0.3      
    #>  [79] scales_1.4.0        globals_0.19.1      xtable_1.8-8       
    #>  [82] class_7.3-23        glue_1.8.1          toastui_0.4.0      
    #>  [85] tools_4.6.1         data.table_1.18.4   reactable_0.4.5    
    #>  [88] fs_2.1.0            grid_4.6.1          tidyr_1.3.2        
    #>  [91] bibtex_0.5.2        nlme_3.1-169        cli_3.6.6          
    #>  [94] textshaping_1.0.5   doFuture_1.2.2      dplyr_1.2.1        
    #>  [97] gtable_0.3.6        sass_0.4.10         digest_0.6.39      
    #> [100] classInt_0.4-11     htmlwidgets_1.6.4   farver_2.1.2       
    #> [103] htmltools_0.5.9     pkgdown_2.2.1       lifecycle_1.0.5    
    #> [106] httr_1.4.8          shinyWidgets_0.9.1  mime_0.13          
    #> [109] MASS_7.3-65         bit64_4.8.2

## Bibliography

[\[1\]](#cite-ferrarisamueleefficient) B. S. Ferrari Samuele Jacob
Aurelien. “Efficient gene editing of human long-term hematopoietic stem
cells validated by clonal tracking”. In: *Nat Biotechnol 38, 1298–1308*
(Nov. 2020). DOI:
[https://doi.org/10.1038/s41587-020-0551-y](https://doi.org/https://doi.org/10.1038/s41587-020-0551-y).
