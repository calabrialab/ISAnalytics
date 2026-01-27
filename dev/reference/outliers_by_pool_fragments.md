# Identify and flag outliers based on pool fragments.

**\[stable\]** Identify and flag outliers based on expected number of
raw reads per pool.

## Usage

``` r
outliers_by_pool_fragments(
  metadata,
  key = "BARCODE_MUX",
  outlier_p_value_threshold = 0.01,
  normality_test = FALSE,
  normality_p_value_threshold = 0.05,
  transform_log2 = TRUE,
  per_pool_test = TRUE,
  pool_col = "PoolID",
  min_samples_per_pool = 5,
  flag_logic = "AND",
  keep_calc_cols = TRUE,
  report_path = default_report_path()
)
```

## Arguments

- metadata:

  The metadata data frame

- key:

  A character vector of numeric column names

- outlier_p_value_threshold:

  The p value threshold for a read to be considered an outlier

- normality_test:

  Perform normality test? Normality is assessed for each column in the
  key using Shapiro-Wilk test and if the values do not follow a normal
  distribution, other calculations are skipped

- normality_p_value_threshold:

  Normality threshold

- transform_log2:

  Perform a log2 trasformation on values prior the actual calculations?

- per_pool_test:

  Perform the test for each pool?

- pool_col:

  A character vector of the names of the columns that uniquely identify
  a pool

- min_samples_per_pool:

  The minimum number of samples that a pool needs to contain in order to
  be processed - relevant only if `per_pool_test = TRUE`

- flag_logic:

  A character vector of logic operators to obtain a global flag
  formula - only relevant if the key is longer than one. All operators
  must be chosen between: AND, OR, XOR, NAND, NOR, XNOR

- keep_calc_cols:

  Keep the calculation columns in the output data frame?

- report_path:

  The path where the report file should be saved. Can be a folder, a
  file or NULL if no report should be produced. Defaults to
  `{user_home}/ISAnalytics_reports`.

## Value

A data frame of metadata with the column `to_remove`

## Details

### Modular structure

The outlier filtering functions are structured in a modular fashion.
There are 2 kind of functions:

- Outlier tests - Functions that perform some kind of calculation based
  on inputs and flags metadata

- Outlier filter - A function that takes one or more outlier tests,
  combines all the flags with a given logic and filters out rows that
  are flagged as outliers

This function is an outlier test, and calculates for each column in the
key

- The zscore of the values

- The tstudent of the values

- The the associated p-value (tdist)

Optionally the test can be performed for each pool and a normality test
can be run prior the actual calculations. Samples are flagged if this
condition is respected:

- tdist \< outlier_p_value_threshold & zscore \< 0

If the key contains more than one column an additional flag logic can be
specified for combining the results. Example: let's suppose the key
contains the names of two columns, X and Y `key = c("X", "Y")` if we
specify the the argument `flag_logic = "AND"` then the reads will be
flagged based on this global condition: (tdist_X \<
outlier_p_value_threshold & zscore_X \< 0) AND (tdist_Y \<
outlier_p_value_threshold & zscore_Y \< 0)

The user can specify one or more logical operators that will be applied
in sequence.

## See also

Other Data cleaning and pre-processing:
[`aggregate_metadata()`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_metadata.md),
[`aggregate_values_by_key()`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_values_by_key.md),
[`compute_near_integrations()`](https://calabrialab.github.io/ISAnalytics/dev/reference/compute_near_integrations.md),
[`default_meta_agg()`](https://calabrialab.github.io/ISAnalytics/dev/reference/default_meta_agg.md),
[`outlier_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/outlier_filter.md),
[`purity_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/purity_filter.md),
[`realign_after_collisions()`](https://calabrialab.github.io/ISAnalytics/dev/reference/realign_after_collisions.md),
[`remove_collisions()`](https://calabrialab.github.io/ISAnalytics/dev/reference/remove_collisions.md),
[`threshold_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/threshold_filter.md)

## Examples

``` r
data("association_file", package = "ISAnalytics")
flagged <- outliers_by_pool_fragments(association_file,
    report_path = NULL
)
#> Removing NAs from data...
#> Log2 transformation, removing values <= 0
head(flagged)
#> # A tibble: 6 × 91
#>   ProjectID FUSIONID  PoolID TagSequence SubjectID VectorType VectorID
#>   <chr>     <chr>     <chr>  <chr>       <chr>     <chr>      <chr>   
#> 1 PJ01      ET#382.46 POOL01 LTR75LC38   PT001     lenti      GLOBE   
#> 2 PJ01      ET#381.40 POOL01 LTR53LC32   PT001     lenti      GLOBE   
#> 3 PJ01      ET#381.9  POOL01 LTR83LC66   PT001     lenti      GLOBE   
#> 4 PJ01      ET#381.71 POOL01 LTR27LC94   PT001     lenti      GLOBE   
#> 5 PJ01      ET#381.2  POOL01 LTR69LC52   PT001     lenti      GLOBE   
#> 6 PJ01      ET#382.28 POOL01 LTR37LC2    PT001     lenti      GLOBE   
#> # ℹ 84 more variables: ExperimentID <chr>, Tissue <chr>, TimePoint <chr>,
#> #   DNAFragmentation <chr>, PCRMethod <chr>, TagIDextended <chr>,
#> #   Keywords <chr>, CellMarker <chr>, TagID <chr>, NGSProvider <chr>,
#> #   NGSTechnology <chr>, ConverrtedFilesDir <chr>, ConverrtedFilesName <chr>,
#> #   SourceFileFolder <chr>, SourceFileNameR1 <chr>, SourceFileNameR2 <chr>,
#> #   DNAnumber <chr>, ReplicateNumber <int>, DNAextractionDate <date>,
#> #   DNAngUsed <dbl>, LinearPCRID <chr>, LinearPCRDate <date>, …
```
