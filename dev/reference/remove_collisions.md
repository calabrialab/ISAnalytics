# Identifies and removes collisions

**\[stable\]** A collision is an integration (aka a unique combination
of the provided
[`mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md))
which is observed in more than one independent sample. The function
tries to decide to which independent sample should an integration event
be assigned to, and if no decision can be taken, the integration is
completely removed from the data frame. For more details refer to the
vignette "Collision removal functionality":
[`vignette("workflow_start", package = "ISAnalytics")`](https://calabrialab.github.io/ISAnalytics/dev/articles/workflow_start.md)

## Usage

``` r
remove_collisions(
  x,
  association_file,
  independent_sample_id = c("ProjectID", "SubjectID"),
  date_col = "SequencingDate",
  reads_ratio = 10,
  quant_cols = c(seqCount = "seqCount", fragmentEstimate = "fragmentEstimate"),
  fold_threshold = 10,
  report_path = default_report_path(),
  max_workers = NULL
)
```

## Arguments

- x:

  Either a multi-quantification matrix (recommended) or a named list of
  matrices (names must be quantification types)

- association_file:

  The association file imported via
  [`import_association_file()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_association_file.md)

- independent_sample_id:

  A character vector of column names that identify independent samples

- date_col:

  The date column that should be considered.

- reads_ratio:

  Deprecated alias for `fold_threshold`, kept for backward
  compatibility. If both parameters are supplied they must have the same
  value.

- quant_cols:

  A named character vector where names are quantification types and
  values are the names of the corresponding columns. The quantification
  `seqCount` MUST be included in the vector.

- fold_threshold:

  A single numeric value greater than 1. For each collision, the
  sequence count values are summed within independent samples and
  compared to the maximum summed sequence count observed for the same
  integration. Observations from independent samples with
  `max(seqCount) / seqCount >= fold_threshold` are removed before
  applying the temporal rule. The default is 10, and a ratio exactly
  equal to the threshold is considered sufficient for removal. If more
  than one independent sample is not clearly separated by this fold
  rule, only those remaining observations are passed to the temporal
  rule. Zero sequence counts are allowed: if the maximum abundance is
  positive, zero-abundance observations are removed by the fold rule; if
  all abundances are zero, the fold rule is unresolved and the temporal
  rule is used. Missing, non-finite or negative sequence counts are
  rejected.

- report_path:

  The path where the report file should be saved. Can be a folder or
  `NULL` if no report should be produced. Defaults to
  `{user_home}/ISAnalytics_reports`.

- max_workers:

  Maximum number of parallel workers to distribute the workload. If
  `NULL` (default) produces the maximum amount of workers allowed, a
  numeric value is requested otherwise. WARNING: a higher number of
  workers speeds up computation at the cost of memory consumption! Tune
  this parameter accordingly.

## Value

Either a multi-quantification matrix or a list of data frames

## Required tags

The function will explicitly check for the presence of these tags:

- project_id

- pool_id

- pcr_replicate

## See also

Other Data cleaning and pre-processing:
[`aggregate_metadata()`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_metadata.md),
[`aggregate_values_by_key()`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_values_by_key.md),
[`compute_near_integrations()`](https://calabrialab.github.io/ISAnalytics/dev/reference/compute_near_integrations.md),
[`default_meta_agg()`](https://calabrialab.github.io/ISAnalytics/dev/reference/default_meta_agg.md),
[`outlier_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/outlier_filter.md),
[`outliers_by_pool_fragments()`](https://calabrialab.github.io/ISAnalytics/dev/reference/outliers_by_pool_fragments.md),
[`purity_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/purity_filter.md),
[`realign_after_collisions()`](https://calabrialab.github.io/ISAnalytics/dev/reference/realign_after_collisions.md),
[`threshold_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/threshold_filter.md)

## Examples

``` r
data("integration_matrices", package = "ISAnalytics")
data("association_file", package = "ISAnalytics")
no_coll <- remove_collisions(
    x = integration_matrices,
    association_file = association_file,
    fold_threshold = 10,
    report_path = NULL
)
#> Identifying collisions...
#> Processing collisions...
#> Finished!
head(no_coll)
#> # A tibble: 6 × 8
#>   chr   integration_locus strand GeneName GeneStrand CompleteAmplificationID    
#>   <chr>             <dbl> <chr>  <chr>    <chr>      <chr>                      
#> 1 1              16602483 +      FBXO42   -          PJ01_POOL01_LTR83LC46_PT00…
#> 2 1              16602483 +      FBXO42   -          PJ01_POOL01_LTR37LC2_PT001…
#> 3 1              16602483 +      FBXO42   -          PJ01_POOL01_LTR85LC54_PT00…
#> 4 1              26446899 +      PDIK1L   +          PJ01_POOL01_LTR85LC54_PT00…
#> 5 1              26446899 +      PDIK1L   +          PJ01_POOL01_LTR83LC46_PT00…
#> 6 1              26446899 +      PDIK1L   +          PJ01_POOL01_LTR69LC52_PT00…
#> # ℹ 2 more variables: seqCount <dbl>, fragmentEstimate <dbl>
```
