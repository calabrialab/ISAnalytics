# Performs aggregation on metadata contained in the association file.

**\[stable\]** Groups metadata by the specified grouping keys and
returns a summary of info for each group. For more details on how to use
this function:
[`vignette("workflow_start", package = "ISAnalytics")`](https://calabrialab.github.io/ISAnalytics/dev/articles/workflow_start.md)

## Usage

``` r
aggregate_metadata(
  association_file,
  grouping_keys = c("SubjectID", "CellMarker", "Tissue", "TimePoint"),
  aggregating_functions = default_meta_agg(),
  import_stats = lifecycle::deprecated()
)
```

## Arguments

- association_file:

  The imported association file (via
  [import_association_file](https://calabrialab.github.io/ISAnalytics/dev/reference/import_association_file.md))

- grouping_keys:

  A character vector of column names to form a grouping operation

- aggregating_functions:

  A data frame containing specifications of the functions to be applied
  to columns in the association file during aggregation. It defaults to
  [default_meta_agg](https://calabrialab.github.io/ISAnalytics/dev/reference/default_meta_agg.md).
  The structure of this data frame should be maintained if the user
  wishes to change the defaults.

- import_stats:

  **\[deprecated\]** The import of VISPA2 stats has been moved to its
  dedicated function, see
  [import_Vispa2_stats](https://calabrialab.github.io/ISAnalytics/dev/reference/import_Vispa2_stats.md).

## Value

An aggregated data frame

## See also

Other Data cleaning and pre-processing:
[`aggregate_values_by_key()`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_values_by_key.md),
[`compute_near_integrations()`](https://calabrialab.github.io/ISAnalytics/dev/reference/compute_near_integrations.md),
[`default_meta_agg()`](https://calabrialab.github.io/ISAnalytics/dev/reference/default_meta_agg.md),
[`outlier_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/outlier_filter.md),
[`outliers_by_pool_fragments()`](https://calabrialab.github.io/ISAnalytics/dev/reference/outliers_by_pool_fragments.md),
[`purity_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/purity_filter.md),
[`realign_after_collisions()`](https://calabrialab.github.io/ISAnalytics/dev/reference/realign_after_collisions.md),
[`remove_collisions()`](https://calabrialab.github.io/ISAnalytics/dev/reference/remove_collisions.md),
[`threshold_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/threshold_filter.md)

## Examples

``` r
data("association_file", package = "ISAnalytics")
aggreg_meta <- aggregate_metadata(
    association_file = association_file
)
head(aggreg_meta)
#> # A tibble: 6 × 19
#>   SubjectID CellMarker Tissue TimePoint FusionPrimerPCRDate_…¹ LinearPCRDate_min
#>   <chr>     <chr>      <chr>  <chr>     <date>                 <date>           
#> 1 PT001     MNC        BM     0030      2016-11-03             Inf              
#> 2 PT001     MNC        BM     0060      2016-11-03             Inf              
#> 3 PT001     MNC        BM     0090      2016-11-03             Inf              
#> 4 PT001     MNC        BM     0180      2016-11-03             Inf              
#> 5 PT001     MNC        BM     0360      2017-04-21             Inf              
#> 6 PT001     MNC        PB     0030      2016-11-03             Inf              
#> # ℹ abbreviated name: ¹​FusionPrimerPCRDate_min
#> # ℹ 13 more variables: VCN_avg <dbl>, `ng DNA corrected_avg` <dbl>,
#> #   Kapa_avg <dbl>, `ng DNA corrected_sum` <dbl>, ulForPool_sum <dbl>,
#> #   BARCODE_MUX_sum <int>, TRIMMING_FINAL_LTRLC_sum <int>, LV_MAPPED_sum <int>,
#> #   BWA_MAPPED_OVERALL_sum <int>, ISS_MAPPED_OVERALL_sum <int>,
#> #   PCRMethod <chr>, NGSTechnology <chr>, DNAnumber <chr>
```
