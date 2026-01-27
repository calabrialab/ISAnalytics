# Filter integration sites based on purity.

**\[stable\]** Filter that targets possible contamination between cell
lines based on a numeric quantification (likely abundance or sequence
count).

## Usage

``` r
purity_filter(
  x,
  lineages = blood_lineages_default(),
  aggregation_key = c("SubjectID", "CellMarker", "Tissue", "TimePoint"),
  group_key = c("CellMarker", "Tissue"),
  selected_groups = NULL,
  join_on = "CellMarker",
  min_value = 3,
  impurity_threshold = 10,
  by_timepoint = TRUE,
  timepoint_column = "TimePoint",
  value_column = "seqCount_sum"
)
```

## Arguments

- x:

  An aggregated integration matrix, obtained via
  [`aggregate_values_by_key()`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_values_by_key.md)

- lineages:

  A data frame containing cell lineages information

- aggregation_key:

  The key used for aggregating `x`

- group_key:

  A character vector of column names for re-aggregation. Column names
  must be either in `x` or in `lineages`. See details.

- selected_groups:

  Either NULL, a character vector or a data frame for group selection.
  See details.

- join_on:

  Common columns to perform a join operation on

- min_value:

  A minimum value to filter the input matrix. Integrations with a value
  strictly lower than `min_value` are excluded (dropped) from the
  output.

- impurity_threshold:

  The ratio threshold for impurity in groups

- by_timepoint:

  Should filtering be applied on each time point? If `FALSE`, all time
  points are merged together

- timepoint_column:

  Column in `x` containing the time point

- value_column:

  Column in `x` containing the numeric quantification of interest

## Value

A data frame

## Details

### Setting input arguments

The input matrix can be re-aggregated with the provided `group_key`
argument. This key contains the names of the columns to group on
(besides the columns holding genomic coordinates of the integration
sites) and must be contained in at least one of `x` or `lineages` data
frames. If the key is not found only in `x`, then a join operation with
the `lineages` data frame is performed on the common column(s)
`join_on`.

### Group selection

It is possible for the user to specify on which groups the logic of the
filter should be applied to. For example: if we have
`group_key = c("HematoLineage")` and we set
`selected_groups = c("CD34", "Myeloid","Lymphoid")` it means that a
single integration will be evaluated for the filter only for groups that
have the values of "CD34", "Myeloid" and "Lymphoid" in the
"HematoLineage" column. If the same integration is present in other
groups it is kept as it is. `selected_groups` can be set to `NULL` if we
want the logic to apply to every group present in the data frame, it can
be set as a simple character vector as the example above if the group
key has length 1 (and there is no need to filter on time point). If the
group key is longer than 1 then the filter is applied only on the first
element of the key.

If a more refined selection on groups is needed, a data frame can be
provided instead:

    group_key = c("CellMarker", "Tissue")
    selected_groups = tibble::tribble(
    ~ CellMarker, ~ Tissue,
    "CD34", "BM",
    "CD14", "BM",
    "CD14", "PB"
    )

Columns in the data frame should be the same as group key (plus,
eventually, the time point column). In this example only those groups
identified by the rows in the provided data frame are processed.

## See also

Other Data cleaning and pre-processing:
[`aggregate_metadata()`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_metadata.md),
[`aggregate_values_by_key()`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_values_by_key.md),
[`compute_near_integrations()`](https://calabrialab.github.io/ISAnalytics/dev/reference/compute_near_integrations.md),
[`default_meta_agg()`](https://calabrialab.github.io/ISAnalytics/dev/reference/default_meta_agg.md),
[`outlier_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/outlier_filter.md),
[`outliers_by_pool_fragments()`](https://calabrialab.github.io/ISAnalytics/dev/reference/outliers_by_pool_fragments.md),
[`realign_after_collisions()`](https://calabrialab.github.io/ISAnalytics/dev/reference/realign_after_collisions.md),
[`remove_collisions()`](https://calabrialab.github.io/ISAnalytics/dev/reference/remove_collisions.md),
[`threshold_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/threshold_filter.md)

## Examples

``` r
data("integration_matrices", package = "ISAnalytics")
data("association_file", package = "ISAnalytics")
aggreg <- aggregate_values_by_key(
    x = integration_matrices,
    association_file = association_file,
    value_cols = c("seqCount", "fragmentEstimate")
)
filtered_by_purity <- purity_filter(
    x = aggreg,
    value_column = "seqCount_sum"
)
head(filtered_by_purity)
#> # A tibble: 6 × 9
#>   chr   integration_locus strand GeneName GeneStrand TimePoint CellMarker Tissue
#>   <chr>             <dbl> <chr>  <chr>    <chr>      <chr>     <chr>      <chr> 
#> 1 1              27058748 -      ARID1A   +          0060      MNC        BM    
#> 2 1              89337032 +      GTF2B    -          0360      MNC        BM    
#> 3 1              89337032 +      GTF2B    -          0360      MNC        PB    
#> 4 10             74765793 -      P4HA1    -          0360      MNC        BM    
#> 5 10             74765793 -      P4HA1    -          0360      MNC        PB    
#> 6 11              3899973 +      STIM1    +          0090      MNC        BM    
#> # ℹ 1 more variable: Value <dbl>
```
