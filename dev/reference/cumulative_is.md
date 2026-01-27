# Expands integration matrix with the cumulative IS union over time.

**\[experimental\]** Given an input integration matrix that can be
grouped over time, this function adds integrations in groups assuming
that if an integration is observed at time point "t" then it is also
observed in time point "t+1".

## Usage

``` r
cumulative_is(
  x,
  key = c("SubjectID", "CellMarker", "Tissue", "TimePoint"),
  timepoint_col = "TimePoint",
  include_tp_zero = FALSE,
  counts = TRUE,
  keep_og_is = FALSE,
  expand = TRUE
)
```

## Arguments

- x:

  An integration matrix, ideally aggregated via
  [`aggregate_values_by_key()`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_values_by_key.md)

- key:

  The aggregation key used

- timepoint_col:

  The name of the time point column

- include_tp_zero:

  Should time point 0 be included?

- counts:

  Add cumulative counts? Logical

- keep_og_is:

  Keep original set of integrations as a separate column?

- expand:

  If `FALSE`, for each group, the set of integration sites is returned
  in a separate column as a nested table, otherwise the resulting column
  is unnested.

## Value

A data frame

## Required tags

The function will explicitly check for the presence of these tags:

- All columns declared in
  [`mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md)

- Checks if the matrix is annotated by assessing presence of
  [`annotation_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md)

## See also

Other Analysis functions:
[`CIS_grubbs()`](https://calabrialab.github.io/ISAnalytics/dev/reference/CIS_grubbs.md),
[`HSC_population_size_estimate()`](https://calabrialab.github.io/ISAnalytics/dev/reference/HSC_population_size_estimate.md),
[`compute_abundance()`](https://calabrialab.github.io/ISAnalytics/dev/reference/compute_abundance.md),
[`gene_frequency_fisher()`](https://calabrialab.github.io/ISAnalytics/dev/reference/gene_frequency_fisher.md),
[`is_sharing()`](https://calabrialab.github.io/ISAnalytics/dev/reference/is_sharing.md),
[`iss_source()`](https://calabrialab.github.io/ISAnalytics/dev/reference/iss_source.md),
[`sample_statistics()`](https://calabrialab.github.io/ISAnalytics/dev/reference/sample_statistics.md),
[`top_integrations()`](https://calabrialab.github.io/ISAnalytics/dev/reference/top_integrations.md),
[`top_targeted_genes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/top_targeted_genes.md)

## Examples

``` r
data("integration_matrices", package = "ISAnalytics")
data("association_file", package = "ISAnalytics")
aggreg <- aggregate_values_by_key(
    x = integration_matrices,
    association_file = association_file,
    value_cols = c("seqCount", "fragmentEstimate")
)
cumulated_is <- cumulative_is(aggreg)
cumulated_is
#> $coordinates
#> # A tibble: 2,375 × 9
#>    SubjectID CellMarker Tissue TimePoint chr   integration_locus strand GeneName
#>    <chr>     <chr>      <chr>      <dbl> <chr>             <dbl> <chr>  <chr>   
#>  1 PT001     MNC        BM            30 1               8464757 -      RERE    
#>  2 PT001     MNC        BM            30 1              16186297 -      SPEN    
#>  3 PT001     MNC        BM            30 1              40689188 +      RLF     
#>  4 PT001     MNC        BM            30 1             157759338 -      FCRL1   
#>  5 PT001     MNC        BM            30 1             234596545 -      TARBP1  
#>  6 PT001     MNC        BM            30 10            122533902 -      WDR11-A…
#>  7 PT001     MNC        BM            30 11              5306480 +      HBE1    
#>  8 PT001     MNC        BM            30 11             64633964 +      EHD1    
#>  9 PT001     MNC        BM            30 11             65949729 -      PACS1   
#> 10 PT001     MNC        BM            30 11             72097513 +      CLPB    
#> # ℹ 2,365 more rows
#> # ℹ 1 more variable: GeneStrand <chr>
#> 
#> $counts
#> # A tibble: 20 × 5
#>    SubjectID CellMarker Tissue TimePoint is_n_cumulative
#>    <chr>     <chr>      <chr>      <dbl>           <int>
#>  1 PT001     MNC        BM            30              54
#>  2 PT001     MNC        BM            60             147
#>  3 PT001     MNC        BM            90             179
#>  4 PT001     MNC        BM           180             240
#>  5 PT001     MNC        BM           360             240
#>  6 PT001     MNC        PB            30              28
#>  7 PT001     MNC        PB            60              77
#>  8 PT001     MNC        PB            90             104
#>  9 PT001     MNC        PB           180             121
#> 10 PT001     MNC        PB           360             121
#> 11 PT002     MNC        BM            30              98
#> 12 PT002     MNC        BM            60             126
#> 13 PT002     MNC        BM            90             141
#> 14 PT002     MNC        BM           180             184
#> 15 PT002     MNC        BM           360             265
#> 16 PT002     MNC        PB            30              15
#> 17 PT002     MNC        PB            60              26
#> 18 PT002     MNC        PB            90              38
#> 19 PT002     MNC        PB           180              62
#> 20 PT002     MNC        PB           360             109
#> 
```
