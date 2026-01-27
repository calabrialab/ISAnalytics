# Default metadata aggregation function table

A default columns-function specifications for
[aggregate_metadata](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_metadata.md)

## Usage

``` r
default_meta_agg()
```

## Value

A data frame

## Details

This data frame contains four columns:

- `Column`: holds the name of the column in the association file that
  should be processed

- `Function`: contains either the name of a function (e.g. mean) or a
  purrr-style lambda (e.g. `~ mean(.x, na.rm = TRUE)`). This function
  will be applied to the corresponding column specified in `Column`

- `Args`: optional additional arguments to pass to the corresponding
  function. This is relevant ONLY if the corresponding `Function` is a
  simple function and not a purrr-style lambda.

- `Output_colname`: a `glue` specification that will be used to
  determine a unique output column name. See
  [glue](https://glue.tidyverse.org/reference/glue.html) for more
  details.

## See also

Other Data cleaning and pre-processing:
[`aggregate_metadata()`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_metadata.md),
[`aggregate_values_by_key()`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_values_by_key.md),
[`compute_near_integrations()`](https://calabrialab.github.io/ISAnalytics/dev/reference/compute_near_integrations.md),
[`outlier_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/outlier_filter.md),
[`outliers_by_pool_fragments()`](https://calabrialab.github.io/ISAnalytics/dev/reference/outliers_by_pool_fragments.md),
[`purity_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/purity_filter.md),
[`realign_after_collisions()`](https://calabrialab.github.io/ISAnalytics/dev/reference/realign_after_collisions.md),
[`remove_collisions()`](https://calabrialab.github.io/ISAnalytics/dev/reference/remove_collisions.md),
[`threshold_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/threshold_filter.md)

## Examples

``` r
default_meta_agg()
#> # A tibble: 15 × 4
#>    Column               Function  Args  Output_colname
#>    <chr>                <list>    <lgl> <chr>         
#>  1 FusionPrimerPCRDate  <formula> NA    {.col}_min    
#>  2 LinearPCRDate        <formula> NA    {.col}_min    
#>  3 VCN                  <formula> NA    {.col}_avg    
#>  4 ng DNA corrected     <formula> NA    {.col}_avg    
#>  5 Kapa                 <formula> NA    {.col}_avg    
#>  6 ng DNA corrected     <formula> NA    {.col}_sum    
#>  7 ulForPool            <formula> NA    {.col}_sum    
#>  8 BARCODE_MUX          <formula> NA    {.col}_sum    
#>  9 TRIMMING_FINAL_LTRLC <formula> NA    {.col}_sum    
#> 10 LV_MAPPED            <formula> NA    {.col}_sum    
#> 11 BWA_MAPPED_OVERALL   <formula> NA    {.col}_sum    
#> 12 ISS_MAPPED_OVERALL   <formula> NA    {.col}_sum    
#> 13 PCRMethod            <formula> NA    {.col}        
#> 14 NGSTechnology        <formula> NA    {.col}        
#> 15 DNAnumber            <formula> NA    {.col}        
```
