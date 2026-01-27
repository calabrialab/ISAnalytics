# Aggregates matrices values based on specified key.

**\[stable\]** Performs aggregation on values contained in the
integration matrices based on the key and the specified lambda. For more
details on how to use this function:
[`vignette("workflow_start", package = "ISAnalytics")`](https://calabrialab.github.io/ISAnalytics/dev/articles/workflow_start.md)

## Usage

``` r
aggregate_values_by_key(
  x,
  association_file,
  value_cols = "Value",
  key = c("SubjectID", "CellMarker", "Tissue", "TimePoint"),
  lambda = list(sum = ~sum(.x, na.rm = TRUE)),
  group = c(mandatory_IS_vars(), annotation_IS_vars()),
  join_af_by = "CompleteAmplificationID"
)
```

## Arguments

- x:

  A single integration matrix or a list of imported integration matrices

- association_file:

  The imported association file

- value_cols:

  A character vector containing the names of the columns to apply the
  given lambdas. Must be numeric or integer columns.

- key:

  A string or a character vector with column names of the association
  file to take as key

- lambda:

  A named list of functions or purrr-style lambdas. See details section.

- group:

  Other variables to include in the grouping besides `key`, can be set
  to NULL

- join_af_by:

  A character vector representing the joining key between the matrix and
  the metadata. Useful to re-aggregate already aggregated matrices.

## Value

A list of data frames or a single data frame aggregated according to the
specified arguments

## Details

### Setting the lambda parameter

The lambda parameter should always contain a named list of either
functions or purrr-style lambdas. It is also possible to specify the
namespace of the function in both ways, for example:

    lambda = list(sum = sum, desc = psych::describe)

Using purrr-style lambdas allows to specify arguments for the functions,
keeping in mind that the first parameter should always be `.x`:

    lambda = list(sum = ~sum(.x, na.rm = TRUE))

It is also possible to use custom user-defined functions, keeping in
mind that the symbol will be evaluated in the calling environment, for
example if the function is called in the global environment and lambda
contains "foo" as a function, "foo" will be evaluated in the global
environment.

    foo <- function(x) {
      sum(x)
    }

    lambda = list(sum = ~sum(.x, na.rm = TRUE), foo = foo)

    # Or with lambda notation
    lambda = list(sum = ~sum(.x, na.rm = TRUE), foo = ~foo(.x))

### Constraints on aggregation functions

Functions passed in the lambda parameters must respect a few constraints
to properly work and it's the user responsibility to ensure this.

- Functions have to accept as input a numeric or integer vector

- Function should return a single value or a list/data frame: if a list
  or a data frame is returned as a result, all the columns will be added
  to the final data frame.

## See also

Other Data cleaning and pre-processing:
[`aggregate_metadata()`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_metadata.md),
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
data("integration_matrices", package = "ISAnalytics")
data("association_file", package = "ISAnalytics")
aggreg <- aggregate_values_by_key(
    x = integration_matrices,
    association_file = association_file,
    value_cols = c("seqCount", "fragmentEstimate")
)
head(aggreg)
#> # A tibble: 6 × 11
#>   chr   integration_locus strand GeneName GeneStrand SubjectID CellMarker Tissue
#>   <chr>             <dbl> <chr>  <chr>    <chr>      <chr>     <chr>      <chr> 
#> 1 1               8464757 -      RERE     -          PT001     MNC        BM    
#> 2 1               8464757 -      RERE     -          PT001     MNC        BM    
#> 3 1               8607357 +      RERE     -          PT001     MNC        BM    
#> 4 1               8607357 +      RERE     -          PT001     MNC        BM    
#> 5 1               8607357 +      RERE     -          PT001     MNC        BM    
#> 6 1               8607362 -      RERE     -          PT001     MNC        BM    
#> # ℹ 3 more variables: TimePoint <chr>, seqCount_sum <dbl>,
#> #   fragmentEstimate_sum <dbl>
```
