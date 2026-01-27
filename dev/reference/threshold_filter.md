# Filter data frames with custom predicates

**\[deprecated\]** This function is deprecated and it's likely going to
be dropped in the next release cycle.

Filter a single data frame or a list of data frames with custom
predicates assembled from the function parameters.

## Usage

``` r
threshold_filter(x, threshold, cols_to_compare = "Value", comparators = ">")
```

## Arguments

- x:

  A data frame or a list of data frames

- threshold:

  A numeric/integer vector or a named list of numeric/integer vectors

- cols_to_compare:

  A character vector or a named list of character vectors

- comparators:

  A character vector or a named list of character vectors. Must be one
  of the allowed values between `c("<", ">", "==", "!=", ">=", "<=")`

## Value

A data frame or a list of data frames

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
[`remove_collisions()`](https://calabrialab.github.io/ISAnalytics/dev/reference/remove_collisions.md)

## Examples

``` r
if (FALSE) { # \dontrun{
example_df <- tibble::tibble(
    a = c(20, 30, 40),
    b = c(40, 50, 60),
    c = c("a", "b", "c"),
    d = c(3L, 4L, 5L)
)
example_list <- list(
    first = example_df,
    second = example_df,
    third = example_df
)

filtered <- threshold_filter(example_list,
    threshold = list(
        first = c(20, 60),
        third = c(25)
    ),
    cols_to_compare = list(
        first = c("a", "b"),
        third = c("a")
    ),
    comparators = list(
        first = c(">", "<"),
        third = c(">=")
    )
)
filtered
} # }
```
