# Integrations cumulative count in time by sample

**\[defunct\]** This function was deprecated in favour of a single
function, please use `cumulative_is` instead.

## Usage

``` r
cumulative_count_union(
  x,
  association_file = NULL,
  timepoint_column = "TimePoint",
  key = c("SubjectID", "CellMarker", "Tissue", "TimePoint"),
  include_tp_zero = FALSE,
  zero = "0000",
  aggregate = FALSE,
  ...
)
```

## Arguments

- x:

  A simple integration matrix or an aggregated matrix (see details)

- association_file:

  NULL or the association file for x if `aggregate` is set to TRUE

- timepoint_column:

  What is the name of the time point column?

- key:

  The aggregation key - must always contain the `timepoint_column`

- include_tp_zero:

  Include timepoint 0?

- zero:

  How is 0 coded in the data frame?

- aggregate:

  Should x be aggregated?

- ...:

  Additional parameters to pass to `aggregate_values_by_key`

## Value

A data frame

## Examples

``` r
if (FALSE) { # \dontrun{
data("integration_matrices", package = "ISAnalytics")
data("association_file", package = "ISAnalytics")
aggreg <- aggregate_values_by_key(
    x = integration_matrices,
    association_file = association_file,
    value_cols = c("seqCount", "fragmentEstimate")
)
cumulative_count <- cumulative_count_union(aggreg)
cumulative_count
} # }
```
