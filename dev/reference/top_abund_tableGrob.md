# Summary top abundant tableGrobs for plots.

Produce summary tableGrobs as R graphics. For this functionality the
suggested package
[gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html)
is required. To visualize the resulting object:

    gridExtra::grid.arrange(tableGrob)

## Usage

``` r
top_abund_tableGrob(
  df,
  id_cols = mandatory_IS_vars(),
  quant_col = "fragmentEstimate_sum_PercAbundance",
  by = "TimePoint",
  alluvial_plot = NULL,
  top_n = 10,
  tbl_cols = "GeneName",
  include_id_cols = FALSE,
  digits = 2,
  perc_symbol = TRUE,
  transform_by = NULL
)
```

## Arguments

- df:

  A data frame

- id_cols:

  Character vector of id column names. To plot after alluvial, these
  columns must be the same as the `alluvia` argument of
  [integration_alluvial_plot](https://calabrialab.github.io/ISAnalytics/dev/reference/integration_alluvial_plot.md).

- quant_col:

  Column name holding the quantification value. To plot after alluvial,
  these columns must be the same as the `plot_y` argument of
  [integration_alluvial_plot](https://calabrialab.github.io/ISAnalytics/dev/reference/integration_alluvial_plot.md).

- by:

  The column name to subdivide tables for. The function will produce one
  table for each distinct value in `by`. To plot after alluvial, these
  columns must be the same as the `plot_x` argument of
  [integration_alluvial_plot](https://calabrialab.github.io/ISAnalytics/dev/reference/integration_alluvial_plot.md).

- alluvial_plot:

  Either NULL or an alluvial plot for color mapping between values of y.

- top_n:

  Integer. How many rows should the table contain at most?

- tbl_cols:

  Table columns to show in the final output besides `quant_col`.

- include_id_cols:

  Logical. Include `id_cols` in the output?

- digits:

  Integer. Digits to show for the quantification column

- perc_symbol:

  Logical. Show percentage symbol in the quantification column?

- transform_by:

  Either a function or a purrr-style lambda. This function is applied to
  the column `by` before separating columns. If `NULL` no function is
  applied. Useful to modify column order in final table.

## Value

A tableGrob object

## See also

Other Plotting functions:
[`CIS_volcano_plot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/CIS_volcano_plot.md),
[`HSC_population_plot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/HSC_population_plot.md),
[`circos_genomic_density()`](https://calabrialab.github.io/ISAnalytics/dev/reference/circos_genomic_density.md),
[`fisher_scatterplot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/fisher_scatterplot.md),
[`integration_alluvial_plot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/integration_alluvial_plot.md),
[`sharing_heatmap()`](https://calabrialab.github.io/ISAnalytics/dev/reference/sharing_heatmap.md),
[`sharing_venn()`](https://calabrialab.github.io/ISAnalytics/dev/reference/sharing_venn.md),
[`top_cis_overtime_heatmap()`](https://calabrialab.github.io/ISAnalytics/dev/reference/top_cis_overtime_heatmap.md)

## Examples

``` r
data("integration_matrices", package = "ISAnalytics")
data("association_file", package = "ISAnalytics")
aggreg <- aggregate_values_by_key(
    x = integration_matrices,
    association_file = association_file,
    value_cols = c("seqCount", "fragmentEstimate")
)
abund <- compute_abundance(x = aggreg)
grob <- top_abund_tableGrob(abund)
gridExtra::grid.arrange(grob)


# with transform
grob <- top_abund_tableGrob(abund, transform_by = ~ as.numeric(.x))
```
