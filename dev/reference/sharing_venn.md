# Produce tables to plot sharing venn or euler diagrams.

**\[stable\]** This function processes a sharing data frame obtained via
[`is_sharing()`](https://calabrialab.github.io/ISAnalytics/dev/reference/is_sharing.md)
with the option `table_for_venn = TRUE` to obtain a list of objects that
can be plotted as venn or euler diagrams.

## Usage

``` r
sharing_venn(sharing_df, row_range = NULL, euler = TRUE)
```

## Arguments

- sharing_df:

  The sharing data frame

- row_range:

  Either `NULL` or a numeric vector of row indexes (e.g. `c(1, 4, 5)`
  will produce tables only for rows 1, 4 and 5)

- euler:

  If `TRUE` will produce tables for euler diagrams, otherwise will
  produce tables for venn diagrams

## Value

A list of data frames

## Details

The functions requires the package
[eulerr](https://jolars.github.io/eulerr/index.html). Each row of the
input data frame is representable as a venn/euler diagram. The function
allows to specify a range of row indexes to obtain a list of plottable
objects all at once, leave it to NULL to process all rows.

To actually plot the data it is sufficient to call the function
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) and specify
optional customization arguments. See [eulerr
docs](https://jolars.github.io/eulerr/reference/plot.euler.html) for
more detail on this.

## See also

Other Plotting functions:
[`CIS_volcano_plot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/CIS_volcano_plot.md),
[`HSC_population_plot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/HSC_population_plot.md),
[`circos_genomic_density()`](https://calabrialab.github.io/ISAnalytics/dev/reference/circos_genomic_density.md),
[`fisher_scatterplot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/fisher_scatterplot.md),
[`integration_alluvial_plot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/integration_alluvial_plot.md),
[`sharing_heatmap()`](https://calabrialab.github.io/ISAnalytics/dev/reference/sharing_heatmap.md),
[`top_abund_tableGrob()`](https://calabrialab.github.io/ISAnalytics/dev/reference/top_abund_tableGrob.md),
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
sharing <- is_sharing(aggreg, n_comp = 3, table_for_venn = TRUE)
venn_tbls <- sharing_venn(sharing, row_range = 1:3, euler = FALSE)
venn_tbls
#> [[1]]
#> 3 set Venn diagram 
#> 
#>                       h     k    a    b  phi
#> PT001_MNC_BM_0030 -0.42 -0.36 1.05 1.05 3.76
#> PT001_MNC_BM_0060  0.42 -0.36 1.05 1.05 3.76
#> PT001_MNC_BM_0090  0.00  0.36 1.05 1.05 3.76
#> 
#> [[2]]
#> 3 set Venn diagram 
#> 
#>                       h     k    a    b  phi
#> PT001_MNC_BM_0030 -0.42 -0.36 1.05 1.05 3.76
#> PT001_MNC_BM_0060  0.42 -0.36 1.05 1.05 3.76
#> PT001_MNC_BM_0180  0.00  0.36 1.05 1.05 3.76
#> 
#> [[3]]
#> 3 set Venn diagram 
#> 
#>                       h     k    a    b  phi
#> PT001_MNC_BM_0030 -0.42 -0.36 1.05 1.05 3.76
#> PT001_MNC_BM_0060  0.42 -0.36 1.05 1.05 3.76
#> PT001_MNC_BM_0360  0.00  0.36 1.05 1.05 3.76
#> 
plot(venn_tbls[[1]])
```
