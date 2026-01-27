# Plot IS sharing heatmaps.

**\[stable\]** Displays the IS sharing calculated via
[is_sharing](https://calabrialab.github.io/ISAnalytics/dev/reference/is_sharing.md)
as heatmaps.

## Usage

``` r
sharing_heatmap(
  sharing_df,
  show_on_x = "g1",
  show_on_y = "g2",
  absolute_sharing_col = "shared",
  title_annot = NULL,
  plot_relative_sharing = TRUE,
  rel_sharing_col = c("on_g1", "on_union"),
  show_perc_symbol_rel = TRUE,
  interactive = FALSE
)
```

## Arguments

- sharing_df:

  The data frame containing the IS sharing data

- show_on_x:

  Name of the column to plot on the x axis

- show_on_y:

  Name of the column to plot on the y axis

- absolute_sharing_col:

  Name of the column that contains the absolute values of IS sharing

- title_annot:

  Additional text to display in the title

- plot_relative_sharing:

  Logical. Compute heatmaps also for relative sharing?

- rel_sharing_col:

  Names of the columns to consider as relative sharing. The function is
  going to plot one heatmap per column in this argument.

- show_perc_symbol_rel:

  Logical. Only relevant if `plot_relative_sharing` is set to TRUE,
  should the percentage symbol be displayed in relative heatmaps?

- interactive:

  Logical. Requires the package
  [plotly](https://plotly.com/r/getting-started/) is required for this
  functionality. Returns the heatmaps as interactive HTML widgets.

## Value

A list of plots or widgets

## See also

[is_sharing](https://calabrialab.github.io/ISAnalytics/dev/reference/is_sharing.md)

Other Plotting functions:
[`CIS_volcano_plot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/CIS_volcano_plot.md),
[`HSC_population_plot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/HSC_population_plot.md),
[`circos_genomic_density()`](https://calabrialab.github.io/ISAnalytics/dev/reference/circos_genomic_density.md),
[`fisher_scatterplot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/fisher_scatterplot.md),
[`integration_alluvial_plot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/integration_alluvial_plot.md),
[`sharing_venn()`](https://calabrialab.github.io/ISAnalytics/dev/reference/sharing_venn.md),
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
sharing <- is_sharing(aggreg,
    minimal = FALSE,
    include_self_comp = TRUE
)
sharing_heatmaps <- sharing_heatmap(sharing_df = sharing)
sharing_heatmaps$absolute

sharing_heatmaps$on_g1

sharing_heatmaps$on_union
```
