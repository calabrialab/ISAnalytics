# Plot of the estimated HSC population size for each patient.

Plot of the estimated HSC population size for each patient.

## Usage

``` r
HSC_population_plot(
  estimates,
  project_name,
  timepoints = "Consecutive",
  models = "Mth Chao (LB)"
)
```

## Arguments

- estimates:

  The estimates data frame, obtained via
  [`HSC_population_size_estimate`](https://calabrialab.github.io/ISAnalytics/dev/reference/HSC_population_size_estimate.md)

- project_name:

  The project name, will be included in the plot title

- timepoints:

  Which time points to plot? One between "All", "Stable" and
  "Consecutive"

- models:

  Name of the models to plot (as they appear in the column of the
  estimates)

## Value

A plot

## See also

Other Plotting functions:
[`CIS_volcano_plot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/CIS_volcano_plot.md),
[`circos_genomic_density()`](https://calabrialab.github.io/ISAnalytics/dev/reference/circos_genomic_density.md),
[`fisher_scatterplot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/fisher_scatterplot.md),
[`integration_alluvial_plot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/integration_alluvial_plot.md),
[`sharing_heatmap()`](https://calabrialab.github.io/ISAnalytics/dev/reference/sharing_heatmap.md),
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
aggreg_meta <- aggregate_metadata(
    association_file = association_file
)
estimate <- HSC_population_size_estimate(
    x = aggreg,
    metadata = aggreg_meta,
    stable_timepoints = c(90, 180, 360),
    cell_type = "Other"
)
#> Calculating number of IS for each group...
p <- HSC_population_plot(estimate$est, "PJ01")
#> Warning: Arguments in `...` must be used.
#> ✖ Problematic arguments:
#> • na.rm = TRUE
#> • se = TRUE
#> ℹ Did you misspell an argument name?
p
```
