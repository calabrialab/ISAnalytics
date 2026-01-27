# Trace a circos plot of genomic densities.

**\[stable\]** For this functionality the suggested package
[circlize](https://cran.r-project.org/web/packages/circlize/index.html)
is required. Please note that this function is a simple wrapper of basic
`circlize` functions, for an in-depth explanation on how the functions
work and additional arguments please refer to the official documentation
[Circular Visualization in
R](https://jokergoo.github.io/circlize_book/book/)

## Usage

``` r
circos_genomic_density(
  data,
  gene_labels = NULL,
  label_col = NULL,
  cytoband_specie = "hg19",
  track_colors = "navyblue",
  grDevice = c("png", "pdf", "svg", "jpeg", "bmp", "tiff", "default"),
  file_path = getwd(),
  ...
)
```

## Arguments

- data:

  Either a single integration matrix or a list of integration matrices.
  If a list is provided, a separate density track for each data frame is
  plotted.

- gene_labels:

  Either `NULL` or a data frame in bed format. See details.

- label_col:

  Numeric index of the column of `gene_labels` that contains the actual
  labels. Relevant only if `gene_labels` is not set to `NULL`.

- cytoband_specie:

  Specie for initializing the cytoband

- track_colors:

  Colors to give to density tracks. If more than one integration matrix
  is provided as `data` should be of the same length. Values are
  recycled if length of `track_colors` is smaller than the length of the
  input data.

- grDevice:

  The graphical device where the plot should be traced. `default`, if
  executing from RStudio is the viewer.

- file_path:

  If a device other than `default` is chosen, the path on disk where the
  file should be saved. Defaults to
  `{current directory}/circos_plot.{device}`.

- ...:

  Additional named arguments to pass on to chosen device,
  [`circlize::circos.par()`](https://rdrr.io/pkg/circlize/man/circos.par.html),
  [`circlize::circos.genomicDensity()`](https://rdrr.io/pkg/circlize/man/circos.genomicDensity.html)
  and
  [`circlize::circos.genomicLabels()`](https://rdrr.io/pkg/circlize/man/circos.genomicLabels.html)

## Value

`NULL`

## Details

### Providing genomic labels

If genomic labels should be plotted alongside genomic density tracks,
the user should provide them as a simple data frame in standard bed
format, namely `chr`, `start`, `end` plus a column containing the
labels. NOTE: if the user decides to plot on the default device (viewer
in RStudio), he must ensure there is enough space for all elements to be
plotted, otherwise an error message is thrown.

## See also

Other Plotting functions:
[`CIS_volcano_plot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/CIS_volcano_plot.md),
[`HSC_population_plot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/HSC_population_plot.md),
[`fisher_scatterplot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/fisher_scatterplot.md),
[`integration_alluvial_plot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/integration_alluvial_plot.md),
[`sharing_heatmap()`](https://calabrialab.github.io/ISAnalytics/dev/reference/sharing_heatmap.md),
[`sharing_venn()`](https://calabrialab.github.io/ISAnalytics/dev/reference/sharing_venn.md),
[`top_abund_tableGrob()`](https://calabrialab.github.io/ISAnalytics/dev/reference/top_abund_tableGrob.md),
[`top_cis_overtime_heatmap()`](https://calabrialab.github.io/ISAnalytics/dev/reference/top_cis_overtime_heatmap.md)

## Examples

``` r
# \donttest{
data("integration_matrices", package = "ISAnalytics")
data("association_file", package = "ISAnalytics")
aggreg <- aggregate_values_by_key(
    x = integration_matrices,
    association_file = association_file,
    value_cols = c("seqCount", "fragmentEstimate")
)
by_subj <- aggreg |>
    dplyr::group_by(.data$SubjectID) |>
    dplyr::group_split()
circos_genomic_density(by_subj,
    track_colors = c("navyblue", "gold"),
    grDevice = "default", track.height = 0.1
)

# }
```
