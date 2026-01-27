# Separate a multiple-quantification matrix into single quantification matrices.

**\[stable\]** The function separates a single multi-quantification
integration matrix, obtained via
[comparison_matrix](https://calabrialab.github.io/ISAnalytics/dev/reference/comparison_matrix.md),
into single quantification matrices as a named list of tibbles.

## Usage

``` r
separate_quant_matrices(
  x,
  fragmentEstimate = "fragmentEstimate",
  seqCount = "seqCount",
  barcodeCount = "barcodeCount",
  cellCount = "cellCount",
  ShsCount = "ShsCount",
  key = c(mandatory_IS_vars(), annotation_IS_vars(), "CompleteAmplificationID")
)
```

## Arguments

- x:

  Single integration matrix with multiple quantification value columns,
  obtained via
  [comparison_matrix](https://calabrialab.github.io/ISAnalytics/dev/reference/comparison_matrix.md).

- fragmentEstimate:

  Name of the fragment estimate values column in input

- seqCount:

  Name of the sequence count values column in input

- barcodeCount:

  Name of the barcode count values column in input

- cellCount:

  Name of the cell count values column in input

- ShsCount:

  Name of the shs count values column in input

- key:

  Key columns to perform the joining operation

## Value

A named list of data frames, where names are quantification types

## See also

[quantification_types](https://calabrialab.github.io/ISAnalytics/dev/reference/quantification_types.md)

Other Utilities:
[`as_sparse_matrix()`](https://calabrialab.github.io/ISAnalytics/dev/reference/as_sparse_matrix.md),
[`comparison_matrix()`](https://calabrialab.github.io/ISAnalytics/dev/reference/comparison_matrix.md),
[`enable_progress_bars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/enable_progress_bars.md),
[`export_ISA_settings()`](https://calabrialab.github.io/ISAnalytics/dev/reference/export_ISA_settings.md),
[`generate_Vispa2_launch_AF()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_Vispa2_launch_AF.md),
[`generate_blank_association_file()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_blank_association_file.md),
[`generate_default_folder_structure()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_default_folder_structure.md),
[`import_ISA_settings()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_ISA_settings.md),
[`transform_columns()`](https://calabrialab.github.io/ISAnalytics/dev/reference/transform_columns.md)

## Examples

``` r
data("integration_matrices", package = "ISAnalytics")
separated <- separate_quant_matrices(
    integration_matrices
)
```
