# Converts tidy integration matrices in the original sparse matrix form.

**\[stable\]** This function is particularly useful when a sparse matrix
structure is needed by a specific function (mainly from other packages).

## Usage

``` r
as_sparse_matrix(
  x,
  single_value_col = "Value",
  fragmentEstimate = "fragmentEstimate",
  seqCount = "seqCount",
  barcodeCount = "barcodeCount",
  cellCount = "cellCount",
  ShsCount = "ShsCount",
  key = pcr_id_column()
)
```

## Arguments

- x:

  A single tidy integration matrix or a list of integration matrices.
  Supports also multi-quantification matrices obtained via
  [comparison_matrix](https://calabrialab.github.io/ISAnalytics/dev/reference/comparison_matrix.md)

- single_value_col:

  Name of the column containing the values when providing a
  single-quantification matrix

- fragmentEstimate:

  For multi-quantification matrix support: the name of the fragment
  estimate values column

- seqCount:

  For multi-quantification matrix support: the name of the sequence
  count values column

- barcodeCount:

  For multi-quantification matrix support: the name of the barcode count
  values column

- cellCount:

  For multi-quantification matrix support: the name of the cell count
  values column

- ShsCount:

  For multi-quantification matrix support: the name of the Shs Count
  values column

- key:

  The name of the sample identifier fields (for aggregated matrices can
  be a vector with more than 1 element)

## Value

Depending on input, 2 possible outputs:

- A single sparse matrix (data frame) if input is a single
  quantification matrix

- A list of sparse matrices divided by quantification if input is a
  single multi-quantification matrix or a list of matrices

## See also

Other Utilities:
[`comparison_matrix()`](https://calabrialab.github.io/ISAnalytics/dev/reference/comparison_matrix.md),
[`enable_progress_bars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/enable_progress_bars.md),
[`export_ISA_settings()`](https://calabrialab.github.io/ISAnalytics/dev/reference/export_ISA_settings.md),
[`generate_Vispa2_launch_AF()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_Vispa2_launch_AF.md),
[`generate_blank_association_file()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_blank_association_file.md),
[`generate_default_folder_structure()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_default_folder_structure.md),
[`import_ISA_settings()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_ISA_settings.md),
[`separate_quant_matrices()`](https://calabrialab.github.io/ISAnalytics/dev/reference/separate_quant_matrices.md),
[`transform_columns()`](https://calabrialab.github.io/ISAnalytics/dev/reference/transform_columns.md)

## Examples

``` r
data("integration_matrices", package = "ISAnalytics")
sparse <- as_sparse_matrix(integration_matrices)
```
