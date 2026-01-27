# Obtain a single integration matrix from individual quantification matrices.

**\[stable\]** Takes a list of integration matrices referring to
different quantification types and merges them into a single data frame
with multiple value columns, each renamed according to their
quantification type of reference.

## Usage

``` r
comparison_matrix(
  x,
  fragmentEstimate = "fragmentEstimate",
  seqCount = "seqCount",
  barcodeCount = "barcodeCount",
  cellCount = "cellCount",
  ShsCount = "ShsCount",
  value_col_name = "Value"
)
```

## Arguments

- x:

  A named list of integration matrices, ideally obtained via
  [import_parallel_Vispa2Matrices](https://calabrialab.github.io/ISAnalytics/dev/reference/import_parallel_Vispa2Matrices.md).
  Names must be quantification types in
  [`quantification_types()`](https://calabrialab.github.io/ISAnalytics/dev/reference/quantification_types.md).

- fragmentEstimate:

  The name of the output column for fragment estimate values

- seqCount:

  The name of the output column for sequence count values

- barcodeCount:

  The name of the output column for barcode count values

- cellCount:

  The name of the output column for cell count values

- ShsCount:

  The name of the output column for Shs count values

- value_col_name:

  Name of the column containing the corresponding values in the single
  matrices

## Value

A single data frame

## See also

[quantification_types](https://calabrialab.github.io/ISAnalytics/dev/reference/quantification_types.md)

Other Utilities:
[`as_sparse_matrix()`](https://calabrialab.github.io/ISAnalytics/dev/reference/as_sparse_matrix.md),
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
sc <- tibble::tribble(
    ~chr, ~integration_locus, ~strand, ~CompleteAmplificationID, ~Value,
    "1", 45324, "+", "ID1", 543,
    "2", 52423, "-", "ID1", 42,
    "6", 54623, "-", "ID2", 67,
    "X", 12314, "+", "ID3", 8
)
fe <- tibble::tribble(
    ~chr, ~integration_locus, ~strand, ~CompleteAmplificationID, ~Value,
    "1", 45324, "+", "ID1", 56.76,
    "2", 52423, "-", "ID1", 78.32,
    "6", 54623, "-", "ID2", 123.45,
    "X", 12314, "+", "ID3", 5.34
)
comparison_matrix(list(
    fragmentEstimate = fe,
    seqCount = sc
))
#> # A tibble: 4 × 6
#>   chr   integration_locus strand CompleteAmplificationID fragmentEstimate
#>   <chr>             <dbl> <chr>  <chr>                              <dbl>
#> 1 1                 45324 +      ID1                                56.8 
#> 2 2                 52423 -      ID1                                78.3 
#> 3 6                 54623 -      ID2                               123.  
#> 4 X                 12314 +      ID3                                 5.34
#> # ℹ 1 more variable: seqCount <dbl>
```
