# Import integration matrices from association file.

**\[defunct\]** This function was deprecated to avoid redundancy. Please
refer to
[`import_parallel_Vispa2Matrices`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_parallel_Vispa2Matrices.md).

## Usage

``` r
import_parallel_Vispa2Matrices_auto(
  association_file,
  quantification_type,
  matrix_type = "annotated",
  workers = 2,
  multi_quant_matrix = TRUE,
  patterns = NULL,
  matching_opt = matching_options(),
  export_report_path = NULL,
  ...
)
```

## Value

A data frame or a list
