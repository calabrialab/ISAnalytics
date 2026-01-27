# Generate a default folder structure, following VISPA2 standards

The function produces a folder structure in the file system at the
provided path that respects VISPA2 standards, with package-included
data.

## Usage

``` r
generate_default_folder_structure(
  type = "correct",
  dir = tempdir(),
  af = "default",
  matrices = "default"
)
```

## Arguments

- type:

  One value between `"correct"`, `"incorrect"` and `"both"`. Tells the
  function wheter to produce a correct structure or introduce some
  errors (mainly for testing purposes).

- dir:

  Path to the folder in which the structure will be produced

- af:

  Either `"default"` for the association file provided as example in the
  package or a custom association file as a data frame

- matrices:

  Either `"default"` for integration matrices provided as example in the
  package or a custom multi-quantification matrix

## Value

A named list containing the path to the association file and the path to
the top level folder(s) of the structure

## Required tags

The function will explicitly check for the presence of these tags:

- project_id

- tag_seq

- vispa_concatenate

## See also

Other Utilities:
[`as_sparse_matrix()`](https://calabrialab.github.io/ISAnalytics/dev/reference/as_sparse_matrix.md),
[`comparison_matrix()`](https://calabrialab.github.io/ISAnalytics/dev/reference/comparison_matrix.md),
[`enable_progress_bars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/enable_progress_bars.md),
[`export_ISA_settings()`](https://calabrialab.github.io/ISAnalytics/dev/reference/export_ISA_settings.md),
[`generate_Vispa2_launch_AF()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_Vispa2_launch_AF.md),
[`generate_blank_association_file()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_blank_association_file.md),
[`import_ISA_settings()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_ISA_settings.md),
[`separate_quant_matrices()`](https://calabrialab.github.io/ISAnalytics/dev/reference/separate_quant_matrices.md),
[`transform_columns()`](https://calabrialab.github.io/ISAnalytics/dev/reference/transform_columns.md)

## Examples

``` r
fs_path <- generate_default_folder_structure(type = "correct")
fs_path
#> $af
#> /tmp/RtmpaYX0p0/asso_file.tsv
#> 
#> $root
#> /tmp/RtmpaYX0p0/fs
#> 
```
