# Creates a reduced association file for a VISPA2 run, given project and pool

The function selects the appropriate columns and prepares a file for the
launch of VISPA2 pipeline for each project/pool pair specified.

## Usage

``` r
generate_Vispa2_launch_AF(association_file, project, pool, path)
```

## Arguments

- association_file:

  The imported association file (via
  [`import_association_file()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_association_file.md))

- project:

  A vector of characters containing project names

- pool:

  A vector of characters containing pool names

- path:

  A single string representing the path to the folder where files should
  be written. If the folder doesn't exist it will be created.

## Value

`NULL`

## Details

Note: the function is vectorized, meaning you can specify more than one
project and more than one pool as vectors of characters, but you must
ensure that:

- Both `project` and `pool` vectors have the same length

- You correclty type names in corresponding positions, for example
  c("PJ01", "PJ01") - c("POOL01", "POOL02"). If you type a pool in the
  position of a corresponding project that doesn't match no file will be
  produced since that pool doesn't exist in the corresponding project.

## Required tags

The function will explicitly check for the presence of these tags:

- cell_marker

- fusion_id

- pcr_repl_id

- pool_id

- project_id

- subject

- tag_id

- tissue

- tp_days

- vector_id

The names of the pools in the `pool` argument is checked against the
column corresponding to the `pool_id` tag.

## See also

Other Utilities:
[`as_sparse_matrix()`](https://calabrialab.github.io/ISAnalytics/dev/reference/as_sparse_matrix.md),
[`comparison_matrix()`](https://calabrialab.github.io/ISAnalytics/dev/reference/comparison_matrix.md),
[`enable_progress_bars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/enable_progress_bars.md),
[`export_ISA_settings()`](https://calabrialab.github.io/ISAnalytics/dev/reference/export_ISA_settings.md),
[`generate_blank_association_file()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_blank_association_file.md),
[`generate_default_folder_structure()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_default_folder_structure.md),
[`import_ISA_settings()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_ISA_settings.md),
[`separate_quant_matrices()`](https://calabrialab.github.io/ISAnalytics/dev/reference/separate_quant_matrices.md),
[`transform_columns()`](https://calabrialab.github.io/ISAnalytics/dev/reference/transform_columns.md)

## Examples

``` r
temp <- tempdir()
data("association_file", package = "ISAnalytics")
generate_Vispa2_launch_AF(association_file, "PJ01", "POOL01", temp)
```
