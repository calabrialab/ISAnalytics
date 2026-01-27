# Import a single integration matrix from file

**\[stable\]** This function allows to read and import an integration
matrix (ideally produced by VISPA2) and converts it to a tidy format.

## Usage

``` r
import_single_Vispa2Matrix(
  path,
  separator = "\t",
  additional_cols = NULL,
  transformations = NULL,
  sample_names_to = pcr_id_column(),
  values_to = "Value",
  to_exclude = lifecycle::deprecated(),
  keep_excluded = lifecycle::deprecated()
)
```

## Arguments

- path:

  The path to the file on disk

- separator:

  The column delimiter used, defaults to `\t`

- additional_cols:

  Either `NULL`, a named character vector or a named list. See details.

- transformations:

  Either `NULL` or a named list of purrr-style lambdas where names are
  column names the function should be applied to.

- sample_names_to:

  Name of the output column holding the sample identifier. Defaults to
  [`pcr_id_column()`](https://calabrialab.github.io/ISAnalytics/dev/reference/pcr_id_column.md)

- values_to:

  Name of the output column holding the quantification values. Defaults
  to `Value`.

- to_exclude:

  **\[deprecated\]** Deprecated. Use `additonal_cols` instead

- keep_excluded:

  **\[deprecated\]** Deprecated. Use `additonal_cols` instead

## Value

A data frame object in tidy format

## Details

### Additional columns

Additional columns are annotation columns present in the integration
matrix to import that are not

- part of the mandatory IS vars (see
  [`mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md))

- part of the annotation IS vars (see
  [`annotation_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md))

- the sample identifier column

- the quantification column

When specified they tell the function how to treat those columns in the
import phase, by providing a named character vector, where names
correspond to the additional column names and values are a choice of the
following:

- `"char"` for character (strings)

- `"int"` for integers

- `"logi"` for logical values (TRUE / FALSE)

- `"numeric"` for numeric values

- `"factor"` for factors

- `"date"` for generic date format - note that functions that need to
  read and parse files will try to guess the format and parsing may fail

- One of the accepted date/datetime formats by `lubridate`, you can use
  [`ISAnalytics::date_formats()`](https://calabrialab.github.io/ISAnalytics/dev/reference/date_formats.md)
  to view the accepted formats

- `"_"` to drop the column

For more details see the "How to use import functions" vignette:
[`vignette("workflow_start", package = "ISAnalytics")`](https://calabrialab.github.io/ISAnalytics/dev/articles/workflow_start.md)

### Transformations

Lambdas provided in input in the `transformations` argument, must be
transformations, aka functions that take in input a vector and return a
vector of the same length as the input.

If the transformation list contains column names that are not present in
the data frame, they are simply ignored.

## Required tags

The function will explicitly check for the presence of these tags:

- All columns declared in
  [`mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md)

## See also

[`transform_columns`](https://calabrialab.github.io/ISAnalytics/dev/reference/transform_columns.md)

Other Import functions:
[`import_Vispa2_stats()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_Vispa2_stats.md),
[`import_association_file()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_association_file.md),
[`import_parallel_Vispa2Matrices()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_parallel_Vispa2Matrices.md)

## Examples

``` r
fs_path <- generate_default_folder_structure(type = "correct")
matrix_path <- fs::path(
    fs_path$root, "PJ01", "quantification",
    "POOL01-1", "PJ01_POOL01-1_seqCount_matrix.no0.annotated.tsv.gz"
)
matrix <- import_single_Vispa2Matrix(matrix_path)
#> Reading file...
#> ℹ Mode: fread
#> Reshaping...
#> *** File info *** 
#> • --- Annotated: TRUE
#> • --- Dimensions: 261 x 29
#> • --- Read mode: fread
#> • --- Sample count: 24
head(matrix)
#> # A tibble: 6 × 7
#>   chr   integration_locus strand GeneName     GeneStrand CompleteAmplificationID
#>   <chr>             <int> <chr>  <chr>        <chr>      <chr>                  
#> 1 16             68164148 +      NFATC3       +          PJ01_POOL01_LTR75LC38_…
#> 2 4             129390130 +      LOC100507487 +          PJ01_POOL01_LTR75LC38_…
#> 3 5              84009671 -      EDIL3        -          PJ01_POOL01_LTR75LC38_…
#> 4 12             54635693 -      CBX5         -          PJ01_POOL01_LTR75LC38_…
#> 5 2             181930711 +      UBE2E3       +          PJ01_POOL01_LTR75LC38_…
#> 6 20             35920986 +      MANBAL       +          PJ01_POOL01_LTR75LC38_…
#> # ℹ 1 more variable: Value <int>
```
