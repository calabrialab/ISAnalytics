# Sets the look-up table for matrix file suffixes.

The function automatically produces and sets a look-up table of matrix
file suffixes based on user input.

## Usage

``` r
set_matrix_file_suffixes(
  quantification_suffix = list(seqCount = "seqCount", fragmentEstimate =
    "fragmentEstimate", barcodeCount = "barcodeCount", cellCount = "cellCount", ShsCount
    = "ShsCount"),
  annotation_suffix = list(annotated = ".no0.annotated", not_annotated = ""),
  file_ext = "tsv.gz",
  glue_file_spec = "{quantification_suffix}_matrix{annotation_suffix}.{file_ext}"
)
```

## Arguments

- quantification_suffix:

  A named list - names must be quantification types in
  [`quantification_types()`](https://calabrialab.github.io/ISAnalytics/dev/reference/quantification_types.md),
  and values must be single strings, containing the associated suffix.
  Please note that ALL quantification types must be specified or the
  function will produce an error.

- annotation_suffix:

  A named list - names must be `annotated` and `not_annotated`, values
  must be single strings, containing the associated suffix. Please note
  that both names must be present in the list or the function will
  produce an error.

- file_ext:

  The file extension (e.g. `tsv`, `tsv.gz`)

- glue_file_spec:

  A string specifying the pattern used to form the entire suffix, as per
  [`glue::glue()`](https://glue.tidyverse.org/reference/glue.html)
  requirements. The string should contain the reference to
  `quantification_suffix`, `annotation_suffix` and `file_ext`.

## Value

`NULL`

## See also

Other dynamic vars:
[`inspect_tags()`](https://calabrialab.github.io/ISAnalytics/dev/reference/inspect_tags.md),
[`mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md),
[`pcr_id_column()`](https://calabrialab.github.io/ISAnalytics/dev/reference/pcr_id_column.md),
[`reset_mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/reset_mandatory_IS_vars.md),
[`set_mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/set_mandatory_IS_vars.md)

## Examples

``` r
set_matrix_file_suffixes(
    quantification_suffix = list(
        seqCount = "sc",
        fragmentEstimate = "fe",
        barcodeCount = "barcodeCount",
        cellCount = "cellCount",
        ShsCount = "ShsCount"
    ),
    annotation_suffix = list(annotated = "annot", not_annotated = "")
)
#> Matrix suffixes specs successfully changed
matrix_file_suffixes()
#> # A tibble: 10 × 3
#>    quantification   matrix_type   file_suffix                    
#>    <chr>            <chr>         <chr>                          
#>  1 seqCount         annotated     sc_matrixannot.tsv.gz          
#>  2 seqCount         not_annotated sc_matrix.tsv.gz               
#>  3 fragmentEstimate annotated     fe_matrixannot.tsv.gz          
#>  4 fragmentEstimate not_annotated fe_matrix.tsv.gz               
#>  5 barcodeCount     annotated     barcodeCount_matrixannot.tsv.gz
#>  6 barcodeCount     not_annotated barcodeCount_matrix.tsv.gz     
#>  7 cellCount        annotated     cellCount_matrixannot.tsv.gz   
#>  8 cellCount        not_annotated cellCount_matrix.tsv.gz        
#>  9 ShsCount         annotated     ShsCount_matrixannot.tsv.gz    
#> 10 ShsCount         not_annotated ShsCount_matrix.tsv.gz         
reset_matrix_file_suffixes()
#> Matrix suffixes specs reset to default
```
