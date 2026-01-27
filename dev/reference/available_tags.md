# All available tags for dynamic vars look-up tables.

Contains all information associated with critical tags used in the
dynamic vars system. To know more see
[`vignette("workflow_start", package="ISAnalytics")`](https://calabrialab.github.io/ISAnalytics/dev/articles/workflow_start.md).

## Usage

``` r
available_tags()
```

## Value

A data frame

## Examples

``` r
available_tags()
#> # A tibble: 27 × 4
#>    tag         needed_in description                                dyn_vars_tbl
#>    <chr>       <list>    <chr>                                      <chr>       
#>  1 chromosome  <chr [3]> Number of the chromosome                   mand_vars   
#>  2 locus       <chr [3]> The locus at which the integration occurs  mand_vars   
#>  3 is_strand   <chr [2]> The DNA strand in which the integration o… mand_vars   
#>  4 gene_symbol <chr [5]> The symbol of the gene                     annot_vars  
#>  5 gene_strand <chr [2]> The strand of the gene                     annot_vars  
#>  6 project_id  <chr [6]> Unique identifier of a project             af_vars     
#>  7 pool_id     <chr [3]> Unique identifier of a sequencing pool     af_vars     
#>  8 fusion_id   <chr [1]> Identification code/number of the barcode… af_vars     
#>  9 tag_seq     <chr [3]> The barcode tag sequence                   af_vars     
#> 10 subject     <chr [2]> Unique identifier of a study subject (usu… af_vars     
#> # ℹ 17 more rows
```
