# Default blood lineages info

A default table with info relative to different blood lineages
associated with cell markers that can be supplied as a parameter to
[`HSC_population_size_estimate`](https://calabrialab.github.io/ISAnalytics/dev/reference/HSC_population_size_estimate.md)

## Usage

``` r
blood_lineages_default()
```

## Value

A data frame

## Examples

``` r
blood_lineages_default()
#> # A tibble: 47 × 6
#>    CellMarker Keywords CellType  HematoLineage SuperGroup LineageByPurity
#>    <chr>      <chr>    <chr>     <chr>         <chr>      <chr>          
#>  1 CD13       MYELO    Myeloid   Myeloid       CD13       Myeloid        
#>  2 CD14       MYELO    Myeloid   Myeloid       CD14       Myeloid        
#>  3 CD15       MYELO    Myeloid   Myeloid       CD15       Myeloid        
#>  4 CD19       B        B         Lymphoid      CD19       Lymphoid       
#>  5 CD3        T        T         Lymphoid      CD3        Lymphoid       
#>  6 CD34       CD34     CD34      CD34          CD34       CD34           
#>  7 CD34Molmed CD34     Other     Other         Other      Other          
#>  8 CD34NEG    CD34     Other     Other         Other      Other          
#>  9 CD36       TE       Erythroid Erythroid     CD36       Erythroid      
#> 10 CD4        T        T         Lymphoid      CD3        Lymphoid       
#> # ℹ 37 more rows
```
