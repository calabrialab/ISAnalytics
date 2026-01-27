# Names of the columns of the association file to consider for Vispa2 launch.

Selection of column names from the association file to be considered for
Vispa2 launch. NOTE: the `TagID` column appears only once but needs to
be repeated twice for generating the launch file. Use the appropriate
function to generate the file automatically.

## Usage

``` r
reduced_AF_columns()
```

## Value

A character vector

## Examples

``` r
reduced_AF_columns()
#> # A tibble: 10 × 2
#>    names                   tag        
#>    <chr>                   <chr>      
#>  1 CellMarker              cell_marker
#>  2 FUSIONID                fusion_id  
#>  3 CompleteAmplificationID pcr_repl_id
#>  4 PoolID                  pool_id    
#>  5 ProjectID               project_id 
#>  6 SubjectID               subject    
#>  7 TagID                   tag_id     
#>  8 Tissue                  tissue     
#>  9 TimePoint               tp_days    
#> 10 VectorID                vector_id  
```
