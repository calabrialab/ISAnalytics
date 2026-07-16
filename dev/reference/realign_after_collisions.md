# Re-aligns matrices of other quantification types based on the processed sequence count matrix.

**\[stable\]** This function should be used to keep data consistent
among the same analysis: if for some reason you removed the collisions
by passing only the sequence count matrix to
[`remove_collisions()`](https://calabrialab.github.io/ISAnalytics/dev/reference/remove_collisions.md),
you should call this function afterwards, providing a list of other
quantification matrices. NOTE: if you provided a list of several
quantification types to
[`remove_collisions()`](https://calabrialab.github.io/ISAnalytics/dev/reference/remove_collisions.md)
before, there is no need to call this function.

## Usage

``` r
realign_after_collisions(
  sc_matrix,
  other_matrices,
  sample_column = pcr_id_column()
)
```

## Arguments

- sc_matrix:

  The sequence count matrix already processed for collisions via
  [`remove_collisions()`](https://calabrialab.github.io/ISAnalytics/dev/reference/remove_collisions.md)

- other_matrices:

  A named list of matrices to re-align. Names in the list must be
  quantification types
  ([`quantification_types()`](https://calabrialab.github.io/ISAnalytics/dev/reference/quantification_types.md))
  except "seqCount".

- sample_column:

  The name of the column containing the sample identifier

## Value

A named list with re-aligned matrices

## Details

For more details on how to use collision removal functionality:
[`vignette("workflow_start", package = "ISAnalytics")`](https://calabrialab.github.io/ISAnalytics/dev/articles/workflow_start.md)

## See also

[`remove_collisions`](https://calabrialab.github.io/ISAnalytics/dev/reference/remove_collisions.md)

Other Data cleaning and pre-processing:
[`aggregate_metadata()`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_metadata.md),
[`aggregate_values_by_key()`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_values_by_key.md),
[`compute_near_integrations()`](https://calabrialab.github.io/ISAnalytics/dev/reference/compute_near_integrations.md),
[`default_meta_agg()`](https://calabrialab.github.io/ISAnalytics/dev/reference/default_meta_agg.md),
[`outlier_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/outlier_filter.md),
[`outliers_by_pool_fragments()`](https://calabrialab.github.io/ISAnalytics/dev/reference/outliers_by_pool_fragments.md),
[`purity_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/purity_filter.md),
[`remove_collisions()`](https://calabrialab.github.io/ISAnalytics/dev/reference/remove_collisions.md),
[`threshold_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/threshold_filter.md)

## Examples

``` r
data("integration_matrices", package = "ISAnalytics")
data("association_file", package = "ISAnalytics")
separated <- separate_quant_matrices(
    integration_matrices
)
no_coll <- remove_collisions(
    x = separated$seqCount,
    association_file = association_file,
    quant_cols = c(seqCount = "Value"),
    report_path = NULL
)
#> Identifying collisions...
#> Processing collisions...
#> Finished!
realigned <- realign_after_collisions(
    sc_matrix = no_coll,
    other_matrices = list(fragmentEstimate = separated$fragmentEstimate)
)
realigned
#> $fragmentEstimate
#>          chr integration_locus strand     GeneName GeneStrand
#>       <char>             <num> <char>       <char>     <char>
#>    1:     16          68164148      +       NFATC3          +
#>    2:      4         129390130      + LOC100507487          +
#>    3:      5          84009671      -        EDIL3          -
#>    4:     12          54635693      -         CBX5          -
#>    5:      5          84009671      -        EDIL3          -
#>   ---                                                        
#> 1662:      6           3388625      -         HTR4          +
#> 1663:     16           3207754      +       UBE2D2          +
#> 1664:     19          13631664      -    LINC01133          +
#> 1665:     12          48119158      -        KMT2D          -
#> 1666:     14          30726412      -     PLEKHG4B          -
#>                                                    CompleteAmplificationID
#>                                                                     <char>
#>    1: PJ01_POOL01_LTR75LC38_PT001_PT001-103_lenti_GLOBE_PB_1_SLiM_0060_MNC
#>    2:  PJ01_POOL01_LTR53LC32_PT001_PT001-81_lenti_GLOBE_BM_1_SLiM_0180_MNC
#>    3:  PJ01_POOL01_LTR53LC32_PT001_PT001-81_lenti_GLOBE_BM_1_SLiM_0180_MNC
#>    4:  PJ01_POOL01_LTR83LC66_PT001_PT001-81_lenti_GLOBE_BM_1_SLiM_0180_MNC
#>    5:  PJ01_POOL01_LTR83LC66_PT001_PT001-81_lenti_GLOBE_BM_1_SLiM_0180_MNC
#>   ---                                                                     
#> 1662: PJ01_POOL03_LTR93LC90_PT002_PT002-464_lenti_GLOBE_PB_1_SLiM_0360_MNC
#> 1663: PJ01_POOL03_LTR93LC90_PT002_PT002-464_lenti_GLOBE_PB_1_SLiM_0360_MNC
#> 1664: PJ01_POOL03_LTR93LC90_PT002_PT002-464_lenti_GLOBE_PB_1_SLiM_0360_MNC
#> 1665: PJ01_POOL03_LTR93LC90_PT002_PT002-464_lenti_GLOBE_PB_1_SLiM_0360_MNC
#> 1666: PJ01_POOL03_LTR93LC90_PT002_PT002-464_lenti_GLOBE_PB_1_SLiM_0360_MNC
#>            Value
#>            <num>
#>    1: 102.945718
#>    2:  68.737467
#>    3:  67.123486
#>    4:  65.157600
#>    5:  61.469810
#>   ---           
#> 1662:  11.761600
#> 1663:  12.600475
#> 1664:   1.704548
#> 1665:  11.366729
#> 1666:   6.047534
#> 
```
