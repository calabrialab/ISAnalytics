# Computes the abundance for every integration event in the input data frame.

**\[stable\]** Abundance is obtained for every integration event by
calculating the ratio between the single value and the total value for
the given group.

## Usage

``` r
compute_abundance(
  x,
  columns = c("fragmentEstimate_sum"),
  percentage = TRUE,
  key = c("SubjectID", "CellMarker", "Tissue", "TimePoint"),
  keep_totals = FALSE
)
```

## Arguments

- x:

  An integration matrix - aka a data frame that includes the
  [`mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md)
  as columns. The matrix can either be aggregated (via
  [`aggregate_values_by_key()`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_values_by_key.md))
  or not.

- columns:

  A character vector of column names to process, must be numeric or
  integer columns

- percentage:

  Add abundance as percentage?

- key:

  The key to group by when calculating totals

- keep_totals:

  A value between `TRUE`, `FALSE` or `df`. If `TRUE`, the intermediate
  totals for each group will be kept in the output data frame as a
  dedicated column with a trailing "\_tot". If `FALSE`, totals won't be
  included in the output data frame. If `df`, the totals are returned to
  the user as a separate data frame, together with the abundance data
  frame.

## Value

Either a single data frame with computed abundance values or a list of 2
data frames (abundance_df, quant_totals)

## Details

Abundance will be computed upon the user selected columns in the
`columns` parameter. For each column a corresponding relative abundance
column (and optionally a percentage abundance column) will be produced.

## Required tags

The function will explicitly check for the presence of these tags:

- All columns declared in
  [`mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md)

## See also

Other Analysis functions:
[`CIS_grubbs()`](https://calabrialab.github.io/ISAnalytics/dev/reference/CIS_grubbs.md),
[`HSC_population_size_estimate()`](https://calabrialab.github.io/ISAnalytics/dev/reference/HSC_population_size_estimate.md),
[`cumulative_is()`](https://calabrialab.github.io/ISAnalytics/dev/reference/cumulative_is.md),
[`gene_frequency_fisher()`](https://calabrialab.github.io/ISAnalytics/dev/reference/gene_frequency_fisher.md),
[`is_sharing()`](https://calabrialab.github.io/ISAnalytics/dev/reference/is_sharing.md),
[`iss_source()`](https://calabrialab.github.io/ISAnalytics/dev/reference/iss_source.md),
[`sample_statistics()`](https://calabrialab.github.io/ISAnalytics/dev/reference/sample_statistics.md),
[`top_integrations()`](https://calabrialab.github.io/ISAnalytics/dev/reference/top_integrations.md),
[`top_targeted_genes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/top_targeted_genes.md)

## Examples

``` r
data("integration_matrices", package = "ISAnalytics")
abund <- compute_abundance(
    x = integration_matrices,
    columns = "fragmentEstimate",
    key = "CompleteAmplificationID"
)
head(abund)
#>       chr integration_locus strand     GeneName GeneStrand
#>    <char>             <num> <char>       <char>     <char>
#> 1:     16          68164148      +       NFATC3          +
#> 2:      4         129390130      + LOC100507487          +
#> 3:      5          84009671      -        EDIL3          -
#> 4:     12          54635693      -         CBX5          -
#> 5:      5          84009671      -        EDIL3          -
#> 6:     12          54635693      -         CBX5          -
#>                                                 CompleteAmplificationID
#>                                                                  <char>
#> 1: PJ01_POOL01_LTR75LC38_PT001_PT001-103_lenti_GLOBE_PB_1_SLiM_0060_MNC
#> 2:  PJ01_POOL01_LTR53LC32_PT001_PT001-81_lenti_GLOBE_BM_1_SLiM_0180_MNC
#> 3:  PJ01_POOL01_LTR53LC32_PT001_PT001-81_lenti_GLOBE_BM_1_SLiM_0180_MNC
#> 4:  PJ01_POOL01_LTR83LC66_PT001_PT001-81_lenti_GLOBE_BM_1_SLiM_0180_MNC
#> 5:  PJ01_POOL01_LTR83LC66_PT001_PT001-81_lenti_GLOBE_BM_1_SLiM_0180_MNC
#> 6:  PJ01_POOL01_LTR27LC94_PT001_PT001-81_lenti_GLOBE_BM_1_SLiM_0180_MNC
#>    seqCount fragmentEstimate fragmentEstimate_RelAbundance
#>       <num>            <num>                         <num>
#> 1:      182        102.94572                     0.5120030
#> 2:    23219         68.73747                     0.1472510
#> 3:    20205         67.12349                     0.1437935
#> 4:    13269         65.15760                     0.1461826
#> 5:    14748         61.46981                     0.1379090
#> 6:    12588         60.84781                     0.1718091
#>    fragmentEstimate_PercAbundance
#>                             <num>
#> 1:                       51.20030
#> 2:                       14.72510
#> 3:                       14.37935
#> 4:                       14.61826
#> 5:                       13.79090
#> 6:                       17.18091
```
