# Computes user specified functions on numerical columns and updates the metadata data frame accordingly.

**\[stable\]** The function operates on a data frame by grouping the
content by the sample key and computing every function specified on
every column in the `value_columns` parameter. After that the metadata
data frame is updated by including the computed results as columns for
the corresponding key. For this reason it's required that both `x` and
`metadata` have the same sample key, and it's particularly important if
the user is working with previously aggregated data. For example:

    data("integration_matrices", package = "ISAnalytics")
    data("association_file", package = "ISAnalytics")
    aggreg <- aggregate_values_by_key(
     x = integration_matrices,
     association_file = association_file,
     value_cols = c("seqCount", "fragmentEstimate")
    )
    aggreg_meta <- aggregate_metadata(association_file = association_file)

    sample_stats <- sample_statistics(x = aggreg,
    metadata = aggreg_meta,
    value_columns = c("seqCount", "fragmentEstimate"),
    sample_key = c("SubjectID", "CellMarker","Tissue", "TimePoint"))

## Usage

``` r
sample_statistics(
  x,
  metadata,
  sample_key = "CompleteAmplificationID",
  value_columns = "Value",
  functions = default_stats(),
  add_integrations_count = TRUE
)
```

## Arguments

- x:

  A data frame

- metadata:

  The metadata data frame

- sample_key:

  Character vector representing the key for identifying a sample

- value_columns:

  The name of the columns to be computed, must be numeric or integer

- functions:

  A named list of function or purrr-style lambdas

- add_integrations_count:

  Add the count of distinct integration sites for each group? Can be
  computed only if `x` contains the mandatory columns
  [`mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md)

## Value

A list with modified x and metadata data frames

## Required tags

The function will explicitly check for the presence of these tags:

- All columns declared in
  [`mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md)

These are checked only if `add_integrations_count = TRUE`.

## See also

Other Analysis functions:
[`CIS_grubbs()`](https://calabrialab.github.io/ISAnalytics/dev/reference/CIS_grubbs.md),
[`HSC_population_size_estimate()`](https://calabrialab.github.io/ISAnalytics/dev/reference/HSC_population_size_estimate.md),
[`compute_abundance()`](https://calabrialab.github.io/ISAnalytics/dev/reference/compute_abundance.md),
[`cumulative_is()`](https://calabrialab.github.io/ISAnalytics/dev/reference/cumulative_is.md),
[`gene_frequency_fisher()`](https://calabrialab.github.io/ISAnalytics/dev/reference/gene_frequency_fisher.md),
[`is_sharing()`](https://calabrialab.github.io/ISAnalytics/dev/reference/is_sharing.md),
[`iss_source()`](https://calabrialab.github.io/ISAnalytics/dev/reference/iss_source.md),
[`top_integrations()`](https://calabrialab.github.io/ISAnalytics/dev/reference/top_integrations.md),
[`top_targeted_genes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/top_targeted_genes.md)

## Examples

``` r
data("integration_matrices", package = "ISAnalytics")
data("association_file", package = "ISAnalytics")
stats <- sample_statistics(
    x = integration_matrices,
    metadata = association_file,
    value_columns = c("seqCount", "fragmentEstimate")
)
stats
#> $x
#> # A tibble: 53 × 38
#>    CompleteAmplificationID          seqCount_sum seqCount_count seqCount_shannon
#>    <chr>                                   <dbl>          <int>            <dbl>
#>  1 PJ01_POOL01_LTR19LC2_PT001_PT00…         3029             22             1.82
#>  2 PJ01_POOL01_LTR27LC94_PT001_PT0…        80610             67             3.02
#>  3 PJ01_POOL01_LTR29LC90_PT001_PT0…           30              8             1.89
#>  4 PJ01_POOL01_LTR37LC18_PT001_PT0…           44             10             2.00
#>  5 PJ01_POOL01_LTR37LC2_PT001_PT00…        50859             82             3.76
#>  6 PJ01_POOL01_LTR49LC12_PT001_PT0…           33             10             1.93
#>  7 PJ01_POOL01_LTR53LC22_PT001_PT0…            1              1             0   
#>  8 PJ01_POOL01_LTR53LC32_PT001_PT0…       129867             72             3.24
#>  9 PJ01_POOL01_LTR57LC20_PT001_PT0…        57443             20             2.41
#> 10 PJ01_POOL01_LTR5LC64_PT001_PT00…        41488             39             2.89
#> # ℹ 43 more rows
#> # ℹ 34 more variables: seqCount_simpson <dbl>, seqCount_invsimpson <dbl>,
#> #   fragmentEstimate_sum <dbl>, fragmentEstimate_count <int>,
#> #   fragmentEstimate_shannon <dbl>, fragmentEstimate_simpson <dbl>,
#> #   fragmentEstimate_invsimpson <dbl>, seqCount_describe_vars <dbl>,
#> #   seqCount_describe_n <dbl>, seqCount_describe_mean <dbl>,
#> #   seqCount_describe_sd <dbl>, seqCount_describe_median <dbl>, …
#> 
#> $metadata
#>     ProjectID  FUSIONID PoolID TagSequence SubjectID VectorType VectorID
#>        <char>    <char> <char>      <char>    <char>     <char>   <char>
#>  1:      PJ01 ET#382.46 POOL01   LTR75LC38     PT001      lenti    GLOBE
#>  2:      PJ01 ET#381.40 POOL01   LTR53LC32     PT001      lenti    GLOBE
#>  3:      PJ01  ET#381.9 POOL01   LTR83LC66     PT001      lenti    GLOBE
#>  4:      PJ01 ET#381.71 POOL01   LTR27LC94     PT001      lenti    GLOBE
#>  5:      PJ01  ET#381.2 POOL01   LTR69LC52     PT001      lenti    GLOBE
#>  6:      PJ01 ET#382.28 POOL01    LTR37LC2     PT001      lenti    GLOBE
#>  7:      PJ01  ET#382.2 POOL01   LTR77LC46     PT001      lenti    GLOBE
#>  8:      PJ01 ET#382.50 POOL01   LTR83LC46     PT001      lenti    GLOBE
#>  9:      PJ01 ET#381.21 POOL01    LTR9LC90     PT001      lenti    GLOBE
#> 10:      PJ01 ET#407.28 POOL02   LTR65LC56     PT001      lenti    GLOBE
#> 11:      PJ01 ET#381.33 POOL01   LTR37LC18     PT001      lenti    GLOBE
#> 12:      PJ01 ET#382.15 POOL01    LTR7LC72     PT001      lenti    GLOBE
#> 13:      PJ01  ET#382.6 POOL01   LTR85LC54     PT001      lenti    GLOBE
#> 14:      PJ01 ET#381.52 POOL01   LTR77LC56     PT001      lenti    GLOBE
#> 15:      PJ01 ET#382.55 POOL01   LTR93LC56     PT001      lenti    GLOBE
#> 16:      PJ01 ET#381.25 POOL01    LTR19LC2     PT001      lenti    GLOBE
#> 17:      PJ01 ET#382.11 POOL01   LTR95LC64     PT001      lenti    GLOBE
#> 18:      PJ01 ET#382.33 POOL01   LTR49LC12     PT001      lenti    GLOBE
#> 19:      PJ01 ET#382.37 POOL01   LTR57LC20     PT001      lenti    GLOBE
#> 20:      PJ01 ET#382.59 POOL01    LTR5LC64     PT001      lenti    GLOBE
#> 21:      PJ01 ET#381.56 POOL01   LTR85LC64     PT001      lenti    GLOBE
#> 22:      PJ01 ET#381.87 POOL01   LTR61LC30     PT001      lenti    GLOBE
#> 23:      PJ01 ET#408.37 POOL02   LTR87LC74     PT001      lenti    GLOBE
#> 24:      PJ01 ET#382.24 POOL01   LTR29LC90     PT001      lenti    GLOBE
#> 25:      PJ01 ET#381.83 POOL01   LTR53LC22     PT001      lenti    GLOBE
#> 26:      PJ01 ET#381.64 POOL01    LTR9LC80     PT001      lenti    GLOBE
#> 27:      PJ01  FB585.26 POOL03   LTR73LC86     PT002      lenti    GLOBE
#> 28:      PJ01 ET#408.58 POOL04   LTR41LC20     PT002      lenti    GLOBE
#> 29:      PJ01 ET#408.27 POOL04   LTR65LC54     PT002      lenti    GLOBE
#> 30:      PJ01   FB585.4 POOL03   LTR51LC86     PT002      lenti    GLOBE
#> 31:      PJ01 ET#414.31 POOL04   LTR85LC62     PT002      lenti    GLOBE
#> 32:      PJ01 ET#408.89 POOL04   LTR11LC82     PT002      lenti    GLOBE
#> 33:      PJ01  FB585.24 POOL03   LTR71LC90     PT002      lenti    GLOBE
#> 34:      PJ01 ET#414.16 POOL04   LTR51LC32     PT002      lenti    GLOBE
#> 35:      PJ01  FB585.46 POOL03   LTR93LC90     PT002      lenti    GLOBE
#> 36:      PJ01 ET#411.42 POOL04    LTR5LC84     PT002      lenti    GLOBE
#> 37:      PJ01  ET#414.1 POOL04    LTR13LC2     PT002      lenti    GLOBE
#> 38:      PJ01 ET#412.41 POOL04    LTR5LC82     PT002      lenti    GLOBE
#> 39:      PJ01 ET#411.49 POOL04    LTR25LC2     PT002      lenti    GLOBE
#> 40:      PJ01 ET#411.18 POOL04   LTR49LC36     PT002      lenti    GLOBE
#> 41:      PJ01 ET#412.23 POOL04   LTR63LC46     PT002      lenti    GLOBE
#> 42:      PJ01 ET#411.80 POOL04   LTR93LC64     PT002      lenti    GLOBE
#> 43:      PJ01   FB585.2 POOL03   LTR49LC90     PT002      lenti    GLOBE
#> 44:      PJ01  FB585.48 POOL03   LTR95LC86     PT002      lenti    GLOBE
#> 45:      PJ01 ET#412.72 POOL04   LTR79LC48     PT002      lenti    GLOBE
#> 46:      PJ01 ET#412.65 POOL04   LTR63LC34     PT002      lenti    GLOBE
#> 47:      PJ01 ET#412.85 POOL04    LTR9LC74     PT002      lenti    GLOBE
#> 48:      PJ01 ET#412.10 POOL04   LTR33LC20     PT002      lenti    GLOBE
#> 49:      PJ01 ET#408.82 POOL04   LTR93LC68     PT002      lenti    GLOBE
#> 50:      PJ01  ET#412.3 POOL04    LTR13LC6     PT002      lenti    GLOBE
#> 51:      PJ01 ET#412.54 POOL04   LTR37LC12     PT002      lenti    GLOBE
#> 52:      PJ01 ET#412.34 POOL04   LTR87LC68     PT002      lenti    GLOBE
#> 53:      PJ01 ET#408.51 POOL04    LTR25LC6     PT002      lenti    GLOBE
#>     ProjectID  FUSIONID PoolID TagSequence SubjectID VectorType VectorID
#>        <char>    <char> <char>      <char>    <char>     <char>   <char>
#>     ExperimentID Tissue TimePoint DNAFragmentation PCRMethod TagIDextended
#>           <char> <char>    <char>           <char>    <char>        <char>
#>  1:         <NA>     PB      0060            SONIC      SLiM     LTR75LC38
#>  2:         <NA>     BM      0180            SONIC      SLiM     LTR53LC32
#>  3:         <NA>     BM      0180            SONIC      SLiM     LTR83LC66
#>  4:         <NA>     BM      0180            SONIC      SLiM     LTR27LC94
#>  5:         <NA>     PB      0180            SONIC      SLiM     LTR69LC52
#>  6:         <NA>     BM      0060            SONIC      SLiM      LTR37LC2
#>  7:         <NA>     PB      0060            SONIC      SLiM     LTR77LC46
#>  8:         <NA>     BM      0060            SONIC      SLiM     LTR83LC46
#>  9:         <NA>     PB      0030            SONIC      SLiM      LTR9LC90
#> 10:         <NA>     PB      0360            SONIC      SLiM     LTR65LC56
#> 11:         <NA>     PB      0180            SONIC      SLiM     LTR37LC18
#> 12:         <NA>     BM      0090            SONIC      SLiM      LTR7LC72
#> 13:         <NA>     BM      0060            SONIC      SLiM     LTR85LC54
#> 14:         <NA>     PB      0030            SONIC      SLiM     LTR77LC56
#> 15:         <NA>     PB      0090            SONIC      SLiM     LTR93LC56
#> 16:         <NA>     BM      0030            SONIC      SLiM      LTR19LC2
#> 17:         <NA>     PB      0090            SONIC      SLiM     LTR95LC64
#> 18:         <NA>     PB      0090            SONIC      SLiM     LTR49LC12
#> 19:         <NA>     BM      0090            SONIC      SLiM     LTR57LC20
#> 20:         <NA>     BM      0090            SONIC      SLiM      LTR5LC64
#> 21:         <NA>     BM      0030            SONIC      SLiM     LTR85LC64
#> 22:         <NA>     BM      0030            SONIC      SLiM     LTR61LC30
#> 23:         <NA>     BM      0360            SONIC      SLiM     LTR87LC74
#> 24:         <NA>     PB      0060            SONIC      SLiM     LTR29LC90
#> 25:         <NA>     PB      0030            SONIC      SLiM     LTR53LC22
#> 26:         <NA>     PB      0180            SONIC      SLiM      LTR9LC80
#> 27:         <NA>     BM      0360            SONIC      SLiM     LTR73LC86
#> 28:         <NA>     BM      0030            SONIC      SLiM     LTR41LC20
#> 29:         <NA>     BM      0030            SONIC      SLiM     LTR65LC54
#> 30:         <NA>     BM      0360            SONIC      SLiM     LTR51LC86
#> 31:         <NA>     BM      0180            SONIC      SLiM     LTR85LC62
#> 32:         <NA>     BM      0030            SONIC      SLiM     LTR11LC82
#> 33:         <NA>     PB      0360            SONIC      SLiM     LTR71LC90
#> 34:         <NA>     BM      0180            SONIC      SLiM     LTR51LC32
#> 35:         <NA>     PB      0360            SONIC      SLiM     LTR93LC90
#> 36:         <NA>     PB      0060            SONIC      SLiM      LTR5LC84
#> 37:         <NA>     BM      0180            SONIC      SLiM      LTR13LC2
#> 38:         <NA>     BM      0090            SONIC      SLiM      LTR5LC82
#> 39:         <NA>     BM      0060            SONIC      SLiM      LTR25LC2
#> 40:         <NA>     BM      0060            SONIC      SLiM     LTR49LC36
#> 41:         <NA>     PB      0180            SONIC      SLiM     LTR63LC46
#> 42:         <NA>     BM      0060            SONIC      SLiM     LTR93LC64
#> 43:         <NA>     PB      0360            SONIC      SLiM     LTR49LC90
#> 44:         <NA>     BM      0360            SONIC      SLiM     LTR95LC86
#> 45:         <NA>     BM      0090            SONIC      SLiM     LTR79LC48
#> 46:         <NA>     PB      0090            SONIC      SLiM     LTR63LC34
#> 47:         <NA>     PB      0180            SONIC      SLiM      LTR9LC74
#> 48:         <NA>     BM      0090            SONIC      SLiM     LTR33LC20
#> 49:         <NA>     PB      0030            SONIC      SLiM     LTR93LC68
#> 50:         <NA>     PB      0090            SONIC      SLiM      LTR13LC6
#> 51:         <NA>     PB      0180            SONIC      SLiM     LTR37LC12
#> 52:         <NA>     PB      0090            SONIC      SLiM     LTR87LC68
#> 53:         <NA>     PB      0030            SONIC      SLiM      LTR25LC6
#>     ExperimentID Tissue TimePoint DNAFragmentation PCRMethod TagIDextended
#>           <char> <char>    <char>           <char>    <char>        <char>
#>     Keywords CellMarker      TagID NGSProvider NGSTechnology ConverrtedFilesDir
#>       <char>     <char>     <char>      <char>        <char>             <char>
#>  1:     <NA>        MNC LTR75.LC38        <NA>         HiSeq               <NA>
#>  2:     <NA>        MNC LTR53.LC32        <NA>         HiSeq               <NA>
#>  3:     <NA>        MNC LTR83.LC66        <NA>         HiSeq               <NA>
#>  4:     <NA>        MNC LTR27.LC94        <NA>         HiSeq               <NA>
#>  5:     <NA>        MNC LTR69.LC52        <NA>         HiSeq               <NA>
#>  6:     <NA>        MNC  LTR37.LC2        <NA>         HiSeq               <NA>
#>  7:     <NA>        MNC LTR77.LC46        <NA>         HiSeq               <NA>
#>  8:     <NA>        MNC LTR83.LC46        <NA>         HiSeq               <NA>
#>  9:     <NA>        MNC  LTR9L.C90        <NA>         HiSeq               <NA>
#> 10:     <NA>        MNC LTR65.LC56        <NA>         HiSeq               <NA>
#> 11:     <NA>        MNC LTR37.LC18        <NA>         HiSeq               <NA>
#> 12:     <NA>        MNC  LTR7L.C72        <NA>         HiSeq               <NA>
#> 13:     <NA>        MNC LTR85.LC54        <NA>         HiSeq               <NA>
#> 14:     <NA>        MNC LTR77.LC56        <NA>         HiSeq               <NA>
#> 15:     <NA>        MNC LTR93.LC56        <NA>         HiSeq               <NA>
#> 16:     <NA>        MNC  LTR19.LC2        <NA>         HiSeq               <NA>
#> 17:     <NA>        MNC LTR95.LC64        <NA>         HiSeq               <NA>
#> 18:     <NA>        MNC LTR49.LC12        <NA>         HiSeq               <NA>
#> 19:     <NA>        MNC LTR57.LC20        <NA>         HiSeq               <NA>
#> 20:     <NA>        MNC  LTR5L.C64        <NA>         HiSeq               <NA>
#> 21:     <NA>        MNC LTR85.LC64        <NA>         HiSeq               <NA>
#> 22:     <NA>        MNC LTR61.LC30        <NA>         HiSeq               <NA>
#> 23:     <NA>        MNC LTR87.LC74        <NA>         HiSeq               <NA>
#> 24:     <NA>        MNC LTR29.LC90        <NA>         HiSeq               <NA>
#> 25:     <NA>        MNC LTR53.LC22        <NA>         HiSeq               <NA>
#> 26:     <NA>        MNC  LTR9L.C80        <NA>         HiSeq               <NA>
#> 27:     <NA>        MNC LTR73.LC86        <NA>         HiSeq               <NA>
#> 28:     <NA>        MNC LTR41.LC20        <NA>         HiSeq               <NA>
#> 29:     <NA>        MNC LTR65.LC54        <NA>         HiSeq               <NA>
#> 30:     <NA>        MNC LTR51.LC86        <NA>         HiSeq               <NA>
#> 31:     <NA>        MNC LTR85.LC62        <NA>         HiSeq               <NA>
#> 32:     <NA>        MNC LTR11.LC82        <NA>         HiSeq               <NA>
#> 33:     <NA>        MNC LTR71.LC90        <NA>         HiSeq               <NA>
#> 34:     <NA>        MNC LTR51.LC32        <NA>         HiSeq               <NA>
#> 35:     <NA>        MNC LTR93.LC90        <NA>         HiSeq               <NA>
#> 36:     <NA>        MNC  LTR5L.C84        <NA>         HiSeq               <NA>
#> 37:     <NA>        MNC  LTR13.LC2        <NA>         HiSeq               <NA>
#> 38:     <NA>        MNC  LTR5L.C82        <NA>         HiSeq               <NA>
#> 39:     <NA>        MNC  LTR25.LC2        <NA>         HiSeq               <NA>
#> 40:     <NA>        MNC LTR49.LC36        <NA>         HiSeq               <NA>
#> 41:     <NA>        MNC LTR63.LC46        <NA>         HiSeq               <NA>
#> 42:     <NA>        MNC LTR93.LC64        <NA>         HiSeq               <NA>
#> 43:     <NA>        MNC LTR49.LC90        <NA>         HiSeq               <NA>
#> 44:     <NA>        MNC LTR95.LC86        <NA>         HiSeq               <NA>
#> 45:     <NA>        MNC LTR79.LC48        <NA>         HiSeq               <NA>
#> 46:     <NA>        MNC LTR63.LC34        <NA>         HiSeq               <NA>
#> 47:     <NA>        MNC  LTR9L.C74        <NA>         HiSeq               <NA>
#> 48:     <NA>        MNC LTR33.LC20        <NA>         HiSeq               <NA>
#> 49:     <NA>        MNC LTR93.LC68        <NA>         HiSeq               <NA>
#> 50:     <NA>        MNC  LTR13.LC6        <NA>         HiSeq               <NA>
#> 51:     <NA>        MNC LTR37.LC12        <NA>         HiSeq               <NA>
#> 52:     <NA>        MNC LTR87.LC68        <NA>         HiSeq               <NA>
#> 53:     <NA>        MNC  LTR25.LC6        <NA>         HiSeq               <NA>
#>     Keywords CellMarker      TagID NGSProvider NGSTechnology ConverrtedFilesDir
#>       <char>     <char>     <char>      <char>        <char>             <char>
#>     ConverrtedFilesName SourceFileFolder SourceFileNameR1 SourceFileNameR2
#>                  <char>           <char>           <char>           <char>
#>  1:                <NA>             <NA>             <NA>             <NA>
#>  2:                <NA>             <NA>             <NA>             <NA>
#>  3:                <NA>             <NA>             <NA>             <NA>
#>  4:                <NA>             <NA>             <NA>             <NA>
#>  5:                <NA>             <NA>             <NA>             <NA>
#>  6:                <NA>             <NA>             <NA>             <NA>
#>  7:                <NA>             <NA>             <NA>             <NA>
#>  8:                <NA>             <NA>             <NA>             <NA>
#>  9:                <NA>             <NA>             <NA>             <NA>
#> 10:                <NA>             <NA>             <NA>             <NA>
#> 11:                <NA>             <NA>             <NA>             <NA>
#> 12:                <NA>             <NA>             <NA>             <NA>
#> 13:                <NA>             <NA>             <NA>             <NA>
#> 14:                <NA>             <NA>             <NA>             <NA>
#> 15:                <NA>             <NA>             <NA>             <NA>
#> 16:                <NA>             <NA>             <NA>             <NA>
#> 17:                <NA>             <NA>             <NA>             <NA>
#> 18:                <NA>             <NA>             <NA>             <NA>
#> 19:                <NA>             <NA>             <NA>             <NA>
#> 20:                <NA>             <NA>             <NA>             <NA>
#> 21:                <NA>             <NA>             <NA>             <NA>
#> 22:                <NA>             <NA>             <NA>             <NA>
#> 23:                <NA>             <NA>             <NA>             <NA>
#> 24:                <NA>             <NA>             <NA>             <NA>
#> 25:                <NA>             <NA>             <NA>             <NA>
#> 26:                <NA>             <NA>             <NA>             <NA>
#> 27:                <NA>             <NA>             <NA>             <NA>
#> 28:                <NA>             <NA>             <NA>             <NA>
#> 29:                <NA>             <NA>             <NA>             <NA>
#> 30:                <NA>             <NA>             <NA>             <NA>
#> 31:                <NA>             <NA>             <NA>             <NA>
#> 32:                <NA>             <NA>             <NA>             <NA>
#> 33:                <NA>             <NA>             <NA>             <NA>
#> 34:                <NA>             <NA>             <NA>             <NA>
#> 35:                <NA>             <NA>             <NA>             <NA>
#> 36:                <NA>             <NA>             <NA>             <NA>
#> 37:                <NA>             <NA>             <NA>             <NA>
#> 38:                <NA>             <NA>             <NA>             <NA>
#> 39:                <NA>             <NA>             <NA>             <NA>
#> 40:                <NA>             <NA>             <NA>             <NA>
#> 41:                <NA>             <NA>             <NA>             <NA>
#> 42:                <NA>             <NA>             <NA>             <NA>
#> 43:                <NA>             <NA>             <NA>             <NA>
#> 44:                <NA>             <NA>             <NA>             <NA>
#> 45:                <NA>             <NA>             <NA>             <NA>
#> 46:                <NA>             <NA>             <NA>             <NA>
#> 47:                <NA>             <NA>             <NA>             <NA>
#> 48:                <NA>             <NA>             <NA>             <NA>
#> 49:                <NA>             <NA>             <NA>             <NA>
#> 50:                <NA>             <NA>             <NA>             <NA>
#> 51:                <NA>             <NA>             <NA>             <NA>
#> 52:                <NA>             <NA>             <NA>             <NA>
#> 53:                <NA>             <NA>             <NA>             <NA>
#>     ConverrtedFilesName SourceFileFolder SourceFileNameR1 SourceFileNameR2
#>                  <char>           <char>           <char>           <char>
#>     DNAnumber ReplicateNumber DNAextractionDate DNAngUsed LinearPCRID
#>        <char>           <int>            <Date>     <num>      <char>
#>  1: PT001-103               3        2016-03-16   23.1840        <NA>
#>  2:  PT001-81               2        2016-07-15  181.4400        <NA>
#>  3:  PT001-81               1        2016-07-15  181.4400        <NA>
#>  4:  PT001-81               3        2016-07-15  181.4400        <NA>
#>  5:  PT001-74               1        2016-07-15   23.0580        <NA>
#>  6: PT001-107               2        2016-03-16  171.3600        <NA>
#>  7: PT001-103               1        2016-03-16   23.1840        <NA>
#>  8: PT001-107               3        2016-03-16  171.3600        <NA>
#>  9:  PT001-93               1        2016-02-16   23.8500        <NA>
#> 10: PT001-149               1        2017-01-13   45.1500        <NA>
#> 11:  PT001-74               2        2016-07-15   23.0580        <NA>
#> 12: PT001-116               1        2016-04-15   89.2080        <NA>
#> 13: PT001-107               1        2016-03-16  171.3600        <NA>
#> 14:  PT001-93               2        2016-02-16   23.8500        <NA>
#> 15: PT001-112               3        2016-04-15   26.5680        <NA>
#> 16:  PT001-97               1        2016-02-16  300.0327        <NA>
#> 17: PT001-112               1        2016-04-15   26.5680        <NA>
#> 18: PT001-112               2        2016-04-15   26.5680        <NA>
#> 19: PT001-116               2        2016-04-15   89.2080        <NA>
#> 20: PT001-116               3        2016-04-15   89.2080        <NA>
#> 21:  PT001-97               2        2016-02-16  300.0327        <NA>
#> 22:  PT001-97               3        2016-02-16  300.0327        <NA>
#> 23: PT001-156               2        2017-01-13   42.0000        <NA>
#> 24: PT001-103               2        2016-03-16   23.1840        <NA>
#> 25:  PT001-93               3        2016-02-16   23.8500        <NA>
#> 26:  PT001-74               3        2016-07-15   23.0580        <NA>
#> 27: PT002-466               2        2017-05-24  299.3220        <NA>
#> 28: PT002-177               2        2016-07-11  181.5000        <NA>
#> 29: PT002-177               1        2016-07-11  181.5000        <NA>
#> 30: PT002-466               1        2017-05-24  299.3220        <NA>
#> 31: PT002-238               3        2016-12-09  196.5000        <NA>
#> 32: PT002-177               3        2016-07-11  181.5000        <NA>
#> 33: PT002-464               2        2017-05-24  300.4800        <NA>
#> 34: PT002-238               2        2016-12-09  196.5000        <NA>
#> 35: PT002-464               3        2017-05-24  300.4800        <NA>
#> 36: PT002-190               2        2016-08-10   77.4000        <NA>
#> 37: PT002-238               1        2016-12-09  196.5000        <NA>
#> 38: PT002-218               2        2016-09-09  135.1500        <NA>
#> 39: PT002-197               2        2016-08-10  172.5000        <NA>
#> 40: PT002-197               1        2016-08-10  172.5000        <NA>
#> 41: PT002-231               1        2016-12-14  122.5500        <NA>
#> 42: PT002-197               3        2016-08-10  172.5000        <NA>
#> 43: PT002-464               1        2017-05-24  300.4800        <NA>
#> 44: PT002-466               3        2017-05-24  299.3220        <NA>
#> 45: PT002-218               3        2016-09-09  135.1500        <NA>
#> 46: PT002-211               3        2016-09-09  123.0000        <NA>
#> 47: PT002-231               3        2016-12-14  122.5500        <NA>
#> 48: PT002-218               1        2016-09-09  135.1500        <NA>
#> 49: PT002-170               3        2016-07-12  106.3500        <NA>
#> 50: PT002-211               1        2016-09-09  123.0000        <NA>
#> 51: PT002-231               2        2016-12-14  122.5500        <NA>
#> 52: PT002-211               2        2016-09-09  123.0000        <NA>
#> 53: PT002-170               2        2016-07-12  106.3500        <NA>
#>     DNAnumber ReplicateNumber DNAextractionDate DNAngUsed LinearPCRID
#>        <char>           <int>            <Date>     <num>      <char>
#>     LinearPCRDate SonicationDate LigationDate 1stExpoPCRID 1stExpoPCRDate
#>            <Date>         <Date>       <Date>       <char>         <Date>
#>  1:          <NA>     2016-11-02   2016-11-02    ET#380.46     2016-11-02
#>  2:          <NA>     2016-11-02   2016-11-02    ET#379.40     2016-11-02
#>  3:          <NA>     2016-11-02   2016-11-02     ET#379.9     2016-11-02
#>  4:          <NA>     2016-11-02   2016-11-02    ET#379.71     2016-11-02
#>  5:          <NA>     2016-11-02   2016-11-02     ET#379.2     2016-11-02
#>  6:          <NA>     2016-11-02   2016-11-02    ET#380.28     2016-11-02
#>  7:          <NA>     2016-11-02   2016-11-02     ET#380.2     2016-11-02
#>  8:          <NA>     2016-11-02   2016-11-02    ET#380.50     2016-11-02
#>  9:          <NA>     2016-11-02   2016-11-02    ET#379.21     2016-11-02
#> 10:          <NA>     2017-04-19   2017-04-19    ET#405.28     2017-04-20
#> 11:          <NA>     2016-11-02   2016-11-02    ET#379.33     2016-11-02
#> 12:          <NA>     2016-11-02   2016-11-02    ET#380.15     2016-11-02
#> 13:          <NA>     2016-11-02   2016-11-02     ET#380.6     2016-11-02
#> 14:          <NA>     2016-11-02   2016-11-02    ET#379.52     2016-11-02
#> 15:          <NA>     2016-11-02   2016-11-02    ET#380.55     2016-11-02
#> 16:          <NA>     2016-11-02   2016-11-02    ET#379.25     2016-11-02
#> 17:          <NA>     2016-11-02   2016-11-02    ET#380.11     2016-11-02
#> 18:          <NA>     2016-11-02   2016-11-02    ET#380.33     2016-11-02
#> 19:          <NA>     2016-11-02   2016-11-02    ET#380.37     2016-11-02
#> 20:          <NA>     2016-11-02   2016-11-02    ET#380.59     2016-11-02
#> 21:          <NA>     2016-11-02   2016-11-02    ET#379.56     2016-11-02
#> 22:          <NA>     2016-11-02   2016-11-02    ET#379.87     2016-11-02
#> 23:          <NA>     2017-04-19   2017-04-19    ET#406.37     2017-04-20
#> 24:          <NA>     2016-11-02   2016-11-02    ET#380.24     2016-11-02
#> 25:          <NA>     2016-11-02   2016-11-02    ET#379.83     2016-11-02
#> 26:          <NA>     2016-11-02   2016-11-02    ET#379.64     2016-11-02
#> 27:          <NA>     2018-03-12   2018-03-12     FB584.26     2018-03-12
#> 28:          <NA>     2017-04-19   2017-04-19    ET#406.58     2017-04-20
#> 29:          <NA>     2017-04-19   2017-04-19    ET#406.27     2017-04-20
#> 30:          <NA>     2018-03-12   2018-03-12      FB584.4     2018-03-12
#> 31:          <NA>     2017-05-12   2017-05-12    ET#413.31     2017-05-12
#> 32:          <NA>     2017-04-19   2017-04-19    ET#406.89     2017-04-20
#> 33:          <NA>     2018-03-12   2018-03-12     FB584.24     2018-03-12
#> 34:          <NA>     2017-05-12   2017-05-12    ET#413.16     2017-05-12
#> 35:          <NA>     2018-03-12   2018-03-12     FB584.46     2018-03-12
#> 36:          <NA>     2017-05-03   2017-05-03    ET#409.42     2017-05-04
#> 37:          <NA>     2017-05-12   2017-05-12     ET#413.1     2017-05-12
#> 38:          <NA>     2017-05-03   2017-05-03    ET#410.41     2017-05-04
#> 39:          <NA>     2017-05-03   2017-05-03    ET#409.49     2017-05-04
#> 40:          <NA>     2017-05-03   2017-05-03    ET#409.18     2017-05-04
#> 41:          <NA>     2017-05-03   2017-05-03    ET#410.23     2017-05-04
#> 42:          <NA>     2017-05-03   2017-05-03    ET#409.80     2017-05-04
#> 43:          <NA>     2018-03-12   2018-03-12      FB584.2     2018-03-12
#> 44:          <NA>     2018-03-12   2018-03-12     FB584.48     2018-03-12
#> 45:          <NA>     2017-05-03   2017-05-03    ET#410.72     2017-05-04
#> 46:          <NA>     2017-05-03   2017-05-03    ET#410.65     2017-05-04
#> 47:          <NA>     2017-05-03   2017-05-03    ET#410.85     2017-05-04
#> 48:          <NA>     2017-05-03   2017-05-03    ET#410.10     2017-05-04
#> 49:          <NA>     2017-04-19   2017-04-19    ET#406.82     2017-04-20
#> 50:          <NA>     2017-05-03   2017-05-03     ET#410.3     2017-05-04
#> 51:          <NA>     2017-05-03   2017-05-03    ET#410.54     2017-05-04
#> 52:          <NA>     2017-05-03   2017-05-03    ET#410.34     2017-05-04
#> 53:          <NA>     2017-04-19   2017-04-19    ET#406.51     2017-04-20
#>     LinearPCRDate SonicationDate LigationDate 1stExpoPCRID 1stExpoPCRDate
#>            <Date>         <Date>       <Date>       <char>         <Date>
#>     2ndExpoID 2ndExpoDate FusionPrimerPCRID FusionPrimerPCRDate   PoolDate
#>        <char>      <Date>            <char>              <Date>     <Date>
#>  1:      <NA>        <NA>         ET#382.46          2016-11-03 2016-11-07
#>  2:      <NA>        <NA>         ET#381.40          2016-11-03 2016-11-07
#>  3:      <NA>        <NA>          ET#381.9          2016-11-03 2016-11-07
#>  4:      <NA>        <NA>         ET#381.71          2016-11-03 2016-11-07
#>  5:      <NA>        <NA>          ET#381.2          2016-11-03 2016-11-07
#>  6:      <NA>        <NA>         ET#382.28          2016-11-03 2016-11-07
#>  7:      <NA>        <NA>          ET#382.2          2016-11-03 2016-11-07
#>  8:      <NA>        <NA>         ET#382.50          2016-11-03 2016-11-07
#>  9:      <NA>        <NA>         ET#381.21          2016-11-03 2016-11-07
#> 10:      <NA>        <NA>         ET#407.28          2017-04-21 2017-04-27
#> 11:      <NA>        <NA>         ET#381.33          2016-11-03 2016-11-07
#> 12:      <NA>        <NA>         ET#382.15          2016-11-03 2016-11-07
#> 13:      <NA>        <NA>          ET#382.6          2016-11-03 2016-11-07
#> 14:      <NA>        <NA>         ET#381.52          2016-11-03 2016-11-07
#> 15:      <NA>        <NA>         ET#382.55          2016-11-03 2016-11-07
#> 16:      <NA>        <NA>         ET#381.25          2016-11-03 2016-11-07
#> 17:      <NA>        <NA>         ET#382.11          2016-11-03 2016-11-07
#> 18:      <NA>        <NA>         ET#382.33          2016-11-03 2016-11-07
#> 19:      <NA>        <NA>         ET#382.37          2016-11-03 2016-11-07
#> 20:      <NA>        <NA>         ET#382.59          2016-11-03 2016-11-07
#> 21:      <NA>        <NA>         ET#381.56          2016-11-03 2016-11-07
#> 22:      <NA>        <NA>         ET#381.87          2016-11-03 2016-11-07
#> 23:      <NA>        <NA>         ET#408.37          2017-04-21 2017-04-27
#> 24:      <NA>        <NA>         ET#382.24          2016-11-03 2016-11-07
#> 25:      <NA>        <NA>         ET#381.83          2016-11-03 2016-11-07
#> 26:      <NA>        <NA>         ET#381.64          2016-11-03 2016-11-07
#> 27:      <NA>        <NA>          FB585.26          2018-03-12 2018-03-13
#> 28:      <NA>        <NA>         ET#408.58          2017-04-21 2017-05-17
#> 29:      <NA>        <NA>         ET#408.27          2017-04-21 2017-05-17
#> 30:      <NA>        <NA>           FB585.4          2018-03-12 2018-03-13
#> 31:      <NA>        <NA>         ET#414.31          2017-05-16 2017-05-17
#> 32:      <NA>        <NA>         ET#408.89          2017-04-21 2017-05-17
#> 33:      <NA>        <NA>          FB585.24          2018-03-12 2018-03-13
#> 34:      <NA>        <NA>         ET#414.16          2017-05-16 2017-05-17
#> 35:      <NA>        <NA>          FB585.46          2018-03-12 2018-03-13
#> 36:      <NA>        <NA>         ET#411.42          2017-05-05 2017-05-17
#> 37:      <NA>        <NA>          ET#414.1          2017-05-16 2017-05-17
#> 38:      <NA>        <NA>         ET#412.41          2017-05-05 2017-05-17
#> 39:      <NA>        <NA>         ET#411.49          2017-05-05 2017-05-17
#> 40:      <NA>        <NA>         ET#411.18          2017-05-05 2017-05-17
#> 41:      <NA>        <NA>         ET#412.23          2017-05-05 2017-05-17
#> 42:      <NA>        <NA>         ET#411.80          2017-05-05 2017-05-17
#> 43:      <NA>        <NA>           FB585.2          2018-03-12 2018-03-13
#> 44:      <NA>        <NA>          FB585.48          2018-03-12 2018-03-13
#> 45:      <NA>        <NA>         ET#412.72          2017-05-05 2017-05-17
#> 46:      <NA>        <NA>         ET#412.65          2017-05-05 2017-05-17
#> 47:      <NA>        <NA>         ET#412.85          2017-05-05 2017-05-17
#> 48:      <NA>        <NA>         ET#412.10          2017-05-05 2017-05-17
#> 49:      <NA>        <NA>         ET#408.82          2017-04-21 2017-05-17
#> 50:      <NA>        <NA>          ET#412.3          2017-05-05 2017-05-17
#> 51:      <NA>        <NA>         ET#412.54          2017-05-05 2017-05-17
#> 52:      <NA>        <NA>         ET#412.34          2017-05-05 2017-05-17
#> 53:      <NA>        <NA>         ET#408.51          2017-04-21 2017-05-17
#>     2ndExpoID 2ndExpoDate FusionPrimerPCRID FusionPrimerPCRDate   PoolDate
#>        <char>      <Date>            <char>              <Date>     <Date>
#>     SequencingDate   VCN Genome SequencingRound Genotype TestGroup    MOI
#>             <Date> <num> <char>           <int>   <char>    <char> <char>
#>  1:     2016-11-15  0.30   hg19               1     <NA>      <NA>   <NA>
#>  2:     2016-11-15  0.27   hg19               1     <NA>      <NA>   <NA>
#>  3:     2016-11-15  0.27   hg19               1     <NA>      <NA>   <NA>
#>  4:     2016-11-15  0.27   hg19               1     <NA>      <NA>   <NA>
#>  5:     2016-11-15  0.24   hg19               1     <NA>      <NA>   <NA>
#>  6:     2016-11-15  0.42   hg19               1     <NA>      <NA>   <NA>
#>  7:     2016-11-15  0.30   hg19               1     <NA>      <NA>   <NA>
#>  8:     2016-11-15  0.42   hg19               1     <NA>      <NA>   <NA>
#>  9:     2016-11-15  0.23   hg19               1     <NA>      <NA>   <NA>
#> 10:     2017-06-23  0.19   hg19               1     <NA>      <NA>   <NA>
#> 11:     2016-11-15  0.24   hg19               1     <NA>      <NA>   <NA>
#> 12:     2016-11-15  0.35   hg19               1     <NA>      <NA>   <NA>
#> 13:     2016-11-15  0.42   hg19               1     <NA>      <NA>   <NA>
#> 14:     2016-11-15  0.23   hg19               1     <NA>      <NA>   <NA>
#> 15:     2016-11-15  0.26   hg19               1     <NA>      <NA>   <NA>
#> 16:     2016-11-15  0.26   hg19               1     <NA>      <NA>   <NA>
#> 17:     2016-11-15  0.26   hg19               1     <NA>      <NA>   <NA>
#> 18:     2016-11-15  0.26   hg19               1     <NA>      <NA>   <NA>
#> 19:     2016-11-15  0.35   hg19               1     <NA>      <NA>   <NA>
#> 20:     2016-11-15  0.35   hg19               1     <NA>      <NA>   <NA>
#> 21:     2016-11-15  0.26   hg19               1     <NA>      <NA>   <NA>
#> 22:     2016-11-15  0.26   hg19               1     <NA>      <NA>   <NA>
#> 23:     2017-06-23  0.18   hg19               1     <NA>      <NA>   <NA>
#> 24:     2016-11-15  0.30   hg19               1     <NA>      <NA>   <NA>
#> 25:     2016-11-15  0.23   hg19               1     <NA>      <NA>   <NA>
#> 26:     2016-11-15  0.24   hg19               1     <NA>      <NA>   <NA>
#> 27:     2018-03-15  2.52   hg19               1     <NA>      <NA>   <NA>
#> 28:     2017-06-23  1.01   hg19               1     <NA>      <NA>   <NA>
#> 29:     2017-06-23  1.01   hg19               1     <NA>      <NA>   <NA>
#> 30:     2018-03-15  2.52   hg19               1     <NA>      <NA>   <NA>
#> 31:     2017-06-23  2.30   hg19               1     <NA>      <NA>   <NA>
#> 32:     2017-06-23  1.01   hg19               1     <NA>      <NA>   <NA>
#> 33:     2018-03-15  2.07   hg19               1     <NA>      <NA>   <NA>
#> 34:     2017-06-23  2.30   hg19               1     <NA>      <NA>   <NA>
#> 35:     2018-03-15  2.07   hg19               1     <NA>      <NA>   <NA>
#> 36:     2017-06-23  1.07   hg19               1     <NA>      <NA>   <NA>
#> 37:     2017-06-23  2.30   hg19               1     <NA>      <NA>   <NA>
#> 38:     2017-06-23  2.05   hg19               1     <NA>      <NA>   <NA>
#> 39:     2017-06-23  1.24   hg19               1     <NA>      <NA>   <NA>
#> 40:     2017-06-23  1.24   hg19               1     <NA>      <NA>   <NA>
#> 41:     2017-06-23  1.43   hg19               1     <NA>      <NA>   <NA>
#> 42:     2017-06-23  1.24   hg19               1     <NA>      <NA>   <NA>
#> 43:     2018-03-15  2.07   hg19               1     <NA>      <NA>   <NA>
#> 44:     2018-03-15  2.52   hg19               1     <NA>      <NA>   <NA>
#> 45:     2017-06-23  2.05   hg19               1     <NA>      <NA>   <NA>
#> 46:     2017-06-23  1.09   hg19               1     <NA>      <NA>   <NA>
#> 47:     2017-06-23  1.43   hg19               1     <NA>      <NA>   <NA>
#> 48:     2017-06-23  2.05   hg19               1     <NA>      <NA>   <NA>
#> 49:     2017-06-23  0.77   hg19               1     <NA>      <NA>   <NA>
#> 50:     2017-06-23  1.09   hg19               1     <NA>      <NA>   <NA>
#> 51:     2017-06-23  1.43   hg19               1     <NA>      <NA>   <NA>
#> 52:     2017-06-23  1.09   hg19               1     <NA>      <NA>   <NA>
#> 53:     2017-06-23  0.77   hg19               1     <NA>      <NA>   <NA>
#>     SequencingDate   VCN Genome SequencingRound Genotype TestGroup    MOI
#>             <Date> <num> <char>           <int>   <char>    <char> <char>
#>     Engraftment Transduction  Notes AddedField1 AddedField2 AddedField3
#>           <num>        <num> <char>      <char>      <char>      <char>
#>  1:          NA           NA   <NA>        <NA>        <NA>        <NA>
#>  2:          NA           NA   <NA>        <NA>        <NA>        <NA>
#>  3:          NA           NA   <NA>        <NA>        <NA>        <NA>
#>  4:          NA           NA   <NA>        <NA>        <NA>        <NA>
#>  5:          NA           NA   <NA>        <NA>        <NA>        <NA>
#>  6:          NA           NA   <NA>        <NA>        <NA>        <NA>
#>  7:          NA           NA   <NA>        <NA>        <NA>        <NA>
#>  8:          NA           NA   <NA>        <NA>        <NA>        <NA>
#>  9:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 10:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 11:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 12:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 13:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 14:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 15:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 16:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 17:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 18:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 19:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 20:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 21:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 22:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 23:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 24:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 25:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 26:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 27:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 28:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 29:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 30:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 31:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 32:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 33:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 34:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 35:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 36:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 37:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 38:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 39:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 40:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 41:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 42:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 43:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 44:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 45:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 46:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 47:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 48:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 49:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 50:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 51:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 52:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 53:          NA           NA   <NA>        <NA>        <NA>        <NA>
#>     Engraftment Transduction  Notes AddedField1 AddedField2 AddedField3
#>           <num>        <num> <char>      <char>      <char>      <char>
#>     AddedField4 concatenatePoolIDSeqRun AddedField6_RelativeBloodPercentage
#>          <char>                  <char>                              <char>
#>  1:        <NA>                POOL01-1                                <NA>
#>  2:        <NA>                POOL01-1                                <NA>
#>  3:        <NA>                POOL01-1                                <NA>
#>  4:        <NA>                POOL01-1                                <NA>
#>  5:        <NA>                POOL01-1                                <NA>
#>  6:        <NA>                POOL01-1                                <NA>
#>  7:        <NA>                POOL01-1                                <NA>
#>  8:        <NA>                POOL01-1                                <NA>
#>  9:        <NA>                POOL01-1                                <NA>
#> 10:        <NA>                POOL02-1                                <NA>
#> 11:        <NA>                POOL01-1                                <NA>
#> 12:        <NA>                POOL01-1                                <NA>
#> 13:        <NA>                POOL01-1                                <NA>
#> 14:        <NA>                POOL01-1                                <NA>
#> 15:        <NA>                POOL01-1                                <NA>
#> 16:        <NA>                POOL01-1                                <NA>
#> 17:        <NA>                POOL01-1                                <NA>
#> 18:        <NA>                POOL01-1                                <NA>
#> 19:        <NA>                POOL01-1                                <NA>
#> 20:        <NA>                POOL01-1                                <NA>
#> 21:        <NA>                POOL01-1                                <NA>
#> 22:        <NA>                POOL01-1                                <NA>
#> 23:        <NA>                POOL02-1                                <NA>
#> 24:        <NA>                POOL01-1                                <NA>
#> 25:        <NA>                POOL01-1                                <NA>
#> 26:        <NA>                POOL01-1                                <NA>
#> 27:        <NA>                POOL03-1                                <NA>
#> 28:        <NA>                POOL04-1                                <NA>
#> 29:        <NA>                POOL04-1                                <NA>
#> 30:        <NA>                POOL03-1                                <NA>
#> 31:        <NA>                POOL04-1                                <NA>
#> 32:        <NA>                POOL04-1                                <NA>
#> 33:        <NA>                POOL03-1                                <NA>
#> 34:        <NA>                POOL04-1                                <NA>
#> 35:        <NA>                POOL03-1                                <NA>
#> 36:        <NA>                POOL04-1                                <NA>
#> 37:        <NA>                POOL04-1                                <NA>
#> 38:        <NA>                POOL04-1                                <NA>
#> 39:        <NA>                POOL04-1                                <NA>
#> 40:        <NA>                POOL04-1                                <NA>
#> 41:        <NA>                POOL04-1                                <NA>
#> 42:        <NA>                POOL04-1                                <NA>
#> 43:        <NA>                POOL03-1                                <NA>
#> 44:        <NA>                POOL03-1                                <NA>
#> 45:        <NA>                POOL04-1                                <NA>
#> 46:        <NA>                POOL04-1                                <NA>
#> 47:        <NA>                POOL04-1                                <NA>
#> 48:        <NA>                POOL04-1                                <NA>
#> 49:        <NA>                POOL04-1                                <NA>
#> 50:        <NA>                POOL04-1                                <NA>
#> 51:        <NA>                POOL04-1                                <NA>
#> 52:        <NA>                POOL04-1                                <NA>
#> 53:        <NA>                POOL04-1                                <NA>
#>     AddedField4 concatenatePoolIDSeqRun AddedField6_RelativeBloodPercentage
#>          <char>                  <char>                              <char>
#>     AddedField7_PurityTestFeasibility AddedField8_FacsSeparationPurity
#>                                 <num>                            <num>
#>  1:                                NA                               NA
#>  2:                                NA                               NA
#>  3:                                NA                               NA
#>  4:                                NA                               NA
#>  5:                                NA                               NA
#>  6:                                NA                               NA
#>  7:                                NA                               NA
#>  8:                                NA                               NA
#>  9:                                NA                               NA
#> 10:                                NA                               NA
#> 11:                                NA                               NA
#> 12:                                NA                               NA
#> 13:                                NA                               NA
#> 14:                                NA                               NA
#> 15:                                NA                               NA
#> 16:                                NA                               NA
#> 17:                                NA                               NA
#> 18:                                NA                               NA
#> 19:                                NA                               NA
#> 20:                                NA                               NA
#> 21:                                NA                               NA
#> 22:                                NA                               NA
#> 23:                                NA                               NA
#> 24:                                NA                               NA
#> 25:                                NA                               NA
#> 26:                                NA                               NA
#> 27:                                NA                               NA
#> 28:                                NA                               NA
#> 29:                                NA                               NA
#> 30:                                NA                               NA
#> 31:                                NA                               NA
#> 32:                                NA                               NA
#> 33:                                NA                               NA
#> 34:                                NA                               NA
#> 35:                                NA                               NA
#> 36:                                NA                               NA
#> 37:                                NA                               NA
#> 38:                                NA                               NA
#> 39:                                NA                               NA
#> 40:                                NA                               NA
#> 41:                                NA                               NA
#> 42:                                NA                               NA
#> 43:                                NA                               NA
#> 44:                                NA                               NA
#> 45:                                NA                               NA
#> 46:                                NA                               NA
#> 47:                                NA                               NA
#> 48:                                NA                               NA
#> 49:                                NA                               NA
#> 50:                                NA                               NA
#> 51:                                NA                               NA
#> 52:                                NA                               NA
#> 53:                                NA                               NA
#>     AddedField7_PurityTestFeasibility AddedField8_FacsSeparationPurity
#>                                 <num>                            <num>
#>          Kapa ulForPool
#>         <num>     <num>
#>  1:        NA        NA
#>  2:        NA        NA
#>  3:        NA        NA
#>  4:        NA        NA
#>  5:        NA        NA
#>  6:        NA        NA
#>  7:        NA        NA
#>  8:        NA        NA
#>  9:        NA        NA
#> 10:        NA        NA
#> 11:        NA        NA
#> 12:        NA        NA
#> 13:        NA        NA
#> 14:        NA        NA
#> 15:        NA        NA
#> 16:        NA        NA
#> 17:        NA        NA
#> 18:        NA        NA
#> 19:        NA        NA
#> 20:        NA        NA
#> 21:        NA        NA
#> 22:        NA        NA
#> 23:        NA        NA
#> 24:        NA        NA
#> 25:        NA        NA
#> 26:        NA        NA
#> 27: 57.128501  1.000000
#> 28:        NA        NA
#> 29:        NA        NA
#> 30: 41.796937  1.000000
#> 31:        NA        NA
#> 32:        NA        NA
#> 33:  9.813033  1.528579
#> 34:        NA        NA
#> 35: 52.331328  1.000000
#> 36:        NA        NA
#> 37:        NA        NA
#> 38:        NA        NA
#> 39:        NA        NA
#> 40:        NA        NA
#> 41:        NA        NA
#> 42:        NA        NA
#> 43:  3.820188  3.926508
#> 44: 16.312071  1.000000
#> 45:        NA        NA
#> 46:        NA        NA
#> 47:        NA        NA
#> 48:        NA        NA
#> 49:        NA        NA
#> 50:        NA        NA
#> 51:        NA        NA
#> 52:        NA        NA
#> 53:        NA        NA
#>          Kapa ulForPool
#>         <num>     <num>
#>                                                  CompleteAmplificationID
#>                                                                   <char>
#>  1: PJ01_POOL01_LTR75LC38_PT001_PT001-103_lenti_GLOBE_PB_1_SLiM_0060_MNC
#>  2:  PJ01_POOL01_LTR53LC32_PT001_PT001-81_lenti_GLOBE_BM_1_SLiM_0180_MNC
#>  3:  PJ01_POOL01_LTR83LC66_PT001_PT001-81_lenti_GLOBE_BM_1_SLiM_0180_MNC
#>  4:  PJ01_POOL01_LTR27LC94_PT001_PT001-81_lenti_GLOBE_BM_1_SLiM_0180_MNC
#>  5:  PJ01_POOL01_LTR69LC52_PT001_PT001-74_lenti_GLOBE_PB_1_SLiM_0180_MNC
#>  6:  PJ01_POOL01_LTR37LC2_PT001_PT001-107_lenti_GLOBE_BM_1_SLiM_0060_MNC
#>  7: PJ01_POOL01_LTR77LC46_PT001_PT001-103_lenti_GLOBE_PB_1_SLiM_0060_MNC
#>  8: PJ01_POOL01_LTR83LC46_PT001_PT001-107_lenti_GLOBE_BM_1_SLiM_0060_MNC
#>  9:   PJ01_POOL01_LTR9LC90_PT001_PT001-93_lenti_GLOBE_PB_1_SLiM_0030_MNC
#> 10: PJ01_POOL02_LTR65LC56_PT001_PT001-149_lenti_GLOBE_PB_1_SLiM_0360_MNC
#> 11:  PJ01_POOL01_LTR37LC18_PT001_PT001-74_lenti_GLOBE_PB_1_SLiM_0180_MNC
#> 12:  PJ01_POOL01_LTR7LC72_PT001_PT001-116_lenti_GLOBE_BM_1_SLiM_0090_MNC
#> 13: PJ01_POOL01_LTR85LC54_PT001_PT001-107_lenti_GLOBE_BM_1_SLiM_0060_MNC
#> 14:  PJ01_POOL01_LTR77LC56_PT001_PT001-93_lenti_GLOBE_PB_1_SLiM_0030_MNC
#> 15: PJ01_POOL01_LTR93LC56_PT001_PT001-112_lenti_GLOBE_PB_1_SLiM_0090_MNC
#> 16:   PJ01_POOL01_LTR19LC2_PT001_PT001-97_lenti_GLOBE_BM_1_SLiM_0030_MNC
#> 17: PJ01_POOL01_LTR95LC64_PT001_PT001-112_lenti_GLOBE_PB_1_SLiM_0090_MNC
#> 18: PJ01_POOL01_LTR49LC12_PT001_PT001-112_lenti_GLOBE_PB_1_SLiM_0090_MNC
#> 19: PJ01_POOL01_LTR57LC20_PT001_PT001-116_lenti_GLOBE_BM_1_SLiM_0090_MNC
#> 20:  PJ01_POOL01_LTR5LC64_PT001_PT001-116_lenti_GLOBE_BM_1_SLiM_0090_MNC
#> 21:  PJ01_POOL01_LTR85LC64_PT001_PT001-97_lenti_GLOBE_BM_1_SLiM_0030_MNC
#> 22:  PJ01_POOL01_LTR61LC30_PT001_PT001-97_lenti_GLOBE_BM_1_SLiM_0030_MNC
#> 23: PJ01_POOL02_LTR87LC74_PT001_PT001-156_lenti_GLOBE_BM_1_SLiM_0360_MNC
#> 24: PJ01_POOL01_LTR29LC90_PT001_PT001-103_lenti_GLOBE_PB_1_SLiM_0060_MNC
#> 25:  PJ01_POOL01_LTR53LC22_PT001_PT001-93_lenti_GLOBE_PB_1_SLiM_0030_MNC
#> 26:   PJ01_POOL01_LTR9LC80_PT001_PT001-74_lenti_GLOBE_PB_1_SLiM_0180_MNC
#> 27: PJ01_POOL03_LTR73LC86_PT002_PT002-466_lenti_GLOBE_BM_1_SLiM_0360_MNC
#> 28: PJ01_POOL04_LTR41LC20_PT002_PT002-177_lenti_GLOBE_BM_1_SLiM_0030_MNC
#> 29: PJ01_POOL04_LTR65LC54_PT002_PT002-177_lenti_GLOBE_BM_1_SLiM_0030_MNC
#> 30: PJ01_POOL03_LTR51LC86_PT002_PT002-466_lenti_GLOBE_BM_1_SLiM_0360_MNC
#> 31: PJ01_POOL04_LTR85LC62_PT002_PT002-238_lenti_GLOBE_BM_1_SLiM_0180_MNC
#> 32: PJ01_POOL04_LTR11LC82_PT002_PT002-177_lenti_GLOBE_BM_1_SLiM_0030_MNC
#> 33: PJ01_POOL03_LTR71LC90_PT002_PT002-464_lenti_GLOBE_PB_1_SLiM_0360_MNC
#> 34: PJ01_POOL04_LTR51LC32_PT002_PT002-238_lenti_GLOBE_BM_1_SLiM_0180_MNC
#> 35: PJ01_POOL03_LTR93LC90_PT002_PT002-464_lenti_GLOBE_PB_1_SLiM_0360_MNC
#> 36:  PJ01_POOL04_LTR5LC84_PT002_PT002-190_lenti_GLOBE_PB_1_SLiM_0060_MNC
#> 37:  PJ01_POOL04_LTR13LC2_PT002_PT002-238_lenti_GLOBE_BM_1_SLiM_0180_MNC
#> 38:  PJ01_POOL04_LTR5LC82_PT002_PT002-218_lenti_GLOBE_BM_1_SLiM_0090_MNC
#> 39:  PJ01_POOL04_LTR25LC2_PT002_PT002-197_lenti_GLOBE_BM_1_SLiM_0060_MNC
#> 40: PJ01_POOL04_LTR49LC36_PT002_PT002-197_lenti_GLOBE_BM_1_SLiM_0060_MNC
#> 41: PJ01_POOL04_LTR63LC46_PT002_PT002-231_lenti_GLOBE_PB_1_SLiM_0180_MNC
#> 42: PJ01_POOL04_LTR93LC64_PT002_PT002-197_lenti_GLOBE_BM_1_SLiM_0060_MNC
#> 43: PJ01_POOL03_LTR49LC90_PT002_PT002-464_lenti_GLOBE_PB_1_SLiM_0360_MNC
#> 44: PJ01_POOL03_LTR95LC86_PT002_PT002-466_lenti_GLOBE_BM_1_SLiM_0360_MNC
#> 45: PJ01_POOL04_LTR79LC48_PT002_PT002-218_lenti_GLOBE_BM_1_SLiM_0090_MNC
#> 46: PJ01_POOL04_LTR63LC34_PT002_PT002-211_lenti_GLOBE_PB_1_SLiM_0090_MNC
#> 47:  PJ01_POOL04_LTR9LC74_PT002_PT002-231_lenti_GLOBE_PB_1_SLiM_0180_MNC
#> 48: PJ01_POOL04_LTR33LC20_PT002_PT002-218_lenti_GLOBE_BM_1_SLiM_0090_MNC
#> 49: PJ01_POOL04_LTR93LC68_PT002_PT002-170_lenti_GLOBE_PB_1_SLiM_0030_MNC
#> 50:  PJ01_POOL04_LTR13LC6_PT002_PT002-211_lenti_GLOBE_PB_1_SLiM_0090_MNC
#> 51: PJ01_POOL04_LTR37LC12_PT002_PT002-231_lenti_GLOBE_PB_1_SLiM_0180_MNC
#> 52: PJ01_POOL04_LTR87LC68_PT002_PT002-211_lenti_GLOBE_PB_1_SLiM_0090_MNC
#> 53:  PJ01_POOL04_LTR25LC6_PT002_PT002-170_lenti_GLOBE_PB_1_SLiM_0030_MNC
#>                                                  CompleteAmplificationID
#>                                                                   <char>
#>                   UniqueID StudyTestID StudyTestGroup MouseID Tigroup Tisource
#>                     <char>      <char>          <int>   <int>  <char>   <char>
#>  1: ID00000000000000007433        <NA>             NA      NA    <NA>     <NA>
#>  2: ID00000000000000007340        <NA>             NA      NA    <NA>     <NA>
#>  3: ID00000000000000007310        <NA>             NA      NA    <NA>     <NA>
#>  4: ID00000000000000007370        <NA>             NA      NA    <NA>     <NA>
#>  5: ID00000000000000007303        <NA>             NA      NA    <NA>     <NA>
#>  6: ID00000000000000007417        <NA>             NA      NA    <NA>     <NA>
#>  7: ID00000000000000007393        <NA>             NA      NA    <NA>     <NA>
#>  8: ID00000000000000007437        <NA>             NA      NA    <NA>     <NA>
#>  9: ID00000000000000007322        <NA>             NA      NA    <NA>     <NA>
#> 10: ID00000000000000009697        <NA>             NA      NA    <NA>     <NA>
#> 11: ID00000000000000007333        <NA>             NA      NA    <NA>     <NA>
#> 12: ID00000000000000007406        <NA>             NA      NA    <NA>     <NA>
#> 13: ID00000000000000007397        <NA>             NA      NA    <NA>     <NA>
#> 14: ID00000000000000007352        <NA>             NA      NA    <NA>     <NA>
#> 15: ID00000000000000007442        <NA>             NA      NA    <NA>     <NA>
#> 16: ID00000000000000007326        <NA>             NA      NA    <NA>     <NA>
#> 17: ID00000000000000007402        <NA>             NA      NA    <NA>     <NA>
#> 18: ID00000000000000007422        <NA>             NA      NA    <NA>     <NA>
#> 19: ID00000000000000007426        <NA>             NA      NA    <NA>     <NA>
#> 20: ID00000000000000007446        <NA>             NA      NA    <NA>     <NA>
#> 21: ID00000000000000007356        <NA>             NA      NA    <NA>     <NA>
#> 22: ID00000000000000007386        <NA>             NA      NA    <NA>     <NA>
#> 23: ID00000000000000009748        <NA>             NA      NA    <NA>     <NA>
#> 24: ID00000000000000007413        <NA>             NA      NA    <NA>     <NA>
#> 25: ID00000000000000007382        <NA>             NA      NA    <NA>     <NA>
#> 26: ID00000000000000007363        <NA>             NA      NA    <NA>     <NA>
#> 27: ID00000000000000013890        <NA>             NA      NA    <NA>     <NA>
#> 28: ID00000000000000009785        <NA>             NA      NA    <NA>     <NA>
#> 29: ID00000000000000009775        <NA>             NA      NA    <NA>     <NA>
#> 30: ID00000000000000013868        <NA>             NA      NA    <NA>     <NA>
#> 31: ID00000000000000009953        <NA>             NA      NA    <NA>     <NA>
#> 32: ID00000000000000009793        <NA>             NA      NA    <NA>     <NA>
#> 33: ID00000000000000013888        <NA>             NA      NA    <NA>     <NA>
#> 34: ID00000000000000009942        <NA>             NA      NA    <NA>     <NA>
#> 35: ID00000000000000013910        <NA>             NA      NA    <NA>     <NA>
#> 36: ID00000000000000009826        <NA>             NA      NA    <NA>     <NA>
#> 37: ID00000000000000009934        <NA>             NA      NA    <NA>     <NA>
#> 38: ID00000000000000009890        <NA>             NA      NA    <NA>     <NA>
#> 39: ID00000000000000009833        <NA>             NA      NA    <NA>     <NA>
#> 40: ID00000000000000009807        <NA>             NA      NA    <NA>     <NA>
#> 41: ID00000000000000009875        <NA>             NA      NA    <NA>     <NA>
#> 42: ID00000000000000009850        <NA>             NA      NA    <NA>     <NA>
#> 43: ID00000000000000013866        <NA>             NA      NA    <NA>     <NA>
#> 44: ID00000000000000013912        <NA>             NA      NA    <NA>     <NA>
#> 45: ID00000000000000009917        <NA>             NA      NA    <NA>     <NA>
#> 46: ID00000000000000009910        <NA>             NA      NA    <NA>     <NA>
#> 47: ID00000000000000009928        <NA>             NA      NA    <NA>     <NA>
#> 48: ID00000000000000009866        <NA>             NA      NA    <NA>     <NA>
#> 49: ID00000000000000009789        <NA>             NA      NA    <NA>     <NA>
#> 50: ID00000000000000009859        <NA>             NA      NA    <NA>     <NA>
#> 51: ID00000000000000009903        <NA>             NA      NA    <NA>     <NA>
#> 52: ID00000000000000009884        <NA>             NA      NA    <NA>     <NA>
#> 53: ID00000000000000009780        <NA>             NA      NA    <NA>     <NA>
#>                   UniqueID StudyTestID StudyTestGroup MouseID Tigroup Tisource
#>                     <char>      <char>          <int>   <int>  <char>   <char>
#>     PathToFolderProjectID SamplesNameCheck TimepointDays TimepointMonths
#>                    <char>           <char>        <char>          <char>
#>  1:                 /PJ01             <NA>          0060              02
#>  2:                 /PJ01             <NA>          0180              06
#>  3:                 /PJ01             <NA>          0180              06
#>  4:                 /PJ01             <NA>          0180              06
#>  5:                 /PJ01             <NA>          0180              06
#>  6:                 /PJ01             <NA>          0060              02
#>  7:                 /PJ01             <NA>          0060              02
#>  8:                 /PJ01             <NA>          0060              02
#>  9:                 /PJ01             <NA>          0030              01
#> 10:                 /PJ01             <NA>          0360              12
#> 11:                 /PJ01             <NA>          0180              06
#> 12:                 /PJ01             <NA>          0090              03
#> 13:                 /PJ01             <NA>          0060              02
#> 14:                 /PJ01             <NA>          0030              01
#> 15:                 /PJ01             <NA>          0090              03
#> 16:                 /PJ01             <NA>          0030              01
#> 17:                 /PJ01             <NA>          0090              03
#> 18:                 /PJ01             <NA>          0090              03
#> 19:                 /PJ01             <NA>          0090              03
#> 20:                 /PJ01             <NA>          0090              03
#> 21:                 /PJ01             <NA>          0030              01
#> 22:                 /PJ01             <NA>          0030              01
#> 23:                 /PJ01             <NA>          0360              12
#> 24:                 /PJ01             <NA>          0060              02
#> 25:                 /PJ01             <NA>          0030              01
#> 26:                 /PJ01             <NA>          0180              06
#> 27:                 /PJ01             <NA>          0360              12
#> 28:                 /PJ01             <NA>          0030              01
#> 29:                 /PJ01             <NA>          0030              01
#> 30:                 /PJ01             <NA>          0360              12
#> 31:                 /PJ01             <NA>          0180              06
#> 32:                 /PJ01             <NA>          0030              01
#> 33:                 /PJ01             <NA>          0360              12
#> 34:                 /PJ01             <NA>          0180              06
#> 35:                 /PJ01             <NA>          0360              12
#> 36:                 /PJ01             <NA>          0060              02
#> 37:                 /PJ01             <NA>          0180              06
#> 38:                 /PJ01             <NA>          0090              03
#> 39:                 /PJ01             <NA>          0060              02
#> 40:                 /PJ01             <NA>          0060              02
#> 41:                 /PJ01             <NA>          0180              06
#> 42:                 /PJ01             <NA>          0060              02
#> 43:                 /PJ01             <NA>          0360              12
#> 44:                 /PJ01             <NA>          0360              12
#> 45:                 /PJ01             <NA>          0090              03
#> 46:                 /PJ01             <NA>          0090              03
#> 47:                 /PJ01             <NA>          0180              06
#> 48:                 /PJ01             <NA>          0090              03
#> 49:                 /PJ01             <NA>          0030              01
#> 50:                 /PJ01             <NA>          0090              03
#> 51:                 /PJ01             <NA>          0180              06
#> 52:                 /PJ01             <NA>          0090              03
#> 53:                 /PJ01             <NA>          0030              01
#>     PathToFolderProjectID SamplesNameCheck TimepointDays TimepointMonths
#>                    <char>           <char>        <char>          <char>
#>     TimepointYears ng DNA corrected      RUN_NAME PHIX_MAPPING
#>             <char>            <num>        <char>        <int>
#>  1:             01            23.18 PJ01|POOL01-1     43586699
#>  2:             01           181.44 PJ01|POOL01-1     43586699
#>  3:             01           181.44 PJ01|POOL01-1     43586699
#>  4:             01           181.44 PJ01|POOL01-1     43586699
#>  5:             01            23.06 PJ01|POOL01-1     43586699
#>  6:             01           171.36 PJ01|POOL01-1     43586699
#>  7:             01            23.18 PJ01|POOL01-1     43586699
#>  8:             01           171.36 PJ01|POOL01-1     43586699
#>  9:             01            23.85 PJ01|POOL01-1     43586699
#> 10:             01            45.15 PJ01|POOL02-1     20409455
#> 11:             01            23.06 PJ01|POOL01-1     43586699
#> 12:             01            89.21 PJ01|POOL01-1     43586699
#> 13:             01           171.36 PJ01|POOL01-1     43586699
#> 14:             01            23.85 PJ01|POOL01-1     43586699
#> 15:             01            26.57 PJ01|POOL01-1     43586699
#> 16:             01           300.03 PJ01|POOL01-1     43586699
#> 17:             01            26.57 PJ01|POOL01-1     43586699
#> 18:             01            26.57 PJ01|POOL01-1     43586699
#> 19:             01            89.21 PJ01|POOL01-1     43586699
#> 20:             01            89.21 PJ01|POOL01-1     43586699
#> 21:             01           300.03 PJ01|POOL01-1     43586699
#> 22:             01           300.03 PJ01|POOL01-1     43586699
#> 23:             01            42.00 PJ01|POOL02-1     20409455
#> 24:             01            23.18 PJ01|POOL01-1     43586699
#> 25:             01            23.85 PJ01|POOL01-1     43586699
#> 26:             01            23.06 PJ01|POOL01-1     43586699
#> 27:             01           299.32 PJ01|POOL03-1     51183662
#> 28:             01           181.50 PJ01|POOL04-1     18979629
#> 29:             01           181.50 PJ01|POOL04-1     18979629
#> 30:             01           299.32 PJ01|POOL03-1     51183662
#> 31:             01           196.50 PJ01|POOL04-1     18979629
#> 32:             01           181.50 PJ01|POOL04-1     18979629
#> 33:             01           300.48 PJ01|POOL03-1     51183662
#> 34:             01           196.50 PJ01|POOL04-1     18979629
#> 35:             01           300.48 PJ01|POOL03-1     51183662
#> 36:             01            77.40 PJ01|POOL04-1     18979629
#> 37:             01           196.50 PJ01|POOL04-1     18979629
#> 38:             01           135.15 PJ01|POOL04-1     18979629
#> 39:             01           172.50 PJ01|POOL04-1     18979629
#> 40:             01           172.50 PJ01|POOL04-1     18979629
#> 41:             01           122.55 PJ01|POOL04-1     18979629
#> 42:             01           172.50 PJ01|POOL04-1     18979629
#> 43:             01           300.48 PJ01|POOL03-1     51183662
#> 44:             01           299.32 PJ01|POOL03-1     51183662
#> 45:             01           135.15 PJ01|POOL04-1     18979629
#> 46:             01           123.00 PJ01|POOL04-1     18979629
#> 47:             01           122.55 PJ01|POOL04-1     18979629
#> 48:             01           135.15 PJ01|POOL04-1     18979629
#> 49:             01           106.35 PJ01|POOL04-1     18979629
#> 50:             01           123.00 PJ01|POOL04-1     18979629
#> 51:             01           122.55 PJ01|POOL04-1     18979629
#> 52:             01           123.00 PJ01|POOL04-1     18979629
#> 53:             01           106.35 PJ01|POOL04-1     18979629
#>     TimepointYears ng DNA corrected      RUN_NAME PHIX_MAPPING
#>             <char>            <num>        <char>        <int>
#>     PLASMID_MAPPED_BYPOOL BARCODE_MUX LTR_IDENTIFIED TRIMMING_FINAL_LTRLC
#>                     <int>       <int>          <int>                <int>
#>  1:               2256176      645026         645026               630965
#>  2:               2256176      652208         652177               649044
#>  3:               2256176      451519         451512               449669
#>  4:               2256176      426500         426499               425666
#>  5:               2256176       18300          18300                18290
#>  6:               2256176      729327         729327               727219
#>  7:               2256176       14999          14999                14996
#>  8:               2256176      748869         748865               746189
#>  9:               2256176      445754         445754               445717
#> 10:               1560035     1178016        1178016              1170116
#> 11:               2256176       36760          36760                36759
#> 12:               2256176      738518         738518               738512
#> 13:               2256176      971374         971339               968244
#> 14:               2256176      430438         430438               430423
#> 15:               2256176       33023          33023                33022
#> 16:               2256176      513783         513365               511701
#> 17:               2256176      635010         635010               634929
#> 18:               2256176       84541          84541                84540
#> 19:               2256176      633931         633931               631334
#> 20:               2256176      408400         408400               408393
#> 21:               2256176      553160         553160               551962
#> 22:               2256176      454734         454734               453088
#> 23:               1560035      800640         800640               800594
#> 24:               2256176       44888          44888                44876
#> 25:               2256176        7477           7477                 7474
#> 26:               2256176       71185          71185                71181
#> 27:                392602     6120512        6120512              6104882
#> 28:               1074651      456215         456210               452788
#> 29:               1074651      484852         484843               480224
#> 30:                392602     3536618        3536575              3532982
#> 31:               1074651      641193         640603               636224
#> 32:               1074651      660233         660117               654590
#> 33:                392602      963093         963017               961733
#> 34:               1074651      790616         790360               783702
#> 35:                392602     5308967        5308966              5280083
#> 36:               1074651      342636         342485               338616
#> 37:               1074651      490058         489719               487186
#> 38:               1074651      524854         524774               517187
#> 39:               1074651      308320         308114               306191
#> 40:               1074651      508548         508425               501959
#> 41:               1074651      995886         995864               980033
#> 42:               1074651      342626         342539               340481
#> 43:                392602      594361         594361               593948
#> 44:                392602     1827545        1827545              1825772
#> 45:               1074651      347728         347697               343609
#> 46:               1074651      346363         346362               342932
#> 47:               1074651      354092         353983               350767
#> 48:               1074651      405322         405276               402861
#> 49:               1074651      631307         631307               622869
#> 50:               1074651      361156         361138               358344
#> 51:               1074651      339082         338995               335745
#> 52:               1074651      331492         331477               328046
#> 53:               1074651      392967         392765               390451
#>     PLASMID_MAPPED_BYPOOL BARCODE_MUX LTR_IDENTIFIED TRIMMING_FINAL_LTRLC
#>                     <int>       <int>          <int>                <int>
#>     LV_MAPPED BWA_MAPPED_OVERALL ISS_MAPPED_OVERALL RAW_READS QUALITY_PASSED
#>         <int>              <int>              <int>     <int>          <int>
#>  1:    211757             402477             219452        NA             NA
#>  2:    303300             322086             222646        NA             NA
#>  3:    204810             227275             149385        NA             NA
#>  4:    185752             223915             143283        NA             NA
#>  5:      6962              10487               5907        NA             NA
#>  6:    318653             369117             235640        NA             NA
#>  7:      6006               8591               5829        NA             NA
#>  8:    339742             376479             247960        NA             NA
#>  9:    217780             227190             192502        NA             NA
#> 10:    567207             533620             310663        NA             NA
#> 11:     31750               4654               4028        NA             NA
#> 12:    345886             372292             251017        NA             NA
#> 13:    410946             512137             320550        NA             NA
#> 14:    195339             215676             176573        NA             NA
#> 15:     19255              12545               8953        NA             NA
#> 16:    213205             274453             178588        NA             NA
#> 17:    273550             312464             173183        NA             NA
#> 18:     55966              26795              24950        NA             NA
#> 19:    273665             331492             207297        NA             NA
#> 20:    193851             205558             142060        NA             NA
#> 21:    249875             283399             189193        NA             NA
#> 22:    201674             234982             160477        NA             NA
#> 23:    396738             370781             234541        NA             NA
#> 24:     18502              23817              16606        NA             NA
#> 25:      3182               3944               2824        NA             NA
#> 26:     37973              30617              19250        NA             NA
#> 27:   2926769            2900368            1966725        NA             NA
#> 28:    200981             210935             119095        NA             NA
#> 29:    209373             214620             114170        NA             NA
#> 30:   1663771            1787575            1334185        NA             NA
#> 31:    297970             300563             180099        NA             NA
#> 32:    278720             310685             171359        NA             NA
#> 33:    449611             496042             374392        NA             NA
#> 34:    348153             378898             222706        NA             NA
#> 35:   2528864            2455083            1633331        NA             NA
#> 36:    156154             162504              94468        NA             NA
#> 37:    218398             241840             150078        NA             NA
#> 38:    223785             235568             125591        NA             NA
#> 39:    139167             142056              83601        NA             NA
#> 40:    215787             234853             128363        NA             NA
#> 41:    448587             368752             173206        NA             NA
#> 42:    150283             161629              93421        NA             NA
#> 43:    291931             293238             221993        NA             NA
#> 44:    885302             903768             673344        NA             NA
#> 45:    155701             155486              84035        NA             NA
#> 46:    149781             167400              94120        NA             NA
#> 47:    153595             163215              92043        NA             NA
#> 48:    177016             187252             107526        NA             NA
#> 49:    290076             271014             133798        NA             NA
#> 50:    156230             165088              90776        NA             NA
#> 51:    149265             158919              88033        NA             NA
#> 52:    146490             147277              75439        NA             NA
#> 53:    184167             172582              96911        NA             NA
#>     LV_MAPPED BWA_MAPPED_OVERALL ISS_MAPPED_OVERALL RAW_READS QUALITY_PASSED
#>         <int>              <int>              <int>     <int>          <int>
#>     ISS_MAPPED_PP seqCount_sum seqCount_count seqCount_shannon seqCount_simpson
#>             <int>        <num>          <int>            <num>            <num>
#>  1:            NA          357             15        1.4928134        0.6693815
#>  2:            NA       129867             72        3.2386944        0.9183413
#>  3:            NA        77474             76        3.1390767        0.9088955
#>  4:            NA        80610             67        3.0212128        0.9041340
#>  5:            NA          147             26        2.4915716        0.8738951
#>  6:            NA        50859             82        3.7551750        0.9623687
#>  7:            NA          111             41        3.1703844        0.9263858
#>  8:            NA        53887             90        3.8384970        0.9721518
#>  9:            NA       121304             11        1.7647042        0.8132592
#> 10:            NA        16723             43        3.4154154        0.9626384
#> 11:            NA           44             10        1.9999913        0.8326446
#> 12:            NA        75905             31        2.7312148        0.9128491
#> 13:            NA        57684             77        3.8446041        0.9736556
#> 14:            NA       120472             19        2.3275229        0.8971868
#> 15:            NA           81             13        2.3180241        0.8876696
#> 16:            NA         3029             22        1.8161366        0.8042655
#> 17:            NA       116284             29        1.8520688        0.8373618
#> 18:            NA           33             10        1.9298441        0.8191001
#> 19:            NA        57443             20        2.4093193        0.8977445
#> 20:            NA        41488             39        2.8907905        0.9382894
#> 21:            NA         3697             32        1.9732310        0.8417839
#> 22:            NA          559             10        0.8013509        0.4811877
#> 23:            NA        42122             78        4.0692422        0.9797955
#> 24:            NA           30              8        1.8858896        0.8288889
#> 25:            NA            1              1        0.0000000        0.0000000
#> 26:            NA            1              1        0.0000000        0.0000000
#> 27:            NA        51114             48        3.5182831        0.9645795
#> 28:            NA         7516             65        3.9064574        0.9751075
#> 29:            NA         5978             74        3.7693672        0.9683265
#> 30:            NA        47816             46        3.3455043        0.9558975
#> 31:            NA         5014             48        3.4524970        0.9589726
#> 32:            NA         9150             87        4.0217856        0.9787792
#> 33:            NA         4287             12        2.1443631        0.8693443
#> 34:            NA         5238             52        3.5545359        0.9636314
#> 35:            NA        45416             35        3.2141718        0.9518653
#> 36:            NA         1468             18        1.9669258        0.8235138
#> 37:            NA         3280             52        3.6004539        0.9671765
#> 38:            NA         1029             36        2.2300233        0.8582969
#> 39:            NA         1851             20        2.5605560        0.9056275
#> 40:            NA         2089             19        2.4844667        0.9017753
#> 41:            NA           69             24        2.7403346        0.9103130
#> 42:            NA         1271             25        2.5482273        0.8953097
#> 43:            NA         1291              9        1.8695596        0.8200513
#> 44:            NA        16360             26        2.8263876        0.9281503
#> 45:            NA          565             10        1.9839406        0.8496076
#> 46:            NA         1538             14        1.8294052        0.8298129
#> 47:            NA         1958              8        1.7666363        0.8149891
#> 48:            NA          661              9        1.9948568        0.8502498
#> 49:            NA         1388             12        1.6183566        0.7803891
#> 50:            NA          117             11        0.5982805        0.2084886
#> 51:            NA            1              1        0.0000000        0.0000000
#> 52:            NA           40              2        0.1169068        0.0487500
#> 53:            NA            3              3        1.0986123        0.6666667
#>     ISS_MAPPED_PP seqCount_sum seqCount_count seqCount_shannon seqCount_simpson
#>             <int>        <num>          <int>            <num>            <num>
#>     seqCount_invsimpson fragmentEstimate_sum fragmentEstimate_count
#>                   <num>                <num>                  <int>
#>  1:            3.024634           201.064659                     15
#>  2:           12.246099           466.804667                     72
#>  3:           10.976410           445.727398                     76
#>  4:           10.431224           354.159459                     67
#>  5:            7.929908           134.338758                     26
#>  6:           26.573601           232.251224                     82
#>  7:           13.584344            86.913928                     41
#>  8:           35.908927           279.440177                     90
#>  9:            5.355015            54.445672                     11
#> 10:           26.765430          1812.911554                     43
#> 11:            5.975309            34.333137                     10
#> 12:           11.474357            98.541329                     31
#> 13:           37.958724           221.918390                     77
#> 14:            9.726378            63.390419                     19
#> 15:            8.902307            20.097692                     13
#> 16:            5.108961            53.231043                     22
#> 17:            6.148615            66.365138                     29
#> 18:            5.527919            26.162592                     10
#> 19:            9.779427            66.290122                     20
#> 20:           16.204669            91.324857                     39
#> 21:            6.320471            46.100979                     32
#> 22:            1.927479            15.029208                     10
#> 23:           49.493862          3547.631728                     78
#> 24:            5.844156             8.007481                      8
#> 25:            1.000000             1.000452                      1
#> 26:            1.000000             1.000310                      1
#> 27:           28.232219           172.660748                     48
#> 28:           40.172734           206.011551                     65
#> 29:           31.572121           218.425224                     74
#> 30:           22.674426           147.531730                     46
#> 31:           24.373976           130.581517                     48
#> 32:           47.123641           256.236461                     87
#> 33:            7.653703            40.176384                     12
#> 34:           27.496256           154.626814                     52
#> 35:           20.775038           166.082193                     35
#> 36:            5.666165            94.906121                     18
#> 37:           30.466007           134.490805                     52
#> 38:            7.057011            57.190086                     36
#> 39:           10.596312            46.149863                     20
#> 40:           10.180735            37.113766                     19
#> 41:           11.149883            37.174678                     24
#> 42:            9.551983            51.164592                     25
#> 43:            5.557141            26.094255                      9
#> 44:           13.917941            55.138395                     26
#> 45:            6.649274            27.081549                     10
#> 46:            5.875887            94.719278                     14
#> 47:            5.405088           106.417161                      8
#> 48:            6.677788            22.061603                      9
#> 49:            4.553509            75.778512                     12
#> 50:            1.263406            12.021148                     11
#> 51:            1.000000             1.001736                      1
#> 52:            1.051248             2.002185                      2
#> 53:            3.000000             3.004764                      3
#>     seqCount_invsimpson fragmentEstimate_sum fragmentEstimate_count
#>                   <num>                <num>                  <int>
#>     fragmentEstimate_shannon fragmentEstimate_simpson
#>                        <num>                    <num>
#>  1:                1.5158967                0.6644034
#>  2:                3.4957593                0.9358864
#>  3:                3.4547233                0.9325986
#>  4:                3.3077507                0.9214312
#>  5:                2.4694808                0.8728520
#>  6:                4.1255479                0.9773086
#>  7:                3.2043690                0.9288182
#>  8:                4.2520301                0.9826514
#>  9:                2.0916345                0.8525232
#> 10:                3.4978977                0.9656616
#> 11:                1.8788583                0.7963169
#> 12:                3.2110253                0.9515061
#> 13:                4.1785406                0.9823708
#> 14:                2.7062808                0.9234798
#> 15:                2.1595245                0.8085066
#> 16:                2.8040916                0.9232304
#> 17:                3.1435438                0.9474470
#> 18:                2.0115316                0.8337413
#> 19:                2.8796549                0.9388584
#> 20:                3.4951040                0.9657862
#> 21:                3.3264870                0.9583635
#> 22:                2.1757058                0.8709946
#> 23:                4.1321436                0.9821639
#> 24:                2.0794415                0.8750000
#> 25:                0.0000000                0.0000000
#> 26:                0.0000000                0.0000000
#> 27:                3.7466139                0.9724395
#> 28:                4.0347118                0.9794325
#> 29:                4.1017170                0.9801725
#> 30:                3.7422290                0.9739783
#> 31:                3.6532389                0.9690120
#> 32:                4.2743192                0.9837312
#> 33:                2.3602980                0.8907893
#> 34:                3.7796517                0.9736284
#> 35:                3.3921932                0.9600307
#> 36:                2.5886973                0.9080593
#> 37:                3.7994565                0.9747936
#> 38:                3.4030063                0.9589507
#> 39:                2.8516684                0.9347119
#> 40:                2.7383206                0.9217140
#> 41:                2.9976077                0.9391220
#> 42:                3.0478044                0.9449613
#> 43:                2.0881376                0.8668208
#> 44:                3.1241984                0.9513700
#> 45:                2.2105119                0.8833528
#> 46:                2.0462669                0.8345178
#> 47:                1.9118719                0.8417363
#> 48:                2.1289383                0.8759944
#> 49:                2.0395769                0.8360567
#> 50:                2.3693447                0.9027688
#> 51:                0.0000000                0.0000000
#> 52:                0.6931471                0.5000000
#> 53:                1.0986117                0.6666663
#>     fragmentEstimate_shannon fragmentEstimate_simpson
#>                        <num>                    <num>
#>     fragmentEstimate_invsimpson seqCount_describe_vars seqCount_describe_n
#>                           <num>                  <num>               <num>
#>  1:                    2.979768                      1                  15
#>  2:                   15.597315                      1                  72
#>  3:                   14.836491                      1                  76
#>  4:                   12.727698                      1                  67
#>  5:                    7.864852                      1                  26
#>  6:                   44.069549                      1                  82
#>  7:                   14.048533                      1                  41
#>  8:                   57.641549                      1                  90
#>  9:                    6.780730                      1                  11
#> 10:                   29.121899                      1                  43
#> 11:                    4.909587                      1                  10
#> 12:                   20.621135                      1                  31
#> 13:                   56.724054                      1                  77
#> 14:                   13.068445                      1                  19
#> 15:                    5.222113                      1                  13
#> 16:                   13.025989                      1                  22
#> 17:                   19.028416                      1                  29
#> 18:                    6.014722                      1                  10
#> 19:                   16.355473                      1                  20
#> 20:                   29.227978                      1                  39
#> 21:                   24.017379                      1                  32
#> 22:                    7.751614                      1                  10
#> 23:                   56.066071                      1                  78
#> 24:                    7.999999                      1                   8
#> 25:                    1.000000                      1                   1
#> 26:                    1.000000                      1                   1
#> 27:                   36.283780                      1                  48
#> 28:                   48.620400                      1                  65
#> 29:                   50.434894                      1                  74
#> 30:                   38.429519                      1                  46
#> 31:                   32.270520                      1                  48
#> 32:                   61.467277                      1                  87
#> 33:                    9.156610                      1                  12
#> 34:                   37.919638                      1                  52
#> 35:                   25.019190                      1                  35
#> 36:                   10.876572                      1                  18
#> 37:                   39.672460                      1                  52
#> 38:                   24.360939                      1                  36
#> 39:                   15.316731                      1                  20
#> 40:                   12.773677                      1                  19
#> 41:                   16.426308                      1                  24
#> 42:                   18.169039                      1                  25
#> 43:                    7.508679                      1                   9
#> 44:                   20.563431                      1                  26
#> 45:                    8.572855                      1                  10
#> 46:                    6.042945                      1                  14
#> 47:                    6.318567                      1                   8
#> 48:                    8.064152                      1                   9
#> 49:                    6.099669                      1                  12
#> 50:                   10.284764                      1                  11
#> 51:                    1.000000                      1                   1
#> 52:                    2.000000                      1                   2
#> 53:                    2.999996                      1                   3
#>     fragmentEstimate_invsimpson seqCount_describe_vars seqCount_describe_n
#>                           <num>                  <num>               <num>
#>     seqCount_describe_mean seqCount_describe_sd seqCount_describe_median
#>                      <num>                <num>                    <num>
#>  1:              23.800000            49.019238                      4.0
#>  2:            1803.708333          4012.246817                    903.0
#>  3:            1019.394737          2497.605611                    455.0
#>  4:            1203.134328          2822.931234                    420.0
#>  5:               5.653846             8.703757                      2.0
#>  6:             620.231707           901.262931                    442.0
#>  7:               2.707317             3.893866                      1.0
#>  8:             598.744444           738.974290                    411.0
#>  9:           11027.636364         11874.895623                  11559.0
#> 10:             388.906977           306.470446                    334.0
#> 11:               4.400000             3.806427                      3.0
#> 12:            2448.548387          3246.888632                   1489.0
#> 13:             749.142857           764.732250                    547.0
#> 14:            6340.631579          6360.953101                   6249.0
#> 15:               6.230769             4.399883                      4.0
#> 16:             137.681818           256.236012                      5.0
#> 17:            4009.793103          7867.007184                      2.0
#> 18:               3.300000             3.128720                      1.5
#> 19:            2872.150000          3012.494250                   2149.5
#> 20:            1063.794872          1278.207239                    275.0
#> 21:             115.531250           236.598746                      1.0
#> 22:              55.900000           120.586944                      1.0
#> 23:             540.025641           412.486599                    565.5
#> 24:               3.750000             2.434866                      4.0
#> 25:               1.000000                   NA                      1.0
#> 26:               1.000000                   NA                      1.0
#> 27:            1064.875000           900.485605                    860.5
#> 28:             115.630769            91.609151                    102.0
#> 29:              80.783784            94.287126                     47.5
#> 30:            1039.478261          1065.948169                    763.5
#> 31:             104.458333           103.931434                     73.5
#> 32:             105.172414            97.308453                     79.0
#> 33:             357.250000           281.183901                    368.5
#> 34:             100.730769            96.019221                     73.5
#> 35:            1297.600000          1089.406698                   1061.0
#> 36:              81.555556           123.814008                     16.5
#> 37:              63.076923            53.547797                     52.0
#> 38:              28.583333            58.707203                      1.0
#> 39:              92.550000            89.451236                     73.5
#> 40:             109.947368           105.136141                     92.0
#> 41:               2.875000             3.152811                      1.0
#> 42:              50.840000            65.987170                     24.0
#> 43:             143.444444           119.755074                    111.0
#> 44:             629.230769           597.873920                    425.0
#> 45:              56.500000            42.277522                     59.5
#> 46:             109.857143           134.051573                      4.0
#> 47:             244.750000           181.291990                    317.0
#> 48:              73.444444            45.937760                     58.0
#> 49:             115.666667           154.491738                     12.5
#> 50:              10.636364            30.968606                      1.0
#> 51:               1.000000                   NA                      1.0
#> 52:              20.000000            26.870058                     20.0
#> 53:               1.000000             0.000000                      1.0
#>     seqCount_describe_mean seqCount_describe_sd seqCount_describe_median
#>                      <num>                <num>                    <num>
#>     seqCount_describe_trimmed seqCount_describe_mad seqCount_describe_min
#>                         <num>                 <num>                 <num>
#>  1:                 13.384615                4.4478                     1
#>  2:                948.258621              876.2166                     1
#>  3:                485.274194              582.6618                     1
#>  4:                556.254545              517.4274                     1
#>  5:                  3.772727                1.4826                     1
#>  6:                469.787879              544.1142                     1
#>  7:                  1.787879                0.0000                     1
#>  8:                456.513889              583.4031                     1
#>  9:               9547.666667            15758.5554                     1
#> 10:                368.628571              314.3112                     2
#> 11:                  3.875000                2.2239                     1
#> 12:               1813.800000             2204.6262                     1
#> 13:                639.126984              621.2094                     1
#> 14:               6057.764706             9263.2848                     1
#> 15:                  6.090909                4.4478                     1
#> 16:                 82.111111                5.9304                     1
#> 17:               2912.120000                1.4826                     1
#> 18:                  2.875000                0.7413                     1
#> 19:               2500.312500             3169.0575                     1
#> 20:                930.606061              406.2324                     1
#> 21:                 61.615385                0.0000                     1
#> 22:                 25.125000                0.0000                     1
#> 23:                520.734375              440.3322                     1
#> 24:                  3.750000                2.2239                     1
#> 25:                  1.000000                0.0000                     1
#> 26:                  1.000000                0.0000                     1
#> 27:                969.150000              873.2514                     3
#> 28:                102.905660               72.6474                     1
#> 29:                 66.366667               57.8214                     1
#> 30:                884.184211              960.7248                     2
#> 31:                 89.400000               72.6474                     1
#> 32:                 94.169014              100.8168                     1
#> 33:                344.200000              366.9435                    10
#> 34:                 87.047619               71.1648                     1
#> 35:               1187.586207             1068.9546                    38
#> 36:                 64.500000               22.9803                     1
#> 37:                 56.666667               53.3736                     1
#> 38:                 15.133333                0.0000                     1
#> 39:                 78.062500               74.8713                     1
#> 40:                102.941176              114.1602                     1
#> 41:                  2.200000                0.0000                     1
#> 42:                 40.190476               29.6520                     1
#> 43:                143.444444              100.8168                     1
#> 44:                575.000000              510.7557                     2
#> 45:                 56.250000               58.5627                     1
#> 46:                 99.916667                4.4478                     1
#> 47:                244.750000              174.2055                    12
#> 48:                 73.444444               50.4084                     1
#> 49:                101.200000               17.0499                     1
#> 50:                  1.333333                0.0000                     1
#> 51:                  1.000000                0.0000                     1
#> 52:                 20.000000               28.1694                     1
#> 53:                  1.000000                0.0000                     1
#>     seqCount_describe_trimmed seqCount_describe_mad seqCount_describe_min
#>                         <num>                 <num>                 <num>
#>     seqCount_describe_max seqCount_describe_range seqCount_describe_skew
#>                     <num>                   <num>                  <num>
#>  1:                   182                     181            2.301241908
#>  2:                 23219                   23218            4.241584903
#>  3:                 14748                   14747            4.426622167
#>  4:                 17030                   17029            4.090477930
#>  5:                    32                      31            2.045164619
#>  6:                  6831                    6830            4.440115101
#>  7:                    21                      20            3.479912926
#>  8:                  3420                    3419            1.802565558
#>  9:                 35374                   35373            0.542170566
#> 10:                   980                     978            0.445493428
#> 11:                    12                      11            0.919078343
#> 12:                 15387                   15386            2.203411709
#> 13:                  3099                    3098            1.189911727
#> 14:                 17489                   17488            0.258600919
#> 15:                    13                      12            0.214947146
#> 16:                   899                     898            1.695818613
#> 17:                 23099                   23098            1.433389491
#> 18:                     9                       8            0.729559680
#> 19:                  9255                    9254            0.665725699
#> 20:                  4081                    4080            0.764076544
#> 21:                   931                     930            1.904838107
#> 22:                   357                     356            1.610560306
#> 23:                  2776                    2775            1.867297249
#> 24:                     8                       7            0.305241830
#> 25:                     1                       0                     NA
#> 26:                     1                       0                     NA
#> 27:                  3447                    3444            0.872205813
#> 28:                   485                     484            1.803069328
#> 29:                   602                     601            2.602506887
#> 30:                  4633                    4631            1.365036835
#> 31:                   492                     491            1.598850712
#> 32:                   399                     398            0.920512611
#> 33:                   835                     825            0.130786772
#> 34:                   485                     484            1.638509198
#> 35:                  4600                    4562            0.970808861
#> 36:                   435                     434            1.452104290
#> 37:                   224                     223            0.984500650
#> 38:                   263                     262            2.380294965
#> 39:                   300                     299            1.081245294
#> 40:                   338                     337            0.715416659
#> 41:                    12                      11            1.721437804
#> 42:                   230                     229            1.460662654
#> 43:                   379                     378            0.609430398
#> 44:                  1847                    1845            0.814842416
#> 45:                   114                     113           -0.003001335
#> 46:                   338                     337            0.446900149
#> 47:                   446                     434           -0.275019641
#> 48:                   163                     162            0.381570108
#> 49:                   375                     374            0.656807751
#> 50:                   104                     103            2.465923327
#> 51:                     1                       0                     NA
#> 52:                    39                      38            0.000000000
#> 53:                     1                       0                    NaN
#>     seqCount_describe_max seqCount_describe_range seqCount_describe_skew
#>                     <num>                   <num>                  <num>
#>     seqCount_describe_kurtosis seqCount_describe_se
#>                          <num>                <num>
#>  1:                  4.4006968           12.6567129
#>  2:                 17.5610094          472.8478220
#>  3:                 19.2823139          286.4950121
#>  4:                 17.1580681          344.8759404
#>  5:                  3.0002933            1.7069472
#>  6:                 26.1076318           99.5278409
#>  7:                 12.2312876            0.6081198
#>  8:                  3.3451795           77.8947296
#>  9:                 -1.0023479         3580.4157460
#> 10:                 -1.0232860           46.7363048
#> 11:                 -0.7765928            1.2036980
#> 12:                  5.7210093          583.1584134
#> 13:                  0.6675981           87.1493280
#> 14:                 -1.6185786         1459.3027238
#> 15:                 -1.7809891            1.2203081
#> 16:                  1.6582823           54.6297014
#> 17:                  0.1905641         1460.8665594
#> 18:                 -1.3581940            0.9893881
#> 19:                 -0.9237911          673.6141925
#> 20:                 -0.9036271          204.6769654
#> 21:                  2.7235555           41.8251445
#> 22:                  1.0549935           38.1329400
#> 23:                  8.7391643           46.7049103
#> 24:                 -1.3723939            0.8608551
#> 25:                         NA                   NA
#> 26:                         NA                   NA
#> 27:                 -0.1202778          129.9739016
#> 28:                  4.1475510           11.3627167
#> 29:                 10.5808065           10.9606558
#> 30:                  1.6176863          157.1654832
#> 31:                  2.5227001           15.0012103
#> 32:                  0.2614462           10.4325612
#> 33:                 -1.5703281           81.1708006
#> 34:                  3.3180379           13.3154702
#> 35:                  0.4571943          184.1433411
#> 36:                  1.1279457           29.1832415
#> 37:                  0.3744349            7.4257433
#> 38:                  5.4638897            9.7845339
#> 39:                  0.2435027           20.0019045
#> 40:                 -0.6542499           24.1198850
#> 41:                  1.9755420            0.6435649
#> 42:                  0.8237687           13.1974341
#> 43:                 -0.9585728           39.9183580
#> 44:                 -0.7629889          117.2527226
#> 45:                 -1.6485903           13.3693264
#> 46:                 -1.7003497           35.8267898
#> 47:                 -1.9433976           64.0963978
#> 48:                 -0.6556674           15.3125866
#> 49:                 -1.5340806           44.5979231
#> 50:                  4.5174735            9.3373860
#> 51:                         NA                   NA
#> 52:                 -2.7500000           19.0000000
#> 53:                        NaN            0.0000000
#>     seqCount_describe_kurtosis seqCount_describe_se
#>                          <num>                <num>
#>     fragmentEstimate_describe_vars fragmentEstimate_describe_n
#>                              <num>                       <num>
#>  1:                              1                          15
#>  2:                              1                          72
#>  3:                              1                          76
#>  4:                              1                          67
#>  5:                              1                          26
#>  6:                              1                          82
#>  7:                              1                          41
#>  8:                              1                          90
#>  9:                              1                          11
#> 10:                              1                          43
#> 11:                              1                          10
#> 12:                              1                          31
#> 13:                              1                          77
#> 14:                              1                          19
#> 15:                              1                          13
#> 16:                              1                          22
#> 17:                              1                          29
#> 18:                              1                          10
#> 19:                              1                          20
#> 20:                              1                          39
#> 21:                              1                          32
#> 22:                              1                          10
#> 23:                              1                          78
#> 24:                              1                           8
#> 25:                              1                           1
#> 26:                              1                           1
#> 27:                              1                          48
#> 28:                              1                          65
#> 29:                              1                          74
#> 30:                              1                          46
#> 31:                              1                          48
#> 32:                              1                          87
#> 33:                              1                          12
#> 34:                              1                          52
#> 35:                              1                          35
#> 36:                              1                          18
#> 37:                              1                          52
#> 38:                              1                          36
#> 39:                              1                          20
#> 40:                              1                          19
#> 41:                              1                          24
#> 42:                              1                          25
#> 43:                              1                           9
#> 44:                              1                          26
#> 45:                              1                          10
#> 46:                              1                          14
#> 47:                              1                           8
#> 48:                              1                           9
#> 49:                              1                          12
#> 50:                              1                          11
#> 51:                              1                           1
#> 52:                              1                           2
#> 53:                              1                           3
#>     fragmentEstimate_describe_vars fragmentEstimate_describe_n
#>                              <num>                       <num>
#>     fragmentEstimate_describe_mean fragmentEstimate_describe_sd
#>                              <num>                        <num>
#>  1:                      13.404311                 2.786707e+01
#>  2:                       6.483398                 1.241552e+01
#>  3:                       5.864834                 1.198705e+01
#>  4:                       5.285962                 1.099775e+01
#>  5:                       5.166875                 8.001288e+00
#>  6:                       2.832332                 2.643828e+00
#>  7:                       2.119852                 2.972646e+00
#>  8:                       3.104891                 2.339368e+00
#>  9:                       4.949607                 4.094939e+00
#> 10:                      42.160734                 2.944916e+01
#> 11:                       3.433314                 3.685074e+00
#> 12:                       3.178753                 2.292427e+00
#> 13:                       2.882057                 1.734395e+00
#> 14:                       3.336338                 2.309314e+00
#> 15:                       1.545976                 1.963774e+00
#> 16:                       2.419593                 2.055568e+00
#> 17:                       2.288453                 1.685942e+00
#> 18:                       2.616259                 2.244818e+00
#> 19:                       3.314506                 1.605262e+00
#> 20:                       2.341663                 1.371696e+00
#> 21:                       1.440656                 8.438481e-01
#> 22:                       1.502921                 8.532066e-01
#> 23:                      45.482458                 2.863215e+01
#> 24:                       1.000935                 4.346800e-04
#> 25:                       1.000452                           NA
#> 26:                       1.000310                           NA
#> 27:                       3.597099                 2.065673e+00
#> 28:                       3.169408                 1.853904e+00
#> 29:                       2.951692                 2.031396e+00
#> 30:                       3.207212                 1.439228e+00
#> 31:                       2.720448                 1.919404e+00
#> 32:                       2.945247                 1.909232e+00
#> 33:                       3.348032                 1.948654e+00
#> 34:                       2.973593                 1.829670e+00
#> 35:                       4.745206                 3.040856e+00
#> 36:                       5.272562                 4.390682e+00
#> 37:                       2.586362                 1.455792e+00
#> 38:                       1.588614                 1.113647e+00
#> 39:                       2.307493                 1.309092e+00
#> 40:                       1.953356                 1.401134e+00
#> 41:                       1.548945                 1.074389e+00
#> 42:                       2.046584                 1.280763e+00
#> 43:                       2.899362                 1.370511e+00
#> 44:                       2.120707                 1.112019e+00
#> 45:                       2.708155                 1.164725e+00
#> 46:                       6.765663                 8.056651e+00
#> 47:                      13.302145                 7.335810e+00
#> 48:                       2.451289                 8.857151e-01
#> 49:                       6.314876                 6.487004e+00
#> 50:                       1.092832                 3.022575e-01
#> 51:                       1.001736                           NA
#> 52:                       1.001093                 3.673758e-04
#> 53:                       1.001588                 1.338704e-03
#>     fragmentEstimate_describe_mean fragmentEstimate_describe_sd
#>                              <num>                        <num>
#>     fragmentEstimate_describe_median fragmentEstimate_describe_trimmed
#>                                <num>                             <num>
#>  1:                         3.012512                          7.470644
#>  2:                         3.012057                          3.811294
#>  3:                         3.009877                          3.142067
#>  4:                         3.008092                          2.737114
#>  5:                         1.003097                          3.513520
#>  6:                         2.505033                          2.324454
#>  7:                         1.002015                          1.398289
#>  8:                         3.007780                          2.718686
#>  9:                         4.015796                          4.581450
#> 10:                        39.036122                         40.562169
#> 11:                         2.004800                          2.772710
#> 12:                         2.005189                          2.811093
#> 13:                         3.007980                          2.659478
#> 14:                         4.008330                          3.194304
#> 15:                         1.001018                          1.001403
#> 16:                         1.001317                          2.006402
#> 17:                         2.004555                          2.048540
#> 18:                         1.503943                          2.262289
#> 19:                         3.010762                          3.262810
#> 20:                         2.005080                          2.249665
#> 21:                         1.001360                          1.310164
#> 22:                         1.001153                          1.377210
#> 23:                        42.000196                         44.901123
#> 24:                         1.001025                          1.000935
#> 25:                         1.000452                          1.000452
#> 26:                         1.000310                          1.000310
#> 27:                         3.010014                          3.335589
#> 28:                         3.009442                          2.896971
#> 29:                         2.505801                          2.611686
#> 30:                         3.009387                          3.113838
#> 31:                         2.005851                          2.433007
#> 32:                         3.008802                          2.671993
#> 33:                         3.007868                          3.007499
#> 34:                         3.007741                          2.723031
#> 35:                         3.011246                          4.200535
#> 36:                         4.032253                          4.801168
#> 37:                         2.506984                          2.436239
#> 38:                         1.001915                          1.369925
#> 39:                         2.004245                          2.131404
#> 40:                         1.001656                          1.769694
#> 41:                         1.002144                          1.354908
#> 42:                         2.004350                          1.862934
#> 43:                         3.012131                          2.899362
#> 44:                         2.003264                          2.050231
#> 45:                         3.010009                          2.758162
#> 46:                         2.009365                          5.674497
#> 47:                        15.161805                         13.302145
#> 48:                         3.006917                          2.451289
#> 49:                         4.079891                          5.237292
#> 50:                         1.001625                          1.001750
#> 51:                         1.001736                          1.001736
#> 52:                         1.001093                          1.001093
#> 53:                         1.001111                          1.001588
#>     fragmentEstimate_describe_median fragmentEstimate_describe_trimmed
#>                                <num>                             <num>
#>     fragmentEstimate_describe_mad fragmentEstimate_describe_min
#>                             <num>                         <num>
#>  1:                  2.981298e+00                      1.000572
#>  2:                  2.239587e+00                      1.000302
#>  3:                  2.976246e+00                      1.000270
#>  4:                  1.499789e+00                      1.000251
#>  5:                  3.649738e-03                      1.000438
#>  6:                  7.510140e-01                      1.000188
#>  7:                  1.265862e-03                      1.000363
#>  8:                  2.974732e+00                      1.000237
#>  9:                  4.468958e+00                      1.000326
#> 10:                  3.926041e+01                      2.005726
#> 11:                  1.488527e+00                      1.000361
#> 12:                  1.488346e+00                      1.001245
#> 13:                  1.492398e+00                      1.000395
#> 14:                  4.457778e+00                      1.000797
#> 15:                  8.451267e-04                      1.000448
#> 16:                  1.496358e-03                      1.000189
#> 17:                  1.487282e+00                      1.000170
#> 18:                  7.462303e-01                      1.000481
#> 19:                  1.490373e+00                      1.000878
#> 20:                  1.488259e+00                      1.001013
#> 21:                  2.548275e-04                      1.000364
#> 22:                  3.453061e-04                      1.000484
#> 23:                  3.203641e+01                      1.000095
#> 24:                  5.625647e-04                      1.000271
#> 25:                  0.000000e+00                      1.000452
#> 26:                  0.000000e+00                      1.000310
#> 27:                  1.204654e-02                      1.000891
#> 28:                  1.487735e+00                      1.000853
#> 29:                  2.229894e+00                      1.000748
#> 30:                  9.548069e-03                      1.000647
#> 31:                  1.489289e+00                      1.000233
#> 32:                  1.498433e+00                      1.000451
#> 33:                  3.832020e-03                      1.001297
#> 34:                  1.489516e+00                      1.000130
#> 35:                  7.114164e-03                      1.704548
#> 36:                  4.291451e+00                      1.000211
#> 37:                  2.231648e+00                      1.000260
#> 38:                  9.476613e-04                      1.000654
#> 39:                  1.487938e+00                      1.000463
#> 40:                  1.351624e-03                      1.000429
#> 41:                  1.444463e-03                      1.000319
#> 42:                  1.487306e+00                      1.000216
#> 43:                  1.492358e+00                      1.000102
#> 44:                  1.486179e+00                      1.000286
#> 45:                  1.488739e+00                      1.000356
#> 46:                  1.495089e+00                      1.000901
#> 47:                  3.344285e+00                      2.004955
#> 48:                  1.028171e-02                      1.000461
#> 49:                  4.564592e+00                      1.000782
#> 50:                  5.953127e-04                      1.001224
#> 51:                  0.000000e+00                      1.001736
#> 52:                  3.851408e-04                      1.000833
#> 53:                  8.263112e-04                      1.000553
#>     fragmentEstimate_describe_mad fragmentEstimate_describe_min
#>                             <num>                         <num>
#>     fragmentEstimate_describe_max fragmentEstimate_describe_range
#>                             <num>                           <num>
#>  1:                    102.945718                    1.019451e+02
#>  2:                     68.737467                    6.773717e+01
#>  3:                     65.157600                    6.415733e+01
#>  4:                     60.847813                    5.984756e+01
#>  5:                     28.628174                    2.762774e+01
#>  6:                     18.345126                    1.734494e+01
#>  7:                     15.439133                    1.443877e+01
#>  8:                     13.173494                    1.217326e+01
#>  9:                     12.212294                    1.121197e+01
#> 10:                    102.491428                    1.004857e+02
#> 11:                     11.151093                    1.015073e+01
#> 12:                      9.105808                    8.104562e+00
#> 13:                      9.105112                    8.104718e+00
#> 14:                      8.086461                    7.085664e+00
#> 15:                      8.081811                    7.081363e+00
#> 16:                      8.065660                    7.065471e+00
#> 17:                      7.078347                    6.078178e+00
#> 18:                      7.063803                    6.063322e+00
#> 19:                      6.042502                    5.041624e+00
#> 20:                      5.034588                    4.033575e+00
#> 21:                      3.013021                    2.012657e+00
#> 22:                      3.011042                    2.010558e+00
#> 23:                    100.161988                    9.916189e+01
#> 24:                      1.001434                    1.163127e-03
#> 25:                      1.000452                    0.000000e+00
#> 26:                      1.000310                    0.000000e+00
#> 27:                     14.083418                    1.308253e+01
#> 28:                     11.158685                    1.015783e+01
#> 29:                     10.158522                    9.157774e+00
#> 30:                     10.115103                    9.114456e+00
#> 31:                      9.119998                    8.119765e+00
#> 32:                      9.101815                    8.101365e+00
#> 33:                      9.100094                    8.098797e+00
#> 34:                      9.084035                    8.083906e+00
#> 35:                     12.600475                    1.089593e+01
#> 36:                     17.087221                    1.608701e+01
#> 37:                      7.049681                    6.049421e+00
#> 38:                      6.046701                    5.046047e+00
#> 39:                      6.029366                    5.028903e+00
#> 40:                      6.028539                    5.028110e+00
#> 41:                      5.054662                    4.054343e+00
#> 42:                      5.022865                    4.022649e+00
#> 43:                      5.015339                    4.015236e+00
#> 44:                      4.016451                    3.016165e+00
#> 45:                      4.015899                    3.015543e+00
#> 46:                     25.624419                    2.462352e+01
#> 47:                     23.579866                    2.157491e+01
#> 48:                      3.014156                    2.013695e+00
#> 49:                     22.404813                    2.140403e+01
#> 50:                      2.004172                    1.002948e+00
#> 51:                      1.001736                    0.000000e+00
#> 52:                      1.001352                    5.195478e-04
#> 53:                      1.003100                    2.546580e-03
#>     fragmentEstimate_describe_max fragmentEstimate_describe_range
#>                             <num>                           <num>
#>     fragmentEstimate_describe_skew fragmentEstimate_describe_kurtosis
#>                              <num>                              <num>
#>  1:                   2.316325e+00                          4.3379571
#>  2:                   4.200284e+00                         17.0134959
#>  3:                   4.164724e+00                         16.6881143
#>  4:                   4.174675e+00                         17.0548285
#>  5:                   1.953187e+00                          2.5068058
#>  6:                   3.564088e+00                         15.7237949
#>  7:                   3.465180e+00                         11.6912367
#>  8:                   1.472620e+00                          2.7446859
#>  9:                   6.186531e-01                         -1.2808229
#> 10:                   3.158471e-01                         -1.1336562
#> 11:                   1.157340e+00                         -0.4195354
#> 12:                   1.196926e+00                          0.5337682
#> 13:                   1.209707e+00                          1.7311584
#> 14:                   3.529692e-01                         -1.1813391
#> 15:                   2.816170e+00                          6.4437855
#> 16:                   1.401263e+00                          1.0168970
#> 17:                   1.327987e+00                          0.8065693
#> 18:                   9.507548e-01                         -0.7724153
#> 19:                   2.800850e-01                         -0.9652756
#> 20:                   3.757505e-01                         -1.3774092
#> 21:                   1.297396e+00                         -0.3239794
#> 22:                   9.780605e-01                         -0.9618423
#> 23:                   1.839668e-01                         -1.0980930
#> 24:                  -2.284244e-01                         -1.6953538
#> 25:                             NA                                 NA
#> 26:                             NA                                 NA
#> 27:                   2.904116e+00                         11.6588565
#> 28:                   2.210633e+00                          6.5758780
#> 29:                   1.537537e+00                          2.3768799
#> 30:                   2.228271e+00                          9.3082986
#> 31:                   1.320622e+00                          1.7397660
#> 32:                   1.171759e+00                          1.2772169
#> 33:                   1.970219e+00                          3.5225017
#> 34:                   1.188657e+00                          1.4823160
#> 35:                   1.543861e+00                          0.9973010
#> 36:                   1.138082e+00                          0.6679724
#> 37:                   7.519567e-01                          0.2093692
#> 38:                   2.211101e+00                          5.1257353
#> 39:                   9.645585e-01                          0.7892475
#> 40:                   1.373334e+00                          1.2086107
#> 41:                   1.779517e+00                          2.3609818
#> 42:                   9.728921e-01                         -0.1698843
#> 43:                  -9.170249e-02                         -1.4314296
#> 44:                   2.941909e-01                         -1.4700960
#> 45:                  -2.468475e-01                         -1.5748596
#> 46:                   1.065480e+00                         -0.2502699
#> 47:                  -4.122360e-01                         -1.3069937
#> 48:                  -8.236812e-01                         -1.2891641
#> 49:                   1.107043e+00                          0.3805447
#> 50:                   2.466904e+00                          4.5206378
#> 51:                             NA                                 NA
#> 52:                  -4.533028e-13                         -2.7500000
#> 53:                   3.112158e-01                         -2.3333333
#>     fragmentEstimate_describe_skew fragmentEstimate_describe_kurtosis
#>                              <num>                              <num>
#>     fragmentEstimate_describe_se   nIS
#>                            <num> <int>
#>  1:                 7.1952466279    15
#>  2:                 1.4631824759    72
#>  3:                 1.3750095083    76
#>  4:                 1.3435884856    67
#>  5:                 1.5691817431    26
#>  6:                 0.2919620099    82
#>  7:                 0.4642492918    41
#>  8:                 0.2465910846    90
#>  9:                 1.2346706770    11
#> 10:                 4.4909549219    43
#> 11:                 1.1653227353    10
#> 12:                 0.4117321131    31
#> 13:                 0.1976526781    77
#> 14:                 0.5297928780    19
#> 15:                 0.5446529515    13
#> 16:                 0.4382486380    22
#> 17:                 0.3130716475    29
#> 18:                 0.7098736619    10
#> 19:                 0.3589475809    20
#> 20:                 0.2196470704    39
#> 21:                 0.1491726757    32
#> 22:                 0.2698076248    10
#> 23:                 3.2419530711    78
#> 24:                 0.0001536826     8
#> 25:                           NA     1
#> 26:                           NA     1
#> 27:                 0.2981542158    48
#> 28:                 0.2299485144    65
#> 29:                 0.2361450154    74
#> 30:                 0.2122025482    46
#> 31:                 0.2770420812    48
#> 32:                 0.2046911329    87
#> 33:                 0.5625279552    12
#> 34:                 0.2537296307    52
#> 35:                 0.5139984901    35
#> 36:                 1.0348937038    18
#> 37:                 0.2018820919    52
#> 38:                 0.1856077872    36
#> 39:                 0.2927218410    20
#> 40:                 0.3214422872    19
#> 41:                 0.2193087193    24
#> 42:                 0.2561526253    25
#> 43:                 0.4568368515     9
#> 44:                 0.2180848137    26
#> 45:                 0.3683184518    10
#> 46:                 2.1532305518    14
#> 47:                 2.5936005456     8
#> 48:                 0.2952383720     9
#> 49:                 1.8726366147    12
#> 50:                 0.0911340761    11
#> 51:                           NA     1
#> 52:                 0.0002597739     2
#> 53:                 0.0007729012     3
#>     fragmentEstimate_describe_se   nIS
#>                            <num> <int>
#> 
```
