# Filter out outliers in metadata, identified by the chosen outlier test.

**\[experimental\]** Filter out outliers in metadata by using
appropriate outlier tests.

## Usage

``` r
outlier_filter(
  metadata,
  pcr_id_col = pcr_id_column(),
  outlier_test = c(outliers_by_pool_fragments),
  outlier_test_outputs = NULL,
  combination_logic = c("AND"),
  negate = FALSE,
  report_path = default_report_path(),
  ...
)
```

## Arguments

- metadata:

  The metadata data frame

- pcr_id_col:

  The name of the pcr identifier column

- outlier_test:

  One or more outlier tests. Must be functions, either from
  [`available_outlier_tests()`](https://calabrialab.github.io/ISAnalytics/dev/reference/available_outlier_tests.md)
  or custom functions that produce an appropriate output format (see
  details).

- outlier_test_outputs:

  `NULL`, a data frame or a list of data frames. See details.

- combination_logic:

  One or more logical operators ("AND", "OR", "XOR", "NAND", "NOR",
  "XNOR"). See datails.

- negate:

  If `TRUE` will return only the metadata that was flagged to be
  removed. If `FALSE` will return only the metadata that wasn't flagged
  to be removed.

- report_path:

  The path where the report file should be saved. Can be a folder or
  `NULL` if no report should be produced. Defaults to
  `{user_home}/ISAnalytics_reports`.

- ...:

  Additional named arguments passed to `outliers_test`

## Value

A data frame of metadata which has less or the same amount of rows

## Details

### Modular structure

The outlier filtering functions are structured in a modular fashion.
There are 2 kind of functions:

- Outlier tests - Functions that perform some kind of calculation based
  on inputs and flags metadata

- Outlier filter - A function that takes one or more outlier tests,
  combines all the flags with a given logic and filters out rows that
  are flagged as outliers

This function acts as the filter. It can either take one or more outlier
tests as functions and call them through the argument `outlier_test`, or
it can take directly outputs produced by individual tests in the
argument `outlier_test_outputs` - if both are provided the second one
has priority. The second method offers a bit more freedom, since single
tests can be run independently and intermediate results saved and
examined more in detail. If more than one test is to be performed, the
argument `combination_logic` tells the function how to combine the
flags: you can specify 1 logical operator or more than 1, provided it is
compatible with the number of tests.

### Writing custom outlier tests

You have the freedom to provide your own functions as outlier tests. For
this purpose, functions provided must respect this guidelines:

- Must take as input the whole metadata df

- Must return a df containing AT LEAST the `pcr_id_col` and a logical
  column `"to_remove"` that contains the flag

- The `pcr_id_col` must contain all the values originally present in the
  metadata df

## See also

Other Data cleaning and pre-processing:
[`aggregate_metadata()`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_metadata.md),
[`aggregate_values_by_key()`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_values_by_key.md),
[`compute_near_integrations()`](https://calabrialab.github.io/ISAnalytics/dev/reference/compute_near_integrations.md),
[`default_meta_agg()`](https://calabrialab.github.io/ISAnalytics/dev/reference/default_meta_agg.md),
[`outliers_by_pool_fragments()`](https://calabrialab.github.io/ISAnalytics/dev/reference/outliers_by_pool_fragments.md),
[`purity_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/purity_filter.md),
[`realign_after_collisions()`](https://calabrialab.github.io/ISAnalytics/dev/reference/realign_after_collisions.md),
[`remove_collisions()`](https://calabrialab.github.io/ISAnalytics/dev/reference/remove_collisions.md),
[`threshold_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/threshold_filter.md)

## Examples

``` r
data("association_file", package = "ISAnalytics")
filtered_af <- outlier_filter(association_file,
    key = "BARCODE_MUX",
    report_path = NULL
)
#> Removing NAs from data...
#> Log2 transformation, removing values <= 0
head(filtered_af)
#>    ProjectID  FUSIONID PoolID TagSequence SubjectID VectorType VectorID
#>       <char>    <char> <char>      <char>    <char>     <char>   <char>
#> 1:      PJ01 ET#382.46 POOL01   LTR75LC38     PT001      lenti    GLOBE
#> 2:      PJ01 ET#381.40 POOL01   LTR53LC32     PT001      lenti    GLOBE
#> 3:      PJ01  ET#381.9 POOL01   LTR83LC66     PT001      lenti    GLOBE
#> 4:      PJ01 ET#381.71 POOL01   LTR27LC94     PT001      lenti    GLOBE
#> 5:      PJ01  ET#381.2 POOL01   LTR69LC52     PT001      lenti    GLOBE
#> 6:      PJ01 ET#382.28 POOL01    LTR37LC2     PT001      lenti    GLOBE
#>    ExperimentID Tissue TimePoint DNAFragmentation PCRMethod TagIDextended
#>          <char> <char>    <char>           <char>    <char>        <char>
#> 1:         <NA>     PB      0060            SONIC      SLiM     LTR75LC38
#> 2:         <NA>     BM      0180            SONIC      SLiM     LTR53LC32
#> 3:         <NA>     BM      0180            SONIC      SLiM     LTR83LC66
#> 4:         <NA>     BM      0180            SONIC      SLiM     LTR27LC94
#> 5:         <NA>     PB      0180            SONIC      SLiM     LTR69LC52
#> 6:         <NA>     BM      0060            SONIC      SLiM      LTR37LC2
#>    Keywords CellMarker      TagID NGSProvider NGSTechnology ConverrtedFilesDir
#>      <char>     <char>     <char>      <char>        <char>             <char>
#> 1:     <NA>        MNC LTR75.LC38        <NA>         HiSeq               <NA>
#> 2:     <NA>        MNC LTR53.LC32        <NA>         HiSeq               <NA>
#> 3:     <NA>        MNC LTR83.LC66        <NA>         HiSeq               <NA>
#> 4:     <NA>        MNC LTR27.LC94        <NA>         HiSeq               <NA>
#> 5:     <NA>        MNC LTR69.LC52        <NA>         HiSeq               <NA>
#> 6:     <NA>        MNC  LTR37.LC2        <NA>         HiSeq               <NA>
#>    ConverrtedFilesName SourceFileFolder SourceFileNameR1 SourceFileNameR2
#>                 <char>           <char>           <char>           <char>
#> 1:                <NA>             <NA>             <NA>             <NA>
#> 2:                <NA>             <NA>             <NA>             <NA>
#> 3:                <NA>             <NA>             <NA>             <NA>
#> 4:                <NA>             <NA>             <NA>             <NA>
#> 5:                <NA>             <NA>             <NA>             <NA>
#> 6:                <NA>             <NA>             <NA>             <NA>
#>    DNAnumber ReplicateNumber DNAextractionDate DNAngUsed LinearPCRID
#>       <char>           <int>            <Date>     <num>      <char>
#> 1: PT001-103               3        2016-03-16    23.184        <NA>
#> 2:  PT001-81               2        2016-07-15   181.440        <NA>
#> 3:  PT001-81               1        2016-07-15   181.440        <NA>
#> 4:  PT001-81               3        2016-07-15   181.440        <NA>
#> 5:  PT001-74               1        2016-07-15    23.058        <NA>
#> 6: PT001-107               2        2016-03-16   171.360        <NA>
#>    LinearPCRDate SonicationDate LigationDate 1stExpoPCRID 1stExpoPCRDate
#>           <Date>         <Date>       <Date>       <char>         <Date>
#> 1:          <NA>     2016-11-02   2016-11-02    ET#380.46     2016-11-02
#> 2:          <NA>     2016-11-02   2016-11-02    ET#379.40     2016-11-02
#> 3:          <NA>     2016-11-02   2016-11-02     ET#379.9     2016-11-02
#> 4:          <NA>     2016-11-02   2016-11-02    ET#379.71     2016-11-02
#> 5:          <NA>     2016-11-02   2016-11-02     ET#379.2     2016-11-02
#> 6:          <NA>     2016-11-02   2016-11-02    ET#380.28     2016-11-02
#>    2ndExpoID 2ndExpoDate FusionPrimerPCRID FusionPrimerPCRDate   PoolDate
#>       <char>      <Date>            <char>              <Date>     <Date>
#> 1:      <NA>        <NA>         ET#382.46          2016-11-03 2016-11-07
#> 2:      <NA>        <NA>         ET#381.40          2016-11-03 2016-11-07
#> 3:      <NA>        <NA>          ET#381.9          2016-11-03 2016-11-07
#> 4:      <NA>        <NA>         ET#381.71          2016-11-03 2016-11-07
#> 5:      <NA>        <NA>          ET#381.2          2016-11-03 2016-11-07
#> 6:      <NA>        <NA>         ET#382.28          2016-11-03 2016-11-07
#>    SequencingDate   VCN Genome SequencingRound Genotype TestGroup    MOI
#>            <Date> <num> <char>           <int>   <char>    <char> <char>
#> 1:     2016-11-15  0.30   hg19               1     <NA>      <NA>   <NA>
#> 2:     2016-11-15  0.27   hg19               1     <NA>      <NA>   <NA>
#> 3:     2016-11-15  0.27   hg19               1     <NA>      <NA>   <NA>
#> 4:     2016-11-15  0.27   hg19               1     <NA>      <NA>   <NA>
#> 5:     2016-11-15  0.24   hg19               1     <NA>      <NA>   <NA>
#> 6:     2016-11-15  0.42   hg19               1     <NA>      <NA>   <NA>
#>    Engraftment Transduction  Notes AddedField1 AddedField2 AddedField3
#>          <num>        <num> <char>      <char>      <char>      <char>
#> 1:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 2:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 3:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 4:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 5:          NA           NA   <NA>        <NA>        <NA>        <NA>
#> 6:          NA           NA   <NA>        <NA>        <NA>        <NA>
#>    AddedField4 concatenatePoolIDSeqRun AddedField6_RelativeBloodPercentage
#>         <char>                  <char>                              <char>
#> 1:        <NA>                POOL01-1                                <NA>
#> 2:        <NA>                POOL01-1                                <NA>
#> 3:        <NA>                POOL01-1                                <NA>
#> 4:        <NA>                POOL01-1                                <NA>
#> 5:        <NA>                POOL01-1                                <NA>
#> 6:        <NA>                POOL01-1                                <NA>
#>    AddedField7_PurityTestFeasibility AddedField8_FacsSeparationPurity  Kapa
#>                                <num>                            <num> <num>
#> 1:                                NA                               NA    NA
#> 2:                                NA                               NA    NA
#> 3:                                NA                               NA    NA
#> 4:                                NA                               NA    NA
#> 5:                                NA                               NA    NA
#> 6:                                NA                               NA    NA
#>    ulForPool
#>        <num>
#> 1:        NA
#> 2:        NA
#> 3:        NA
#> 4:        NA
#> 5:        NA
#> 6:        NA
#>                                                 CompleteAmplificationID
#>                                                                  <char>
#> 1: PJ01_POOL01_LTR75LC38_PT001_PT001-103_lenti_GLOBE_PB_1_SLiM_0060_MNC
#> 2:  PJ01_POOL01_LTR53LC32_PT001_PT001-81_lenti_GLOBE_BM_1_SLiM_0180_MNC
#> 3:  PJ01_POOL01_LTR83LC66_PT001_PT001-81_lenti_GLOBE_BM_1_SLiM_0180_MNC
#> 4:  PJ01_POOL01_LTR27LC94_PT001_PT001-81_lenti_GLOBE_BM_1_SLiM_0180_MNC
#> 5:  PJ01_POOL01_LTR69LC52_PT001_PT001-74_lenti_GLOBE_PB_1_SLiM_0180_MNC
#> 6:  PJ01_POOL01_LTR37LC2_PT001_PT001-107_lenti_GLOBE_BM_1_SLiM_0060_MNC
#>                  UniqueID StudyTestID StudyTestGroup MouseID Tigroup Tisource
#>                    <char>      <char>          <int>   <int>  <char>   <char>
#> 1: ID00000000000000007433        <NA>             NA      NA    <NA>     <NA>
#> 2: ID00000000000000007340        <NA>             NA      NA    <NA>     <NA>
#> 3: ID00000000000000007310        <NA>             NA      NA    <NA>     <NA>
#> 4: ID00000000000000007370        <NA>             NA      NA    <NA>     <NA>
#> 5: ID00000000000000007303        <NA>             NA      NA    <NA>     <NA>
#> 6: ID00000000000000007417        <NA>             NA      NA    <NA>     <NA>
#>    PathToFolderProjectID SamplesNameCheck TimepointDays TimepointMonths
#>                   <char>           <char>        <char>          <char>
#> 1:                 /PJ01             <NA>          0060              02
#> 2:                 /PJ01             <NA>          0180              06
#> 3:                 /PJ01             <NA>          0180              06
#> 4:                 /PJ01             <NA>          0180              06
#> 5:                 /PJ01             <NA>          0180              06
#> 6:                 /PJ01             <NA>          0060              02
#>    TimepointYears ng DNA corrected      RUN_NAME PHIX_MAPPING
#>            <char>            <num>        <char>        <int>
#> 1:             01            23.18 PJ01|POOL01-1     43586699
#> 2:             01           181.44 PJ01|POOL01-1     43586699
#> 3:             01           181.44 PJ01|POOL01-1     43586699
#> 4:             01           181.44 PJ01|POOL01-1     43586699
#> 5:             01            23.06 PJ01|POOL01-1     43586699
#> 6:             01           171.36 PJ01|POOL01-1     43586699
#>    PLASMID_MAPPED_BYPOOL BARCODE_MUX LTR_IDENTIFIED TRIMMING_FINAL_LTRLC
#>                    <int>       <int>          <int>                <int>
#> 1:               2256176      645026         645026               630965
#> 2:               2256176      652208         652177               649044
#> 3:               2256176      451519         451512               449669
#> 4:               2256176      426500         426499               425666
#> 5:               2256176       18300          18300                18290
#> 6:               2256176      729327         729327               727219
#>    LV_MAPPED BWA_MAPPED_OVERALL ISS_MAPPED_OVERALL RAW_READS QUALITY_PASSED
#>        <int>              <int>              <int>     <int>          <int>
#> 1:    211757             402477             219452        NA             NA
#> 2:    303300             322086             222646        NA             NA
#> 3:    204810             227275             149385        NA             NA
#> 4:    185752             223915             143283        NA             NA
#> 5:      6962              10487               5907        NA             NA
#> 6:    318653             369117             235640        NA             NA
#>    ISS_MAPPED_PP
#>            <int>
#> 1:            NA
#> 2:            NA
#> 3:            NA
#> 4:            NA
#> 5:            NA
#> 6:            NA
```
