# Import the association file from disk

**\[stable\]** Imports the association file and optionally performs a
check on the file system starting from the root to assess the alignment
between the two.

## Usage

``` r
import_association_file(
  path,
  root = NULL,
  dates_format = "ymd",
  separator = "\t",
  filter_for = NULL,
  import_iss = FALSE,
  convert_tp = TRUE,
  report_path = default_report_path(),
  transformations = default_af_transform(convert_tp),
  tp_padding = lifecycle::deprecated(),
  ...
)
```

## Arguments

- path:

  The path on disk to the association file.

- root:

  The path on disk of the root folder of VISPA2 output or `NULL`. See
  details.

- dates_format:

  A single string indicating how dates should be parsed. Must be a value
  in:
  [`date_formats()`](https://calabrialab.github.io/ISAnalytics/dev/reference/date_formats.md)

- separator:

  The column separator used in the file

- filter_for:

  A named list where names represent column names that must be filtered.
  For example: `list(ProjectID = c("PROJECT1", "PROJECT2))` will filter
  the association file so that it contains only those rows for which the
  value of the column "ProjectID" is one of the specified values. If
  multiple columns are present in the list all filtering conditions are
  applied as a logical AND.

- import_iss:

  Import VISPA2 pool stats and merge them with the association file?
  Logical value

- convert_tp:

  Should be time points be converted into months and years? Logical
  value

- report_path:

  The path where the report file should be saved. Can be a folder or
  `NULL` if no report should be produced. Defaults to
  `{user_home}/ISAnalytics_reports`.

- transformations:

  Either `NULL` or a named list of purrr-style lambdas where names are
  column names the function should be applied to.

- tp_padding:

  **\[deprecated\]** Deprecated. Use `transformations` instead.

- ...:

  Additional arguments to pass to
  [`import_Vispa2_stats`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_Vispa2_stats.md)

## Value

The data frame containing metadata

## Details

### Transformations

Lambdas provided in input in the `transformations` argument, must be
transformations, aka functions that take in input a vector and return a
vector of the same length as the input.

If the transformation list contains column names that are not present in
the data frame, they are simply ignored.

### File system alignment

If the `root` argument is set to `NULL` no file system alignment is
performed. This allows to import the basic file but it won't be possible
to perform automated matrix and stats import. For more details see the
"How to use import functions" vignette:
[`vignette("workflow_start", package = "ISAnalytics")`](https://calabrialab.github.io/ISAnalytics/dev/articles/workflow_start.md)

### Time point conversion

The time point conversion is based on the following logic, given `TPD`
is the column containing the time point expressed in days and `TPM` and
`TPY` are respectively the time points expressed as month and years

- If `TPD` is `NA` –\> `NA` (for both months and years)

- `TPM` = 0, `TPY` = 0 if and only if `TPD` = 0

For conversion in months:

- `TPM` = ceiling(`TPD`/30) if `TPD` \< 30 otherwise `TPM` =
  round(`TPD`/30)

For conversion in years:

- `TPY` = ceiling(`TPD`/360)

## Required tags

The function will explicitly check for the presence of these tags:

- project_id

- pool_id

- tag_seq

- subject

- tissue

- tp_days

- cell_marker

- pcr_replicate

- vispa_concatenate

- pcr_repl_id

- proj_folder

The function will use all the available specifications contained in
`association_file_columns(TRUE)` to read and parse the file. If the
specifications contain columns with a type `"date"`, the function will
parse the generic date with the format in the `dates_format` argument.

## See also

[`transform_columns`](https://calabrialab.github.io/ISAnalytics/dev/reference/transform_columns.md)

[`date_formats`](https://calabrialab.github.io/ISAnalytics/dev/reference/date_formats.md)

Other Import functions:
[`import_Vispa2_stats()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_Vispa2_stats.md),
[`import_parallel_Vispa2Matrices()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_parallel_Vispa2Matrices.md),
[`import_single_Vispa2Matrix()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_single_Vispa2Matrix.md)

## Examples

``` r
fs_path <- generate_default_folder_structure(type = "correct")
af <- import_association_file(fs_path$af,
    root = fs_path$root,
    report_path = NULL
)
head(af)
#> # A tibble: 6 × 74
#>   ProjectID FUSIONID  PoolID TagSequence SubjectID VectorType VectorID
#>   <chr>     <chr>     <chr>  <chr>       <chr>     <chr>      <chr>   
#> 1 PJ01      ET#382.46 POOL01 LTR75LC38   PT001     lenti      GLOBE   
#> 2 PJ01      ET#381.40 POOL01 LTR53LC32   PT001     lenti      GLOBE   
#> 3 PJ01      ET#381.9  POOL01 LTR83LC66   PT001     lenti      GLOBE   
#> 4 PJ01      ET#381.71 POOL01 LTR27LC94   PT001     lenti      GLOBE   
#> 5 PJ01      ET#381.2  POOL01 LTR69LC52   PT001     lenti      GLOBE   
#> 6 PJ01      ET#382.28 POOL01 LTR37LC2    PT001     lenti      GLOBE   
#> # ℹ 67 more variables: ExperimentID <chr>, Tissue <chr>, TimePoint <chr>,
#> #   DNAFragmentation <chr>, PCRMethod <chr>, TagIDextended <chr>,
#> #   Keywords <chr>, CellMarker <chr>, TagID <chr>, NGSProvider <chr>,
#> #   NGSTechnology <chr>, ConverrtedFilesDir <chr>, ConverrtedFilesName <chr>,
#> #   SourceFileFolder <chr>, SourceFileNameR1 <chr>, SourceFileNameR2 <chr>,
#> #   DNAnumber <chr>, ReplicateNumber <int>, DNAextractionDate <date>,
#> #   DNAngUsed <dbl>, LinearPCRID <chr>, LinearPCRDate <date>, …
```
