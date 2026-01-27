# Import Vispa2 stats given the aligned association file.

**\[stable\]** Imports all the Vispa2 stats files for each pool provided
the association file has been aligned with the file system (see
[`import_association_file`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_association_file.md)).

## Usage

``` r
import_Vispa2_stats(
  association_file,
  file_prefixes = default_iss_file_prefixes(),
  join_with_af = TRUE,
  pool_col = "concatenatePoolIDSeqRun",
  report_path = default_report_path()
)
```

## Arguments

- association_file:

  The file system aligned association file (contains columns with
  absolute paths to the 'iss' folder)

- file_prefixes:

  A character vector with known file prefixes to match on file names.
  NOTE: the elements represent regular expressions. For defaults see
  [default_iss_file_prefixes](https://calabrialab.github.io/ISAnalytics/dev/reference/default_iss_file_prefixes.md).

- join_with_af:

  Logical, if `TRUE` the imported stats files will be merged with the
  association file, if `FALSE` a single data frame holding only the
  stats will be returned.

- pool_col:

  A single string. What is the name of the pool column used in the
  Vispa2 run? This will be used as a key to perform a join operation
  with the stats files `POOL` column.

- report_path:

  The path where the report file should be saved. Can be a folder or
  `NULL` if no report should be produced. Defaults to
  `{user_home}/ISAnalytics_reports`.

## Value

A data frame

## Required tags

The function will explicitly check for the presence of these tags:

- project_id

- tag_seq

- vispa_concatenate

- pcr_repl_id

## See also

Other Import functions:
[`import_association_file()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_association_file.md),
[`import_parallel_Vispa2Matrices()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_parallel_Vispa2Matrices.md),
[`import_single_Vispa2Matrix()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_single_Vispa2Matrix.md)

## Examples

``` r
fs_path <- generate_default_folder_structure(type = "correct")
af <- import_association_file(fs_path$af,
    root = fs_path$root,
    import_iss = FALSE,
    report_path = NULL
)
stats_files <- import_Vispa2_stats(af,
    join_with_af = FALSE,
    report_path = NULL
)
head(stats_files)
#> # A tibble: 6 × 14
#>   POOL     TAG       RUN_NAME     PHIX_MAPPING PLASMID_MAPPED_BYPOOL BARCODE_MUX
#>   <chr>    <chr>     <chr>               <dbl>                 <dbl>       <dbl>
#> 1 POOL01-1 LTR75LC38 PJ01|POOL01…     43586699               2256176      645026
#> 2 POOL01-1 LTR53LC32 PJ01|POOL01…     43586699               2256176      652208
#> 3 POOL01-1 LTR83LC66 PJ01|POOL01…     43586699               2256176      451519
#> 4 POOL01-1 LTR27LC94 PJ01|POOL01…     43586699               2256176      426500
#> 5 POOL01-1 LTR69LC52 PJ01|POOL01…     43586699               2256176       18300
#> 6 POOL01-1 LTR37LC2  PJ01|POOL01…     43586699               2256176      729327
#> # ℹ 8 more variables: LTR_IDENTIFIED <dbl>, TRIMMING_FINAL_LTRLC <dbl>,
#> #   LV_MAPPED <dbl>, BWA_MAPPED_OVERALL <dbl>, ISS_MAPPED_OVERALL <dbl>,
#> #   RAW_READS <lgl>, QUALITY_PASSED <lgl>, ISS_MAPPED_PP <lgl>
```
