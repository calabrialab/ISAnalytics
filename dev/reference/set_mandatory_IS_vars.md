# Define custom dynamic vars.

This set of function allows users to specify custom look-up tables for
dynamic variables. For more details, refer to the dedicated vignette
[`vignette("workflow_start", package="ISAnalytics")`](https://calabrialab.github.io/ISAnalytics/dev/articles/workflow_start.md).

- `set_mandatory_IS_vars()` sets the look-up table for mandatory IS
  vars.

&nbsp;

- `set_annotation_IS_vars()` sets the look-up table for genomic
  annotation IS vars.

&nbsp;

- `set_af_columns_def()` sets the look-up table for association file
  columns vars

&nbsp;

- `set_iss_stats_specs()` sets the look-up table for VISPA2 pool
  statistics vars

## Usage

``` r
set_mandatory_IS_vars(specs)

set_annotation_IS_vars(specs)

set_af_columns_def(specs)

set_iss_stats_specs(specs)
```

## Arguments

- specs:

  Either a named vector or a data frame with specific format. See
  details.

## Value

`NULL`

## Details

The user can supply specifications in the form of a named vector or a
data frame.

### Named vector

When using a named vector, names should be the names of the columns,
values should be the type associated with each column in the form of a
string. The vector gets automatically converted into a data frame with
the right format (default values for the columns `transform` and `flag`
are `NULL` and `required` respectively). Use of this method is however
discouraged: data frame inputs are preferred since they offer more
control.

### Look-up table structure

The look-up table for dynamic vars should always follow this structure:

|                        |          |                      |          |         |
|------------------------|----------|----------------------|----------|---------|
| names                  | types    | transform            | flag     | tag     |
| `<name of the column>` | `<type>` | `<a lambda or NULL>` | `<flag>` | `<tag>` |

where

- `names` contains the name of the column as a character

- `types` contains the type of the column. Type should be expressed as a
  string and should be in one of the allowed types

- `char` for character (strings)

- `int` for integers

- `logi` for logical values (TRUE / FALSE)

- `numeric` for numeric values

- `factor` for factors

- `date` for generic date format - note that functions that need to read
  and parse files will try to guess the format and parsing may fail

- One of the accepted date/datetime formats by `lubridate`, you can use
  [`ISAnalytics::date_formats()`](https://calabrialab.github.io/ISAnalytics/dev/reference/date_formats.md)
  to view the accepted formats

- `transform`: a purrr-style lambda that is applied immediately after
  importing. This is useful to operate simple transformations like
  removing unwanted characters or rounding to a certain precision.
  Please note that these lambdas need to be functions that accept a
  vector as input and only operate a transformation, aka they output a
  vector of the same length as the input. For more complicated
  applications that may require the value of other columns, appropriate
  functions should be manually applied post-import.

- `flag`: as of now, it should be set either to `required` or
  `optional` - some functions internally check for only required tags
  presence and if those are missing from inputs they fail, signaling
  failure to the user

- `tag`: a specific tag expressed as a string

### Column types:

Type should be expressed as a string and should be in one of the allowed
types

- `char` for character (strings)

- `int` for integers

- `logi` for logical values (TRUE / FALSE)

- `numeric` for numeric values

- `factor` for factors

- `date` for generic date format - note that functions that need to read
  and parse files will try to guess the format and parsing may fail

- One of the accepted date/datetime formats by `lubridate`, you can use
  [`ISAnalytics::date_formats()`](https://calabrialab.github.io/ISAnalytics/dev/reference/date_formats.md)
  to view the accepted formats

## See also

Other dynamic vars:
[`inspect_tags()`](https://calabrialab.github.io/ISAnalytics/dev/reference/inspect_tags.md),
[`mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md),
[`pcr_id_column()`](https://calabrialab.github.io/ISAnalytics/dev/reference/pcr_id_column.md),
[`reset_mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/reset_mandatory_IS_vars.md),
[`set_matrix_file_suffixes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/set_matrix_file_suffixes.md)

## Examples

``` r
tmp_mand_vars <- tibble::tribble(
    ~names, ~types, ~transform, ~flag, ~tag,
    "chrom", "char", ~ stringr::str_replace_all(.x, "chr", ""), "required",
    "chromosome",
    "position", "int", NULL, "required", "locus",
    "strand", "char", NULL, "required", "is_strand",
    "gap", "int", NULL, "required", NA_character_,
    "junction", "int", NULL, "required", NA_character_
)
set_mandatory_IS_vars(tmp_mand_vars)
#> Mandatory IS vars successfully changed
print(mandatory_IS_vars(TRUE))
#> # A tibble: 5 × 5
#>   names    types transform flag     tag       
#>   <chr>    <chr> <list>    <chr>    <chr>     
#> 1 chrom    char  <formula> required chromosome
#> 2 position int   <NULL>    required locus     
#> 3 strand   char  <NULL>    required is_strand 
#> 4 gap      int   <NULL>    required NA        
#> 5 junction int   <NULL>    required NA        
reset_mandatory_IS_vars()
#> Mandatory IS vars reset to default

tmp_annot_vars <- tibble::tribble(
    ~names, ~types, ~transform, ~flag, ~tag,
    "gene", "char", NULL, "required",
    "gene_symbol",
    "gene_strand", "char", NULL, "required", "gene_strand"
)
print(annotation_IS_vars(TRUE))
#> # A tibble: 2 × 5
#>   names      types transform flag     tag        
#>   <chr>      <chr> <list>    <chr>    <chr>      
#> 1 GeneName   char  <NULL>    required gene_symbol
#> 2 GeneStrand char  <NULL>    required gene_strand
reset_annotation_IS_vars()
#> Annotation IS vars reset to default

temp_af_cols <- tibble::tribble(
    ~names, ~types, ~transform, ~flag, ~tag,
    "project", "char", NULL, "required",
    "project_id",
    "pcr_id", "char", NULL, "required", "pcr_repl_id",
    "subject", "char", NULL, "required", "subject"
)
set_af_columns_def(temp_af_cols)
#> Warning: Warning: important tags missing
#> ℹ Some tags are required for proper execution of some functions. If these tags are not provided, execution of dependent functions might fail. Review your inputs carefully.
#> ℹ Missing tags: pool_id, fusion_id, tag_seq, vector_id, tissue, tp_days, cell_marker, tag_id, pcr_replicate, vispa_concatenate, proj_folder
#> ℹ To see where these are involved type `inspect_tags(c('pool_id','fusion_id','tag_seq','vector_id','tissue','tp_days','cell_marker','tag_id','pcr_replicate','vispa_concatenate','proj_folder'))`
#> Association file columns specs successfully changed
print(association_file_columns(TRUE))
#> # A tibble: 3 × 5
#>   names   types transform flag     tag        
#>   <chr>   <chr> <list>    <chr>    <chr>      
#> 1 project char  <NULL>    required project_id 
#> 2 pcr_id  char  <NULL>    required pcr_repl_id
#> 3 subject char  <NULL>    required subject    
reset_af_columns_def()
#> Association file columns specs reset to default

tmp_iss_vars <- tibble::tribble(
    ~names, ~types, ~transform, ~flag, ~tag,
    "pool", "char", NULL, "required",
    "vispa_concatenate",
    "tag", "char", NULL, "required", "tag_seq",
    "barcode", "int", NULL, "required", NA_character_
)
set_iss_stats_specs(tmp_iss_vars)
#> ISS stats specs successfully changed
iss_stats_specs(TRUE)
#> # A tibble: 3 × 5
#>   names   types transform flag     tag              
#>   <chr>   <chr> <list>    <chr>    <chr>            
#> 1 pool    char  <NULL>    required vispa_concatenate
#> 2 tag     char  <NULL>    required tag_seq          
#> 3 barcode int   <NULL>    required NA               
reset_iss_stats_specs()
#> ISS stats specs reset to default
```
