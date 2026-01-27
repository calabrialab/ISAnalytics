# Current dynamic vars specifications getters.

Fetches the look-up tables for different categories of dynamic vars. For
more details, refer to the dedicated vignette
[`vignette("workflow_start", package="ISAnalytics")`](https://calabrialab.github.io/ISAnalytics/dev/articles/workflow_start.md).

- `mandatory_IS_vars` returns the look-up table of variables that are
  used to uniquely identify integration events

&nbsp;

- `annotation_IS_vars()` returns the look-up table of variables that
  contain genomic annotations

&nbsp;

- `association_file_columns()` returns the look-up table of variables
  that contains information on how metadata is structured

&nbsp;

- `iss_stats_specs()` returns the look-up table of variables that
  contains information on the format of pool statistics files produced
  automatically by VISPA2

&nbsp;

- `matrix_file_suffixes()` returns the look-up table of variables that
  contains all default file names for each quantification type and it is
  used by automated import functions

## Usage

``` r
mandatory_IS_vars(include_types = FALSE)

annotation_IS_vars(include_types = FALSE)

association_file_columns(include_types = FALSE)

iss_stats_specs(include_types = FALSE)

matrix_file_suffixes()
```

## Arguments

- include_types:

  If set to `TRUE` returns both the names and the types associated,
  otherwise returns only a character vector of names

## Value

A character vector or a data frame

## See also

Other dynamic vars:
[`inspect_tags()`](https://calabrialab.github.io/ISAnalytics/dev/reference/inspect_tags.md),
[`pcr_id_column()`](https://calabrialab.github.io/ISAnalytics/dev/reference/pcr_id_column.md),
[`reset_mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/reset_mandatory_IS_vars.md),
[`set_mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/set_mandatory_IS_vars.md),
[`set_matrix_file_suffixes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/set_matrix_file_suffixes.md)

## Examples

``` r
# Names only
mandatory_IS_vars()
#> [1] "chr"               "integration_locus" "strand"           

# Names and types
mandatory_IS_vars(TRUE)
#> # A tibble: 3 × 5
#>   names             types transform flag     tag       
#>   <chr>             <chr> <list>    <chr>    <chr>     
#> 1 chr               char  <NULL>    required chromosome
#> 2 integration_locus int   <NULL>    required locus     
#> 3 strand            char  <NULL>    required is_strand 

# Names only
annotation_IS_vars()
#> [1] "GeneName"   "GeneStrand"

# Names and types
annotation_IS_vars(TRUE)
#> # A tibble: 2 × 5
#>   names      types transform flag     tag        
#>   <chr>      <chr> <list>    <chr>    <chr>      
#> 1 GeneName   char  <NULL>    required gene_symbol
#> 2 GeneStrand char  <NULL>    required gene_strand

# Names only
association_file_columns()
#>  [1] "ProjectID"                           "FUSIONID"                           
#>  [3] "PoolID"                              "TagSequence"                        
#>  [5] "SubjectID"                           "VectorType"                         
#>  [7] "VectorID"                            "ExperimentID"                       
#>  [9] "Tissue"                              "TimePoint"                          
#> [11] "DNAFragmentation"                    "PCRMethod"                          
#> [13] "TagIDextended"                       "Keywords"                           
#> [15] "CellMarker"                          "TagID"                              
#> [17] "NGSProvider"                         "NGSTechnology"                      
#> [19] "ConverrtedFilesDir"                  "ConverrtedFilesName"                
#> [21] "SourceFileFolder"                    "SourceFileNameR1"                   
#> [23] "SourceFileNameR2"                    "DNAnumber"                          
#> [25] "ReplicateNumber"                     "DNAextractionDate"                  
#> [27] "DNAngUsed"                           "LinearPCRID"                        
#> [29] "LinearPCRDate"                       "SonicationDate"                     
#> [31] "LigationDate"                        "1stExpoPCRID"                       
#> [33] "1stExpoPCRDate"                      "2ndExpoID"                          
#> [35] "2ndExpoDate"                         "FusionPrimerPCRID"                  
#> [37] "FusionPrimerPCRDate"                 "PoolDate"                           
#> [39] "SequencingDate"                      "VCN"                                
#> [41] "Genome"                              "SequencingRound"                    
#> [43] "Genotype"                            "TestGroup"                          
#> [45] "MOI"                                 "Engraftment"                        
#> [47] "Transduction"                        "Notes"                              
#> [49] "AddedField1"                         "AddedField2"                        
#> [51] "AddedField3"                         "AddedField4"                        
#> [53] "concatenatePoolIDSeqRun"             "AddedField6_RelativeBloodPercentage"
#> [55] "AddedField7_PurityTestFeasibility"   "AddedField8_FacsSeparationPurity"   
#> [57] "Kapa"                                "ulForPool"                          
#> [59] "CompleteAmplificationID"             "UniqueID"                           
#> [61] "StudyTestID"                         "StudyTestGroup"                     
#> [63] "MouseID"                             "Tigroup"                            
#> [65] "Tisource"                            "PathToFolderProjectID"              
#> [67] "SamplesNameCheck"                    "TimepointDays"                      
#> [69] "TimepointMonths"                     "TimepointYears"                     
#> [71] "ng DNA corrected"                   

# Names and types
association_file_columns(TRUE)
#> # A tibble: 71 × 5
#>    names        types transform flag     tag       
#>    <chr>        <chr> <list>    <chr>    <chr>     
#>  1 ProjectID    char  <NULL>    required project_id
#>  2 FUSIONID     char  <NULL>    optional fusion_id 
#>  3 PoolID       char  <NULL>    required pool_id   
#>  4 TagSequence  char  <NULL>    required tag_seq   
#>  5 SubjectID    char  <NULL>    required subject   
#>  6 VectorType   char  <NULL>    optional NA        
#>  7 VectorID     char  <NULL>    required vector_id 
#>  8 ExperimentID char  <NULL>    optional NA        
#>  9 Tissue       char  <NULL>    required tissue    
#> 10 TimePoint    char  <formula> required tp_days   
#> # ℹ 61 more rows

# Names only
iss_stats_specs()
#>  [1] "RUN_NAME"              "POOL"                  "TAG"                  
#>  [4] "RAW_READS"             "QUALITY_PASSED"        "PHIX_MAPPING"         
#>  [7] "PLASMID_MAPPED_BYPOOL" "BARCODE_MUX"           "LTR_IDENTIFIED"       
#> [10] "TRIMMING_FINAL_LTRLC"  "LV_MAPPED"             "BWA_MAPPED_OVERALL"   
#> [13] "ISS_MAPPED_OVERALL"    "ISS_MAPPED_PP"        

# Names and types
iss_stats_specs(TRUE)
#> # A tibble: 14 × 5
#>    names                 types transform flag     tag              
#>    <chr>                 <chr> <list>    <chr>    <chr>            
#>  1 RUN_NAME              char  <NULL>    required NA               
#>  2 POOL                  char  <NULL>    required vispa_concatenate
#>  3 TAG                   char  <formula> required tag_seq          
#>  4 RAW_READS             int   <NULL>    optional NA               
#>  5 QUALITY_PASSED        int   <NULL>    optional NA               
#>  6 PHIX_MAPPING          int   <NULL>    optional NA               
#>  7 PLASMID_MAPPED_BYPOOL int   <NULL>    optional NA               
#>  8 BARCODE_MUX           int   <NULL>    required NA               
#>  9 LTR_IDENTIFIED        int   <NULL>    optional NA               
#> 10 TRIMMING_FINAL_LTRLC  int   <NULL>    optional NA               
#> 11 LV_MAPPED             int   <NULL>    optional NA               
#> 12 BWA_MAPPED_OVERALL    int   <NULL>    optional NA               
#> 13 ISS_MAPPED_OVERALL    int   <NULL>    optional NA               
#> 14 ISS_MAPPED_PP         int   <NULL>    optional NA               

# Names only
matrix_file_suffixes()
#> # A tibble: 10 × 3
#>    quantification   matrix_type   file_suffix                                 
#>    <chr>            <chr>         <chr>                                       
#>  1 seqCount         annotated     seqCount_matrix.no0.annotated.tsv.gz        
#>  2 seqCount         not_annotated seqCount_matrix.tsv.gz                      
#>  3 fragmentEstimate annotated     fragmentEstimate_matrix.no0.annotated.tsv.gz
#>  4 fragmentEstimate not_annotated fragmentEstimate_matrix.tsv.gz              
#>  5 barcodeCount     annotated     barcodeCount_matrix.no0.annotated.tsv.gz    
#>  6 barcodeCount     not_annotated barcodeCount_matrix.tsv.gz                  
#>  7 cellCount        annotated     cellCount_matrix.no0.annotated.tsv.gz       
#>  8 cellCount        not_annotated cellCount_matrix.tsv.gz                     
#>  9 ShsCount         annotated     ShsCount_matrix.no0.annotated.tsv.gz        
#> 10 ShsCount         not_annotated ShsCount_matrix.tsv.gz                      
```
