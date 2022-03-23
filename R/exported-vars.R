#------------------------------------------------------------------------------#
# Exported/Internal variables
#------------------------------------------------------------------------------#

# Internal: default mandatory IS vars and associated column types.
# The combination of these fields defines a unique integration site.
.default_mandatory_IS_vars <- function() {
  tibble::tribble(
    ~ names, ~ types, ~ transform, ~ flag, ~ tag,
    "chr", "char", NULL, "required", "chromosome",
    "integration_locus", "int", NULL, "required", "locus",
    "strand", "char", NULL, "required", "is_strand"
  )
}

# Internal: default genomic annotation IS vars and associated column types.
.default_annotation_IS_vars <- function() {
  tibble::tribble(
    ~ names, ~ types, ~ transform, ~ flag, ~ tag,
    "GeneName", "char", NULL, "required", "gene_symbol",
    "GeneStrand", "char", NULL, "required", "gene_strand"
  )
}

# Internal: default association file columns and types
.default_af_cols <- function() {
  tibble::tribble(
    ~ names, ~ types, ~ transform, ~ flag, ~ tag,
    "ProjectID", "char", NULL, "required", "project_id",
    "FUSIONID", "char", NULL, "optional", NA_character_,
    "PoolID", "char", NULL, "required", "pool_id",
    "TagSequence", "char", NULL, "required", "tag_seq",
    "SubjectID", "char", NULL, "required", "subject",
    "VectorType", "char", NULL, "optional", NA_character_,
    "VectorID", "char", NULL, "required", NA_character_,
    "ExperimentID", "char", NULL, "optional", NA_character_,
    "Tissue", "char", NULL, "required", "tissue",
    "TimePoint", "char", ~ stringr::str_pad(.x, 4, side = "left", pad = "0"),
    "required", "tp_days",
    "DNAFragmentation", "char", NULL, "optional", NA_character_,
    "PCRMethod", "char", NULL, "required", "pcr_method",
    "TagIDextended", "char", NULL, "optional", NA_character_,
    "Keywords","char", NULL, "optional", NA_character_,
    "CellMarker", "char", NULL, "required", "cell_marker",
    "TagID", "char", NULL, "required", NA_character_,
    "NGSProvider", "char", NULL, "optional", NA_character_,
    "NGSTechnology", "char", NULL, "required", "ngs_tech",
    "ConverrtedFilesDir", "char", NULL, "optional", NA_character_,
    "ConverrtedFilesName", "char", NULL, "optional", NA_character_,
    "SourceFileFolder", "char", NULL, "optional", NA_character_,
    "SourceFileNameR1", "char", NULL, "optional", NA_character_,
    "SourceFileNameR2", "char", NULL, "optional", NA_character_,
    "DNAnumber", "char", NULL, "required", "dna_num",
    "ReplicateNumber", "int", NULL, "required", "pcr_replicate",
    "DNAextractionDate", "date", NULL, "optional", NA_character_,
    "DNAngUsed", "numeric", NULL, "required", NA_character_,
    "LinearPCRID", "char", NULL, "optional", NA_character_,
    "LinearPCRDate", "date", NULL, "optional", NA_character_,
    "SonicationDate", "date", NULL, "optional", NA_character_,
    "LigationDate", "date", NULL, "optional", NA_character_,
    "1stExpoPCRID", "char", NULL, "optional", NA_character_,
    "1stExpoPCRDate", "date", NULL, "optional", NA_character_,
    "2ndExpoID", "char", NULL, "optional", NA_character_,
    "2ndExpoDate", "date", NULL, "optional", NA_character_,
    "FusionPrimerPCRID", "char", NULL, "optional", NA_character_,
    "FusionPrimerPCRDate", "date", NULL, "optional", NA_character_,
    "PoolDate", "date", NULL, "optional", NA_character_,
    "SequencingDate", "date", NULL, "required", NA_character_,
    "VCN", "numeric", NULL, "required", "vcn",
    "Genome", "char", NULL, "required", "genome",
    "SequencingRound", "int", NULL, "optional", NA_character_,
    "Genotype", "char", NULL, "optional", NA_character_,
    "TestGroup", "char", NULL, "optional", NA_character_,
    "MOI", "char", NULL, "optional", NA_character_,
    "Engraftment", "numeric", NULL, "optional", NA_character_,
    "Transduction", "numeric", NULL, "optional", NA_character_,
    "Notes", "char", NULL, "optional", NA_character_,
    "AddedField1", "char", NULL, "optional", NA_character_,
    "AddedField2", "char", NULL, "optional", NA_character_,
    "AddedField3", "char", NULL, "optional", NA_character_,
    "AddedField4", "char", NULL, "optional", NA_character_,
    "concatenatePoolIDSeqRun", "char", NULL,"required", "vispa_concatenate",
    "AddedField6_RelativeBloodPercentage", "char", NULL, "optional",
    NA_character_,
    "AddedField7_PurityTestFeasibility", "char", NULL, "optional",
    NA_character_,
    "AddedField8_FacsSeparationPurity", "char", NULL, "optional", NA_character_,
    "Kapa", "numeric", NULL, "required", NA_character_,
    "ulForPool", "numeric", NULL, "required", NA_character_,
    "CompleteAmplificationID", "char", NULL, "required", "pcr_repl_id",
    "UniqueID", "char", NULL, "required", NA_character_,
    "StudyTestID", "char", NULL, "optional", NA_character_,
    "StudyTestGroup", "char", NULL, "optional", NA_character_,
    "MouseID", "char", NULL, "optional", NA_character_,
    "Tigroup", "char", NULL, "optional", NA_character_,
    "Tisource", "char", NULL, "optional", NA_character_,
    "PathToFolderProjectID", "char", NULL, "required", "proj_folder",
    "SamplesNameCheck", "char", NULL, "optional", NA_character_,
    "TimepointDays", "char", NULL, "optional", NA_character_,
    "TimepointMonths", "char", NULL, "optional", NA_character_,
    "TimepointYears", "char", NULL, "optional", NA_character_,
    "ng DNA corrected", "numeric", NULL, "optional", NA_character_
  )
}

# Internal: default columns and types of vispa2 stats cols
.default_iss_stats_specs <- function() {
  tibble::tribble(
    ~ names, ~ types, ~ transform, ~ flag, ~ tag,
    "RUN_NAME", "char", NULL, "required", NA_character_,
    "POOL", "char", NULL, "required", "vispa_concatenate",
    "TAG", "char", ~ stringr::str_replace_all(.x, pattern = "\\.",
                                              replacement = ""), "required",
    "tag_seq",
    "RAW_READS", "int", NULL, "optional", NA_character_,
    "QUALITY_PASSED", "int", NULL, "optional", NA_character_,
    "PHIX_MAPPING", "int", NULL, "optional", NA_character_,
    "PLASMID_MAPPED_BYPOOL", "int", NULL, "optional", NA_character_,
    "BARCODE_MUX", "int", NULL, "required", NA_character_,
    "LTR_IDENTIFIED", "int", NULL, "optional", NA_character_,
    "TRIMMING_FINAL_LTRLC", "int", NULL, "optional", NA_character_,
    "LV_MAPPED", "int", NULL, "optional", NA_character_,
    "BWA_MAPPED_OVERALL", "int", NULL, "optional", NA_character_,
    "ISS_MAPPED_OVERALL", "int", NULL, "optional", NA_character_,
    "ISS_MAPPED_PP", "int", NULL, "optional", NA_character_
  )
}

# Mappings between input format and formats requested by external parsing
# functions
.types_mapping <- function() {
  tibble::tribble(
    ~ types, ~ mapping, ~ fread,
    "char", "c", "character",
    "int", "i", "integer",
    "logi", "l", "logical",
    "numeric", "d", "numeric",
    "factor", "f", "factor",
    "date", "c", "charcter",
    "ymd", "c", "character",
    "ydm", "c", "character",
    "mdy", "c", "character",
    "myd", "c", "character",
    "dmy", "c", "character",
    "yq", "c", "character",
    "ym", "c", "character",
    "my", "c", "character",
    "ymd_hms", "c", "character",
    "ymd_hm", "c", "character",
    "ymd_h", "c", "character",
    "dmy_hms", "c", "character",
    "dmy_hm", "c", "character",
    "dmy_h", "c", "character",
    "mdy_hms", "c", "character",
    "mdy_hm", "c", "character",
    "mdy_h", "c", "character",
    "ydm_hms", "c", "character",
    "ydm_hm", "c", "character",
    "ydm_h", "c", "character"
  )
}

# Internal: associates column types with column names for a more precise
# import
.mandatory_IS_types <- function(mode) {
  specs <- mandatory_IS_vars(include_types = TRUE)
  specs_mappings <- specs %>%
    dplyr::left_join(.types_mapping(), by = "types")
  if (mode == "fread") {
    specs_mappings <- specs_mappings %>%
      dplyr::select(.data$names, .data$fread) %>%
      dplyr::group_by(.data$fread)
    types <- specs_mappings %>%
      dplyr::group_keys() %>%
      dplyr::pull(.data$fread)
    specs_mappings <- specs_mappings %>%
      dplyr::group_split(.keep = FALSE)
    names(specs_mappings) <- types
    types <- purrr::map(specs_mappings, ~ .x$names)
    return(types)
  }
  types <- as.list(setNames(specs_mappings$mapping, specs_mappings$names))
  return(types)
}

# Internal: associates column types with column names for a more precise
# import
.annotation_IS_types <- function(mode) {
  specs <- annotation_IS_vars(include_types = TRUE)
  specs_mappings <- specs %>%
    dplyr::left_join(.types_mapping(), by = "types")
  if (mode == "fread") {
    specs_mappings <- specs_mappings %>%
      dplyr::select(.data$names, .data$fread) %>%
      dplyr::group_by(.data$fread)
    types <- specs_mappings %>%
      dplyr::group_keys() %>%
      dplyr::pull(.data$fread)
    specs_mappings <- specs_mappings %>%
      dplyr::group_split(.keep = FALSE)
    names(specs_mappings) <- types
    types <- purrr::map(specs_mappings, ~ .x$names)
    return(types)
  }
  types <- as.list(setNames(specs_mappings$mapping, specs_mappings$names))
  return(types)
}

# Internal: associates column types with column names for a more precise
# import
.af_col_types <- function(mode) {
  specs <- association_file_columns(include_types = TRUE)
  specs_mappings <- specs %>%
    dplyr::left_join(.types_mapping(), by = "types")
  if (mode == "fread") {
    specs_mappings <- specs_mappings %>%
      dplyr::select(.data$names, .data$fread) %>%
      dplyr::group_by(.data$fread)
    types <- specs_mappings %>%
      dplyr::group_keys() %>%
      dplyr::pull(.data$fread)
    specs_mappings <- specs_mappings %>%
      dplyr::group_split(.keep = FALSE)
    names(specs_mappings) <- types
    types <- purrr::map(specs_mappings, ~ .x$names)
    return(types)
  }
  types <- as.list(setNames(specs_mappings$mapping, specs_mappings$names))
}


# Internal: used for file system alignment in import_association_file,
# gives the names of the columns that respectively contain:
# - the absolute path on disk of the project
# - the path to the quantification folder
# - the path to the iss folder
.path_cols_names <- function() {
    list(project = "Path", quant = "Path_quant", iss = "Path_iss")
}


.default_matrix_suffixes <- function() {

}

.matrix_annotated_suffixes <- function() {
    c(".no0.annotated")
}

.matrix_not_annotated_suffixes <- function() {
    c("")
}

#' Names of the columns of the association file to consider for
#' Vispa2 launch.
#'
#' Selection of column names from the association file to be considered for
#' Vispa2 launch.
#' NOTE: the `TagID` column appears only once but needs to be
#' repeated twice for generating the launch file. Use the appropriate
#' function to generate the file automatically.
#'
#' @return A character vector
#' @export
#'
#' @examples
#' reduced_AF_columns()
reduced_AF_columns <- function() {
    c(
        "TagID", "Tissue", "SubjectID", "TimePoint", "FUSIONID",
        "CompleteAmplificationID", "CellMarker", "ProjectID", "VectorID",
        "PoolID"
    )
}

# Names of the columns of iss stats considered for aggregation
# USED IN:
# - .join_and_aggregate
.agg_iss_cols <- function() {
    c(
        "BARCODE_MUX", "TRIMMING_FINAL_LTRLC",
        "LV_MAPPED",
        "BWA_MAPPED_OVERALL",
        "ISS_MAPPED_PP"
    )
}

.compressed_formats <- function() {
    c("gz", "bz2", "xz", "zip")
}

.supported_fread_compression_formats <- function() {
    c("gz", "bz2")
}


flag_logics <- function() {
    c("AND", "OR", "XOR", "NAND", "NOR", "XNOR")
}

### Articles links
.lentiviral_CIS_paper <- function() {
    paste0(
        "https://ashpublications.org/blood/article/117/20/5332/21206/",
        "Lentiviral-vector-common-integration-sites-in"
    )
}

.vispa2_paper_link <- function() {
    paste0("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5702242/")
}

.vispa2_paper_title <- function() {
    paste0(
        "VISPA2:A Scalable Pipeline for High-Throughput ",
        "Identification and Annotation of Vector Integration Sites"
    )
}


#' Required columns for refGene file.
#'
#' @return Character vector of column names
#' @export
#'
#' @examples
#' refGene_table_cols()
refGene_table_cols <- function() {
    c(
        "name2", "chrom", "strand", "min_txStart", "max_txEnd",
        "minmax_TxLen", "average_TxLen", "name", "min_cdsStart",
        "max_cdsEnd", "minmax_CdsLen", "average_CdsLen"
    )
}


available_column_tags <- function() {
    list(
        critical = list(af = c(
            "project_id", "pool_id", "tag_seq", "subject", "tissue",
            "cell_marker", "pcr_replicate", "vispa_concatenate",
            "pcr_repl_id", "proj_folder"
        ),
        matrix = c("chromosome", "locus", "is_strand", "gene_symbol"),
        stats = c("vispa_concatenate", "tag_seq")),
        optional = list(af = c(
            "tp_days", "pcr_method", "ngs_tech", "dna_num",
            "vcn", "genome"
        ),
        matrix = c("gene_strand"),
        stats = c())
    )
}
