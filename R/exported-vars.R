#------------------------------------------------------------------------------#
# Exported/Internal variables
#------------------------------------------------------------------------------#

#' Names of mandatory variables for an integration matrix.
#'
#' Contains the names of the columns that need to be present in order for a
#' tibble to be considered an integration matrix.
#'
#' @return A character vector
#' @export
#'
#' @examples
#' mandatory_IS_vars()
mandatory_IS_vars <- function() {
    c("chr", "integration_locus", "strand")
}

# Internal: associates column types with column names for a more precise
# import
.mandatory_IS_types <- function(mode) {
    if (mode == "fread") {
        return(list(
            character = c("chr", "strand"),
            integer = "integration_locus"
        ))
    } else {
        return(
            list(
                chr = "c",
                integration_locus = "i",
                strand = "c"
            )
        )
    }
}

#' Names of the annotation variables for an integration matrix.
#'
#' Contains the names of the columns that are present if the integration matrix
#' is annotated.
#'
#' @return A character vector
#' @export
#'
#' @examples
#' annotation_IS_vars()
annotation_IS_vars <- function() {
    c("GeneName", "GeneStrand")
}

# Internal: associates column types with column names for a more precise
# import
.annotation_IS_types <- function(mode) {
    if (mode == "fread") {
        return(list(character = c("GeneName", "GeneStrand")))
    } else {
        return(list(
            GeneName = "c",
            GeneStrand = "c"
        ))
    }
}

#' Names of the columns in the association file.
#'
#' All the names of the columns present in the association file.
#'
#' @return A character vector
#' @export
#'
#' @examples
#' association_file_columns()
association_file_columns <- function() {
    c(
        "ProjectID", "FUSIONID", "PoolID", "TagSequence", "SubjectID",
        "VectorType", "VectorID", "ExperimentID", "Tissue", "TimePoint",
        "DNAFragmentation", "PCRMethod", "TagIDextended", "Keywords",
        "CellMarker",
        "TagID", "NGSProvider", "NGSTechnology", "ConverrtedFilesDir",
        "ConverrtedFilesName", "SourceFileFolder", "SourceFileNameR1",
        "SourceFileNameR2", "DNAnumber", "ReplicateNumber", "DNAextractionDate",
        "DNAngUsed", "LinearPCRID", "LinearPCRDate", "SonicationDate",
        "LigationDate", "1stExpoPCRID", "1stExpoPCRDate", "2ndExpoID",
        "2ndExpoDate", "FusionPrimerPCRID", "FusionPrimerPCRDate",
        "PoolDate", "SequencingDate", "VCN", "Genome", "SequencingRound",
        "Genotype", "TestGroup", "MOI", "Engraftment", "Transduction", "Notes",
        "AddedField1", "AddedField2", "AddedField3", "AddedField4",
        "concatenatePoolIDSeqRun", "AddedField6_RelativeBloodPercentage",
        "AddedField7_PurityTestFeasibility", "AddedField8_FacsSeparationPurity",
        "Kapa", "ulForPool", "CompleteAmplificationID", "UniqueID",
        "StudyTestID",
        "StudyTestGroup", "MouseID", "Tigroup", "Tisource",
        "PathToFolderProjectID",
        "SamplesNameCheck",
        "TimepointDays", "TimepointMonths",
        "TimepointYears", "ng DNA corrected"
    )
}

# Internal: associates column types with column names for a more precise
# import
.af_col_types <- function(mode) {
    if (mode == "fread") {
        types <- list(
            character = c(
                "ProjectID", "FUSIONID", "PoolID", "TagSequence",
                "SubjectID", "VectorType", "VectorID", "ExperimentID",
                "Tissue", "TimePoint", "DNAFragmentation",
                "PCRMethod", "TagIDextended", "Keywords",
                "CellMarker", "TagID", "NGSProvider", "NGSTechnology",
                "ConverrtedFilesDir", "ConverrtedFilesName",
                "SourceFileFolder", "SourceFileNameR1",
                "SourceFileNameR2", "DNAnumber", "LinearPCRID",
                "1stExpoPCRID", "2ndExpoID", "FusionPrimerPCRID",
                "Genome", "Genotype", "Notes", "AddedField1",
                "AddedField2", "AddedField3", "AddedField4",
                "concatenatePoolIDSeqRun", "CompleteAmplificationID",
                "UniqueID", "StudyTestID", "Tigroup", "Tisource",
                "PathToFolderProjectID", "SamplesNameCheck",
                "DNAextractionDate", "LinearPCRDate",
                "SonicationDate", "LigationDate",
                "FusionPrimerPCRDate", "PoolDate", "SequencingDate",
                "MOI", "AddedField6_RelativeBloodPercentage",
                "TestGroup"
            ),
            double = c(
                "DNAngUsed", "VCN", "Engraftment", "Transduction",
                "AddedField7_PurityTestFeasibility",
                "AddedField8_FacsSeparationPurity", "Kapa",
                "ulForPool", "TimepointMonths", "TimepointYears",
                "ng DNA corrected"
            ),
            integer = c(
                "ReplicateNumber", "SequencingRound",
                "StudyTestGroup", "MouseID", "TimepointDays"
            )
        )
        return(types)
    }
    if (mode == "readr") {
        types <- list(
            ProjectID = "c", FUSIONID = "c", PoolID = "c", TagSequence = "c",
            SubjectID = "c", VectorType = "c", VectorID = "c",
            ExperimentID = "c", Tissue = "c", TimePoint = "c",
            DNAFragmentation = "c", PCRMethod = "c", TagIDextended = "c",
            Keywords = "c", CellMarker = "c", TagID = "c",
            NGSProvider = "c", NGSTechnology = "c",
            ConverrtedFilesDir = "c", ConverrtedFilesName = "c",
            SourceFileFolder = "c", SourceFileNameR1 = "c",
            SourceFileNameR2 = "c", DNAnumber = "c", LinearPCRID = "c",
            `1stExpoPCRID` = "c", `2ndExpoID` = "c",
            FusionPrimerPCRID = "c", Genome = "c", Genotype = "c",
            Notes = "c", AddedField1 = "c",
            AddedField2 = "c", AddedField3 = "c", AddedField4 = "c",
            concatenatePoolIDSeqRun = "c", CompleteAmplificationID = "c",
            UniqueID = "c", StudyTestID = "c", Tigroup = "c", Tisource = "c",
            PathToFolderProjectID = "c", SamplesNameCheck = "c",
            DNAextractionDate = "c",
            LinearPCRDate = "c",
            SonicationDate = "c",
            LigationDate = "c",
            FusionPrimerPCRDate = "c",
            PoolDate = "c",
            SequencingDate = "c",
            MOI = "c", AddedField6_RelativeBloodPercentage = "c",
            DNAngUsed = "d", VCN = "d", Engraftment = "d", Transduction = "d",
            AddedField7_PurityTestFeasibility = "d",
            AddedField8_FacsSeparationPurity = "d", Kapa = "d",
            ulForPool = "d", TimepointMonths = "d", TimepointYears = "d",
            `ng DNA corrected` = "d",
            ReplicateNumber = "i", SequencingRound = "i", TestGroup = "c",
            MouseID = "i", TimepointDays = "i",
            `1stExpoPCRDate` = "c",
            `2ndExpoDate` = "c",
            StudyTestGroup = "i"
        )
        return(types)
    }
}

# Internal: used for file system alignment in import_association_file,
# gives the names of the columns that respectively contain:
# - the absolute path on disk of the project
# - the path to the quantification folder
# - the path to the iss folder
.path_cols_names <- function() {
    list(project = "Path", quant = "Path_quant", iss = "Path_iss")
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
