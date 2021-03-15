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
        "PathToFolderProjectID"
    )
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
