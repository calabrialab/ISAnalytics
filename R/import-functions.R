#------------------------------------------------------------------------------#
# Importing functions
#------------------------------------------------------------------------------#

#' Imports a single integration matrix from file.
#'
#' @description This function allows to read and import an integration matrix
#' produced as the output of Vispa2 pipeline and converts it to a tidy
#' ISADataFrame.
#'
#' @param path the path to the file on disk
#'
#' @return a tidy ISADataFrame
#' @importFrom tidyr separate
#' @importFrom utils read.csv
#' @export
#'
#' @examples
#' path_to_file <- system.file("extdata", "ex_annotated_ISMatrix.tsv.xz",
#' package = "ISAnalytics")
#' isa_dataframe <- import_single_Vispa2Matrix(path_to_file)
import_single_Vispa2Matrix <- function(path) {
  stopifnot(!missing(path) & is.character(path))
  if (!file.exists(path)) {
    stop(paste("File not found at", path))
  }
  df <- read.csv(path, sep = "\t", header = TRUE, fill = TRUE,
                 check.names = FALSE, stringsAsFactors = FALSE)

  df_type <- .auto_detect_type(df)
  switch(df_type,
         "OLD" = {
           df <- df %>% tidyr::separate(col = .data$IS_genomicID,
                                        into = c("chr", "integration_locus",
                                                 "strand"),
                                        sep = "_", remove = TRUE,
                                        convert = TRUE)
           isadf <- ISADataFrame(df)
         },
         "NEW_ANNOTATED" = {
           isadf <- ISADataFrame(df, metadata = c("GeneName", "GeneStrand"))
         },
         "NEW_NOTANN" = {
           isadf <- ISADataFrame(df)
         },
         "MALFORMED" = stop("Error - the IS matrix you're trying to import is
                            malformed, aborting")
         )
  isadf <- .messy_to_tidy(isadf)
  isadf
}


#' Internal function to convert a messy matrix to a tidy data frame
#'
#' @description The function uses the suite of functions provided by the
#' tidyverse to have a more dense and ordered structure. This function is not
#' exported and should be called in other importing functions.
#'
#' @param df ISADataFrame to convert to tidy
#'
#' @return a tidy ISADataFrame
#' @importFrom rlang .data
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr arrange all_of
#' @importFrom forcats fct_inseq as_factor
.messy_to_tidy <- function(df) {
  stopifnot(is.ISADataFrame(df))
  exp_cols <- which(!(colnames(df) %in% c(mandatoryVars(df), metadata(df))))
  isadf_tidy <- df %>% tidyr::pivot_longer(cols = dplyr::all_of(exp_cols),
                                           names_to = "ExperimentID",
                                           values_to = "Value",
                                           values_drop_na = TRUE) %>%
        dplyr::arrange(forcats::fct_inseq(forcats::as_factor(.data$chr)))
  new_meta <- c(metadata(df), "ExperimentID")
  attr(isadf_tidy, "metadata") <- new_meta
  isadf_tidy
}

#' Internal function to auto-detect the type of IS based on the headers.
#'
#' @param df the data frame to inspect
#'
#' @return one value among:
#' * "OLD" : for old-style matrices that had only one column holding all genomic
#' coordinates
#' * "NEW_ANNOTATED" :  for the classic Vispa2 annotated matrices
#' * "NEW_NOTANN" : for Vispa2 not annotated matrices
#' * "MALFORMED" : in any other case
.auto_detect_type <- function(df) {
  if ("IS_genomicID" %in% colnames(df) &
      all(!(c("chr", "integration_locus", "strand") %in% colnames(df)))) {
    return("OLD")
  }
  if (all(c("chr", "integration_locus", "strand") %in% colnames(df))) {
    if (all(c("GeneName", "GeneStrand") %in% colnames(df))) {
      return("NEW_ANNOTATED")
    } else {
      return("NEW_NOTANN")
    }
  }
  return("MALFORMED")
}
