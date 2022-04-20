#' @param onco_db_file Uniprot file for proto-oncogenes (see details).
#' If different from default, should be supplied as a path to a file.
#' @param tumor_suppressors_db_file Uniprot file for tumor-suppressor genes.
#' If different from default, should be supplied as a path to a file.
#' @param species One between `"human"`, `"mouse"` and `"all"`
#' @param known_onco Data frame with known oncogenes. See details.
#' @param suspicious_genes Data frame with clinical relevant suspicious
#' genes. See details.
#'
#' @details
#' ## Oncogene and tumor suppressor genes files
#' These files are included in the package for user convenience and are
#' simply UniProt files with gene annotations for human and mouse.
#' For more details on how this files were generated use the help
#' `?tumor_suppressors`, `?proto_oncogenes`
#'
#' ## Known oncogenes
#' The default values are included in this package and
#' it can be accessed by doing:
#'
#' ```{r eval=FALSE}
#' known_clinical_oncogenes()
#'
#' ```
#' If the user wants to change this parameter the input data frame must
#' preserve the column structure. The same goes for the `suspicious_genes`
#' parameter (DOIReference column is optional):
#'
#' ```{r eval=FALSE}
#' clinical_relevant_suspicious_genes()
#'
#' ```
