#----------------------------------------------------------------------------------------------------------------------------#
# ISADataFrame S3 class specs and generics
#----------------------------------------------------------------------------------------------------------------------------#

#' Low-level efficient constructor for ISADataFrame objects.
#'
#' @description This is a devel function and should **not** be used interactively.\cr\cr
#' ISADataFrame is a sub-class of tbl_df (tibble): it supports all of tibble functionality and adds two
#' attributes, `mandatoryVars` and `metadata` which represent respectively the columns which are mandatory
#' in an ISADataFrame (chr, integration_locus, strand) and various annotations which are not experimental
#' data (for example GeneName, GeneStrand...).\cr
#' **NOTE**: the aim of this function is to be efficent, therefore no formal correctness check of the produced
#' data frame is performed and this is why this function should never be called interactively.\cr\cr
#' From a devel perspective you could directly use this function in those pieces of code where there is a certainty
#' of having the correct input parameters.
#' For more insight on this topic take a look at [Hadley Wickham - Advanced R](https://adv-r.hadley.nz/s3.html#s3-classes).
#' @param x a named list
#' @param mandVars a character vector containing the names of the mandatory vars that must be present in the data frame
#' @param meta a character vector containing the names of the variables representing metadata or annotations (optional)
#'
#' @return a new object of S3 class ISADataFrame
#' @importFrom tibble new_tibble
#' @seealso [new_tibble], [Hadley Wickham - Advanced R](https://adv-r.hadley.nz/s3.html#s3-classes)
#' @aliases ISADataFrame
#'
#' @examples
#' \dontrun{
#'  # Specifing the named list only returns an ISAdf where mandatoryVars are as defaults (chr, integration_locus, strand),
#'  # and empty metadata
#'  isaDf <- new_ISADataFrame(list(a=1:10, b=10:1))
#'
#'  # You can change the mandatory variables by explicitly specifying the names
#'  isaDf <- new_ISADataFrame(list(a=1:10, b=10:1), mandVars = c("myvar1", "myvar2"))
#'
#'  # You can specify the metadata columns also
#'  isaDf <- new_ISADataFrame(list(a=1:10, b=10:1), mandVars = c("myvar1", "myvar2"), meta = c("m1","m2", "m3"))
#' }
new_ISADataFrame <- function(x, mandVars = c("chr", "integration_locus", "strand"), meta = character()) {
  stopifnot(is.list(x))

  tibble::new_tibble(x, mandatoryVars = mandVars, metadata = meta, nrow = length(x[[1]]), class="ISADataFrame")
}

#' Validator for ISADataFrame objects.
#' @description This is a devel function and should **not** be used interactively.
#' The validator takes an ISADataFrame as input to check if the object was built correctly. More specifically:\cr
#' * Checks if the data frame contains the mandatory variables specified by the `mandatoryVars` attribute
#' * Checks if the data frame contains the metadata variables specified by the `metadata` attribute (if not, generates only a warning)
#' * Checks if there is at least one experimental data column: a column is considered experimental data if it's name is not contained
#' both in `mandatoryVars` and `metadata` and the column contains numeric values.
#'
#' If all checks pass the function returns TRUE, otherwise some kind of error is shown.
#'
#' @param x the ISADataFrame object to validate
#'
#' @return TRUE if all checks pass, error otherwise
#'
#' @examples
validate_ISADataFrame <- function(x) {
  stopifnot(is.ISADataFrame(x))
  # checks if ISAdf contains the mandatory vars columns
  if(!all(sapply(X=attr(x, "mandatoryVars"), FUN = is.element, set = colnames(x)))) {
    stop("Validation of ISADataFrame failed: the input data doesn't contain the mandatory variables")
  }
  # checks if the specified metadata are present in the data frame
  if(!all(sapply(X=attr(x, "metadata"), FUN = is.element, set = colnames(x)))) {
    warning("Validation of ISADataFrame - warning: the input data doesn't contain the specified metadata columns")
  }
  # checks if there is at least one experimental data column (column type must be numeric)
  mandAndMeta <- c(attr(x, "mandatoryVars"), attr(x, "metadata"))
  if(length(colnames(x))<= length(mandAndMeta)) {
    stop("Validation of ISADataFrame failed: no experimental variables found")
  }
  remCols <- colnames(x)[which(!(colnames(x) %in% mandAndMeta))]
  remColsTypes <- sapply(x[remCols], class)
  if(any(remColsTypes != "numeric")) {
    warning("Validation of ISADataFrame - warning: found experimental columns with non numeric type")
  }
  return(TRUE)
}

ISADataFrame <- function(x) {
  stopifnot(#x is not data.frame or a tibble or a named list
    )
}

`[.ISADataFrame`<- function(x, i, j, drop = TRUE) {
  new_ISADataFrame(NextMethod())
}

attributes.ISADataFrame <- function(x) {
  new_ISADataFrame(NextMethod())
}

is.ISADataFrame <- function(x) {
  any(class(x) == "ISADataFrame")
}


