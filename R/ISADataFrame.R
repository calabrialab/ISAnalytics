#------------------------------------------------------------------------------#
# ISADataFrame S3 class
#------------------------------------------------------------------------------#

#' Low-level efficient constructor for ISADataFrame objects.
#'
#' @description This is a devel function and should **not** be used
#' interactively.\cr\cr
#' ISADataFrame is a sub-class of tbl_df (tibble): it supports all of tibble
#' functionality and adds two attributes, `mandatoryVars` and `metadata` which
#' represent respectively the columns which are mandatory in an ISADataFrame
#' (chr, integration_locus, strand) and various annotations which are not
#' experimental data (for example GeneName, GeneStrand...).\cr
#' **NOTE**: the aim of this function is to be efficent, therefore no formal
#' correctness check of the produced data frame is performed and this is why
#' this function should never be called interactively.\cr\cr
#' From a devel perspective you could directly use this function in those
#' pieces of code where there is a certainty of having the correct input
#' parameters.
#' For more insight on this topic take a look at
#' [Hadley Wickham - Advanced R](https://adv-r.hadley.nz/s3.html#s3-classes).
#' @param x a named list, a tibble or a data.frame
#' @param mandVars a character vector containing the names of the mandatory
#' vars that must be present in the data frame
#' @param meta a character vector containing the names of the variables
#' representing metadata or annotations (optional)
#' @param ... optional arguments, to be used for those who want to extend
#' ISADataFrame
#' @param class character vector representing all the classes
#'
#' @return a new object of S3 class ISADataFrame
#' @importFrom tibble new_tibble
#' @details Note that if the constructor is supplied with a named list with
#' elements having different length, the resulting
#' ISADataFrame will have truncated length equal to the minimum of the lenghts
#' of the elements in the list. Appropriate checks
#' should be performed in validators and/or helpers.
#' @seealso [new_tibble],
#' [Hadley Wickham - Advanced R](https://adv-r.hadley.nz/s3.html#s3-classes)
#' @export
#' @examples
#' # Specifing the named list only returns an ISAdf where mandatoryVars
#' # are as defaults (chr, integration_locus, strand),
#' # and empty metadata
#' isaDf <- new_ISADataFrame(list(a = 1:10, b = 10:1))
#'
#' # You can change the mandatory variables by explicitly specifying the names
#' isaDf <- new_ISADataFrame(list(a = 1:10, b = 10:1),
#' mandVars = c("myvar1", "myvar2"))
#'
#' # You can specify the metadata columns also
#' isaDf <- new_ISADataFrame(list(a = 1:10, b = 10:1),
#' mandVars = c("myvar1", "myvar2"),
#' meta = c("m1", "m2", "m3"))
new_ISADataFrame <- function(x,
                             mandVars = c("chr", "integration_locus", "strand"),
                             meta = character(), ..., class = character()) {
    stopifnot(is.list(x))
    minLength <- min(vapply(x, length, FUN.VALUE = numeric(1)))

    tibble::new_tibble(x,
                       mandatoryVars = mandVars,
                       metadata = meta, ...,
                       nrow = minLength, class = c(class, "ISADataFrame"))
}

#' Validator for ISADataFrame objects.
#' @description This is a devel function and should **not** be used
#' interactively.
#' The validator takes an ISADataFrame as input to check if the object was built
#' correctly. More specifically:\cr
#' * Checks if the data frame contains the mandatory variables specified by
#' the `mandatoryVars` attribute
#' * Checks if the data frame contains the metadata variables specified by
#' the `metadata` attribute (if not, generates only a warning)
#' * Checks if there is at least one experimental data column: a column is
#' considered experimental data if it's name is not contained
#' both in `mandatoryVars` and `metadata` and the column contains numeric
#' values.
#'
#' If all checks pass the function returns TRUE, otherwise some kind of error
#' is shown.
#'
#' @param x the ISADataFrame object to validate
#'
#' @return 'TRUE' if all checks pass, error otherwise
#' @export
#'
#' @examples
#' isadf <- new_ISADataFrame(list(chr = c(as.character(1:10)),
#' integration_locus = runif(10, min = 100, max = 10000),
#' strand = sample(c("+", "-"), 10, replace = TRUE),
#' exp_1 = runif(10, min = 0, max = 10000),
#' exp_2 = runif(10, min = 0, max = 10000),
#' exp_3 = runif(10, min = 0, max = 10000)
#' ))
#'
#' validate_ISADataFrame(isadf)
#'
validate_ISADataFrame <- function(x) {
    stopifnot(is.ISADataFrame(x))
    # checks if ISAdf contains the mandatory vars columns
    if (!all(vapply(X = mandatoryVars(x), FUN = is.element,
                    set = colnames(x), FUN.VALUE = logical(1)))) {
        stop(paste("Validation of ISADataFrame failed: the input data",
        "doesn't contain the mandatory variables"))
    }
    # checks if there is at least one experimental data column (column type
    # must be numeric)
    checknonnum <- .check_nonNumdata(x)
    if (checknonnum == FALSE) {
        stop(paste("Validation of ISADataFrame failed: no experimental",
                       "variables found"))
    }
    if (checknonnum == "Warning") {
        warning(paste("Validation of ISADataFrame - warning: found",
                      "experimental columns with non numeric type"))
    }
    # checks if the specified metadata are present in the data frame
    checkM <- .check_metadata(x)
    if (checkM == FALSE) {
        warning(paste("Validation of ISADataFrame - warning: the input",
                      "data doesn't contain the specified metadata columns"))
    }

    return(TRUE)
}




#' Helper function to obtain ISADataFrame object.
#'
#' @description This function is intended to be used interactively and should be
#' used to build correct ISADataFrames.
#' If called with parameter `try.correct = TRUE` the function is able to catch
#' and correct minor issues such as:
#' * When provided with a named list as a parameter, if the elements do not have
#' the same length the shortest are filled with NAs to match the longest element
#' * When there are metadata attributes declared that are not present in the
#' data frame they're removed
#' * When non-numeric columns are detected but are not declared as metadata
#' they're added to the metadata attribute.
#'
#' Errors will be thrown in at least 2 cases:
#' * The mandatoryVars are not included in the data frame
#' * There are no experimental data columns. Note that experimental data columns
#' must be numeric.
#' @param x a named list, a tibble or a data.frame
#' @param metadata the metadata fields that are present in the table (should be
#' all variables that are not mandatory and are not experimental data)
#' @param try.correct if set to TRUE is able to fix minor issues
#'
#' @return a properly built ISADataFrame
#' @export
#' @importFrom rlang env_bind env_parent
#' @importFrom tibble is_tibble
#' @examples
#' aListWithSomeIssues <- list(chr = c(as.character(1:10)),
#' integration_locus = runif(10, min = 100, max = 10000),
#' strand = sample(c("+", "-"), 10, replace = TRUE),
#' meta1 = rep_len("m1", 10),
#' nonNumericdata = rep_len("random", 10),
#' exp_1 = runif(5, min = 0, max = 10000),
#' exp_2 = runif(10, min = 0, max = 10000),
#' exp_3 = runif(8, min = 0, max = 10000))
#'
#' isadf <- ISADataFrame(aListWithSomeIssues, metadata = c("meta1"),
#' try.correct = TRUE)
#' head(isadf)
ISADataFrame <- function(x, metadata = character(), try.correct = TRUE) {
    stopifnot(is.list(x) | is.data.frame(x) | is_tibble(x) | is.ISADataFrame(x))
    if (is.list(x)) {
        lengths <- vapply(x, length, FUN.VALUE = numeric(1))
        equalLeng <- (lengths == lengths[[1]])
        if (!all(equalLeng)) {
            if (!try.correct) {
                stop(paste("Error in ISADataFrame(): list provided as input",
                     "has elements with different lengths.",
                     "Try try.correct = TRUE"))
            } else {
                max <- max(lengths)
                x <- lapply(x, FUN = function(x) {
                    if (length(x) < max) {
                        append(x, rep_len(NA, max - length(x)))
                    } else {
                        x
                    }
                })
                message(paste("Warning - introduced NAs to fix issues in",
                        "provided list"))
            }
        }
    }
    isaDf <- new_ISADataFrame(x, meta = metadata)
    resultValidation <- withCallingHandlers({
            validate_ISADataFrame(isaDf)
        },
        error = function(cond) {
            stop(paste("Couldn't build ISADataFrame from provided input.",
                       "Aborting.", conditionMessage(cond)))
        },
        warning = function(cond) {
            if (try.correct) {
                if (conditionMessage(cond)==paste("Validation of ISADataFrame",
                    "- warning: the input data doesn't contain the specified",
                    "metadata columns")) {
                    present <- which(vapply(X = metadata(isaDf),
                                            FUN = is.element,
                                            set = colnames(isaDf),
                                            FUN.VALUE = logical(1)))
                    attr(isaDf, "metadata") <- names(present)
                    rlang::env_bind(env_parent(), isaDf = isaDf)
                    message(paste("Auto-corrected: the input data does not",
                            "contain the specified metadata columns"))
                    invokeRestart("muffleWarning")
                }
                if (conditionMessage(cond)==paste("Validation of ISADataFrame",
                    "- warning: found experimental columns with non numeric",
                    "type")) {
                    nnum <- .find_nonNumData(isaDf)
                    attr(isaDf, "metadata") <- c(metadata(isaDf), names(nnum))
                    if (.check_atLeastOneExp(isaDf)) {
                        rlang::env_bind(env_parent(), isaDf = isaDf)
                        message(paste("Auto-corrected: found experimental",
                        "columns with non numeric type"))
                        invokeRestart("muffleWarning")
                    } else {
                        stop(paste("Validation of ISADataFrame failed:",
                        "no experimental variables found"))
                    }
                }
            } else {
                stop(paste("Could not build ISADataFrame - warnings thrown:",
                           conditionMessage(cond),
                           paste("\n", "Try auto-correct function by using",
                               "try.correct = TRUE")))
            }
        }
    )
    isaDf
}


#' Gets the value of the attribute mandatoryVars.
#'
#' @param x an ISADataFrame object
#'
#' @return a character vector
#' @export
#'
#' @examples
#' expList <- list(chr = c(as.character(1:10)),
#' integration_locus = runif(10, min = 100, max = 10000),
#' strand = sample(c("+","-"), 10, replace = TRUE),
#' meta1 = rep_len("m1", 10),
#' exp_1 = runif(10, min = 0, max = 10000),
#' exp_2 = runif(10, min = 0, max = 10000),
#' exp_3 = runif(10, min = 0, max = 10000))
#' isadf <- ISADataFrame(expList, metadata = c("meta1"))
#' mandatory <- mandatoryVars(isadf)
mandatoryVars <- function(x) {
    stopifnot(is.ISADataFrame(x))
    attr(x, "mandatoryVars")
}

#' Gets the value of the attribute metadata.
#'
#' @param x an ISADataFrame object
#'
#' @return a character vector
#' @export
#'
#' @examples
#' expList <- list(chr = c(as.character(1:10)),
#' integration_locus = runif(10, min = 100, max = 10000),
#' strand = sample(c("+","-"), 10, replace = TRUE),
#' meta1 = rep_len("m1", 10),
#' exp_1 = runif(10, min = 0, max = 10000),
#' exp_2 = runif(10, min = 0, max = 10000),
#' exp_3 = runif(10, min = 0, max = 10000))
#' isadf <- ISADataFrame(expList, metadata = c("meta1"))
#' meta <- metadata(isadf)
metadata <- function(x) {
    stopifnot(is.ISADataFrame(x))
    attr(x, "metadata")
}

