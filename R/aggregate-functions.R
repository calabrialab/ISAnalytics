#------------------------------------------------------------------------------#
# Aggregate functions
#------------------------------------------------------------------------------#

#' Performs aggregation on metadata contained in the association file.
#'
#' \lifecycle{experimental}
#' Groups metadata by grouping_keys and returns a summary of info for each
#' group. For more details on how to use this function:
#' \code{vignette("Working with aggregate functions", package = "ISAnalytics")}
#'
#' @param association_file The imported association file
#' (via `import_association_file`)
#' @param grouping_keys A character vector of column names to form a group
#' @param import_stats Should Vispa2 stats files be imported and included?
#' @family Aggregate functions
#' @importFrom purrr is_empty
#' @importFrom tibble is_tibble
#'
#' @return A tibble
#' @export
#'
#' @examples
#' op <- options("ISAnalytics.widgets" = FALSE)
#' path_AF <- system.file("extdata", "ex_association_file.tsv",
#'     package = "ISAnalytics"
#' )
#' root_correct <- system.file("extdata", "fs.zip", package = "ISAnalytics")
#' root_correct <- unzip_file_system(root_correct, "fs")
#' association_file <- import_association_file(path_AF, root_correct)
#' aggregated_meta <- aggregate_metadata(association_file, import_stats = FALSE)
#' options(op)
aggregate_metadata <- function(association_file,
    grouping_keys = c(
        "SubjectID",
        "CellMarker",
        "Tissue",
        "TimePoint"
    ),
    import_stats = TRUE) {
    # Check parameters
    stopifnot(tibble::is_tibble(association_file))
    min_missing <- setdiff(.min_var_set(), colnames(association_file))
    if (!purrr::is_empty(min_missing)) {
        stop(paste(c(
            "Association file is missing some of the mandatory columns:",
            min_missing
        ), collapse = "\n"))
    }
    stopifnot(!is.null(grouping_keys))
    stopifnot(is.character(grouping_keys))
    keys_missing <- setdiff(grouping_keys, colnames(association_file))
    if (!purrr::is_empty(keys_missing)) {
        stop(paste(c(
            "Some of the grouping keys you provided were not found:",
            keys_missing
        ), collapse = "\n"))
    }
    stopifnot(is.logical(import_stats) & length(import_stats) == 1)
    # Import if true
    stats <- NULL
    if (import_stats == TRUE) {
        stats <- .import_stats_iss(association_file)
        if (is.null(stats)) {
            if (getOption("ISAnalytics.verbose") == TRUE) {
                message(paste("No Vispa2 stats files found for import,
                    ignoring this step"))
            }
        } else {
            if (getOption("ISAnalytics.widgets") == TRUE) {
                withCallingHandlers(
                    {
                        report <- stats[[2]]
                        stats <- stats[[1]]
                        widg <- .iss_import_widget(report)
                        print(widg)
                    },
                    error = function(cnd) {
                        message(conditionMessage(cnd))
                        message(.widgets_error())
                        if (getOption("ISAnalytics.verbose") == TRUE) {
                            print(paste0(
                                "--- REPORT IMPORT VISPA2",
                                "STATS: FILES IMPORTED ---"
                            ))
                            print(stats[[2]],
                                width = Inf,
                                n = nrow(stats[[2]])
                            )
                        }
                    }
                )
            } else {
                if (getOption("ISAnalytics.verbose") == TRUE) {
                    print(paste0(
                        "--- REPORT IMPORT VISPA2",
                        "STATS: FILES IMPORTED ---"
                    ))
                    print(stats[[2]],
                        width = Inf,
                        n = nrow(stats[[2]])
                    )
                }
                stats <- stats[[1]]
            }
        }
    }
    aggregated <- .join_and_aggregate(association_file, stats, grouping_keys)
    aggregated
}

#' Aggregates matrices values based on specified key.
#'
#' \lifecycle{experimental}
#' Performs aggregation on values contained in the integration matrices based
#' on the key and the specified lambda. For more details on how to use this
#' function:
#' \code{vignette("Working with aggregate functions", package = "ISAnalytics")}
#'
#' @details
#' ## Setting the lambda parameter
#' The lambda parameter should always contain a named list of either
#' functions or purrr-style lambdas.
#' It is also possible to specify the namespace of the function in both
#' ways, for example:
#'
#' ```{r}
#' lambda = list(sum = sum, desc = psych::describe)
#' ```
#' Using purrr-style lambdas allows to specify arguments for the functions,
#' keeping in mind that the first parameter should always be `.x`:
#'
#' ```{r}
#' lambda = list(sum = ~sum(.x, na.rm = TRUE))
#' ```
#' It is also possible to use custom user-defined functions, keeping in
#' mind that the symbol will be evaluated in the calling environment,
#' for example if the function is called in the global environment
#' and lambda contains "foo" as a function, "foo" will be evaluated in
#' the global environment.
#'
#' ```{r}
#' foo <- function(x) {
#'   sum(x)
#' }
#'
#' lambda = list(sum = ~sum(.x, na.rm = TRUE), foo = foo)
#'
#' # Or with lambda notation
#' lambda = list(sum = ~sum(.x, na.rm = TRUE), foo = ~foo(.x))
#' ```
#' ## Constraints on aggregation functions
#' Functions passed in the lambda parameters must respect a few constraints
#' to properly work and it's the user responsibility to ensure this.
#' * Functions have to accept as input a numeric or integer vector
#' * Function should return a single value or a list/data frame:
#' if a list or a data frame is returned as a result, all the columns
#' will be added to the final data frame.
#' @param x A single integration matrix (tibble) or a list of imported
#' integration matrices (tibble)
#' @param association_file The imported association file
#' @param value_cols A character vector containing the names of the
#' columns to apply the given lambdas. Must be numeric or integer
#' columns.
#' @param key A string or a character vector with column names of the
#' association file to take as key
#' @param lambda A named list of functions or purrr-style lambdas.
#' See details section.
#' @param group Other variables to include in the grouping besides `key`,
#' can be set to NULL
#' @family Aggregate functions
#'
#' @importFrom purrr walk
#' @importFrom rlang expr eval_tidy
#'
#' @return A list of tibbles or a single tibble according to input
#' @export
#'
#' @examples
#' op <- options("ISAnalytics.widgets" = FALSE)
#' path_AF <- system.file("extdata", "ex_association_file.tsv",
#'     package = "ISAnalytics"
#' )
#' root_correct <- system.file("extdata", "fs.zip", package = "ISAnalytics")
#' root_correct <- unzip_file_system(root_correct, "fs")
#' association_file <- import_association_file(path_AF, root_correct)
#' matrices <- import_parallel_Vispa2Matrices_auto(
#'     association_file = association_file, root = NULL,
#'     quantification_type = c("fragmentEstimate", "seqCount"),
#'     matrix_type = "annotated", workers = 2, matching_opt = "ANY"
#' )
#' agg <- aggregate_values_by_key(
#'     x = matrices$seqCount,
#'     association_file = association_file
#' )
#' options(op)
aggregate_values_by_key <- function(
    x,
    association_file,
    value_cols = "Value",
    key = c(
        "SubjectID",
        "CellMarker",
        "Tissue",
        "TimePoint"
    ),
    lambda = list(sum = ~ sum(.x, na.rm = TRUE)),
    group = c(
        mandatory_IS_vars(),
        annotation_IS_vars()
    )) {
    stopifnot(is.data.frame(x) || is.list(x))
    if (!is.data.frame(x)) {
        purrr::walk(x, function(df) {
            stopifnot(is.data.frame(df))
            if (.check_mandatory_vars(df) == FALSE) {
                stop(.non_ISM_error())
            }
            if (.check_complAmpID(df) == FALSE) {
                stop(.missing_complAmpID_error())
            }
            if (!all(value_cols %in% colnames(df))) {
                stop(.missing_user_cols_error())
            }
            purrr::walk(value_cols, function(col) {
                expr <- rlang::expr(`$`(df, !!col))
                if (!is.numeric(rlang::eval_tidy(expr)) &&
                    !is.integer(rlang::eval_tidy(expr))) {
                    stop(.non_num_user_cols_error())
                }
            })
        })
    } else {
        if (.check_mandatory_vars(x) == FALSE) {
            stop(.non_ISM_error())
        }
        if (.check_complAmpID(x) == FALSE) {
            stop(.missing_complAmpID_error())
        }
        if (!all(value_cols %in% colnames(x))) {
            stop(.missing_user_cols_error())
        }
        purrr::walk(value_cols, function(col) {
            expr <- rlang::expr(`$`(x, !!col))
            if (!is.numeric(rlang::eval_tidy(expr)) &&
                !is.integer(rlang::eval_tidy(expr))) {
                stop(.non_num_user_cols_error())
            }
        })
    }
    # Check association file
    stopifnot(is.data.frame(association_file))
    # Check key
    stopifnot(is.character(key))
    if (!all(key %in% colnames(association_file))) {
        stop("Key fields are missing from association file")
    }
    # Check lambda
    stopifnot(is.list(lambda))
    # Check group
    stopifnot(is.character(group) | is.null(group))
    if (is.data.frame(x)) {
        if (!all(group %in% c(colnames(association_file), colnames(x)))) {
            stop(paste("Grouping variables not found"))
        }
    } else {
        purrr::walk(x, function(df) {
            if (!all(group %in% c(colnames(association_file), colnames(df)))) {
                stop(paste("Grouping variables not found"))
            }
        })
    }
    if (is.data.frame(x)) {
        x <- list(x)
        agg_matrix <- .aggregate_lambda(
            x, association_file, key, value_cols, lambda, group
        )
        return(agg_matrix[[1]])
    }
    agg_matrix <- .aggregate_lambda(
        x, association_file, key, value_cols, lambda, group
    )
    agg_matrix
}
