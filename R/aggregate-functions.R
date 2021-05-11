#------------------------------------------------------------------------------#
# Aggregate functions
#------------------------------------------------------------------------------#

#' Performs aggregation on metadata contained in the association file.
#'
#' \lifecycle{maturing}
#' Groups metadata by the specified grouping keys and returns a
#' summary of info for each group. For more details on how to use this function:
#' \code{vignette("Working with aggregate functions", package = "ISAnalytics")}
#'
#' @param association_file The imported association file
#' (via link{import_association_file})
#' @param grouping_keys A character vector of column names to form a group
#' @param aggregating_functions A data frame containing specifications
#' of the functions to be applied to columns in the association file during
#' aggregation. It defaults to \link{default_meta_agg}. The structure of
#' this data frame should be maintained if the user wishes to change the
#' defaults.
#' @param import_stats `r lifecycle::badge("deprecated")` The import
#' of VISPA2 stats has been moved to its dedicated function,
#' see \link{import_Vispa2_stats}.
#' @family Aggregate functions
#' @importFrom rlang abort inform
#' @importFrom purrr is_empty
#' @import lifecycle
#'
#' @return An aggregated data frame
#' @export
#'
#' @examples
#' op <- options("ISAnalytics.widgets" = FALSE, "ISAnalytics.verbose" = FALSE)
#' path_AF <- system.file("extdata", "ex_association_file.tsv",
#'     package = "ISAnalytics"
#' )
#' root_correct <- system.file("extdata", "fs.zip", package = "ISAnalytics")
#' root_correct <- unzip_file_system(root_correct, "fs")
#' association_file <- import_association_file(path_AF, root_correct,
#'     dates_format = "dmy"
#' )
#' aggregated_meta <- aggregate_metadata(association_file)
#' options(op)
aggregate_metadata <- function(association_file,
    grouping_keys = c(
        "SubjectID",
        "CellMarker",
        "Tissue",
        "TimePoint"
    ),
    aggregating_functions = default_meta_agg(),
    import_stats = lifecycle::deprecated()
    ) {
    # Check parameters
    stopifnot(is.data.frame(association_file))
    stopifnot(!is.null(grouping_keys))
    stopifnot(is.character(grouping_keys))
    keys_missing <- grouping_keys[!grouping_keys %in%
                                      colnames(association_file)]
    if (!purrr::is_empty(keys_missing)) {
        rlang::abort(.missing_user_cols_error(keys_missing))
    }
    if (lifecycle::is_present(import_stats)) {
        lifecycle::deprecate_warn(
            when = "1.1.11",
            what = "aggregate_metadata(import_stats)",
            details = c("Import Vispa2 stats functionality moved",
                        i = paste("Please use `import_Vispa2_stats()`",
                                  "or",
                                  "`import_association_file(import_iss = TRUE)`",
                                  "instead."))
        )
    }
    aggregated <- .aggregate_meta(association_file = association_file,
                                  grouping_keys = grouping_keys,
                                  function_tbl = aggregating_functions)
    if (is.null(aggregated)) {
        rlang::inform(paste("No columns in `aggregating_functions$Column`",
                            "was found in column names of the association",
                            "file. Nothing to return."))
    }
    aggregated
}


#' Default metadata aggregation function table
#'
#' A default columns-function specifications for \link{aggregate_metadata}
#'
#' @details
#' This data frame contains four columns:
#'
#' * `Column`: holds the name of the column in the association file that
#' should be processed
#' * `Function`: contains either the name of a function (e.g. mean)
#' or a purrr-style lambda (e.g. `~ mean(.x, na.rm = TRUE)`). This function
#' will be applied to the corresponding column specified in `Column`
#' * `Args`: optional additional arguments to pass to the corresponding
#' function. This is relevant ONLY if the corresponding `Function` is a
#' simple function and not a purrr-style lambda.
#' * `Output_colname`: a `glue` specification that will be used to determine
#' a unique output column name. See \link[glue]{glue} for more details.
#'
#' @importFrom tibble tribble
#' @return A data frame
#' @family Aggregate functions
#' @export
#'
#' @examples
#' default_meta_agg()
default_meta_agg <- function() {
    tibble::tribble(
        ~ Column, ~ Function, ~ Args, ~ Output_colname,
        "FusionPrimerPCRDate", ~ suppressWarnings(min(.x, na.rm = TRUE)),
        NA, "{.col}_min",
        "LinearPCRDate", ~ suppressWarnings(min(.x, na.rm = TRUE)),
        NA, "{.col}_min",
        "VCN", ~ suppressWarnings(mean(.x, na.rm = TRUE)),
        NA, "{.col}_avg",
        "ng DNA corrected", ~ suppressWarnings(mean(.x, na.rm = TRUE)),
        NA, "{.col}_avg",
        "Kapa", ~ suppressWarnings(mean(.x, na.rm = TRUE)),
        NA, "{.col}_avg",
        "ng DNA corrected", ~ sum(.x, na.rm = TRUE),
        NA, "{.col}_sum",
        "ulForPool", ~ sum(.x, na.rm = TRUE),
        NA, "{.col}_sum",
        "BARCODE_MUX", ~ sum(.x, na.rm = TRUE),
        NA, "{.col}_sum",
        "TRIMMING_FINAL_LTRLC", ~ sum(.x, na.rm = TRUE),
        NA, "{.col}_sum",
        "LV_MAPPED", ~ sum(.x, na.rm = TRUE),
        NA, "{.col}_sum",
        "BWA_MAPPED_OVERALL", ~ sum(.x, na.rm = TRUE),
        NA, "{.col}_sum",
        "ISS_MAPPED_PP", ~ sum(.x, na.rm = TRUE),
        NA, "{.col}_sum"
    )
}

#' Aggregates matrices values based on specified key.
#'
#' \lifecycle{maturing}
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
#' @param join_af_by A character vector representing the joining key
#' between the matrix and the metadata. Useful to re-aggregate already
#' aggregated matrices.
#' @family Aggregate functions
#'
#' @importFrom purrr walk set_names map_lgl
#' @importFrom rlang expr eval_tidy abort
#'
#' @return A list of tibbles or a single tibble aggregated according to
#' the specified arguments
#' @export
#'
#' @examples
#' op <- options("ISAnalytics.widgets" = FALSE, "ISAnalytics.verbose" = FALSE)
#' path_AF <- system.file("extdata", "ex_association_file.tsv",
#'     package = "ISAnalytics"
#' )
#' root_correct <- system.file("extdata", "fs.zip", package = "ISAnalytics")
#' root_correct <- unzip_file_system(root_correct, "fs")
#' association_file <- import_association_file(path_AF, root_correct,
#'     dates_format = "dmy"
#' )
#' matrices <- import_parallel_Vispa2Matrices_auto(
#'     association_file = association_file, root = NULL,
#'     quantification_type = c("fragmentEstimate", "seqCount"),
#'     matrix_type = "annotated", workers = 2, matching_opt = "ANY"
#' )
#' agg <- aggregate_values_by_key(
#'     x = matrices,
#'     association_file = association_file,
#'     value_cols = c("fragmentEstimate", "seqCount")
#' )
#' options(op)
aggregate_values_by_key <- function(x,
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
    ),
    join_af_by = "CompleteAmplificationID") {
    stopifnot(is.data.frame(x) || is.list(x))
    if (!is.data.frame(x)) {
        purrr::walk(x, function(df) {
            stopifnot(is.data.frame(df))
            if (.check_mandatory_vars(df) == FALSE) {
                rlang::abort(.non_ISM_error())
            }
            if (!all(join_af_by %in% colnames(df))) {
                rlang::abort(c(x = paste("Missing common columns",
                                         "to join metadata"),
                               i = paste("Missing: ",
                                         paste0(join_af_by[!join_af_by %in%
                                                               colnames(df)],
                                                collapse = ", "))))
            }
            if (!all(value_cols %in% colnames(df))) {
                rlang::abort(.missing_user_cols_error(
                    value_cols[!value_cols %in% colnames(df)]
                ))
            }
            is_numeric_col <- purrr::map_lgl(value_cols, function(col) {
                if (!is.double(df[[col]]) &&
                    !is.integer(df[[col]])) {
                    FALSE
                } else {
                    TRUE
                }
            }) %>% purrr::set_names(value_cols)
            if (any(!is_numeric_col)) {
                rlang::abort(.non_num_user_cols_error(
                    names(is_numeric_col)[!is_numeric_col]
                ))
            }
        })
    } else {
        if (.check_mandatory_vars(x) == FALSE) {
            rlang::abort(.non_ISM_error())
        }
        if (!all(join_af_by %in% colnames(x))) {
            rlang::abort(c(x = paste("Missing common columns",
                                     "to join metadata"),
                           i = paste("Missing: ",
                                     paste0(join_af_by[!join_af_by %in%
                                                           colnames(x)],
                                            collapse = ", "))))
        }
        if (!all(value_cols %in% colnames(x))) {
            rlang::abort(.missing_user_cols_error(
                value_cols[!value_cols %in% colnames(x)]
            ))
        }
        is_numeric_col <- purrr::map_lgl(value_cols, function(col) {
            if (!is.double(x[[col]]) &&
                !is.integer(x[[col]])) {
                FALSE
            } else {
                TRUE
            }
        }) %>% purrr::set_names(value_cols)
        if (any(!is_numeric_col)) {
            rlang::abort(.non_num_user_cols_error(
                names(is_numeric_col)[!is_numeric_col]
            ))
        }
    }
    # Check association file
    stopifnot(is.data.frame(association_file))
    # Check key
    stopifnot(is.character(key))
    if (!all(key %in% colnames(association_file))) {
        rlang::abort(c(x = "Key fields are missing from association file"))
    }
    # Check lambda
    stopifnot(is.list(lambda))
    # Check group
    stopifnot(is.character(group) | is.null(group))
    if (is.data.frame(x)) {
        if (!all(group %in% c(colnames(association_file), colnames(x)))) {
            rlang::abort(paste("Grouping variables not found"))
        }
    } else {
        purrr::walk(x, function(df) {
            if (!all(group %in% c(colnames(association_file), colnames(df)))) {
                rlang::abort(paste("Grouping variables not found"))
            }
        })
    }
    if (is.data.frame(x)) {
        x <- list(x)
        agg_matrix <- .aggregate_lambda(
            x, association_file, key, value_cols, lambda, group, join_af_by
        )
        return(agg_matrix[[1]])
    }
    agg_matrix <- .aggregate_lambda(
        x, association_file, key, value_cols, lambda, group, join_af_by
    )
    agg_matrix
}
