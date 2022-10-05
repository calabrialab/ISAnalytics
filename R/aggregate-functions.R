#------------------------------------------------------------------------------#
# Aggregate functions
#------------------------------------------------------------------------------#
#' Performs aggregation on metadata contained in the association file.
#'
#' @description
#' `r lifecycle::badge("stable")`
#' Groups metadata by the specified grouping keys and returns a
#' summary of info for each group. For more details on how to use this function:
#' \code{vignette("workflow_start", package = "ISAnalytics")}
#'
#' @param association_file The imported association file
#' (via \link{import_association_file})
#' @param grouping_keys A character vector of column names to form a grouping
#' operation
#' @param aggregating_functions A data frame containing specifications
#' of the functions to be applied to columns in the association file during
#' aggregation. It defaults to \link{default_meta_agg}. The structure of
#' this data frame should be maintained if the user wishes to change the
#' defaults.
#' @param import_stats `r lifecycle::badge("deprecated")` The import
#' of VISPA2 stats has been moved to its dedicated function,
#' see \link{import_Vispa2_stats}.
#'
#' @family Data cleaning and pre-processing
#'
#' @return An aggregated data frame
#' @export
#'
#' @examples
#' data("association_file", package = "ISAnalytics")
#' aggreg_meta <- aggregate_metadata(
#'     association_file = association_file
#' )
#' head(aggreg_meta)
aggregate_metadata <- function(association_file,
    grouping_keys = c(
        "SubjectID",
        "CellMarker",
        "Tissue",
        "TimePoint"
    ),
    aggregating_functions = default_meta_agg(),
    import_stats = lifecycle::deprecated()) {
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
                i = paste(
                    "Please use `import_Vispa2_stats()`",
                    "or",
                    "`import_association_file(import_iss = TRUE)`",
                    "instead."
                )
            )
        )
    }
    aggregated <- .aggregate_meta(
        association_file = association_file,
        grouping_keys = grouping_keys,
        function_tbl = aggregating_functions
    )
    if (is.null(aggregated)) {
        msg <- paste(
            "No columns in `aggregating_functions$Column`",
            "was found in column names of the association",
            "file. Nothing to return."
        )
        rlang::inform(msg)
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
#' @return A data frame
#' @family Data cleaning and pre-processing
#' @export
#'
#' @examples
#' default_meta_agg()
default_meta_agg <- function() {
    tibble::tribble(
        ~Column, ~Function, ~Args, ~Output_colname,
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
        "ISS_MAPPED_OVERALL", ~ sum(.x, na.rm = TRUE),
        NA, "{.col}_sum",
        "PCRMethod", ~ paste0(unique(.x), collapse = "|"),
        NA, "{.col}",
        "NGSTechnology", ~ paste0(unique(.x), collapse = "|"),
        NA, "{.col}",
        "DNAnumber", ~ paste0(unique(.x), collapse = "|"),
        NA, "{.col}"
    )
}

#' Aggregates matrices values based on specified key.
#'
#' @description
#' `r lifecycle::badge("stable")`
#' Performs aggregation on values contained in the integration matrices based
#' on the key and the specified lambda. For more details on how to use this
#' function:
#' \code{vignette("workflow_start", package = "ISAnalytics")}
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
#' @param x A single integration matrix or a list of imported
#' integration matrices
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
#'
#' @family Data cleaning and pre-processing
#'
#' @importFrom purrr walk set_names map_lgl
#' @importFrom rlang expr eval_tidy abort
#'
#' @return A list of data frames or a single data frame aggregated according to
#' the specified arguments
#' @export
#'
#' @examples
#' data("integration_matrices", package = "ISAnalytics")
#' data("association_file", package = "ISAnalytics")
#' aggreg <- aggregate_values_by_key(
#'     x = integration_matrices,
#'     association_file = association_file,
#'     value_cols = c("seqCount", "fragmentEstimate")
#' )
#' head(aggreg)
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
    stopifnot(is.character(value_cols))
    stopifnot(is.character(key))
    stopifnot(is.null(group) || is.character(group))
    stopifnot(is.character(join_af_by))
    stopifnot(is.data.frame(association_file))
    data.table::setDT(association_file)
    stopifnot(is.list(lambda) && !is.null(names(lambda)))
    join_key_err <- c("Join key not present in in both data frames",
        x = paste(
            "Fields specified in the argument",
            "`join_af_by` must appear in both",
            "the association file and the matrix(es)"
        )
    )
    check_single_matrix <- function(df) {
        stopifnot(is.data.frame(df))
        is_numeric_col <- purrr::map_lgl(
            value_cols,
            ~ is.numeric(df[[.x]]) ||
                is.double(df[[.x]]) ||
                is.integer(df[[.x]])
        ) %>% purrr::set_names(value_cols)
        if (any(!is_numeric_col)) {
            rlang::abort(.non_num_user_cols_error(
                names(is_numeric_col)[!is_numeric_col]
            ))
        }
        if (!all(join_af_by %in% colnames(df))) {
            rlang::abort(join_key_err, class = "join_key_err_agg")
        }
    }
    if (!is.data.frame(x)) {
        purrr::walk(x, check_single_matrix)
    } else {
        check_single_matrix(x)
    }
    join_key_in_af <- all(join_af_by %in% colnames(association_file))
    if (!join_key_in_af) {
        rlang::abort(join_key_err, class = "join_key_err_agg")
    }
    agg_matrix <- if (is.data.frame(x)) {
        .aggregate_lambda(
            x, association_file, key, value_cols, lambda, group, join_af_by
        )
    } else {
        agg_matrix <- purrr::map(x, ~ .aggregate_lambda(
            .x, association_file, key, value_cols, lambda, group, join_af_by
        ))
    }
    return(agg_matrix)
}
