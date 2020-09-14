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
        if (is.null(stats) & getOption("ISAnalytics.verbose") == TRUE) {
            message(paste("No Vispa2 stats files found for import,
                    ignoring this step"))
        } else {
            if (getOption("ISAnalytics.widgets") == TRUE) {
                report <- stats[[2]]
                stats <- stats[[1]]
                widg <- .iss_import_widget(report)
                print(widg)
            } else {
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
#' @param x A single integration matrix (tibble) or a list of imported
#' integration matrices (tibble)
#' @param association_file The imported association file
#' @param key A string or a character vector with column names of the
#' association file to take as keys
#' @param lambda The name of the aggregating function to apply to values
#' @param args Other arguments to pass to lambda, can be set to NULL
#' @param group Other variables to include in the grouping besides `key`,
#' can be set to NULL
#' @param namespace The namespace that exports `lambda`. Can be set to NULL if
#' lambda is not an exported object but rather a user-defined function in some
#' environment.
#' @param env The environment in which the function should look for symbols
#' @family Aggregate functions
#'
#' @importFrom rlang is_installed is_function eval_tidy env_get
#' @importFrom tibble is_tibble
#' @importFrom purrr map_lgl
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
#'     association_file = association_file,
#'     key = "SubjectID", args = list(na.rm = TRUE)
#' )
#' options(op)
aggregate_values_by_key <- function(
    x,
    association_file,
    key = "SubjectID",
    lambda = "sum",
    args = NULL,
    group = c(
        mandatory_IS_vars(),
        annotation_IS_vars()
    ),
    namespace = "base",
    env = .GlobalEnv) {
    stopifnot(tibble::is_tibble(x) || is.list(x))
    if (!tibble::is_tibble(x)) {
        all_int_mat <- purrr::map_lgl(x, function(y) {
            if (tibble::is_tibble(y)) {
                mand <- .check_mandatory_vars(y)
                mand
            } else {
                FALSE
            }
        })
        value_present <- purrr::map_lgl(x, .check_value_col)
    } else {
        all_int_mat <- .check_mandatory_vars(x)
        value_present <- .check_value_col(x)
    }
    if (!all(all_int_mat == TRUE)) {
        stop(.non_ISM_error())
    }
    if (!all(value_present == TRUE)) {
        stop(.missing_value_col_error())
    }
    # Check association file
    stopifnot(tibble::is_tibble(association_file))
    # Check key
    stopifnot(is.character(key))
    if (!all(key %in% colnames(association_file))) {
        stop("Key fields are missing from association file")
    }
    # Check lambda
    stopifnot(is.character(lambda) & length(lambda) == 1)
    # Check group
    stopifnot(is.character(group) | is.null(group))
    matrix_cols <- if (tibble::is_tibble(x)) {
        colnames(x)
    } else {
        colnames(x[[1]])
    }
    if (!all(group %in% c(colnames(association_file), matrix_cols))) {
        stop(paste("Grouping variables not found"))
    }
    # Check args
    stopifnot((is.list(args) & !is.null(names(args)) | is.null(args)))
    # Check namespace
    stopifnot((is.character(namespace) & length(namespace) == 1) |
        is.null(namespace))
    # Check env
    stopifnot(is.environment(env))
    # Check if lambda is an exported function of namespace (only if namespace is
    # not null)
    if (!is.null(namespace)) {
        if (!rlang::is_installed(namespace)) {
            stop(paste("Namespace", namespace, " is not installed"))
        }
        tryCatch(
            expr = {
                fn <- getExportedValue(
                    rlang::eval_tidy(namespace),
                    rlang::eval_tidy(lambda)
                )
                if (!rlang::is_function(fn)) {
                    stop("Provided lambda is not a function")
                }
            },
            error = function(cond) {
                stop(paste(
                    lambda, "is not an exported object from namespace",
                    namespace
                ))
            }
        )
    } else {
        # If no namespace check if the lambda is correctly defined in the
        # environment and it's actually a function
        fn <- rlang::env_get(env = env, nm = lambda, default = NULL)
        if (is.null(fn)) {
            stop(paste(
                "No binding found for", lambda, "in the specified",
                "environment"
            ))
        }
        if (!rlang::is_function(fn)) {
            stop(paste(
                lambda, "defined in the specified environment is not a",
                "function"
            ))
        }
    }
    if (tibble::is_tibble(x)) {
        x <- list(x)
        agg_matrix <- .aggregate_lambda(
            x, association_file, key, lambda, group,
            args, namespace, env
        )
        return(agg_matrix[[1]])
    }
    agg_matrix <- .aggregate_lambda(
        x, association_file, key, lambda, group,
        args, namespace, env
    )
    agg_matrix
}
