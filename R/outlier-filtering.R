#------------------------------------------------------------------------------#
# Outliers filtering functions
#------------------------------------------------------------------------------#

#---- Filter outliers ---------------------------------------------------------#
#' Filter out outliers in metadata, identified by the chosen outlier test.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#' Filter out outliers in metadata by using appropriate outlier tests.
#'
#' @details
#' ## Modular structure
#' The outlier filtering functions are structured in a modular fashion.
#' There are 2 kind of functions:
#' * Outlier tests - Functions that perform some kind of calculation based
#' on inputs and flags metadata
#' * Outlier filter - A function that takes one or more outlier tests,
#' combines all the flags with a given logic and filters out rows that
#' are flagged as outliers
#'
#' This function acts as the filter. It can either take one or more outlier
#' tests as functions and call them through the argument `outlier_test`,
#' or it can take directly outputs produced by individual tests in
#' the argument `outlier_test_outputs` - if both are provided the second one
#' has priority. The second method offers a bit more freedom, since single
#' tests can be run independently and intermediate results saved and examined
#' more in detail. If more than one test is to be performed, the argument
#' `combination_logic` tells the function how to combine the flags: you can
#' specify 1 logical operator or more than 1, provided it is compatible
#' with the number of tests.
#'
#' ## Writing custom outlier tests
#' You have the freedom to provide your own functions as outlier tests. For
#' this purpose, functions provided must respect this guidelines:
#' * Must take as input the whole metadata df
#' * Must return a df containing AT LEAST the `pcr_id_col` and a logical column
#' `"to_remove"` that contains the flag
#' * The `pcr_id_col` must contain all the values originally present in the
#' metadata df
#'
#' @param metadata The metadata data frame
#' @param pcr_id_col The name of the pcr identifier column
#' @param outlier_test One or more outlier tests. Must be functions,
#' either from `available_outlier_tests()` or custom functions that
#' produce an appropriate output format (see details).
#' @param outlier_test_outputs `NULL`, a data frame or a list of data frames.
#' See details.
#' @param combination_logic One or more logical operators
#' ("AND", "OR", "XOR", "NAND", "NOR", "XNOR"). See datails.
#' @param negate If `TRUE` will return only the metadata that was flagged to
#' be removed. If `FALSE` will return only the metadata that wasn't flagged
#' to be removed.
#' @param ... Additional named arguments passed to `outliers_test`
#'
#' @template report_path_param
#'
#' @family Data cleaning and pre-processing
#' @return A data frame of metadata which has less or the same amount of rows
#'
#' @export
#' @examples
#' data("association_file", package = "ISAnalytics")
#' filtered_af <- outlier_filter(association_file,
#'     key = "BARCODE_MUX",
#'     report_path = NULL
#' )
#' head(filtered_af)
outlier_filter <- function(metadata,
    pcr_id_col = pcr_id_column(),
    outlier_test = c(outliers_by_pool_fragments),
    outlier_test_outputs = NULL,
    combination_logic = c("AND"),
    negate = FALSE,
    report_path = default_report_path(),
    ...) {
    stopifnot(is.data.frame(metadata))
    data.table::setDT(metadata)
    stopifnot(is.character(pcr_id_col))
    pcr_id_col <- pcr_id_col[1]
    if (!pcr_id_col %in% colnames(metadata)) {
        rlang::abort(.missing_af_needed_cols(pcr_id_col))
    }
    stopifnot(is.null(outlier_test) ||
        all(purrr::map_lgl(outlier_test, is.function)))
    stopifnot(is.null(outlier_test_outputs) ||
        is.data.frame(outlier_test_outputs) ||
        all(purrr::map_lgl(outlier_test_outputs, is.data.frame)))
    stopifnot(is.logical(negate))
    params_for_report <- list(
        dyn_vars = list(pcr_id = pcr_id_col),
        input_stats = list(
            nrow = nrow(metadata),
            n_samples = length(
                unique(
                    metadata[[pcr_id_col]]
                )
            )
        )
    )
    mode <- "CALL"
    test_num <- 1
    if (is.null(outlier_test) & is.null(outlier_test_outputs)) {
        err_msg <- c(paste(
            "One between `outlier_test` and",
            "`outlier_test_outputs` should not be NULL"
        ),
        i = "See documentation with `?outlier_filter`"
        )
        rlang::abort(err_msg, class = "no_outlier_tests")
    } else if (!is.null(outlier_test) & !is.null(outlier_test_outputs)) {
        mode <- "RES"
        if (!is.data.frame(outlier_test_outputs)) {
            test_num <- length(outlier_test_outputs)
        } else {
            outlier_test_outputs <- list(outlier_test_outputs)
        }
    }
    if (mode == "CALL") {
        test_num <- length(outlier_test)
        .outlier_test_verify_logiop(
            outlier_test, combination_logic,
            "combination_logic"
        )
    } else {
        .outlier_test_verify_logiop(
            seq_len(test_num),
            combination_logic,
            "combination_logic"
        )
    }
    params_for_report$mode <- mode
    params_for_report$dyn_vars$operators <- combination_logic
    if (mode == "CALL") {
        ## Get varargs
        dots <- rlang::dots_list(..., .named = TRUE, .homonyms = "first")
        call_test <- function(test, vargs, df, test_name) {
            test_args_names <- do.call(rlang::fn_fmls_names, args = list(test))
            test_args_names <- test_args_names[test_args_names != "metadata"]
            test_args <- if ("report_path" %in% test_args_names) {
                append(
                    vargs[names(vargs) %in% test_args_names],
                    list(report_path = report_path)
                )
            } else {
                vargs[names(vargs) %in% test_args_names]
            }
            f_call <- rlang::call2(test, metadata = metadata, !!!test_args)
            result <- rlang::eval_tidy(f_call)
            .validate_outlier_output_format(
                result,
                unique(metadata[[pcr_id_col]]),
                pcr_id_col
            )
            data.table::setDT(result)
            return(list(
                result = result,
                call_args = test_args
            ))
        }
        if (is.null(names(outlier_test))) {
            names(outlier_test) <- paste("test", seq_along(outlier_test),
                sep = "_"
            )
        }
        outlier_test_outputs <- purrr::map2(
            outlier_test,
            names(outlier_test),
            ~ call_test(
                .x, dots, metadata,
                .y
            )
        )
        params_for_report$call_args <- purrr::map(
            outlier_test_outputs,
            ~ .x$call_args
        )
        outlier_test_outputs <- purrr::map(outlier_test_outputs, ~ .x$result)
    } else {
        purrr::walk(outlier_test_outputs, ~ {
            data.table::setDT(.x)
            .validate_outlier_output_format(
                .x, unique(metadata[[pcr_id_col]]),
                pcr_id_col
            )
        })
        if (is.null(names(outlier_test_outputs))) {
            outlier_test_outputs <- outlier_test_outputs %>%
                purrr::set_names(paste("test", seq(
                    1,
                    length(outlier_test_outputs)
                ),
                sep = "_"
                ))
        }
    }
    params_for_report$test_results <- outlier_test_outputs
    to_rem_cols <- purrr::map(
        outlier_test_outputs,
        ~ .x[, mget(c(pcr_id_col, "to_remove"))]
    )
    to_rem_total <- if (test_num == 1) {
        to_rem_cols[[1]]
        params_for_report$joint <- to_rem_cols[[1]]
    } else {
        purrr::walk2(to_rem_cols, names(to_rem_cols), ~ {
            setnames(.x, "to_remove", paste0("to_remove_", .y))
        })
        joint <- purrr::reduce(to_rem_cols, ~ {
            .y[.x, on = pcr_id_col]
        })
        column_names_to_rem <- paste0("to_remove_", names(to_rem_cols))
        joint[, c("to_remove") := list(
            .apply_flag_logic(!!!mget(column_names_to_rem),
                logic = combination_logic
            )
        )]
        params_for_report$joint <- data.table::copy(joint)
        joint[, c(column_names_to_rem) := NULL]
    }
    metadata <- to_rem_total[metadata, on = pcr_id_col]
    if (negate) {
        metadata_filtered <- metadata[
            eval(sym("to_remove")) == TRUE,
            !c("to_remove")
        ]
    } else {
        metadata_filtered <- metadata[
            eval(sym("to_remove")) == FALSE,
            !c("to_remove")
        ]
    }
    if (getOption("ISAnalytics.reports", TRUE) == TRUE &&
        !is.null(report_path)) {
        withCallingHandlers(
            {
                .produce_report(
                    report_type = "outlier_filter",
                    params = params_for_report,
                    path = report_path
                )
            },
            error = function(cnd) {
                rest <- findRestart("report_fail")
                invokeRestart(rest, cnd)
            }
        )
    }
    return(metadata_filtered)
}


#---- Outlier tests -----------------------------------------------------------#
### Outlier tests perform some calculation on data and flag outliers according
### to parameters

#' Identify and flag outliers based on pool fragments.
#'
#' @description
#' `r lifecycle::badge("stable")`
#' Identify and flag outliers based on expected number of raw reads per pool.
#'
#' @details
#' ## Modular structure
#' The outlier filtering functions are structured in a modular fashion.
#' There are 2 kind of functions:
#' * Outlier tests - Functions that perform some kind of calculation based
#' on inputs and flags metadata
#' * Outlier filter - A function that takes one or more outlier tests,
#' combines all the flags with a given logic and filters out rows that
#' are flagged as outliers
#'
#' This function is an outlier test, and calculates for each column in the key
#'
#' * The zscore of the values
#' * The tstudent of the values
#' * The the associated p-value (tdist)
#'
#' Optionally the test can be performed for each pool and a normality test
#' can be run prior the actual calculations.
#' Samples are flagged if this condition is respected:
#'
#' * tdist < outlier_p_value_threshold & zscore < 0
#'
#' If the key contains more than one column an additional flag logic can be
#' specified for combining the results.
#' Example:
#' let's suppose the key contains the names of two columns, X and Y
#' `key = c("X", "Y")`
#' if we specify the the argument `flag_logic = "AND"` then the reads will
#' be flagged based on this global condition:
#' (tdist_X < outlier_p_value_threshold & zscore_X < 0) AND
#' (tdist_Y < outlier_p_value_threshold & zscore_Y < 0)
#'
#' The user can specify one or more logical operators that will be applied
#' in sequence.
#'
#' @param metadata The metadata data frame
#' @param key A character vector of numeric column names
#' @param outlier_p_value_threshold The p value threshold for a read to be
#' considered an outlier
#' @param normality_test Perform normality test? Normality is assessed for
#' each column in the key using Shapiro-Wilk test and if the values do not
#' follow a normal distribution, other calculations are skipped
#' @param normality_p_value_threshold Normality threshold
#' @param transform_log2 Perform a log2 trasformation on values prior the
#' actual calculations?
#' @param per_pool_test Perform the test for each pool?
#' @param pool_col A character vector of the names of the columns that
#' uniquely identify a pool
#' @param min_samples_per_pool The minimum number of samples that a pool
#' needs to contain in order to be processed - relevant only if
#' `per_pool_test = TRUE`
#' @param flag_logic A character vector of logic operators to obtain a
#' global flag formula - only relevant if the key is longer than one.
#' All operators must be chosen between:
#' `r paste0(flag_logics(), collapse = ", ")`
#' @param keep_calc_cols Keep the calculation columns in the output data frame?
#' @param report_path The path where the report file should be saved.
#' Can be a folder, a file or NULL if no report should be produced.
#' Defaults to `{user_home}/ISAnalytics_reports`.
#'
#' @return A data frame of metadata with the column `to_remove`
#'
#' @importFrom rlang .data
#'
#' @family Data cleaning and pre-processing
#' @export
#' @examples
#' data("association_file", package = "ISAnalytics")
#' flagged <- outliers_by_pool_fragments(association_file,
#'     report_path = NULL
#' )
#' head(flagged)
outliers_by_pool_fragments <- function(metadata,
    key = "BARCODE_MUX",
    outlier_p_value_threshold = 0.01,
    normality_test = FALSE,
    normality_p_value_threshold = 0.05,
    transform_log2 = TRUE,
    per_pool_test = TRUE,
    pool_col = "PoolID",
    min_samples_per_pool = 5,
    flag_logic = "AND",
    keep_calc_cols = TRUE,
    report_path = default_report_path()) {
    ## Check
    pcr_id_col <- pcr_id_column()
    .outlier_pool_frag_base_checks(
        metadata, key, outlier_p_value_threshold,
        normality_test, normality_p_value_threshold,
        transform_log2, min_samples_per_pool,
        per_pool_test, pool_col,
        pcr_id_col
    )
    .outlier_test_verify_logiop(key, flag_logic, "flag_logic")
    ## Cleaning
    if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
        rlang::inform("Removing NAs from data...")
    }
    metadata <- data.table::setDT(metadata)
    predicate_na <- purrr::map(key, ~ {
        rlang::expr(!is.na(!!sym(.x)))
    }) %>% purrr::reduce(~ rlang::expr(!!.x & !!.y))
    metadata_filtered <- metadata[eval(predicate_na), ]
    if (nrow(metadata_filtered) == 0) {
        empty_meta_msg <- c("Metadata empty",
            paste(
                "Metadata has 0 rows after filtering,",
                "nothing to do"
            ),
            i = "Returning input data frame"
        )
        rlang::inform(empty_meta_msg)
        return(metadata)
    }
    removed_nas <- metadata[!metadata_filtered, on = pcr_id_col]
    ## Calculation
    calc_res <- .pool_frag_calc(
        meta = metadata_filtered,
        key = key,
        by_pool = per_pool_test,
        normality_test = normality_test,
        normality_threshold = normality_p_value_threshold,
        pool_col = pool_col,
        min_samples_per_pool = min_samples_per_pool,
        log2 = transform_log2, pcr_id_col = pcr_id_col
    )
    ## Extract dfs for report
    non_proc_samples <- if ("non_proc_samples" %in% names(calc_res)) {
        unique(calc_res$non_proc_samples[, mget(c(pool_col, pcr_id_col))])
    } else {
        NULL
    }
    removed_zeros <- calc_res$removed_zeros
    calc_res <- calc_res$metadata

    ## Obtain final df
    calc_res <- data.table::rbindlist(list(
        calc_res, removed_nas[, c("processed") := FALSE]
    ), fill = TRUE)

    ## Flag outliers
    for (k in key) {
        col <- paste("flag", k, sep = "_")
        calc_res[, c(col) := list(
            .flag_cond_outliers_pool_frag(
                proc = .SD$processed,
                tdist = .SD[[paste0("tdist_", k)]],
                zscore = .SD[[paste0("zscore_", k)]],
                outlier_threshold = outlier_p_value_threshold
            )
        )]
    }
    to_rem <- if (length(key) == 1) {
        calc_res[[paste0("flag_", key)]]
    } else {
        purrr::pmap_lgl(calc_res[, seq(
            from = length(calc_res) - length(key) + 1,
            to = length(calc_res)
        ), with = FALSE], .apply_flag_logic, logic = flag_logic)
    }
    calc_res[, c("to_remove") := to_rem]
    ## Produce report
    if (getOption("ISAnalytics.reports", TRUE) == TRUE &&
        !is.null(report_path)) {
        withCallingHandlers(
            {
                af_cols_specs <- data.table::setDT(
                    association_file_columns(TRUE)
                )
                proj <- af_cols_specs[eval(sym("tag")) == "project_id"]
                if (nrow(proj) == 0) {
                    proj <- NULL
                } else {
                    if (proj$names[1] %in% colnames(calc_res)) {
                        proj <- proj$names[1]
                    } else {
                        proj <- NULL
                    }
                }
                columns_to_get <- c(
                    "processed", "to_remove", proj,
                    pool_col, pcr_id_col,
                    unique(colnames(calc_res)[
                        dplyr::contains(match = key, vars = colnames(calc_res))
                    ])
                )
                flagged <- calc_res[, mget(columns_to_get)]
                .produce_report(
                    report_type = "outlier_flag",
                    params = list(
                        dyn_vars = list(project_id = proj, pcr_id = pcr_id_col),
                        by_pool = per_pool_test,
                        pool_col = pool_col,
                        norm_test = normality_test,
                        key = key,
                        flag_logic = flag_logic,
                        outlier_thresh = outlier_p_value_threshold,
                        log2_req = transform_log2,
                        removed_nas = removed_nas,
                        removed_zeros = removed_zeros,
                        non_proc_pools = non_proc_samples,
                        flag_df = flagged
                    ), path = report_path
                )
            },
            error = function(cnd) {
                rest <- findRestart("report_fail")
                invokeRestart(rest, cnd)
            }
        )
    }
    if (!keep_calc_cols) {
        calc_res <- calc_res[, mget(c(colnames(metadata), "to_remove"))]
    }
    calc_res
}

#' A character vector containing all the names of the currently supported
#' outliers tests that can be called in the function \link{outlier_filter}.
#'
#' @return A character vector
#'
#' @family Outlier tests
#' @export
#' @examples
#' available_outlier_tests()
available_outlier_tests <- function() {
    c("outliers_by_pool_fragments")
}
