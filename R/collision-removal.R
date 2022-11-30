#------------------------------------------------------------------------------#
# Collision removal functions
#------------------------------------------------------------------------------#

#' Identifies and removes collisions.
#'
#' @description
#' `r lifecycle::badge("stable")`
#' A collision is an integration (aka a unique combination of the provided
#' `mandatory_IS_vars()`) which is observed in more than one
#' independent sample.
#' The function tries to decide to which independent sample should
#' an integration event be assigned to, and if no
#' decision can be taken, the integration is completely removed from the data
#' frame.
#' For more details refer to the vignette "Collision removal functionality":
#' \code{vignette("workflow_start", package = "ISAnalytics")}
#'
#' @param x Either a multi-quantification matrix (recommended) or a
#' named list of matrices (names must be quantification types)
#' @param association_file The association file imported via
#' `import_association_file()`
#' @param independent_sample_id A character vector of column names that
#' identify independent samples
#' @param date_col The date column that should be considered.
#' @param reads_ratio A single numeric value that represents the ratio that has
#' to be considered when deciding between `seqCount` value.
#' @param quant_cols A named character vector where names are
#' quantification types and
#' values are the names of the corresponding columns. The quantification
#' `seqCount` MUST be included in the vector.
#' @param max_workers Maximum number of parallel workers to distribute the
#' workload. If `NULL` (default) produces the maximum amount of workers allowed,
#' a numeric value is requested otherwise. WARNING: a higher number of workers
#' speeds up computation at the cost of memory consumption! Tune this parameter
#' accordingly.
#'
#' @section Required tags:
#' The function will explicitly check for the presence of these tags:
#'
#' ```{r echo=FALSE, results="asis"}
#' all_tags <- available_tags()
#' needed <- unique(all_tags[purrr::map_lgl(eval(rlang::sym("needed_in")),
#'  ~ "remove_collisions" %in% .x)][["tag"]])
#'  cat(paste0("* ", needed, collapse="\n"))
#' ```
#'
#' @template report_path_param
#'
#' @family Data cleaning and pre-processing
#' @importFrom rlang .data sym
#'
#' @return Either a multi-quantification matrix or a list of data frames
#' @export
#'
#' @examples
#' data("integration_matrices", package = "ISAnalytics")
#' data("association_file", package = "ISAnalytics")
#' no_coll <- remove_collisions(
#'     x = integration_matrices,
#'     association_file = association_file,
#'     report_path = NULL
#' )
#' head(no_coll)
remove_collisions <- function(x,
    association_file,
    independent_sample_id = c("ProjectID", "SubjectID"),
    date_col = "SequencingDate",
    reads_ratio = 10,
    quant_cols = c(
        seqCount = "seqCount",
        fragmentEstimate = "fragmentEstimate"
    ),
    report_path = default_report_path(),
    max_workers = NULL) {
    # Check basic parameter correctness
    association_file <- data.table::setDT(association_file)
    ### --- For matrices
    mode <- .collisions_check_input_matrices(x, quant_cols)
    ### --- For AF
    stopifnot(is.character(date_col))
    date_col <- date_col[1]
    req_tag_cols <- .collisions_check_input_af(
        association_file,
        date_col,
        independent_sample_id
    )
    req_tag_cols <- data.table::setDT(req_tag_cols)
    ### --- Other checks
    stopifnot(is.null(max_workers[1]) || is.numeric(max_workers[1]))
    max_workers <- max_workers[1]
    stopifnot(is.integer(reads_ratio) || is.numeric(reads_ratio))
    reads_ratio <- reads_ratio[1]
    seq_count_col <- quant_cols["seqCount"]
    ## Check if association file contains all info relative to content the of
    ## the matrix
    verbose <- getOption("ISAnalytics.verbose", TRUE)
    pcr_col <- pcr_id_column()
    ### - Are all sample in the matrix present in the AF?
    missing_ind <- if (!all(x[[pcr_col]] %in%
        association_file[[pcr_col]])) {
        which(!x[[pcr_col]] %in%
            association_file[[pcr_col]])
    } else {
        NULL
    }
    pre_process <- if (!is.null(missing_ind)) {
        rlang::inform(.missing_af_samples_msg(
            length(unique(x[[pcr_col]][missing_ind]))
        ))
        x[-missing_ind, ]
    } else {
        x
    }
    ## If after removing missing the matrix is empty, stop
    if (nrow(pre_process) == 0) {
        rlang::inform("Matrix is empty after removing missing samples")
        to_return <- if (mode == "LIST") {
            args <- append(list(x = x), as.list(quant_cols))
            rlang::exec(separate_quant_matrices, !!!args)
        } else {
            x
        }
        return(to_return)
    }

    ## Check if association file contains more info than matrix and
    ## keep only metadata that concerns projectIDs in the matrix
    add_samples <- .check_same_info(association_file, pre_process,
        req_tag_cols = req_tag_cols,
        indep_sample_id = independent_sample_id
    )
    association_file <- add_samples$reduced_af
    add_samples <- add_samples$miss
    if (nrow(add_samples) > 0) {
        if (verbose) {
            rlang::inform(.additional_ad_samples_msg())
        }
    } else {
        add_samples <- NULL
    }

    # Begin workflow
    ## - Join based on pcr identifier
    replicate_n_col <- req_tag_cols[eval(sym("tag")) ==
        "pcr_replicate"][["names"]]
    pool_col <- req_tag_cols[eval(sym("tag")) == "pool_id"][["names"]]
    joined <- pre_process %>%
        dplyr::left_join(association_file, by = pcr_col) %>%
        dplyr::select(dplyr::all_of(
            c(
                colnames(pre_process), date_col,
                replicate_n_col, independent_sample_id, pool_col
            )
        ))
    if (verbose) {
        rlang::inform("Identifying collisions...")
    }
    ## - Separate collisions from non-collisions
    split_df <- .identify_independent_samples(joined,
        indep_sample_id = independent_sample_id
    )

    if (nrow(split_df$collisions) == 0) {
        rlang::inform(.no_collisions_msg())
        to_return <- if (mode == "LIST") {
            args <- append(list(x = x), as.list(quant_cols))
            rlang::exec(separate_quant_matrices, !!!args)
        } else {
            x
        }
        return(to_return)
    }

    # Remove collisions
    if (verbose) {
        rlang::inform("Processing collisions...")
    }
    fixed_collisions <- .process_collisions(
        collisions = split_df$collisions,
        date_col = date_col,
        reads_ratio = reads_ratio,
        seqCount_col = seq_count_col,
        repl_col = replicate_n_col,
        ind_sample_key = independent_sample_id,
        max_workers = max_workers
    )
    reassigned <- fixed_collisions$reassigned
    removed <- fixed_collisions$removed
    fixed_collisions <- fixed_collisions$coll
    final_matr <- data.table::rbindlist(list(
        fixed_collisions,
        split_df$non_collisions
    ))
    final_matr <- final_matr[, mget(colnames(pre_process))]

    ## ---- REPORT PRODUCTION
    if (getOption("ISAnalytics.reports", TRUE) == TRUE &
        !is.null(report_path)) {
        post_joined <- association_file[final_matr, on = pcr_col]
        post_joined <- post_joined[, mget(c(
            colnames(final_matr), date_col,
            replicate_n_col, independent_sample_id, pool_col
        ))]
        summaries <- .collisions_obtain_report_summaries(
            x = x, association_file = association_file,
            quant_cols = quant_cols, missing_ind = missing_ind,
            pcr_col = pcr_col, pre_process = pre_process,
            collisions = split_df$collisions, removed = removed,
            reassigned = reassigned,
            joined = joined,
            final_matr = final_matr,
            post_joined = post_joined,
            pool_col = pool_col, replicate_n_col = replicate_n_col,
            independent_sample_id = independent_sample_id,
            seq_count_col = seq_count_col
        )
        sharing_heatmaps <- .collisions_obtain_sharing_heatmaps(
            joined = joined, independent_sample_id = independent_sample_id,
            post_joined = post_joined, report_path = report_path
        )
        report_params <- append(summaries, sharing_heatmaps)
        report_params[["additional_info"]] <- add_samples
        report_params[["sample_key"]] <- independent_sample_id
        report_params[["dynamic_cols"]] <- list(
            pcr_id = pcr_col,
            pool = pool_col
        )
        withCallingHandlers(
            {
                .produce_report("collisions",
                    params = report_params,
                    path = report_path
                )
            },
            error = function(cnd) {
                rest <- findRestart("report_fail")
                invokeRestart(rest, cnd)
            }
        )
    }

    ## If input was list, return the result in list form
    final <- if (mode == "LIST") {
        args <- append(list(x = final_matr), as.list(quant_cols))
        rlang::exec(separate_quant_matrices, !!!args)
    } else {
        final_matr
    }
    if (verbose) {
        rlang::inform("Finished!")
    }
    return(final)
}

#' Re-aligns matrices of other quantification types based on the processed
#' sequence count matrix.
#'
#' @description
#' `r lifecycle::badge("stable")`
#' This function should be used to keep data consistent among the same analysis:
#' if for some reason you removed the collisions by passing only the sequence
#' count matrix to `remove_collisions()`, you should call this
#' function afterwards, providing a list of other quantification matrices.
#' NOTE: if you provided a list of several quantification types to
#' `remove_collisions()` before, there is no need to call this function.
#'
#' @details For more details on how to use collision removal functionality:
#' \code{vignette("workflow_start", package = "ISAnalytics")}
#'
#' @param sc_matrix The sequence count matrix already processed for collisions
#' via `remove_collisions()`
#' @param other_matrices A named list of matrices to re-align. Names in the list
#' must be quantification types (\code{quantification_types()}) except
#' "seqCount".
#' @param sample_column The name of the column containing the sample identifier
#'
#' @family Data cleaning and pre-processing
#' @seealso \code{\link{remove_collisions}}
#'
#' @return A named list with re-aligned matrices
#' @export
#'
#' @examples
#' data("integration_matrices", package = "ISAnalytics")
#' data("association_file", package = "ISAnalytics")
#' separated <- separate_quant_matrices(
#'     integration_matrices
#' )
#' no_coll <- remove_collisions(
#'     x = separated$seqCount,
#'     association_file = association_file,
#'     quant_cols = c(seqCount = "Value"),
#'     report_path = NULL
#' )
#' realigned <- realign_after_collisions(
#'     sc_matrix = no_coll,
#'     other_matrices = list(fragmentEstimate = separated$fragmentEstimate)
#' )
#' realigned
realign_after_collisions <- function(sc_matrix,
    other_matrices,
    sample_column = pcr_id_column()) {
    stopifnot(is.list(other_matrices) & !is.null(names(other_matrices)))
    stopifnot(all(names(other_matrices) %in% quantification_types()))
    stopifnot(is.character(sample_column))
    sample_column <- sample_column[1]
    all_ISm <- purrr::map_lgl(other_matrices, .check_mandatory_vars)
    if (!all(all_ISm)) {
        rlang::abort(.non_ISM_error())
    }
    all_campid <- purrr::map_lgl(
        other_matrices,
        ~ sample_column %in% colnames(.x)
    )
    if (!all(all_campid)) {
        rlang::abort(.missing_needed_cols(sample_column))
    }
    realigned <- purrr::map(other_matrices, function(x) {
        x %>% dplyr::semi_join(sc_matrix,
            by = c(
                mandatory_IS_vars(),
                sample_column
            )
        )
    })
    realigned
}
