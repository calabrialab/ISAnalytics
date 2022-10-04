#------------------------------------------------------------------------------#
# Re-calibration functions
#------------------------------------------------------------------------------#
#' Scans input matrix to find and merge near integration sites.
#'
#' @description
#' `r lifecycle::badge("stable")`
#' This function scans the input integration matrix to detect eventual
#' integration sites that are too "near" to each other and merges them
#' into single integration sites adjusting their values if needed.
#'
#' @details
#' ## The concept of "near"
#' An integration event is uniquely identified by all fields specified in
#' the `mandatory_IS_vars()` look-up table. It can happen to find IS that
#' are formally distinct (different combination of values in the fields),
#' but that should not considered distinct in practice,
#' since they represent the same integration event - this may be due
#' to artefacts at the putative locus of the IS in the merging of multiple
#' sequencing libraries.
#'
#' We say that an integration event IS1 is near to another integration event
#' IS2 if the absolute difference of their loci is strictly lower than the
#' set `threshold`.
#'
#' ## The IS identity
#' There is also another aspect to be considered. Since the algorithm
#' is based on a sliding window mechanism, on which groups of IS should we
#' set and slide the window?
#'
#' By default, we have 3 fields in the `mandatory_IS_vars()`:
#' `r paste(mandatory_IS_vars(), sep = ", ")`, and we assume that all the fields
#' contribute to the identity of the IS. This means that IS1 and IS2 can be
#' compared only if they have the same chromosome and the same strand.
#' However, if we would like to exclude the strand of the integration from
#' our considerations then IS1 and IS2 can be selected from all the events
#' that fall on the same chromosome. A practical example:
#'
#' IS1 = `(chr = "1", strand = "+", integration_locus = 14568)`
#'
#' IS2 = `(chr = "1", strand = "-", integration_locus = 14567)`
#'
#' if `is_identity_tags = c("chromosome", "is_strand")` IS1 and IS2 are
#' considered distinct because they differ in strand, therefore no correction
#' will be applied to loci of either of the 2.
#' If `is_identity_tags = c("chromosome")` then IS1 and IS2 are considered
#' near, because the strand is irrelevant, hence one of the 2 IS will change
#' locus.
#'
#' ## Aggregating near IS
#' IS that fall in the same interval are evaluated according to the
#' criterion selected - if recalibration is necessary, rows with the same
#' sample ID are aggregated in a single row with a quantification value that
#' is the sum of all the merged rows.
#'
#' If the input integration matrix contains annotation columns, that is
#' additional columns that are not
#' * part of the mandatory IS vars (see `mandatory_IS_vars()`)
#' * part of the annotation IS vars (see `annotation_IS_vars()`)
#' * the sample identifier column
#' * the quantification column
#'
#' it is possible to specify how they should be aggregated.
#' Defaults are provided for each column type (character, integer, numeric...),
#' but custom functions can be specified as a named list, where names are
#' column names in `x` and values are functions to be applied.
#' NOTE: functions must be purrr-style lambdas and they must perform some kind
#' of aggregating operation, aka they must take a vector as input and return
#' a single value. The type of the output should match the type of the
#' target column. If you specify custom lambdas, provide defaults in the
#' special element `.defaults`.
#' Example:
#' ```r
#' list(
#'   numeric_col = ~ sum(.x),
#'   char_col = ~ paste0(.x, collapse = ", "),
#'   .defaults = default_rec_agg_lambdas()
#' )
#' ```
#' @section Required tags:
#' The function will explicitly check for the presence of these tags:
#'
#' ```{r echo=FALSE, results="asis"}
#' all_tags <- available_tags()
#' needed <- unique(all_tags[purrr::map_lgl(eval(rlang::sym("needed_in")),
#'  ~ "compute_near_integrations" %in% .x)][["tag"]])
#'  cat(paste0("* ", needed, collapse="\n"))
#' ```
#'
#' @note We do recommend to use this function in combination with
#' \link{comparison_matrix} to automatically perform re-calibration on
#' all quantification matrices.
#'
#' @param x An integration matrix
#' @param threshold A single integer that represents an absolute
#' number of bases for which two integrations are considered distinct.
#' If the threshold is set to 3 it means, provided fields `chr`
#' and `strand` are the same, integrations sites
#' which have at least 3 bases in between them are
#' considered distinct.
#' @param is_identity_tags Character vector of tags that identify the
#' integration event as distinct (except for `"locus"`). See details.
#' @param keep_criteria While scanning, which integration should be kept?
#' The 2 possible choices for this parameter are:
#' * "max_value": keep the integration site which has the highest value
#' (and collapse other values on that integration).
#' * "keep_first": keeps the first integration
#' @param value_columns Character vector, contains the names of the numeric
#' experimental columns
#' @param max_value_column The column that has to be considered for
#' searching the maximum value
#' @param sample_id_column The name of the column containing the sample
#' identifier
#' @param additional_agg_lambda A named list containing aggregating functions
#' for additional columns. See details.
#' @param max_workers Maximum parallel workers allowed
#' @param map_as_file Produce recalibration map as a .tsv file?
#' @param file_path String representing the path were the file will be
#' saved. Must be a folder. Relevant only if `map_as_file` is
#' `TRUE`.
#' @param strand_specific `r lifecycle::badge("deprecated")`
#' Deprecated, use `is_identity_tags`
#'
#' @importFrom rlang .data sym
#'
#' @family Data cleaning and pre-processing
#'
#' @return An integration matrix with same or less number of rows
#' @export
#'
#' @examples
#' data("integration_matrices", package = "ISAnalytics")
#' rec <- compute_near_integrations(
#'     x = integration_matrices, map_as_file = FALSE
#' )
#' head(rec)
compute_near_integrations <- function(x,
    threshold = 4,
    is_identity_tags = c("chromosome", "is_strand"),
    keep_criteria = c("max_value", "keep_first"),
    value_columns = c("seqCount", "fragmentEstimate"),
    max_value_column = "seqCount",
    sample_id_column = pcr_id_column(),
    additional_agg_lambda = list(.default = default_rec_agg_lambdas()),
    max_workers = 4,
    map_as_file = TRUE,
    file_path = default_report_path(),
    strand_specific = lifecycle::deprecated()) {
    ## --- Check parameters
    stopifnot(is.data.frame(x))
    stopifnot(is.numeric(threshold) || is.integer(threshold))
    threshold <- threshold[1]
    stopifnot(is.character(value_columns))
    stopifnot(is.character(max_value_column))
    num_cols <- unique(c(value_columns, max_value_column[1]))
    if (!all(num_cols %in% colnames(x))) {
        rlang::abort(.missing_needed_cols(
            num_cols[!num_cols %in% colnames(x)]
        ))
    }
    criteria <- rlang::arg_match(keep_criteria)
    stopifnot(is.character(sample_id_column))
    sample_id_column <- sample_id_column[1]
    stopifnot(is.logical(map_as_file))
    map_as_file <- map_as_file[1]
    if (map_as_file == TRUE) {
        stopifnot(is.character(file_path))
        file_path <- file_path[1]
    }
    if (lifecycle::is_present(strand_specific)) {
        lifecycle::deprecate_warn(
            when = "1.5.4",
            what = "compute_near_integrations(strand_specific)",
            with = "compute_near_integrations(is_identity_tags)",
            details = paste(
                "See the documentation for details,",
                "the argument will be likely dropped in the next",
                "release cycle."
            )
        )
        stopifnot(is.logical(strand_specific))
        strand_specific <- strand_specific[1]
        if (!"is_strand" %in% is_identity_tags) {
            is_identity_tags <- c(is_identity_tags, "is_strand")
        }
    }
    ## --- Check tags
    stopifnot(is.null(is_identity_tags) || is.character(is_identity_tags))
    required_tags <- list("locus" = c("int", "numeric"))
    if (!is.null(is_identity_tags)) {
        is_identity_tags <- is_identity_tags[is_identity_tags != "locus"]
        id_tags_types <- purrr::map(is_identity_tags, ~NULL) %>%
            purrr::set_names(is_identity_tags)
        required_tags <- append(required_tags, id_tags_types)
    }
    required_tag_cols <- .check_required_cols(
        required_tags = required_tags,
        vars_df = mandatory_IS_vars(TRUE),
        duplicate_politic = "error"
    )
    if (!all(required_tag_cols$names %in% colnames(x))) {
        rlang::abort(.missing_needed_cols(required_tag_cols$names[
            !required_tag_cols$names %in% colnames(x)
        ]))
    }
    if (!sample_id_column %in% colnames(x)) {
        rlang::abort(.missing_needed_cols(sample_id_column))
    }
    ## Is x annotated?
    annotated <- .is_annotated(x)
    ## Any additional columns present?
    additional_cols <- colnames(x)[!colnames(x) %in% c(
        required_tag_cols$names, sample_id_column, value_columns,
        annotation_IS_vars()
    )]
    find_lambda <- function(.x) {
        if (.x %in% names(additional_agg_lambda)) {
            return(additional_agg_lambda[[.x]])
        }
        defaults <- additional_agg_lambda[[".defaults"]]
        if (is.null(defaults)) {
            return(NULL)
        }
        col_type <- typeof(x[[.x]])
        return(defaults[[col_type]])
    }
    add_cols_lambdas <- purrr::map(additional_cols, find_lambda) %>%
        purrr::set_names(additional_cols)
    add_cols_lambdas <- add_cols_lambdas[purrr::map_lgl(
        add_cols_lambdas,
        ~ !is.null(.x)
    )]
    # Process
    if (is.null(is_identity_tags) || purrr::is_empty(is_identity_tags)) {
        result <- .sliding_window(
            x = x,
            threshold = threshold,
            keep_criteria = criteria,
            annotated = annotated,
            num_cols = num_cols,
            max_val_col = max_value_column,
            sample_col = sample_id_column,
            req_tags = required_tag_cols,
            add_col_lambdas = add_cols_lambdas,
            produce_map = map_as_file
        )
        recalibr_m <- result$recalibrated_matrix
        maps <- result$map
    } else {
        stopifnot(is.numeric(max_workers))
        max_workers <- max_workers[1]
        is_identity_names <- required_tag_cols %>%
            dplyr::filter(.data$tag != "locus") %>%
            dplyr::pull(.data$names)
        # Split data for parallel execution
        split <- x %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(is_identity_names))) %>%
            dplyr::group_split()
        ## Select only groups with 2 or more rows
        split_to_process <- split[purrr::map_lgl(split, ~ nrow(.x) > 1)]
        result <- .execute_map_job(
            data_list = split_to_process,
            fun_to_apply = .sliding_window,
            fun_args = list(
                threshold = threshold,
                keep_criteria = criteria,
                annotated = annotated,
                num_cols = num_cols,
                max_val_col = max_value_column,
                sample_col = sample_id_column,
                req_tags = required_tag_cols,
                add_col_lambdas = add_cols_lambdas,
                produce_map = map_as_file
            ),
            stop_on_error = TRUE,
            max_workers = max_workers
        )
        ## Obtain single list
        recalibr_m <- purrr::map(result$res, ~ .x$recalibrated_matrix)
        recalibr_m <- data.table::rbindlist(recalibr_m)
        maps <- purrr::map(result$res, ~ .x$map)
        maps <- data.table::rbindlist(maps)
        ## Add all rows that were not part of recalibration
        split_fine <- split[purrr::map_lgl(split, ~ !nrow(.x) > 1)]
        if (length(split_fine) > 0) {
            split_fine <- purrr::reduce(
                split_fine,
                function(g1, g2) {
                    g1 <- data.table::setDT(g1)
                    g2 <- data.table::setDT(g2)
                    data.table::rbindlist(list(g1, g2))
                }
            )
            data.table::setDT(split_fine)
            map_fine <- unique(split_fine[, mget(required_tag_cols$names)])
            data.table::setnames(
                x = map_fine,
                old = required_tag_cols$names,
                new = paste0(required_tag_cols$names, "_before")
            )
            map_fine[, c(
                paste0(required_tag_cols$names, "_after")
            ) :=
                purrr::map(
                    required_tag_cols$names,
                    ~ eval(sym(paste0(.x, "_before")))
                )]
            recalibr_m <- data.table::rbindlist(list(recalibr_m, split_fine))
            maps <- data.table::rbindlist(list(maps, map_fine))
        }
    }
    if (map_as_file & !is.null(file_path)) {
        ### Manage file
        withCallingHandlers(
            {
                .write_recalibr_map(maps, file_path)
            },
            error = function(e) {
                r <- findRestart("skip_write")
                if (is.null(r)) {
                    return()
                }
                invokeRestart(r)
            }
        )
    }
    return(recalibr_m)
}

#' Defaults for column aggregations in `compute_near_integrations()`.
#'
#' @return A named list of lambdas
#' @export
#'
#' @examples
#' default_rec_agg_lambdas()
default_rec_agg_lambdas <- function() {
    list(
        character = ~ paste0(.x, collapse = ";"),
        integer = ~ sum(.x, na.rm = TRUE),
        double = ~ sum(.x, na.rm = TRUE),
        logical = ~ all(.x)
    )
}
