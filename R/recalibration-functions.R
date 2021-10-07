#------------------------------------------------------------------------------#
# Re-calibration functions
#------------------------------------------------------------------------------#
#' Scans input matrix to find and merge near integration sites.
#'
#' \lifecycle{stable}
#' This function scans the input integration matrix to detect eventual
#' integration sites that are too "near" to each other and merges them
#' into single integration sites adjusting their values if needed.
#'
#' @details The whole matrix is scanned with a sliding window mechanism:
#' for each row in the integration matrix an interval is calculated
#' based on the `threshold` value, then a "look ahead" operation is
#' performed to detect subsequent rows which integration locuses fall
#' in the interval. If `CompleteAmplificationID`s of the near integrations
#' are different only the locus value (and optionally
#' `GeneName` and `GeneStrand` if the matrix is annotated) is modified,
#' otherwise rows with the same id are aggregated and values are summed.
#' The function will also
#' produce a re-calibration map: this data frame contains the reference
#' of pre-recalibration values for `chr`, `strand` and `integration_locus` and
#' the value to which that integration was changed to.
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
#' considered distinct (e.g. (1, 14576, +) and (1, 14580, +) are
#' considered distinct)
#' @param keep_criteria While scanning, which integration should be kept?
#' The 2 possible choices for this parameter are:
#' * "max_value": keep the integration site which has the highest value
#' (and collapse other values on that integration).
#' * "keep_first": keeps the first integration
#' @param strand_specific Should strand be considered? If yes,
#' for example these two integration sites
#' `(chr = "1", strand = "+", integration_locus = 14568)` and
#' `(chr = "1", strand = "-", integration_locus = 14568)` are considered
#' different and not grouped together.
#' @param value_columns Character vector, contains the names of the numeric
#' experimental columns
#' @param max_value_column The column that has to be considered for
#' searching the maximum value
#' @param map_as_file Produce recalibration map as a .tsv file?
#' @param file_path String representing the path were the file will be
#' saved. Can be either a folder or a file. Relevant only if `map_as_file` is
#' `TRUE`.
#' @importFrom magrittr `%>%`
#' @importFrom data.table rbindlist setDT setnames
#' @importFrom dplyr group_by group_split
#' @importFrom purrr reduce map_lgl map_dfr
#' @importFrom rlang .data abort arg_match
#' @importFrom BiocParallel SnowParam MulticoreParam bptry bplapply bpstop
#'
#' @family Recalibration functions
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
    keep_criteria = "max_value",
    strand_specific = TRUE,
    value_columns = c("seqCount", "fragmentEstimate"),
    max_value_column = "seqCount",
    map_as_file = TRUE,
    file_path = default_report_path()) {
    # Check parameters
    stopifnot(is.data.frame(x))
    ## Check tibble is an integration matrix
    if (!.check_mandatory_vars(x)) {
        rlang::abort(.missing_mand_vars())
    }
    ## Check it contains CompleteAmplificationID column
    if (!.check_complAmpID(x)) {
        rlang::abort(.missing_cAmp_sub_msg())
    }
    ## Check numeric columns
    stopifnot(is.character(value_columns))
    stopifnot(is.character(max_value_column))
    num_cols <- unique(c(value_columns, max_value_column[1]))
    if (!all(num_cols %in% colnames(x))) {
        rlang::abort(.missing_user_cols_error(
            num_cols[!num_cols %in% colnames(x)]
        ))
    }
    stopifnot(is.numeric(threshold) || is.integer(threshold))
    threshold <- threshold[1]
    criteria <- rlang::arg_match(keep_criteria, values = c(
        "max_value",
        "keep_first"
    ))
    stopifnot(is.logical(strand_specific))
    strand_specific <- strand_specific[1]
    stopifnot(is.logical(map_as_file))
    map_as_file <- map_as_file[1]
    if (map_as_file == TRUE) {
        stopifnot(is.character(file_path))
        file_path <- file_path[1]
    }
    ## Is x annotated?
    annotated <- .is_annotated(x)
    # Split data for parallel execution
    split <- if (strand_specific == TRUE) {
        x %>%
            dplyr::group_by(.data$chr, .data$strand) %>%
            dplyr::group_split()
    } else {
        x %>%
            dplyr::group_by(.data$chr) %>%
            dplyr::group_split()
    }
    ## Select only groups with 2 or more rows
    split_to_process <- split[purrr::map_lgl(split, ~ nrow(.x) > 1)]
    # # Set up parallel workers
    if (.Platform$OS.type == "windows") {
        p <- BiocParallel::SnowParam(
            stop.on.error = FALSE,
            progressbar = TRUE,
            tasks = length(split_to_process),
            exportglobals = FALSE
        )
    } else {
        p <- BiocParallel::MulticoreParam(
            stop.on.error = FALSE,
            progressbar = TRUE,
            tasks = length(split_to_process),
            exportglobals = FALSE
        )
    }
    FUN <- function(x) {
        .sliding_window(x,
            threshold = threshold,
            keep_criteria = criteria,
            num_cols = num_cols,
            annotated = annotated,
            max_val_col = max_value_column,
            produce_map = map_as_file
        )
    }
    ## Obtain result: list of lists
    result <- BiocParallel::bptry(
        BiocParallel::bplapply(
            split_to_process,
            FUN = FUN,
            BPPARAM = p
        )
    )
    BiocParallel::bpstop(p)
    ## Obtain single list
    recalibr_m <- purrr::map(result, ~ .x$recalibrated_matrix)
    recalibr_m <- data.table::rbindlist(recalibr_m)
    maps <- purrr::map(result, ~ .x$map)
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
        ) %>% data.table::setDT()
        mand_vars <- mandatory_IS_vars()
        map_fine <- unique(split_fine[, ..mand_vars])
        data.table::setnames(
            x = map_fine,
            old = mandatory_IS_vars(),
            new = paste0(mandatory_IS_vars(), "_before")
        )
        map_fine[, c("chr_after", "integration_locus_after", "strand_after") :=
            list(chr_before, integration_locus_before, strand_before)]
        recalibr_m <- data.table::rbindlist(list(recalibr_m, split_fine))
        maps <- data.table::rbindlist(list(maps, map_fine))
    }
    if (map_as_file) {
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
