#------------------------------------------------------------------------------#
# Re-calibration functions
#------------------------------------------------------------------------------#
#' Scans input matrix to find and merge near integration sites.
#'
#' This function scans the input integration matrix to detect eventual
#' integration sites that are too "near" to each other and merges them
#' into single integration sites that have a value equal to the sum
#' of values.
#'
#' @details The whole matrix is scanned with a sliding window mechanism,
#' a window being a group of 3 integration sites (or rows which share
#' the same chromosome and optionally strand). The distance between
#' two integration sites that share the same chr (and strand optionally),
#' is calculated in the window by comparing the 2 external integrations
#' with the central one.
#'
#' @param x A single integration matrix
#' @param threshold A single numeric value, represents the number of bases
#' under which two integrations are deemed "near"
#' @param keep_criteria While scanning, which integration should be kept?
#' The 3 possible choices for this parameter are:
#' * "max_value": keep the integration site which has the highest value
#' (and collapse other values on that integration). It is highly reccomended,
#' if you choose this criteria, to express a second preference in case
#' values can't be compared, for example c("max_value", "keep_first")
#' * "keep_first": keeps the first integration in the considered
#'  window of 3
#' * "keep_central": keeps the central integration in the considered
#' window of 3
#' @param strand_specific Should strand be considered? If yes,
#' for example these two integration sites
#' c(chr = "1", strand = "+", integration_locus = 14568) and
#' c(chr = "1", strand = "-", integration_locus = 14568) are considered
#' different and not grouped together.
#' @importFrom magrittr `%>%`
#' @importFrom tibble is_tibble
#' @importFrom dplyr group_by group_split bind_rows summarise
#' @importFrom purrr reduce map_lgl
#' @importFrom rlang .data
#' @import BiocParallel
#'
#' @family Recalibration functions
#'
#' @return An integration matrix with same or less number of rows
#' @export
#'
#' @examples
#' path <- system.file("extdata", "ex_annotated_ISMatrix.tsv.xz",
#' package = "ISAnalytics")
#' matrix <- import_single_Vispa2Matrix(path)
#' near <- compute_near_integrations(matrix)
compute_near_integrations <- function(x,
    threshold = 4,
    keep_criteria = c(
        "max_value",
        "keep_central",
        "keep_first"
    ),
    strand_specific = TRUE) {
    # Check parameters
    stopifnot(tibble::is_tibble(x))
    if (.check_mandatory_vars(x) == FALSE) {
        stop(.non_ISM_error())
    }
    if (.check_value_col(x) == FALSE) {
        stop(.missing_value_col_error())
    }
    stopifnot(is.numeric(threshold) | is.integer(threshold))
    stopifnot(length(threshold) == 1)
    stopifnot(is.character(keep_criteria))
    stopifnot(all(keep_criteria %in% c(
        "max_value", "keep_central",
        "keep_first"
    )))
    if (length(keep_criteria) == 1 && keep_criteria == "max_value") {
        keep_criteria <- c(keep_criteria, "keep_first")
    }
    stopifnot(is.logical(strand_specific) &
        length(strand_specific) == 1)
    # Pre-process data
    x <- x %>%
        dplyr::group_by(dplyr::across(mandatory_IS_vars())) %>%
        dplyr::summarise(Value = sum(.data$Value), .groups = "drop")
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
    split_to_process <- split[purrr::map_lgl(split, function(x) {
        nrow(x) > 1
    })]
    # # Set up parallel workers
    if (.Platform$OS.type == "windows") {
        p <- BiocParallel::SnowParam(stop.on.error = FALSE)
    } else {
        p <- BiocParallel::MulticoreParam(
            stop.on.error = FALSE
        )
    }
    FUN <- function(x) {
        .sliding_window(x, threshold = threshold, keep_criteria = keep_criteria)
    }
    suppressWarnings({
        result <- BiocParallel::bptry(
            BiocParallel::bplapply(split_to_process,
                FUN = FUN,
                BPPARAM = p
            )
        )
    })
    BiocParallel::bpstop(p)
    result <- purrr::reduce(result, dplyr::bind_rows)
    split_fine <- split[purrr::map_lgl(split, function(x) {
        !nrow(x) > 1
    })]
    result <- result %>% dplyr::bind_rows(split_fine)
    return(result)
}
