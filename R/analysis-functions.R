#------------------------------------------------------------------------------#
# Analysis functions
#------------------------------------------------------------------------------#
#' Computes the abundance of every integration in the sample.
#'
#' \lifecycle{maturing}
#' Abundance is obtained for every row by calculating the ratio
#' between the single value and the total value for the sample.
#'
#' @param x An integration matrix
#' @param percentage Add abundance as percentage?
#' @family Analysis functions
#'
#' @importFrom magrittr `%>%`
#' @importFrom tibble is_tibble
#' @importFrom dplyr group_by summarise left_join mutate select
#' @importFrom rlang .data
#' @return An integration matrix
#' @export
#'
#' @examples
#' path <- system.file("extdata", "ex_annotated_ISMatrix.tsv.xz",
#'     package = "ISAnalytics"
#' )
#' matrix <- import_single_Vispa2Matrix(path)
#' abundance <- compute_abundance(matrix)
compute_abundance <- function(x, percentage = TRUE) {
    ## Check parameters
    stopifnot(tibble::is_tibble(x))
    if (.check_mandatory_vars(x) == FALSE) {
        stop(.non_ISM_error())
    }
    if (.check_complAmpID(x) == FALSE) {
        stop(.missing_complAmpID_error())
    }
    if (.check_value_col(x) == FALSE) {
        stop(.missing_value_col_error())
    }
    stopifnot(is.logical(percentage) & length(percentage) == 1)
    ## Computation
    totals <- x %>%
        dplyr::group_by(.data$CompleteAmplificationID) %>%
        dplyr::summarise(
            QuantificationSum = sum(.data$Value)
        )
    abundance_df <- x %>%
        dplyr::left_join(totals, by = "CompleteAmplificationID") %>%
        dplyr::mutate(AbsAbundance = .data$Value / .data$QuantificationSum) %>%
        dplyr::select(-c(.data$QuantificationSum))
    if (percentage == TRUE) {
        abundance_df <- abundance_df %>%
            dplyr::mutate(
                PercAbundance = .data$AbsAbundance * 100
            )
    }
    abundance_df
}


#' obtain a single integration matrix from individual quantification
#' matrices.
#'
#' \lifecycle{maturing}
#' Takes a list of integration matrices referring to different qunatification
#' types and merges them in a single data frame that has multiple
#' value columns, each renamed according to their quantification type
#' of reference.
#'
#' @param x A named list of integration matrices, ideally obtained via
#' \link{import_parallel_Vispa2Matrices_interactive} or
#' \link{import_parallel_Vispa2Matrices_auto}. Names must be
#' quantification types.
#' @param fragmentEstimate The name of the output column for fragment
#' estimate values
#' @param seqCount The name of the output column for sequence
#' count values
#' @param barcodeCount The name of the output column for barcode count
#' values
#' @param cellCount The name of the output column for cell count values
#' @param ShsCount The name of the output column for Shs count values
#'
#' @importFrom purrr walk map2 reduce
#' @importFrom dplyr rename full_join intersect
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data `:=`
#'
#' @family Analysis functions
#'
#' @seealso \link{quantification_types}
#'
#' @return A tibble
#' @export
#'
#' @examples
#' op <- options("ISAnalytics.widgets" = FALSE)
#' path <- system.file("extdata", "ex_association_file.tsv",
#'     package = "ISAnalytics"
#' )
#' root_pth <- system.file("extdata", "fs.zip", package = "ISAnalytics")
#' root <- unzip_file_system(root_pth, "fs")
#' matrices <- import_parallel_Vispa2Matrices_auto(
#'     path, root,
#'     c("fragmentEstimate", "seqCount"), "annotated", 2, NULL, "ANY"
#' )
#' total_matrix <- comparison_matrix(matrices)
#' options(op)
comparison_matrix <- function(x,
    fragmentEstimate = "fragmentEstimate",
    seqCount = "seqCount",
    barcodeCount = "barcodeCount",
    cellCount = "cellCount",
    ShsCount = "ShsCount") {
    stopifnot(is.list(x) & !is.data.frame(x))
    stopifnot(all(names(x) %in% quantification_types()))
    purrr::walk(x, function(m) {
        mand <- .check_mandatory_vars(m)
        amp <- .check_complAmpID(m)
        val <- .check_value_col(m)
        if (any(c(mand, amp, val) == FALSE)) {
            stop(.non_ISM_error())
        }
    })
    stopifnot(is.character(fragmentEstimate) & length(fragmentEstimate) == 1)
    stopifnot(is.character(seqCount) & length(seqCount) == 1)
    stopifnot(is.character(barcodeCount) & length(barcodeCount) == 1)
    stopifnot(is.character(cellCount) & length(cellCount) == 1)
    stopifnot(is.character(ShsCount) & length(ShsCount) == 1)
    x <- purrr::map2(x, names(x), function(matrix, quant_type) {
        matrix %>% dplyr::rename({{ quant_type }} := .data$Value)
    })
    result <- purrr::reduce(x, function(matrix1, matrix2) {
        commoncols <- dplyr::intersect(colnames(matrix1), colnames(matrix2))
        matrix1 %>%
            dplyr::full_join(matrix2, by = commoncols)
    })
    if (any(is.na(result)) & getOption("ISAnalytics.verbose") == TRUE) {
        message(.nas_introduced_msg())
    }
    result
}


#' Separate a multiple-quantification matrix into single quantification
#' matrices.
#'
#' \lifecycle{maturing}
#' The function separates a single multi-quantification integration
#' matrix, obtained via \link{comparison_matrix}, into single
#' quantification matrices as a named list of tibbles.
#'
#' @param x Single integration matrix with multiple quantification
#' value columns, likely obtained via \link{comparison_matrix}.
#' @param fragmentEstimate Name of the fragment estimate values column
#' in input
#' @param seqCount Name of the sequence count values column
#' in input
#' @param barcodeCount Name of the barcode count values column
#' in input
#' @param cellCount Name of the cell count values column
#' in input
#' @param ShsCount Name of the shs count values column
#' in input
#'
#' @importFrom purrr is_empty map set_names
#' @importFrom dplyr rename
#' @importFrom magrittr `%>%`
#'
#' @family Analysis functions
#'
#' @return A named list of tibbles, where names are quantification types
#' @seealso \link{quantification_types}
#' @export
#'
#' @examples
#' op <- options("ISAnalytics.widgets" = FALSE)
#' path <- system.file("extdata", "ex_association_file.tsv",
#'     package = "ISAnalytics"
#' )
#' root_pth <- system.file("extdata", "fs.zip", package = "ISAnalytics")
#' root <- unzip_file_system(root_pth, "fs")
#' matrices <- import_parallel_Vispa2Matrices_auto(
#'     path, root,
#'     c("fragmentEstimate", "seqCount"), "annotated", 2, NULL, "ANY"
#' )
#' total_matrix <- comparison_matrix(matrices)
#' separated_matrix <- separate_quant_matrices(total_matrix)
#' options(op)
separate_quant_matrices <- function(x, fragmentEstimate = "fragmentEstimate",
    seqCount = "seqCount",
    barcodeCount = "barcodeCount",
    cellCount = "cellCount",
    ShsCount = "ShsCount") {
    stopifnot(is.data.frame(x))
    if (.check_mandatory_vars(x) == FALSE) {
        stop(.non_ISM_error())
    }
    if (.check_complAmpID(x) == FALSE) {
        stop(.missing_complAmpID_error())
    }
    num_cols <- .find_exp_cols(x)
    if (purrr::is_empty(num_cols)) {
        stop(.missing_num_cols_error())
    }
    stopifnot(is.character(fragmentEstimate) & length(fragmentEstimate) == 1)
    stopifnot(is.character(seqCount) & length(seqCount) == 1)
    stopifnot(is.character(barcodeCount) & length(barcodeCount) == 1)
    stopifnot(is.character(cellCount) & length(cellCount) == 1)
    stopifnot(is.character(ShsCount) & length(ShsCount) == 1)
    param_col <- c(
        fragmentEstimate, seqCount, barcodeCount, cellCount,
        ShsCount
    )
    to_copy <- if (any(!num_cols %in% param_col)) {
        if (all(!num_cols %in% param_col)) {
            stop(.non_quant_cols_error())
        }
        num_cols[!num_cols %in% param_col]
    }
    num_cols <- num_cols[num_cols %in% param_col]
    annot <- if (.is_annotated(x)) {
        annotation_IS_vars()
    } else {
        character(0)
    }
    if (!purrr::is_empty(to_copy) & getOption("ISAnalytics.verbose") == TRUE) {
        message(.non_quant_cols_msg(to_copy))
    }
    separated <- purrr::map(num_cols, function(quant) {
        x[c(
            mandatory_IS_vars(), annot, "CompleteAmplificationID",
            to_copy, quant
        )] %>% dplyr::rename(Value = quant)
    }) %>% purrr::set_names(num_cols)
    separated
}
