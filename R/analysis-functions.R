#------------------------------------------------------------------------------#
# Analysis functions
#------------------------------------------------------------------------------#
#' Computes the abundance of every integration in the sample.
#'
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
#' package = "ISAnalytics")
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
        dplyr::mutate(AbsAbundance = .data$Value/.data$QuantificationSum) %>%
        dplyr::select(-c(.data$QuantificationSum))
    if (percentage == TRUE) {
        abundance_df <- abundance_df %>%
            dplyr::mutate(
                PercAbundance = .data$AbsAbundance * 100
            )
    }
    abundance_df
}
