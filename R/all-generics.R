
#' Pivot_longer implementation for ISADataFrame.
#' @inheritParams tidyr::pivot_longer
#' @export
#' @importFrom tidyr pivot_longer
#' @return An ISADataFrame
#' @examples
#' path <- system.file("extdata", "ex_annotated_ISMatrix.tsv.xz",
#'     package = "ISAnalytics"
#' )
#' df <- read.csv(path,
#'     sep = "\t", check.names = FALSE,
#'     stringsAsFactors = FALSE
#' )
#' isadf <- ISADataFrame(df, metadata = c("GeneName", "GeneStrand"))
#' isadf <- tidyr::pivot_longer(isadf,
#'     cols = 6:10, names_to = "ExpID",
#'     values_to = "Val", values_drop_na = TRUE
#' )
pivot_longer.ISADataFrame <- function(data,
    cols,
    names_to = "name",
    names_prefix = NULL,
    names_sep = NULL,
    names_pattern = NULL,
    names_ptypes = list(),
    names_transform = list(),
    names_repair = "check_unique",
    values_to = "value",
    values_drop_na = FALSE,
    values_ptypes = list(),
    values_transform = list(),
    ...) {
    df <- vctrs::vec_restore(NextMethod(), data)
    df <- .fix_metadata(df)
    df
}

#' inner_join implementation for ISADataFrame.
#' @inheritParams dplyr::inner_join
#' @importFrom dplyr inner_join
#'
#' @return An ISADataFrame
#' @export
#'
#' @examples
#' path <- system.file("extdata", "ex_annotated_ISMatrix.tsv.xz",
#'     package = "ISAnalytics"
#' )
#' isadf1 <- import_single_Vispa2Matrix(path)
#' isadf2 <- import_single_Vispa2Matrix(path)
#' isadf1 %>% dplyr::inner_join(isadf2)
inner_join.ISADataFrame <- function(x,
    y,
    by = NULL,
    copy = FALSE,
    suffix = c(".x", ".y"), ...) {
    df <- vctrs::vec_restore(NextMethod(), x)
    df <- .fix_metadata(df)
    df
}

#' left_join implementation for ISADataFrame.
#' @inheritParams dplyr::left_join
#' @importFrom dplyr left_join
#'
#' @return An ISADataFrame
#' @export
#'
#' @examples
#' path <- system.file("extdata", "ex_annotated_ISMatrix.tsv.xz",
#'     package = "ISAnalytics"
#' )
#' isadf1 <- import_single_Vispa2Matrix(path)
#' isadf2 <- import_single_Vispa2Matrix(path)
#' isadf1 %>% dplyr::left_join(isadf2)
left_join.ISADataFrame <- function(x,
    y,
    by = NULL,
    copy = FALSE,
    suffix = c(".x", ".y"), ...,
    keep = FALSE) {
    df <- vctrs::vec_restore(NextMethod(), x)
    df <- .fix_metadata(df)
    df
}
#' semi_join implementation for ISADataFrame.
#' @inheritParams dplyr::semi_join
#' @importFrom dplyr semi_join
#'
#' @return An ISADataFrame
#' @export
#'
#' @examples
#' path <- system.file("extdata", "ex_annotated_ISMatrix.tsv.xz",
#'     package = "ISAnalytics"
#' )
#' isadf1 <- import_single_Vispa2Matrix(path)
#' isadf2 <- import_single_Vispa2Matrix(path)
#' isadf1 %>% dplyr::semi_join(isadf2)
semi_join.ISADataFrame <- function(x,
    y,
    by = NULL,
    copy = FALSE, ...) {
    df <- vctrs::vec_restore(NextMethod(), x)
    df <- .fix_metadata(df)
    df
}

#' anti_join implementation for ISADataFrame.
#' @inheritParams dplyr::anti_join
#' @importFrom dplyr anti_join
#'
#' @return An ISADataFrame
#' @export
#'
#' @examples
#' path <- system.file("extdata", "ex_annotated_ISMatrix.tsv.xz",
#'     package = "ISAnalytics"
#' )
#' isadf1 <- import_single_Vispa2Matrix(path)
#' isadf2 <- import_single_Vispa2Matrix(path)
#' isadf1 %>% dplyr::anti_join(isadf2)
anti_join.ISADataFrame <- function(x, y, by = NULL, copy = FALSE, ...) {
    df <- vctrs::vec_restore(NextMethod(), x)
    df <- .fix_metadata(df)
    df
}

#' select implementation for ISADataFrame.
#' @inheritParams dplyr::select
#' @importFrom dplyr select
#'
#' @return An ISADataFrame
#' @export
#'
#' @examples
#' path <- system.file("extdata", "ex_annotated_ISMatrix.tsv.xz",
#'     package = "ISAnalytics"
#' )
#' isadf1 <- import_single_Vispa2Matrix(path)
#' isadf1 %>% dplyr::select(chr)
select.ISADataFrame <- function(.data, ...) {
    df <- vctrs::vec_restore(NextMethod(), .data)
    df <- .fix_metadata(df)
    df
}

#' Is the object an ISADataFrame?
#'
#' @param x an object
#'
#' @importFrom methods is
#' @return TRUE or FALSE
#' @export
#'
#' @examples
#' is.ISADataFrame(1:10)
is.ISADataFrame <- function(x) {
    is(x, "ISADataFrame")
}

#' Printing ISADataFrames
#'
#' @param x an ISADataFrame object
#' @param ... optional arguments to print
#'
#' @export
#' @return Nothing
#'
#' @examples
#' expList <- list(
#'     chr = c(as.character(1:10)),
#'     integration_locus = runif(10, min = 100, max = 10000),
#'     strand = sample(c("+", "-"), 10, replace = TRUE),
#'     meta1 = rep_len("m1", 10),
#'     exp_1 = runif(10, min = 0, max = 10000),
#'     exp_2 = runif(10, min = 0, max = 10000),
#'     exp_3 = runif(10, min = 0, max = 10000)
#' )
#' isadf <- ISADataFrame(expList, metadata = c("meta1"))
#' print(isadf)
print.ISADataFrame <- function(x, ...) {
    cat(
        "mandatoryVars: ",
        paste0(mandatoryVars(x)[seq_len(length(mandatoryVars(x)) - 1)], ", "),
        mandatoryVars(x)[length(mandatoryVars(x))],
        "\n"
    )
    if (length(metadata(x)) > 1) {
        cat("metadata: ", paste0(
            metadata(x)[seq_len(length(metadata(x)) - 1)],
            ", "
        ), metadata(x)[length(metadata(x))], "\n")
    } else {
        cat("metadata: ", metadata(x), "\n")
    }
    NextMethod(...)
}


#' Implementation of vec_restore for ISADataFrame.
#' @inheritParams vctrs::vec_restore
#' @return See official documentation at \code{
#' \link[vctrs:vec_proxy]{vec_restore}}
#' @export
vec_restore.ISADataFrame <- function(x, to, ...) {
    new_ISADataFrame(x, meta = attr(to, "metadata"))
}
