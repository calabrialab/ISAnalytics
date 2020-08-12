#------------------------------------------------------------------------------#
# Importing functions
#------------------------------------------------------------------------------#

#' Import a single integration matrix from file
#'
#' @description This function allows to read and import an integration matrix
#' produced as the output of Vispa2 pipeline and converts it to a tidy
#' ISADataFrame.
#'
#' @param path The path to the file on disk
#'
#' @return A tidy ISADataFrame
#' @family Import functions
#' @importFrom tidyr separate
#' @importFrom utils read.csv
#' @details The import series of functions is designed to work in combination
#' with the use of Vispa2 pipeline, please refer to this article for more
#' details: \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5702242/}{VISPA2:
#' A Scalable Pipeline for High-Throughput Identification and Annotation of
#' Vector Integration Sites}.
#' For more details on how to properly use these functions, refer to the
#' vignette - vignette("how_to_import_functions")
#' @export
#'
#' @examples
#' path_to_file <- system.file("extdata", "ex_annotated_ISMatrix.tsv.xz",
#'     package = "ISAnalytics"
#' )
#' isa_dataframe <- import_single_Vispa2Matrix(path_to_file)
import_single_Vispa2Matrix <- function(path) {
    stopifnot(!missing(path) & is.character(path))
    if (!file.exists(path)) {
        stop(paste("File not found at", path))
    }
    df <- read.csv(path,
        sep = "\t", header = TRUE, fill = TRUE,
        check.names = FALSE, stringsAsFactors = FALSE
    )

    df_type <- .auto_detect_type(df)
    switch(df_type,
        "OLD" = {
            df <- df %>% tidyr::separate(
                col = .data$IS_genomicID,
                into = c(
                    "chr", "integration_locus",
                    "strand"
                ),
                sep = "_", remove = TRUE,
                convert = TRUE
            )
            isadf <- ISADataFrame(df)
        },
        "NEW_ANNOTATED" = {
            isadf <- ISADataFrame(df, metadata = c("GeneName", "GeneStrand"))
        },
        "NEW_NOTANN" = {
            isadf <- ISADataFrame(df)
        },
        "MALFORMED" = stop(paste(
            "Error - the IS matrix you're trying to",
            "import is malformed, aborting"
        ))
    )
    isadf <- .messy_to_tidy(isadf)
    isadf
}


#' Import the association file from disk
#'
#' @description Imports the association file and immediately performs a check on
#' the file system starting from the root to assess the alignment between the
#' two.
#' @param path The path on disk to the association file.
#' @param root The path on disk of the root folder of Vispa2 output.
#' See details.
#' @param tp_padding Timepoint padding, indicates the number of digits of the
#' "Timepoint" column once imported. Fills the content with 0s up to the length
#' specified (ex: 1 becomes 0001 with a tp_padding of 4)
#' @param dates_format A single string indicating how dates should be parsed.
#' Must be a value in: \code{date_formats()}
#' @family Import functions
#' @return A tibble with the contents of the association file plus a column
#' containing the path in the file system for every project and pool if found.
#' @details The import series of functions is designed to work in combination
#' with the use of Vispa2 pipeline, please refer to this article for more
#' details: \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5702242/}{VISPA2:
#' A Scalable Pipeline for High-Throughput Identification
#' and Annotation of Vector Integration Sites}.\cr
#' The pipeline automatically produces an hierarchical structure in the file
#' system which follows this schema:
#' * /root_folder
#'   * Optional intermediate folders
#'     * ProjectID\cr
#'       |_bam\cr
#'       |_bcmuxall\cr
#'       |_bed\cr
#'       |_iss\cr
#'       |_quality\cr
#'       |_report\cr
#'       |_quantification\cr
#'        *|___concatenatePoolIDSeqRun\cr
#'
#' For each ProjectID there may be several nested PoolIDs. The alignment
#' function only looks for PoolIDs in the quantification folder, since it's
#' the location of the matrices to import.
#' For more details on how to properly use these functions, refer to the
#' vignette - vignette("how_to_import_functions")
#' @export
#'
#' @seealso \code{\link{date_formats}}
#' @examples
#' path <- system.file("extdata", "ex_association_file.tsv",
#'     package = "ISAnalytics"
#' )
#' root_pth <- system.file("extdata", "fs.zip", package = "ISAnalytics")
#' root <- unzip_file_system(root_pth, "fs")
#' association_file <- import_association_file(path, root)
import_association_file <- function(path,
    root, tp_padding = 4, dates_format = "dmy") {
    # Check parameters
    stopifnot(is.character(path) & length(path) == 1)
    stopifnot(is.character(root) & length(root) == 1)
    stopifnot(file.exists(path))
    stopifnot(file.exists(root))
    stopifnot((is.numeric(tp_padding) |
                   is.integer(tp_padding)) & length(tp_padding) == 1)
    stopifnot(length(dates_format) == 1 & dates_format %in% date_formats())

    # Read file and check the correctness
    as_file <- .read_and_correctness_af(path, tp_padding, dates_format)

    # Checks if the association file and the file system are aligned
    checks <- .check_file_system_alignment(as_file, root_folder = root)
    widg <- .checker_widget(checks)
    print(widg)

    as_file <- .update_af_after_alignment(as_file, checks, root)
    as_file
}


#' Import integration matrices based on the association file.
#'
#' @description These functions are designed to import the appropriate
#' integration matrix files given the association file and the root folder of
#' the file system where Vispa2 matrices are generated.
#' @details Import family functions are designed to work in combination with
#' Vispa2, for more details on this take a look here:
#' \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5702242/}{VISPA2:
#' A Scalable Pipeline for High-Throughput Identification
#' and Annotation of Vector Integration Sites}.\cr
#' For more details on how to properly use these functions, refer to the
#' vignette - vignette("how_to_import_functions")
#' @rdname import_parallel_Vispa2Matrices
#' @section Interactive version:
#' The interactive version of import_parallel_Vispa2Matrices asks user for input
#' and allows a more detailed choice of projects to import, pools to import and,
#' if necessary, duplicate files. During the execution, a series of reports is
#' shown in html format.
#' @param association_file A single string containing the path to the
#' association file on disk, or a tibble resulting from the previous call of
#' `import_association_file`
#' @param root A single string containing the path to the root folder containing
#' Vispa2 output. Can be NULL if association_file parameter is a tibble
#' @param quantification_type A vector of requested quantification_types. Must
#' be one in `quantification_types()`
#' @param matrix_type A single string representing the type of matrices to
#' be imported. Can only be one in `"annotated"` or `"not_annotated"`
#' @param workers A single number representing the number of parallel workers to
#' use for the import
#' @param tp_padding Timepoint padding, indicates the number of digits of the
#' "Timepoint" column once imported. Fills the content with 0s up to the length
#' specified (ex: 1 becomes 0001 with a tp_padding of 4)
#' @param dates_format A single string indicating how dates should be parsed.
#' Must be a value in: \code{date_formats()}
#' @importFrom htmltools browsable tagList
#' @importFrom tibble is_tibble
#' @family Import functions
#'
#' @return A named list of tibbles containing data from all imported integration
#' matrices, divided by quantification type
#' @export
#'
#' @examples
#' \dontrun{
#' # Can't run because it's interactive and requires user input
#' matrices <- import_parallel_Vispa2Matrices_interactive(
#'     association_file,
#'     root, quantification_type, matrix_type, workers
#' )
#' }
import_parallel_Vispa2Matrices_interactive <- function(association_file,
    root,
    quantification_type,
    matrix_type = "annotated",
    workers = 2,
    tp_padding = 4,
    dates_format = "dmy") {
    # Check parameters
    stopifnot(!missing(association_file))
    stopifnot(is.character(association_file) |
        tibble::is_tibble(association_file))
    if (is.character(association_file)) {
        stopifnot(length(association_file) == 1)
        stopifnot(!missing(root) && is.character(root) && length(root) == 1)
    } else {
        root <- NULL
    }
    stopifnot(is.numeric(workers) & length(workers) == 1)
    workers <- floor(workers)
    stopifnot(!missing(quantification_type))
    stopifnot(all(quantification_type %in% quantification_types()))
    stopifnot(is.character(matrix_type) & matrix_type %in% c(
        "annotated",
        "not_annotated"
    ))
    stopifnot((is.numeric(tp_padding) |
                   is.integer(tp_padding)) & length(tp_padding) == 1)
    stopifnot(length(dates_format) == 1 & dates_format %in% date_formats())

    # Manage association file
    association_file <- .manage_association_file(association_file, root,
                                                 tp_padding, dates_format)
    checker_widg <- association_file[[2]]
    association_file <- association_file[[1]]
    if (!is.null(checker_widg)) {
        print(checker_widg)
    }
    ## User selects projects to keep
    association_file <- .interactive_select_projects_import(association_file)
    ## User selects pools to keep
    association_file <- .interactive_select_pools_import(association_file)
    ## Scan the appropriate file system paths and look for files
    files_found <- .lookup_matrices(
        association_file, quantification_type,
        matrix_type
    )
    ff_widget <- .files_found_widget(files_found)
    if (!is.null(checker_widg)) {
        print(htmltools::browsable(htmltools::tagList(
            ff_widget,
            checker_widg
        )))
    } else {
        print(ff_widget)
    }
    ## Manage missing files and duplicates
    files_to_import <- .manage_anomalies_interactive(files_found)
    fimp_widget <- .files_to_import_widget(files_to_import)

    if (!is.null(checker_widg)) {
        print(htmltools::browsable(htmltools::tagList(
            fimp_widget,
            ff_widget,
            checker_widg
        )))
    } else {
        print(htmltools::browsable(htmltools::tagList(
            fimp_widget,
            ff_widget
        )))
    }

    ## If files to import are 0 just terminate
    if (nrow(files_to_import) == 0) {
        stop("No files to import")
    }

    ## Import
    matrices <- .parallel_import_merge(files_to_import, workers)
    fimported_widg <- matrices[[2]]
    fimported_widg <- .files_imported_widget(fimported_widg)
    if (!is.null(checker_widg)) {
        print(htmltools::browsable(htmltools::tagList(
            fimported_widg,
            fimp_widget,
            ff_widget,
            checker_widg
        )))
    } else {
        print(htmltools::browsable(htmltools::tagList(
            fimported_widg,
            fimp_widget,
            ff_widget
        )))
    }
    matrices <- matrices[[1]]
    matrices
}


#' @describeIn import_parallel_Vispa2Matrices Automatic version
#'
#' @section Automatic version:
#' The automatic version of import_parallel_Vispa2Matrices doesn't interact with
#' the user directly, for this reason options in this modality are more limited
#' compared to the interactive version. In automatic version you can't:
#' * Choose single projects or pools: to have a selection import the association
#' file first and filter it according to your needs before calling the function
#' (more details on this in the vignette)
#' * Choose duplicates: if, after filtering by the specified patterns,
#' duplicates are found they are automatically ignored
#'
#' @param patterns A character vector of additional patterns to match on file
#' names. Please note that patterns must be regular expressions. Can be NULL if
#' no patterns needs to be matched.
#' @param matching_opt A single value between `matching_options`
#' @seealso \code{\link{matching_options}},
#' \url{https://stringr.tidyverse.org/articles/regular-expressions.html}
#' @family Import functions
#' @importFrom htmltools browsable tagList
#' @importFrom tibble is_tibble
#'
#' @export
#'
#' @examples
#' path <- system.file("extdata", "ex_association_file.tsv",
#'     package = "ISAnalytics"
#' )
#' root_pth <- system.file("extdata", "fs.zip", package = "ISAnalytics")
#' root <- unzip_file_system(root_pth, "fs")
#' matrices <- import_parallel_Vispa2Matrices_auto(
#'     path, root,
#'     c("fragmentEstimate", "seqCount"), "annotated", 2, NULL, "ANY"
#' )
import_parallel_Vispa2Matrices_auto <- function(association_file,
    root,
    quantification_type,
    matrix_type = "annotated",
    workers = 2,
    patterns = NULL,
    matching_opt = matching_options(),
    tp_padding = 4,
    dates_format = "dmy") {
    # Check parameters
    stopifnot(!missing(association_file))
    stopifnot(is.character(association_file) |
        tibble::is_tibble(association_file))
    if (is.character(association_file)) {
        stopifnot(length(association_file) == 1)
        stopifnot(!missing(root) && is.character(root) && length(root) == 1)
    } else {
        root <- NULL
    }
    stopifnot(is.numeric(workers) & length(workers) == 1)
    workers <- floor(workers)
    stopifnot(!missing(quantification_type))
    stopifnot(all(quantification_type %in% quantification_types()))
    stopifnot(is.character(matrix_type) & matrix_type %in% c(
        "annotated",
        "not_annotated"
    ))
    if (!is.null(patterns)) {
        stopifnot(is.character(patterns))
    }
    stopifnot((is.numeric(tp_padding) |
                   is.integer(tp_padding)) & length(tp_padding) == 1)
    stopifnot(length(dates_format) == 1 & dates_format %in% date_formats())

    ### Evaluate matching_opt
    matching_option <- match.arg(matching_opt)
    stopifnot(is.character(matching_option))

    # Manage association file
    association_file <- .manage_association_file(association_file, root,
                                                 tp_padding, dates_format)
    checker_widg <- association_file[[2]]
    association_file <- association_file[[1]]
    if (!is.null(checker_widg)) {
        print(checker_widg)
    }

    # Automatic workflow - limited options
    ## In automatic workflow all projects and pools contained in the association
    ## file are considered. If more precise selection is needed on this, user
    ## should use the interactive version or filter the association file
    ## appropriately before calling the function.

    ## Scan the appropriate file system paths and look for files
    files_found <- .lookup_matrices_auto(
        association_file, quantification_type,
        matrix_type, patterns, matching_option
    )

    ff_widget <- .files_found_widget(files_found)
    if (!is.null(checker_widg)) {
        print(htmltools::browsable(htmltools::tagList(
            ff_widget,
            checker_widg
        )))
    } else {
        print(ff_widget)
    }

    ## Manage missing files and duplicates
    files_to_import <- .manage_anomalies_auto(files_found)
    fimp_widget <- .files_to_import_widget(files_to_import)

    if (!is.null(checker_widg)) {
        print(htmltools::browsable(htmltools::tagList(
            fimp_widget,
            ff_widget,
            checker_widg
        )))
    } else {
        print(htmltools::browsable(htmltools::tagList(
            fimp_widget,
            ff_widget
        )))
    }
    ## If files to import are 0 just terminate
    if (nrow(files_to_import) == 0) {
        stop("No files to import")
    }

    ## Import
    matrices <- .parallel_import_merge(files_to_import, workers)
    fimported_widg <- matrices[[2]]
    fimported_widg <- .files_imported_widget(fimported_widg)
    if (!is.null(checker_widg)) {
        print(htmltools::browsable(htmltools::tagList(
            fimported_widg,
            fimp_widget,
            ff_widget,
            checker_widg
        )))
    } else {
        print(htmltools::browsable(htmltools::tagList(
            fimported_widg,
            fimp_widget,
            ff_widget
        )))
    }
    matrices <- matrices[[1]]
    matrices
}

#' Possible choices for the `quantification_type` parameter.
#'
#' These are all the possible values for the `quantification_type` parameter in
#' `import_parallel_vispa2Matrices_interactive` and
#' `import_parallel_vispa2Matrices_auto`.
#'
#' @details The possible values are:
#' * fragmentEstimate
#' * seqCount
#' * barcodeCount
#' * cellCount
#' * ShsCount
#' @return A vector of characters for quantification types
#' @export
#' @seealso \code{
#' \link[import_parallel_vispa2Matrices]{
#' import_parallel_vispa2Matrices_interactive}},
#' \code{
#' \link[import_parallel_vispa2Matrices]{
#' import_parallel_vispa2Matrices_auto}}
#'
#' @family Import functions helpers
#'
#' @examples
#' quant_types <- quantification_types()
quantification_types <- function() {
    c(
        "fragmentEstimate", "seqCount",
        "barcodeCount", "cellCount",
        "ShsCount"
    )
}


#' Possible choices for the `matching_opt` parameter.
#'
#' These are all the possible values for the `matching_opt` parameter in
#' `import_parallel_vispa2Matrices_auto`.
#' @details The values "ANY", "ALL" and "OPTIONAL", represent how the patterns
#' should be matched, more specifically
#' * ANY = look only for files that match AT LEAST one of the patterns specified
#' * ALL = look only for files that match ALL of the patterns specified
#' * OPTIONAL = look preferentially for files that match, in order, all
#' patterns or any pattern and if no match is found return what is found (keep
#' in mind that duplicates are discarded in automatic mode)
#' @return A vector of characters for matching_opt
#' @export
#' @family Import functions helpers
#' @seealso \code{
#' \link[import_parallel_vispa2Matrices]{import_parallel_vispa2Matrices_auto}}
#'
#' @examples
#' opts <- matching_options()
matching_options <- function() {
    c("ANY", "ALL", "OPTIONAL")
}


#' Possible choices for the `dates_format` parameter in
#' `import_association_file`, `import_parallel_vispa2Matrices_interactive` and
#' `import_parallel_vispa2Matrices_auto`.
#'
#' All options correspond to `lubridate` functions:
#' * ymd: year, month, date
#' * ydm: year, day, month
#' * mdy: month, day, year
#' * myd: month, year, day
#' * dmy: day, month, year - default value
#' * dym: day, year, month
#' * yq: year quantile
#'
#' **NOTE: use the same date format across the association file.**
#'
#' @return A character vector
#' @export
#' @seealso \code{\link{import_association_file}}, \code{
#' \link[import_parallel_vispa2Matrices]{import_parallel_vispa2Matrices_auto}}
#'
#' @examples
#' date_formats()
date_formats <- function() {
    c("ymd", "ydm", "mdy", "myd", "dmy", "dym", "yq")
}
