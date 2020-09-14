#------------------------------------------------------------------------------#
# Utility functions
#------------------------------------------------------------------------------#

#' Creates a blank association file.
#'
#' This function is useful if you want a blank association file to start using
#' both Vispa2 and this package or simply if you want a correct framework to
#' fix a malformed association file you have already.
#'
#' @param path The path on disk where the file should be written
#' @importFrom fs path_dir as_fs_path dir_create dir_exists
#' @importFrom readr write_tsv
#'
#' @return returns NULL
#' @export
#' @family Utility functions
#'
#' @examples
#' temp <- tempfile()
#' generate_blank_association_file(temp)
generate_blank_association_file <- function(path) {
    stopifnot(is.character(path))
    af <- data.frame(matrix(
        ncol = length(association_file_columns()),
        nrow = 0
    ))
    colnames(af) <- association_file_columns()
    path <- fs::as_fs_path(path)
    dir <- fs::path_dir(path)
    if (!fs::dir_exists(dir)) {
        fs::dir_create(dir, recurse = TRUE)
    }
    readr::write_tsv(af, path = path)
}

#' Creates a reduced association file for Vispa2 run,
#' given project and pool
#'
#' The function selects the appropriate columns and prepares a file for the
#' launch of Vispa2 pipeline for each project/pool pair specified.
#'
#' @details Note: the function is vectorized, meaning you can specify more than
#' one project and more than one pool as vectors of characters, but you must
#' ensure that:
#' * Both `project` and `pool` vectors have the same length
#' * You correclty type names in corresponding positions, for example
#' c("CLOEXP", "PROJECT1100", "PROJECT1100") - c("POOL6", "ABX-LR-PL5-POOL14-1",
#' "ABX-LR-PL6-POOL15-1"). If you type a pool in the position of a corresponding
#' project that doesn't match no file will be produced since that pool doesn't
#' exist in the corresponding project.
#'
#' @param association_file The imported association file (via
#' `import_association_file`)
#' @param project A vector of characters containing project names
#' @param pool A vector of characters containing pool names.
#' **NOTE: the names should refer to the values contained in the
#' PoolID column of the association file and NOT the
#' concatenatePoolIDSeqRun column!**
#' @param path A single string representing the path to the folder where files
#' should be written. If the folder doesn't exist it will be created.
#' @importFrom fs as_fs_path file_exists dir_create path
#' @importFrom purrr map2 set_names walk2
#' @importFrom dplyr select filter mutate all_of
#' @importFrom readr write_tsv
#' @importFrom magrittr `%>%`
#' @family Utility functions
#'
#' @return returns NULL
#' @export
#'
#' @examples
#' op <- options("ISAnalytics.widgets" = FALSE)
#' temp <- tempdir()
#' path_af <- system.file("extdata", "ex_association_file.tsv",
#'     package = "ISAnalytics"
#' )
#' root_pth <- system.file("extdata", "fs.zip", package = "ISAnalytics")
#' root <- unzip_file_system(root_pth, "fs")
#' association_file <- import_association_file(path_af, root)
#' generate_Vispa2_launch_AF(association_file, "CLOEXP", "POOL6", temp)
#' options(op)
generate_Vispa2_launch_AF <- function(association_file, project, pool, path) {
    stopifnot(is.data.frame(association_file))
    stopifnot(is.character(project))
    stopifnot(is.character(pool))
    stopifnot(length(project) == length(pool))
    stopifnot(.check_af_correctness(association_file))
    stopifnot(is.character(path) & length(path) == 1)
    path <- fs::as_fs_path(path)
    if (!fs::file_exists(path)) {
        fs::dir_create(path)
    }
    files <- purrr::map2(project, pool, function(x, y) {
        selected_cols <- association_file %>%
            dplyr::select(dplyr::all_of(reduced_AF_columns())) %>%
            dplyr::filter(.data$ProjectID == x, .data$PoolID == y) %>%
            dplyr::mutate(TagID2 = .data$TagID, .before = .data$TagID)
    }) %>% purrr::set_names(paste0(project, "-", pool, "_AF.tsv"))
    purrr::walk2(files, names(files), function(x, y) {
        complete_path <- fs::path(path, y)
        if (nrow(x) > 0) {
            readr::write_tsv(x, complete_path, col_names = FALSE)
        } else {
            message(paste("Nothing to write for ", y, ", skipping."))
        }
    })
}


#' Converts tidy integration matrices in the original sparse matrix
#' form.
#'
#' \lifecycle{maturing}
#' This function is particularly useful when a sparce matrix structure
#' is needed by a specific function (mainly from other packages).
#'
#' @param x A single tidy integration matrix or a list of integration
#' matrices. Supports also multi-quantification matrices
#' obtained via \link{comparison_matrix}
#' @param fragmentEstimate For multi-quantification matrix support:
#' the name of the fragment estimate values column
#' @param seqCount For multi-quantification matrix support:
#' the name of the sequence count values column
#' @param barcodeCount For multi-quantification matrix support:
#' the name of the barcode count values column
#' @param cellCount For multi-quantification matrix support:
#' the name of the cell count values column
#' @param ShsCount For multi-quantification matrix support:
#' the name of the Shs Count values column
#'
#' @importFrom tidyr pivot_wider
#' @importFrom purrr is_empty set_names walk map
#' @importFrom dplyr select
#'
#' @family Utility functions
#'
#' @return Depending on input, 2 possible outputs:
#' * A single sparce matrix (tibble) if input is a single quantification
#' matrix
#' * A list of sparce matrices divided by quantification if input
#' is a single multi-quantification matrix or a list of matrices
#' @export
#'
#' @examples
#' path <- system.file("extdata", "ex_annotated_ISMatrix.tsv.xz",
#'     package = "ISAnalytics"
#' )
#' matrix <- import_single_Vispa2Matrix(path)
#' sparse <- as_sparse_matrix(matrix)
as_sparse_matrix <- function(x, fragmentEstimate = "fragmentEstimate",
    seqCount = "seqCount",
    barcodeCount = "barcodeCount",
    cellCount = "cellCount",
    ShsCount = "ShsCount") {
    stopifnot(is.list(x))
    if (is.data.frame(x)) {
        ## SINGLE DATA FRAME
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
        ### SINGLE QUANT
        if (all(num_cols == "Value")) {
            sparse_m <- tidyr::pivot_wider(x,
                names_from =
                    .data$CompleteAmplificationID,
                values_from =
                    .data$Value
            )
            return(sparse_m)
        }
        ### MULTI QUANT
        param_cols <- c(
            fragmentEstimate, seqCount, barcodeCount, cellCount,
            ShsCount
        )
        found <- param_cols[param_cols %in% num_cols]
        if (purrr::is_empty(found)) {
            stop(.non_quant_cols_error())
        }
        annot <- if (.is_annotated(x)) {
            annotation_IS_vars()
        } else {
            character(0)
        }
        sparse_m <- purrr::map(found, function(quant) {
            temp <- x %>% dplyr::select(
                mandatory_IS_vars(),
                .data$CompleteAmplificationID,
                annot, quant
            )
            tidyr::pivot_wider(temp,
                names_from = .data$CompleteAmplificationID,
                values_from = quant
            )
        }) %>% purrr::set_names(found)
        return(sparse_m)
    } else {
        ## LIST
        purrr::walk(x, function(m) {
            mand <- .check_mandatory_vars(m)
            amp <- .check_complAmpID(m)
            val <- .check_value_col(m)
            if (any(c(mand, amp, val) == FALSE)) {
                stop(.non_ISM_error())
            }
        })
        sparse_m <- purrr::map(x, function(data) {
            tidyr::pivot_wider(data,
                names_from = .data$CompleteAmplificationID,
                values_from = .data$Value
            )
        })
        return(sparse_m)
    }
}

#### ---- Utilities for tests and examples ----####
#' A utility function to unzip and use example file systems included in the
#' package
#'
#' This utility function is a simple shortcut to create a temporary directory,
#' unzip and reference the examples file systems included in the package for
#' testing purposes.
#'
#' @param zipfile The zipped file to decompress
#' @param name The name of the folder in the zipped archive ("fs" or "fserr")
#' @importFrom zip unzip
#' @family Utility functions
#' @export
#'
#' @return A path to reference
#' @examples
#' root_pth <- system.file("extdata", "fs.zip", package = "ISAnalytics")
#' root <- unzip_file_system(root_pth, "fs")
unzip_file_system <- function(zipfile, name) {
    root_folder <- tempdir()
    zip::unzip(zipfile, exdir = root_folder)
    root_folder <- file.path(root_folder, name, fsep = "\\")
    root_folder <- gsub('"', "", gsub("\\\\", "/", root_folder))
    root_folder
}
