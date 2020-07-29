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
#'
#' @examples
#' temp <- tempfile()
#' generate_blank_association_file(temp)
generate_blank_association_file <- function(path) {
  stopifnot(is.character(path))
  af <- data.frame(matrix(
    ncol = length(association_file_columns),
    nrow = 0)
    )
  colnames(af) <- association_file_columns
  path <- fs::as_fs_path(path)
  dir <- fs::path_dir(path)
  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir, recurse = TRUE)
  }
  readr::write_tsv(af, path = path)
}

#' Creates a reduced association file for Vispa2 run, given project and pool
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
#' @param pool A vector of characters containing pool names. **NOTE: the names
#' should refer to the values contained in the PoolID column of the association
#' file and NOT the concatenatePoolIDSeqRun column!**
#' @param path A single string representing the path to the folder where files
#' should be written. If the folder doesn't exist it will be created.
#' @importFrom fs as_fs_path file_exists dir_create path
#' @importFrom purrr map2 set_names walk2
#' @importFrom dplyr select filter mutate all_of
#' @importFrom readr write_tsv
#'
#' @return returns NULL
#' @export
#'
#' @examples
#' temp <- tempdir()
#' path_af <- system.file("extdata", "ex_association_file.tsv",
#' package = "ISAnalytics")
#' root_pth <- system.file("extdata", "fs.zip", package = "ISAnalytics")
#' root <- tempdir()
#' zip::unzip(root_pth, exdir = root)
#' root <- file.path(root, "fs")
#' root <- gsub('"', "", gsub("\\\\", "/", root))
#' association_file <- import_association_file(path_af, root)
#' generate_Vispa2_launch_AF(association_file, "CLOEXP", "POOL6", temp)
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
      dplyr::select(dplyr::all_of(reduced_AF_columns)) %>%
      dplyr::filter(.data$ProjectID == x, .data$PoolID == y) %>%
      dplyr::mutate(TagID2 = .data$TagID, .before = .data$TagID)
  }) %>% purrr::set_names(paste0(project,"-", pool,"_AF.tsv"))
  purrr::walk2(files, names(files), function(x, y) {
    complete_path <- fs::path(path, y)
    if (nrow(x) > 0) {
      readr::write_tsv(x, complete_path, col_names = FALSE)
    } else {
      message(paste("Nothing to write for ", y, ", skipping."))
    }
  })
}
