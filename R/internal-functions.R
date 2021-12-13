#------------------------------------------------------------------------------#
# Internal functions
#------------------------------------------------------------------------------#
## All functions in this file are NOT exported, to be used internally only.

#### ---- Internals for multiple functions/general purpose ----####
# Returns the file format for each of the file paths passed as a parameter.
#' @importFrom fs path_ext
#' @importFrom tools file_path_sans_ext
#' @importFrom purrr is_empty
.check_file_extension <- function(file_path) {
    ## Get last portion of file name
    last <- fs::path_ext(file_path)
    compressed <- which(last %in% .compressed_formats())
    if (!purrr::is_empty(compressed)) {
        ## for compressed files try extracting the true extension
        file_path[compressed] <- tools::file_path_sans_ext(
            file_path[compressed]
        )
        last <- fs::path_ext(file_path)
    }
    return(last)
}

#### ---- Internals for checks on integration matrices----####

# Internal helper function for checking mandatory vars presence in x.
#
# Checks if the elements of `mandatory_IS_vars` are present as column names
# in the data frame.
# @param x A data.frame object (or any extending class)
# @keywords internal
#
# @return FALSE if all or some elements are not found in the data frame, TRUE
# otherwise
.check_mandatory_vars <- function(x) {
    stopifnot(is.data.frame(x))
    res <- if (all(mandatory_IS_vars() %in% colnames(x))) {
        TRUE
    } else {
        FALSE
    }
    return(res)
}

# Internal helper function for checking `Value` column presence in x.
#
# Checks if the column `Value` is present in the data frame and also
# checks if the column is numeric or integer.
# @param x A data.frame object (or any extending class)
# @keywords internal
#
# @return FALSE if not found or contains non-numeric data, TRUE otherwise
.check_value_col <- function(x) {
    stopifnot(is.data.frame(x))
    present <- if ("Value" %in% colnames(x)) {
        TRUE
    } else {
        FALSE
    }
    if (present == TRUE) {
        return(is.numeric(x$Value) | is.integer(x$Value))
    } else {
        return(FALSE)
    }
}

# Internal helper function for checking `CompleteAmplifcationID`
# column presence in x.
#
# Checks if the column `CompleteAmplifcationID` is present in the data frame.
#
# @param x A data.frame object (or any extending class)
# @keywords internal
#
# @return FALSE if not found, TRUE otherwise
.check_complAmpID <- function(x) {
    stopifnot(is.data.frame(x))
    if ("CompleteAmplificationID" %in% colnames(x)) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

# Finds experimental columns in an integration matrix.
#
# The function checks if there are numeric columns which are not
# standard integration matrix columns, if there are returns their names.
#
# @param x A data.frame
#' @importFrom purrr map_lgl
#' @importFrom rlang expr eval_tidy
#' @keywords internal
#
# @return A character vector of column names
.find_exp_cols <- function(x, to_exclude) {
    stopifnot(is.data.frame(x))
    remaining <- colnames(x)[!colnames(x) %in% to_exclude]
    remaining_numeric <- purrr::map_lgl(remaining, function(y) {
        exp <- rlang::expr(`$`(x, !!y))
        exp <- rlang::eval_tidy(exp)
        is.numeric(rlang::eval_tidy(exp)) | is.integer(rlang::eval_tidy(exp))
    })
    remaining[remaining_numeric]
}

# Checks if the integration matrix is annotated or not.
#
# @param x A data.frame
# @keywords internal
#
# @return A logical value
.is_annotated <- function(x) {
    stopifnot(is.data.frame(x))
    if (all(annotation_IS_vars() %in% colnames(x))) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

#### ---- Internals for matrix import ----####

#---- USED IN : import_single_Vispa2Matrix ----

# Internal function to auto-detect the type of IS based on the headers.
#
# @param df the data frame to inspect
# @keywords internal
#
# @return one value among:
# * "OLD" : for old-style matrices that had only one column holding
# all genomic coordinates
# * "NEW" :  for the classic Vispa2 annotated/not annotated matrices
# * "MALFORMED" : in any other case
.auto_detect_type <- function(df) {
    if ("IS_genomicID" %in% colnames(df) &
        all(!mandatory_IS_vars() %in% colnames(df))) {
        return("OLD")
    }
    if (all(mandatory_IS_vars() %in% colnames(df)) &
        !"IS_genomicID" %in% colnames(df)) {
        return("NEW")
    }
    return("MALFORMED")
}

# Reads an integration matrix using data.table::fread
#' @importFrom data.table fread
.read_with_fread <- function(path, to_drop, df_type, annotated, sep) {
    df <- if (df_type == "OLD") {
        data.table::fread(
            file = path,
            sep = sep,
            na.strings = c("NONE", "NA", "NULL", "NaN", ""),
            verbose = FALSE,
            drop = to_drop,
            colClasses = list(
                character = "IS_genomicID"
            ),
            showProgress = getOption("ISAnalytics.verbose"),
            data.table = TRUE
        )
    } else {
        col_types <- .mandatory_IS_types("fread")
        if (annotated) {
            col_types$character <- append(
                col_types$character,
                .annotation_IS_types("fread")$character
            )
        }
        data.table::fread(
            file = path,
            sep = sep,
            na.strings = c("NONE", "NA", "NULL", "NaN", ""),
            verbose = FALSE,
            drop = to_drop,
            colClasses = col_types,
            showProgress = getOption("ISAnalytics.verbose"),
            data.table = TRUE
        )
    }
    return(df)
}

# Reads an integration matrix using readr::read_delim
#' @importFrom readr read_delim cols
#' @importFrom data.table setDT
.read_with_readr <- function(path, to_drop, df_type, annotated, sep) {
    col_types <- if (df_type == "NEW") {
        .mandatory_IS_types("classic")
    } else {
        list(IS_genomicID = "c")
    }
    if (annotated) {
        col_types <- append(
            col_types,
            .annotation_IS_types("classic")
        )
    }
    if (!is.null(to_drop)) {
        for (x in to_drop) {
            col_types[[x]] <- "_"
        }
    }
    col_types[[".default"]] <- "n"
    df <- readr::read_delim(
        file = path,
        delim = sep,
        col_types = do.call(readr::cols, col_types),
        na = c("NONE", "NA", "NULL", "NaN", ""),
        trim_ws = TRUE,
        progress = getOption("ISAnalytics.verbose")
    )
    df <- data.table::setDT(df)
    return(df)
}

#---- USED IN : import_association_file ----

# Checks if the association file has the right format (standard headers).
#
# @param df The imported association file
# @keywords internal
#
# @return TRUE if the check passes, FALSE otherwise
.check_af_correctness <- function(df) {
    if (all(association_file_columns() %in% colnames(df))) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}


# Imports association file from disk. Converts dates and pads timepoints,
# reporting parsing problems.
#' @importFrom rlang inform
#' @importFrom readr read_delim cols problems
#' @importFrom purrr map2_lgl is_empty map_chr map_dfr
#' @importFrom lubridate parse_date_time
#' @importFrom tibble tibble
#' @importFrom dplyr mutate across select all_of
.read_af <- function(path, padding, date_format, delimiter) {
    mode <- "readr"
    ## - Check file extension
    file_ext <- .check_file_extension(path)
    if (file_ext %in% c("xls", "xlsx")) {
        mode <- "readxl"
        if (!requireNamespace("readxl", quietly = TRUE)) {
            rlang::abort(.missing_pkg_error("readxl"))
        }
        if (getOption("ISAnalytics.verbose") == TRUE) {
            rlang::inform(.xls_file_warning(),
                class = "xls_file"
            )
        }
    }
    ## - Peek headers to know column types
    col_types <- if (mode != "readxl") {
        headers_peek <- readr::read_delim(path,
            delim = delimiter,
            n_max = 0,
            col_types = readr::cols()
        )
        t_readr <- .af_col_types("readr")
        t_readr[names(t_readr) %in% colnames(headers_peek)]
    } else {
        headers_peek <- readxl::read_excel(path, n_max = 0)
        t_readr <- .af_col_types("readr")
        ordered <- purrr::map_chr(
            colnames(headers_peek),
            function(x) {
                el <- getElement(t_readr, x)
                el <- if (is.null(el)) {
                    "guess"
                } else if (el == "c") {
                    "text"
                } else if (el %in% c("i", "d")) {
                    "numeric"
                } else {
                    "guess"
                }
            }
        )
        ordered
    }
    suppressWarnings({
        as_file <- if (mode == "readr") {
            df <- readr::read_delim(path,
                delim = delimiter,
                na = c("NONE", "NA", "NULL", "NaN", ""),
                col_types = do.call(readr::cols, col_types),
                trim_ws = TRUE,
                progress = getOption("ISAnalytics.verbose")
            )
            problems <- readr::problems(df)
            df
        } else {
            df <- readxl::read_excel(path,
                col_types = col_types,
                na = c("NONE", "NA", "NULL", "NaN", ""),
                trim_ws = TRUE,
                progress = getOption("ISAnalytics.verbose")
            )
            problems <- NULL
            df
        }
        if ("TimePoint" %in% colnames(as_file)) {
            as_file <- as_file %>%
                dplyr::mutate(TimePoint = stringr::str_pad(
                    as.character(.data$TimePoint),
                    padding,
                    side = "left",
                    pad = "0"
                ))
        }
        date_cols <- colnames(as_file)[stringr::str_detect(
            colnames(as_file), "Date"
        )]
        before <- as_file %>% dplyr::select(dplyr::all_of(date_cols))
        as_file <- as_file %>%
            dplyr::mutate(dplyr::across(date_cols,
                .fns = ~ lubridate::as_date(lubridate::parse_date_time(.x,
                    orders = date_format
                ))
            ))
        date_failures <- purrr::map_dfr(date_cols, function(col) {
            before_col <- purrr::pluck(before, col)
            row_failed <- which(
                purrr::map2_lgl(before_col, as_file[[col]], function(d1, d2) {
                    if (!is.na(d1) & is.na(d2)) {
                        TRUE
                    } else {
                        FALSE
                    }
                })
            )
            if (!is.null(row_failed) && !purrr::is_empty(row_failed)) {
                tibble::tibble(
                    row = row_failed, col = col,
                    original = before_col[row_failed],
                    parsed = NA, format = date_format
                )
            } else {
                return(NULL)
            }
        })
    })
    return(list(af = as_file, probs = problems, date_fail = date_failures))
}

# Internal function to check alignment between association file and file
# system starting from the root. The alignment is checked at folder level.
#
# @param df The imported association file (data.frame or tibble)
# @param root_folder Path to the root folder
# @keywords internal
#' @importFrom rlang .data `:=`
#
# @return A data frame containing, for each ProjectID and
# concatenatePoolIDSeqRun the
# corresponding path on disk if found, NA otherwise.
.check_file_system_alignment <- function(df, root_folder) {
    if (!"PathToFolderProjectID" %in% colnames(df)) {
        rlang::abort(.af_missing_pathfolder_error())
    }
    temp_df <- df %>%
        dplyr::select(
            .data$ProjectID,
            .data$concatenatePoolIDSeqRun,
            .data$PathToFolderProjectID
        ) %>%
        dplyr::distinct()
    path_cols <- .path_cols_names()
    proj_folders_exist <- temp_df %>%
        dplyr::select(.data$PathToFolderProjectID) %>%
        dplyr::distinct() %>%
        dplyr::mutate(Found := !is.na(
            .data$PathToFolderProjectID
        ) &
            unname(fs::dir_exists(
                fs::path(
                    fs::path(root_folder),
                    .data$PathToFolderProjectID
                )
            )))
    partial_check <- temp_df %>%
        dplyr::left_join(proj_folders_exist, by = "PathToFolderProjectID")
    temp_df <- partial_check %>%
        dplyr::filter(Found == TRUE)
    partial_check <- partial_check %>%
        dplyr::filter(Found == FALSE) %>%
        dplyr::mutate(
            !!path_cols$project := NA_character_,
            !!path_cols$quant := NA_character_,
            !!path_cols$iss := NA_character_
        )
    FUN <- function(...) {
        cur <- tibble::tibble(...)
        project_folder <- fs::path(
            fs::path(root_folder),
            cur$PathToFolderProjectID
        )
        quant_folder <- if (!is.na(cur$concatenatePoolIDSeqRun)) {
            paste0(fs::path(
                "quantification",
                fs::path(cur$concatenatePoolIDSeqRun)
            ), "$")
        } else {
            if (getOption("ISAnalytics.verbose") == TRUE) {
                rlang::inform(c(paste(
                    "Warning: found NA",
                    "concatenatePoolIDSeqRun field"
                ),
                i = "Check association file for possible issues"
                ))
            }
            NA_character_
        }
        iss_folder <- if (!is.na(cur$concatenatePoolIDSeqRun)) {
            paste0(fs::path(
                "iss",
                fs::path(cur$concatenatePoolIDSeqRun)
            ), "$")
        } else {
            NA_character_
        }
        quant_found <- unique(fs::dir_ls(
            path = project_folder, recurse = TRUE,
            type = "directory", fail = FALSE,
            regexp = quant_folder
        ))
        if (length(quant_found) == 0 || all(is.na(quant_found))) {
            quant_found <- NA_character_
        }
        iss_found <- unique(fs::dir_ls(
            path = project_folder, recurse = TRUE,
            type = "directory", fail = FALSE,
            regexp = iss_folder
        ))
        if (length(iss_found) == 0 || all(is.na(iss_found))) {
            iss_found <- NA_character_
        }
        return(
            cur %>%
                dplyr::mutate(
                    !!path_cols$project := project_folder,
                    !!path_cols$quant := quant_found,
                    !!path_cols$iss := iss_found
                )
        )
    }
    results_df <- purrr::pmap_dfr(
        temp_df,
        FUN
    )
    checker_df <- results_df %>%
        dplyr::bind_rows(partial_check)
    checker_df
}

# Updates the association file after the alignment check.
#
# The function updates the association
# file by adding a column `Path` where the absolute path on disk for the
# project and pool is found, if no path is found NA is inserted instead.
#
# @param as_file The tibble representing the read association_file
# @param checks The tibble representing the results
# of `.check_file_system_alignment`
# @param root The root folder
# @keywords internal
#
# @return An updated association file with absolute paths
#' @importFrom dplyr left_join select
.update_af_after_alignment <- function(as_file, checks, root) {
    as_file <- as_file %>%
        dplyr::left_join(checks %>%
            dplyr::select(
                -.data$PathToFolderProjectID,
                -.data$Found
            ),
        by = c("ProjectID", "concatenatePoolIDSeqRun")
        )
    as_file
}

#---- USED IN : import_Vispa2_stats ----

# Finds automatically the path on disk to each stats file.
#' @importFrom purrr pmap_dfr detect_index
#' @importFrom tibble tibble
#' @importFrom fs dir_ls
#' @importFrom rlang .data
# @keywords internal
# @return A tibble with columns: ProjectID, concatenatePoolIDSeqRun,
# Path_iss (or designated dynamic name), stats_files, info
.stats_report <- function(association_file, prefixes) {
    path_col_names <- .path_cols_names()
    temp <- association_file %>%
        dplyr::select(
            .data$ProjectID,
            .data$concatenatePoolIDSeqRun,
            .data[[path_col_names$iss]]
        ) %>%
        dplyr::distinct() %>%
        dplyr::mutate(
            ProjectID = as.factor(.data$ProjectID),
            concatenatePoolIDSeqRun = as.factor(
                .data$concatenatePoolIDSeqRun
            )
        )
    # If paths are all NA return
    if (all(is.na(temp[[path_col_names$iss]]))) {
        return(temp %>% dplyr::mutate(
            stats_files = NA_character_,
            info = list("NO FOLDER FOUND")
        ))
    }
    match_pattern <- function(pattern, temp_row) {
        if (is.na(temp_row[[path_col_names$iss]])) {
            return(tibble::tibble(pattern = pattern, file = NA_character_))
        }
        # For each prefix pattern search in the iss folder
        # Note: there can be
        # - a single file matching the prefix (ideal)
        # - multiple files matching
        # - no file matching
        files <- fs::dir_ls(temp_row[[path_col_names$iss]],
            type = "file", fail = FALSE,
            regexp = pattern
        )
        if (length(files) == 0) {
            files <- NA_character_
        }
        tibble::tibble(pattern = pattern, file = files)
    }
    stats_paths <- purrr::pmap_dfr(temp, function(...) {
        temp_row <- tibble::tibble(...)
        matches <- purrr::map_dfr(prefixes, ~ match_pattern(.x, temp_row))
        if (all(is.na(matches$file))) {
            # No stats files found for all prefixes
            return(temp_row %>%
                dplyr::mutate(
                    stats_files = NA_character_,
                    info = list("NOT FOUND")
                ))
        }
        if (length(which(!is.na(matches$file))) == 1) {
            return(temp_row %>%
                dplyr::mutate(
                    stats_files = (matches %>%
                        dplyr::filter(!is.na(.data$file)) %>%
                        dplyr::pull(.data$file)),
                    info = list("NULL")
                ))
        }
        # Get the count of files found for each pattern and order them
        # as preference - preferred pattern is the first in the prefixes
        # parameter
        pattern_counts <- matches %>%
            dplyr::mutate(pattern = factor(.data$pattern,
                levels = prefixes
            )) %>%
            dplyr::filter(!is.na(.data$file)) %>%
            dplyr::group_by(.data$pattern) %>%
            dplyr::tally() %>%
            dplyr::arrange(.data$pattern)
        one_index <- purrr::detect_index(pattern_counts$n, ~ .x == 1)
        if (one_index == 0) {
            ## Means there are duplicates for all patterns, gets reported
            ## but no files will be imported
            return(temp_row %>%
                dplyr::mutate(
                    stats_files = NA_character_,
                    info = list(matches %>%
                        dplyr::filter(!is.na(.data$file)) %>%
                        dplyr::mutate(info = "DUPLICATE"))
                ))
        } else {
            ## Pick the file with the pattern that has a count of 1 (no
            ## duplicates)
            chosen_pattern <- pattern_counts$pattern[one_index]
            return(temp_row %>%
                dplyr::mutate(
                    stats_files = matches %>%
                        dplyr::filter(.data$pattern == chosen_pattern) %>%
                        dplyr::pull(.data$file),
                    info = list("NULL")
                ))
        }
    })
    stats_paths
}

# Imports all found Vispa2 stats files.
#
#' @importFrom purrr map_chr pmap_dfr reduce is_empty
#' @importFrom tibble tibble
#' @importFrom data.table fread
# @keywords internal
#
# @return A list with the imported stats and a report of imported files. If no
# files were imported returns NULL instead
.import_stats_iss <- function(association_file, prefixes) {
    # Obtain paths
    stats_paths <- .stats_report(association_file, prefixes)
    stats_paths <- stats_paths %>%
        tidyr::unnest(.data$info)
    if (all(is.na(stats_paths$stats_files))) {
        stats_paths <- stats_paths %>%
            dplyr::mutate(
                Imported = FALSE,
                reason = as.factor(.data$info)
            ) %>%
            dplyr::mutate(-.data$info)
        return(list(stats = NULL, report = stats_paths))
    }
    # Setup parallel workers and import
    # Register backend according to platform
    if (.Platform$OS.type == "windows") {
        p <- BiocParallel::SnowParam(
            stop.on.error = FALSE,
            tasks = length(stats_paths$stats_files),
            progressbar = getOption("ISAnalytics.verbose"),
            exportglobals = TRUE
        )
    } else {
        p <- BiocParallel::MulticoreParam(
            stop.on.error = FALSE,
            tasks = length(stats_paths$stats_files),
            progressbar = getOption("ISAnalytics.verbose"),
            exportglobals = FALSE
        )
    }
    FUN <- function(x) {
        if (is.na(x)) {
            return(NULL)
        }
        stats <- data.table::fread(
            file = x, sep = "\t",
            na.strings = c("", "NA", "na", "NONE"),
            data.table = TRUE
        )
        ok <- .check_stats(stats)
        if (ok == TRUE) {
            return(stats %>%
                dplyr::mutate(TAG = stringr::str_replace_all(
                    .data$TAG,
                    pattern = "\\.", replacement = ""
                )))
        } else {
            return("MALFORMED")
        }
    }
    stats_dfs <- BiocParallel::bptry(
        BiocParallel::bplapply(stats_paths$stats_files,
            FUN,
            BPPARAM = p
        )
    )
    BiocParallel::bpstop(p)
    correct <- purrr::map_chr(stats_dfs, function(x) {
        if (all(is.data.frame(x))) {
            "TRUE"
        } else if (is.null(x)) {
            "FALSE"
        } else {
            x
        }
    })
    stats_paths <- stats_paths %>%
        dplyr::mutate(Imported = correct)
    stats_paths <- purrr::pmap_dfr(stats_paths, function(...) {
        row <- tibble::tibble(...)
        condition <- if (row$Imported == "MALFORMED") {
            "MALFORMED"
        } else if (row$Imported == "TRUE") {
            NA_character_
        } else {
            row$info
        }
        row %>%
            dplyr::mutate(reason = as.factor(condition)) %>%
            dplyr::mutate(Imported = dplyr::if_else(
                condition = (.data$Imported == "MALFORMED"),
                true = FALSE,
                false = as.logical(.data$Imported)
            ))
    })
    stats_paths <- stats_paths %>%
        dplyr::select(-.data$info)
    stats_dfs <- stats_dfs[stats_paths$Imported]
    # Bind rows in single tibble for all files
    if (purrr::is_empty(stats_dfs)) {
        return(list(stats = NULL, report = stats_paths))
    }
    stats_dfs <- purrr::reduce(stats_dfs, function(x, y) {
        x %>%
            dplyr::bind_rows(y) %>%
            dplyr::distinct()
    })
    list(stats = stats_dfs, report = stats_paths)
}

#---- USED IN : import_parallel_Vispa2Matrices_interactive ----
# Base arg check (valid for both versions)
.base_param_check <- function(association_file,
    quantification_type,
    matrix_type,
    workers,
    multi_quant_matrix) {
    # Check parameters
    stopifnot((is.character(association_file) &
        length(association_file) == 1) ||
        is.data.frame(association_file))
    stopifnot(is.numeric(workers) & length(workers) == 1)
    stopifnot(!missing(quantification_type))
    stopifnot(all(quantification_type %in% quantification_types()))
    stopifnot(is.character(matrix_type) & matrix_type %in% c(
        "annotated",
        "not_annotated"
    ))
    stopifnot(is.logical(multi_quant_matrix) & length(multi_quant_matrix) == 1)
}

.pre_manage_af <- function(association_file, import_af_args) {
    ## Import association file if provided a path
    if (is.character(association_file)) {
        association_file <- rlang::eval_tidy(
            rlang::call2("import_association_file",
                path = association_file,
                !!!import_af_args
            )
        )
    }
    ## Check there are the appropriate columns
    if (!.path_cols_names()$quant %in% colnames(association_file)) {
        rlang::abort(.af_missing_path_error(.path_cols_names()$quant),
            class = "missing_path_col"
        )
    }
    association_file <- association_file %>%
        dplyr::filter(!is.na(.data[[.path_cols_names()$quant]]))
    return(association_file)
}

# Helper function to be used internally to treat association file.
#' @importFrom rlang inform enexpr eval_tidy parse_expr
#' @importFrom purrr is_empty map2_chr
.manage_association_file <- function(association_file,
    root, padding, format,
    delimiter, filter) {
    # If it's a path to file import the association file
    association_file <- .read_af(
        association_file,
        padding, format,
        delimiter
    )
    parsing_probs <- association_file$probs
    date_probs <- association_file$date_fail
    association_file <- association_file$af
    if (!is.null(filter)) {
        # Pre-filtering of association file
        if (!all(names(filter) %in% colnames(association_file))) {
            if (getOption("ISAnalytics.verbose") == TRUE) {
                rlang::inform(
                    c("Some or all the names in the filter not found",
                        i = paste(
                            "Columns not found: ",
                            paste0(
                                names(filter)[names(filter) %in%
                                    colnames(association_file)],
                                collapse = ", "
                            )
                        ),
                        "Ignoring the missing columns"
                    ),
                    class = "filter_warn"
                )
            }
            filter <- filter[names(filter) %in% colnames(association_file)]
        }
        if (!purrr::is_empty(filter)) {
            predicate <- purrr::map2_chr(
                filter, names(filter),
                function(x, y) {
                    if (is.character(x)) {
                        paste0(
                            y, " %in% c(", paste0("'",
                                rlang::enexpr(x),
                                "'",
                                collapse = ", "
                            ),
                            ")"
                        )
                    } else {
                        paste0(
                            y, " %in% c(", paste0(rlang::enexpr(x),
                                collapse = ", "
                            ),
                            ")"
                        )
                    }
                }
            )
            predicate <- paste0(c(
                "dplyr::filter(association_file, ",
                paste0(predicate, collapse = ", "),
                ")"
            ), collapse = "")
            association_file <- rlang::eval_tidy(
                rlang::parse_expr(predicate)
            )
        }
    }
    checks <- NULL
    if (!is.null(root) & nrow(association_file) > 0) {
        checks <- .check_file_system_alignment(association_file, root)
        association_file <- .update_af_after_alignment(
            association_file,
            checks, root
        )
    }
    res <- list(
        af = association_file, check = checks,
        parsing_probs = parsing_probs, date_probs = date_probs
    )
    return(res)
}

# Allows the user to choose interactively the projects
# to consider for import.
#
# @param association_file The tibble representing the imported association file
# @keywords internal
#' @importFrom dplyr distinct select filter
#' @importFrom rlang .data
#' @importFrom stringr str_split
#' @importFrom purrr map_dbl
#
# @return A modified version of the association file where only selected
# projects are present
.interactive_select_projects_import <- function(association_file) {
    repeat {
        cat("Which projects would you like to import?\n")
        cat("[1] ALL", "\n", "[2] ONLY SOME\n", "[0] QUIT\n", sep = "")
        # Keep asking user until the choice is valid
        repeat {
            cat("Your choice: ")
            connection <- getOption("ISAnalytics.connection")
            if (length(connection) > 1) {
                con <- connection[1]
            } else {
                con <- connection
            }
            n_projects_to_import <- readLines(con = con, n = 1)
            n_projects_to_import <- as.numeric(n_projects_to_import)
            if (!is.na(n_projects_to_import) &
                is.numeric(n_projects_to_import) &
                n_projects_to_import %in% c(0, 1, 2)) {
                if (n_projects_to_import == 0) {
                    stop("Quitting", call. = FALSE)
                } else {
                    break
                }
            } else {
                message("Invalid choice, please choose between the above.")
            }
        }
        # If only some projects selected
        if (n_projects_to_import == 2) {
            cat(
                "Here is the list of available projects, type the indexes",
                "separated by a comma:\n",
                sep = ""
            )
            project_list <- dplyr::distinct(
                dplyr::select(association_file, .data$ProjectID)
            )$ProjectID
            cat(paste0("[", seq_along(project_list), "] ", project_list),
                "[0] QUIT\n",
                sep = "\n"
            )
            # Keep asking user until a valid choice
            repeat {
                cat("Your choice: ")
                if (length(connection) > 1) {
                    con <- connection[2]
                } else {
                    con <- connection
                }
                projects_to_import <- readLines(con = con, n = 1)
                projects_to_import <- purrr::map_dbl(unlist(
                    stringr::str_split(projects_to_import, ",")
                ), as.numeric)
                if (!any(is.na(projects_to_import)) &
                    all(is.numeric(projects_to_import)) &
                    all(projects_to_import %in% c(
                        seq_along(project_list),
                        0
                    ))) {
                    if (any(projects_to_import == 0)) {
                        stop("Quitting", call. = FALSE)
                    } else {
                        break
                    }
                } else {
                    message(paste(
                        "Invalid choice, please select valid indexes separated",
                        "by a comma (example: 1, 2, 3) "
                    ))
                }
            }
            # Ask confirm
            cat("\nYour choices: ",
                paste(project_list[projects_to_import], collapse = ","),
                "\n",
                "Confirm your choices? [y/n]",
                sep = ""
            )
            repeat {
                if (length(connection) > 1) {
                    con <- connection[3]
                } else {
                    con <- connection
                }
                confirm <- readLines(con = con, n = 1)
                if (!confirm %in% c("y", "n")) {
                    message(paste(
                        "Invalid choice, please select one between y",
                        "or n"
                    ))
                } else {
                    break
                }
            }
            if (confirm == "y") {
                result <- association_file %>%
                    dplyr::filter(.data$ProjectID %in%
                        project_list[projects_to_import])
                return(result)
            } else {
                next
            }
        } else {
            # Ask confirm
            cat("\nYour choices: ", "import all projects", "\n",
                "Confirm your choices? [y/n]",
                sep = ""
            )
            repeat {
                if (length(connection) > 1) {
                    con <- connection[3]
                } else {
                    con <- connection
                }
                confirm <- readLines(con = con, n = 1)
                if (!confirm %in% c("y", "n")) {
                    message(paste(
                        "Invalid choice, please select one between y",
                        "or n"
                    ))
                } else {
                    break
                }
            }
            if (confirm == "y") {
                return(association_file)
            } else {
                next
            }
        }
    }
}

# Simple internal helper function to handle user input for selection
# of number of pools.
# @keywords internal
#
# @return Numeric representing user selection (1 for all pools, 2 for
# only some pools, 0 to exit)
.pool_number_IN <- function() {
    cat("Which pools for each project would you like to import?\n")
    cat("[1] ALL", "[2] ONLY SOME", "[0] QUIT", sep = "\n")
    repeat {
        cat("Your choice: ")
        connection <- getOption("ISAnalytics.connection")
        if (length(connection) > 1) {
            connection <- connection[4]
        }
        n_pools_to_import <- readLines(con = connection, n = 1)
        n_pools_to_import <- as.numeric(n_pools_to_import)
        if (!is.na(n_pools_to_import) & is.numeric(n_pools_to_import) &
            n_pools_to_import %in% c(0, 1, 2)) {
            if (n_pools_to_import == 0) {
                stop("Quitting", call. = FALSE)
            } else {
                break
            }
        } else {
            message("Invalid choice, please choose between the above.")
        }
    }
    n_pools_to_import
}

# Simple helper interal function to handle user input for actual
# pool choices.
#
# @param indexes A vector of integer indexes available
# @keywords internal
#' @importFrom stringr str_split
#' @importFrom purrr map_dbl
#
# @return The user selection as a numeric vector
.pool_choices_IN <- function(indexes) {
    repeat {
        cat("\nYour choice: ")
        connection <- getOption("ISAnalytics.connection")
        if (length(connection) > 1) {
            connection <- connection[5]
        }
        to_imp <- readLines(
            con = connection,
            n = 1
        )
        to_imp <- purrr::map_dbl(unlist(
            stringr::str_split(to_imp, ",")
        ), as.numeric)
        if (!any(is.na(to_imp)) & all(is.numeric(to_imp))) {
            if (all(to_imp %in%
                c(indexes, 0))) {
                if (any(to_imp == 0)) {
                    stop("Quitting", call. = FALSE)
                } else {
                    break
                }
            } else {
                message(paste(
                    "Invalid choice, please select valid indexes",
                    "separated by a comma (example: 1, 2, 3) "
                ))
                next
            }
        } else {
            message(paste(
                "Invalid choice, please select valid indexes separated",
                "by a comma (example: 1, 2, 3) "
            ))
        }
    }
    to_imp
}

# Allows the user to choose interactively
# the pools to consider for import.
#
# @param association_file The tibble representing the imported
# association file
#' @importFrom dplyr select distinct group_by bind_rows inner_join
#' @importFrom tibble tibble
#' @importFrom tidyr nest
#' @importFrom purrr map pmap reduce
#' @importFrom rlang .data
# @keywords internal
#
# @return A modified version of the association file where only selected
# pools for each project are present
.interactive_select_pools_import <- function(association_file) {
    repeat {
        n_pools_to_import <- .pool_number_IN()
        if (n_pools_to_import == 2) {
            cat("Here is the list of available pools for each project,",
                "type the indexes separated by a comma or 0 to quit: \n\n",
                sep = ""
            )
            available <- association_file %>%
                dplyr::select(
                    .data$ProjectID,
                    .data$concatenatePoolIDSeqRun
                ) %>%
                dplyr::distinct() %>%
                tidyr::nest(data = c(.data$concatenatePoolIDSeqRun))
            pools_to_import <- purrr::pmap(available, function(...) {
                l <- list(...)
                current <- l["ProjectID"]
                current_data <- tibble::as_tibble(
                    purrr::flatten(l["data"])
                )
                indexes <- seq_along(current_data$concatenatePoolIDSeqRun)
                cat("ProjectID: ", current$ProjectID, "\n\n")
                cat("[0] QUIT\n")
                purrr::walk2(
                    current_data$concatenatePoolIDSeqRun,
                    indexes,
                    function(x, y) {
                        cat(paste0("[", y, "] ", x),
                            sep = "\n"
                        )
                    }
                )
                to_imp <- .pool_choices_IN(indexes)
                tibble::tibble(
                    ProjectID = current$ProjectID,
                    concatenatePoolIDSeqRun =
                        current_data$concatenatePoolIDSeqRun[to_imp]
                )
            })
            pools_to_import <- purrr::reduce(pools_to_import, function(x, y) {
                dplyr::bind_rows(x, y)
            })
            cat("\nYour choices: ", sep = "\n")
            print(pools_to_import)
            cat("\nConfirm your choices? [y/n]", sep = "")
            repeat {
                connection <- getOption("ISAnalytics.connection")
                if (length(connection) > 1) {
                    connection <- connection[6]
                }
                confirm <- readLines(con = connection, n = 1)
                if (!confirm %in% c("y", "n")) {
                    message(paste(
                        "Invalid choice, please select one between y",
                        "or n"
                    ))
                } else {
                    break
                }
            }
            if (confirm == "y") {
                association_file <- association_file %>%
                    dplyr::inner_join(pools_to_import,
                        by = c("ProjectID", "concatenatePoolIDSeqRun")
                    )
                return(association_file)
            } else {
                next
            }
        } else {
            cat("\nYour choices: ", sep = "\n")
            print(association_file %>%
                dplyr::select(
                    .data$ProjectID,
                    .data$concatenatePoolIDSeqRun
                ) %>%
                dplyr::distinct())
            cat("\nConfirm your choices? [y/n]", sep = "")
            repeat {
                connection <- getOption("ISAnalytics.connection")
                if (length(connection) > 1) {
                    connection <- connection[6]
                }
                confirm <- readLines(con = connection, n = 1)
                if (!confirm %in% c("y", "n")) {
                    message(paste(
                        "Invalid choice, please select one between y",
                        "or n"
                    ))
                } else {
                    break
                }
            }
            if (confirm == "y") {
                return(association_file)
            } else {
                next
            }
        }
    }
}

# Updates a files_found tibble with Files_count and Anomalies column.
#
# @param lups A files_found tibble obtained in a lookup function. Must contain
# the Files column (nested table quantification type and files)
# @keywords internal
#
# @return Updated files_found with Anomalies and Files_count columns
.trace_anomalies <- function(lups) {
    files_count <- purrr::pmap(lups, function(...) {
        temp <- tibble::tibble(...)
        an <- temp$Files
        an <- an %>%
            dplyr::group_by(.data$Quantification_type) %>%
            dplyr::summarise(
                Found = sum(!is.na(.data$Files_found)),
                .groups = "drop_last"
            )
        an
    })

    lups <- lups %>%
        dplyr::mutate(Files_count = files_count) %>%
        dplyr::mutate(
            Anomalies = purrr::map_lgl(
                .data$Files_count,
                ~ any(.x["Found"] != 1)
            ),
            .before = .data$Files
        )
    lups
}

# Looks up matrices to import given the association file and the
# root of the file system.
#
# @param association_file Tibble representing the association file
# @param quantification_type The type of quantification matrices to look for
# (one in `quantification_types()`)
# @param matrix_type The matrix_type to lookup (one between "annotated" or
# "not_annotated")
# @keywords internal
#' @importFrom tibble tibble
#' @importFrom fs dir_ls as_fs_path
#' @importFrom purrr pmap_dfr cross_df
#' @importFrom stringr str_replace_all
#' @importFrom tidyr nest
#
# @return A tibble containing all found files, including duplicates and missing
.lookup_matrices <- function(association_file,
    quantification_type,
    matrix_type) {
    path_col_names <- .path_cols_names()
    temp <- association_file %>%
        dplyr::select(
            .data$ProjectID,
            .data$concatenatePoolIDSeqRun,
            .data[[path_col_names$quant]]
        ) %>%
        dplyr::distinct()
    ## Obtain a df with all possible combination of suffixes for each
    ## quantification
    ms <- if (matrix_type == "annotated") {
        .matrix_annotated_suffixes()
    } else {
        .matrix_not_annotated_suffixes()
    }
    cross <- purrr::cross_df(list(
        quant = quantification_type,
        ms = ms
    ))
    cross <- cross %>%
        dplyr::mutate(suffix = paste0(
            .data$quant,
            "_matrix",
            .data$ms,
            ".tsv"
        )) %>%
        dplyr::mutate(suffix = stringr::str_replace_all(
            .data$suffix,
            "\\.",
            "\\\\."
        ))
    ## For each row in temp (aka for each ProjectID and
    ## concatenatePoolIDSeqRun) scan the quantification folder
    lups <- purrr::pmap_dfr(temp, function(...) {
        temp_row <- tibble::tibble(...)
        found <- purrr::pmap_dfr(cross, function(...) {
            ## For each quantification scan the folder for the
            ## corresponding suffixes
            cross_row <- tibble::tibble(...)
            matches <- fs::dir_ls(temp_row[[path_col_names$quant]],
                type = "file", fail = FALSE,
                regexp = cross_row$suffix
            )
            if (length(matches) == 0) {
                matches <- NA_character_
            }
            tibble::tibble(
                Quantification_type = cross_row$quant,
                Files_found = matches
            )
        })
        found <- found %>%
            dplyr::group_by(.data$Quantification_type) %>%
            dplyr::distinct() %>%
            dplyr::group_modify(~ {
                if (nrow(.x) > 1) {
                    if (any(is.na(.x$Files_found)) &
                        !all(is.na(.x$Files_found))) {
                        .x %>%
                            dplyr::filter(!is.na(.data$Files_found))
                    } else {
                        .x
                    }
                } else {
                    .x
                }
            })
        tibble::tibble(
            ProjectID = temp_row$ProjectID,
            concatenatePoolIDSeqRun = temp_row$concatenatePoolIDSeqRun,
            found
        )
    })
    lups <- lups %>%
        tidyr::nest(Files = c(.data$Quantification_type, .data$Files_found))
    lups <- .trace_anomalies(lups)
    lups
}


# Simple function to manage user input for duplicate file choices for each
# quantification type in a single project/pool pair (use internally in
# `.manage_anomalies_interactive`).
#
# @param q_types Vector of characters containing the unique quantification
# types that are detected as duplicates
# @param dupl The tibble containing quantification types and path to the files
# found for a single project/pool pair
# @keywords internal
#' @importFrom dplyr filter slice bind_rows
#' @importFrom purrr map reduce
#' @importFrom rlang .data
#
# @return An updated tibble containing files chosen for each type
.choose_duplicates_files_interactive <- function(q_types, dupl) {
    non_dup <- dupl %>% dplyr::filter(!.data$Quantification_type %in% q_types)
    type_choice <- purrr::map(q_types, function(x) {
        cat("---Quantification type: ", x, "---\n\n")
        temp <- dupl %>% dplyr::filter(.data$Quantification_type == x)
        cat("[0] QUIT\n",
            paste0("[", seq_along(temp$Files_found), "] ", temp$Files_found,
                collapse = "\n\n"
            ),
            sep = "\n"
        )
        repeat {
            cat("\n\nType the index number: ")
            connection <- getOption("ISAnalytics.connection")
            if (length(connection) > 1) {
                connection <- connection[7]
            }
            ch <- readLines(con = connection, n = 1)
            ch <- as.numeric(ch)
            if (!is.na(ch) &
                is.numeric(ch) &
                ch >= 0 & ch <= length(temp$Files_found)) {
                if (ch == 0) {
                    stop("Quitting", call. = FALSE)
                } else {
                    break
                }
            }
        }
        temp %>%
            dplyr::slice(ch) %>%
            dplyr::bind_rows(non_dup)
    })
    type_choice <- purrr::reduce(type_choice, dplyr::bind_rows)
}

# Manages anomalies for files found.
#
# The function manages anomalies found for files after scanning appropriate
# folders by:
# * Removing files not found (files for which Files_count$Found == 0 and Path
# is NA) and printing a message to notify user
# * Removing duplicates by asking the user which files to keep for each
# quantification type, pool and project
#
# @param files_found The tibble obtained via calling `.lookup_matrices`
# @keywords internal
#' @importFrom dplyr filter select rename bind_rows arrange
#' @importFrom tidyr unnest
#' @importFrom tibble as_tibble tibble
#' @importFrom purrr pmap flatten reduce
#
# @return A tibble containing for each project, pool and quantification type
# the files chosen (ideally 1 for each quantification type if found, no more
# than 1 per type)
.manage_anomalies_interactive <- function(files_found) {
    # Isolate anomalies in files found
    anomalies <- files_found %>% dplyr::filter(.data$Anomalies == TRUE)
    files_found <- files_found %>% dplyr::filter(.data$Anomalies == FALSE)
    # If there are no anomalies to fix, return chosen files
    if (nrow(anomalies) == 0) {
        files_to_import <- files_found %>%
            dplyr::select(
                .data$ProjectID,
                .data$concatenatePoolIDSeqRun, .data$Files
            ) %>%
            tidyr::unnest(.data$Files) %>%
            dplyr::rename(Files_chosen = "Files_found")
        return(files_to_import)
    }
    repeat {
        # For each ProjectID and PoolID
        files_to_import <-
            purrr::pmap(anomalies, function(...) {
                l <- list(...)
                current <- tibble::as_tibble(l[c(
                    "ProjectID",
                    "concatenatePoolIDSeqRun"
                )])
                current_files <- tibble::as_tibble(purrr::flatten(l["Files"]))
                current_count <- tibble::as_tibble(
                    purrr::flatten(l["Files_count"])
                )
                # Find missing files first
                missing <- current_count %>% dplyr::filter(.data$Found == 0)
                if (nrow(missing) > 0) {
                    message("Some files are missing and will be ignored")
                    current_files <- current_files %>%
                        dplyr::filter(!.data$Quantification_type %in%
                            missing$Quantification_type)
                }
                # Manage duplicates
                duplicate <- current_count %>% dplyr::filter(.data$Found > 1)
                if (nrow(duplicate) > 0) {
                    message("Duplicates found for some files")
                    cat("Plese select one file for each group, type the index ",
                        "as requested\n\n",
                        sep = ""
                    )
                    cat("#### ProjectID: ", current$ProjectID, "####\n\n")
                    cat(
                        "#### Pool: ", current$concatenatePoolIDSeqRun,
                        "####\n\n"
                    )
                    current_files <- .choose_duplicates_files_interactive(
                        duplicate$Quantification_type, current_files
                    )
                }
                to_import <- tibble::tibble(
                    ProjectID = current$ProjectID,
                    concatenatePoolIDSeqRun =
                        current$concatenatePoolIDSeqRun,
                    current_files
                )
                to_import <- to_import %>%
                    dplyr::rename(Files_chosen = "Files_found")
            })
        # Obtain single tibble for all anomalies
        files_to_import <- purrr::reduce(files_to_import, dplyr::bind_rows)
        cat("\nYour choices: ", sep = "\n")
        print(files_to_import)
        cat("\nFiles chosen details:\n")
        cat(paste0(
            seq_along(files_to_import$Files_chosen), "-",
            files_to_import$Files_chosen
        ), sep = "\n\n")
        cat("\nConfirm your choices? [y/n]", sep = "")
        repeat {
            connection <- getOption("ISAnalytics.connection")
            if (length(connection) > 1) {
                connection <- connection[8]
            }
            confirm <- readLines(con = connection, n = 1)
            if (!confirm %in% c("y", "n")) {
                message(paste(
                    "Invalid choice, please select one between",
                    "y or n"
                ))
            } else {
                break
            }
        }
        if (confirm == "y") {
            # Add non-anomalies
            if (nrow(files_found) > 0) {
                files_found <- files_found %>%
                    dplyr::select(
                        .data$ProjectID,
                        .data$concatenatePoolIDSeqRun, .data$Files
                    ) %>%
                    tidyr::unnest(.data$Files) %>%
                    dplyr::rename(Files_chosen = "Files_found")
                files_to_import <- files_to_import %>%
                    dplyr::bind_rows(files_found) %>%
                    dplyr::arrange(.data$ProjectID)
            }
            return(files_to_import)
        } else {
            next
        }
    }
}

# A single threaded and simplified version of import_single_Vispa2Matrix
# to use for parallel import. Preferred instead of calling the exported one
# because:
# - Better control over process forking (don't risk too many threads at once)
# - It was proved that on smaller matrices it is more efficient NOT to split
# and to reshape the entire matrix directly
# - No need for info msgs since there is a final report
#' @importFrom rlang abort
#' @importFrom fs path_ext
#' @importFrom readr read_delim cols
#' @importFrom tidyr separate
#' @importFrom magrittr `%>%`
#' @importFrom dplyr mutate
#' @importFrom stringr str_replace
#' @importFrom data.table melt.data.table
.import_single_matrix <- function(path, to_exclude = NULL, separator = "\t") {
    stopifnot(!missing(path) & is.character(path))
    stopifnot(is.null(to_exclude) || is.character(to_exclude))
    if (!file.exists(path)) {
        rlang::abort(paste("File not found at", path))
    }
    if (!fs::is_file(path)) {
        rlang::abort(paste("Path exists but is not a file"))
    }
    mode <- "fread"
    ## Is the file compressed?
    is_compressed <- fs::path_ext(path) %in% .compressed_formats()
    if (is_compressed) {
        ## The compression type is supported by data.table::fread?
        compression_type <- fs::path_ext(path)
        if (!compression_type %in% .supported_fread_compression_formats()) {
            ### If not, switch to classic for reading
            mode <- "classic"
        }
    }
    ### Peak headers
    peek_headers <- readr::read_delim(path,
        delim = separator, n_max = 1,
        col_types = readr::cols()
    )
    ## - Detect type
    df_type <- .auto_detect_type(peek_headers)
    if (df_type == "MALFORMED") {
        rlang::abort(.malformed_ISmatrix_error(),
            class = "malformed_ism"
        )
    }
    is_annotated <- .is_annotated(peek_headers)
    ## - Start reading
    df <- if (mode == "fread") {
        .read_with_fread(
            path = path, to_drop = to_exclude,
            df_type = df_type, annotated = is_annotated,
            sep = separator
        )
    } else {
        .read_with_readr(
            path = path, to_drop = to_exclude,
            df_type = df_type, annotated = is_annotated,
            sep = separator
        )
    }
    if (df_type == "OLD") {
        df <- df %>%
            tidyr::separate(
                col = .data$IS_genomicID,
                into = mandatory_IS_vars(),
                sep = "_", remove = TRUE,
                convert = TRUE
            ) %>%
            dplyr::mutate(chr = stringr::str_replace(
                .data$chr, "chr", ""
            ))
    }
    ## - Melt
    mt <- function(data, annot) {
        id_vars <- if (annot) {
            c(
                mandatory_IS_vars(),
                annotation_IS_vars()
            )
        } else {
            mandatory_IS_vars()
        }
        data.table::melt.data.table(data,
            id.vars = id_vars,
            variable.name = "CompleteAmplificationID",
            value.name = "Value",
            na.rm = TRUE,
            verbose = FALSE
        )
    }
    tidy <- mt(df, is_annotated)
    tidy <- tidy["Value" > 0]
    return(tidy)
}

# Internal function for parallel import of a single quantification
# type files.
#
# @param q_type The quantification type (single string)
# @param files Files_found table were absolute paths of chosen files
# are stored
# @param workers Number of parallel workers
# @keywords internal
#' @importFrom dplyr filter mutate bind_rows distinct
#' @importFrom BiocParallel SnowParam MulticoreParam bptry bplapply bpstop bpok
#' @importFrom purrr is_empty reduce
#
# @return A single tibble with all data from matrices of same quantification
# type in tidy format
.import_type <- function(q_type, files, workers) {
    files <- files %>% dplyr::filter(.data$Quantification_type == q_type)
    # Register backend according to platform
    if (.Platform$OS.type == "windows") {
        p <- BiocParallel::SnowParam(workers = workers, stop.on.error = FALSE)
    } else {
        p <- BiocParallel::MulticoreParam(
            workers = workers,
            stop.on.error = FALSE
        )
    }
    # Import every file
    FUN <- function(x) {
        matrix <- .import_single_matrix(x)
    }
    suppressMessages(suppressWarnings({
        matrices <- BiocParallel::bptry(
            BiocParallel::bplapply(files$Files_chosen, FUN, BPPARAM = p)
        )
    }))
    BiocParallel::bpstop(p)
    correct <- BiocParallel::bpok(matrices)
    imported_files <- files %>% dplyr::mutate(Imported = correct)
    matrices <- matrices[correct]
    # Bind rows in single tibble for all files
    if (purrr::is_empty(matrices)) {
        return(NULL)
    }
    matrices <- purrr::reduce(matrices, function(x, y) {
        x %>%
            dplyr::bind_rows(y) %>%
            dplyr::distinct()
    })
    list(matrices, imported_files)
}

# Internal function for importing all files for each quantification type.
#
# @param files_to_import The tibble containing the files to import
# @param workers Number of parallel workers
# @keywords internal
#' @importFrom dplyr select distinct bind_rows
#' @importFrom purrr map set_names reduce flatten
#' @importFrom tibble as_tibble
#
# @return A named list of tibbles
.parallel_import_merge <- function(files_to_import, workers) {
    # Find the actual quantification types included
    q_types <- files_to_import %>%
        dplyr::select(.data$Quantification_type) %>%
        dplyr::distinct()
    q_types <- q_types$Quantification_type
    # Import and merge for every quantification type
    imported_matrices <- purrr::map(q_types,
        .f = ~ .import_type(
            .x,
            files_to_import,
            workers
        )
    ) %>%
        purrr::set_names(q_types)
    summary_files <- purrr::map(imported_matrices, function(x) {
        tibble::as_tibble(purrr::flatten(x[2]))
    }) %>% purrr::reduce(dplyr::bind_rows)
    imported_matrices <- purrr::map(imported_matrices, function(x) {
        tibble::as_tibble(purrr::flatten(x[1]))
    })

    list(imported_matrices, summary_files)
}

#---- USED IN : import_parallel_Vispa2Matrices_auto ----

# Internal function to match user defined patterns on a vector
# of file names.
#
# For each pattern specified by the user, the function tries to find a match
# on all the file names and combines the results as a tibble in which column
# names are the patterns and the values are TRUE if the pattern matched on the
# element at that index or FALSE otherwise.
#
# @param filenames A character vector of file names
# @param patterns A character vector of patterns to be matched
# @keywords internal
#' @importFrom tibble as_tibble_col
#' @importFrom stringr str_detect
#' @importFrom purrr map reduce
#' @importFrom dplyr bind_cols
#
# @return A tibble
.pattern_matching <- function(filenames, patterns) {
    p_matches <- purrr::map(patterns, function(x) {
        mtc <- stringr::str_detect(filenames, x)
        tibble::as_tibble_col(mtc, column_name = x)
    })
    p_matches <- purrr::reduce(p_matches, function(x, y) {
        x %>% dplyr::bind_cols(y)
    })
    p_matches
}

# Helper function for checking if any of the elements of the list is true.
#
# @param ... A list of logical values
# @keywords internal
#
# @return TRUE if any of the parameters is true
.any_match <- function(...) {
    l <- unlist(list(...))
    any(l)
}

# Helper function for checking if all of the elements of the list is true.
#
# @param ... A list of logical values
# @keywords internal
#
# @return TRUE if all of the parameters is true
.all_match <- function(...) {
    l <- unlist(list(...))
    all(l)
}

# Updates files_found tibble according to pattern and matching options.
#
# @param files_nested The tibble containing Quantification_type and Files_found
# columns relative to a single project/pool pair
# @param p_matches The tibble representing the pattern matchings resulting from
# `pattern_matching`
# @param matching_opt The matching option
# @keywords internal
#' @importFrom purrr map reduce
#' @importFrom rlang .data
#' @importFrom fs as_fs_path
#
# @return An updated files_found tibble according to the matching option
.update_as_option <- function(files_nested, p_matches, matching_opt) {
    patterns <- colnames(p_matches)
    # Bind columns of Files nested tbl with pattern matches results
    p_matches <- files_nested %>% dplyr::bind_cols(p_matches)
    # Find the different quantification types in the files_found
    types <- p_matches %>%
        dplyr::select(.data$Quantification_type) %>%
        dplyr::distinct()
    # For each quantification type
    to_keep <- purrr::map(types$Quantification_type, function(x) {
        # Obtain a summary of matches (matches ALL patterns or ANY pattern)
        temp <- p_matches %>%
            dplyr::filter(.data$Quantification_type == x) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(
                ALL = .all_match(
                    dplyr::c_across(dplyr::all_of(patterns))
                ),
                ANY = .any_match(dplyr::c_across(dplyr::all_of(patterns)))
            ) %>%
            dplyr::ungroup()
        # If the matching option is "ANY" - preserve files that match any of the
        # given patterns
        if (matching_opt == "ANY") {
            if (!all(is.na(temp$ANY)) & any(temp$ANY)) {
                # If there is at least 1 match in the group, files that match
                # are kept
                files_keep <- temp %>%
                    dplyr::filter(.data$ANY == TRUE) %>%
                    dplyr::select(.data$Quantification_type, .data$Files_found)
            } else {
                # If none of the files match none is preserved, file is
                # converted to NA
                files_keep <- temp %>%
                    dplyr::slice(1) %>%
                    dplyr::mutate(Files_found = fs::as_fs_path(
                        NA_character_
                    )) %>%
                    dplyr::select(.data$Quantification_type, .data$Files_found)
            }
            files_keep
            # If the matching option is "ALL" - preserve only files that match
            # all the given patterns
        } else if (matching_opt == "ALL") {
            if (!all(is.na(temp$ALL)) & any(temp$ALL)) {
                # If there is at least 1 match in the group, files that match
                # are kept
                files_keep <- temp %>%
                    dplyr::filter(.data$ALL == TRUE) %>%
                    dplyr::select(.data$Quantification_type, .data$Files_found)
            } else {
                # If none of the files match none is preserved, file is
                # converted to NA
                files_keep <- temp %>%
                    dplyr::slice(1) %>%
                    dplyr::mutate(Files_found = fs::as_fs_path(
                        NA_character_
                    )) %>%
                    dplyr::select(.data$Quantification_type, .data$Files_found)
            }
            files_keep
            # If the matching option is "OPTIONAL" - preserve preferentially
            # files that match all or any of the given patterns, if none is
            # matching preserves file anyway
        } else {
            # Check first if some files match all the patterns
            if (!all(is.na(temp$ALL)) & any(temp$ALL)) {
                # Preserve only the ones that match
                files_keep <- temp %>%
                    dplyr::filter(.data$ALL == TRUE) %>%
                    dplyr::select(.data$Quantification_type, .data$Files_found)
            } else if (!all(is.na(temp$ANY)) & any(temp$ANY)) {
                # If none match all the patterns check files that match any of
                # the patterns and preserve only those
                files_keep <- temp %>%
                    dplyr::filter(.data$ANY == TRUE) %>%
                    dplyr::select(.data$Quantification_type, .data$Files_found)
            } else {
                # If there are no files that match any of the patterns simply
                # preserve the files found
                files_keep <- temp %>%
                    dplyr::select(.data$Quantification_type, .data$Files_found)
            }
            files_keep
        }
    })
    to_keep <- purrr::reduce(to_keep, dplyr::bind_rows)
    to_keep
}

# Looks up matrices to import given the association file and the
# root of the file system.
#
# @inheritParams .lookup_matrices
# @param patterns A character vector of patterns to be matched
# @param matching_opt A single character representing the matching option (one
# of "ANY", "ALL" or "OPTIONAL")
# @keywords internal
#' @importFrom purrr pmap flatten map
#' @importFrom tibble as_tibble
#' @importFrom stringr str_split
#' @importFrom dplyr mutate select
#' @importFrom utils tail
#
# @return A tibble containing all found files, including duplicates and missing
.lookup_matrices_auto <- function(association_file,
    quantification_type,
    matrix_type,
    patterns,
    matching_opt) {
    files_found <- .lookup_matrices(
        association_file,
        quantification_type,
        matrix_type
    )
    if (nrow(files_found) > 0 & !is.null(patterns)) {
        matching_patterns <- purrr::pmap(files_found, function(...) {
            l <- list(...)
            files_nested <- tibble::as_tibble(purrr::flatten(l["Files"]))
            splitted <- stringr::str_split(files_nested$Files_found, "\\/")
            filenames <- unlist(purrr::map(splitted, function(x) {
                utils::tail(x, n = 1)
            }))
            p_matches <- .pattern_matching(filenames, patterns)
            to_keep <- .update_as_option(files_nested, p_matches, matching_opt)
        })
        files_found <- files_found %>%
            dplyr::mutate(Files = matching_patterns) %>%
            dplyr::select(
                .data$ProjectID, .data$concatenatePoolIDSeqRun,
                .data$Files
            )
        files_found <- .trace_anomalies(files_found)
    }
    files_found
}

# Manages anomalies for files found.
#
# The function manages anomalies found for files after scanning appropriate
# folders by:
# * Removing files not found (files for which Files_count$Found == 0 and Path
# is NA) and printing a message to notify user
# * Removing duplicates
#
# @param files_found The tibble obtained via calling `.lookup_matrices_auto`
# @keywords internal
#' @importFrom dplyr filter select rename bind_rows arrange
#' @importFrom rlang .data
#' @importFrom tidyr unnest
#' @importFrom tibble as_tibble tibble
#' @importFrom purrr pmap flatten
#
# @return A tibble containing for each project, pool and quantification type
# the files chosen
.manage_anomalies_auto <- function(files_found) {
    # Isolate anomalies in files found
    anomalies <- files_found %>% dplyr::filter(.data$Anomalies == TRUE)
    files_found <- files_found %>% dplyr::filter(.data$Anomalies == FALSE)
    # If there are no anomalies to fix, return chosen files
    if (nrow(anomalies) == 0) {
        files_to_import <- files_found %>%
            dplyr::select(
                .data$ProjectID,
                .data$concatenatePoolIDSeqRun, .data$Files
            ) %>%
            tidyr::unnest(.data$Files) %>%
            dplyr::rename(Files_chosen = "Files_found")
        return(files_to_import)
    }
    # For each ProjectID and PoolID
    files_to_import <- purrr::pmap(anomalies, function(...) {
        l <- list(...)
        current <- tibble::as_tibble(l[c(
            "ProjectID",
            "concatenatePoolIDSeqRun"
        )])
        current_files <- tibble::as_tibble(purrr::flatten(l["Files"]))
        current_count <- tibble::as_tibble(purrr::flatten(l["Files_count"]))
        # Find missing files first
        missing <- current_count %>% dplyr::filter(.data$Found == 0)
        if (nrow(missing) > 0) {
            message("Some files are missing and will be ignored")
            current_files <- current_files %>%
                dplyr::filter(!.data$Quantification_type %in%
                    missing$Quantification_type)
        }
        # Manage duplicates
        duplicate <- current_count %>% dplyr::filter(.data$Found > 1)
        if (nrow(duplicate) > 0) {
            message(paste(
                "Duplicates found for some files: in automatic mode",
                "duplicates are not preserved - use interactive mode for more",
                "accurate file selection"
            ))
            current_files <- current_files %>%
                dplyr::filter(!.data$Quantification_type %in%
                    duplicate$Quantification_type)
        }
        to_import <- tibble::tibble(
            ProjectID = current$ProjectID,
            concatenatePoolIDSeqRun =
                current$concatenatePoolIDSeqRun,
            current_files
        )
        to_import <- to_import %>% dplyr::rename(Files_chosen = "Files_found")
    })
    # Obtain single tibble for all anomalies
    files_to_import <- purrr::reduce(files_to_import, dplyr::bind_rows)
    # Add non-anomalies
    if (nrow(files_found) > 0) {
        files_found <- files_found %>%
            dplyr::select(
                .data$ProjectID,
                .data$concatenatePoolIDSeqRun, .data$Files
            ) %>%
            tidyr::unnest(.data$Files) %>%
            dplyr::rename(Files_chosen = "Files_found")
        files_to_import <- files_to_import %>%
            dplyr::bind_rows(files_found) %>%
            dplyr::arrange(.data$ProjectID)
    }
    files_to_import
}

#### ---- Internals for collision removal ----####

#---- USED IN : remove_collisions ----

# Checks if association file contains more information than the matrix.
#
# Used to notify the user that wants to know if for the projects contained in
# the examined matrix there are additional CompleteAmplificationIDs contained
# in the association file that weren't included in the integration matrix (for
# example failed runs).
#' @importFrom rlang .data
.check_same_info <- function(association_file, df) {
    joined <- association_file %>%
        dplyr::left_join(df, by = "CompleteAmplificationID")
    projects <- joined %>%
        dplyr::filter(!is.na(.data$chr)) %>%
        dplyr::distinct(.data$ProjectID) %>%
        dplyr::pull(.data$ProjectID)
    af_cols <- colnames(association_file)
    wanted <- c(
        "ProjectID", "PoolID", "CompleteAmplificationID",
        "SubjectID", "Tissue", "CellMarker"
    )
    cols <- wanted[wanted %in% af_cols]
    missing_in_df <- joined %>%
        dplyr::filter(is.na(.data$chr), .data$ProjectID %in% projects) %>%
        dplyr::select(dplyr::all_of(cols)) %>%
        dplyr::distinct()
    reduced_af <- association_file %>%
        dplyr::filter(.data$ProjectID %in% projects)
    list(miss = missing_in_df, reduced_af = reduced_af)
}

# Produces a joined tibble between the sequence count matrix and the
# association file
#
# @param seq_count_df The sequence count tibble
# @param association_file The association file tibble
# @param date_col The date column chosen
# @keywords internal
#' @importFrom rlang .data
#
# @return A tibble
.join_matrix_af <- function(df, association_file, date_col) {
    joined <- df %>%
        dplyr::left_join(association_file, by = "CompleteAmplificationID") %>%
        dplyr::select(
            dplyr::all_of(colnames(df)), all_of(date_col),
            .data$ProjectID, .data$PoolID, .data$SubjectID,
            .data$ReplicateNumber
        )
    joined
}

# Identifies independent samples and separates the joined_df in
# collisions and non-collisions
#
# @param joined_df The joined tibble obtained via `.join_matrix_af`
#' @importFrom rlang .data
# @keywords internal
#
# @return A named list containing the splitted joined_df for collisions and
# non-collisions
.identify_independent_samples <- function(joined_df) {
    temp <- joined_df %>%
        dplyr::select(
            dplyr::all_of(mandatory_IS_vars()),
            .data$ProjectID, .data$SubjectID
        ) %>%
        dplyr::group_by(dplyr::across(mandatory_IS_vars())) %>%
        dplyr::distinct(.data$ProjectID, .data$SubjectID, .keep_all = TRUE) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
        dplyr::filter(.data$n > 1) %>%
        dplyr::select(-c(.data$n))

    non_collisions <- joined_df %>%
        dplyr::anti_join(temp, by = mandatory_IS_vars())
    collisions <- joined_df %>%
        dplyr::right_join(temp, by = mandatory_IS_vars())
    list(collisions = collisions, non_collisions = non_collisions)
}

# Returns a polished version of collisions table for analysis.
#
# Polishing consists in nesting all information relative to a single
# integration, in other words, obtaining a linked data frame for each
# integration.
#
# @param collisions The collisions table obtained via
# `.identify_independent_samples`
# @keywords internal
#' @importFrom dplyr group_by across
#' @importFrom rlang .data
#' @importFrom tidyr nest
#
# @return A nested tibble
.obtain_nested <- function(collisions) {
    collisions %>%
        dplyr::group_by(dplyr::across(mandatory_IS_vars())) %>%
        tidyr::nest()
}


# Internal for date discrimination in collision removal.
#
# It's the first of 4 steps in the algorithm for collision removal: it tries to
# find a single sample who has an associated date which is earlier than any
# other. If comparison is not possible the analysis fails and returns the
# orginal input data.
#
# @details NOTE: this function is meant to be used inside a mapping function
# such as `purrr::pmap`. The function only works on data regarding a SINGLE
# integration (triplet chr, integration_locus, strand).
#
# @param nest The nested table associated with a single integration
# @param date_col The name of the date column chosen for collision removal,
# one in \code{date_columns_coll()}
#' @importFrom rlang expr eval_tidy .data
#' @importFrom dplyr filter arrange
# @keywords internal
#
# @return A named list with:
# * $data: a tibble, containing the data (unmodified or modified)
# * $check: a logical value indicating whether the analysis was successful or
# not (and therefore there is the need to perform the next step)
.discriminate_by_date <- function(nest, date_col) {
    dates <- rlang::expr(`$`(nest, !!date_col))
    dates <- rlang::eval_tidy(dates)
    # Test for all dates equality
    all_equal_dates <- length(unique(dates)) == 1
    if (all_equal_dates == TRUE) {
        return(list(data = nest, check = FALSE))
    }
    # If not all are equal sort them asc
    nest <- nest %>% dplyr::arrange(`$`(.data, !!date_col))
    # Test if first date is unique
    dates <- rlang::expr(`$`(nest, !!date_col))
    dates <- rlang::eval_tidy(dates)
    if (length(dates[dates %in% dates[1]]) > 1) {
        return(list(data = nest, check = FALSE))
    }
    winning_subj <- nest$SubjectID[1]
    # Filter the winning rows
    nest <- nest %>% dplyr::filter(.data$SubjectID == winning_subj)
    return(list(data = nest, check = TRUE))
}

# Internal for replicate discrimination in collision removal.
#
# It's the second of 4 steps in the algorithm for collision removal:
# grouping by independent sample (ProjectID, SubjectID) it counts the number of
# rows (replicates) found for each group, orders them from biggest to smallest
# and, if a single group has more rows than any other group the integration is
# assigned to that sample, otherwise the analysis fails and returns the
# original input to be submitted to the next step.
#
# @details NOTE: this function is meant to be used inside a mapping function
# such as `purrr::pmap`. The function only works on data regarding a SINGLE
# integration (triplet chr, integration_locus, strand).
#
# @param nest The nested table associated with a single integration
#' @importFrom rlang .data
# @keywords internal
#
# @return A named list with:
# * $data: a tibble, containing the data (unmodified or modified)
# * $check: a logical value indicating whether the analysis was successful or
# not (and therefore there is the need to perform the next step)
.discriminate_by_replicate <- function(nest) {
    temp <- nest %>%
        dplyr::group_by(.data$ProjectID, .data$SubjectID) %>%
        dplyr::summarise(n_rows = dplyr::n(), .groups = "keep") %>%
        dplyr::arrange(dplyr::desc(.data$n_rows))
    if (length(temp$n_rows) != 1 & !temp$n_rows[1] > temp$n_rows[2]) {
        return(list(data = nest, check = FALSE))
    }
    nest <- nest %>%
        dplyr::filter(
            .data$ProjectID == temp$ProjectID[1],
            .data$SubjectID == temp$SubjectID[1]
        )
    return(list(data = nest, check = TRUE))
}

# Internal for sequence count discrimination in collision removal.
#
# It's the third of 4 steps in the algorithm for collision removal:
# grouping by independent sample (ProjectID, SubjectID), it sums the value of
# the sequence count for each group, sorts the groups from highest to lowest
# cumulative value and then checks the ratio between the first element and the
# second: if the ratio is > `reads_ratio` the integration is assigned to
# the first group, otherwise the analysis fails and returns the original input.
#
# @details NOTE: this function is meant to be used inside a mapping function
# such as `purrr::pmap`. The function only works on data regarding a SINGLE
# integration (triplet chr, integration_locus, strand).
#
# @param nest The nested table associated with a single integration
# @param reads_ratio The value of the ratio between sequence count values to
# check
# @param seqCount_col The name of the sequence count column (support
# for multi quantification matrix)
#' @importFrom rlang .data eval_tidy expr
# @keywords internal
#
# @return A named list with:
# * $data: a tibble, containing the data (unmodified or modified)
# * $check: a logical value indicating whether the analysis was successful or
# not (and therefore there is the need to perform the next step)
.discriminate_by_seqCount <- function(nest, reads_ratio, seqCount_col) {
    temp <- nest %>%
        dplyr::group_by(.data$ProjectID, .data$SubjectID) %>%
        dplyr::summarise(sum_seqCount = sum(
            rlang::eval_tidy(rlang::expr(`$`(.data, !!seqCount_col)))
        ), .groups = "keep") %>%
        dplyr::arrange(dplyr::desc(.data$sum_seqCount))

    ratio <- temp$sum_seqCount[1] / temp$sum_seqCount[2]
    if (!ratio > reads_ratio) {
        return(list(data = nest, check = FALSE))
    }
    nest <- nest %>%
        dplyr::filter(
            .data$ProjectID == temp$ProjectID[1],
            .data$SubjectID == temp$SubjectID[1]
        )
    return(list(data = nest, check = TRUE))
}

# Internal function that performs four-step-check of collisions for
# a single integration.
#
# @details NOTE: this function is meant to be used inside a mapping function
# such as `purrr::pmap`. The function only works on data regarding a SINGLE
# integration (triplet chr, integration_locus, strand).
#
# @param ... Represents a single row of a tibble obtained via
# `obtain_nested`. One row contains 4 variables: chr, integration_locus,
# strand and data, where data is a nested table that contains all rows
# that share that integration (collisions).
# @param date_col The date column to consider
# @param reads_ratio The value of the ratio between sequence count values to
# check
# @param seqCount_col The name of the sequence count column (support
# for multi quantification matrix)
# @keywords internal
#
#' @importFrom tibble as_tibble tibble
#' @importFrom purrr flatten
#' @importFrom tidyr unnest
#' @importFrom rlang .data env_bind
# @return A list with:
# * $data: an updated tibble with processed collisions or NULL if no
# criteria was sufficient
# * $reassigned: 1 if the integration was successfully reassigned, 0 otherwise
# * $removed: 1 if the integration is removed entirely because no criteria was
# met, 0 otherwise
.four_step_check <- function(..., date_col, reads_ratio, seqCount_col) {
    l <- list(...)
    current <- tibble::as_tibble(l[mandatory_IS_vars()])
    current_data <- tibble::as_tibble(purrr::flatten(l["data"]))
    # Try to discriminate by date
    result <- .discriminate_by_date(current_data, date_col)
    if (result$check == TRUE) {
        current_data <- result$data
        res <- tibble::tibble(
            chr = current$chr,
            integration_locus = current$integration_locus,
            strand = current$strand,
            data = list(current_data)
        )
        res <- res %>% tidyr::unnest(.data$data)
        return(list(data = res, reassigned = 1, removed = 0))
    }
    current_data <- result$data
    # If first check fails try to discriminate by replicate
    result <- .discriminate_by_replicate(current_data)
    if (result$check == TRUE) {
        current_data <- result$data
        res <- tibble::tibble(
            chr = current$chr,
            integration_locus = current$integration_locus,
            strand = current$strand,
            data = list(current_data)
        )
        res <- res %>% tidyr::unnest(.data$data)
        return(list(data = res, reassigned = 1, removed = 0))
    }
    current_data <- result$data
    # If second check fails try to discriminate by seqCount
    result <- .discriminate_by_seqCount(current_data, reads_ratio, seqCount_col)
    if (result$check == TRUE) {
        current_data <- result$data
        res <- tibble::tibble(
            chr = current$chr,
            integration_locus = current$integration_locus,
            strand = current$strand,
            data = list(current_data)
        )
        res <- res %>% tidyr::unnest(.data$data)
        return(list(data = res, reassigned = 1, removed = 0))
    }
    # If all check fails remove the integration from all subjects
    return(list(data = NULL, reassigned = 0, removed = 1))
}

# Internal function to process collisions on multiple integrations.
#
# @param x The nested collision data frame (or chunks for parallel execution)
# @param date_col The date column to consider
# @param reads_ratio The value of the ratio between sequence count values to
# check
# @param seqCount_col The name of the sequence count column (support
# for multi quantification matrix)
# @keywords internal
#' @importFrom purrr pmap reduce
#' @importFrom dplyr bind_rows
# @return A list containing the updated collisions, a numeric value
# representing the number of integrations removed and a numeric value
# representing the number of integrations reassigned
.coll_mapping <- function(x, date_col, reads_ratio, seqCount_col) {
    result <- purrr::pmap(x,
        .f = .four_step_check,
        date_col = date_col,
        reads_ratio = reads_ratio,
        seqCount_col = seqCount_col
    )
    proc_collisions <- purrr::map(result, function(x) {
        x$data
    })
    proc_collisions <- purrr::reduce(proc_collisions, dplyr::bind_rows)
    removed_tot <- purrr::map(result, function(x) {
        x$removed
    })
    removed_tot <- purrr::reduce(removed_tot, sum)
    reassigned_tot <- purrr::map(result, function(x) {
        x$reassigned
    })
    reassigned_tot <- purrr::reduce(reassigned_tot, sum)
    list(
        coll = proc_collisions,
        removed = removed_tot,
        reassigned = reassigned_tot
    )
}

# Internal function to process collisions on multiple integrations,
# parallelized.
#
# @param collisions The collision data frame
# @param date_col The date column to consider
# @param reads_ratio The value of the ratio between sequence count values to
# check
# @param seqCount_col The name of the sequence count column (support
# for multi quantification matrix)
# @keywords internal
#' @importFrom purrr map reduce
#' @importFrom tibble as_tibble_col add_column
#' @importFrom dplyr group_by group_split
# @return A list containing the updated collisions, a numeric value
# representing the number of integrations removed and a numeric value
# representing the number of integrations reassigned
.process_collisions <- function(collisions, date_col, reads_ratio,
    seqCount_col, max_workers) {
    # Obtain nested version of collisions
    nested <- .obtain_nested(collisions)
    # Manage workers
    if (.Platform$OS.type == "windows" & is.null(max_workers)) {
        max_workers <- BiocParallel::snowWorkers()
    } else if (.Platform$OS.type != "windows" & is.null(max_workers)) {
        max_workers <- BiocParallel::multicoreWorkers()
    }
    # Split the data frame in chunks according to number of workers
    exceeding <- nrow(nested) %% max_workers
    ampl <- trunc(nrow(nested) / max_workers)
    chunks <- unlist(lapply(seq_len(max_workers),
        FUN = rep_len,
        length.out = ampl
    ))
    if (exceeding > 0) {
        chunks <- c(chunks, rep_len(
            x = tail(max_workers, n = 1),
            length.out = exceeding
        ))
    }
    chunks <- tibble::as_tibble_col(chunks, column_name = "chunk")
    nested <- nested %>% tibble::add_column(chunks)
    split_data <- nested %>%
        dplyr::group_by(.data$chunk) %>%
        dplyr::group_split(.keep = FALSE)

    # Register backend according to platform
    if (.Platform$OS.type == "windows") {
        p <- BiocParallel::SnowParam(
            stop.on.error = FALSE,
            progressbar = getOption("ISAnalytics.verbose", default = TRUE),
            tasks = nrow(nested),
            workers = max_workers
        )
    } else {
        p <- BiocParallel::MulticoreParam(
            stop.on.error = FALSE,
            progressbar = getOption("ISAnalytics.verbose", default = TRUE),
            tasks = nrow(nested),
            workers = max_workers
        )
    }

    # For each chunk process collisions in parallel
    suppressMessages(suppressWarnings({
        ### result is a list with n elements, each element is a list of 3,
        ### containing processed collisions, number of removed collisions and
        ### number of reassigned collisions
        result <- BiocParallel::bptry(
            BiocParallel::bplapply(split_data,
                FUN = .coll_mapping,
                date_col = date_col,
                reads_ratio = reads_ratio,
                seqCount_col = seqCount_col,
                BPPARAM = p
            )
        )
    }))
    BiocParallel::bpstop(p)
    # For each element of result extract and bind the 3 components
    processed_collisions_df <- purrr::map(result, function(x) {
        x$coll
    })
    processed_collisions_df <- purrr::reduce(
        processed_collisions_df,
        dplyr::bind_rows
    )
    removed_total <- purrr::map(result, function(x) {
        x$removed
    })
    removed_total <- purrr::reduce(removed_total, sum)
    reassigned_total <- purrr::map(result, function(x) {
        x$reassigned
    })
    reassigned_total <- purrr::reduce(reassigned_total, sum)
    # Return the summary
    list(
        coll = processed_collisions_df, removed = removed_total,
        reassigned = reassigned_total
    )
}

# Internal function to obtain a summary table after collision processing.
#
# @param before The joined table before collision removing
# (obtained via `.join_matrix_af`)
# @param after The final matrix obtained after collision processing
# @param association_file The association file
#' @importFrom rlang .data
# @keywords internal
#
# @return A tibble with a summary containing for each SubjectID the number of
# integrations found before and after, the sum of the value of the sequence
# count for each subject before and after and the corresponding deltas.
.summary_table <- function(before, after, seqCount_col) {
    after_matr <- after %>%
        dplyr::select(dplyr::all_of(colnames(after)), .data$SubjectID)

    n_int_after <- after_matr %>%
        dplyr::group_by(.data$SubjectID) %>%
        dplyr::tally(name = "count_integrations_after")
    sumReads_after <- after_matr %>%
        dplyr::group_by(.data$SubjectID) %>%
        dplyr::summarise(
            sumSeqReads_after = sum(.data[[seqCount_col]]),
            .groups = "drop_last"
        )
    temp_aft <- n_int_after %>% dplyr::left_join(sumReads_after,
        by = c("SubjectID")
    )

    before_matr <- before %>%
        dplyr::select(dplyr::all_of(colnames(after)), .data$SubjectID)
    n_int_before <- before_matr %>%
        dplyr::group_by(.data$SubjectID) %>%
        dplyr::tally(name = "count_integrations_before")
    sumReads_before <- before_matr %>%
        dplyr::group_by(.data$SubjectID) %>%
        dplyr::summarise(
            sumSeqReads_before = sum(.data[[seqCount_col]]),
            .groups = "drop_last"
        )
    temp_bef <- n_int_before %>% dplyr::left_join(sumReads_before,
        by = c("SubjectID")
    )
    summary <- temp_bef %>%
        dplyr::left_join(temp_aft, by = c("SubjectID")) %>%
        dplyr::mutate(
            delta_integrations =
                .data$count_integrations_before -
                    .data$count_integrations_after,
            delta_seqReads = .data$sumSeqReads_before - .data$sumSeqReads_after
        )
    summary
}

# Internal for obtaining summary info about the input sequence count matrix
# as is (no pre-processing).
# Info contain
# - Number of distinct integration sites
# - Total for each quantification present (if multi-quant) or seqCount
.summary_input <- function(input_df, quant_cols) {
    names(quant_cols) <- NULL
    # Distinct integration sites
    n_IS <- input_df %>%
        dplyr::distinct(dplyr::across(mandatory_IS_vars())) %>%
        nrow()
    quant_totals <- input_df %>%
        dplyr::select(dplyr::all_of(quant_cols)) %>%
        dplyr::summarise(dplyr::across(
            .cols = dplyr::all_of(quant_cols),
            .fns = list(
                sum = ~ sum(.x, na.rm = TRUE)
            ),
            .names = "{.col}"
        ),
        .groups = "drop"
        ) %>%
        as.list()
    list(total_iss = n_IS, quant_totals = quant_totals)
}

.per_pool_stats <- function(joined, quant_cols) {
    ## Joined is the matrix already joined with metadata
    df_by_pool <- joined %>%
        dplyr::group_by(.data$PoolID) %>%
        dplyr::summarise(dplyr::across(
            .cols = dplyr::all_of(quant_cols),
            .fns = list(
                sum = ~ sum(.x, na.rm = TRUE),
                count = length,
                describe = ~ psych::describe(.x)
            )
        ))
    df_cols <- purrr::map_lgl(df_by_pool, ~ is.data.frame(.x))
    if (any(df_cols == TRUE)) {
        df_cols <- df_cols[df_cols]
        df_cols <- names(df_cols)
        dfss <- purrr::map(df_cols, function(col) {
            exp <- rlang::expr(`$`(df_by_pool, !!col))
            df <- rlang::eval_tidy(exp)
            df <- df %>% dplyr::rename_with(.fn = ~ paste0(col, "_", .x))
            df
        }) %>% purrr::set_names(df_cols)
        for (dfc in df_cols) {
            df_by_pool <- df_by_pool %>%
                dplyr::select(-dplyr::all_of(dfc)) %>%
                dplyr::bind_cols(dfss[[dfc]])
        }
    }
    df_by_pool
}

#### ---- Internals for aggregate functions ----####

#---- USED IN : aggregate_metadata ----

# Minimal stats column set.
#
# Contains the name of the columns that are a minimum requirement for
# aggregation.
# @keywords internal
#
# @return A character vector
.stats_columns_min <- function() {
    c(
        "POOL", "TAG", "BARCODE_MUX"
    )
}

# Checks if the stats file contains the minimal set of columns.
#
# @param x The stats df
# @keywords internal
#
# @return TRUE or FALSE
.check_stats <- function(x) {
    if (all(.stats_columns_min() %in% colnames(x))) {
        TRUE
    } else {
        FALSE
    }
}

# Aggregates the association file based on the function table.
#' @importFrom rlang .data is_function is_formula
#' @importFrom tibble tibble
#' @importFrom purrr pmap_df
# @keywords internal
#
# @return A tibble
.aggregate_meta <- function(association_file, grouping_keys, function_tbl) {
    ## Discard columns that are not present in the af
    function_tbl <- function_tbl %>%
        dplyr::filter(.data$Column %in% colnames(association_file))
    ## If no columns are left return
    if (nrow(function_tbl) == 0) {
        return(NULL)
    }
    function_tbl <- function_tbl %>%
        tidyr::nest(cols = .data$Column)
    apply_function <- function(Function, Args, Output_colname, cols) {
        Function <- list(Function)
        Args <- list(Args)
        row <- tibble::tibble(Function, Args,
            Output_colname,
            cols = list(cols)
        )
        if (rlang::is_formula(row$Function[[1]])) {
            res <- association_file %>%
                dplyr::group_by(dplyr::across(dplyr::all_of(grouping_keys))) %>%
                dplyr::summarise(
                    dplyr::across(
                        dplyr::all_of(row$cols[[1]]$Column),
                        .fns = row$Function[[1]],
                        .names = row$Output_colname
                    ),
                    .groups = "drop"
                )
            return(res)
        }
        if (rlang::is_function(row$Function[[1]])) {
            args <- if (is.na(row$Args[[1]]) | is.null(is.na(row$Args[[1]]))) {
                list()
            } else {
                row$Args[[1]]
            }
            res <- association_file %>%
                dplyr::group_by(dplyr::across(dplyr::all_of(grouping_keys))) %>%
                dplyr::summarise(
                    dplyr::across(
                        .cols = dplyr::all_of(row$cols[[1]]$Column),
                        .fns = row$Function[[1]],
                        .names = row$Output_colname,
                        !!!args
                    ),
                    .groups = "drop"
                )
            return(res)
        }
        return(NULL)
    }
    agg <- purrr::pmap(function_tbl, apply_function)
    agg <- purrr::reduce(agg, ~ dplyr::full_join(.x, .y, by = grouping_keys))
    return(agg)
}

#---- USED IN : aggregate_values_by_key ----

# Internal function that performs aggregation on values with a lambda
# operation.
#
# @param x The list of matrices to aggregate. If a single matrix has to be
# supplied it must be enclosed in a list. For example `x = list(matrix)`.
# @param af The association file
# @param key A string or character vector to use as key
# @param lambda The aggregating operation to apply to values. Must take as
# input a numeric/integer vector and return a single value
# @param group The additional variables to add to grouping
# @param args Additional arguments passed on to lambda (named list)
# @param namespace The namespace from which the lambda is exported
# @param envir The environment in which symbols must be evaluated
#' @importFrom rlang .data
# @keywords internal
#
# @return A list of tibbles with aggregated values
.aggregate_lambda <- function(x, af, key, value_cols, lambda, group,
    join_af_by) {
    # Vectorize
    aggregated_matrices <- purrr::map(x, function(y) {
        cols <- c(colnames(y), key)
        agg <- y %>%
            dplyr::left_join(af, by = dplyr::all_of(join_af_by)) %>%
            dplyr::select(dplyr::all_of(cols)) %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(c(group, key)))) %>%
            dplyr::summarise(dplyr::across(
                .cols = dplyr::all_of(value_cols),
                .fns = lambda
            ), .groups = "drop")
        agg
    })
    aggregated_matrices
}


#### ---- Internals for re-calibration functions ----####

#---- USED IN : compute_near_integrations ----

# Finds a unique maximum value in a vector of numbers.
#
# Unlike the standard `max` function, this function returns a single maximum
# value if and only if a maximum exists, more in details that means:
# * It ignores NA values by default
# * If the identified maximum value is not unique in the vector, it returns
# an empty vector instead
#
# @param x A numeric or integer vector
#' @importFrom purrr is_empty
#
# @return A single numeric value or an empty numeric vector
# @keywords internal
.find_unique_max <- function(x) {
    if (any(is.na(x))) {
        x <- x[!is.na(x)]
    }
    uniques <- unique(x)
    if (purrr::is_empty(uniques)) {
        return(uniques)
    }
    if (length(uniques) == 1) {
        if (length(x[x == uniques]) > 1) {
            return(numeric(0))
        } else {
            return(uniques)
        }
    }
    max_val <- max(uniques)
    if (length(x[x == max_val]) > 1) {
        return(numeric(0))
    }
    max_val
}

# Internal function that implements the sliding window algorithm.
#
# **NOTE: this function is meant to be called on a SINGLE GROUP,
# meaning a subset of an integration matrix in which all rows
# share the same chr (and optionally same strand).**
#
# @param x An integration matrix subset (see description)
# @param threshold The numeric value representing an absolute number
# of bases for which two integrations are considered distinct
# @param keep_criteria The string with selecting criteria
# @param annotated Is `x` annotated? Logical value
# @param num_cols A character vector with the names of the numeric columns
# @param max_val_col The column to consider if criteria is max_value
# @param produce_map Produce recalibration map?
#' @importFrom dplyr arrange
#' @importFrom purrr is_empty map_dfr
#' @importFrom rlang expr eval_tidy .data
#' @importFrom data.table setDT setnames `%chin%`
# @keywords internal
#
# @return A named list with recalibrated matrix and recalibration map.
.sliding_window <- function(x,
    threshold,
    keep_criteria,
    annotated,
    num_cols,
    max_val_col,
    produce_map) {
    ## Order by integration_locus
    x <- x %>% dplyr::arrange(.data$integration_locus)
    x <- data.table::setDT(x)
    map_recalibr <- if (produce_map == TRUE) {
        mand_vars <- mandatory_IS_vars()
        rec_map <- unique(x[, ..mand_vars])
        data.table::setnames(
            x = rec_map,
            old = mandatory_IS_vars(),
            new = paste0(mandatory_IS_vars(), "_before")
        )
        rec_map[, c("chr_after", "integration_locus_after", "strand_after") :=
            list(NA_character_, NA_real_, NA_character_)]
    } else {
        NULL
    }
    index <- 1
    repeat {
        # If index is the last row in the data frame return
        if (index == nrow(x)) {
            if (produce_map == TRUE) {
                map_recalibr[
                    chr_before == x[index, ]$chr &
                        integration_locus_before ==
                            x[index, ]$integration_locus &
                        strand_before == x[index, ]$strand,
                    c(
                        "chr_after",
                        "integration_locus_after",
                        "strand_after"
                    ) :=
                        list(
                            chr_before,
                            integration_locus_before,
                            strand_before
                        )
                ]
            }
            return(list(recalibrated_matrix = x, map = map_recalibr))
        }
        ## Compute interval for row
        interval <- x[index, ]$integration_locus + 1 + threshold
        ## Look ahead for every integration that falls in the interval
        near <- numeric()
        k <- index
        repeat {
            if (k == nrow(x)) {
                break
            }
            k <- k + 1
            if (x[k, ]$integration_locus < interval) {
                # Saves the indexes of the rows that are in the interval
                near <- append(near, k)
            } else {
                break
            }
        }
        window <- c(index, near)
        if (!purrr::is_empty(near)) {
            ## Change loci according to criteria
            ######## CRITERIA PROCESSING
            row_to_keep <- index
            if (keep_criteria == "max_value") {
                expr <- rlang::expr(`$`(x[window, ], !!max_val_col))
                to_check <- rlang::eval_tidy(expr)
                max <- .find_unique_max(to_check)
                if (!purrr::is_empty(max)) {
                    row_to_keep <- window[which(to_check == max)]
                    near <- window[!window == row_to_keep]
                }
            }
            # Fill map if needed
            if (produce_map == TRUE) {
                map_recalibr[
                    chr_before %chin% x[window, ]$chr &
                        integration_locus_before %in%
                            x[window, ]$integration_locus &
                        strand_before %chin% x[window, ]$strand,
                    c(
                        "chr_after",
                        "integration_locus_after",
                        "strand_after"
                    ) :=
                        list(
                            x[row_to_keep, ]$chr,
                            x[row_to_keep, ]$integration_locus,
                            x[row_to_keep, ]$strand
                        )
                ]
            }
            # Change loci and strand of near integrations
            if (annotated) {
                x[near, c(
                    "integration_locus",
                    "strand",
                    annotation_IS_vars()
                ) := list(
                    x[row_to_keep, ]$integration_locus,
                    x[row_to_keep, ]$strand,
                    x[row_to_keep, ]$GeneName,
                    x[row_to_keep, ]$GeneStrand
                )]
            } else {
                x[near, c("integration_locus", "strand") := list(
                    x[row_to_keep, ]$integration_locus,
                    x[row_to_keep, ]$strand
                )]
            }
            ## Aggregate same IDs
            starting_rows <- nrow(x)
            repeat {
                t <- x$CompleteAmplificationID[window]
                d <- unique(t[duplicated(t)])
                if (purrr::is_empty(d)) {
                    break
                }
                dupl_indexes <- which(t == d[1])
                values_sum <- colSums(x[window[dupl_indexes], ..num_cols],
                    na.rm = TRUE
                )
                x[window[dupl_indexes[1]], c(num_cols) := as.list(values_sum)]
                x <- x[-window[dupl_indexes[-1]], ]
                to_drop <- seq(
                    from = length(window),
                    length.out = length(dupl_indexes[-1]),
                    by = -1
                )
                window <- window[-to_drop]
            }
            ## Increment index
            if (nrow(x) == starting_rows) {
                index <- k
            } else {
                index <- tail(window, n = 1) + 1
            }
        } else {
            if (produce_map == TRUE) {
                map_recalibr[
                    chr_before == x[index, ]$chr &
                        integration_locus_before ==
                            x[index, ]$integration_locus &
                        strand_before == x[index, ]$strand,
                    c(
                        "chr_after",
                        "integration_locus_after",
                        "strand_after"
                    ) :=
                        list(
                            chr_before,
                            integration_locus_before,
                            strand_before
                        )
                ]
            }
            index <- index + 1
        }
    }
    list(recalibrated_matrix = x, map = map_recalibr)
}


# Generates a file name for recalibration maps.
#
# Unique file names include current date and timestamp.
# @keywords internal
#' @importFrom stringr str_replace_all
#
# @return A string
.generate_rec_map_filename <- function() {
    def <- "recalibration_map.tsv.gz"
    date <- lubridate::today()
    return(paste0(date, "_", def))
}

# Tries to write the recalibration map to a tsv file.
#
# @param map The recalibration map
# @param file_path The file path as a string
#' @importFrom fs dir_create path_wd path
#' @importFrom readr write_tsv
# @keywords internal
#
# @return Nothing
.write_recalibr_map <- function(map, file_path) {
    if (fs::is_dir(file_path)) {
        if (!fs::dir_exists(file_path)) {
            fs::dir_create(file_path)
        }
        gen_filename <- .generate_rec_map_filename()
        file_path <- fs::path(file_path, gen_filename)
    } else {
        if (fs::path_ext(file_path) == "") {
            file_path <- fs::path_ext_set(file_path, "tsv.gz")
        }
    }
    withRestarts(
        {
            data.table::fwrite(map, file = file_path, sep = "\t", na = "")
            saved_msg <- paste(
                "Recalibration map saved to: ",
                file_path
            )
            if (getOption("ISAnalytics.verbose") == TRUE) {
                rlang::inform(saved_msg)
            }
        },
        skip_write = function() {
            skip_msg <- paste(
                "Could not write recalibration map file.",
                "Skipping."
            )
            rlang::inform(skip_msg)
        }
    )
}

#### ---- Internals for analysis functions ----####


#---- USED IN : threshold_filter ----

# @keywords internal
.check_threshold_param <- function(threshold) {
    if (is.list(threshold)) {
        # List must have names
        if (is.null(names(threshold))) {
            stop(.names_list_param_err())
        }
        all_nums <- purrr::map_lgl(threshold, function(t) {
            is.numeric(t) || is.integer(t)
        })
        if (!all(all_nums)) {
            stop(.threshold_err(), call. = FALSE)
        }
    } else {
        if (!(is.numeric(threshold) || is.integer(threshold))) {
            stop(.threshold_err(), call. = FALSE)
        }
    }
}

# @keywords internal
.check_comparators_param <- function(comparators) {
    if (is.list(comparators)) {
        # List must have names
        if (is.null(names(comparators))) {
            stop(.names_list_param_err())
        }
        chk <- purrr::map_lgl(comparators, function(comp) {
            is.character(comp) &
                all(comp %in% c("<", ">", "==", "!=", ">=", "<="))
        })
        if (!all(chk)) {
            stop(.comparators_err(), call. = FALSE)
        }
    } else {
        if (!(is.character(comparators) &
            all(comparators %in% c("<", ">", "==", "!=", ">=", "<=")))) {
            stop(.comparators_err(), call. = FALSE)
        }
    }
}

# @keywords internal
.check_cols_to_compare_param <- function(cols_to_compare) {
    if (is.list(cols_to_compare)) {
        # List must have names
        if (is.null(names(cols_to_compare))) {
            stop(.names_list_param_err())
        }
        chk <- purrr::map_lgl(cols_to_compare, function(col) {
            is.character(col)
        })
        if (!all(chk)) {
            stop(paste("`cols_to_compare` elements must be character"),
                call. = FALSE
            )
        }
    } else {
        if (!is.character(cols_to_compare)) {
            stop(paste("`cols_to_compare` must be character"), call. = FALSE)
        }
    }
}

# @keywords internal
.list_names_to_mod <- function(x, list_params) {
    ### list_params is a list with 1-3 elements with names as parameters
    ### single elements are also lists
    names_equal_list <- purrr::map_lgl(list_params, function(p) {
        all(names(p) %in% names(x))
    })
    if (!all(names_equal_list)) {
        stop(.param_names_not_in_list_err())
    }
    names_same <- purrr::reduce(list_params, function(p1, p2) {
        names_same <- all(names(p1) == names(p2))
        if (names_same == FALSE) {
            stop(.param_names_not_equal_err())
        }
    })
    return(names(list_params[[1]]))
}

# @keywords internal
.check_col_correct <- function(cols_to_compare, x, names_to_mod) {
    if (is.list(cols_to_compare)) {
        purrr::walk2(
            cols_to_compare,
            names(cols_to_compare),
            function(set, name) {
                df <- x[[name]]
                if (!all(set %in% colnames(df))) {
                    rlang::abort(.missing_user_cols_list_error(set, name))
                }
                purrr::walk(set, function(col) {
                    expr <- rlang::expr(`$`(df, !!col))
                    if (!is.numeric(rlang::eval_tidy(expr)) &&
                        !is.integer(rlang::eval_tidy(expr))) {
                        rlang::abort(.non_num_user_cols_error(
                            rlang::eval_tidy(expr)
                        ))
                    }
                })
            }
        )
    } else {
        if (is.data.frame(x)) {
            if (!all(cols_to_compare %in% colnames(x))) {
                rlang::abort(.missing_user_cols_error(
                    cols_to_compare[!cols_to_compare %in% colnames(x)]
                ))
            }
            ### Check they're numeric
            all_numeric <- purrr::map_lgl(cols_to_compare, function(col) {
                expr <- rlang::expr(`$`(x, !!col))
                is.numeric(rlang::eval_tidy(expr)) ||
                    is.integer(rlang::eval_tidy(expr))
            })
            if (!all(all_numeric)) {
                rlang::abort(
                    .non_num_user_cols_error(
                        cols_to_compare[!all_numeric]
                    )
                )
            }
        } else {
            purrr::walk(x[names_to_mod], function(ddf) {
                if (!all(cols_to_compare %in% colnames(ddf))) {
                    rlang::abort(
                        .missing_user_cols_error(
                            cols_to_compare[!cols_to_compare %in% colnames(ddf)]
                        )
                    )
                }
                purrr::walk(cols_to_compare, function(col) {
                    expr <- rlang::expr(`$`(ddf, !!col))
                    if (!is.numeric(rlang::eval_tidy(expr)) &&
                        !is.integer(rlang::eval_tidy(expr))) {
                        rlang::abort(
                            .non_num_user_cols_error(rlang::eval_tidy(expr))
                        )
                    }
                })
            })
        }
    }
}

# @keywords internal
.check_length_list_params <- function(list_params) {
    purrr::pwalk(list(
        list_params$threshold, list_params$cols_to_compare,
        list_params$comparators
    ), function(x, y, z) {
        if (length(x) != length(y) || length(x) != length(z) ||
            length(y) != length(z)) {
            stop(.diff_leng_args())
        }
    })
}

# @keywords internal
.check_length_vector_params <- function(vec_params, list_params) {
    ### Check lengths between vector params first
    vec_lengths <- purrr::map_dbl(vec_params, length)
    max <- max(vec_lengths)
    vec_lengths_ok <- purrr::map_lgl(vec_lengths, function(l) {
        l == 1 || l == max
    })
    if (any(vec_lengths_ok == FALSE)) {
        stop(.diff_leng_args(), call. = FALSE)
    }
    if (!is.null(list_params) || !purrr::is_empty(list_params)) {
        list_lengths_ok <- purrr::map_lgl(list_params, function(p) {
            l <- purrr::map_dbl(p, length)
            all(l == max)
        })
        if (any(list_lengths_ok == FALSE)) {
            stop(.diff_leng_args(), call. = FALSE)
        }
    }
}

# Threshold filter for data frames.
#
# @param x A data frame
# @param threshold Vector of threshold values
# @param cols_to_compare Vector of columns to compare
# @param comparators Vector of comparators
#
#' @importFrom purrr pmap_chr reduce
#' @importFrom rlang eval_tidy parse_expr .data
#' @importFrom magrittr `%>%`
#
# @keywords internal
# @return A data frame
.tf_data_frame <- function(x, threshold, cols_to_compare, comparators) {
    ### Check all columns in input are present
    if (is.list(threshold) || is.list(comparators) ||
        is.list(cols_to_compare)) {
        stop(.list_params_err())
    }
    .check_threshold_param(threshold)
    .check_comparators_param(comparators)
    .check_cols_to_compare_param(cols_to_compare)
    .check_col_correct(cols_to_compare, x, NULL)
    .check_length_vector_params(
        list(threshold, cols_to_compare, comparators),
        NULL
    )
    ### Obtain single expression for all columns
    conditions <- purrr::pmap_chr(
        list(
            cols_to_compare,
            comparators,
            threshold
        ),
        function(q, c, t) {
            temp <- paste(q, c, t)
            paste0(".data$", temp)
        }
    )
    single_predicate_string <- purrr::reduce(conditions, function(c1, c2) {
        paste0(c1, ", ", c2)
    })
    single_predicate_string <- paste0(
        "x %>% dplyr::filter(",
        single_predicate_string,
        ")"
    )
    filtered <- rlang::eval_tidy(rlang::parse_expr(single_predicate_string))
    return(filtered)
}

# Threshold filter for lists.
#
# @param x The list of data frames
# @param threshold Vector or named list for threshold values
# @param cols_to_compare Vector or named list for columns to
# compare
# @param comparators Vector or named list for comparators
#
#' @importFrom purrr map_lgl is_empty map pmap_chr reduce map2
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data parse_expr eval_tidy
# @keywords internal
#
# @return A list of data frames
.tf_list <- function(x, threshold, cols_to_compare, comparators) {
    ### Check list contains data frames
    contain_df <- purrr::map_lgl(x, function(el) {
        is.data.frame(el)
    })
    if (!all(contain_df)) {
        stop(paste("Elements of the list must be data frames"))
    }
    ### Check the nature of the parameters (list or vector)
    list_has_names <- !is.null(names(x))
    params <- list(
        threshold = threshold, comparators = comparators,
        cols_to_compare = cols_to_compare
    )
    param_are_lists <- purrr::map_lgl(params, is.list)
    ### With unnamed lists only vector parameters allowed
    if (any(param_are_lists) & list_has_names == FALSE) {
        stop(.unnamed_list_err())
    }
    ### Check single params
    .check_threshold_param(threshold)
    .check_cols_to_compare_param(cols_to_compare)
    .check_comparators_param(comparators)
    list_params <- params[param_are_lists]
    names_to_mod <- if (!purrr::is_empty(list_params)) {
        .list_names_to_mod(x, list_params)
    } else {
        names(x)
    }
    if (all(param_are_lists)) {
        .check_length_list_params(params)
    }
    if (any(param_are_lists == FALSE)) {
        .check_length_vector_params(params[!param_are_lists], list_params)
    }
    .check_col_correct(cols_to_compare, x, names_to_mod)
    ### Filter each
    if (!list_has_names) {
        filtered <- purrr::map(x, function(df) {
            conditions <- purrr::pmap_chr(
                list(
                    cols_to_compare,
                    comparators,
                    threshold
                ),
                function(q, c, t) {
                    temp <- paste(q, c, t)
                    paste0(".data$", temp)
                }
            )
            single_predicate_string <- purrr::reduce(
                conditions,
                function(c1, c2) {
                    paste0(c1, ", ", c2)
                }
            )
            single_predicate_string <- paste0(
                "df %>% dplyr::filter(",
                single_predicate_string,
                ")"
            )
            rlang::eval_tidy(
                rlang::parse_expr(single_predicate_string)
            )
        })
        return(filtered)
    }
    filtered <- purrr::map2(x, names(x), function(df, nm) {
        if (!nm %in% names_to_mod) {
            df
        } else {
            col_set <- if (is.list(cols_to_compare)) {
                cols_to_compare[[nm]]
            } else {
                cols_to_compare
            }
            comp <- if (is.list(comparators)) {
                comparators[[nm]]
            } else {
                comparators
            }
            thr <- if (is.list(threshold)) {
                threshold[[nm]]
            } else {
                threshold
            }
            conditions <- purrr::pmap_chr(
                list(
                    col_set,
                    comp,
                    thr
                ),
                function(q, c, t) {
                    temp <- paste(q, c, t)
                    paste0(".data$", temp)
                }
            )
            single_predicate_string <- purrr::reduce(
                conditions,
                function(c1, c2) {
                    paste0(c1, ", ", c2)
                }
            )
            single_predicate_string <- paste0(
                "df %>% dplyr::filter(",
                single_predicate_string,
                ")"
            )
            rlang::eval_tidy(
                rlang::parse_expr(single_predicate_string)
            )
        }
    })
    return(filtered)
}

#---- USED IN : CIS_grubbs ----
# Internal computation of CIS_grubbs test
.cis_grubb_calc <- function(x,
    refgenes,
    grubbs_flanking_gene_bp,
    threshold_alpha) {
    df_by_gene <- if (requireNamespace("psych", quietly = TRUE)) {
        x %>%
            dplyr::group_by(
                .data$GeneName,
                .data$GeneStrand,
                .data$chr
            ) %>%
            dplyr::summarise(
                n_IS_perGene = dplyr::n_distinct(
                    .data$integration_locus
                ),
                min_bp_integration_locus =
                    min(.data$integration_locus),
                max_bp_integration_locus =
                    max(.data$integration_locus),
                IS_span_bp = (max(.data$integration_locus) -
                    min(.data$integration_locus)),
                avg_bp_integration_locus =
                    mean(.data$integration_locus),
                median_bp_integration_locus =
                    stats::median(.data$integration_locus),
                distinct_orientations = dplyr::n_distinct(.data$strand),
                describe = list(tibble::as_tibble(
                    psych::describe(.data$integration_locus)
                )),
                .groups = "drop"
            ) %>%
            tidyr::unnest(cols = .data$describe)
    } else {
        x %>%
            dplyr::group_by(
                .data$GeneName,
                .data$GeneStrand,
                .data$chr
            ) %>%
            dplyr::summarise(
                n_IS_perGene = dplyr::n_distinct(
                    .data$integration_locus
                ),
                min_bp_integration_locus =
                    min(.data$integration_locus),
                max_bp_integration_locus =
                    max(.data$integration_locus),
                IS_span_bp = (max(.data$integration_locus) -
                    min(.data$integration_locus)),
                avg_bp_integration_locus =
                    mean(.data$integration_locus),
                median_bp_integration_locus =
                    stats::median(.data$integration_locus),
                distinct_orientations = dplyr::n_distinct(.data$strand),
                .groups = "drop"
            )
    }
    df_bygene_withannotation <- df_by_gene %>%
        dplyr::inner_join(refgenes, by = c(
            "chr" = "chrom",
            "GeneStrand" = "strand",
            "GeneName" = "name2"
        )) %>%
        dplyr::select(c(
            dplyr::all_of(colnames(df_by_gene)),
            .data$average_TxLen
        ))
    n_elements <- nrow(df_bygene_withannotation)

    df_bygene_withannotation <- df_bygene_withannotation %>%
        dplyr::mutate(
            geneIS_frequency_byHitIS = .data$n_IS_perGene / n_elements
        )

    ### Grubbs test
    ### --- Gene Frequency
    df_bygene_withannotation <- df_bygene_withannotation %>%
        dplyr::mutate(
            raw_gene_integration_frequency =
                .data$n_IS_perGene / .data$average_TxLen,
            integration_frequency_withtolerance = .data$n_IS_perGene /
                (.data$average_TxLen + grubbs_flanking_gene_bp) * 1000,
            minus_log2_integration_freq_withtolerance =
                -log(x = .data$integration_frequency_withtolerance, base = 2)
        )
    ### --- z score
    z_mlif <- function(x) {
        sqrt((n_elements * (n_elements - 2) * x^2) /
            (((n_elements - 1)^2) - (n_elements * x^2)))
    }
    df_bygene_withannotation <- df_bygene_withannotation %>%
        dplyr::mutate(
            zscore_minus_log2_int_freq_tolerance =
                scale(-log(
                    x = .data$integration_frequency_withtolerance,
                    base = 2
                ))[, 1],
            neg_zscore_minus_log2_int_freq_tolerance =
                -.data$zscore_minus_log2_int_freq_tolerance,
            t_z_mlif = z_mlif(
                .data$neg_zscore_minus_log2_int_freq_tolerance
            )
        )
    ### --- tdist
    t_dist_2t <- function(x, deg) {
        return((1 - stats::pt(x, deg)) * 2)
    }
    df_bygene_withannotation <- df_bygene_withannotation %>%
        dplyr::mutate(
            tdist2t = t_dist_2t(.data$t_z_mlif, n_elements - 2),
            tdist_pt = stats::pt(
                q = .data$t_z_mlif,
                df = n_elements - 2
            ),
            tdist_bonferroni_default = ifelse(
                .data$tdist2t * n_elements > 1, 1,
                .data$tdist2t * n_elements
            ),
            tdist_bonferroni = stats::p.adjust(
                .data$tdist2t,
                method = "bonferroni",
                n = length(.data$tdist2t)
            ),
            tdist_fdr = stats::p.adjust(
                .data$tdist2t,
                method = "fdr",
                n = length(.data$tdist2t)
            ),
            tdist_benjamini = stats::p.adjust(
                .data$tdist2t,
                method = "BY",
                n = length(.data$tdist2t)
            )
        )
    df_bygene_withannotation <- df_bygene_withannotation %>%
        dplyr::mutate(
            tdist_positive_and_corrected =
                ifelse(
                    (.data$tdist_bonferroni_default < threshold_alpha &
                        .data$neg_zscore_minus_log2_int_freq_tolerance > 0),
                    .data$tdist_bonferroni_default,
                    NA
                ),
            tdist_positive = ifelse(
                (.data$tdist2t < threshold_alpha &
                    .data$neg_zscore_minus_log2_int_freq_tolerance > 0),
                .data$tdist2t,
                NA
            )
        )
    EM_correction_N <- length(
        df_bygene_withannotation$tdist_positive[
            !is.na(df_bygene_withannotation$tdist_positive)
        ]
    )
    df_bygene_withannotation <- df_bygene_withannotation %>%
        dplyr::mutate(
            tdist_positive_and_correctedEM =
                ifelse(
                    (.data$tdist2t * EM_correction_N <
                        threshold_alpha &
                        .data$neg_zscore_minus_log2_int_freq_tolerance > 0),
                    .data$tdist2t * EM_correction_N,
                    NA
                )
        )
    return(df_bygene_withannotation)
}

#---- USED IN : CIS_volcano_plot ----
#' @importFrom rlang arg_match abort inform current_env
#' @importFrom utils read.delim
.load_onco_ts_genes <- function(onco_db_file,
    tumor_suppressors_db_file,
    species) {
    onco_db <- if (onco_db_file == "proto_oncogenes") {
        utils::data("proto_oncogenes", envir = rlang::current_env())
        rlang::current_env()$proto_oncogenes
    } else {
        if (!file.exists(onco_db_file)) {
            onco_db_not_found <- c(paste("`onco_db_file` was not found"),
                x = "Did you provide the correct path?"
            )
            rlang::abort(onco_db_not_found)
        }
        utils::read.delim(
            file = onco_db_file, header = TRUE,
            fill = TRUE, sep = "\t",
            check.names = FALSE,
            na.strings = c("NONE", "NA", "NULL", "NaN", "")
        )
    }
    tumsup_db <- if (tumor_suppressors_db_file == "tumor_suppressors") {
        utils::data("tumor_suppressors", envir = rlang::current_env())
        rlang::current_env()$tumor_suppressors
    } else {
        if (!file.exists(tumor_suppressors_db_file)) {
            tum_db_not_found <- c("`tumor_suppressors_db_file` was not found",
                x = "Did you provide the correct path?"
            )
            rlang::abort(tum_db_not_found)
        }
        utils::read.delim(
            file = tumor_suppressors_db_file,
            header = TRUE,
            fill = TRUE, sep = "\t",
            check.names = FALSE,
            na.strings = c("NONE", "NA", "NULL", "NaN", "")
        )
    }
    specie <- rlang::arg_match(species, values = c("human", "mouse", "all"))
    specie_name <- switch(specie,
        "human" = "Homo sapiens (Human)",
        "mouse" = "Mus musculus (Mouse)",
        "all" = c(
            "Homo sapiens (Human)",
            "Mus musculus (Mouse)"
        )
    )
    if (getOption("ISAnalytics.verbose") == TRUE) {
        load_onco_msg <- c(
            "Loading annotated genes -  species selected: ",
            paste(c(specie_name), collapse = ", ")
        )
        rlang::inform(load_onco_msg)
    }
    # Filter and merge
    onco_df <- .filter_db(onco_db, specie_name, "OncoGene")
    tumsup_df <- .filter_db(tumsup_db, specie_name, "TumorSuppressor")
    oncots_df_to_use <- .merge_onco_tumsup(onco_df, tumsup_df)
    if (getOption("ISAnalytics.verbose") == TRUE) {
        rlang::inform(c(
            "Loading annotated genes -  done"
        ))
    }
    return(oncots_df_to_use)
}

#' @importFrom rlang .data
#' @importFrom magrittr `%>%`
.filter_db <- function(onco_db, species, output_col_name) {
    filtered_db <- onco_db %>%
        dplyr::filter(.data$Organism %in% species) %>%
        dplyr::filter(.data$Status == "reviewed" &
            !is.na(.data$`Gene names`)) %>%
        dplyr::select(.data$`Gene names`) %>%
        dplyr::distinct()

    filtered_db <- filtered_db %>%
        dplyr::mutate(`Gene names` = stringr::str_split(
            .data$`Gene names`, " "
        )) %>%
        tidyr::unnest("Gene names") %>%
        dplyr::mutate({{ output_col_name }} := .data$`Gene names`) %>%
        dplyr::rename("GeneName" = "Gene names")

    return(filtered_db)
}

#' @importFrom rlang .data
#' @importFrom magrittr `%>%`
.merge_onco_tumsup <- function(onco_df, tum_df) {
    onco_tumsup <- onco_df %>%
        dplyr::full_join(tum_df, by = "GeneName") %>%
        dplyr::mutate(Onco1_TS2 = ifelse(
            !is.na(.data$OncoGene) & !is.na(.data$TumorSuppressor),
            yes = 3,
            no = ifelse(!is.na(.data$OncoGene),
                yes = 1,
                no = ifelse(!is.na(.data$TumorSuppressor),
                    yes = 2,
                    no = NA
                )
            )
        )) %>%
        dplyr::select(-c("OncoGene", "TumorSuppressor")) %>%
        dplyr::distinct()
    return(onco_tumsup)
}


#---- USED IN : outliers_by_pool_fragments ----
# Actual computation of statistical test on pre-filtered metadata (no NA values)
#' @importFrom rlang inform
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply bpstop
#' @importFrom tibble add_column
#' @importFrom purrr reduce
#' @importFrom magrittr `%>%`
.pool_frag_calc <- function(meta,
    key,
    by_pool,
    normality_test,
    normality_threshold,
    pool_col,
    min_samples_per_pool,
    log2) {
    log2_removed_report <- NULL
    if (log2) {
        ## discard values < = 0 and report removed
        if (getOption("ISAnalytics.verbose") == TRUE) {
            rlang::inform("Log2 transformation, removing values <= 0")
        }
        old_meta <- meta
        meta <- meta %>%
            dplyr::filter(dplyr::across(dplyr::all_of(key), ~ .x > 0))
        log2_removed_report <- old_meta %>%
            dplyr::anti_join(meta, by = "CompleteAmplificationID") %>%
            dplyr::select(
                dplyr::all_of(pool_col),
                .data$CompleteAmplificationID,
                dplyr::all_of(key)
            ) %>%
            dplyr::distinct()
    }
    if (by_pool) {
        # Group by pool
        grouped <- meta %>%
            dplyr::group_by(dplyr::across({{ pool_col }})) %>%
            dplyr::add_tally(name = "sample_count")
        split <- grouped %>%
            dplyr::filter(.data$sample_count >= min_samples_per_pool) %>%
            dplyr::select(-.data$sample_count) %>%
            dplyr::group_split()
        p <- if (.Platform$OS.type == "windows") {
            BiocParallel::SnowParam(
                tasks = length(split),
                progressbar = getOption("ISAnalytics.verbose"),
                exportglobals = FALSE
            )
        } else {
            BiocParallel::MulticoreParam(
                tasks = length(split),
                progressbar = getOption("ISAnalytics.verbose"),
                exportglobals = FALSE
            )
        }
        test_res <- BiocParallel::bplapply(
            X = split,
            FUN = .process_pool_frag,
            key = key,
            normality_test = normality_test,
            normality_threshold = normality_threshold,
            log2 = log2
        )
        BiocParallel::bpstop(p)
        test_res <- purrr::reduce(test_res, dplyr::bind_rows)
        test_res <- test_res %>% tibble::add_column(processed = TRUE)
        # Groups not processed because not enough samples
        non_proc <- grouped %>%
            dplyr::ungroup() %>%
            dplyr::filter(.data$sample_count < min_samples_per_pool) %>%
            dplyr::select(-.data$sample_count) %>%
            dplyr::mutate(processed = FALSE)
        final <- test_res %>% dplyr::bind_rows(non_proc)
        # Rows not processed because log2 requested and value <= 0
        if (exists("old_meta")) {
            log2_removed <- old_meta %>%
                dplyr::anti_join(meta, by = "CompleteAmplificationID") %>%
                dplyr::mutate(processed = FALSE)
            final <- final %>% dplyr::bind_rows(log2_removed)
        }
        return(list(
            metadata = final,
            removed_zeros = log2_removed_report,
            non_proc_samples = grouped %>%
                dplyr::ungroup() %>%
                dplyr::filter(.data$sample_count < min_samples_per_pool)
        ))
    } else {
        test_res <- .process_pool_frag(
            chunk = meta,
            key = key,
            normality_test = normality_test,
            normality_threshold = normality_threshold,
            log2 = log2
        )
        final <- test_res %>% dplyr::mutate(processed = TRUE)
        # Rows not processed because log2 requested and value <= 0
        if (exists("old_meta")) {
            log2_removed <- old_meta %>%
                dplyr::filter(dplyr::across(dplyr::all_of(key), ~ .x <= 0)) %>%
                dplyr::mutate(processed = FALSE)
            final <- final %>% dplyr::bind_rows(log2_removed)
        }
        return(list(
            metadata = final,
            removed_zeros = log2_removed_report
        ))
    }
}

# Only computes for each key the actual calculations (either on pool chunk
# or entire df)
#' @importFrom purrr map reduce flatten
#' @importFrom dplyr bind_cols
.process_pool_frag <- function(chunk, key, normality_test,
    normality_threshold,
    log2) {
    res <- purrr::map(key, ~ .tests_pool_frag(
        x = chunk[[.x]],
        suffix = .x,
        normality_test = normality_test,
        normality_threshold = normality_threshold,
        log2 = log2
    ))
    res <- purrr::reduce(
        list(og = chunk, purrr::flatten(res)),
        dplyr::bind_cols
    )
    res
}

# Computation on single vectors (columns)
#' @importFrom stats dt shapiro.test
#' @importFrom tibble tibble
.tests_pool_frag <- function(x,
    suffix, normality_test,
    normality_threshold,
    log2) {
    if (log2) {
        x <- log2(x)
    }
    norm <- NA
    zscore <- NA
    tstudent <- NA
    tdist <- NA
    if (normality_test) {
        withCallingHandlers(
            {
                withRestarts(
                    {
                        shapiro_test <- stats::shapiro.test(x)
                        norm <- shapiro_test$p.value >= normality_threshold
                        if (norm) {
                            zscore <- scale(x)
                            zscore <- zscore[, 1]
                            count <- length(x)
                            tstudent <- sqrt(
                                ((count * (count - 2) * (zscore)^2) /
                                    ((count - 1)^2 - count * (zscore)^2))
                            )
                            tdist <- stats::dt(tstudent, df = count - 2)
                        }
                    },
                    test_err = function() {
                        rlang::inform(
                            c("Unable to perform normality test, skipping")
                        )
                    }
                )
            },
            error = function(cnd) {
                rlang::inform(cnd$message)
                invokeRestart("test_err")
            }
        )
    } else {
        zscore <- scale(x)
        zscore <- zscore[, 1]
        count <- length(x)
        tstudent <- sqrt(
            ((count * (count - 2) * (zscore)^2) /
                ((count - 1)^2 - count * (zscore)^2))
        )
        tdist <- stats::dt(tstudent, df = count - 2)
    }
    results <- if (log2) {
        tibble::tibble(
            "log2_{suffix}" := x,
            "normality_{suffix}" := norm,
            "zscore_{suffix}" := zscore,
            "tstudent_{suffix}" := tstudent,
            "tdist_{suffix}" := tdist
        )
    } else {
        tibble::tibble(
            "normality_{suffix}" := norm,
            "zscore_{suffix}" := zscore,
            "tstudent_{suffix}" := tstudent,
            "tdist_{suffix}" := tdist
        )
    }
    return(results)
}

# Flag logic to use on a single key component
# returns a logical vector
#' @importFrom purrr pmap_lgl
.flag_cond_outliers_pool_frag <- function(proc,
    tdist,
    zscore,
    outlier_threshold) {
    # basic flag condition is tdist < out.thresh AND zscore < 0
    purrr::pmap_lgl(list(proc, tdist, zscore), function(p, t, z) {
        if (p == FALSE) {
            return(FALSE)
        }
        if (is.na(t) | is.na(z)) {
            return(FALSE)
        }
        return(t < outlier_threshold & z < 0)
    })
}

# Obtain and apply global logic on multiple columns based on the single
# key logic
#' @importFrom tibble tibble
.apply_flag_logic <- function(..., logic) {
    all_flags <- tibble::tibble(...)
    partial <- NULL
    index <- 1
    for (op in logic) {
        switch(op,
            "AND" = {
                if (is.null(partial)) {
                    partial <- all_flags[[1]] & all_flags[[2]]
                    index <- 3
                } else {
                    partial <- partial & all_flags[[index]]
                    index <- index + 1
                }
            },
            "OR" = {
                if (is.null(partial)) {
                    partial <- all_flags[[1]] | all_flags[[2]]
                    index <- 3
                } else {
                    partial <- partial | all_flags[[index]]
                    index <- index + 1
                }
            },
            "XOR" = {
                if (is.null(partial)) {
                    partial <- xor(all_flags[[1]], all_flags[[2]])
                    index <- 3
                } else {
                    partial <- xor(partial, all_flags[[index]])
                    index <- index + 1
                }
            },
            "NAND" = {
                if (is.null(partial)) {
                    partial <- !(all_flags[[1]] & all_flags[[2]])
                    index <- 3
                } else {
                    partial <- !(partial & all_flags[[index]])
                    index <- index + 1
                }
            },
            "NOR" = {
                if (is.null(partial)) {
                    partial <- !(all_flags[[1]] | all_flags[[2]])
                    index <- 3
                } else {
                    partial <- !(partial | all_flags[[index]])
                    index <- index + 1
                }
            },
            "XNOR" = {
                if (is.null(partial)) {
                    partial <- !xor(all_flags[[1]], all_flags[[2]])
                    index <- 3
                } else {
                    partial <- !xor(partial, all_flags[[index]])
                    index <- index + 1
                }
            }
        )
    }
    return(partial)
}

#---- USED IN : HSC_population_size_estimate ----
# Calculates population estimates (all)
#' @importFrom purrr map_lgl detect_index
#' @importFrom tidyr pivot_wider
.estimate_pop <- function(df,
    seqCount_column,
    timepoint_column,
    annotation_cols,
    stable_timepoints,
    subject) {
    # --- STABLE TIMEPOINTS?
    found_stable <- purrr::map_lgl(
        stable_timepoints,
        ~ .x %in%
            as.numeric(df[[timepoint_column]])
    )
    first_stable_index <- if (length(found_stable) > 0) {
        purrr::detect_index(found_stable, ~ .x == TRUE)
    } else {
        0
    }
    stable_tps <- if (first_stable_index > 0) {
        TRUE
    } else {
        FALSE
    }
    # --- OBTAIN MATRIX (ALL TPs)
    matrix_desc <- df %>%
        dplyr::mutate(bin = 1) %>%
        dplyr::select(-.data[[seqCount_column]]) %>%
        tidyr::pivot_wider(
            names_from = c(
                "CellType",
                "Tissue",
                timepoint_column
            ),
            values_from = .data$bin,
            names_sort = TRUE,
            values_fill = 0
        ) %>%
        dplyr::select(-c(mandatory_IS_vars(), annotation_cols)) %>%
        as.matrix()
    # --- OBTAIN MATRIX (STABLE TPs)
    patient_slice_stable <- if (first_stable_index > 0) {
        first_stable <- stable_timepoints[first_stable_index]
        df %>%
            dplyr::filter(as.numeric(.data[[timepoint_column]]) >=
                first_stable) %>%
            dplyr::mutate(bin = 1) %>%
            dplyr::select(-.data[[seqCount_column]]) %>%
            tidyr::pivot_wider(
                names_from = c(
                    "CellType",
                    "Tissue",
                    timepoint_column
                ),
                values_from = .data$bin,
                names_sort = TRUE,
                values_fill = 0
            ) %>%
            dplyr::select(-c(mandatory_IS_vars(), annotation_cols)) %>%
            as.matrix()
    } else {
        NULL
    }

    timecaptures <- length(colnames(matrix_desc))
    cols_estimate_mcm <- c("Model", "abundance", "stderr")
    # --- ESTIMATES
    ## models0
    ### Estimate abundance without bias correction nor het
    #### For all tps
    models0_all <- .closed_m0_est(
        matrix_desc = matrix_desc,
        timecaptures = timecaptures,
        cols_estimate_mcm = cols_estimate_mcm,
        subject = subject, stable = FALSE
    )
    #### For the slice (first stable - last tp)
    models0_stable <- if (!is.null(patient_slice_stable) &&
        ncol(patient_slice_stable) > 1) {
        .closed_m0_est(
            matrix_desc = patient_slice_stable,
            timecaptures = length(colnames(patient_slice_stable)),
            cols_estimate_mcm = cols_estimate_mcm,
            subject = subject, stable = TRUE
        )
    } else {
        NULL
    }
    ## BC models
    ### Estimate abundance without bias correction nor het
    models_bc_all <- .closed_mthchaobc_est(
        matrix_desc = matrix_desc,
        timecaptures = timecaptures,
        cols_estimate_mcm = cols_estimate_mcm,
        subject = subject, stable = FALSE
    )
    models_bc_stable <- if (!is.null(patient_slice_stable) &&
        ncol(patient_slice_stable) > 1) {
        .closed_mthchaobc_est(
            matrix_desc = patient_slice_stable,
            timecaptures = length(colnames(patient_slice_stable)),
            cols_estimate_mcm = cols_estimate_mcm,
            subject = subject, stable = TRUE
        )
    } else {
        NULL
    }
    ## Consecutive timepoints
    ### Model m0 bc
    estimate_consecutive_m0 <- if (ncol(matrix_desc) > 1) {
        .consecutive_m0bc_est(
            matrix_desc = matrix_desc,
            cols_estimate_mcm = cols_estimate_mcm,
            subject = subject
        )
    } else {
        NULL
    }
    ### Model Mth
    estimate_consecutive_mth <- if (stable_tps & ncol(matrix_desc) > 2) {
        # - Note: pass the whole matrix, not only stable slice
        .consecutive_mth_est(
            matrix_desc = matrix_desc,
            cols_estimate_mcm = cols_estimate_mcm,
            subject = subject
        )
    } else {
        NULL
    }
    results <- models0_all %>%
        dplyr::bind_rows(models0_stable) %>%
        dplyr::bind_rows(models_bc_all) %>%
        dplyr::bind_rows(models_bc_stable) %>%
        dplyr::bind_rows(estimate_consecutive_m0) %>%
        dplyr::bind_rows(estimate_consecutive_mth)
    results
}

#' @importFrom Rcapture closedp.0
#' @importFrom purrr flatten_chr
#' @importFrom stringr str_split
#' @importFrom tibble as_tibble add_column
#' @importFrom dplyr select all_of
.closed_m0_est <- function(matrix_desc, timecaptures, cols_estimate_mcm,
    subject, stable) {
    ### Estimate abundance without bias correction nor het
    models0 <- Rcapture::closedp.0(matrix_desc,
        dfreq = FALSE,
        dtype = "hist",
        t = timecaptures,
        t0 = timecaptures,
        neg = FALSE
    )
    splitted_col_nms1 <- purrr::flatten_chr(stringr::str_split(
        colnames(matrix_desc)[1],
        pattern = "_"
    ))
    splitted_col_nmsn <- purrr::flatten_chr(stringr::str_split(
        colnames(matrix_desc)[ncol(matrix_desc)],
        pattern = "_"
    ))
    tp_string <- if (stable) {
        "Stable"
    } else {
        "All"
    }
    estimate_results <- tibble::as_tibble(models0$results,
        rownames = "Model"
    ) %>%
        dplyr::select(dplyr::all_of(cols_estimate_mcm)) %>%
        tibble::add_column(
            SubjectID = subject,
            Timepoints = tp_string,
            CellType = splitted_col_nms1[1],
            Tissue = splitted_col_nms1[2],
            TimePoint_from = as.numeric(splitted_col_nms1[3]),
            TimePoint_to = as.numeric(splitted_col_nmsn[3]),
            ModelType = "ClosedPopulation",
            ModelSetUp = "models0"
        )
    estimate_results
}

#' @importFrom Rcapture closedp.bc
#' @importFrom purrr flatten_chr
#' @importFrom stringr str_split
#' @importFrom tibble as_tibble add_column
#' @importFrom dplyr select all_of
.closed_mthchaobc_est <- function(matrix_desc, timecaptures, cols_estimate_mcm,
    subject, stable) {
    mthchaobc <- Rcapture::closedp.bc(matrix_desc,
        dfreq = FALSE,
        dtype = "hist",
        t = timecaptures,
        t0 = timecaptures,
        m = c("M0", "Mt", "Mth")
    )
    splitted_col_nms1 <- purrr::flatten_chr(stringr::str_split(
        colnames(matrix_desc)[1],
        pattern = "_"
    ))
    splitted_col_nmsn <- purrr::flatten_chr(stringr::str_split(
        colnames(matrix_desc)[ncol(matrix_desc)],
        pattern = "_"
    ))
    tp_string <- if (stable) {
        "Stable"
    } else {
        "All"
    }
    estimate_results_chaobc <- tibble::as_tibble(mthchaobc$results,
        rownames = "Model"
    ) %>%
        dplyr::select(dplyr::all_of(cols_estimate_mcm)) %>%
        tibble::add_column(
            SubjectID = subject,
            Timepoints = tp_string,
            CellType = splitted_col_nms1[1],
            Tissue = splitted_col_nms1[2],
            TimePoint_from = as.numeric(splitted_col_nms1[3]),
            TimePoint_to = as.numeric(splitted_col_nmsn[3]),
            ModelType = "ClosedPopulation",
            ModelSetUp = "mthchaobc"
        )
    estimate_results_chaobc
}

#' @importFrom Rcapture closedp.bc
#' @importFrom purrr flatten_chr map2_dfr
#' @importFrom stringr str_split
#' @importFrom tibble as_tibble add_column
#' @importFrom dplyr select all_of
.consecutive_m0bc_est <- function(matrix_desc, cols_estimate_mcm, subject) {
    ## Consecutive timepoints
    indexes <- seq(from = 1, to = ncol(matrix_desc) - 1, by = 1)
    purrr::map2_dfr(
        indexes,
        indexes + 1,
        function(t1, t2) {
            sub_matrix <- matrix_desc[, seq(from = t1, to = t2, by = 1)]
            splitted_col_nms1 <- purrr::flatten_chr(stringr::str_split(
                colnames(sub_matrix)[1],
                pattern = "_"
            ))
            splitted_col_nmsn <- purrr::flatten_chr(stringr::str_split(
                colnames(sub_matrix)[ncol(sub_matrix)],
                pattern = "_"
            ))
            patient_trend_M0 <- Rcapture::closedp.bc(sub_matrix,
                dfreq = FALSE,
                dtype = "hist",
                t = 2,
                t0 = 2,
                m = c("M0")
            )
            results <- tibble::as_tibble(patient_trend_M0$results,
                rownames = "Model"
            ) %>%
                dplyr::select(dplyr::all_of(cols_estimate_mcm)) %>%
                tibble::add_column(
                    SubjectID = subject,
                    Timepoints = "Consecutive",
                    CellType = splitted_col_nms1[1],
                    Tissue = splitted_col_nms1[2],
                    TimePoint_from = as.numeric(splitted_col_nms1[3]),
                    TimePoint_to = as.numeric(splitted_col_nmsn[3]),
                    ModelType = "ClosedPopulation",
                    ModelSetUp = "models0BC"
                )
            results
        }
    )
}

#' @importFrom Rcapture closedp.bc
#' @importFrom purrr flatten_chr map2_dfr
#' @importFrom stringr str_split
#' @importFrom tibble as_tibble add_column
#' @importFrom dplyr select all_of
.consecutive_mth_est <- function(matrix_desc, cols_estimate_mcm, subject) {
    indexes_s <- seq(from = 1, to = ncol(matrix_desc) - 2, by = 1)
    purrr::map2_dfr(
        indexes_s,
        indexes_s + 2,
        function(t1, t2) {
            sub_matrix <- matrix_desc[, seq(from = t1, to = t2, by = 1)]
            splitted_col_nms1 <- purrr::flatten_chr(stringr::str_split(
                colnames(sub_matrix)[1],
                pattern = "_"
            ))
            splitted_col_nmsn <- purrr::flatten_chr(stringr::str_split(
                colnames(sub_matrix)[ncol(sub_matrix)],
                pattern = "_"
            ))
            patient_trend_Mth <- Rcapture::closedp.bc(sub_matrix,
                dfreq = FALSE,
                dtype = "hist",
                t = 3, t0 = 3,
                m = c("Mth"),
                h = "Chao"
            )
            result <- tibble::as_tibble(patient_trend_Mth$results,
                rownames = "Model"
            ) %>%
                dplyr::select(dplyr::all_of(cols_estimate_mcm)) %>%
                tibble::add_column(
                    SubjectID = subject,
                    Timepoints = "Consecutive",
                    CellType = splitted_col_nms1[1],
                    Tissue = splitted_col_nms1[2],
                    TimePoint_from = as.numeric(splitted_col_nms1[3]),
                    TimePoint_to = as.numeric(splitted_col_nmsn[3]),
                    ModelType = "ClosedPopulation",
                    ModelSetUp = "modelMTHBC"
                )
            result
        }
    )
}

#---- USED IN : is_sharing ----
## Internal to find absolute shared number of is between an arbitrary
## number of groups
## Dots are group names in sharing df (actually a data.table)
.find_in_common <- function(..., lookup_tbl, keep_genomic_coord) {
    groups <- as.list(...)
    in_common <- purrr::pmap(groups, function(...) {
        grps <- list(...)
        filt <- lookup_tbl[group_id %chin% grps, ]
        common <- purrr::reduce(filt$is, function(l, r) {
            l[r, on = mandatory_IS_vars(), nomatch = 0]
        })
        common
    })
    if (!keep_genomic_coord) {
        return(purrr::map(in_common, ~ nrow(.x)))
    }
    return(list(purrr::map(in_common, ~ nrow(.x)), in_common))
}

## Internal, to use on each row of the combinations df.
## Expands the row with all its permutations keeping same absolute shared is
## and counts if present
## - ... : row passed as a list
## - g_names: names of the groups (g1, g2...)
## - counts: are counts present? TRUE/FALSE
.sh_row_permut <- function(..., g_names, counts) {
    og_row <- list(...)
    coord <- if ("is_coord" %in% names(og_row)) {
        og_row$is_coord <- list(og_row$is_coord)
        TRUE
    } else {
        FALSE
    }
    ids <- unlist(og_row[g_names])
    og_row <- data.table::setDT(og_row)
    # If row of all equal elements no need for permutations
    if (length(unique(ids)) == 1) {
        return(og_row)
    }
    # If elements are different
    shared_is <- og_row$shared
    if (counts) {
        count_union <- og_row$count_union
    }
    perm <- gtools::permutations(
        n = length(g_names),
        r = length(g_names),
        v = g_names,
        set = TRUE,
        repeats.allowed = FALSE
    )
    colnames(perm) <- g_names
    perm <- data.table::setDT(as.data.frame(perm))
    perm[, shared := shared_is]
    if (coord) {
        perm[, is_coord := list(og_row$is_coord)]
    }
    if (counts) {
        for (g in g_names) {
            count_col <- paste0("count_", g)
            perm[, c(count_col) := paste0("count_", get(g))]
        }
        perm[, count_union := count_union]
    }
    sub_with_val <- function(val) {
        unlist(purrr::map(val, ~ og_row[[.x]]))
    }
    for (g in g_names) {
        perm[, c(g) := .(sub_with_val(get(g)))]
        if (counts) {
            count_col <- paste0("count_", g)
            perm[, c(count_col) := .(sub_with_val(get(count_col)))]
        }
    }
    return(perm)
}

# Counts the number of integrations in the union of all groups of a row
.count_group_union <- function(..., col_groups, lookup_tbl) {
    dots <- list(...)
    if ("is_coord" %in% names(dots)) {
        dots$is_coord <- list(dots$is_coord)
    }
    row <- data.table::setDT(dots)
    groups_in_row <- row[1, ..col_groups]
    sub_lookup <- lookup_tbl[group_id %chin% groups_in_row]
    count_union <- nrow(purrr::reduce(sub_lookup$is, data.table::funion))
    row[, count_union := count_union]
    return(row)
}

# Obtains lookup table for groups
.sh_obtain_lookup <- function(key, df) {
    temp <- df %>%
        dplyr::select(dplyr::all_of(c(key, mandatory_IS_vars()))) %>%
        dplyr::group_by(dplyr::across({{ key }})) %>%
        dplyr::distinct(dplyr::across(dplyr::all_of(mandatory_IS_vars())),
            .keep_all = TRUE
        ) %>%
        tidyr::unite(col = "group_id", dplyr::all_of(key))
    temp <- data.table::setDT(temp)
    temp <- temp[, .(is = list(.SD)), by = "group_id"]
    return(temp)
}

# Obtains truth table for venn diagrams. To apply to each row with pmap.
# - groups : g1, g2, g3...
# - lookup : either list (mult/mult) or df
.sh_truth_tbl_venn <- function(..., lookup, groups) {
    row <- list(...)
    if (!is.data.frame(lookup)) {
        # If lookup is a list and not single df
        retrieve_is <- function(group_name) {
            label <- row[[group_name]] # a string
            retrieved <- (lookup[[group_name]])[group_id == label]
            retrieved[, group_id := paste0(group_id, "(", group_name, ")")]
            retrieved <- retrieved %>%
                tidyr::unnest(.data$is) %>%
                tidyr::unite(
                    col = "int_id",
                    dplyr::all_of(mandatory_IS_vars())
                ) %>%
                data.table::setDT()
        }
        retrieved_iss <- purrr::map(groups, retrieve_is)
        retrieved_iss <- data.table::rbindlist(retrieved_iss)
        retrieved_iss[, observed := TRUE]
        truth_tbl <- data.table::dcast(retrieved_iss, int_id ~ group_id,
            value.var = "observed",
            fill = FALSE
        )
        return(truth_tbl)
    } else {
        labels <- unlist(row[groups])
        retrieved <- lookup[group_id %in% labels]
        retrieved <- retrieved %>%
            tidyr::unnest(.data$is) %>%
            tidyr::unite(
                col = "int_id",
                dplyr::all_of(mandatory_IS_vars())
            ) %>%
            data.table::setDT()
        retrieved[, observed := TRUE]
        truth_tbl <- data.table::dcast(retrieved, int_id ~ group_id,
            value.var = "observed",
            fill = FALSE
        )
        return(truth_tbl)
    }
}

## Computes sharing table for single input df and single key
## -key: name of the columns to do a group by
## -minimal: true or false, if true computes ONLY combinations and not
## all permutations
## -n_comp: number of comparisons (2-way sharing, 3-way sharing...).
## Should be less or equal to the distinct number of groups
## -is_count: keep the counts in the table?
## -rel_sharing: compute relative sharing? (on g_i & on_union)
## -include_self_comp: include rows with the same group for each comparison?
## useful for heatmaps. Sharing for this rows is always 100%
.sharing_singledf_single_key <- function(df, key, minimal, n_comp,
    is_count, rel_sharing,
    include_self_comp,
    keep_genomic_coord,
    venn) {
    temp <- .sh_obtain_lookup(key, df)
    # Check n and k
    if (nrow(temp) < n_comp) {
        n_comp <- nrow(temp)
        if (getOption("ISAnalytics.verbose") == TRUE) {
            warn_msg <- c("Number of requested comparisons too big",
                i = paste(
                    "The number of requested comparisons",
                    "is greater than the number of groups.",
                    "Reducing comparisons to the biggest value",
                    "allowed"
                )
            )
            rlang::inform(warn_msg)
        }
    }
    if (getOption("ISAnalytics.verbose") == TRUE) {
        rlang::inform("Calculating combinations...")
    }
    group_comb <- gtools::combinations(
        n = length(temp$group_id),
        r = n_comp,
        v = temp$group_id,
        set = TRUE,
        repeats.allowed = FALSE
    )
    cols <- paste0("g", seq_len(ncol(group_comb)))
    colnames(group_comb) <- cols
    group_comb <- data.table::setDT(as.data.frame(group_comb))
    is_counts <- temp %>%
        dplyr::mutate(count = purrr::map_int(is, ~ nrow(.x))) %>%
        dplyr::select(-.data$is)
    sharing_df <- if (!keep_genomic_coord) {
        group_comb[, shared := .find_in_common(.SD,
            lookup_tbl = temp,
            keep_genomic_coord = keep_genomic_coord
        ),
        .SDcols = cols
        ]
    } else {
        group_comb[, c("shared", "is_coord") := .find_in_common(.SD,
            lookup_tbl = temp,
            keep_genomic_coord = keep_genomic_coord
        ),
        .SDcols = cols
        ]
    }

    if (include_self_comp) {
        ## Calculate groups with equal components
        if (getOption("ISAnalytics.verbose") == TRUE) {
            rlang::inform("Calculating self groups (requested)...")
        }
        self_rows <- purrr::map2_df(
            is_counts$group_id, is_counts$count,
            function(x, y) {
                row_ls <- as.list(setNames(rep_len(
                    x,
                    length.out = n_comp
                ), nm = cols))
                row <- data.table::setDT(row_ls)
                row[, shared := y]
                if (keep_genomic_coord) {
                    row[, is_coord := list(temp[group_id == x]$is)]
                }
                return(row)
            }
        )
        sharing_df <- data.table::rbindlist(list(
            self_rows,
            sharing_df
        ))
    }
    group_cols <- colnames(sharing_df)[!colnames(sharing_df) %in% c(
        "shared",
        "is_coord"
    )]
    if (is_count || rel_sharing) {
        ## Add counts -groups
        for (col in group_cols) {
            count_col_name <- paste0("count_", col)
            sharing_df <- sharing_df[
                dplyr::rename(
                    is_counts,
                    !!col := .data$group_id,
                    !!count_col_name := .data$count
                ),
                on = col,
                nomatch = 0
            ]
        }
        ## Add counts -union
        sharing_df <- purrr::pmap_df(sharing_df, .count_group_union,
            col_groups = group_cols,
            lookup_tbl = temp
        )
    }
    if (!minimal) {
        if (getOption("ISAnalytics.verbose") == TRUE) {
            rlang::inform("Calculating permutations (requested)...")
        }
        sharing_df <- purrr::pmap_df(sharing_df, .sh_row_permut,
            g_names = group_cols,
            counts = is_count || rel_sharing
        )
    }
    if (rel_sharing) {
        # for groups
        for (col in group_cols) {
            rel_col_name <- paste0("on_", col)
            count_col_name <- paste0("count_", col)
            sharing_df <- sharing_df %>%
                dplyr::mutate(!!rel_col_name := (.data$shared /
                    .data[[count_col_name]]) * 100)
        }
        # for union
        sharing_df <- sharing_df %>%
            dplyr::mutate(on_union = (.data$shared /
                .data$count_union) * 100)
    }
    if (!is_count) {
        sharing_df <- sharing_df %>%
            dplyr::select(!dplyr::contains("count_"))
    }
    if (venn) {
        sharing_df <- sharing_df %>%
            dplyr::mutate(truth_tbl_venn = purrr::pmap(.,
                .sh_truth_tbl_venn,
                lookup = temp,
                groups = group_cols
            ))
    }
    sharing_df
}

## Computes sharing table for single input df and multiple keys
## -keys: name of the columns to do a group by (it is a named list)
## -minimal: true or false, if true computes ONLY combinations and not
## all permutations
## -is_count: keep the counts in the table?
## -rel_sharing: compute relative sharing? (on g_i & on_union)
.sharing_singledf_mult_key <- function(df, keys,
    minimal, is_count,
    rel_sharing, keep_genomic_coord,
    venn) {
    g_names <- names(keys)
    ## Obtain lookup table for each key
    lookup <- purrr::map(keys, ~ .sh_obtain_lookup(.x, df))
    group_labels <- purrr::map(lookup, ~ .x$group_id)
    unique_keys <- names(keys[!duplicated(keys)])
    lookup <- data.table::rbindlist(lookup[unique_keys])
    ## Obtain combinations
    combin <- data.table::setDT(purrr::cross_df(group_labels))
    sharing_df <- if (!keep_genomic_coord) {
        combin[, shared := .find_in_common(.SD,
            lookup_tbl = lookup,
            keep_genomic_coord = keep_genomic_coord
        ),
        .SDcols = g_names
        ]
    } else {
        combin[, c("shared", "is_coord") := .find_in_common(.SD,
            lookup_tbl = lookup,
            keep_genomic_coord = keep_genomic_coord
        ),
        .SDcols = g_names
        ]
    }

    if (is_count || rel_sharing) {
        is_counts <- lookup %>%
            dplyr::mutate(count = purrr::map_int(is, ~ nrow(.x))) %>%
            dplyr::select(-.data$is)
        ## Add counts -groups
        for (col in g_names) {
            count_col_name <- paste0("count_", col)
            sharing_df <- sharing_df[
                dplyr::rename(
                    is_counts,
                    !!col := .data$group_id,
                    !!count_col_name := .data$count
                ),
                on = col,
                nomatch = 0
            ]
        }
        ## Add counts -union
        sharing_df <- purrr::pmap_df(sharing_df, .count_group_union,
            col_groups = g_names,
            lookup_tbl = lookup
        )
    }
    if (!minimal) {
        if (getOption("ISAnalytics.verbose") == TRUE) {
            rlang::inform("Calculating permutations (requested)...")
        }
        sharing_df <- purrr::pmap_df(sharing_df, .sh_row_permut,
            g_names = g_names,
            counts = is_count || rel_sharing
        )
    }
    if (rel_sharing) {
        # for groups
        for (col in g_names) {
            rel_col_name <- paste0("on_", col)
            count_col_name <- paste0("count_", col)
            sharing_df <- sharing_df %>%
                dplyr::mutate(!!rel_col_name := (.data$shared /
                    .data[[count_col_name]]) * 100)
        }
        # for union
        sharing_df <- sharing_df %>%
            dplyr::mutate(on_union = (.data$shared /
                .data$count_union) * 100)
    }
    if (!is_count) {
        sharing_df <- sharing_df %>%
            dplyr::select(!dplyr::contains("count_"))
    }
    if (venn) {
        sharing_df <- sharing_df %>%
            dplyr::mutate(truth_tbl_venn = purrr::pmap(.,
                .sh_truth_tbl_venn,
                lookup = lookup,
                groups = g_names
            ))
    }
    sharing_df
}

## Computes sharing table for mult input dfs and single key
.sharing_multdf_single_key <- function(dfs, key, minimal,
    is_count, rel_sharing, keep_genomic_coord, venn) {
    if (is.null(names(dfs))) {
        g_names <- paste0("g", seq_len(length(dfs)))
        names(dfs) <- g_names
    }
    lookups <- purrr::map(dfs, ~ .sh_obtain_lookup(key, .x))
    group_labels <- purrr::map(lookups, ~ .x$group_id)
    lookup <- data.table::rbindlist(lookups)
    ## Obtain combinations
    combin <- data.table::setDT(purrr::cross_df(group_labels))
    sharing_df <- if (!keep_genomic_coord) {
        combin[, shared := .find_in_common(.SD,
            lookup_tbl = lookup,
            keep_genomic_coord = keep_genomic_coord
        ),
        .SDcols = names(dfs)
        ]
    } else {
        combin[, c("shared", "is_coord") := .find_in_common(.SD,
            lookup_tbl = lookup,
            keep_genomic_coord = keep_genomic_coord
        ),
        .SDcols = names(dfs)
        ]
    }

    if (is_count || rel_sharing) {
        is_counts <- purrr::map2(lookups, names(lookups), function(x, y) {
            count_col_name <- paste0("count_", y)
            x %>%
                dplyr::mutate(!!count_col_name := purrr::map_int(
                    is, ~ nrow(.x)
                )) %>%
                dplyr::select(-.data$is) %>%
                dplyr::rename(!!y := .data$group_id)
        })
        ## Add counts -groups
        for (col in names(dfs)) {
            cnt_tbl <- is_counts[[col]]
            sharing_df <- sharing_df[
                cnt_tbl,
                on = col,
                nomatch = 0
            ]
        }
        ## Add counts -union
        sharing_df <- purrr::pmap_df(sharing_df, .count_group_union,
            col_groups = names(dfs),
            lookup_tbl = lookup
        )
    }
    if (!minimal) {
        if (getOption("ISAnalytics.verbose") == TRUE) {
            rlang::inform("Calculating permutations (requested)...")
        }
        sharing_df <- purrr::pmap_df(sharing_df, .sh_row_permut,
            g_names = names(dfs),
            counts = is_count || rel_sharing
        )
    }
    if (rel_sharing) {
        # for groups
        for (col in names(dfs)) {
            rel_col_name <- paste0("on_", col)
            count_col_name <- paste0("count_", col)
            sharing_df <- sharing_df %>%
                dplyr::mutate(!!rel_col_name := (.data$shared /
                    .data[[count_col_name]]) * 100)
        }
        # for union
        sharing_df <- sharing_df %>%
            dplyr::mutate(on_union = (.data$shared /
                .data$count_union) * 100)
    }
    if (!is_count) {
        sharing_df <- sharing_df %>%
            dplyr::select(!dplyr::contains("count_"))
    }
    if (venn) {
        sharing_df <- sharing_df %>%
            dplyr::mutate(truth_tbl_venn = purrr::pmap(.,
                .sh_truth_tbl_venn,
                lookup = lookups,
                groups = names(dfs)
            ))
    }
    sharing_df
}

## Computes sharing table for mult input dfs and mult key
.sharing_multdf_mult_key <- function(dfs, keys, minimal,
    is_count, rel_sharing, keep_genomic_coord, venn) {
    lookups <- purrr::map2(dfs, keys, ~ .sh_obtain_lookup(.y, .x)) %>%
        purrr::set_names(names(keys))
    group_labels <- purrr::map(lookups, ~ .x$group_id)
    lookup <- data.table::rbindlist(lookups)
    ## Obtain combinations
    combin <- data.table::setDT(purrr::cross_df(group_labels))
    sharing_df <- if (!keep_genomic_coord) {
        combin[, shared := .find_in_common(.SD,
            lookup_tbl = lookup,
            keep_genomic_coord = keep_genomic_coord
        ),
        .SDcols = names(keys)
        ]
    } else {
        combin[, c("shared", "is_coord") := .find_in_common(.SD,
            lookup_tbl = lookup,
            keep_genomic_coord = keep_genomic_coord
        ),
        .SDcols = names(keys)
        ]
    }
    if (is_count || rel_sharing) {
        is_counts <- purrr::map2(lookups, names(lookups), function(x, y) {
            count_col_name <- paste0("count_", y)
            x %>%
                dplyr::mutate(!!count_col_name := purrr::map_int(
                    is, ~ nrow(.x)
                )) %>%
                dplyr::select(-.data$is) %>%
                dplyr::rename(!!y := .data$group_id)
        })
        ## Add counts -groups
        for (col in names(keys)) {
            cnt_tbl <- is_counts[[col]]
            sharing_df <- sharing_df[
                cnt_tbl,
                on = col,
                nomatch = 0
            ]
        }
        ## Add counts -union
        sharing_df <- purrr::pmap_df(sharing_df, .count_group_union,
            col_groups = names(keys),
            lookup_tbl = lookup
        )
    }
    if (!minimal) {
        if (getOption("ISAnalytics.verbose") == TRUE) {
            rlang::inform("Calculating permutations (requested)...")
        }
        sharing_df <- purrr::pmap_df(sharing_df, .sh_row_permut,
            g_names = names(keys),
            counts = is_count || rel_sharing
        )
    }
    if (rel_sharing) {
        # for groups
        for (col in names(keys)) {
            rel_col_name <- paste0("on_", col)
            count_col_name <- paste0("count_", col)
            sharing_df <- sharing_df %>%
                dplyr::mutate(!!rel_col_name := (.data$shared /
                    .data[[count_col_name]]) * 100)
        }
        # for union
        sharing_df <- sharing_df %>%
            dplyr::mutate(on_union = (.data$shared /
                .data$count_union) * 100)
    }
    if (!is_count) {
        sharing_df <- sharing_df %>%
            dplyr::select(!dplyr::contains("count_"))
    }
    if (venn) {
        sharing_df <- sharing_df %>%
            dplyr::mutate(truth_tbl_venn = purrr::pmap(.,
                .sh_truth_tbl_venn,
                lookup = lookups,
                groups = names(keys)
            ))
    }
    sharing_df
}

#---- USED IN : iss_source ----
.assign_iss_by_tp <- function(df, timepoint_column) {
    is_vars <- if (.is_annotated(df)) {
        c(mandatory_IS_vars(), annotation_IS_vars())
    } else {
        mandatory_IS_vars()
    }
    return(df %>%
        dplyr::mutate(!!timepoint_column := as.numeric(
            .data[[timepoint_column]]
        )) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(is_vars))) %>%
        dplyr::group_modify(~ {
            if (nrow(.x) == 1) {
                return(.x)
            }
            min_tp <- min(.x[[timepoint_column]])
            return(.x %>%
                dplyr::filter(.data[[timepoint_column]] == min_tp))
        }) %>%
        dplyr::ungroup())
}

.sharing_for_source <- function(ref,
    sel,
    ref_key,
    sel_key,
    tp_col,
    subj_col) {
    if (!is.data.frame(ref)) {
        ### Workflow by subject
        common_names <- if (!all(names(ref) == names(sel))) {
            intersect(names(ref), names(sel))
        } else {
            names(ref)
        }
        if (length(common_names) == 0) {
            no_common_err <- c("No common subjects",
                x = paste(
                    "No common subjects between",
                    "reference and selection"
                ),
                i = paste(
                    "In reference:",
                    paste0(names(ref), collapse = ", "),
                    "\n",
                    "In selection: ",
                    paste0(names(sel), collapse = ", ")
                )
            )
            rlang::abort(no_common_err)
        }
        if (getOption("ISAnalytics.verbose") == TRUE &&
            (length(common_names) < length(names(ref)) ||
                length(common_names) < length(names(sel)))) {
            all_names <- union(names(ref), names(sel))
            common_warn_msg <- c("Mismatch in subjects found",
                i = paste(
                    "Some subjects were excluded from",
                    "computations because they were",
                    "absent from reference or from",
                    "selection"
                ),
                paste(
                    "Excluded: ",
                    paste0(all_names[!all_names %in%
                        common_names],
                    collapse = ", "
                    )
                )
            )
        }
        if (.Platform$OS.type == "windows") {
            p <- BiocParallel::SnowParam(
                tasks = length(common_names),
                progressbar = getOption("ISAnalytics.verbose"),
                exportglobals = TRUE
            )
        } else {
            p <- BiocParallel::MulticoreParam(
                tasks = length(common_names),
                progressbar = getOption("ISAnalytics.verbose"),
                exportglobals = FALSE
            )
        }
        FUN <- function(subj,
    ref, sel, ref_key, sel_key,
    tp_col) {
            ref_df <- ref[[subj]]
            sel_df <- sel[[subj]]
            ref_key_min <- ref_key[ref_key != tp_col]
            ref_split <- ref_df %>%
                dplyr::group_by(dplyr::across(dplyr::all_of(ref_key_min))) %>%
                dplyr::group_split()
            quiet_sharing <- purrr::quietly(is_sharing)
            sharing <- purrr::map_df(ref_split, ~ {
                assigned_ref <- .assign_iss_by_tp(.x,
                    timepoint_column = tp_col
                )
                sh <- (quiet_sharing(assigned_ref,
                    sel_df,
                    group_keys = list(
                        g1 = ref_key,
                        g2 = sel_key
                    ),
                    keep_genomic_coord = TRUE
                ))$result
                sh <- sh %>%
                    tidyr::separate(
                        col = "g1",
                        into = paste0("g1_", ref_key),
                        remove = FALSE,
                        convert = TRUE
                    ) %>%
                    tidyr::separate(
                        col = "g2",
                        into = paste0("g2_", sel_key),
                        remove = FALSE,
                        convert = TRUE
                    )
            })
        }

        shared <- BiocParallel::bplapply(common_names, FUN,
            ref = ref,
            sel = sel,
            ref_key = ref_key,
            sel_key = sel_key,
            tp_col = tp_col,
            BPPARAM = p
        ) %>%
            purrr::set_names(common_names)
        return(shared)
    }
    ref_key_min <- ref_key[ref_key != tp_col]
    ref_split <- ref %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(ref_key_min))) %>%
        dplyr::group_split()
    quiet_sharing <- purrr::quietly(is_sharing)
    sharing <- purrr::map_df(ref_split, ~ {
        assigned_ref <- .assign_iss_by_tp(.x,
            timepoint_column = tp_col
        )
        (quiet_sharing(assigned_ref,
            sel,
            group_keys = list(
                g1 = ref_key,
                g2 = sel_key
            ),
            keep_genomic_coord = TRUE
        ))$result
    })
    sharing %>%
        tidyr::separate(
            col = "g1",
            into = paste0("g1_", ref_key),
            remove = FALSE,
            convert = TRUE
        ) %>%
        tidyr::separate(
            col = "g2",
            into = paste0("g2_", sel_key),
            remove = FALSE,
            convert = TRUE
        )
}
