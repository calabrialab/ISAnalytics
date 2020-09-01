#------------------------------------------------------------------------------#
# Internal functions
#------------------------------------------------------------------------------#
## All functions in this file are NOT exported, to be used internally only.

### Convenience functions for errors and warnings ###
#' @keywords internal
.malformed_ISmatrix_warning <- function() {
    paste(c(
        "Mandatory integration matrix variables, ", mandatory_IS_vars(),
        ", were not detected"
    ), collapse = " ")
}
#' @keywords internal
.non_ISM_error <- function() {
    paste(
        "One or more elements in x are not integration matrices.",
        "Aborting."
    )
}

#' @keywords internal
.missing_value_col_error <- function() {
    paste(
        "The `Value` column is missing or it contains non-numeric data.",
        "The column is needed for this operation.",
        "Aborting."
    )
}

#' @keywords internal
.missing_complAmpID_error <- function() {
    paste(
        "The `CompleteAmplificationID` column is missing.",
        "The column is needed for this operation.",
        "Aborting."
    )
}

#' @keywords internal
.quant_types_error <- function() {
    paste(
        "The list names must be quantification types",
        ", see quantification_types() for reference"
    )
}


#### ---- Internals for checks on integration matrices----####

#' Internal helper function for checking mandatory vars presence in x.
#'
#' Checks if the elements of `mandatory_IS_vars` are present as column names
#' in the data frame.
#' @param x A data.frame object (or any extending class)
#' @keywords internal
#'
#' @return FALSE if all or some elements are not found in the data frame, TRUE
#' otherwise
.check_mandatory_vars <- function(x) {
    stopifnot(is.data.frame(x))
    res <- if (all(mandatory_IS_vars() %in% colnames(x))) {
        TRUE
    } else {
        FALSE
    }
    return(res)
}

#' Internal helper function for checking `Value` column presence in x.
#'
#' Checks if the column `Value` is present in the data frame and also
#' checks if the column is numeric or integer.
#' @param x A data.frame object (or any extending class)
#' @keywords internal
#'
#' @return FALSE if not found or contains non-numeric data, TRUE otherwise
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

#' Internal helper function for checking `CompleteAmplifcationID`
#' column presence in x.
#'
#' Checks if the column `CompleteAmplifcationID` is present in the data frame.
#'
#' @param x A data.frame object (or any extending class)
#' @keywords internal
#'
#' @return FALSE if not found, TRUE otherwise
.check_complAmpID <- function(x) {
    stopifnot(is.data.frame(x))
    if ("CompleteAmplificationID" %in% colnames(x)) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

#### ---- Internals for matrix import ----####

#---- USED IN : import_single_Vispa2Matrix ----

#' Internal function to convert a messy matrix to a tidy data frame
#'
#' @description Uses the suite of functions provided by the
#' tidyverse to produce a more dense and ordered structure.
#' This function is not exported and should be called in other importing
#' functions.
#'
#' @param df Messy tibble to convert to tidy
#' @keywords internal
#'
#' @return a tidy tibble
#' @importFrom rlang .data
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr arrange all_of
#' @importFrom forcats fct_inseq as_factor
.messy_to_tidy <- function(df) {
    exp_cols <- which(!(colnames(df) %in% c(
        mandatory_IS_vars(),
        annotation_IS_vars()
    )))
    isadf_tidy <- df %>%
        tidyr::pivot_longer(
            cols = dplyr::all_of(exp_cols),
            names_to = "CompleteAmplificationID",
            values_to = "Value",
            values_drop_na = TRUE
        ) %>%
        dplyr::arrange(forcats::fct_inseq(forcats::as_factor(.data$chr)))
    isadf_tidy
}

#' Internal function to auto-detect the type of IS based on the headers.
#'
#' @param df the data frame to inspect
#' @keywords internal
#'
#' @return one value among:
#' * "OLD" : for old-style matrices that had only one column holding
#' all genomic coordinates
#' * "NEW" :  for the classic Vispa2 annotated/not annotated matrices
#' * "MALFORMED" : in any other case
.auto_detect_type <- function(df) {
    if ("IS_genomicID" %in% colnames(df) &
        all(!mandatory_IS_vars() %in% colnames(df))) {
        return("OLD")
    }
    if (all(mandatory_IS_vars() %in% colnames(df))) {
        return("NEW")
    }
    return("MALFORMED")
}

#---- USED IN : import_association_file ----

#' Checks if the association file has the right format (standard headers).
#'
#' @param df The imported association file
#' @keywords internal
#'
#' @return TRUE if the check passes, FALSE otherwise
.check_af_correctness <- function(df) {
    if (all(association_file_columns() %in% colnames(df))) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}


#' Reads association file and checks if it's correct or not.
#'
#' @param path The path to the association file on disk
#' @param padding The padding for TimePoint field
#' @param date_format The date format of date columns
#' @keywords internal
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate across contains
#' @importFrom rlang .data
#' @importFrom stringr str_pad
#' @import lubridate
#'
#' @return A tibble containing the association file.
.read_and_correctness_af <- function(path, padding, date_format) {
    stopifnot(is.character(path))
    stopifnot(file.exists(path))
    as_file <- read.csv(path,
        header = TRUE, check.names = FALSE,
        stringsAsFactors = FALSE, sep = "\t"
    )
    as_file <- tibble::as_tibble(as_file)
    # Checks if association file is correct
    correct <- .check_af_correctness(as_file)
    if (!correct) {
        stop(paste(
            "Malformed association file, could not import.",
            "Check the file or generate a new blank one with",
            "generate_blank_association_file()"
        ), call. = FALSE)
    }
    as_file <- as_file %>%
        dplyr::mutate(TimePoint = stringr::str_pad(
            as.character(.data$TimePoint),
            padding,
            side = "left",
            pad = "0"
        )) %>%
        dplyr::mutate(dplyr::across(dplyr::contains("Date"), ~ do.call(
            getExportedValue("lubridate", date_format),
            list(.x),
            quote = TRUE
        )))
    as_file
}

#' Internal function to check alignment between association file and file
#' system starting from the root. The alignment is checked at folder level.
#'
#' @param df The imported association file (data.frame or tibble)
#' @param root_folder Path to the root folder
#' @keywords internal
#' @importFrom dplyr select distinct mutate bind_rows
#' @importFrom fs dir_ls
#' @importFrom purrr pmap is_empty reduce map_dbl
#' @importFrom stringr str_replace_all str_extract_all
#' @importFrom tibble tibble
#'
#' @return A data frame containing, for each ProjectID and
#' concatenatePoolIDSeqRun the
#' corresponding path on disk if found, NA otherwise.
.check_file_system_alignment <- function(df, root_folder) {
    temp_df <- df %>%
        dplyr::select(
            .data$ProjectID,
            .data$concatenatePoolIDSeqRun,
            .data$PathToFolderProjectID
        ) %>%
        dplyr::distinct()
    tree_struct <- fs::dir_ls(
        path = root_folder, recurse = TRUE,
        type = "directory"
    )
    root_regexp <- stringr::str_replace_all(root_folder, "\\/", "\\\\/")
    results_df <- purrr::pmap(
        temp_df,
        function(ProjectID, concatenatePoolIDSeqRun,
    PathToFolderProjectID) {
            pattern <- paste0(
                "^", root_regexp,
                stringr::str_replace_all(
                    PathToFolderProjectID, "\\/", "\\\\/"
                ),
                "\\/quantification\\/",
                concatenatePoolIDSeqRun, "$"
            )
            pattern <- stringr::str_replace_all(pattern,
                pattern = "\\\\\\/\\\\\\/",
                replacement = "\\\\/"
            )
            found <- stringr::str_extract_all(tree_struct, pattern)
            found <- unlist(found)
            value <- if (purrr::is_empty(found)) {
                NA_character_
            } else {
                found
            }
            tibble::tibble(
                ProjectID = ProjectID,
                concatenatePoolIDSeqRun =
                    concatenatePoolIDSeqRun,
                Path = value
            )
        }
    )
    checker_df <- purrr::reduce(results_df, dplyr::bind_rows) %>%
        dplyr::mutate(
            Found = ifelse(!is.na(.data$Path), TRUE, FALSE),
            .before = .data$Path
        )
    checker_df
}

#' Updates the association file after the alignment check.
#'
#' The function checks if there are missing folders and updates the association
#' file by adding a column `Path` where the absolute path on disk for the
#' project and pool is found, if no path is found NA is inserted instead.
#'
#' @param as_file The tibble representing the read association_file
#' @param checks The tibble representing the results
#' of `.check_file_system_alignment`
#' @param root The root folder
#' @keywords internal
#'
#' @return An updated association file with absolute paths
.update_af_after_alignment <- function(as_file, checks, root) {
    # If some some folders are missing a warning is thrown
    not_found <- checks %>% dplyr::filter(.data$Found == FALSE)
    if (nrow(not_found) > 0) {
        warning(paste0("One or more projects were not found in the file",
            "system starting from ", root, ", please check your",
            "association file for errors and/or your file system.",
            "Until you re-import the association file these ",
            "missing files will be ignored.",
            collapse = ""
        ),
        immediate. = TRUE, call. = FALSE
        )
    }
    # Finally import modified association file
    checks <- checks %>%
        dplyr::select(
            .data$ProjectID, .data$concatenatePoolIDSeqRun,
            .data$Path
        )
    as_file <- as_file %>%
        dplyr::left_join(checks, by = c("ProjectID", "concatenatePoolIDSeqRun"))
    as_file
}

#---- USED IN : import_parallel_Vispa2Matrices_interactive ----

#' Helper function to be used when importing matrices in parallel.
#'
#' @param association_file Either the path to the association file or the tibble
#' representing the imported association file (done via
#' `import_association_file`)
#' @param root If `association_file` is the path to the file, root is a single
#' string holding the path to the root folder, otherwise root is `NULL`
#' @param padding The padding for TimePoint field
#' @param format The date format of date columns
#' @keywords internal
#'
#' @return A list of two elements: the first element is an updated version of
#' the association file with NAs removed, the second element is a widget showing
#' results of alignment checks.
.manage_association_file <- function(association_file, root, padding, format) {
    # Manage association file
    if (is.character(association_file)) {
        # If it's a path to file import the association file
        association_file <- .read_and_correctness_af(
            association_file,
            padding, format
        )
        checks <- .check_file_system_alignment(association_file, root)
        association_file <- .update_af_after_alignment(
            association_file,
            checks, root
        )
        association_file <- association_file %>% dplyr::filter(!is.na(
            .data$Path
        ))
        res <- list(association_file, checks)
        return(res)
    } else {
        # If it's a tibble (file already imported) check the correctness
        correct <- ifelse(.check_af_correctness(association_file),
            ifelse("Path" %in% colnames(association_file), TRUE,
                FALSE
            ), FALSE
        )
        if (isFALSE(correct)) {
            stop(paste("Malformed association file"))
        }
        association_file <- association_file %>% dplyr::filter(!is.na(
            .data$Path
        ))
        return(list(association_file, NULL))
    }
}

#' Allows the user to choose interactively the projects
#' to consider for import.
#'
#' @param association_file The tibble representing the imported association file
#' @keywords internal
#' @importFrom dplyr distinct select filter
#' @importFrom rlang .data
#' @importFrom stringr str_split
#' @importFrom purrr map_dbl
#'
#' @return A modified version of the association file where only selected
#' projects are present
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

#' Simple internal helper function to handle user input for selection
#' of number of pools.
#' @keywords internal
#'
#' @return Numeric representing user selection (1 for all pools, 2 for
#' only some pools, 0 to exit)
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

#' Simple helper interal function to handle user input for actual
#' pool choices.
#'
#' @param indexes A vector of integer indexes available
#' @keywords internal
#' @importFrom stringr str_split
#' @importFrom purrr map_dbl
#'
#' @return The user selection as a numeric vector
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

#' Allows the user to choose interactively
#' the pools to consider for import.
#'
#' @param association_file The tibble representing the imported
#' association file
#' @importFrom dplyr select distinct group_by bind_rows inner_join
#' @importFrom tibble tibble
#' @importFrom tidyr nest
#' @importFrom purrr map pmap reduce
#' @importFrom rlang .data
#' @keywords internal
#'
#' @return A modified version of the association file where only selected
#' pools for each project are present
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

#' Updates a files_found tibble with Files_count and Anomalies column.
#'
#' @param lups A files_found tibble obtained in a lookup function. Must contain
#' the Files column (nested table quantification type and files)
#' @keywords internal
#'
#' @return Updated files_found with Anomalies and Files_count columns
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

#' Looks up matrices to import given the association file and the
#' root of the file system.
#'
#' @param association_file Tibble representing the association file
#' @param quantification_type The type of quantification matrices to look for
#' (one in `quantification_types()`)
#' @param matrix_type The matrix_type to lookup (one between "annotated" or
#' "not_annotated")
#' @keywords internal
#' @importFrom tibble tibble
#' @importFrom fs dir_ls as_fs_path
#' @importFrom purrr map reduce map_dbl
#' @importFrom stringr str_detect
#' @importFrom dplyr select distinct bind_rows mutate
#' @importFrom tidyr nest
#'
#' @return A tibble containing all found files, including duplicates and missing
.lookup_matrices <- function(association_file,
    quantification_type,
    matrix_type) {
    temp <- association_file %>%
        dplyr::select(
            .data$ProjectID,
            .data$concatenatePoolIDSeqRun,
            .data$Path
        ) %>%
        dplyr::distinct()
    # Regex for matrix type matching
    matrix_type_regexp <- paste0("\\.no0\\.annotated")
    # Map, for each row in temp
    lups <- purrr::pmap(temp, function(...) {
        # Obtain a tibble with all the columns in temp but a single row
        current <- tibble::tibble(...)
        # Scan file system starting from Path
        files_found <- fs::dir_ls(path = current$Path)
        # Map for each quantification type requested
        matching <- purrr::map(quantification_type, function(x) {
            # Are there files with the requested quantification type in the
            # scanned folder?
            file_regexp <- paste0("\\_", x, "\\_", "matrix")
            detected <- unlist(stringr::str_detect(files_found,
                pattern = file_regexp
            ))
            found <- if (length(files_found[detected]) > 0) {
                # If yes subset and in the subset find only the ones that match
                # the matrix type regex
                files_found <- files_found[detected]
                detected <- unlist(stringr::str_detect(files_found,
                    pattern = matrix_type_regexp
                ))
                type_match <- if (matrix_type == "annotated") {
                    # If the wanted matrix type is annotated
                    if (length(detected) > 0) {
                        # If some paths were found to contain the type regex
                        # then return the subsetted vector of file paths
                        files_found[detected]
                    } else {
                        # If no paths were found to contain the type regex
                        # return NA
                        NA_character_
                    }
                } else {
                    # If the wanted matrix is NOT annotated
                    if (length(files_found[detected]) > 0) {
                        # If some paths were found to contain the type regex
                        if (length(files_found[!detected]) > 0) {
                            # If the files_found MINUS the detected ones is
                            # not an empty set then return the set
                            files_found[!detected]
                        } else {
                            # If the files_found MINUS the detected ones is an
                            # empty set then return NA
                            NA_character_
                        }
                    } else {
                        # If there were no paths that matched the type regex
                        # simply return the original files_found
                        files_found
                    }
                }
            } else {
                NA_character_
            }
            tibble::tibble(
                Quantification_type = x,
                Files_found = fs::as_fs_path(found)
            )
        })
        matching <- purrr::reduce(matching, dplyr::bind_rows) %>%
            dplyr::mutate(Files_found = fs::as_fs_path(.data$Files_found))

        tibble::tibble(
            ProjectID = current$ProjectID,
            concatenatePoolIDSeqRun = current$concatenatePoolIDSeqRun,
            matching
        )
    })

    lups <- purrr::reduce(lups, dplyr::bind_rows) %>%
        tidyr::nest(Files = c(.data$Quantification_type, .data$Files_found))
    lups <- .trace_anomalies(lups)
    lups
}


#' Simple function to manage user input for duplicate file choices for each
#' quantification type in a single project/pool pair (use internally in
#' `.manage_anomalies_interactive`).
#'
#' @param q_types Vector of characters containing the unique quantification
#' types that are detected as duplicates
#' @param dupl The tibble containing quantification types and path to the files
#' found for a single project/pool pair
#' @keywords internal
#' @importFrom dplyr filter slice bind_rows
#' @importFrom purrr map reduce
#' @importFrom rlang .data
#'
#' @return An updated tibble containing files chosen for each type
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

#' Manages anomalies for files found.
#'
#' The function manages anomalies found for files after scanning appropriate
#' folders by:
#' * Removing files not found (files for which Files_count$Found == 0 and Path
#' is NA) and printing a message to notify user
#' * Removing duplicates by asking the user which files to keep for each
#' quantification type, pool and project
#'
#' @param files_found The tibble obtained via calling `.lookup_matrices`
#' @keywords internal
#' @importFrom dplyr filter select rename bind_rows arrange
#' @importFrom tidyr unnest
#' @importFrom tibble as_tibble tibble
#' @importFrom purrr pmap flatten reduce
#'
#' @return A tibble containing for each project, pool and quantification type
#' the files chosen (ideally 1 for each quantification type if found, no more
#' than 1 per type)
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

#' Internal function for parallel import of a single quantification
#' type files.
#'
#' @param q_type The quantification type (single string)
#' @param files Files_found table were absolute paths of chosen files
#' are stored
#' @param workers Number of parallel workers
#' @keywords internal
#' @importFrom dplyr filter mutate bind_rows distinct
#' @importFrom BiocParallel SnowParam MulticoreParam bptry bplapply bpstop bpok
#' @importFrom purrr is_empty reduce
#'
#' @return A single tibble with all data from matrices of same quantification
#' type in tidy format
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
        matrix <- import_single_Vispa2Matrix(x)
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

#' Internal function for importing all files for each quantification type.
#'
#' @param files_to_import The tibble containing the files to import
#' @param workers Number of parallel workers
#' @keywords internal
#' @importFrom dplyr select distinct bind_rows
#' @importFrom purrr map set_names reduce flatten
#' @importFrom tibble as_tibble
#'
#' @return A named list of tibbles
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

#' Internal function to match user defined patterns on a vector
#' of file names.
#'
#' For each pattern specified by the user, the function tries to find a match
#' on all the file names and combines the results as a tibble in which column
#' names are the patterns and the values are TRUE if the pattern matched on the
#' element at that index or FALSE otherwise.
#'
#' @param filenames A character vector of file names
#' @param patterns A character vector of patterns to be matched
#' @keywords internal
#' @importFrom tibble as_tibble_col
#' @importFrom stringr str_detect
#' @importFrom purrr map reduce
#' @importFrom dplyr bind_cols
#'
#' @return A tibble
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

#' Helper function for checking if any of the elements of the list is true.
#'
#' @param ... A list of logical values
#' @keywords internal
#'
#' @return TRUE if any of the parameters is true
.any_match <- function(...) {
    l <- unlist(list(...))
    any(l)
}

#' Helper function for checking if all of the elements of the list is true.
#'
#' @param ... A list of logical values
#' @keywords internal
#'
#' @return TRUE if all of the parameters is true
.all_match <- function(...) {
    l <- unlist(list(...))
    all(l)
}

#' Updates files_found tibble according to pattern and matching options.
#'
#' @param files_nested The tibble containing Quantification_type and Files_found
#' columns relative to a single project/pool pair
#' @param p_matches The tibble representing the pattern matchings resulting from
#' `pattern_matching`
#' @param matching_opt The matching option
#' @keywords internal
#' @import dplyr
#' @importFrom purrr map reduce
#' @importFrom rlang .data
#' @importFrom fs as_fs_path
#'
#' @return An updated files_found tibble according to the matching option
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

#' Looks up matrices to import given the association file and the
#' root of the file system.
#'
#' @inheritParams .lookup_matrices
#' @param patterns A character vector of patterns to be matched
#' @param matching_opt A single character representing the matching option (one
#' of "ANY", "ALL" or "OPTIONAL")
#' @keywords internal
#' @importFrom purrr pmap flatten map
#' @importFrom tibble as_tibble
#' @importFrom stringr str_split
#' @importFrom dplyr mutate select
#' @importFrom utils tail
#'
#' @return A tibble containing all found files, including duplicates and missing
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

#' Manages anomalies for files found.
#'
#' The function manages anomalies found for files after scanning appropriate
#' folders by:
#' * Removing files not found (files for which Files_count$Found == 0 and Path
#' is NA) and printing a message to notify user
#' * Removing duplicates
#'
#' @param files_found The tibble obtained via calling `.lookup_matrices_auto`
#' @keywords internal
#' @importFrom dplyr filter select rename bind_rows arrange
#' @importFrom rlang .data
#' @importFrom tidyr unnest
#' @importFrom tibble as_tibble tibble
#' @importFrom purrr pmap flatten
#'
#' @return A tibble containing for each project, pool and quantification type
#' the files chosen
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

#### ---- Internals for HTML widgets construction ----####

#' Builds the html widget for the checker table.
#'
#' @param checker_df Tibble obtained via `.check_file_system_alignment`
#' @keywords internal
#' @importFrom reactable reactable reactableTheme colDef
#' @importFrom htmltools div span h2 css browsable
#'
#' @return An html widget
.checker_widget <- function(checker_df) {
    styled_df <- reactable::reactable(
        checker_df,
        striped = TRUE,
        sortable = TRUE,
        showSortable = TRUE,
        bordered = FALSE,
        outlined = TRUE,
        searchable = TRUE,
        pagination = TRUE,
        paginationType = "simple",
        showPageSizeOptions = TRUE,
        pageSizeOptions = c(4, 8, 12),
        defaultPageSize = 5,
        showPagination = TRUE,
        resizable = TRUE,
        defaultSorted = list(Found = "asc"),
        theme = reactable::reactableTheme(
            style = list(
                fontFamily = "Calibri"
            ),
            cellStyle = list(
                display = "flex",
                flexDirection = "column",
                justifyContent = "center"
            )
        ),
        defaultColDef = reactable::colDef(
            headerStyle = list(fontSize = "18px", paddingLeft = "15px"),
            align = "left"
        ),
        columns = list(
            ProjectID = reactable::colDef(
                align = "right",
                style = list(
                    color = "#9e9e9e",
                    fontWeight = "800",
                    borderRight = "2px solid #E6E6E6"
                ),
                minWidth = 60
            ),
            concatenatePoolIDSeqRun = reactable::colDef(
                minWidth = 100,
                style = list(paddingLeft = "15px")
            ),
            Found = reactable::colDef(
                maxWidth = 100,
                align = "center",
                style = function(value) {
                    color <- if (value == TRUE) {
                        "#6afc21"
                    } else {
                        "#d61e1e"
                    }
                    list(
                        color = color, paddingLeft = "15px",
                        fontWeight = "bold"
                    )
                },
                cell = function(value) {
                    if (value == TRUE) "\u2713" else "\u2718"
                }
            ),
            Path = reactable::colDef(
                minWidth = 200,
                style = list(paddingLeft = "15px")
            )
        )
    )
    widget <- htmltools::div(
        htmltools::div(
            htmltools::h2("ALIGNMENT RESULTS"),
            htmltools::span(paste(
                "Results of alignment between file system and",
                "association file. If some folders are not found",
                "they will be ignored until the problem is fixed",
                "and the association file re-imported."
            )),
            style = htmltools::css(font.family = "Calibri"),
            styled_df
        )
    )
    htmltools::browsable(widget)
}

#' Builds the html widget for the files_found table.
#'
#' @param files_found Tibble obtained via `.lookup_matrices` or
#' `.lookup_matrices_auto`
#' @keywords internal
#' @importFrom reactable reactable reactableTheme colDef
#' @importFrom htmltools div span h2 css h3 browsable
#' @importFrom dplyr select
#' @importFrom rlang .data
#'
#' @return An html widget
.files_found_widget <- function(files_found) {
    main_cols <- files_found %>% dplyr::select(
        .data$ProjectID,
        .data$concatenatePoolIDSeqRun,
        .data$Anomalies
    )
    styled_df <- reactable::reactable(
        main_cols,
        striped = TRUE,
        sortable = TRUE,
        showSortable = TRUE,
        bordered = FALSE,
        outlined = TRUE,
        searchable = TRUE,
        pagination = TRUE,
        paginationType = "simple",
        showPageSizeOptions = TRUE,
        pageSizeOptions = c(4, 8, 12),
        defaultPageSize = 5,
        showPagination = TRUE,
        resizable = TRUE,
        theme = reactable::reactableTheme(
            style = list(
                fontFamily = "Calibri"
            ),
            cellStyle = list(
                display = "flex",
                flexDirection = "column",
                justifyContent = "center"
            )
        ),
        defaultColDef = reactable::colDef(
            headerStyle = list(fontSize = "18px", paddingLeft = "15px"),
            align = "left"
        ),
        columns = list(
            ProjectID = reactable::colDef(
                align = "right",
                style = list(
                    color = "#9e9e9e",
                    fontWeight = "800",
                    borderRight = "2px solid #E6E6E6"
                ),
                minWidth = 60
            ),
            concatenatePoolIDSeqRun = reactable::colDef(
                minWidth = 100,
                style = list(paddingLeft = "15px")
            ),
            Anomalies = reactable::colDef(
                minWidth = 150,
                align = "center",
                style = function(value) {
                    color <- if (value == TRUE) {
                        "#f2cd29"
                    } else {
                        "#6afc21"
                    }
                    list(
                        color = color, paddingLeft = "15px",
                        fontWeight = "bold",
                        fontSize = "20px"
                    )
                },
                cell = function(value) {
                    if (value == TRUE) "\u26A0" else "\u2713"
                }
            )
        ),
        details = function(index) {
            files <- files_found$Files[[index]]
            count <- files_found$Files_count[[index]]
            styled_files <- reactable::reactable(
                files,
                bordered = FALSE,
                outlined = TRUE,
                resizable = TRUE,
                striped = TRUE,
                pagination = TRUE,
                defaultPageSize = 4,
                showPagination = TRUE,
                paginationType = "simple",
                defaultColDef = reactable::colDef(
                    headerStyle = list(fontSize = "18px", paddingLeft = "15px"),
                    align = "left",
                    style = list(paddingLeft = "15px"),
                    header = function(value) gsub("_", " ", value, fixed = TRUE)
                ),
                columns = list(
                    Quantification_type = reactable::colDef(
                        minWidth = 200,
                        maxWidth = 200,
                        filterable = TRUE
                    )
                )
            )
            styled_count <- reactable::reactable(
                count,
                bordered = FALSE,
                outlined = TRUE,
                resizable = TRUE,
                striped = TRUE,
                defaultColDef = reactable::colDef(
                    headerStyle = list(fontSize = "18px", paddingLeft = "15px"),
                    align = "left",
                    style = list(paddingLeft = "15px")
                ),
                columns = list(
                    Quantification_type = reactable::colDef(
                        header = function(value) {
                            gsub("_", " ", value,
                                fixed = TRUE
                            )
                        }
                    ),
                    Found = reactable::colDef(
                        cell = function(value) {
                            if (value == 1) {
                                value
                            } else {
                                if (value > 1) {
                                    paste(value, "\u2691")
                                } else {
                                    paste(value, "\u2716")
                                }
                            }
                        },
                        style = function(value) {
                            if (value == 1) {
                                color <- "black"
                                weight <- "normal"
                            } else {
                                if (value > 1) {
                                    color <- "#f2cd29"
                                    weight <- "bold"
                                } else {
                                    color <- "#d61e1e"
                                    weight <- "bold"
                                }
                            }
                            list(
                                color = color, fontWeight = weight,
                                paddingLeft = "15px"
                            )
                        }
                    )
                )
            )
            htmltools::div(
                style = paste(
                    "padding-left: 40px;",
                    "padding-right: 40px;",
                    "padding-bottom: 20px"
                ),
                htmltools::h3(paste(
                    "Summary of files count for each",
                    "quantification type"
                )),
                styled_count,
                htmltools::h3(paste(
                    "Summary of files found for each",
                    "quantification type"
                )),
                styled_files
            )
        }
    )
    widget <- htmltools::div(
        htmltools::h2("INTEGRATION MATRICES FOUND REPORT"),
        htmltools::span(paste(
            "Report of all files found for each quantification",
            "type. Click on the arrow on the left side of each",
            "row to see details."
        )),
        style = htmltools::css(font.family = "Calibri"),
        styled_df
    )
    htmltools::browsable(widget)
}

#' Builds the html widget for the files_to_import table.
#'
#' @param files_to_import Tibble obtained via
#' `.manage_anomalies_interactive` or
#' `.manage_anomalies_auto`
#' @keywords internal
#' @importFrom reactable reactable reactableTheme colDef
#' @importFrom htmltools div span h2 css browsable
#'
#' @return An html widget
.files_to_import_widget <- function(files_to_import) {
    styled_df <- reactable::reactable(
        files_to_import,
        striped = TRUE,
        sortable = TRUE,
        showSortable = TRUE,
        bordered = FALSE,
        outlined = TRUE,
        searchable = TRUE,
        pagination = TRUE,
        paginationType = "simple",
        showPageSizeOptions = TRUE,
        pageSizeOptions = c(4, 8, 12),
        defaultPageSize = 5,
        showPagination = TRUE,
        resizable = TRUE,
        theme = reactable::reactableTheme(
            style = list(
                fontFamily = "Calibri"
            ),
            cellStyle = list(
                display = "flex",
                flexDirection = "column",
                justifyContent = "center"
            )
        ),
        defaultColDef = reactable::colDef(
            headerStyle = list(
                fontSize = "18px", paddingLeft = "15px",
                display = "flex",
                flexDirection = "column",
                justifyContent = "center"
            ),
            style = list(paddingLeft = "15px"),
            align = "left",
            header = function(value) gsub("_", " ", value, fixed = TRUE)
        ),
        columns = list(
            Files_chosen = reactable::colDef(
                minWidth = 250
            ),
            ProjectID = reactable::colDef(
                filterable = TRUE
            ),
            concatenatePoolIDSeqRun = reactable::colDef(
                filterable = TRUE
            ),
            Quantification_type = reactable::colDef(
                filterable = TRUE,
                align = "center"
            )
        )
    )
    widget <- htmltools::div(
        style = htmltools::css(font.family = "Calibri"),
        htmltools::h2("SUMMARY OF FILES CHOSEN FOR IMPORT"),
        htmltools::span("Here is a summary of all files chosen for import"),
        styled_df
    )
    htmltools::browsable(widget)
}

#' Builds the html widget for the files_imported table.
#'
#' @param files_imported Tibble obtained via `.parallel_import_merge`
#' @keywords internal
#' @importFrom reactable reactable reactableTheme colDef
#' @importFrom htmltools div span h2 css browsable
#'
#' @return An html widget
.files_imported_widget <- function(files_imported) {
    styled_df <- reactable::reactable(
        files_imported,
        striped = TRUE,
        sortable = TRUE,
        showSortable = TRUE,
        bordered = FALSE,
        outlined = TRUE,
        searchable = TRUE,
        pagination = TRUE,
        paginationType = "simple",
        showPageSizeOptions = TRUE,
        pageSizeOptions = c(4, 8, 12),
        defaultPageSize = 5,
        showPagination = TRUE,
        resizable = TRUE,
        theme = reactable::reactableTheme(
            style = list(
                fontFamily = "Calibri"
            ),
            cellStyle = list(
                display = "flex",
                flexDirection = "column",
                justifyContent = "center"
            )
        ),
        defaultColDef = reactable::colDef(
            headerStyle = list(
                fontSize = "18px", paddingLeft = "15px",
                display = "flex",
                flexDirection = "column",
                justifyContent = "center"
            ),
            style = list(paddingLeft = "15px"),
            align = "left",
            header = function(value) gsub("_", " ", value, fixed = TRUE)
        ),
        columns = list(
            Files_chosen = reactable::colDef(
                minWidth = 250
            ),
            ProjectID = reactable::colDef(
                filterable = TRUE
            ),
            concatenatePoolIDSeqRun = reactable::colDef(
                filterable = TRUE
            ),
            Quantification_type = reactable::colDef(
                filterable = TRUE,
                align = "center"
            ),
            Imported = reactable::colDef(
                style = function(value) {
                    color <- if (value == TRUE) {
                        "#6afc21"
                    } else {
                        "#d61e1e"
                    }
                    list(
                        paddingLeft = "15px",
                        textTransform = "uppercase",
                        color = color,
                        fontWeight = "bold"
                    )
                },
                align = "center"
            )
        )
    )
    widget <- htmltools::div(
        style = htmltools::css(font.family = "Calibri"),
        htmltools::h2("REPORT: FILES IMPORTED"),
        htmltools::span("Here is a summary of all files actually imported for
        each quantification type. If you see 'false' in the column Imported,
        some errors might have occurred and the function was unable to import
                        that matrix."),
        styled_df
    )
    htmltools::browsable(widget)
}

#' Builds the html widget for the summary table.
#'
#' @param removed Number of removed collisions
#' @param reassigned Number of re-assigned collisions
#' @param summary Summary table
#' @param tot_rows Total number rows of sequence count matrix before processing
#' @param collision_rows Total number of rows of collisions
#' @keywords internal
#' @importFrom reactable reactable reactableTheme colDef
#' @importFrom htmltools div h2 h4 css browsable
#'
#' @return A widget
.summary_collisions_widget <- function(removed,
    reassigned,
    summary,
    tot_rows,
    collision_rows) {
    theme <- reactable::reactableTheme(
        style = list(
            fontFamily = "Calibri"
        ),
        cellStyle = list(
            display = "flex",
            flexDirection = "column",
            justifyContent = "center"
        )
    )

    styled_df <- reactable::reactable(
        summary,
        striped = TRUE,
        sortable = TRUE,
        showSortable = TRUE,
        bordered = FALSE,
        outlined = TRUE,
        searchable = TRUE,
        pagination = TRUE,
        paginationType = "simple",
        showPageSizeOptions = TRUE,
        pageSizeOptions = c(4, 8, 12),
        defaultPageSize = 5,
        showPagination = TRUE,
        resizable = TRUE,
        theme = theme,
        defaultColDef = reactable::colDef(
            headerStyle = list(
                fontSize = "18px", paddingLeft = "15px",
                display = "flex",
                flexDirection = "column",
                justifyContent = "center"
            ),
            style = list(paddingLeft = "15px"),
            align = "left",
            header = function(value) gsub("_", " ", value, fixed = TRUE),
            filterable = TRUE
        )
    )

    collisions_found <- tibble::tibble(
        abs_number = collision_rows,
        percentage_on_total =
            collision_rows / tot_rows
    )
    collisions_removed <- tibble::tibble(
        abs_number = removed,
        percentage_on_collisions =
            removed / collision_rows,
        percentage_on_total =
            removed / tot_rows
    )
    collisions_reassigned <- tibble::tibble(
        abs_number = reassigned,
        percentage_on_collisions =
            reassigned / collision_rows,
        percentage_on_total =
            reassigned / tot_rows
    )

    collisions_found_styled <- reactable::reactable(
        collisions_found,
        fullWidth = FALSE,
        bordered = FALSE,
        outlined = TRUE,
        theme = theme,
        defaultColDef = reactable::colDef(
            header = function(value) gsub("_", " ", value, fixed = TRUE)
        ),
        columns = list(
            percentage_on_total = reactable::colDef(
                format = reactable::colFormat(
                    percent = TRUE,
                    digits = 2
                )
            )
        )
    )

    collisions_removed_styled <- reactable::reactable(
        collisions_removed,
        fullWidth = FALSE,
        bordered = FALSE,
        outlined = TRUE,
        theme = theme,
        defaultColDef = reactable::colDef(
            header = function(value) gsub("_", " ", value, fixed = TRUE)
        ),
        columns = list(
            percentage_on_total = reactable::colDef(
                format = reactable::colFormat(
                    percent = TRUE,
                    digits = 2
                )
            ),
            percentage_on_collisions = reactable::colDef(
                format = reactable::colFormat(
                    percent = TRUE,
                    digits = 2
                )
            )
        )
    )

    collisions_reassigned_styled <- reactable::reactable(
        collisions_reassigned,
        fullWidth = FALSE,
        bordered = FALSE,
        outlined = TRUE,
        theme = theme,
        defaultColDef = reactable::colDef(
            header = function(value) gsub("_", " ", value, fixed = TRUE)
        ),
        columns = list(
            percentage_on_total = reactable::colDef(
                format = reactable::colFormat(
                    percent = TRUE,
                    digits = 2
                )
            ),
            percentage_on_collisions = reactable::colDef(
                format = reactable::colFormat(
                    percent = TRUE,
                    digits = 2
                )
            )
        )
    )

    widget <- htmltools::div(
        style = htmltools::css(font.family = "Calibri"),
        htmltools::h2("COLLISION REMOVAL SUMMARY"),
        htmltools::h4("TOTAL READS:"),
        htmltools::div(tot_rows),
        htmltools::h4("COLLISIONS FOUND:"),
        collisions_found_styled,
        htmltools::h4("REMOVED:"),
        collisions_removed_styled,
        htmltools::h4("REASSIGNED:"),
        collisions_reassigned_styled,
        htmltools::h4("SUMMARY:"),
        styled_df
    )
    htmltools::browsable(widget)
}

#' Builds the html widget for the iss_import.
#'
#' @param report Table obtained via `import_stats_iss`
#' @keywords internal
#'
#' @importFrom reactable reactable reactableTheme colDef
#' @importFrom htmltools div h2 h4 css browsable
#'
#' @return A widget
.iss_import_widget <- function(report) {
    theme <- reactable::reactableTheme(
        style = list(
            fontFamily = "Calibri"
        ),
        cellStyle = list(
            display = "flex",
            flexDirection = "column",
            justifyContent = "center"
        )
    )

    styled_df <- reactable::reactable(
        report,
        striped = TRUE,
        sortable = TRUE,
        showSortable = TRUE,
        bordered = FALSE,
        outlined = TRUE,
        searchable = TRUE,
        pagination = TRUE,
        paginationType = "simple",
        showPageSizeOptions = TRUE,
        pageSizeOptions = c(4, 8, 12),
        defaultPageSize = 5,
        showPagination = TRUE,
        resizable = TRUE,
        theme = theme,
        defaultColDef = reactable::colDef(
            headerStyle = list(
                fontSize = "18px", paddingLeft = "15px",
                display = "flex",
                flexDirection = "column",
                justifyContent = "center"
            ),
            style = list(paddingLeft = "15px"),
            align = "left",
            header = function(value) gsub("_", " ", value, fixed = TRUE)
        ),
        columns = list(
            ProjectID = reactable::colDef(
                filterable = TRUE
            ),
            Imported = reactable::colDef(
                style = function(value) {
                    color <- if (value == TRUE) {
                        "#6afc21"
                    } else {
                        "#d61e1e"
                    }
                    list(
                        paddingLeft = "15px",
                        textTransform = "uppercase",
                        color = color,
                        fontWeight = "bold"
                    )
                },
                align = "center"
            )
        )
    )
    widget <- htmltools::div(
        style = htmltools::css(font.family = "Calibri"),
        htmltools::h2("REPORT IMPORT VISPA2 STATS: FILES IMPORTED"),
        htmltools::span("Here is a summary of all files actually imported.
        If you see 'FALSE' in the column Imported, some errors might have
        occurred and the function was unable to import the file or simply no
                        path was found for that stats file."),
        styled_df
    )
    htmltools::browsable(widget)
}

#### ---- Internals for collision removal ----####

#---- USED IN : remove_collisions ----

#' Checks if association file contains more information than the matrix.
#'
#' Used to notify the user that wants to know if for the projects contained in
#' the examined matrix there are additional CompleteAmplificationIDs contained
#' in the association file that weren't included in the integration matrix (for
#' example failed runs).
#'
#' @param association_file The imported association file
#' @param df The sequence count matrix to examine
#' @keywords internal
#' @import dplyr
#' @importFrom purrr map reduce
#' @importFrom rlang .data
#'
#' @return A tibble containing ProjectID, PoolID and CompleteAmplificationID
#' only for additional info and only for projects already present in df
.check_same_info <- function(association_file, df) {
    joined <- association_file %>%
        dplyr::left_join(df, by = "CompleteAmplificationID")
    projects <- (joined %>% dplyr::select(.data$ProjectID) %>%
        dplyr::distinct())$ProjectID
    joined1 <- purrr::map(projects, function(x) {
        temp <- joined %>% dplyr::filter(.data$ProjectID == x)
        if (!all(is.na(temp$chr))) {
            temp %>%
                dplyr::filter(is.na(.data$chr)) %>%
                dplyr::select(
                    .data$ProjectID,
                    .data$PoolID,
                    .data$CompleteAmplificationID
                )
        }
    })
    purrr::reduce(joined1, dplyr::bind_rows)
}

#' Produces a joined tibble between the sequence count matrix and the
#' association file
#'
#' @param seq_count_df The sequence count tibble
#' @param association_file The association file tibble
#' @param date_col The date column chosen
#' @keywords internal
#' @import dplyr
#' @importFrom rlang .data
#'
#' @return A tibble
.join_matrix_af <- function(seq_count_df, association_file, date_col) {
    joined <- seq_count_df %>%
        dplyr::left_join(association_file, by = "CompleteAmplificationID") %>%
        dplyr::select(
            dplyr::all_of(colnames(seq_count_df)), all_of(date_col),
            .data$ProjectID, .data$PoolID, .data$SubjectID,
            .data$ReplicateNumber
        )
    joined
}

#' Identifies independent samples and separates the joined_df in
#' collisions and non-collisions
#'
#' @param joined_df The joined tibble obtained via `.join_matrix_af`
#' @import dplyr
#' @importFrom rlang .data
#' @keywords internal
#'
#' @return A named list containing the splitted joined_df for collisions and
#' non-collisions
.identify_independent_samples <- function(joined_df) {
    temp <- joined_df %>%
        dplyr::select(
            dplyr::all_of(mandatory_IS_vars()),
            .data$ProjectID, .data$SubjectID
        ) %>%
        dplyr::group_by(dplyr::across(mandatory_IS_vars())) %>%
        dplyr::distinct(.data$ProjectID, .data$SubjectID, .keep_all = TRUE) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
        dplyr::ungroup() %>%
        dplyr::filter(.data$n > 1) %>%
        dplyr::select(-c(.data$n))

    non_collisions <- joined_df %>%
        dplyr::anti_join(temp, by = mandatory_IS_vars())
    collisions <- joined_df %>%
        dplyr::right_join(temp, by = mandatory_IS_vars())
    list(collisions = collisions, non_collisions = non_collisions)
}

#' Returns a polished version of collisions table for analysis.
#'
#' Polishing consists in nesting all information relative to a single
#' integration, in other words, obtaining a linked data frame for each
#' integration.
#'
#' @param collisions The collisions table obtained via
#' `.identify_independent_samples`
#' @keywords internal
#' @importFrom dplyr group_by across
#' @importFrom rlang .data
#' @importFrom tidyr nest
#'
#' @return A nested tibble
.obtain_nested <- function(collisions) {
    collisions %>%
        dplyr::group_by(dplyr::across(mandatory_IS_vars())) %>%
        tidyr::nest()
}


#' Internal for date discrimination in collision removal.
#'
#' It's the first of 4 steps in the algorithm for collision removal: it tries to
#' find a single sample who has an associated date which is earlier than any
#' other. If comparison is not possible the analysis fails and returns the
#' orginal input data.
#'
#' @details NOTE: this function is meant to be used inside a mapping function
#' such as `purrr::pmap`. The function only works on data regarding a SINGLE
#' integration (triplet chr, integration_locus, strand).
#'
#' @param nest The nested table associated with a single integration
#' @param date_col The name of the date column chosen for collision removal,
#' one in \code{date_columns_coll()}
#' @importFrom rlang expr eval_tidy .data
#' @importFrom dplyr filter arrange
#' @keywords internal
#'
#' @return A named list with:
#' * $data: a tibble, containing the data (unmodified or modified)
#' * $check: a logical value indicating whether the analysis was successful or
#' not (and therefore there is the need to perform the next step)
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

#' Internal for replicate discrimination in collision removal.
#'
#' It's the second of 4 steps in the algorithm for collision removal:
#' grouping by independent sample (ProjectID, SubjectID) it counts the number of
#' rows (replicates) found for each group, orders them from biggest to smallest
#' and, if a single group has more rows than any other group the integration is
#' assigned to that sample, otherwise the analysis fails and returns the
#' original input to be submitted to the next step.
#'
#' @details NOTE: this function is meant to be used inside a mapping function
#' such as `purrr::pmap`. The function only works on data regarding a SINGLE
#' integration (triplet chr, integration_locus, strand).
#'
#' @param nest The nested table associated with a single integration
#' @import dplyr
#' @importFrom rlang .data
#' @keywords internal
#'
#' @return A named list with:
#' * $data: a tibble, containing the data (unmodified or modified)
#' * $check: a logical value indicating whether the analysis was successful or
#' not (and therefore there is the need to perform the next step)
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

#' Internal for sequence count discrimination in collision removal.
#'
#' It's the third of 4 steps in the algorithm for collision removal:
#' grouping by independent sample (ProjectID, SubjectID), it sums the value of
#' the sequence count for each group, sorts the groups from highest to lowest
#' cumulative value and then checks the ratio between the first element and the
#' second: if the ratio is > `reads_ratio` the integration is assigned to
#' the first group, otherwise the analysis fails and returns the original input.
#'
#' @details NOTE: this function is meant to be used inside a mapping function
#' such as `purrr::pmap`. The function only works on data regarding a SINGLE
#' integration (triplet chr, integration_locus, strand).
#'
#' @param nest The nested table associated with a single integration
#' @param reads_ratio The value of the ratio between sequence count values to
#' check
#' @import dplyr
#' @importFrom rlang .data
#' @keywords internal
#'
#' @return A named list with:
#' * $data: a tibble, containing the data (unmodified or modified)
#' * $check: a logical value indicating whether the analysis was successful or
#' not (and therefore there is the need to perform the next step)
.discriminate_by_seqCount <- function(nest, reads_ratio) {
    temp <- nest %>%
        dplyr::group_by(.data$ProjectID, .data$SubjectID) %>%
        dplyr::summarise(sum_seqCount = sum(.data$Value), .groups = "keep") %>%
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

#' Internal function that performs four-step-check of collisions for
#' a single integration.
#'
#' @details NOTE: this function is meant to be used inside a mapping function
#' such as `purrr::pmap`. The function only works on data regarding a SINGLE
#' integration (triplet chr, integration_locus, strand).
#'
#' @param ... Represents a single row of a tibble obtained via
#' `obtain_nested`. One row contains 4 variables: chr, integration_locus,
#' strand and data, where data is a nested table that contains all rows
#' that share that integration (collisions).
#' @param date_col The date column to consider
#' @param reads_ratio The value of the ratio between sequence count values to
#' check
#' @keywords internal
#'
#' @importFrom tibble as_tibble tibble
#' @importFrom purrr flatten
#' @importFrom tidyr unnest
#' @importFrom rlang .data env_bind
#' @return A list with:
#' * $data: an updated tibble with processed collisions or NULL if no
#' criteria was sufficient
#' * $reassigned: 1 if the integration was successfully reassigned, 0 otherwise
#' * $removed: 1 if the integration is removed entirely because no criteria was
#' met, 0 otherwise
.four_step_check <- function(..., date_col, reads_ratio) {
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
    result <- .discriminate_by_seqCount(current_data, reads_ratio)
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

#' Internal function to process collisions on multiple integrations.
#'
#' @param x The nested collision data frame (or chunks for parallel execution)
#' @param date_col The date column to consider
#' @param reads_ratio The value of the ratio between sequence count values to
#' check
#' @keywords internal
#' @importFrom purrr pmap reduce
#' @importFrom dplyr bind_rows
#' @return A list containing the updated collisions, a numeric value
#' representing the number of integrations removed and a numeric value
#' representing the number of integrations reassigned
.coll_mapping <- function(x, date_col, reads_ratio) {
    result <- purrr::pmap(x,
        .f = .four_step_check,
        date_col = date_col,
        reads_ratio = reads_ratio
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

#' Internal function to process collisions on multiple integrations,
#' parallelized.
#'
#' @param collisions The collision data frame
#' @param date_col The date column to consider
#' @param reads_ratio The value of the ratio between sequence count values to
#' check
#' @keywords internal
#'
#' @import BiocParallel
#' @importFrom purrr map reduce
#' @importFrom tibble as_tibble_col add_column
#' @importFrom dplyr group_by group_split
#' @return A list containing the updated collisions, a numeric value
#' representing the number of integrations removed and a numeric value
#' representing the number of integrations reassigned
.process_collisions <- function(collisions, date_col, reads_ratio) {
    # Obtain nested version of collisions
    nested <- .obtain_nested(collisions)
    # Register backend according to platform
    if (.Platform$OS.type == "windows") {
        p <- BiocParallel::SnowParam(stop.on.error = FALSE)
    } else {
        p <- BiocParallel::MulticoreParam(stop.on.error = FALSE)
    }
    # Split the data frame in chunks according to number of workers
    workers <- BiocParallel::bpworkers(p)
    exceeding <- nrow(nested) %% workers
    ampl <- trunc(nrow(nested) / workers)
    chunks <- unlist(lapply(seq_len(workers), FUN = rep_len, length.out = ampl))
    if (exceeding > 0) {
        chunks <- c(chunks, rep_len(
            x = tail(workers, n = 1),
            length.out = exceeding
        ))
    }
    chunks <- tibble::as_tibble_col(chunks, column_name = "chunk")
    nested <- nested %>% tibble::add_column(chunks)
    split_data <- nested %>%
        dplyr::group_by(.data$chunk) %>%
        dplyr::group_split(.keep = FALSE)
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

#' Internal function to obtain a summary table after collision processing.
#'
#' @param before The joined table before collision removing
#' (obtained via `.join_matrix_af`)
#' @param after The final matrix obtained after collision processing
#' @param association_file The association file
#' @import dplyr
#' @importFrom rlang .data
#' @keywords internal
#'
#' @return A tibble with a summary containing for each SubjectID the number of
#' integrations found before and after, the sum of the value of the sequence
#' count for each subject before and after and the corresponding deltas.
.summary_table <- function(before, after, association_file) {
    after_matr <- after %>%
        dplyr::left_join(association_file,
            by = c("CompleteAmplificationID")
        ) %>%
        dplyr::select(dplyr::all_of(colnames(after)), .data$SubjectID)

    n_int_after <- after_matr %>%
        dplyr::group_by(.data$SubjectID) %>%
        dplyr::tally(name = "count_integrations_after")
    sumReads_after <- after_matr %>%
        dplyr::group_by(.data$SubjectID) %>%
        dplyr::summarise(
            sumSeqReads_after = sum(.data$Value),
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
            sumSeqReads_before = sum(.data$Value),
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

#### ---- Internals for aggregate functions ----####

#---- USED IN : aggregate_metadata ----

#' Minimal association_file variable set.
#'
#' Contains the names of the columns of the association file that are a minimum
#' requirement to perform aggregation.
#' @keywords internal
#'
#' @return A character vector
.min_var_set <- function() {
    c(
        "FusionPrimerPCRDate", "LinearPCRDate", "VCN", "DNAngUsed", "Kapa",
        "ulForPool", "Path"
    )
}

#' Minimal stats column set.
#'
#' Contains the name of the columns that are a minimum requirement for
#' aggregation.
#' @keywords internal
#'
#' @return A character vector
.stats_columns_min <- function() {
    c(
        "POOL", "TAG", "BARCODE_MUX", "TRIMMING_FINAL_LTRLC",
        "LV_MAPPED", "BWA_MAPPED_OVERALL", "ISS_MAPPED_PP"
    )
}

#' Checks if the stats file contains the minimal set of columns.
#'
#' @param x The stats df
#' @keywords internal
#'
#' @return TRUE or FALSE
.check_stats <- function(x) {
    if (all(.stats_columns_min() %in% colnames(x))) {
        TRUE
    } else {
        FALSE
    }
}

#' Finds automatically the path on disk to each stats file.
#'
#' @param association_file The association file
#' @import dplyr
#' @importFrom purrr pmap_dfr pmap_df
#' @importFrom tibble as_tibble tibble add_column
#' @importFrom stringr str_extract str_extract_all str_detect
#' @importFrom fs dir_exists dir_ls
#' @importFrom tidyr unnest unite
#' @importFrom rlang .data
#' @keywords internal
#' @return A tibble with ProjectID and the absolute path on disk to each file
#' if found
.stats_report <- function(association_file) {
    # Obtain unique projectID and path
    temp <- association_file %>%
        dplyr::select(.data$ProjectID, .data$Path) %>%
        dplyr::distinct()
    # If paths are all NA return
    if (all(is.na(temp$Path))) {
        return(NULL)
    }
    pattern <- "^stats\\."
    file_reg <- "(?<=iss\\/).*"
    # Obtain paths to iss folder
    stats_paths <- purrr::pmap_dfr(temp, function(...) {
        current <- tibble::tibble(...)
        if (is.na(current$Path)) {
            l <- list(ProjectID = current$ProjectID, stats_path = NA_character_)
            tibble::as_tibble(l)
        } else {
            pj_path <-
                stringr::str_extract(current$Path,
                    pattern = paste0(
                        ".*",
                        current$ProjectID,
                        "(?=\\/quantification)"
                    )
                )
            pj_path <- paste0(pj_path, "/iss")
            l <- list(ProjectID = current$ProjectID, stats_path = pj_path)
            tibble::as_tibble(l)
        }
    })
    stats_paths <- stats_paths %>% dplyr::distinct()
    # Find all stats files in iss folders
    stats_paths <- purrr::pmap_df(stats_paths, function(...) {
        cur <- tibble::tibble(...)
        # Check if folder exists
        # Set to NA the iss folders not found
        if (!fs::dir_exists(cur$stats_path)) {
            cur$stats_path <- NA_character_
            cur %>% tibble::add_column(stats_files = list(NA_character_))
        } else {
            files_in_iss <- unlist(fs::dir_ls(cur$stats_path))
            files_in_iss <- unlist(stringr::str_extract_all(
                files_in_iss, file_reg
            ))
            which_match <- unlist(stringr::str_detect(files_in_iss, pattern))
            files_in_iss <- files_in_iss[which_match]
            cur %>% tibble::add_column(stats_files = list(files_in_iss))
        }
    })
    stats_paths <- stats_paths %>%
        tidyr::unnest(.data$stats_files) %>%
        tidyr::unite(
            col = "files", .data$stats_path, .data$stats_files,
            sep = "/", na.rm = TRUE
        )
}

#' Imports all found Vispa2 stats files.
#'
#' @param association_file The association file
#' @import BiocParallel
#' @importFrom tibble as_tibble
#' @importFrom purrr map2_lgl reduce is_empty
#' @importFrom dplyr mutate bind_rows distinct
#' @keywords internal
#'
#' @return A list with the imported stats and a report of imported files. If no
#' files were imported returns NULL instead
.import_stats_iss <- function(association_file) {
    # Obtain paths
    stats_paths <- .stats_report(association_file)
    if (is.null(stats_paths)) {
        return(NULL)
    }
    # Setup parallel workers and import
    # Register backend according to platform
    if (.Platform$OS.type == "windows") {
        p <- BiocParallel::SnowParam(stop.on.error = FALSE)
    } else {
        p <- BiocParallel::MulticoreParam(
            stop.on.error = FALSE
        )
    }
    FUN <- function(x) {
        stats <- read.csv(
            file = x, sep = "\t", stringsAsFactors = FALSE,
            check.names = FALSE
        )
        stats <- tibble::as_tibble(stats)
        ok <- .check_stats(stats)
        if (ok == TRUE) {
            return(stats)
        } else {
            return(NULL)
        }
    }
    suppressMessages(suppressWarnings({
        stats_dfs <- BiocParallel::bptry(
            BiocParallel::bplapply(stats_paths$files, FUN, BPPARAM = p)
        )
    }))
    BiocParallel::bpstop(p)
    correct <- BiocParallel::bpok(stats_dfs)
    correct <- purrr::map2_lgl(stats_dfs, correct, function(x, y) {
        if (is.null(x)) {
            FALSE
        } else {
            y
        }
    })
    stats_paths <- stats_paths %>% dplyr::mutate(Imported = correct)
    stats_dfs <- stats_dfs[correct]
    # Bind rows in single tibble for all files
    if (purrr::is_empty(stats_dfs)) {
        return(NULL)
    }
    stats_dfs <- purrr::reduce(stats_dfs, function(x, y) {
        x %>%
            dplyr::bind_rows(y) %>%
            dplyr::distinct()
    })
    list(stats_dfs, stats_paths)
}

#' Performs eventual join with stats and aggregation.
#'
#' @param association_file The imported association file
#' @param stats The imported stats df
#' @param grouping_keys The names of the columns to group by
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom stringr str_replace_all
#' @importFrom tidyr unite
#' @keywords internal
#'
#' @return A tibble
.join_and_aggregate <- function(association_file, stats, grouping_keys) {
    # If stats not null join with af
    if (!is.null(stats)) {
        stats <- stats %>% dplyr::mutate(TAG = stringr::str_replace_all(
            .data$TAG,
            pattern = "\\.", replacement = ""
        ))
        association_file <- association_file %>%
            dplyr::left_join(stats,
                by = c(
                    "concatenatePoolIDSeqRun" = "POOL",
                    "TagSequence" = "TAG"
                )
            )
        aggregated <- association_file %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(grouping_keys))) %>%
            dplyr::summarise(
                FusionPrimerPCRDate = suppressWarnings({
                    min(.data$FusionPrimerPCRDate,
                        na.rm = TRUE
                    )
                }),
                LinearPCRDate = suppressWarnings({
                    min(.data$LinearPCRDate, na.rm = TRUE)
                }),
                VCN = suppressWarnings({
                    mean(.data$VCN, na.rm = TRUE)
                }),
                Avg_DNAngUsed = suppressWarnings({
                    mean(.data$DNAngUsed, na.rm = TRUE)
                }),
                Kapa = suppressWarnings({
                    mean(.data$Kapa, na.rm = TRUE)
                }),
                DNAngUsed = suppressWarnings({
                    sum(.data$DNAngUsed, na.rm = TRUE)
                }),
                ulForPool = suppressWarnings({
                    sum(.data$ulForPool, na.rm = TRUE)
                }),
                BARCODE_MUX = suppressWarnings({
                    sum(.data$BARCODE_MUX, na.rm = TRUE)
                }),
                TRIMMING_FINAL_LTRLC = suppressWarnings({
                    sum(.data$TRIMMING_FINAL_LTRLC, na.rm = TRUE)
                }),
                LV_MAPPED = suppressWarnings({
                    sum(.data$LV_MAPPED, na.rm = TRUE)
                }),
                BWA_MAPPED_OVERALL = suppressWarnings({
                    sum(.data$BWA_MAPPED_OVERALL, na.rm = TRUE)
                }),
                ISS_MAPPED_PP = suppressWarnings({
                    sum(.data$ISS_MAPPED_PP, na.rm = TRUE)
                })
            ) %>%
            dplyr::ungroup() %>%
            tidyr::unite(
                col = "AggregateMeta", dplyr::all_of(grouping_keys),
                sep = "_", remove = FALSE
            )
        return(aggregated)
    }
    # Aggregate
    aggregated <- association_file %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(grouping_keys))) %>%
        dplyr::summarise(
            FusionPrimerPCRDate = suppressWarnings({
                min(.data$FusionPrimerPCRDate,
                    na.rm = TRUE
                )
            }),
            LinearPCRDate = suppressWarnings({
                min(.data$LinearPCRDate, na.rm = TRUE)
            }),
            VCN = suppressWarnings({
                mean(.data$VCN, na.rm = TRUE)
            }),
            Avg_DNAngUsed = suppressWarnings({
                mean(.data$DNAngUsed, na.rm = TRUE)
            }),
            Kapa = suppressWarnings({
                mean(.data$Kapa, na.rm = TRUE)
            }),
            DNAngUsed = suppressWarnings({
                sum(.data$DNAngUsed, na.rm = TRUE)
            }),
            ulForPool = suppressWarnings({
                sum(.data$ulForPool, na.rm = TRUE)
            })
        ) %>%
        dplyr::ungroup() %>%
        tidyr::unite(
            col = "AggregateMeta", dplyr::all_of(grouping_keys),
            sep = "_", remove = FALSE
        )
}

#---- USED IN : aggregate_values_by_key ----

#' Internal function that performs aggregation on values with a lambda
#' operation.
#'
#' @param x The list of matrices to aggregate. If a single matrix has to be
#' supplied it must be enclosed in a list. For example `x = list(matrix)`.
#' @param af The association file
#' @param key A string or character vector to use as key
#' @param lambda The aggregating operation to apply to values. Must take as
#' input a numeric/integer vector and return a single value
#' @param group The additional variables to add to grouping
#' @param args Additional arguments passed on to lambda (named list)
#' @param namespace The namespace from which the lambda is exported
#' @param envir The environment in which symbols must be evaluated
#' @import dplyr
#' @importFrom rlang .data
#' @keywords internal
#'
#' @return A list of tibbles with aggregated values
.aggregate_lambda <- function(x, af, key, lambda, group,
    args, namespace, envir) {
    # Vectorize
    aggregated_matrices <- purrr::map(x, function(y) {
        cols <- c(colnames(y), key)
        agg <- y %>%
            dplyr::left_join(af, by = "CompleteAmplificationID") %>%
            dplyr::select(dplyr::all_of(cols)) %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(c(group, key)))) %>%
            dplyr::summarise(Aggregated_value = .wrapper(
                x = .data$Value,
                lambda, args,
                namespace, envir
            )) %>%
            dplyr::ungroup()
        agg
    })
    aggregated_matrices
}


#' Wrapper function to call in summarise.
#'
#' This function tests for the type of result obtained by calling the lambda on
#' the grouped `Value` column: if its an atomic type (and not a list or a
#' vector) simply returns the result, otherwise encapsulates the result in a
#' list to produce a list column that can be unnested with tidyr. Data frames
#' are automatically converted to tibbles to improve performance.
#'
#' This function actually improves both performance and flexibility of the
#' aggregate function since every object of any class can be included in a list
#' and therefore in the aggregated_value column.
#'
#' @param x The numeric vector to which the lambda is applied
#' @param lambda The name of the function to call
#' @param args Additional arguments to pass to lambda
#' @param namespace The namespace that exports lambda (or NULL)
#' @param env The environment where symbols are evaluated
#' @keywords internal
#'
#' @return An atomic type or a list
.wrapper <- function(x, lambda, args, namespace, env) {
    fun <- if (is.null(namespace)) {
        lambda
    } else {
        getExportedValue(ns = namespace, name = lambda)
    }
    res <- do.call(
        what = fun, args = append(list(x), args), quote = TRUE,
        envir = env
    )
    # Check result type, if is anything else than basic types put it in a list
    if (is.atomic(res) & !is.list(res) & length(res) == 1) {
        return(res)
    } else {
        if (is.data.frame(res)) {
            res <- tibble::as_tibble(res)
            return(list(res))
        }
        return(list(res))
    }
}

#### ---- Internals for re-calibration functions ----####

#---- USED IN : compute_near_integrations ----

#' Internal function that computes distance between two loci.
#'
#' @param x First locus
#' @param y Second locus
#' @keywords internal
#'
#' @return The distance as the absolute value of x-y
.locus_distance <- function(x, y) {
    abs((x - y))
}

#' Internal function for keep criteria "keep_first".
#'
#' Returns the first value in a set to keep or returns values to drop if
#' inverted is TRUE.
#'
#' @param x A tibble or a vector of numbers
#' @param inverted Return the values to drop instead of the ones to keep?
#' @importFrom tibble is_tibble
#' @keywords internal
#'
#' @return If x is tibble returns a tibble otherwise returns a vector
.keep_first <- function(x, inverted = TRUE) {
    if (tibble::is_tibble(x)) {
        sum_val <- sum(x$Value)
        result <- x[1, ]
        result$Value <- sum_val
        return(result)
    }
    if (is.integer(x) | is.numeric(x)) {
        if (inverted == TRUE) {
            return(x[c(-1)])
        } else {
            return(x[1])
        }
    }
}

#' Internal function for keep criteria "keep_central".
#'
#' Returns the central value in a set to keep or returns values to drop if
#' inverted is TRUE.
#'
#' @param x A tibble or a vector of numbers
#' @param inverted Return the values to drop instead of the ones to keep?
#' @importFrom tibble is_tibble
#' @keywords internal
#'
#' @return If x is tibble returns a tibble otherwise returns a vector
.keep_central <- function(x, inverted = TRUE) {
    if (tibble::is_tibble(x)) {
        sum_val <- sum(x$Value)
        result <- x[2, ]
        result$Value <- sum_val
        return(result)
    }
    if (is.integer(x) | is.numeric(x)) {
        if (inverted == TRUE) {
            return(x[c(-2)])
        } else {
            return(x[2])
        }
    }
}

#' Internal function for keep criteria "max_value".
#'
#' Returns the value in a set which has (or is) the max "Value" or returns
#' values to drop if inverted is TRUE. If values aren't comparable (because
#' they're equal), a second criteria is used for selection.
#'
#' @param x A tibble or a vector of numbers
#' @param inverted Return the values to drop instead of the ones to keep?
#' @param second_choice Second criteria to use if "max_value" is not applicable
#' @importFrom tibble is_tibble
#' @keywords internal
#'
#' @return If x is tibble returns a tibble otherwise returns a vector
.keep_max_value <- function(x, second_choice, inverted = TRUE) {
    if (tibble::is_tibble(x)) {
        # Check if values are comparable (not equal)
        if (length(unique(x$Value)) == length(x$Value)) {
            sum_val <- sum(x$Value)
            index_to_keep <- which(x$Value == max(x$Value))
            result <- x[index_to_keep, ]
            result$Value <- sum_val
            return(result)
        }
        # If equal use second criteria
        if (second_choice == "keep_central") {
            return(.keep_central(x, inverted = FALSE))
        }
        if (second_choice == "keep_first") {
            return(.keep_first(x, inverted = FALSE))
        }
    }
    if (is.integer(x) | is.numeric(x)) {
        # Check if values are comparable (not equal)
        if (length(unique(x)) == length(x)) {
            if (inverted == TRUE) {
                to_drop <- which(x != max(x))
                return(to_drop)
            } else {
                to_keep <- which(x == max(x))
                return(to_keep)
            }
        }
        # If equal use second criteria
        if (second_choice == "keep_central") {
            return(.keep_central(x, inverted = inverted))
        }
        if (second_choice == "keep_first") {
            return(.keep_first(x, inverted = inverted))
        }
    }
}


#' Wrapper function for criteria evaluation inside `window` function.
#'
#' @param criterias The vector of criterias
#' @param values The numeric vector containing values
#' @param subset Integer vector with indexes of the values vector to pass to
#' individual functions.
#' @keywords internal
#'
#' @return Values to drop
.check_window_criteria <- function(criterias, values, subset) {
    if (criterias[1] == "max_value") {
        if (is.null(subset)) {
            to_drop <- .keep_max_value(values, criterias[2], inverted = TRUE)
        } else {
            to_drop <- .keep_max_value(values[subset], criterias[2],
                inverted = TRUE
            )
        }
        return(to_drop)
    }
    if (criterias[1] == "keep_first") {
        if (is.null(subset)) {
            to_drop <- .keep_first(values, inverted = TRUE)
        } else {
            to_drop <- .keep_first(values[subset], inverted = TRUE)
        }
        return(to_drop)
    }
    if (criterias[1] == "keep_central") {
        if (is.null(subset)) {
            to_drop <- .keep_central(values, inverted = TRUE)
        } else {
            to_drop <- .keep_central(values[subset], inverted = TRUE)
        }
        return(to_drop)
    }
}


#' Internal for the construction of a window of width equal to 3 rows and
#' computation of result on those values.
#'
#' NOTE: to use the function correctly, the parameters
#' `indexes`, `loci` and
#' `values` MUST be named vector with these names:
#' c(first = ..., center = ..., last = ...)
#'
#' @param indexes The indexes of the rows in the window
#' (named vector, see description)
#' @param loci The actual value of the `integration_locus` variable
#' for the rows (named vector, see description)
#' @param values The actual value of the `Value` variable for the rows
#' (named vector, see description)
#' @param criterias The character vector containing the selection criterias
#' @param treshold The numeric value representing the threshold for selection
#' @keywords internal
#'
#' @return A named list with the indexes of rows to drop, the index
#' of the row to collapse on and the value to assign.
#' If no row needs to be dropped returns NULL instead.
.window <- function(indexes, loci, values, criterias, threshold) {
    # Compute distances from center
    dist_x <- .locus_distance(loci["first"], loci["center"])
    dist_y <- .locus_distance(loci["last"], loci["center"])
    # If both above threshold
    if (all(!c(dist_x, dist_y) < threshold)) {
        ## Drop nothing
        return(NULL)
    }
    # If both below threshold
    if (all(c(dist_x, dist_y) < threshold)) {
        ## Check criteria
        to_drop <- .check_window_criteria(criterias, values, subset = NULL)
        indexes_to_drop <- indexes[names(to_drop)]
        collapse_on <- indexes[!names(indexes) %in%
            names(indexes_to_drop)]
        tot_value <- sum(values)
        res <- list(
            drop = indexes_to_drop, collapse_on = collapse_on,
            value = tot_value
        )
        return(res)
    }
    # Only one below
    if (dist_x < threshold) {
        ## Check criteria
        to_drop <- .check_window_criteria(criterias, values, subset = c(1, 2))
        indexes_to_drop <- indexes[names(to_drop)]
        collapse_on <- indexes[c(1, 2)]
        collapse_on <- collapse_on[!names(collapse_on) %in%
            names(indexes_to_drop)]
        tot_value <- sum(values[c(1, 2)])
        res <- list(
            drop = indexes_to_drop, collapse_on = collapse_on,
            value = tot_value
        )
        return(res)
    }
    if (dist_y < threshold) {
        ## Check criteria
        to_drop <- .check_window_criteria(criterias, values, subset = c(2, 3))
        indexes_to_drop <- indexes[names(to_drop)]
        collapse_on <- indexes[c(2, 3)]
        collapse_on <- collapse_on[!names(collapse_on) %in%
            names(indexes_to_drop)]
        tot_value <- sum(values[c(2, 3)])
        res <- list(
            drop = indexes_to_drop, collapse_on = collapse_on,
            value = tot_value
        )
        return(res)
    }
}

#' Internal that represents the window moving and modifying the
#' original table as it progresses.
#'
#' @param x The original tibble (subgroup)
#' @param start The starting index from where the window should be built
#' @param criterias The character vector with selecting criterias
#' @param threshold The numeric value representing the threshold for selection
#' @keywords internal
#'
#' @return A modified tibble with a number of rows less or equal than x.
.window_slide <- function(x, start, criterias, threshold) {
    repeat {
        center <- start + 1
        last <- start + 2
        result <- .window(
            indexes = c(
                first = start,
                center = center,
                last = last
            ),
            loci = c(
                first = x$integration_locus[start],
                center = x$integration_locus[center],
                last = x$integration_locus[last]
            ),
            values = c(
                first = x$Value[start],
                center = x$Value[center],
                last = x$Value[last]
            ),
            criterias = criterias,
            threshold = threshold
        )
        if (last == nrow(x)) {
            if (!is.null(result)) {
                x[result$collapse_on, ]$Value <- result$value
                x <- x[-result$drop, ]
            }
            return(x)
        } else {
            if (!is.null(result)) {
                x[result$collapse_on, ]$Value <- result$value
                x <- x[-result$drop, ]
            }
            start <- if (is.null(result)) {
                start + 1
            } else {
                start
            }
            next
        }
    }
}

#' Internal function that implements the sliding window algorithm.
#'
#' **NOTE: this function is meant to be called on a SINGLE GROUP, meaning a
#' subset of an integration matrix in which all rows share the same chr and same
#' strand.**
#' Also note that is better to call this function on a group that has all
#' distinct integration_locus values (this is ensured by calling function) and
#' also it's not possible to call this function on groups with a single row.
#'
#' @param x An integration matrix subset (see description)
#' @param threshold The numeric value representing the threshold for selection
#' @param keep_criteria The character vector with selecting criterias
#' @importFrom dplyr arrange
#' @keywords internal
#'
#' @return A modified tibble with a number of rows less or equal than x.
.sliding_window <- function(x, threshold, keep_criteria) {
    ## Order by integration_locus
    x <- x %>% dplyr::arrange(.data$integration_locus)
    # ---- If group has only 2 rows ---- #
    if (nrow(x) == 2) {
        ## Compute distance
        distance <- .locus_distance(
            x$integration_locus[1],
            x$integration_locus[2]
        )
        ## Is distance less than the treshold?
        if (distance < threshold) {
            ## If yes examine keep criterias
            if (keep_criteria[1] == "max_value") {
                return(.keep_max_value(x, keep_criteria[2]))
            }
            if (keep_criteria[1] == "keep_first") {
                return(.keep_first(x))
            }
            if (keep_criteria[1] == "keep_central") {
                return(.keep_central(x))
            }
        } else {
            ## If not return the input as is
            return(x)
        }
    }
    # ---- If group has 3 or more rows ---- #
    x <- .window_slide(x,
        start = 1, criterias = keep_criteria,
        threshold = threshold
    )
    return(x)
}
