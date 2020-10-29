#### ---- Internals for HTML widgets construction ----####

# Builds the html widget for the checker table.
#
# @param checker_df Tibble obtained via `.check_file_system_alignment`
# @keywords internal
#' @importFrom reactable reactable reactableTheme colDef
#' @importFrom htmltools div span h2 css browsable
#
# @return An html widget
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

# Builds the html widget for the files_found table.
#
# @param files_found Tibble obtained via `.lookup_matrices` or
# `.lookup_matrices_auto`
# @keywords internal
#' @importFrom reactable reactable reactableTheme colDef
#' @importFrom htmltools div span h2 css h3 browsable
#' @importFrom dplyr select
#' @importFrom rlang .data
#
# @return An html widget
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

# Builds the html widget for the files_to_import table.
#
# @param files_to_import Tibble obtained via
# `.manage_anomalies_interactive` or
# `.manage_anomalies_auto`
# @keywords internal
#' @importFrom reactable reactable reactableTheme colDef
#' @importFrom htmltools div span h2 css browsable
#
# @return An html widget
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

# Builds the html widget for the files_imported table.
#
# @param files_imported Tibble obtained via `.parallel_import_merge`
# @keywords internal
#' @importFrom reactable reactable reactableTheme colDef
#' @importFrom htmltools div span h2 css browsable
#
# @return An html widget
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

# Builds the html widget for the summary table.
#
# @param removed Number of removed collisions
# @param reassigned Number of re-assigned collisions
# @param summary Summary table
# @param tot_rows Total number rows of sequence count matrix before processing
# @param collision_rows Total number of rows of collisions
# @keywords internal
#' @importFrom reactable reactable reactableTheme colDef
#' @importFrom htmltools div h2 h4 css browsable
#
# @return A widget
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

# Builds the html widget for the iss_import.
#
# @param report Table obtained via `import_stats_iss`
# @keywords internal
#
#' @importFrom reactable reactable reactableTheme colDef
#' @importFrom htmltools div h2 h4 css browsable
#
# @return A widget
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

# Builds the html widget for the recalibration map.
#
# @param map_rec Table obtained via `.sliding_window`
# @keywords internal
#
#' @importFrom reactable reactable reactableTheme colDef
#' @importFrom htmltools div h2 css browsable
#
# @return A widget
.recalibr_map_widget <- function(map_rec) {
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
        map_rec,
        striped = TRUE,
        sortable = TRUE,
        showSortable = TRUE,
        bordered = FALSE,
        outlined = TRUE,
        searchable = TRUE,
        pagination = TRUE,
        paginationType = "simple",
        showPageSizeOptions = TRUE,
        pageSizeOptions = c(5, 10, 15),
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
    widget <- htmltools::div(
        style = htmltools::css(font.family = "Calibri"),
        htmltools::h2("RECALIBRATION MAP"),
        styled_df
    )
    htmltools::browsable(widget)
}

# Builds the html widget for the missing info in collisions.
#
# @param map_rec Table obtained in remove_collisions
# @keywords internal
#
#' @importFrom reactable reactable reactableTheme colDef
#' @importFrom htmltools div h2 css browsable span
#
# @return A widget
.missing_info_widget <- function(missing) {
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
        missing,
        striped = TRUE,
        sortable = TRUE,
        showSortable = TRUE,
        bordered = FALSE,
        outlined = TRUE,
        searchable = TRUE,
        pagination = TRUE,
        paginationType = "simple",
        showPageSizeOptions = TRUE,
        pageSizeOptions = c(5, 10, 15),
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
    widget <- htmltools::div(
        style = htmltools::css(font.family = "Calibri"),
        htmltools::h2("MISSING INFORMATION"),
        htmltools::span(paste(
            "All samples in this table are removed from",
            "the matrix"
        )),
        styled_df
    )
    htmltools::browsable(widget)
}

# Builds the html widget for the additional info in collisions.
#
# @param map_rec Table obtained in remove_collisions
# @keywords internal
#
#' @importFrom reactable reactable reactableTheme colDef
#' @importFrom htmltools div h2 css browsable span
#
# @return A widget
.add_info_widget <- function(additional) {
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
        additional,
        striped = TRUE,
        sortable = TRUE,
        showSortable = TRUE,
        bordered = FALSE,
        outlined = TRUE,
        searchable = TRUE,
        pagination = TRUE,
        paginationType = "simple",
        showPageSizeOptions = TRUE,
        pageSizeOptions = c(5, 10, 15),
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
    widget <- htmltools::div(
        style = htmltools::css(font.family = "Calibri"),
        htmltools::h2("ADDITIONAL INFO FOUND"),
        htmltools::span(paste(
            "Additional info found in metadata that",
            "is not included in provided matrix"
        )),
        styled_df
    )
    htmltools::browsable(widget)
}
