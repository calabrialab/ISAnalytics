#### ---- Internals for HTML widgets construction ----####

# Default theme for reactable widgets
#' @importFrom reactable reactableTheme
.default_theme <- function() {
    reactable::reactableTheme(
        style = list(
            fontFamily = "Calibri"
        ),
        cellStyle = list(
            display = "flex",
            flexDirection = "column",
            justifyContent = "center"
        )
    )
}

# Generates a standardized reactable widget for the input df
#' @importFrom reactable reactable colDef
.generate_react_table <- function(df, ...) {
    styled_df <- reactable::reactable(
        df,
        striped = TRUE,
        sortable = TRUE,
        showSortable = TRUE,
        showSortIcon = TRUE,
        filterable = TRUE,
        bordered = FALSE,
        outlined = TRUE,
        searchable = TRUE,
        pagination = TRUE,
        paginationType = "numbers",
        showPageSizeOptions = TRUE,
        pageSizeOptions = c(5, 10, 15),
        defaultPageSize = 5,
        showPagination = TRUE,
        resizable = TRUE,
        theme = .default_theme(),
        defaultColDef = reactable::colDef(
            headerStyle = list(fontSize = "18px", paddingLeft = "15px"),
            align = "left", sortable = TRUE, resizable = TRUE,
            filterable = TRUE, style = list(paddingLeft = "15px"),
            header = function(value) gsub("_", " ", value, fixed = TRUE)
        ),
        ...
    )
}

# Generates a reactable formatted for the sharing matrix (conditional coloring)
#' @importFrom reactable reactable colDef colFormat
#' @importFrom grDevices rgb colorRamp
.generate_react_table_sharing_colored <- function(df, digits, ...) {
    pal <- function(x) {
        grDevices::rgb(grDevices::colorRamp(c("white", "gold", "navyblue"))(x),
            maxColorValue = 255
        )
    }
    styled_df <- reactable::reactable(
        df,
        striped = FALSE,
        sortable = FALSE,
        filterable = FALSE,
        bordered = FALSE,
        outlined = TRUE,
        searchable = FALSE,
        pagination = FALSE,
        resizable = TRUE,
        theme = .default_theme(),
        defaultColDef = reactable::colDef(
            headerStyle = list(fontSize = "18px", paddingLeft = "15px"),
            align = "center", sortable = FALSE, resizable = TRUE,
            filterable = FALSE,
            header = function(value) gsub("_", " ", value, fixed = TRUE),
            format = reactable::colFormat(
                digits = digits,
                locales = "en-US"
            ),
            style = function(value) {
                only_num <- df[colnames(df) != "Project_Subject"]
                if (!is.numeric(value)) {
                    return()
                }
                normalized <- (value - min(only_num)) /
                    (max(only_num) - min(only_num))
                color <- pal(normalized)
                text_color <- if (color == "#000080") {
                    "white"
                } else {
                    "black"
                }
                list(background = color, color = text_color)
            }
        ),
        columns = list(
            Project_Subject = reactable::colDef(
                name = "",
                style = list(fontWeight = "bold", fontSize = "18px")
            )
        ),
        ...
    )
}

# Generates a smaller, more static reactable widget with no pagination
# and no filtering options (ideal for displaying single numeric values like
# quantification sums)
#' @importFrom reactable reactable colDef colFormat
#' @importFrom purrr set_names
.generate_react_table_mini <- function(df, perc_cols, ...) {
    perc_col <- reactable::colDef(
        format = reactable::colFormat(
            percent = TRUE,
            digits = 2
        )
    )

    c_perc <- purrr::map(perc_cols, function(x) {
        perc_col
    }) %>% purrr::set_names(perc_cols)

    styled_df <- reactable::reactable(
        df,
        fullWidth = FALSE,
        bordered = FALSE,
        outlined = TRUE,
        theme = .default_theme(),
        defaultColDef = reactable::colDef(
            headerStyle = list(fontSize = "18px", paddingLeft = "15px"),
            align = "left", sortable = FALSE, resizable = TRUE,
            filterable = FALSE, style = list(paddingLeft = "15px"),
            header = function(value) gsub("_", " ", value, fixed = TRUE)
        ),
        columns = c_perc,
        ...
    )
}

# Contains general css style for html widgets
.widget_css <- function() {
    paste(
        'body {font-family: "Calibri"}',
        "h2 {margin-left: 10px; font-weight: bold; border-bottom-width: 1px;
    border-bottom-style: solid;}",
        "#subtitle {font-size: 18px;}",
        "#section-content {margin-left: 20px;}",
        "#single-numbers {margin-left: 25px;}"
    )
}
### ---- import_association_file ---- ###

# Builds the html widget for summary of import_association_file.
#
# @param checker_df Tibble obtained via `.check_file_system_alignment`
# @keywords internal
#' @importFrom reactable colDef
#' @importFrom htmltools tags h1 h2 div browsable
#
# @return An html widget
.checker_widget <- function(checker_df) {
    columns_def <- list(
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
            minWidth = 100
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
            minWidth = 200
        )
    )

    styled_df <- .generate_react_table(checker_df,
        defaultSorted = list(Found = "asc"),
        columns = columns_def
    )


    widget <- htmltools::tags$html(
        htmltools::tags$head(
            htmltools::tags$style(.widget_css())
        ),
        htmltools::tags$body(
            htmltools::h1("IMPORT ASSOCIATION FILE REPORT"),
            htmltools::h2("ALIGNMENT RESULTS"),
            htmltools::div(
                id = "section-content",
                htmltools::div("Results of alignment between file system and",
                    "association file. If some folders are not found",
                    "they will be ignored until the problem is fixed",
                    "and the association file re-imported.",
                    id = "subtitle"
                )
            ),
            styled_df
        )
    )

    htmltools::browsable(widget)
}

### ---- import_parallel_Vispa2Matrices ---- ###

# Report widgets for parallel import of matrices

#' @importFrom htmltools div h3 browsable
.files_found_details_row <- function(index, files_found) {
    files <- files_found$Files[[index]]
    count <- files_found$Files_count[[index]]
    styled_files <- .files_found_files(files)
    styled_count <- .files_found_count(count)
    w <- htmltools::div(
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
    return(w)
}

#' @importFrom reactable colDef
.files_found_files <- function(files) {
    styled_files <- .generate_react_table(files,
        columns = list(
            Quantification_type = reactable::colDef(
                minWidth = 200,
                maxWidth = 200
            )
        )
    )
    return(styled_files)
}

#' @importFrom reactable colDef
.files_found_count <- function(count) {
    col_def <- list(
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

    styled_count <- .generate_react_table(count, columns = col_def)
    return(styled_count)
}

# Builds the html widget for the files_found table.
#
# @param files_found Tibble obtained via `.lookup_matrices` or
# `.lookup_matrices_auto`
# @keywords internal
#' @importFrom reactable colDef
#' @importFrom htmltools tags h2 div browsable
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

    cols_def <- list(
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
    )

    styled_df <- .generate_react_table(main_cols,
        columns = cols_def,
        details = function(index) {
            .files_found_details_row(
                index,
                files_found
            )
        }
    )

    widget <- htmltools::tags$html(
        htmltools::tags$head(
            htmltools::tags$style(.widget_css())
        ),
        htmltools::tags$body(
            htmltools::h2("INTEGRATION MATRICES FOUND REPORT"),
            htmltools::div(
                id = "section-content",
                htmltools::div("Report of all files found for
                               each quantification",
                    "type. Click on the arrow on the left side of each",
                    "row to see details.",
                    id = "subtitle"
                ),
                styled_df
            )
        )
    )
    htmltools::browsable(widget)
}

# Builds the html widget for the files_to_import table.
#
# @param files_to_import Tibble obtained via
# `.manage_anomalies_interactive` or
# `.manage_anomalies_auto`
# @keywords internal
#' @importFrom reactable colDef
#' @importFrom htmltools tags h2 div browsable
#
# @return An html widget
.files_to_import_widget <- function(files_to_import) {
    cols_def <- list(
        Files_chosen = reactable::colDef(
            minWidth = 250
        ),
        Quantification_type = reactable::colDef(
            align = "center"
        )
    )

    styled_df <- .generate_react_table(files_to_import, columns = cols_def)

    widget <- htmltools::tags$html(
        htmltools::tags$head(
            htmltools::tags$style(.widget_css())
        ),
        htmltools::tags$body(
            htmltools::h2("SUMMARY OF FILES CHOSEN FOR IMPORT"),
            htmltools::div(
                id = "section-content",
                htmltools::div("Here is a summary of all files
                               chosen for import",
                    id = "subtitle"
                ),
                styled_df
            )
        )
    )
    htmltools::browsable(widget)
}

# Builds the html widget for the files_imported table.
#
# @param files_imported Tibble obtained via `.parallel_import_merge`
# @keywords internal
#' @importFrom reactable colDef
#' @importFrom htmltools tags h2 div browsable
#
# @return An html widget
.files_imported_widget <- function(files_imported) {
    cols_def <- list(
        Files_chosen = reactable::colDef(
            minWidth = 250
        ),
        Quantification_type = reactable::colDef(
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

    styled_df <- .generate_react_table(files_imported, columns = cols_def)

    widget <- htmltools::tags$html(
        htmltools::tags$head(
            htmltools::tags$style(.widget_css())
        ),
        htmltools::tags$body(
            htmltools::h2("REPORT: FILES IMPORTED"),
            htmltools::div(
                id = "section-content",
                htmltools::div("Here is a summary of all files
                actually imported for
        each quantification type. If you see 'false' in the column Imported,
        some errors might have occurred and the function was unable to import
        that matrix.", id = "subtitle"),
                styled_df
            )
        )
    )
    htmltools::browsable(widget)
}

#' @importFrom htmltools tags h1 h2 browsable
.import_report_widget <- function(files_found, files_to_import,
    files_imported) {
    files_found_styled <- .files_found_widget(files_found)
    files_toimp_styled <- .files_to_import_widget(files_to_import)
    files_imported_styled <- .files_imported_widget(files_imported)

    widget <- htmltools::tags$html(
        htmltools::tags$head(
            htmltools::tags$style(.widget_css())
        ),
        htmltools::tags$body(
            htmltools::h1("MATRICES IMPORT REPORT"),
            files_found_styled,
            files_toimp_styled,
            files_imported_styled
        )
    )
    htmltools::browsable(widget)
}

### ---- remove_collisions ---- ###

# Section 1 of the report: info on the input seq count or multi quant matrix
# (before collision removal)

## Before joining with metadata
#' @import dplyr
#' @importFrom tidyr pivot_longer everything
#' @importFrom htmltools tags h2 div h4 browsable
.sc_stats_input <- function(input, quant_cols) {
    quant_totals <- input %>%
        dplyr::select(dplyr::all_of(quant_cols)) %>%
        dplyr::summarise(dplyr::across(
            .cols = dplyr::all_of(quant_cols),
            .fns = list(
                sum = ~ sum(.x, na.rm = TRUE)
            )
        ),
        .groups = "drop"
        )
    quant_totals_pivoted <- quant_totals %>%
        tidyr::pivot_longer(tidyr::everything(),
            names_to = c("Quantification", ".value"),
            names_pattern = paste0(
                "(",
                paste0(quant_cols,
                    collapse = "|"
                ),
                ")", "_(.+)"
            ),
            names_transform = list(
                Quantification = ~ readr::parse_factor(
                    .x,
                    levels = quant_cols
                )
            )
        )
    styled_quant_totals <- .generate_react_table(quant_totals_pivoted)
    tot_iss <- input %>%
        dplyr::select(dplyr::all_of(mandatory_IS_vars())) %>%
        dplyr::distinct() %>%
        nrow()
    widget <- htmltools::tags$html(
        htmltools::tags$head(
            htmltools::tags$style(.widget_css())
        ),
        htmltools::tags$body(
            htmltools::h2("INPUT MATRIX INFO"),
            htmltools::div(
                id = "section-content",
                htmltools::div("Information on the matrix provided as input
                        before processing collisions",
                    id = "subtitle"
                ),
                htmltools::h4(paste(
                    "TOTAL ISS (number of distinct",
                    "integration sites):"
                )),
                htmltools::div(tot_iss, id = "single-numbers"),
                htmltools::h4("QUANTIFICATIONS TOTALS:"),
                styled_quant_totals
            )
        )
    )
    htmltools::browsable(widget)
}

## After identifying missing and join with metadata
#' @import dplyr
#' @importFrom htmltools tags h2 div h4 browsable tagList
#' @importFrom psych describe
#' @importFrom stringr str_replace_all
.sc_stats_input_joined <- function(input_joined, quant_cols) {
    pools_and_samples <- input_joined %>%
        dplyr::select(.data$PoolID, .data$CompleteAmplificationID) %>%
        dplyr::distinct()
    pools_and_samples_styled <- .generate_react_table(pools_and_samples,
        groupBy = "PoolID"
    )
    pool_stats <- input_joined %>%
        dplyr::group_by(.data$PoolID) %>%
        dplyr::summarise(dplyr::across(
            .cols = dplyr::all_of(quant_cols),
            .fns = list(
                sum = ~ sum(.x, na.rm = TRUE),
                count = length, describe = psych::describe
            )
        ),
        .groups = "drop"
        )
    single_cells <- c("PoolID", colnames(pool_stats)[grepl(
        "*_sum$|*_count$",
        colnames(pool_stats)
    )])

    desc_cells <- colnames(pool_stats)[grepl(
        "*_describe$",
        colnames(pool_stats)
    )]
    pool_stats_styled <- .generate_react_table(pool_stats[single_cells],
        details = function(index) {
            sub_sect <- lapply(desc_cells, FUN = function(desc) {
                sub_index <- pool_stats[index, ]
                styled_desc <- .generate_react_table(sub_index[[desc]])
                htmltools::tags$html(
                    htmltools::tags$head(
                        htmltools::tags$style(.widget_css())
                    ),
                    htmltools::tags$body(
                        htmltools::h3(stringr::str_replace_all(desc, "_", " ")),
                        styled_desc
                    )
                )
            })
            htmltools::tagList(sub_sect)
        }
    )

    widget <- htmltools::tags$html(
        htmltools::tags$head(
            htmltools::tags$style(.widget_css())
        ),
        htmltools::tags$body(
            htmltools::h2("PRE-PROCESS MATRIX INFO"),
            htmltools::div(
                id = "section-content",
                htmltools::div("Information on the pre-process matrix, aka
                               the input matrix with missing samples removed",
                    id = "subtitle"
                ),
                htmltools::h4("SUMMARY OF POOLS AND SAMPLES PRESENT"),
                pools_and_samples_styled,
                htmltools::h4("PER-POOL STATS"),
                pool_stats_styled
            )
        )
    )
    htmltools::browsable(widget)
}

# Builds the html widget for the missing info in collisions.
#
#' @importFrom htmltools tags h2 h4 div browsable
#' @import upsetjs
#' @import dplyr
#
# @return A widget
.missing_info_widget <- function(missing, input, af) {
    # All samples found missing
    styled_missing <- .generate_react_table(input[missing, ])

    # Venn diagram
    samples_input <- input %>%
        dplyr::select(.data$CompleteAmplificationID) %>%
        dplyr::distinct()
    samples_af <- af %>%
        dplyr::select(.data$CompleteAmplificationID) %>%
        dplyr::distinct()
    venn_data <- list(
        "Input matrix" = samples_input$CompleteAmplificationID,
        "Association file" = samples_af$CompleteAmplificationID
    )
    colors <- list(
        "Input matrix" = "#1f77b4", "Association file" = "#9467bd",
        "Input matrix&Association file" = "#e377c2"
    )

    venn <- upsetjs::upsetjsVennDiagram()

    venn <- upsetjs::fromList(venn, venn_data, colors = colors)
    venn <- upsetjs::chartTheme(venn,
        selection.color = "",
        has.selection.opacity = 0.3
    )
    venn <- upsetjs::chartVennLabels(venn,
        title = paste(
            "Samples shared between input",
            "and association file"
        ),
        description = paste(
            "The number of distinct",
            "CompleteAmplificationIDs",
            "shared between the association",
            "file and the input matrix"
        )
    )
    venn <- upsetjs::interactiveChart(venn)


    widget <- htmltools::tags$html(
        htmltools::tags$head(
            htmltools::tags$style(.widget_css())
        ),
        htmltools::tags$body(
            htmltools::h2("MISSING INFORMATION"),
            htmltools::div(
                id = "section-content",
                htmltools::div("Information that is marked as missing from the
                       association file provided but was found in the
                       input matrix. All missing samples are removed
                       from the matrix in the pre-process phase.",
                    id = "subtitle"
                ),
                htmltools::h4("MISSING SAMPLES"),
                styled_missing,
                venn
            )
        )
    )
    htmltools::browsable(widget)
}

# Builds the html widget for the additional info in collisions.
#
#' @importFrom htmltools tags h2 div h4 browsable
#
# @return A widget
.add_info_widget <- function(additional) {
    styled_add <- .generate_react_table(additional,
        groupBy = c("ProjectID", "PoolID")
    )

    widget <- htmltools::tags$html(
        htmltools::tags$head(
            htmltools::tags$style(.widget_css())
        ),
        htmltools::tags$body(
            htmltools::h2("ADDITIONAL INFORMATION"),
            htmltools::div(
                id = "section-content",
                htmltools::div("Information that is found in
                the association file
                       but not in the input matrix,
                       which is related to the projects of interest",
                    id = "subtitle"
                ),
                htmltools::h4("ADDITIONAL SAMPLES"),
                styled_add
            )
        )
    )
    htmltools::browsable(widget)
}

#' @import dplyr
#' @importFrom purrr map2_dfr map2 set_names
#' @importFrom tibble as_tibble_row
#' @importFrom htmltools tags h2 h4 div browsable
.sharing_widget <- function(input_joined, phase) {
    identifiers <- input_joined %>%
        dplyr::select(.data$ProjectID, .data$SubjectID) %>%
        dplyr::distinct()

    selected_df <- input_joined %>%
        dplyr::select(
            dplyr::all_of(mandatory_IS_vars()),
            .data$ProjectID, .data$SubjectID
        ) %>%
        dplyr::distinct()

    get_shared_cell <- function(p, s, temp) {
        selected2 <- selected_df %>%
            dplyr::filter(
                .data$ProjectID == p,
                .data$SubjectID == s
            ) %>%
            dplyr::select(dplyr::all_of(mandatory_IS_vars()))
        shared <- temp %>%
            dplyr::select(mandatory_IS_vars()) %>%
            dplyr::semi_join(selected2, by = mandatory_IS_vars()) %>%
            nrow()
    }

    sharing_df_abs <- purrr::map2_dfr(
        identifiers$ProjectID, identifiers$SubjectID,
        function(proj, subj) {
            temp <- selected_df %>%
                dplyr::filter(
                    .data$ProjectID == proj,
                    .data$SubjectID == subj
                )
            cols <- purrr::map2(
                identifiers$ProjectID,
                identifiers$SubjectID,
                function(p, s) {
                    get_shared_cell(p, s, temp)
                }
            ) %>% purrr::set_names(~ paste0(
                identifiers$ProjectID,
                "_",
                identifiers$SubjectID
            ))
            tibble::as_tibble_row(
                purrr::prepend(cols, list("Project_Subject" = paste0(
                    proj,
                    "_",
                    subj
                )))
            )
        }
    )

    sharing_df_rel <- purrr::map2_dfr(
        identifiers$ProjectID, identifiers$SubjectID,
        function(proj, subj) {
            temp <- selected_df %>%
                dplyr::filter(
                    .data$ProjectID == proj,
                    .data$SubjectID == subj
                )
            cols <- purrr::map2(
                identifiers$ProjectID,
                identifiers$SubjectID,
                function(p, s) {
                    shared <- get_shared_cell(p, s, temp)
                    shared_ratio <- shared / nrow(temp)
                }
            ) %>% purrr::set_names(~ paste0(
                identifiers$ProjectID,
                "_",
                identifiers$SubjectID
            ))
            tibble::as_tibble_row(
                purrr::prepend(cols, list("Project_Subject" = paste0(
                    proj,
                    "_",
                    subj
                )))
            )
        }
    )

    styled_shared_abs <- .generate_react_table_sharing_colored(
        sharing_df_abs,
        0
    )
    styled_shared_ratio <- .generate_react_table_sharing_colored(
        sharing_df_rel,
        2
    )
    widget <- htmltools::tags$html(
        htmltools::tags$head(
            htmltools::tags$style(.widget_css())
        ),
        htmltools::tags$body(
            htmltools::h2(paste("SHARING", phase)),
            htmltools::div(
                id = "section-content",
                htmltools::div("ISs sharing between independent samples",
                    id = "subtitle"
                ),
                htmltools::h4("ABSOLUTE VALUES:"),
                styled_shared_abs,
                htmltools::h4("RELATIVE VALUES (LEFT TO RIGHT):"),
                styled_shared_ratio
            )
        )
    )
    htmltools::browsable(widget)
}

#' @importFrom tidyr pivot_longer
#' @importFrom readr parse_factor
#' @import dplyr
#' @importFrom psych describe
#' @importFrom htmltools tags h2 h3 h4 div browsable tagList
#' @importFrom stringr str_replace_all
.sc_stats_after <- function(after, quant_cols, meta_joined) {
    quant_totals <- after %>%
        dplyr::select(dplyr::all_of(quant_cols)) %>%
        dplyr::summarise(dplyr::across(
            .cols = dplyr::all_of(quant_cols),
            .fns = list(
                sum = ~ sum(.x, na.rm = TRUE)
            )
        ),
        .groups = "drop"
        )
    quant_totals_pivoted <- quant_totals %>%
        tidyr::pivot_longer(tidyr::everything(),
            names_to = c("Quantification", ".value"),
            names_pattern = paste0(
                "(",
                paste0(quant_cols,
                    collapse = "|"
                ),
                ")", "_(.+)"
            ),
            names_transform = list(
                Quantification = ~ readr::parse_factor(
                    .x,
                    levels = quant_cols
                )
            )
        )
    styled_quant_totals <- .generate_react_table(quant_totals_pivoted)
    tot_iss <- after %>%
        dplyr::select(dplyr::all_of(mandatory_IS_vars())) %>%
        dplyr::distinct() %>%
        nrow()

    pools_and_samples <- meta_joined %>%
        dplyr::select(.data$PoolID, .data$CompleteAmplificationID) %>%
        dplyr::distinct()
    pools_and_samples_styled <- .generate_react_table(pools_and_samples,
        groupBy = "PoolID"
    )

    pool_stats <- meta_joined %>%
        dplyr::group_by(.data$PoolID) %>%
        dplyr::summarise(dplyr::across(
            .cols = dplyr::all_of(quant_cols),
            .fns = list(
                sum = ~ sum(.x, na.rm = TRUE),
                count = length, describe = psych::describe
            )
        ),
        .groups = "drop"
        )
    single_cells <- c("PoolID", colnames(pool_stats)[grepl(
        "*_sum$|*_count$",
        colnames(pool_stats)
    )])
    desc_cells <- colnames(pool_stats)[grepl(
        "*_describe$",
        colnames(pool_stats)
    )]

    pool_stats_styled <- .generate_react_table(pool_stats[single_cells],
        details = function(index) {
            sub_sect <- lapply(desc_cells, FUN = function(desc) {
                sub_index <- pool_stats[index, ]
                styled_desc <- .generate_react_table(sub_index[[desc]])
                htmltools::tags$html(
                    htmltools::tags$head(
                        htmltools::tags$style(.widget_css())
                    ),
                    htmltools::tags$body(
                        htmltools::h3(stringr::str_replace_all(desc, "_", " ")),
                        styled_desc
                    )
                )
            })
            htmltools::tagList(sub_sect)
        }
    )

    widget <- htmltools::tags$html(
        htmltools::tags$head(
            htmltools::tags$style(.widget_css())
        ),
        htmltools::tags$body(
            htmltools::h2("POST-PROCESS FINAL MATRIX INFO"),
            htmltools::div(
                id = "section-content",
                htmltools::div("Information on the post-process matrix",
                    id = "subtitle"
                ),
                htmltools::h4(paste(
                    "TOTAL ISS (number of distinct",
                    "integration sites):"
                )),
                htmltools::div(tot_iss, id = "single-numbers"),
                htmltools::h4("QUANTIFICATIONS TOTALS:"),
                styled_quant_totals,
                htmltools::h4("SUMMARY OF POOLS AND SAMPLES PRESENT"),
                pools_and_samples_styled,
                htmltools::h4("PER-POOL STATS"),
                pool_stats_styled
            )
        )
    )
    htmltools::browsable(widget)
}

#' @importFrom htmltools tags h1 h2 h4 div browsable
#' @import dplyr
#' @importFrom purrr is_empty
#' @importFrom tibble tibble
# @return A widget
.collisions_widget <- function(
    input_df,
    quant_cols,
    input_joined_df,
    association_file,
    missing,
    add_info,
    collision_df,
    after_df,
    removed,
    reassigned,
    summary) {
    report_input <- .sc_stats_input(input_df, quant_cols)
    missing_info <- if (!purrr::is_empty(missing)) {
        .missing_info_widget(missing, input_df, association_file)
    } else {
        htmltools::p()
    }
    additional_info <- if (!is.null(add_info)) {
        .add_info_widget(add_info)
    } else {
        htmltools::p()
    }
    post_join_info <- .sc_stats_input_joined(input_joined_df, quant_cols)
    pre_process_matrix <- .generate_react_table(input_df[-missing, ])
    pre_sharing <- .sharing_widget(input_joined_df, "PRE-PROCESSING")
    summary_w <- .generate_react_table(summary)

    joined_after <- after_df %>%
        dplyr::left_join(association_file, by = "CompleteAmplificationID") %>%
        dplyr::select(
            dplyr::all_of(colnames(after_df)),
            .data$ProjectID, .data$PoolID, .data$SubjectID,
            .data$ReplicateNumber
        )

    report_post <- .sc_stats_after(after_df, quant_cols, joined_after)
    post_sharing <- .sharing_widget(joined_after, "POST-PROCESSING")


    coll_iss <- collision_df %>%
        dplyr::select(.data$chr, .data$integration_locus, .data$strand) %>%
        dplyr::distinct() %>%
        nrow()

    tot_iss_input <- input_joined_df %>%
        dplyr::select(.data$chr, .data$integration_locus, .data$strand) %>%
        dplyr::distinct() %>%
        nrow()

    collisions_found <- tibble::tibble(
        abs_number = coll_iss,
        percentage_on_total =
            coll_iss / tot_iss_input
    )

    collisions_removed <- tibble::tibble(
        abs_number = removed,
        percentage_on_collisions =
            removed / coll_iss,
        percentage_on_total =
            removed / tot_iss_input
    )

    collisions_reassigned <- tibble::tibble(
        abs_number = reassigned,
        percentage_on_collisions =
            reassigned / coll_iss,
        percentage_on_total =
            reassigned / tot_iss_input
    )

    collisions_found_styled <- .generate_react_table_mini(
        collisions_found,
        "percentage_on_total"
    )

    collisions_removed_styled <- .generate_react_table_mini(
        collisions_removed, c("percentage_on_total", "percentage_on_collisions")
    )

    collisions_reassigned_styled <- .generate_react_table_mini(
        collisions_reassigned, c(
            "percentage_on_total",
            "percentage_on_collisions"
        )
    )

    widget <- htmltools::tags$html(
        htmltools::tags$head(
            htmltools::tags$style(.widget_css())
        ),
        htmltools::tags$body(
            htmltools::h1("COLLISION PROCESSING REPORT"),
            report_input,
            missing_info,
            additional_info,
            post_join_info,
            htmltools::h2("PRE-PROCESS FINAL MATRIX"),
            htmltools::div(
                id = "section-content",
                htmltools::div("This is the matrix that will be processed for
                         collisions",
                    id = "subtitle"
                ),
                pre_process_matrix
            ),
            pre_sharing,
            htmltools::h2("POST-PROCESSING SUMMARY"),
            htmltools::div(
                id = "section-content",
                htmltools::div("A summary of the post-processing",
                    id = "subtitle"
                ),
                summary_w
            ),
            htmltools::h2("COLLISIONS INFO"),
            htmltools::div(
                id = "section-content",
                htmltools::h4("COLLISIONS FOUND"),
                collisions_found_styled,
                htmltools::h4("COLLISIONS REMOVED"),
                collisions_removed_styled,
                htmltools::h4("COLLISIONS REASSIGNED"),
                collisions_reassigned_styled
            ),
            report_post,
            post_sharing
        )
    )
    htmltools::browsable(widget)
}

#---- aggregate_metadata ----#

# Builds the html widget for the iss_import.
#
# @param report Table obtained via `import_stats_iss`
# @keywords internal
#
#' @importFrom reactable colDef
#' @importFrom htmltools tags h1 div browsable
#
# @return A widget
.iss_import_widget <- function(report) {
    cols <- list(
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

    styled_df <- .generate_react_table(report, columns = cols)

    widget <- htmltools::tags$html(
        htmltools::tags$head(
            htmltools::tags$style(.widget_css())
        ),
        htmltools::tags$body(
            htmltools::h1("REPORT IMPORT VISPA2 STATS: FILES IMPORTED"),
            htmltools::div(
                id = "section-content",
                htmltools::div("Here is a summary of all files actually
                imported.
        If you see 'FALSE' in the column Imported, some errors might have
        occurred and the function was unable to import the file or simply no
        path was found for that stats file.",
                    id = "subtitle"
                ),
                styled_df
            )
        )
    )
    htmltools::browsable(widget)
}

# Builds the html widget for the recalibration map.
#
# @param map_rec Table obtained via `.sliding_window`
# @keywords internal
#
#' @importFrom htmltools tags h1 div browsable
#
# @return A widget
.recalibr_map_widget <- function(map_rec) {
    styled_df <- .generate_react_table(map_rec)

    widget <- htmltools::tags$html(
        htmltools::tags$head(
            htmltools::tags$style(.widget_css())
        ),
        htmltools::tags$body(
            htmltools::h1("RECALIBRATION REPORT"),
            htmltools::div(
                id = "section-content",
                htmltools::div("Recalibration map", id = "subtitle"),
                styled_df
            )
        )
    )
    htmltools::browsable(widget)
}
