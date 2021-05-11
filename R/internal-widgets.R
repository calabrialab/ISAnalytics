#### ---- Internals for HTML widgets construction ----####

.produce_widget <- function(widget_function, ...) {
    dots <- rlang::list2(...)
    widg <- withCallingHandlers(expr = {
        withRestarts({
            rlang::exec(.fn = widget_function, !!!dots)
        },
        widget_fail = function() {
            rlang::inform(.widgets_error())
            NULL
        })
    }, error = function(cnd) {
        invokeRestart("widget_fail")
    })
    widg
}

.print_widget <- function(widget, else_verbose = NULL) {
    withCallingHandlers(expr = {
        withRestarts({
            print(widget)
        }, print_err = function() {
            rlang::inform(.widgets_print_error())
        })
    }, error = function(cnd) {
        rlang::inform(conditionMessage(cnd))
        if (getOption("ISAnalytics.verbose") == TRUE) {
            else_verbose
        }
        invokeRestart("print_err")
    })
}

#' @importFrom fs is_dir dir_exists dir_create path
#' @importFrom htmltools save_html
.export_widget_file <- function(widget, path, def_file_name) {
    export_widget_path <- if (fs::is_dir(path)) {
        if (!fs::dir_exists(path)) {
            fs::dir_create(path)
        }
        fs::path(
            path,
            def_file_name
        )
    } else {
        path
    }
    withCallingHandlers(
        expr = {
            htmltools::save_html(widget, export_widget_path)
        },
        error = function(cnd) {
            warning(.widgets_save_error())
        }
    )
}

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
    dots <- rlang::list2(...)
    if (!"defaultColDef" %in% names(dots)) {
        def <- list(defaultColDef = reactable::colDef(
            headerStyle = list(fontSize = "18px", paddingLeft = "15px"),
            align = "left", sortable = TRUE, resizable = TRUE,
            filterable = TRUE, style = list(paddingLeft = "15px"),
            header = function(value) gsub("_", " ", value, fixed = TRUE)
        ))
        dots <- append(dots, def)
    }
    styled_df <- rlang::exec(
        reactable::reactable,
        data = df,
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
        !!!dots
    )
    styled_df
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
        "#subtitle {font-size: 18px; margin-bottom: 10px;}",
        "#section-content {margin-left: 20px;}",
        "#single-numbers {margin-left: 25px;}",
        "h3, h4 {margin-left: 10px;}",
        "#simple_txt {color: grey; font-style: italic;}"
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
.checker_widget <- function(parsing_probs, date_probs, checker_df,
    col_probs, critical_nas) {
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

    nothing_to_rep <- htmltools::div(
        id = "section-content",
        htmltools::div("Nothing to report", id = "simple_txt")
    )

    styled_parsing_df <- if (!is.null(parsing_probs) &&
        !purrr::is_empty(parsing_probs)) {
        .generate_react_table(parsing_probs)
    } else {
        nothing_to_rep
    }

    styled_date_df <- if (!is.null(date_probs) &&
        !purrr::is_empty(date_probs)) {
        .generate_react_table(date_probs)
    } else {
        nothing_to_rep
    }

    styled_checker_df <- if (!is.null(checker_df) &&
        !purrr::is_empty(checker_df)) {
        .generate_react_table(checker_df,
            defaultSorted = list(Found = "asc"),
            columns = columns_def
        )
    } else {
        nothing_to_rep
    }


    missing <- if ("missing" %in% names(col_probs)) {
        paste0(col_probs$missing, collapse = ", ")
    } else {
        "no missing columns to report"
    }
    non_st <- if ("non_standard" %in% names(col_probs)) {
        paste0(col_probs$non_standard, collapse = ", ")
    } else {
        "no non-standard columns to report"
    }
    styled_cols_probs <- htmltools::div(
        id = "section-content",
        htmltools::tags$ul(
            htmltools::tags$li(
                "Missing standard columns: ",
                missing
            ),
            htmltools::tags$li(
                "Non standard columns: ",
                non_st
            )
        )
    )
    styled_crit_na <- if (!is.null(critical_nas) &&
        !purrr::is_empty(checker_df)) {
        htmltools::div(
            id = "section-content",
            "NAs found in date columns that can be used for collision removal",
            htmltools::tags$ul(
                purrr::map(critical_nas, ~ htmltools::tags$li(.x))
            )
        )
    } else {
        nothing_to_rep
    }

    widget <- htmltools::tags$html(
        htmltools::tags$head(
            htmltools::tags$style(.widget_css())
        ),
        htmltools::tags$body(
            htmltools::h1("IMPORT ASSOCIATION FILE REPORT"),
            htmltools::h3(lubridate::today()),
            htmltools::h2("PROBLEMS REPORT"),
            htmltools::h3("PARSING PROBLEMS"),
            htmltools::div(
                id = "section-content",
                htmltools::div("Summary of parsing problems when reading from",
                    "file. NOTE: if the input file was in .xls or .xlsx",
                    "format, this section will be empty by default!",
                    id = "subtitle"
                )
            ),
            styled_parsing_df,
            htmltools::h3("DATE CONVERSION PROBLEMS"),
            htmltools::div(
                id = "section-content",
                htmltools::div("Summary of date conversion problems found",
                    id = "subtitle"
                )
            ),
            styled_date_df,
            htmltools::h3("COLUMNS PROBLEMS"),
            htmltools::div(
                id = "section-content",
                htmltools::div("Potential problems in columns - either",
                    "missing or non-standard columns",
                    id = "subtitle"
                )
            ),
            styled_cols_probs,
            htmltools::h3("IMPORTANT MISSING INFO"),
            htmltools::div(
                id = "section-content",
                htmltools::div("Important missing info - this info is",
                    "needed for the correct functioning of other",
                    "operations (eg: collision removal).",
                    "NOTE: this info refers ONLY to data post-",
                    "filtering (if a filter was set in the import",
                    "phase.",
                    id = "subtitle"
                )
            ),
            styled_crit_na,
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
            styled_checker_df
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
.collisions_widget <- function(input_df,
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
    additional_info <- if (is.null(add_info) || nrow(add_info) == 0) {
        htmltools::p()
    } else {
        .add_info_widget(add_info)
    }
    post_join_info <- .sc_stats_input_joined(input_joined_df, quant_cols)
    ppm <- if (length(missing) > 0) {
        input_df[-missing, ]
    } else {
        input_df
    }
    pre_process_matrix <- .generate_react_table(ppm)
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
    nothing_to_rep <- htmltools::div(
        id = "section-content",
        htmltools::div("Nothing to report", id = "simple_txt")
    )
    malf_expl <- paste0("One or more of the minimum required columns ",
                       "were missing. The minimum required columns are: ",
                       paste0(.stats_columns_min(), collapse = ", "),
                       ". Check your stats files for possible problems.")
    dupl_expl <- paste0("More than one file was found with the same prefix",
                        " for all prefixes - impossible to decide which",
                        " file to import")
    not_found_expl <- paste0("The file or the iss folder were not found")
    cols <- list(
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
        ),
        info = reactable::colDef(
            show = FALSE
        ),
        reason = reactable::colDef(
            show = FALSE
        )
    )
    row_details <- function(index, report) {
        content <- if (report$Imported[index] == TRUE) {
            nothing_to_rep
        } else {
            reason <- report$reason[index]
            reason <- if (reason == "MALFORMED") {
                paste0(c(reason, malf_expl), collapse = " - ")
            } else if (reason == "DUPLICATES") {
                paste0(c(reason, dupl_expl), collapse = " - ")
            } else if (reason == "NOT FOUND") {
                paste0(c(reason, not_found_expl), collapse = " - ")
            }
            info <- report$info[index]
            info_cont <- if (is.na(info)) {
                nothing_to_rep
            } else if (is.list(info)) {
                if (is.data.frame(info[[1]])) {
                  .generate_react_table(info[[1]])
                } else if (is.null(info[[1]])) {
                    nothing_to_rep
                } else {
                    if (all(info[[1]] == "NOT FOUND")) {
                        if (is.na((report[[.path_cols_names()$iss]])[index])) {
                            htmltools::p(style = paste(
                                "padding-left: 40px"),
                                "VISPA2 stats folder is missing")
                        } else {
                            htmltools::p(style = paste(
                                "padding-left: 40px"),
                                "No stats files found in the iss folder")
                        }
                    } else {
                        htmltools::p(style = paste(
                            "padding-left: 40px"), info)
                    }
                }
            } else {
                if (info == "NOT FOUND") {
                    if (is.na((report[[.path_cols_names()$iss]])[index])) {
                        htmltools::p(style = paste(
                            "padding-left: 40px"),
                            "VISPA2 stats folder is missing")
                    } else {
                        htmltools::p(style = paste(
                            "padding-left: 40px"),
                            "No stats files found in the iss folder")
                    }
                }
                htmltools::p(style = paste(
                    "padding-left: 40px"), info)
            }
            htmltools::div(
                htmltools::h4("NOT IMPORTED BECAUSE"),
                htmltools::p(reason, style = paste(
                    "padding-left: 40px")),
                htmltools::h4("ADDITIONAL INFO"),
                info_cont
            )
        }
        htmltools::div(
            style = paste(
                "padding-left: 40px;",
                "padding-right: 40px;",
                "padding-bottom: 20px"
            ),
            htmltools::h3("DETAILS"),
            content
        )
    }

    styled_df <- .generate_react_table(report, columns = cols,
                                       details = function(index) {
                                           row_details(index, report)
                                       })

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
        occurred and the function was unable to import the file.
        Expand the rows for more detail.",
                    id = "subtitle"
                ),
                styled_df
            )
        )
    )
    htmltools::browsable(widget)
}

.missing_iss_widget <- function(missing_iss) {
    styled_df <- .generate_react_table(missing_iss)
    widget <- htmltools::tags$html(
        htmltools::tags$head(
            htmltools::tags$style(.widget_css())
        ),
        htmltools::tags$body(
            htmltools::h2("MISSING STATS FOR SAMPLES"),
            htmltools::div(
                id = "section-content",
                htmltools::div("A summary of the samples for which no
                               corresponding stats were found",
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

#---- raw reads ----#
#' @importFrom plotly plot_ly
#' @import dplyr
#' @importFrom stringr str_detect
#' @import htmltools
#' @importFrom reactable colDef colFormat
#' @importFrom purrr map_chr
#' @importFrom tibble tibble add_row
.outliers_report_widg <- function(by_pool,
    pool_col,
    norm_test,
    key,
    flag_logic,
    outlier_thresh,
    log2_req,
    removed_nas,
    removed_zeros,
    non_proc_pools,
    flag_df) {
    li_pool_opt <- if (by_pool) {
        htmltools::tags$li("The test was run for each pool")
    } else {
        htmltools::tags$li("The test was NOT run for each pool")
    }
    li_norm_test <- if (norm_test) {
        htmltools::tags$li("Calculations performed only
                      if distribution found normal")
    } else {
        htmltools::tags$li("Normality test not performed")
    }
    li_key <- htmltools::tags$li(paste(
        "Calculations performed on columns:",
        paste0(key, collapse = ", ")
    ))
    li_thresh <- htmltools::tags$li(paste(
        "Outlier p-value threshold:",
        outlier_thresh
    ))
    li_logic <- if (length(key) > 1) {
        base_flag <- purrr::map_chr(key, function(k) {
            paste0(
                "(tdist_", k, " < ", outlier_thresh,
                " & zscore_", k, " < 0)"
            )
        })
        if (length(flag_logic) == 1) {
            flag_logic <- rep_len(x = flag_logic, length(key) - 1)
        }
        combined <- rbind(base_flag, c(flag_logic, ""))
        htmltools::tags$li(paste(
            "Key length > 1, flagging formula used: ",
            paste(combined, collapse = " ")
        ))
    } else {
        base_flag <- paste0(
                "(tdist_", key, " < ", outlier_thresh,
                " & zscore_", key, " < 0)"
            )
        htmltools::tags$li(paste(
            "Flagging formula used: ",
            paste(base_flag)
        ))
    }
    li_log2 <- if (log2_req) {
        htmltools::tags$li("Log2 transformation prior to calculations")
    } else {
        NULL
    }

    styled_removed_nas <- .generate_react_table(
        removed_nas %>%
        dplyr::select(
            dplyr::all_of(pool_col),
            .data$CompleteAmplificationID,
            dplyr::all_of(key)
        ))
    styled_removed_zeros <- if (log2_req) {
        htmltools::tags$body(
            htmltools::h3("Negative or zero values (log2 transformation)"),
            .generate_react_table(removed_zeros)
        )
    } else {
        NULL
    }

    styled_unproc_pool <- if (by_pool) {
        htmltools::tags$body(
            htmltools::h3("Unprocessed samples (per pool test)"),
            .generate_react_table(non_proc_pools, groupBy = pool_col)
        )
    } else {
        NULL
    }

    unprocessed_perc_tot <- tibble::tibble(
        abs = nrow(flag_df %>%
            dplyr::filter(.data$processed == FALSE)),
        perc = (nrow(flag_df %>%
            dplyr::filter(.data$processed == FALSE)) /
            nrow(flag_df))
    )
    styled_perc_tot <- .generate_react_table_mini(unprocessed_perc_tot,
        perc_cols = "perc"
    )
    unprocessed_perc_diff <- tibble::tibble(
        abs = c(nrow(removed_nas)),
        perc_on_unprocessed = c(nrow(removed_nas) /
            unprocessed_perc_tot$abs[1]),
        perc_on_total = c(nrow(removed_nas) / nrow(flag_df)),
        reason = c("NAs in key")
    )
    if (log2_req) {
        unprocessed_perc_diff <- unprocessed_perc_diff %>%
            tibble::add_row(
                abs = c(nrow(removed_zeros)),
                perc_on_unprocessed = c(nrow(removed_zeros) /
                    unprocessed_perc_tot$abs[1]),
                perc_on_total = c(nrow(removed_zeros) / nrow(flag_df)),
                reason = c("Values <= 0")
            )
    }
    unprocessed_perc_diff <- unprocessed_perc_diff %>%
        tibble::add_row(
            abs = c(unprocessed_perc_tot$abs[1] -
                sum(unprocessed_perc_diff$abs)),
            perc_on_unprocessed = c((unprocessed_perc_tot$abs[1] -
                sum(unprocessed_perc_diff$abs)) /
                unprocessed_perc_tot$abs[1]),
            perc_on_total = c((unprocessed_perc_tot$abs[1] -
                sum(unprocessed_perc_diff$abs)) / nrow(flag_df)),
            reason = c("Pool samples < min samples")
        )
    styled_perc_diff <- .generate_react_table_mini(unprocessed_perc_diff,
        perc_cols = c(
            "perc_on_unprocessed",
            "perc_on_total"
        )
    )

    flagged_perc <- tibble::tibble(
        abs = nrow(flag_df %>%
            dplyr::filter(.data$to_remove == TRUE)),
        perc = (nrow(flag_df %>%
            dplyr::filter(.data$to_remove == TRUE)) /
            nrow(flag_df))
    )
    styled_flagged_perc <- .generate_react_table_mini(flagged_perc,
        perc_cols = "perc"
    )

    col_def <- list(
        processed = reactable::colDef(
            align = "center",
            style = function(value) {
                color <- if (value == TRUE) {
                    "black"
                } else {
                    "#d61e1e"
                }
                list(
                    color = color, paddingLeft = "15px",
                    fontWeight = "bold",
                    fontSize = "20px"
                )
            },
            cell = function(value) {
                toupper(value)
            }
        ),
        to_remove = reactable::colDef(
            align = "center",
            style = function(value) {
                color <- if (value == TRUE) {
                    "#d61e1e"
                } else {
                    "black"
                }
                list(
                    color = color, paddingLeft = "15px",
                    fontWeight = "bold",
                    fontSize = "20px"
                )
            }, cell = function(value) {
                if (value == TRUE) paste("\u2691", value) else value
            }
        )
    )

    styled_final <- .generate_react_table(
        flag_df,
        columns = col_def, defaultSorted = list(to_remove = "desc"),
        defaultColDef = reactable::colDef(
            headerStyle = list(fontSize = "18px", paddingLeft = "15px"),
            align = "left", sortable = TRUE, resizable = TRUE,
            filterable = TRUE, style = list(paddingLeft = "15px"),
            header = function(value) gsub("_", " ", value, fixed = TRUE),
            format = reactable::colFormat(digits = 4, separators = TRUE)
        )
    )

    plot_col_style <- list(
        reactable::colDef(
            cell = function(values) {
                if (all(is.na(values))) {
                    NA
                } else {
                    plotly::plot_ly(
                        x = values[!is.na(values)], type = "histogram",
                        width = 200, height = 150
                    )
                }
            }
        )
    )

    styled_distr <- if (by_pool) {
        if (log2_req) {
            data <- flag_df %>%
                dplyr::select(
                    dplyr::all_of(pool_col),
                    dplyr::contains("log2")
                ) %>%
                dplyr::group_by(dplyr::across(pool_col)) %>%
                dplyr::summarise(dplyr::across(dplyr::contains("log2"),
                    ~ list(.x),
                    .names = "{.col}"
                ))
            log2_cols <- colnames(data)[stringr::str_detect(
                colnames(data),
                "log2"
            )]
            coldef <- rep(plot_col_style, length(log2_cols))
            names(coldef) <- log2_cols
            htmltools::div(
                id = "section-content",
                htmltools::div("Distributions for each pool"),
                .generate_react_table(data,
                    columns = coldef,
                    groupBy = pool_col
                )
            )
        } else {
            data <- flag_df %>%
                dplyr::select(
                    dplyr::all_of(pool_col),
                    dplyr::all_of(key)
                ) %>%
                dplyr::group_by(dplyr::across(pool_col)) %>%
                dplyr::summarise(dplyr::across(key,
                    ~ list(.x),
                    .names = "{.col}"
                ))
            coldef <- rep(plot_col_style, length(key))
            names(coldef) <- key
            htmltools::div(
                id = "section-content",
                htmltools::div("Distributions for each pool"),
                .generate_react_table(data,
                    columns = coldef,
                    groupBy = pool_col
                )
            )
        }
    } else {
        if (log2_req) {
            data <- flag_df %>%
                dplyr::select(dplyr::contains("log2")) %>%
                dplyr::summarise(dplyr::across(
                    .fns = ~ list(.x),
                    .names = "{.col}"
                ))
            log2_cols <- colnames(data)[stringr::str_detect(
                colnames(data),
                "log2"
            )]
            coldef <- rep(plot_col_style, length(log2_cols))
            names(coldef) <- log2_cols
            htmltools::div(
                id = "section-content",
                htmltools::div("Distributions"),
                .generate_react_table(data, columns = coldef)
            )
        } else {
            data <- flag_df %>%
                dplyr::select(dplyr::all_of(key)) %>%
                dplyr::summarise(dplyr::across(key,
                    ~ list(.x),
                    .names = "{.col}"
                ))
            coldef <- rep(plot_col_style, length(key))
            names(coldef) <- key
            htmltools::div(
                id = "section-content",
                htmltools::div("Distributions"),
                .generate_react_table(data, columns = coldef)
            )
        }
    }

    widget <- htmltools::tags$html(
        htmltools::tags$head(
            htmltools::tags$style(.widget_css())
        ),
        htmltools::tags$body(
            htmltools::h1("RAW READS OUTLIERS REPORT"),
            htmltools::h3(lubridate::today()),
            htmltools::h2("PARAMETER CHOICE AND SETTINGS"),
            htmltools::div(
                id = "section-content",
                htmltools::tags$ul(
                    li_pool_opt,
                    li_norm_test,
                    li_key,
                    li_thresh,
                    li_logic,
                    li_log2
                )
            ),
            htmltools::h2("UNPROCESSED READS"),
            htmltools::h3("NAs in key columns"),
            htmltools::div(
                id = "section-content",
                htmltools::div("These reads were not further processed
                                       because NAs were found in one or more
                                       of the key columns",
                    id = "subtitle"
                )
            ),
            styled_removed_nas,
            styled_removed_zeros,
            styled_unproc_pool,
            htmltools::h2("SUMMARY"),
            htmltools::h3("Total unprocessed reads"),
            styled_perc_tot,
            htmltools::h3("Unprocessed reads summary"),
            styled_perc_diff,
            htmltools::h3("Flagged reads"),
            styled_flagged_perc,
            htmltools::h3("Flagged metadata"),
            styled_final,
            htmltools::h3("Data distribution"),
            styled_distr
        )
    )
    htmltools::browsable(widget)
}
