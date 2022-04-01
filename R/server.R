library(shiny)
library(datamods)
library(shinyWidgets)
library(DT)
library(ggplot2)

.modal_warning_data <- function() {
    confirmSweetAlert(
        inputId = "data_override_warn",
        title = "Warning",
        text = paste(
            "Importing another data set will override the current one.",
            "Unsaved plots will be lost. Are you sure?"
        ),
        type = "warning",
        btn_labels = c("No", "Yes")
    )
}

.save_plot_alert <- function() {

}

server <- shinyServer(function(input, output, session) {
    values <- reactiveValues(
        loaded_data = NULL,
        new_data = NULL,
        plotted = NULL
    )
    imported_env <- import_globalenv_server("global_env_input")
    imported_file <- import_file_server("file_input",
        trigger_return = "button",
        return_class = "data.table"
    )
    observeEvent(input$data_override_warn, {
        if (isTRUE(input$data_override_warn)) {
            values$loaded_data <- values$new_data
            sendSweetAlert(
                title = "Success",
                text = "Data changed",
                type = "success"
            )
        }
        closeSweetAlert()
    })
    observeEvent(imported_env$data(), {
        values$new_data <- imported_env$data()
        if (!is.null(values$loaded_data)) {
            .modal_warning_data()
        } else {
            values$loaded_data <- values$new_data
            sendSweetAlert(
                title = "Success",
                text = "Data changed",
                type = "success"
            )
        }
    })
    observeEvent(imported_file$data(), {
        values$new_data <- imported_file$data()
        if (!is.null(values$loaded_data)) {
            .modal_warning_data()
        } else {
            values$loaded_data <- values$new_data
            sendSweetAlert(
                title = "Success",
                text = "Data changed",
                type = "success"
            )
        }
    })
    output$loaded_data <- renderDT({
        datatable(values$loaded_data,
            options = list(
                autoWidth = TRUE,
                class = "stripe",
                scrollX = TRUE
            ),
            filter = "top"
        )
    })
    observeEvent(input$geom, {
        function_name <- paste0("geom_", input$geom)
        values$main_geom <- rlang::as_function(function_name)
    })
    observeEvent(input$theme, {
        fun_name <- paste0("theme_", input$theme)
        values$main_theme <- rlang::as_function(fun_name)
    })
    observeEvent(input$clear_color, {
        updatePickerInput(
            session = session,
            inputId = "color",
            selected = ""
        )
    })
    observeEvent(input$clear_fill, {
        updatePickerInput(
            session = session,
            inputId = "fill",
            selected = ""
        )
    })
    observeEvent(c(
        values$loaded_data,
        input$plot_x,
        input$plot_y,
        values$main_geom,
        input$color,
        input$fill,
        input$alpha,
        values$main_theme,
        input$plot_title,
        input$x_lab,
        input$y_lab
    ), {
        req(values$loaded_data, input$plot_x, input$plot_y)
        color_aes <- if (input$color == "") {
            NULL
        } else {
            rlang::expr(.data[[input$color]])
        }
        fill_aes <- if (input$fill == "") {
            NULL
        } else {
            rlang::expr(.data[[input$fill]])
        }
        alpha_aes <- if (input$alpha == "") {
            NULL
        } else {
            rlang::expr(.data[[input$alpha]])
        }
        base_plot <- ggplot(
            data = values$loaded_data,
            mapping = aes(
                x = .data[[input$plot_x]],
                y = .data[[input$plot_y]],
                color = eval(color_aes),
                fill = eval(fill_aes),
                alpha = eval(alpha_aes)
            )
        )
        values$plotted <- base_plot +
            values$main_geom() +
            values$main_theme() +
            labs(
                title = input$plot_title, x = input$x_lab, y = input$y_lab,
                color = input$color, fill = input$fill, alpha = input$alpha
            )
    })
    observeEvent(c(input$facet_1, input$facet_2), {
        if (!is.null(values$plotted)) {
            non_null1 <- input$facet_1 != ""
            non_null2 <- input$facet_2 != ""
            if (non_null1 & non_null2) {
                values$plotted <- values$plotted +
                    facet_grid(get(input$facet_1) ~ get(input$facet_2))
            } else if (non_null1) {
                values$plotted <- values$plotted +
                    facet_wrap(~ get(input$facet_1))
            } else if (non_null2) {
                values$plotted <- values$plotted +
                    facet_wrap(~ get(input$facet_2))
            }
        }
    })
    output$plotted <- renderPlot({
        values$plotted
    })
    output$plot_tab <- renderUI({
        plot_sidebar <- sidebarPanel(
            pickerInput(
                inputId = "plot_x",
                label = "Plot on x axis",
                choices = colnames(values$loaded_data),
                options = list(
                    `live-search` = TRUE,
                    title = "Select a column"
                )
            ),
            pickerInput(
                inputId = "plot_y",
                label = "Plot on y axis",
                choices = colnames(values$loaded_data),
                options = list(
                    `live-search` = TRUE,
                    title = "Select a column"
                )
            ),
            pickerInput(
                inputId = "geom",
                label = "Geom type",
                choices = c("point", "line", "col")
            ),
            span("Color by"),
            fluidRow(
                column(
                    width = 10,
                    pickerInput(
                        inputId = "color",
                        label = "",
                        choices = colnames(values$loaded_data),
                        options = list(
                            `live-search` = TRUE,
                            title = "Select a column"
                        )
                    )
                ),
                column(
                    width = 2,
                    class = "align-self-center",
                    actionLink(
                        inputId = "clear_color", label = "",
                        icon = icon(name = "times-circle")
                    )
                )
            ),
            span("Fill by"),
            fluidRow(
                column(
                    width = 10,
                    pickerInput(
                        inputId = "fill",
                        label = "",
                        choices = colnames(values$loaded_data),
                        options = list(
                            `live-search` = TRUE,
                            title = "Select a column"
                        )
                    )
                ),
                column(
                    width = 2,
                    class = "align-self-center",
                    actionLink(
                        inputId = "clear_fill", label = "",
                        icon = icon(name = "times-circle")
                    )
                )
            ),
            span("Alpha by"),
            fluidRow(
                column(
                    width = 10,
                    pickerInput(
                        inputId = "alpha",
                        label = "",
                        choices = colnames(values$loaded_data),
                        options = list(
                            `live-search` = TRUE,
                            title = "Select a column"
                        )
                    )
                ),
                column(
                    width = 2,
                    class = "align-self-center",
                    actionLink(
                        inputId = "clear_alpha", label = "",
                        icon = icon(name = "times-circle")
                    )
                )
            ),
            pickerInput(
                inputId = "theme",
                label = "Theme",
                choices = c(
                    "classic", "grey", "bw", "linedraw", "light",
                    "dark", "minimal", "void"
                )
            ),
            span("Faceting"),
            fluidRow(
                column(
                    width = 6,
                    pickerInput(
                        inputId = "facet_1",
                        label = "",
                        choices = colnames(values$loaded_data),
                        options = list(
                            `live-search` = TRUE,
                            title = "Select a column"
                        )
                    )
                ),
                column(
                    width = 6,
                    pickerInput(
                        inputId = "facet_2",
                        label = "",
                        choices = colnames(values$loaded_data),
                        options = list(
                            `live-search` = TRUE,
                            title = "Select a column"
                        )
                    )
                )
            ),
            textInput("plot_title", label = "Plot title"),
            textInput("x_lab", label = "X axis title"),
            textInput("y_lab", label = "Y axis title"),
            fluidRow(
                column(
                    width = 12,
                    dropdownButton(
                        pickerInput("device", "Device",
                            choices = c(
                                "pdf", "jpeg", "tiff", "png",
                                "bmp", "svg", "eps", "ps", "tex"
                            )
                        ),
                        numericInputIcon("plot_file_width",
                            label = "Width",
                            value = 8,
                            min = 1, icon = icon("arrows-alt-h")
                        ),
                        numericInputIcon("plot_file_height",
                            label = "Height",
                            value = 8,
                            min = 1, icon = icon("arrows-alt-v")
                        ),
                        pickerInput("units", "Units",
                            choices = c("", "in", "cm", "mm", "px")
                        ),
                        numericInputIcon("plot_file_res",
                            label = "Resolution (dpi)",
                            value = 300,
                            min = 72, icon = icon("desktop")
                        ),
                        downloadButton("save_plot",
                            label = "Save",
                            icon = NULL
                        ),
                        circle = FALSE,
                        status = "default",
                        icon = icon("download")
                    ),
                    align = "center"
                )
            )
        )
        plot_main <- mainPanel(
            plotOutput("plotted", height = "100%")
        )
        sidebarLayout(
            plot_sidebar,
            plot_main,
            position = "right"
        )
    })
    output$save_plot <- downloadHandler(
        filename = function() {
            paste0(lubridate::today(), "_output-plot.", input$device)
        },
        content = function(file) {
            if (input$units == "") {
                ggsave(
                    plot = values$plotted, filename = file,
                    width = input$plot_file_width,
                    height = input$plot_file_height,
                    dpi = input$plot_file_res,
                    device = input$device
                )
            } else {
                ggsave(
                    plot = values$plotted, filename = file,
                    width = input$plot_file_width,
                    height = input$plot_file_height,
                    dpi = input$plot_file_res,
                    device = input$device,
                    units = input$units
                )
            }
        }
    )
})
