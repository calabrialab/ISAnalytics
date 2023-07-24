#' @import shiny
#' @import shinyWidgets
#' @import datamods
#' @import bslib
ui <- shinyUI(fluidPage(
    useSweetAlert(),
    ### --- Top level navbar
    navbarPage(
        "NGSdataExplorer",
        theme = bs_theme(
            bootswatch = "cosmo",
            primary = "#23687A",
            version = 4,
            "navbar-bg" = "#23687A"
        ),
        collapsible = TRUE,
        ### --- First nav element
        tabPanel("Home", {
            ### --- Lateral bar + tabs
            navlistPanel(
                "Load data",
                tabPanel(
                    "From R environment",
                    fluidRow(
                        column(
                            width = 12,
                            import_globalenv_ui("global_env_input")
                        )
                    )
                ),
                tabPanel(
                    "From file",
                    fluidRow(
                        column(
                            width = 12,
                            import_file_ui("file_input")
                        )
                    )
                )
            )
        }),
        ### --- Second nav element
        tabPanel(
            "Explore",
            navset_pill(
                id = "explore_pills",
                nav_panel(
                    title = "Explore loaded data",
                    DT::DTOutput("loaded_data"),
                    class = "p-3 border rounded"
                ),
                nav_panel(
                    title = "Plotting",
                    class = "p-3 border rounded",
                    uiOutput("plot_tab")
                )
            )
        )
    )
))
