library(shiny)
library(shinyWidgets)
library(datamods)
library(bslib)
library(DT)

# Define UI for application that draws a histogram
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
            navs_pill(
                id = "explore_pills",
                nav(
                    title = "Explore loaded data",
                    DTOutput("loaded_data"),
                    class = "p-3 border rounded"
                ),
                nav(
                    title = "Plotting",
                    class = "p-3 border rounded",
                    uiOutput("plot_tab")
                )
            )
        )
    )
))
