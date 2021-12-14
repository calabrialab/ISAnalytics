#------------------------------------------------------------------------------#
# Plotting functions
#------------------------------------------------------------------------------#

#' Trace volcano plot for computed CIS data.
#'
#' \lifecycle{stable}
#' Traces a volcano plot for IS frequency and CIS results.
#'
#' @details
#' ## Input data frame
#' Users can supply as `x` either a simple integration matrix or a
#' data frame resulting from the call to \link{CIS_grubbs}.
#' In the first case an internal call to
#' the function `CIS_grubbs()` is performed.
#'
#' ## Oncogene and tumor suppressor genes files
#' These files are included in the package for user convenience and are
#' simply UniProt files with gene annotations for human and mouse.
#' For more details on how this files were generated use the help
#' `?tumor_suppressors`, `?proto_oncogenes`
#'
#' ## Known oncogenes
#' The default values are included in this package and
#' it can be accessed by doing:
#'
#' ```{r}
#' head(known_clinical_oncogenes())
#'
#' ```
#' If the user wants to change this parameter the input data frame must
#' preserve the column structure. The same goes for the `suspicious_genes`
#' parameter (DOIReference column is optional):
#'
#' ```{r}
#' head(clinical_relevant_suspicious_genes())
#'
#' ```
#' @family Plotting functions
#'
#' @param x Either a simple integration matrix or a data frame resulting
#' from the call to \link{CIS_grubbs} with `add_standard_padjust = TRUE`
#' @param onco_db_file Uniprot file for proto-oncogenes (see details).
#' If different from default, should be supplied as a path to a file.
#' @param tumor_suppressors_db_file Uniprot file for tumor-suppressor genes.
#' If different from default, should be supplied as a path to a file.
#' @param species One between `"human"`, `"mouse"` and `"all"`
#' @param known_onco Data frame with known oncogenes. See details.
#' @param suspicious_genes Data frame with clinical relevant suspicious
#' genes. See details.
#' @param significance_threshold The significance threshold
#' @param annotation_threshold_ontots Value above which genes are annotated
#' with colorful labels
#' @param highlight_genes Either `NULL` or a character vector of genes to be
#' highlighted in the plot even if they're not above the threshold
#' @param title_prefix A string or character vector to be displayed
#' in the title - usually the
#' project name and other characterizing info. If a vector is supplied,
#' it is concatenated in a single string via `paste()`
#' @param return_df Return the data frame used to generate the plot?
#' This can be useful if the user wants to manually modify the plot with
#' ggplot2. If TRUE the function returns a list containing both the plot
#' and the data frame.
#'
#' @importFrom ggplot2 ggplot aes_ geom_point geom_hline scale_y_continuous
#' @importFrom ggplot2 scale_x_continuous unit theme element_text labs
#' @importFrom ggrepel geom_label_repel
#' @importFrom dplyr left_join mutate filter
#' @importFrom rlang .data inform
#' @importFrom purrr is_empty
#' @importFrom stringr str_to_lower
#'
#' @return A plot or a list containing a plot and a data frame
#' @export
#'
#' @examples
#' data("integration_matrices", package = "ISAnalytics")
#' cis_plot <- CIS_volcano_plot(integration_matrices,
#'     title_prefix = "PJ01"
#' )
#' cis_plot
CIS_volcano_plot <- function(x,
    onco_db_file = "proto_oncogenes",
    tumor_suppressors_db_file = "tumor_suppressors",
    species = "human",
    known_onco = known_clinical_oncogenes(),
    suspicious_genes =
        clinical_relevant_suspicious_genes(),
    significance_threshold = 0.05,
    annotation_threshold_ontots = 0.1,
    highlight_genes = NULL,
    title_prefix = NULL,
    return_df = FALSE) {
    ## Check params
    stopifnot(is.data.frame(x))
    stopifnot(is.character(onco_db_file))
    onco_db_file <- onco_db_file[1]
    stopifnot(is.character(tumor_suppressors_db_file))
    tumor_suppressors_db_file <- tumor_suppressors_db_file[1]
    stopifnot(is.character(species))
    stopifnot(is.data.frame(known_onco))
    stopifnot(is.data.frame(suspicious_genes))
    stopifnot(is.numeric(significance_threshold) |
        is.integer(significance_threshold))
    significance_threshold <- significance_threshold[1]
    stopifnot(is.numeric(annotation_threshold_ontots) |
        is.integer(annotation_threshold_ontots))
    annotation_threshold_ontots <- annotation_threshold_ontots[1]
    stopifnot(is.null(title_prefix) || (is.character(title_prefix)))
    stopifnot(is.null(highlight_genes) || is.character(highlight_genes))
    stopifnot(is.logical(return_df))
    if (is.null(title_prefix)) {
        title_prefix <- ""
    } else {
        title_prefix <- paste(title_prefix, collapse = " ")
    }
    ## Load onco and ts
    oncots_to_use <- .load_onco_ts_genes(
        onco_db_file,
        tumor_suppressors_db_file,
        species
    )
    ## Check if CIS function was already called
    min_cis_col <- c(
        "tdist_bonferroni_default", "tdist_fdr",
        "neg_zscore_minus_log2_int_freq_tolerance"
    )
    cis_grubbs_df <- if (!all(min_cis_col %in% colnames(x))) {
        if (getOption("ISAnalytics.verbose") == TRUE) {
            rlang::inform("Calculating CIS_grubbs for x...")
        }
        CIS_grubbs(x)
    } else {
        x
    }
    ## Join all dfs by gene
    cis_grubbs_df <- cis_grubbs_df %>%
        dplyr::left_join(oncots_to_use, by = "GeneName") %>%
        dplyr::left_join(known_onco, by = "GeneName") %>%
        dplyr::left_join(suspicious_genes, by = "GeneName")
    ## Add info to CIS
    cis_grubbs_df <- cis_grubbs_df %>%
        dplyr::mutate(minus_log_p = -log(.data$tdist_bonferroni_default,
            base = 10
        ))
    cis_grubbs_df <- cis_grubbs_df %>%
        dplyr::mutate(
            minus_log_p_fdr = -log(.data$tdist_fdr, base = 10),
            positive_outlier_and_significant = ifelse(
                test = !is.na(.data$tdist_fdr) &
                    .data$tdist_fdr < significance_threshold,
                yes = TRUE,
                no = FALSE
            )
        )
    cis_grubbs_df <- cis_grubbs_df %>%
        dplyr::mutate(
            KnownGeneClass = ifelse(
                is.na(.data$Onco1_TS2),
                yes = "Other",
                no = ifelse(.data$Onco1_TS2 == 1,
                    yes = "OncoGene",
                    no = "TumSuppressor"
                )
            ),
            CriticalForInsMut = ifelse(!is.na(.data$KnownClonalExpansion),
                yes = TRUE, no = FALSE
            )
        )
    significance_threshold_minus_log_p <- -log(significance_threshold,
        base = 10
    )
    annotation_threshold_ontots_log <- -log(annotation_threshold_ontots,
        base = 10
    )
    ## Trace plot
    plot_cis_fdr_slice <- ggplot2::ggplot(
        data = cis_grubbs_df,
        ggplot2::aes_(
            y = ~minus_log_p_fdr,
            x = ~neg_zscore_minus_log2_int_freq_tolerance,
            color = ~KnownGeneClass,
            fill = ~KnownGeneClass
        ),
        na.rm = TRUE, se = TRUE
    ) +
        ggplot2::geom_point(alpha = .5, size = 3) +
        ggplot2::geom_hline(
            yintercept = significance_threshold_minus_log_p,
            color = "black", size = 1, show.legend = TRUE, linetype = "dashed"
        ) +
        ggplot2::scale_y_continuous(limits = c(0, max(c(
            (significance_threshold_minus_log_p + 0.5),
            max(cis_grubbs_df$minus_log_p_fdr, na.rm = TRUE)
        ), na.rm = TRUE))) +
        ggplot2::scale_x_continuous(breaks = seq(-4, 4, 2)) +
        ggrepel::geom_label_repel(
            data = dplyr::filter(
                cis_grubbs_df,
                .data$tdist_fdr < significance_threshold
            ),
            ggplot2::aes_(label = ~GeneName),
            box.padding = ggplot2::unit(0.35, "lines"),
            point.padding = ggplot2::unit(0.3, "lines"),
            color = "white",
            segment.color = "black",
            max.overlaps = Inf
        ) +
        ggplot2::theme(
            strip.text.y = ggplot2::element_text(
                size = 16,
                colour = "blue",
                angle = 270
            ),
            strip.text.x = ggplot2::element_text(
                size = 16,
                colour = "blue",
                angle = 0
            )
        ) +
        ggplot2::theme(strip.text = ggplot2::element_text(
            face = "bold",
            size = 16
        )) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(size = 16),
            axis.text.y = ggplot2::element_text(size = 16),
            axis.title = ggplot2::element_text(size = 16),
            plot.title = ggplot2::element_text(size = 20)
        ) +
        ggplot2::labs(
            title = paste(
                title_prefix,
                "- Volcano plot of IS gene frequency and",
                "CIS results"
            ),
            y = "P-value Grubbs test (-log10(p))",
            x = "Integration frequency (log2)",
            size = "Avg Transcr. Len",
            color = "Onco TumSupp Genes",
            subtitle = paste0(
                "Significance threshold for annotation",
                " labeling: P-value < ", significance_threshold,
                " (FDR adjusted; ",
                "-log = ", (round(-log(significance_threshold, base = 10), 3)),
                ").\nOnco/TS genes source: UniProt (other genes ",
                "labeled as 'Other'). \nAnnotated if P-value > ",
                round(annotation_threshold_ontots_log, 3), " or in highlighted",
                " genes"
            )
        )
    if (!is.null(highlight_genes) && !purrr::is_empty(highlight_genes)) {
        ## Look for the genes (case insensitive)
        to_highlight <- cis_grubbs_df %>%
            dplyr::filter(
                stringr::str_to_lower(.data$GeneName) %in%
                    stringr::str_to_lower(highlight_genes),
                .data$tdist_fdr >= significance_threshold
            )
        plot_cis_fdr_slice <- plot_cis_fdr_slice +
            ggrepel::geom_label_repel(
                data = to_highlight,
                ggplot2::aes_(label = ~GeneName),
                box.padding = ggplot2::unit(0.35, "lines"),
                point.padding = ggplot2::unit(0.3, "lines"),
                color = "white",
                segment.color = "black",
                max.overlaps = Inf
            )
    }
    if (return_df) {
        return(list(plot = plot_cis_fdr_slice, df = cis_grubbs_df))
    } else {
        return(plot_cis_fdr_slice)
    }
}

#' Known clinical oncogenes (for mouse and human).
#'
#' @return A data frame
#' @importFrom tibble tibble
#'
#' @family Plotting function helpers
#' @export
#'
#' @examples
#' known_clinical_oncogenes()
known_clinical_oncogenes <- function() {
    tibble::tibble(
        GeneName = c("MECOM", "CCND2", "TAL1", "LMO2", "HMGA2"),
        KnownClonalExpansion = TRUE
    )
}

#' Clinical relevant suspicious genes (for mouse and human).
#'
#' @return A data frame
#' @importFrom tibble tibble
#'
#' @family Plotting function helpers
#' @export
#'
#' @examples
#' clinical_relevant_suspicious_genes()
clinical_relevant_suspicious_genes <- function() {
    tibble::tibble(
        GeneName = c(
            "DNMT3A", "TET2", "ASXL1",
            "JAK2", "CBL", "TP53"
        ),
        ClinicalRelevance = TRUE,
        DOIReference =
            "https://doi.org/10.1182/blood-2018-01-829937"
    )
}

#' Plot of the estimated HSC population size for each patient.
#'
#' @param estimates The estimates data frame, obtained via
#' \code{\link{HSC_population_size_estimate}}
#' @param project_name The project name, will be included in the plot title
#' @param timepoints Which time points to plot? One between "All",
#' "Stable" and "Consecutive"
#' @param models Name of the models to plot (as they appear in the column
#' of the estimates)
#'
#' @family Plotting functions
#'
#' @import ggplot2
#'
#' @return A plot
#' @export
#'
#' @examples
#' data("integration_matrices", package = "ISAnalytics")
#' data("association_file", package = "ISAnalytics")
#' aggreg <- aggregate_values_by_key(
#'     x = integration_matrices,
#'     association_file = association_file,
#'     value_cols = c("seqCount", "fragmentEstimate")
#' )
#' aggreg_meta <- aggregate_metadata(
#'     association_file = association_file
#' )
#' estimate <- HSC_population_size_estimate(
#'     x = aggreg,
#'     metadata = aggreg_meta,
#'     stable_timepoints = c(90, 180, 360),
#'     cell_type = "Other"
#' )
#' p <- HSC_population_plot(estimate, "PJ01")
#' p
HSC_population_plot <- function(estimates,
    project_name,
    timepoints = "Consecutive",
    models = "Mth Chao (LB)") {
    if (is.null(estimates)) {
        return(NULL)
    }
    ## Pre-filter
    df <- estimates %>%
        dplyr::filter(
            .data$Timepoints %in% timepoints,
            .data$Model %in% models
        )
    p <- ggplot2::ggplot(
        data = df,
        ggplot2::aes_(
            y = ~PopSize,
            x = ~TimePoint_to,
            color = ~SubjectID
        ),
        na.rm = TRUE, se = TRUE
    ) +
        ggplot2::geom_point(alpha = .5) +
        ggplot2::geom_line(size = 2, alpha = .7) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(size = 14),
            axis.text.y = ggplot2::element_text(size = 14),
            axis.title = ggplot2::element_text(size = 16),
            plot.title = ggplot2::element_text(size = 20),
            strip.text.x = ggplot2::element_text(
                size = 14,
                colour = "darkblue",
                angle = 0,
                face = "bold"
            ),
            strip.text.y = ggplot2::element_text(
                size = 14,
                colour = "darkred",
                angle = 270,
                face = "bold"
            )
        ) +
        ggplot2::labs(
            title = paste(project_name, "- HSC Population size"),
            x = "Time Point (months after GT)",
            y = "HSC size (Chao model with bias correction)",
            colour = "Patient",
            subtitle = "IS from Myeloid PB cells as surrogate of HSC."
        )
    p
}

## --- Alluvial plots --- ##

#' Alluvial plots for IS distribution in time.
#'
#' \lifecycle{experimental}
#' Alluvial plots allow the visualization of integration sites distribution
#' in different points in time in the same group.
#' This functionality requires the suggested package
#' [ggalluvial](https://corybrunson.github.io/ggalluvial/).
#'
#' @details
#' ### Input data frame
#' The input data frame must contain all the columns specified in the
#' arguments `group`, `plot_x`, `plot_y` and `alluvia`. The standard
#' input for this function is the data frame obtained via the
#' \link{compute_abundance} function.
#'
#' ### Plotting threshold on y
#' The plotting threshold on the quantification on the y axis has the
#' function to highlight only relevant information on the plot and reduce
#' computation time. The default value is 1, that acts on the default column
#' plotted on the y axis which holds a percentage value. This translates
#' in natural language roughly as "highlight with colors only those
#' integrations (alluvia) that at least in 1 point in time have an
#' abundance value >= 1 %". The remaining integrations will be plotted
#' as transparent in the strata.
#'
#' @param x A data frame. See details.
#' @param group Character vector containing the column names that identify
#' unique groups.
#' @param plot_x Column name to plot on the x axis
#' @param plot_y Column name to plot on the y axis
#' @param alluvia Character vector of column names that uniquely identify
#' alluvia
#' @param alluvia_plot_y_threshold Numeric value. Everything below this
#' threshold on y will be plotted in grey and aggregated. See details.
#' @param top_abundant_tbl Logical. Produce the summary top abundant tables
#' via \link{top_abund_tableGrob}?
#' @param ... Additional arguments to pass on to \link{top_abund_tableGrob}
#'
#' @family Plotting functions
#' @importFrom rlang abort eval_tidy call2 inform .data fn_fmls_names dots_list
#' @importFrom dplyr group_by across group_keys everything pull group_split
#' @importFrom BiocParallel SnowParam MulticoreParam bplapply bpstop
#' @importFrom purrr set_names
#' @importFrom tidyr unite
#'
#' @return For each group a list with the associated plot and optionally
#' the summary tableGrob
#' @export
#'
#' @examples
#' data("integration_matrices", package = "ISAnalytics")
#' data("association_file", package = "ISAnalytics")
#' aggreg <- aggregate_values_by_key(
#'     x = integration_matrices,
#'     association_file = association_file,
#'     value_cols = c("seqCount", "fragmentEstimate")
#' )
#' abund <- compute_abundance(x = aggreg)
#' alluvial_plots <- integration_alluvial_plot(abund,
#'     alluvia_plot_y_threshold = 0.5
#' )
#' ex_plot <- alluvial_plots[[1]]$plot +
#'     ggplot2::labs(
#'         title = "IS distribution over time",
#'         subtitle = "Patient 1, MNC BM",
#'         y = "Abundance (%)",
#'         x = "Time point (days after GT)"
#'     )
#' print(ex_plot)
integration_alluvial_plot <- function(x,
    group = c(
        "SubjectID",
        "CellMarker",
        "Tissue"
    ),
    plot_x = "TimePoint",
    plot_y = "fragmentEstimate_sum_PercAbundance",
    alluvia = mandatory_IS_vars(),
    alluvia_plot_y_threshold = 1,
    top_abundant_tbl = TRUE,
    ...) {
    if (!requireNamespace("ggalluvial", quietly = TRUE)) {
        rlang::abort(.missing_pkg_error("ggalluvial"))
    }
    stopifnot(is.logical(top_abundant_tbl))
    stopifnot(is.data.frame(x))
    stopifnot(is.character(group))
    stopifnot(is.character(plot_x))
    stopifnot(is.character(plot_y))
    stopifnot(is.numeric(alluvia_plot_y_threshold))
    stopifnot(is.character(alluvia))
    plot_x <- plot_x[1]
    plot_y <- plot_y[1]
    if (any(!c(group, plot_x, plot_y, alluvia) %in% colnames(x))) {
        rlang::abort(.missing_user_cols_error(
            c(group, plot_x, plot_y, alluvia)[!c(group, plot_x, plot_y, alluvia)
            %in% colnames(x)]
        ))
    }
    groups_to_plot <- x %>%
        dplyr::group_by(dplyr::across({{ group }}))
    group_names <- groups_to_plot %>%
        dplyr::group_keys() %>%
        tidyr::unite(col = "id", dplyr::everything()) %>%
        dplyr::pull(.data$id)
    groups_to_plot <- groups_to_plot %>%
        dplyr::group_split() %>%
        purrr::set_names(group_names)

    # Compute plots in parallel
    p <- BiocParallel::MulticoreParam(
        stop.on.error = FALSE, progressbar = TRUE,
        tasks = length(groups_to_plot), exportglobals = FALSE
    )


    FUN <- function(group_df,
    plot_x,
    plot_y,
    alluvia,
    alluvia_plot_y_threshold,
    top_abundant_tbl,
    other_params) {
        cleaned <- .clean_data(
            group_df, plot_x, plot_y,
            alluvia, alluvia_plot_y_threshold
        )
        alluv_plot <- .alluvial_plot(cleaned, plot_x, plot_y)
        if (top_abundant_tbl == TRUE) {
            summary_tbls <- NULL
            withCallingHandlers(
                {
                    withRestarts(
                        {
                            summary_tbls <- rlang::eval_tidy(
                                rlang::call2("top_abund_tableGrob",
                                    group_df,
                                    id_cols = alluvia,
                                    quant_col = plot_y,
                                    by = plot_x,
                                    alluvial_plot = alluv_plot,
                                    !!!other_params
                                )
                            )
                        },
                        missing = function() {
                            rlang::inform(paste(
                                "Failed to produce top",
                                "abundance tables, skipping"
                            ))
                        }
                    )
                },
                error = function(cnd) {
                    rlang::inform(conditionMessage(cnd))
                    invokeRestart("missing")
                }
            )
            return(list(plot = alluv_plot, tables = summary_tbls))
        }
        return(alluv_plot)
    }

    dot_args <- rlang::dots_list(..., .named = TRUE, .homonyms = "first")
    optional_params_names <- rlang::fn_fmls_names(top_abund_tableGrob)
    optional_params_names <- optional_params_names[!optional_params_names %in%
        c(
            "df", "id_cols",
            "quant_col", "by",
            "alluvial_plot"
        )]
    dot_args <- dot_args[names(dot_args) %in% optional_params_names]
    results <- BiocParallel::bplapply(
        X = groups_to_plot,
        FUN = FUN,
        plot_x = plot_x,
        plot_y = plot_y,
        alluvia = alluvia,
        alluvia_plot_y_threshold = alluvia_plot_y_threshold,
        top_abundant_tbl = top_abundant_tbl,
        other_params = dot_args,
        BPPARAM = p
    )
    BiocParallel::bpstop(p)
    return(results)
}

# Internal, used in integration_alluvial_plot to obtain the data frame
# with data to plot
#' @importFrom tidyr unite
#' @importFrom data.table setDT .SD
#' @importFrom dplyr group_by across summarise n filter distinct pull rename
#' @importFrom rlang .data
.clean_data <- function(df,
    plot_x,
    plot_y,
    alluvia,
    alluvia_plot_y_threshold) {
    tbl <- if (length(alluvia) > 1) {
        df %>%
            tidyr::unite(col = "alluvia_id", {{ alluvia }}) %>%
            data.table::setDT()
    } else {
        df %>%
            dplyr::rename(alluvia_id = alluvia) %>%
            data.table::setDT()
    }
    ## Getting counts by x
    counts_by_x <- tbl %>%
        dplyr::group_by(dplyr::across({{ plot_x }})) %>%
        dplyr::summarise(count = dplyr::n(), .groups = "drop")
    # Filtering alluvia to plot
    alluvia_to_plot <- tbl %>%
        dplyr::filter(.data[[plot_y]] >= alluvia_plot_y_threshold[1]) %>%
        dplyr::distinct(.data[["alluvia_id"]]) %>%
        dplyr::pull(.data[["alluvia_id"]])
    # Modify ids - identifiers that are below threshold are converted to
    # NA and the quantities in y are aggregated
    tbl[!alluvia_id %in% alluvia_to_plot, alluvia_id := NA_character_]
    tbl <- tbl[, setNames(list(sum(.SD)), plot_y),
        by = c(plot_x, "alluvia_id"), .SDcols = plot_y
    ]
    # Add counts
    tbl <- merge(tbl, counts_by_x)
    tbl
}

# Internal, used in integration_alluvial_plot to obtain the alluvial plots
# for a single group. NOTE: tbl must contain the column "alluvia_id" and
# "counts"
#' @importFrom ggplot2 ggplot aes_ geom_text scale_fill_viridis_d sym
#' @importFrom dplyr group_by summarise pull across all_of
#' @importFrom rlang .data
.alluvial_plot <- function(tbl, plot_x, plot_y) {
    sums_y <- tbl %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(plot_x))) %>%
        dplyr::summarise(sums = sum(.data[[plot_y]])) %>%
        dplyr::pull(.data$sums)
    max_y <- max(sums_y)
    labels <- tbl %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(plot_x))) %>%
        dplyr::summarise(count = .data$count[1], .groups = "drop")
    alluv <- ggplot2::ggplot(
        tbl,
        ggplot2::aes_(
            x = ggplot2::sym(plot_x),
            y = ggplot2::sym(plot_y),
            alluvium = ~alluvia_id
        )
    ) +
        ggalluvial::geom_stratum(ggplot2::aes_(stratum = ~alluvia_id),
            na.rm = FALSE, decreasing = FALSE,
            fill = NA
        ) +
        ggalluvial::geom_alluvium(ggplot2::aes_(fill = ~alluvia_id),
            na.rm = FALSE,
            decreasing = FALSE,
            alpha = .75,
            aes.bind = "alluvia"
        ) +
        ggplot2::scale_fill_viridis_d() +
        ggplot2::theme(legend.position = "none")
    alluv <- alluv +
        ggplot2::geom_text(
            data = labels,
            ggplot2::aes_(
                x = ggplot2::sym(plot_x),
                y = max_y + 5,
                label = ~count
            ), inherit.aes = FALSE
        )
    return(alluv)
}

#' Summary top abundant tableGrobs for plots.
#'
#' Produce summary tableGrobs as R graphics. For this functionality
#' the suggested package
#' [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html)
#' is required. To visualize the resulting object:
#' ```
#' gridExtra::grid.arrange(tableGrob)
#' ```
#'
#' @param df A data frame
#' @param id_cols Character vector of id column names. To plot after alluvial,
#' these columns must be the same as the `alluvia` argument of
#' \link{integration_alluvial_plot}.
#' @param quant_col Column name holding the quantification value.
#' To plot after alluvial,
#' these columns must be the same as the `plot_y` argument of
#' \link{integration_alluvial_plot}.
#' @param by The column name to subdivide tables for. The function
#' will produce one table for each distinct value in `by`.
#' To plot after alluvial,
#' these columns must be the same as the `plot_x` argument of
#' \link{integration_alluvial_plot}.
#' @param alluvial_plot Either NULL or an alluvial plot for color mapping
#' between values of y.
#' @param top_n Integer. How many rows should the table contain at most?
#' @param tbl_cols Table columns to show in the final output besides
#' `quant_col`.
#' @param include_id_cols Logical. Include `id_cols` in the output?
#' @param digits Integer. Digits to show for the quantification column
#' @param perc_symbol Logical. Show percentage symbol in the quantification
#' column?
#'
#' @family Plotting functions
#' @return A tableGrob object
#' @export
#'
#' @importFrom rlang abort .data
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot_build
#' @importFrom dplyr rename mutate select all_of across group_modify arrange
#' @importFrom dplyr desc slice_head ungroup distinct left_join filter
#' @importFrom dplyr rename_with pull starts_with
#' @importFrom tidyr unite
#' @importFrom purrr map set_names map2 reduce
#' @importFrom stringr str_detect
#'
#' @examples
#' data("integration_matrices", package = "ISAnalytics")
#' data("association_file", package = "ISAnalytics")
#' aggreg <- aggregate_values_by_key(
#'     x = integration_matrices,
#'     association_file = association_file,
#'     value_cols = c("seqCount", "fragmentEstimate")
#' )
#' abund <- compute_abundance(x = aggreg)
#' grob <- top_abund_tableGrob(abund)
#' gridExtra::grid.arrange(grob)
top_abund_tableGrob <- function(df,
    id_cols = mandatory_IS_vars(),
    quant_col = "fragmentEstimate_sum_PercAbundance",
    by = "TimePoint",
    alluvial_plot = NULL,
    top_n = 10,
    tbl_cols = "GeneName",
    include_id_cols = FALSE,
    digits = 2,
    perc_symbol = TRUE) {
    if (!requireNamespace("gridExtra", quietly = TRUE)) {
        rlang::abort(.missing_pkg_error("gridExtra"))
    }
    stopifnot(is.numeric(top_n))
    stopifnot(is.character(tbl_cols))
    stopifnot(all(tbl_cols %in% colnames(df)))
    stopifnot(is.null(alluvial_plot) |
        any(class(alluvial_plot) %in% c("gg", "ggplot")))
    color_map <- if (!is.null(alluvial_plot)) {
        g <- ggplot2::ggplot_build(alluvial_plot)
        tibble::tibble(g$data[[2]]["alluvium"], g$data[[2]]["fill"]) %>%
            dplyr::rename(alluvia_id = "alluvium")
    } else {
        NULL
    }
    top <- df %>%
        dplyr::mutate(dplyr::across({{ by }}, ~ ifelse(is.na(.x),
            yes = "NA",
            no = .x
        ))) %>%
        tidyr::unite({{ id_cols }}, col = "alluvia_id") %>%
        dplyr::select(dplyr::all_of(c(
            "alluvia_id", tbl_cols,
            quant_col, by
        ))) %>%
        dplyr::group_by(dplyr::across({{ by }})) %>%
        dplyr::group_modify(~ {
            .x %>%
                dplyr::arrange(dplyr::desc(.data[[quant_col[1]]])) %>%
                dplyr::slice_head(n = top_n)
        }) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(font_col = "black") %>%
        dplyr::rename(Ab = quant_col)
    if (perc_symbol == TRUE) {
        top <- top %>%
            dplyr::mutate(Ab = paste(format(
                round(.data$Ab, digits = digits),
                nsmall = digits
            ), "%"))
    } else {
        top <- top %>%
            dplyr::mutate(Ab = format(round(.data$Ab,
                digits = digits
            ),
            nsmall = digits
            ))
    }
    if (!is.null(color_map)) {
        top <- top %>%
            dplyr::left_join(color_map, by = c("alluvia_id")) %>%
            dplyr::distinct()
    }
    if (include_id_cols == FALSE) {
        top <- top %>%
            dplyr::select(-.data$alluvia_id)
    }
    distinct_x <- unique(top[[by]])
    tops_by_x <- purrr::map(distinct_x, function(x) {
        tmp <- top %>%
            dplyr::filter(
                dplyr::across(
                    dplyr::all_of(by), ~ .x == x
                )
            ) %>%
            dplyr::select(-dplyr::all_of(by))
        if (nrow(tmp) < top_n) {
            index_1 <- nrow(tmp) + 1
            tmp[seq(from = index_1, to = top_n, by = 1), ] <-
                as.list(c(
                    rep_len(
                        NA,
                        length(colnames(tmp)[colnames(tmp) != "font_col"])
                    ),
                    "transparent"
                ))
        }
        tmp <- tmp %>%
            dplyr::rename_with(.fn = ~ paste0(.x, " ", x))
    }) %>% purrr::set_names(distinct_x)

    obtain_grobs <- function(df, x) {
        fill_var <- colnames(df)[stringr::str_detect(colnames(df), "fill")]
        col_var <- colnames(df)[stringr::str_detect(colnames(df), "font_col")]
        fill_vec <- if (!length(fill_var) == 0) {
            df %>% dplyr::pull(fill_var)
        } else {
            NULL
        }
        col_vec <- df %>% dplyr::pull(col_var)
        df_minus <- df %>%
            dplyr::select(
                -dplyr::starts_with("fill"),
                -dplyr::starts_with("font_col")
            )
        theme <- gridExtra::ttheme_minimal(
            base_size = 8,
            core = list(
                bg_params = list(fill = fill_vec, alpha = 0.75),
                fg_params = list(col = col_vec, parse = FALSE)
            )
        )
        grob_x <- gridExtra::tableGrob(df_minus,
            rows = NULL, theme = theme
        )
        grob_x
    }
    single_grobs <- purrr::map2(tops_by_x, names(tops_by_x), obtain_grobs)
    all_grobs <- purrr::reduce(single_grobs, function(d1, d2) {
        gridExtra::gtable_combine(d1, d2, along = 1)
    })
    all_grobs
}

#' Plot IS sharing heatmaps.
#'
#' \lifecycle{experimental}
#' Displays the IS sharing calculated via \link{is_sharing} as heatmaps.
#'
#' @param sharing_df The data frame containing the IS sharing data
#' @param show_on_x Name of the column to plot on the x axis
#' @param show_on_y Name of the column to plot on the y axis
#' @param absolute_sharing_col Name of the column that contains the absolute
#' values of IS sharing
#' @param title_annot Additional text to display in the title
#' @param plot_relative_sharing Logical. Compute heatmaps also for relative
#' sharing?
#' @param rel_sharing_col Names of the columns to consider as relative sharing.
#' The function is going to plot one heatmap per column in this argument.
#' @param show_perc_symbol_rel Logical. Only relevant if `plot_relative_sharing`
#' is set to TRUE, should the percentage symbol be displayed in relative
#' heatmaps?
#' @param interactive Logical. Requires the package
#' [plotly](https://plotly.com/r/getting-started/) is required for this
#' functionality. Returns the heatmaps as interactive HTML widgets.
#'
#' @family Plotting functions
#' @return A list of plots or widgets
#' @seealso \link{is_sharing}
#' @export
#'
#' @importFrom rlang abort inform .data
#' @importFrom ggplot2 ggplot aes_ geom_raster scale_fill_gradientn geom_text
#' @importFrom ggplot2 scale_alpha_continuous aes theme element_text
#' @importFrom ggplot2 element_blank rel labs
#' @importFrom dplyr mutate across all_of
#' @importFrom purrr map set_names
#'
#' @examples
#' data("integration_matrices", package = "ISAnalytics")
#' data("association_file", package = "ISAnalytics")
#' aggreg <- aggregate_values_by_key(
#'     x = integration_matrices,
#'     association_file = association_file,
#'     value_cols = c("seqCount", "fragmentEstimate")
#' )
#' sharing <- is_sharing(aggreg,
#'     minimal = FALSE,
#'     include_self_comp = TRUE
#' )
#' sharing_heatmaps <- sharing_heatmap(sharing_df = sharing)
#' sharing_heatmaps$absolute
#' sharing_heatmaps$on_g1
#' sharing_heatmaps$on_union
sharing_heatmap <- function(sharing_df,
    show_on_x = "g1",
    show_on_y = "g2",
    absolute_sharing_col = "shared",
    title_annot = NULL,
    plot_relative_sharing = TRUE,
    rel_sharing_col = c("on_g1", "on_union"),
    show_perc_symbol_rel = TRUE,
    interactive = FALSE) {
    ## Check inputs
    stopifnot(is.data.frame(sharing_df))
    stopifnot(is.character(show_on_x))
    stopifnot(is.character(show_on_y))
    stopifnot(is.character(absolute_sharing_col))
    stopifnot(is.logical(plot_relative_sharing))
    if (plot_relative_sharing) {
        stopifnot(is.character(rel_sharing_col))
        if (!all(rel_sharing_col %in% colnames(sharing_df))) {
            rlang::abort(.missing_user_cols_error(
                rel_sharing_col[!rel_sharing_col %in% colnames(sharing_df)]
            ))
        }
    }
    stopifnot(is.null(title_annot) || is.character(title_annot))
    stopifnot(is.logical(show_perc_symbol_rel))
    stopifnot(is.logical(interactive))

    ### --- Absolute sharing
    heatmap_abs <- ggplot2::ggplot(
        sharing_df,
        ggplot2::aes_(
            x = as.name(show_on_x[1]),
            y = as.name(show_on_y[1]),
            fill = as.name(absolute_sharing_col[1]),
            alpha = as.name(absolute_sharing_col[1]),
            label = as.name(absolute_sharing_col[1])
        )
    ) +
        ggplot2::geom_raster() +
        ggplot2::scale_fill_gradientn(colours = c(
            "white", "gold",
            "navyblue"
        )) +
        ggplot2::scale_alpha_continuous(range = c(1, 0.5)) +
        ggplot2::geom_text(ggplot2::aes(alpha = 1)) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
            axis.text.y = ggplot2::element_text(hjust = 1),
            axis.title = ggplot2::element_blank(),
            legend.position = "none",
            plot.title = ggplot2::element_text(size = ggplot2::rel(2))
        ) +
        ggplot2::labs(
            title = paste(c(
                "IS sharing - absolute number",
                title_annot
            ), collapse = " - "),
            subtitle = paste(
                "Absolute number of shared IS",
                "between group pairs"
            )
        )
    ### --- Relative sharing
    heatmap_rel <- NULL
    if (plot_relative_sharing) {
        sharing_df_rounding <- sharing_df %>%
            dplyr::mutate(dplyr::across(
                .cols = dplyr::all_of(rel_sharing_col),
                .fns = ~ round(.x, digits = 2)
            ))
        plot_rel_heat <- function(col, df) {
            plot <- ggplot2::ggplot(
                df,
                ggplot2::aes_(
                    x = as.name(show_on_x[1]),
                    y = as.name(show_on_y[1]),
                    fill = as.name(col),
                    alpha = as.name(col)
                )
            ) +
                ggplot2::geom_raster() +
                ggplot2::scale_alpha_continuous(range = c(1, 0.5)) +
                ggplot2::theme(
                    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                    axis.text.y = ggplot2::element_text(hjust = 1),
                    axis.title = ggplot2::element_blank(),
                    legend.position = "none",
                    plot.title = ggplot2::element_text(size = ggplot2::rel(2))
                ) +
                ggplot2::labs(
                    title = paste(c(
                        "IS sharing - relative",
                        title_annot
                    ), collapse = " - "),
                    subtitle = paste(
                        "Percentage of shared IS",
                        "between group pairs"
                    )
                ) +
                ggplot2::scale_fill_gradientn(colours = c(
                    "white",
                    "gold",
                    "navyblue"
                ))

            if (show_perc_symbol_rel[1]) {
                plot <- plot +
                    ggplot2::geom_text(ggplot2::aes(
                        alpha = 1,
                        label = scales::percent(..fill..,
                            scale = 1,
                            accuracy = 0.1
                        )
                    ),
                    size = 3
                    )
            } else {
                plot <- plot +
                    ggplot2::geom_text(ggplot2::aes(
                        alpha = 1,
                        label = ..fill..
                    ),
                    size = 3
                    )
            }
            plot
        }
    }

    heatmap_rel <- purrr::map(
        unique(rel_sharing_col),
        ~ plot_rel_heat(.x, df = sharing_df_rounding)
    ) %>%
        purrr::set_names(rel_sharing_col)

    result <- list(absolute = heatmap_abs)
    result <- append(result, heatmap_rel)
    if (interactive) {
        if (!requireNamespace("plotly")) {
            rlang::inform(.missing_pkg_error("plotly"))
            rlang::inform("Returning static plots")
            return(result)
        }
        result <- purrr::map(result, plotly::ggplotly)
    }
    return(result)
}


#' Produce tables to plot sharing venn or euler diagrams.
#'
#' @description \lifecycle{experimental}
#' This function processes a sharing data frame obtained via `is_sharing()`
#' with the option `table_for_venn = TRUE` to obtain a list of objects
#' that can be plotted as venn or euler diagrams.
#'
#' @details
#' The functions requires the package
#' [eulerr](https://jolars.github.io/eulerr/index.html). Each row of the
#' input data frame is representable as a venn/euler diagram. The function
#' allows to specify a range of row indexes to obtain a list of plottable
#' objects all at once, leave it to NULL to process all rows.
#'
#' To actually plot the data it is sufficient to call the function `plot()`
#' and specify optional customization arguments. See
#' [eulerr docs](https://jolars.github.io/eulerr/reference/plot.euler.html)
#' for more detail on this.
#'
#' @param sharing_df The sharing data frame
#' @param row_range Either `NULL` or a numeric vector of row indexes (e.g.
#' `c(1, 4, 5)` will produce tables only for rows 1, 4 and 5)
#' @param euler If `TRUE` will produce tables for euler diagrams, otherwise
#' will produce tables for venn diagrams
#'
#' @family Plotting functions
#'
#' @return A list of data frames
#' @export
#'
#' @examples
#' data("integration_matrices", package = "ISAnalytics")
#' data("association_file", package = "ISAnalytics")
#' aggreg <- aggregate_values_by_key(
#'     x = integration_matrices,
#'     association_file = association_file,
#'     value_cols = c("seqCount", "fragmentEstimate")
#' )
#' sharing <- is_sharing(aggreg, n_comp = 3, table_for_venn = TRUE)
#' venn_tbls <- sharing_venn(sharing, row_range = 1:3, euler = FALSE)
#' venn_tbls
#' plot(venn_tbls[[1]])
sharing_venn <- function(sharing_df,
    row_range = NULL,
    euler = TRUE) {
    if (!requireNamespace("eulerr", quietly = TRUE)) {
        rlang::abort(.missing_pkg_error("eulerr"))
    }
    stopifnot(is.data.frame(sharing_df))
    stopifnot(is.null(row_range) ||
        is.numeric(row_range) || is.integer(row_range))
    stopifnot(is.logical(euler))
    # Check row range
    if (is.null(row_range)) {
        row_range <- seq_len(nrow(sharing_df))
    }
    # Check truth table
    if (!"truth_tbl_venn" %in% colnames(sharing_df)) {
        no_truth_tbl_msg <- c("No truth table column",
            x = paste(
                "The column 'truth_tbl_venn'",
                "is required but seems to be missing"
            ),
            i = paste(
                "Did you forget to call",
                "`is_sharing(..., table_for_venn",
                "= TRUE)`?"
            )
        )
        rlang::abort(no_truth_tbl_msg)
    }
    # Filter data
    filtered_df <- sharing_df[row_range]
    if (nrow(filtered_df) == 0) {
        rlang::inform("Empty table, nothing to compute")
        return(NULL)
    }
    fixed_tbls <- if (euler) {
        purrr::map(filtered_df$truth_tbl_venn, function(x) {
            as_matrix <- as.matrix(x, rownames = "int_id")
            eul <- eulerr::euler(as_matrix)
            eul
        })
    } else {
        purrr::map(filtered_df$truth_tbl_venn, function(x) {
            as_matrix <- as.matrix(x, rownames = "int_id")
            eul <- eulerr::venn(as_matrix)
            eul
        })
    }
    fixed_tbls
}

#' Trace a circos plot of genomic densities.
#'
#' @description \lifecycle{experimental}
#' For this functionality
#' the suggested package
#' [circlize](https://cran.r-project.org/web/packages/circlize/index.html)
#' is required.
#' Please note that this function is a simple wrapper of basic `circlize`
#' functions, for an in-depth explanation on how the functions work and
#' additional arguments please refer to the official documentation
#' [Circular Visualization in R](https://jokergoo.github.io/circlize_book/book/)
#'
#' @details
#' ## Providing genomic labels
#' If genomic labels should be plotted alongside genomic density tracks,
#' the user should provide them as a simple data frame in standard bed format,
#' namely `chr`, `start`, `end` plus a column containing the labels.
#' NOTE: if the user decides to plot on the default device (viewer in RStudio),
#' he must ensure there is enough space for all elements to be plotted,
#' otherwise an error message is thrown.
#'
#' @param data Either a single integration matrix or a list of integration
#' matrices. If a list is provided, a separate density track for each
#' data frame is plotted.
#' @param gene_labels Either `NULL` or a data frame in bed format. See details.
#' @param label_col Numeric index of the column of `gene_labels` that contains
#' the actual labels. Relevant only if `gene_labels` is not set to `NULL`.
#' @param cytoband_specie Specie for initializing the cytoband
#' @param track_colors Colors to give to density tracks. If more than one
#' integration matrix is provided as `data` should be of the same length.
#' Values are recycled if length of `track_colors` is smaller than the length
#' of the input data.
#' @param grDevice The graphical device where the plot should be traced.
#' `default`, if executing from RStudio is the viewer.
#' @param file_path If a device other than `default` is chosen, the path on
#' disk where the file should be saved. Defaults to
#' `{current directory}/circos_plot.{device}`.
#' @param ... Additional named arguments to pass on to chosen device,
#' `circlize::circos.par()`,
#' `circlize::circos.genomicDensity()` and `circlize::circos.genomicLabels()`
#'
#' @importFrom rlang abort dots_list .data arg_match fn_fmls_names as_function
#' @importFrom rlang exec
#' @importFrom dplyr select all_of distinct mutate rename
#' @importFrom purrr map walk2
#' @importFrom fs is_dir dir_exists dir_create path path_ext_set path_ext
#' @importFrom lubridate today
#'
#' @family Plotting functions
#' @return `NULL`
#' @export
#'
#' @examples
#' \donttest{
#' data("integration_matrices", package = "ISAnalytics")
#' data("association_file", package = "ISAnalytics")
#' aggreg <- aggregate_values_by_key(
#'     x = integration_matrices,
#'     association_file = association_file,
#'     value_cols = c("seqCount", "fragmentEstimate")
#' )
#' by_subj <- aggreg %>%
#'     dplyr::group_by(.data$SubjectID) %>%
#'     dplyr::group_split()
#' circos_genomic_density(by_subj,
#'     track_colors = c("navyblue", "gold"),
#'     grDevice = "default", track.height = 0.1
#' )
#' }
circos_genomic_density <- function(data,
    gene_labels = NULL,
    label_col = NULL,
    cytoband_specie = "hg19",
    track_colors = "navyblue",
    grDevice = c(
        "png", "pdf", "svg",
        "jpeg", "bmp", "tiff",
        "default"
    ),
    file_path = getwd(),
    ...) {
    if (!requireNamespace("circlize", quietly = TRUE)) {
        rlang::abort(.missing_pkg_error("circlize"))
    }
    stopifnot(is.list(data))
    mode <- if (is.data.frame(data)) {
        "DF"
    } else {
        "LIST"
    }
    stopifnot(is.null(gene_labels) || is.data.frame(gene_labels))
    if (!is.null(gene_labels)) {
        stopifnot(is.numeric(label_col) || is.integer(label_col))
        label_col <- label_col[1]
    }
    stopifnot(is.character(cytoband_specie))
    stopifnot(is.character(track_colors))
    dots <- rlang::dots_list(..., .named = TRUE, .homonyms = "first")
    ## -- Prep data
    .prep_dens_data <- function(df) {
        if (!.check_mandatory_vars(df)) {
            rlang::abort(.missing_mand_vars())
        }
        df %>%
            dplyr::select(dplyr::all_of(mandatory_IS_vars())) %>%
            dplyr::distinct() %>%
            dplyr::mutate(
                chr = paste0("chr", .data$chr),
                end = .data$integration_locus
            ) %>%
            dplyr::rename(start = "integration_locus") %>%
            dplyr::select(.data$chr, .data$start, .data$end, .data$strand)
    }
    data_mod <- if (mode == "DF") {
        .prep_dens_data(data)
    } else {
        purrr::map(data, .prep_dens_data)
    }
    device <- rlang::arg_match(grDevice)
    if (device != "default") {
        ## Open device
        stopifnot(is.character(file_path))
        if (fs::is_dir(file_path)) {
            if (!fs::dir_exists(file_path)) {
                fs::dir_create(file_path)
            }
            def <- paste0("circos_plot.", device)
            date <- lubridate::today()
            gen_filename <- paste0(date, "_", def)
            file_path <- fs::path(file_path, gen_filename)
        } else {
            if (fs::path_ext(file_path) == "") {
                file_path <- fs::path_ext_set(file_path, device)
            }
        }
        device_args <- dots[names(dots) %in%
            rlang::fn_fmls_names(
                rlang::as_function(device)
            )]
        device_args <- device_args[!names(device_args) %in% c("filename")]
        rlang::exec(.fn = device, filename = file_path, !!!device_args)
    }
    ## Draw plot
    par_args <- dots[names(dots) %in%
        rlang::fn_fmls_names(
            circlize::circos.par
        )]
    par_args <- par_args[!names(par_args) %in% c("start.degree")]
    density_args <- dots[names(dots) %in%
        rlang::fn_fmls_names(
            circlize::circos.genomicDensity
        )]
    density_args <- density_args[!names(density_args) %in% c("data", "col")]
    rlang::exec(.fn = circlize::circos.par, "start.degree" = 90, !!!par_args)
    circlize::circos.initializeWithIdeogram(species = cytoband_specie[1])
    if (mode == "DF") {
        rlang::exec(
            .fn = circlize::circos.genomicDensity,
            data = data_mod,
            col = track_colors[1],
            !!!density_args
        )
    } else {
        if (length(track_colors) < length(data_mod)) {
            track_colors <- rep_len(track_colors, length(data_mod))
        } else if (length(track_colors) > length(data_mod)) {
            track_colors <- track_colors[seq_along(data_mod)]
        }
        purrr::walk2(
            data_mod, track_colors,
            ~ rlang::exec(
                .fn = circlize::circos.genomicDensity,
                data = .x,
                col = .y,
                !!!density_args
            )
        )
    }

    if (!is.null(gene_labels)) {
        labs_args <- dots[names(dots) %in%
            rlang::fn_fmls_names(
                circlize::circos.genomicLabels
            )]
        labs_args <- labs_args[!names(labs_args) %in% c("bed", "labels.column")]
        rlang::exec(
            .fn = circlize::circos.genomicLabels,
            bed = gene_labels,
            labels.column = label_col,
            !!!labs_args
        )
    }

    if (device != "default") {
        grDevices::dev.off()
    }
    circlize::circos.clear()
}
