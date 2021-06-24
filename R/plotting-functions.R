#------------------------------------------------------------------------------#
# Plotting functions
#------------------------------------------------------------------------------#

#' Trace volcano plot for computed CIS data.
#'
#' \lifecycle{experimental}
#' Traces a volcano plot for IS frequency and CIS results.
#'
#' @details
#' ## Input data frame
#' Users can supply as `x` either a simple integration matrix or a
#' data frame resulting from the call to \link{CIS_grubbs}
#' with `add_standard_padjust = TRUE`. In the first case an internal call to
#' the function `CIS_grubbs` is performed.
#'
#' ## Oncogene and tumor suppressor genes files
#' These files are included in the package for user convenience and are
#' simply UniProt files with gene annotations for human and mouse.
#' For more details on how this files were generated use the help `?filename`
#' function.
#'
#' ## Known oncogenes
#' The default values are contained in a data frame exported by this package,
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
#' @param onco_db_file Uniprot file for proto-oncogenes (see details)
#' @param tumor_suppressors_db_file Uniprot file for tumor-suppressor genes
#' @param species One between "human", "mouse" and "all"
#' @param known_onco Data frame with known oncogenes. See details.
#' @param suspicious_genes Data frame with clinical relevant suspicious
#' genes. See details.
#' @param significance_threshold The significance threshold
#' @param annotation_threshold_ontots Value above which genes are annotated
#' @param highlight_genes Either NULL or a character vector of genes to be
#' highlighted in the plot even if they're not above the threshold
#' @param title_prefix A string to be displayed in the title - usually the
#' project name and other characterizing info
#' @param return_df Return the data frame used to generate the plot?
#' This can be useful if the user wants to manually modify the plot with
#' ggplot2. If TRUE the function returns a list containing both the plot
#' and the data frame.
#'
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#' @import dplyr
#'
#' @return A plot or a list containing a plot and a data frame
#' @export
#'
#' @examples
#' op <- options(ISAnalytics.widgets = FALSE)
#'
#' path_AF <- system.file("extdata", "ex_association_file.tsv",
#'     package = "ISAnalytics"
#' )
#' root_correct <- system.file("extdata", "fs.zip",
#'     package = "ISAnalytics"
#' )
#' root_correct <- unzip_file_system(root_correct, "fs")
#'
#' matrices <- import_parallel_Vispa2Matrices_auto(
#'     association_file = path_AF, root = root_correct,
#'     quantification_type = c("seqCount", "fragmentEstimate"),
#'     matrix_type = "annotated", workers = 2, patterns = NULL,
#'     matching_opt = "ANY",
#'     dates_format = "dmy"
#' )
#'
#' cis <- CIS_grubbs(matrices)
#' plot <- CIS_volcano_plot(cis)
#' options(op)
CIS_volcano_plot <- function(x,
    onco_db_file = system.file("extdata",
        "201806_uniprot-Proto-oncogene.tsv.xz",
        package = "ISAnalytics"
    ),
    tumor_suppressors_db_file = system.file("extdata",
        "201806_uniprot-Tumor-suppressor.tsv.xz",
        package = "ISAnalytics"
    ),
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
    stopifnot(is.character(onco_db_file) & length(onco_db_file) == 1)
    stopifnot(is.character(tumor_suppressors_db_file) &
        length(tumor_suppressors_db_file) == 1)
    stopifnot(is.character(species))
    stopifnot(is.data.frame(known_onco))
    stopifnot(is.data.frame(suspicious_genes))
    stopifnot(is.numeric(significance_threshold) |
        is.integer(significance_threshold) &
            length(significance_threshold) == 1)
    stopifnot(is.numeric(annotation_threshold_ontots) |
        is.integer(annotation_threshold_ontots) &
            length(annotation_threshold_ontots) == 1)
    stopifnot(is.null(title_prefix) || (is.character(title_prefix) &
        length(title_prefix == 1)))
    stopifnot(is.null(highlight_genes) || is.character(highlight_genes))
    stopifnot(is.logical(return_df))
    if (is.null(title_prefix)) {
        title_prefix <- ""
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
            message(paste("Calculating CIS_grubbs for x..."))
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
                "(FDR adjusted; ",
                "-log = ", (round(-log(significance_threshold, base = 10), 3)),
                ").\nOnco/TS genes source: UniProt (other genes ",
                "labeled as 'Other'). Annotated if P-value > ",
                round(annotation_threshold_ontots_log, 3), "\nexcept ",
                "selected genes to be highlighted"
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
#' @import dplyr
#' @import ggplot2
#'
#' @return A plot
#' @export
#'
#' @examples
#' op <- options("ISAnalytics.widgets" = FALSE, "ISAnalytics.verbose" = FALSE)
#' path_AF <- system.file("extdata", "ex_association_file.tsv",
#'     package = "ISAnalytics"
#' )
#' root_correct <- system.file("extdata", "fs.zip", package = "ISAnalytics")
#' root_correct <- unzip_file_system(root_correct, "fs")
#' association_file <- import_association_file(path_AF, root_correct,
#'     dates_format = "dmy"
#' )
#' aggregated_meta <- aggregate_metadata(association_file)
#' matrices <- import_parallel_Vispa2Matrices_auto(
#'     association_file = association_file, root = NULL,
#'     quantification_type = c("fragmentEstimate", "seqCount"),
#'     matrix_type = "annotated", workers = 2, matching_opt = "ANY"
#' )
#' agg <- aggregate_values_by_key(
#'     x = matrices,
#'     association_file = association_file,
#'     value_cols = "seqCount"
#' )
#' estimate <- HSC_population_size_estimate(
#'     x = agg,
#'     metadata = aggregated_meta,
#'     stable_timepoints = NULL
#' )
#' p <- HSC_population_plot(estimate, "PROJECT1")
#' options(op)
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
#' in natural language roughly as "highlight with colours only those
#' integrations (alluvia) that at least in 1 point in time have an
#' abundance value >= 1 %". The remaining integrations will be plotted in
#' grey in the strata.
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
#' @import dplyr
#' @import BiocParallel
#' @importFrom purrr set_names
#' @importFrom tidyr unite
#'
#' @return For each group a list with the associated plot and optionally
#' the summary tableGrob
#' @export
#'
#' @examples
#' op <- options("ISAnalytics.widgets" = FALSE, "ISAnalytics.verbose" = FALSE)
#' path_AF <- system.file("extdata", "ex_association_file.tsv",
#'     package = "ISAnalytics"
#' )
#' root_correct <- system.file("extdata", "fs.zip", package = "ISAnalytics")
#' root_correct <- unzip_file_system(root_correct, "fs")
#' association_file <- import_association_file(path_AF, root_correct,
#'     dates_format = "dmy"
#' )
#' matrices <- import_parallel_Vispa2Matrices_auto(
#'     association_file = association_file, root = NULL,
#'     quantification_type = c("fragmentEstimate", "seqCount"),
#'     matrix_type = "annotated", workers = 2, matching_opt = "ANY"
#' )
#' agg <- aggregate_values_by_key(
#'     x = matrices,
#'     association_file = association_file,
#'     value_cols = c("fragmentEstimate", "seqCount")
#' )
#' abundance <- compute_abundance(agg,
#'     columns = "fragmentEstimate_sum",
#'     key = c("SubjectID", "CellMarker", "Tissue", "TimePoint")
#' )
#' plots <- integration_alluvial_plot(abundance, top_abundant_tbl = FALSE)
#' options(op)
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
    p <- if (.Platform$OS.type == "windows") {
        BiocParallel::SnowParam(
            stop.on.error = FALSE,
            progressbar = TRUE,
            tasks = length(groups_to_plot),
            exportglobals = FALSE
        )
    } else {
        BiocParallel::MulticoreParam(
            stop.on.error = FALSE, progressbar = TRUE,
            tasks = length(groups_to_plot), exportglobals = FALSE
        )
    }

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
                            NULL
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
#' @import dplyr
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
    tbl <- tbl[, setNames(.(sum(.SD)), plot_y),
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
.alluvial_plot <- function(tbl, plot_x, plot_y) {
    max_y <- max(tbl[[plot_y]])
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
            na.rm = TRUE,
            decreasing = FALSE,
            alpha = .75,
            aes.bind = "alluvia"
        ) +
        ggplot2::geom_text(ggplot2::aes_(
            x = ggplot2::sym(plot_x),
            y = max_y + 5,
            label = ~count
        )) +
        ggplot2::scale_fill_viridis_d() +
        ggplot2::theme(legend.position = "none")
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
#' @return A tableGrob object
#' @export
#'
#' @examples
#' op <- options("ISAnalytics.widgets" = FALSE, "ISAnalytics.verbose" = FALSE)
#' path_AF <- system.file("extdata", "ex_association_file.tsv",
#'     package = "ISAnalytics"
#' )
#' root_correct <- system.file("extdata", "fs.zip", package = "ISAnalytics")
#' root_correct <- unzip_file_system(root_correct, "fs")
#' association_file <- import_association_file(path_AF, root_correct,
#'     dates_format = "dmy"
#' )
#' matrices <- import_parallel_Vispa2Matrices_auto(
#'     association_file = association_file, root = NULL,
#'     quantification_type = c("fragmentEstimate", "seqCount"),
#'     matrix_type = "annotated", workers = 2, matching_opt = "ANY"
#' )
#' agg <- aggregate_values_by_key(
#'     x = matrices,
#'     association_file = association_file,
#'     value_cols = c("fragmentEstimate", "seqCount")
#' )
#' abundance <- compute_abundance(agg,
#'     columns = "fragmentEstimate_sum",
#'     key = c("SubjectID", "CellMarker", "Tissue", "TimePoint")
#' )
#' grob <- top_abund_tableGrob(abundance)
#' gridExtra::grid.arrange(grob)
#' options(op)
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
#' @return A list of plots or widgets
#' @seealso \link{is_sharing}
#' @export
#'
#' @importFrom rlang abort inform .data
#' @import ggplot2
#' @importFrom dplyr mutate across all_of
#' @importFrom purrr map set_names
#'
#' @examples
#' path <- system.file("extdata", "ex_annotated_ISMatrix.tsv.xz",
#'     package = "ISAnalytics"
#' )
#' matrix <- import_single_Vispa2Matrix(path)
#' sharing <- is_sharing(matrix, group_key = "CompleteAmplificationID")
#' heatmaps <- sharing_heatmap(sharing$sharing)
#' heatmaps$absolute
sharing_heatmap <- function(sharing_df,
    show_on_x = "group1",
    show_on_y = "group2",
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
