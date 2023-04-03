#------------------------------------------------------------------------------#
# Plotting functions
#------------------------------------------------------------------------------#

#' Trace volcano plot for computed CIS data.
#'
#' @description
#' `r lifecycle::badge("stable")`
#' Traces a volcano plot for IS frequency and CIS results.
#'
#' @details
#' ## Input data frame
#' Users can supply as `x` either a simple integration matrix or a
#' data frame resulting from the call to \link{CIS_grubbs}.
#' In the first case an internal call to
#' the function `CIS_grubbs()` is performed.
#'
#' @family Plotting functions
#'
#' @param x Either a simple integration matrix or a data frame resulting
#' from the call to \link{CIS_grubbs} with `add_standard_padjust = TRUE`
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
#' @template genes_db
#'
#' @section Required tags:
#' The function will explicitly check for the presence of these tags:
#'
#' ```{r echo=FALSE, results="asis"}
#' all_tags <- available_tags()
#' needed <- all_tags |>
#'    dplyr::mutate(
#'    in_fun = purrr::map_lgl(.data$needed_in,
#'    ~ "CIS_volcano_plot" %in% .x)
#'    ) |>
#'    dplyr::filter(in_fun == TRUE) |>
#'    dplyr::pull("tag")
#'  cat(paste0("* ", needed, collapse="\n"))
#' ```
#'
#' @importFrom rlang .data
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
CIS_volcano_plot <- function(
        x,
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
    ## Check if CIS function was already called
    min_cis_col <- c(
        "tdist_bonferroni_default", "tdist_fdr",
        "neg_zscore_minus_log2_int_freq_tolerance"
    )
    cis_grubbs_df <- if (!all(min_cis_col %in% colnames(x))) {
        if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
            rlang::inform("Calculating CIS_grubbs for x...")
        }
        (CIS_grubbs(x, return_missing_as_df = FALSE))$cis
    } else {
        x
    }
    gene_sym_col <- annotation_IS_vars(TRUE) |>
        dplyr::filter(.data$tag == "gene_symbol") |>
        dplyr::pull(.data$names)
    ## Add info to CIS
    cis_grubbs_df <- .expand_cis_df(
        cis_grubbs_df, gene_sym_col,
        onco_db_file, tumor_suppressors_db_file,
        species, known_onco, suspicious_genes
    )
    cis_grubbs_df <- cis_grubbs_df |>
        dplyr::mutate(minus_log_p = -log(.data$tdist_bonferroni_default,
            base = 10
        ))
    cis_grubbs_df <- cis_grubbs_df |>
        dplyr::mutate(
            minus_log_p_fdr = -log(.data$tdist_fdr, base = 10),
            positive_outlier_and_significant = ifelse(
                test = !is.na(.data$tdist_fdr) &
                    .data$tdist_fdr < significance_threshold,
                yes = TRUE,
                no = FALSE
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
        ggplot2::aes(
            y = .data[["minus_log_p_fdr"]],
            x = .data[["neg_zscore_minus_log2_int_freq_tolerance"]],
            color = .data[["KnownGeneClass"]],
            fill = .data[["KnownGeneClass"]]
        ),
        na.rm = TRUE, se = TRUE
    ) +
        ggplot2::geom_point(alpha = .5, size = 3) +
        ggplot2::geom_hline(
            yintercept = significance_threshold_minus_log_p,
            color = "black", linewidth = 1, show.legend = TRUE,
            linetype = "dashed"
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
            ggplot2::aes(label = .data[[gene_sym_col]]),
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
                "Volcano plot of IS gene frequency and",
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
        to_highlight <- cis_grubbs_df |>
            dplyr::filter(
                stringr::str_to_lower(.data[[gene_sym_col]]) %in%
                    stringr::str_to_lower(highlight_genes),
                .data$tdist_fdr >= significance_threshold
            )
        plot_cis_fdr_slice <- plot_cis_fdr_slice +
            ggrepel::geom_label_repel(
                data = to_highlight,
                ggplot2::aes(label = .data[[gene_sym_col]]),
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

#' Heatmaps for the top N common insertion sites over time.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#' This function computes the visualization of the results of the function
#' `CIS_grubbs_overtime()` in the form of heatmaps for the top N selected
#' genes over time.
#'
#' @template genes_db
#' @family Plotting functions
#'
#' @inheritParams CIS_volcano_plot
#' @param x Output of the function `CIS_grubbs_overtime()`, either in single
#' data frame form or nested lists
#' @param n_genes Number of top genes to consider
#' @param timepoint_col The name of the time point column in `x`
#' @param group_col The name of the group column in `x`
#' @param plot_values Which kind of values should be plotted? Can either be
#' `"p"` for the p-value or `"minus_log_p"` for a scaled p-value of the
#' Grubbs test
#' @param p_value_correction One among `"bonferroni"` and `"fdr"`
#' @param prune_tp_treshold Minimum number of genes to retain a time point.
#' See details.
#' @param gene_selection_param The descriptive statistic measure to decide
#' which genes to plot, possible choices are
#' `"trimmed", "n", "mean", "sd", "median","mad", "min", "max"`. See details.
#' @param fill_0_selection Fill NA values with 0s before computing statistics
#' for each gene? (TRUE/FALSE)
#' @param fill_NA_in_heatmap Fill NA values with 0 when plotting the heatmap?
#' (TRUE/FALSE)
#' @param heatmap_color_palette Colors for values in the heatmaps,
#' either `"default"` or a function producing
#' a color palette, obtainable via `grDevices::colorRampPalette`.
#' @param title_generator Either `NULL` or a function. See details.
#' @param save_as_files Should heatmaps be saved to files on disk? (TRUE/FALSE)
#' @param files_format The extension of the files produced, supported
#' formats are `"pdf", "png", "tiff", "bmp", "jpg"`. Relevant only if
#' `files_format = TRUE`
#' @param folder_path Path to the folder where files will be saved
#' @param ... Other params to pass to `pheatmap::pheatmap`
#'
#' @details
#' ## Top N gene selection
#' Since the genes present in different time point slices are likely different,
#' the decision process to select the final top N genes to represent in the
#' heatmap follows this logic:
#'
#' * Each time point slice is arranged either in ascending order (if we want to
#' plot the p-value) or in descending order (if we want to plot the scaled
#' p-value) and the top n genes are selected
#' * A series of statistics are computed over the union set of genes on ALL
#' time points (min, max, mean, ...)
#' * A decision is taken by considering the ordered `gene_selection_param`
#' (order depends once again if the values are scaled or not), and the first
#' N genes are selected for plotting.
#'
#' ### Filling NA values prior calculations
#' It is possible to fill NA values (aka missing combinations of GENE/TP) with
#' 0s prior computing the descriptive statistics on which gene selection is
#' based. Please keep in mind that this has an impact on the final result,
#' since for computing metrics such as the mean, NA values are usually removed,
#' decreasing the overall number of values considered - this does not hold
#' when NA values are substituted with 0s.
#'
#' ### The statistics
#' Statistics are computed for each gene over all time points of each group.
#' More in detail, `n`: counts the number of instances (rows)
#' in which the genes appears, aka it counts the time points in which the gene
#' is present. NOTE: if
#' `fill_0_selection` option is set to `TRUE` this value will be equal for
#' all genes! All other statistics as per the argument `gene_selection_param`
#' map to the corresponding R functions with the exception of `trimmed` which
#' is a simple call to the `mean` function with the argument `trimmed = 0.1`.
#'
#' ## Aesthetics
#' It is possible to customise the appearence of the plot through different
#' parameters.
#'
#' * `fill_NA_in_heatmap` tells the function whether missing combinations of
#' GENE/TP should be plotted as NA or filled with a value (1 if p-value, 0
#' if scaled p-value)
#' * A title generator function can be provided to dynamically create a title
#' for the plots: the function can accept two positional arguments for
#' the group identifier and the number of selected genes respectively. If one or
#' none of the arguments are of interest, they can be absorbed with `...`.
#' * `heatmap_color_palette` can be used to specify a function from which
#' colors are sampled (refers to the colors of values only)
#' * To change the colors associated with annotations instead, use the
#' argument `annotation_colors` of `pheatmap::pheatmap()` - it must be set to a
#' list with this format:
#' ```
#' list(
#'   KnownGeneClass = c("OncoGene" = color_spec,
#'                      "Other" = color_spec,
#'                      "TumSuppressor" = color_spec),
#'   ClinicalRelevance = c("TRUE" = color_spec,
#'                         "FALSE" = color_spec),
#'   CriticalForInsMut = c("TRUE" = color_spec,
#'                         "FALSE" = color_spec)
#' )
#' ```
#'
#' @return Either a list of graphical objects or a list of paths where
#' plots were saved
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
#' cis_overtime <- CIS_grubbs_overtime(aggreg)
#' hmaps <- top_cis_overtime_heatmap(cis_overtime$cis,
#'     fill_NA_in_heatmap = TRUE
#' )
#'
#' # To re-plot:
#' # grid::grid.newpage()
#' # grid::grid.draw(hmaps$PT001$gtable)
top_cis_overtime_heatmap <- function(
        x,
        n_genes = 20,
        timepoint_col = "TimePoint",
        group_col = "group",
        onco_db_file = "proto_oncogenes",
        tumor_suppressors_db_file = "tumor_suppressors",
        species = "human",
        known_onco = known_clinical_oncogenes(),
        suspicious_genes =
            clinical_relevant_suspicious_genes(),
        significance_threshold = 0.05,
        plot_values = c("minus_log_p", "p"),
        p_value_correction = c("fdr", "bonferroni"),
        prune_tp_treshold = 20,
        gene_selection_param = c(
            "trimmed", "n", "mean", "sd", "median",
            "mad", "min", "max"
        ),
        fill_0_selection = TRUE,
        fill_NA_in_heatmap = FALSE,
        heatmap_color_palette = "default",
        title_generator = NULL,
        save_as_files = FALSE,
        files_format = c("pdf", "png", "tiff", "bmp", "jpg"),
        folder_path = NULL,
        ...) {
    # --- Preliminary checks
    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        rlang::abort(.missing_pkg_error("pheatmap"))
    }
    stopifnot(is.list(x))
    stopifnot(is.numeric(n_genes))
    if (is.list(x) & !is.data.frame(x) & is.null(names(x))) {
        err_named <- c("Input list must have names",
            x = paste(
                "Input should follow the output format of",
                "`CIS_grubbs_overtime()`"
            )
        )
        rlang::abort(err_named)
    } else if (is.data.frame(x) &
        !all(c(timepoint_col, group_col) %in% colnames(x))) {
        err_df <- c("Input df is missing columns",
            x = paste(
                "Input should follow the output format of",
                "`CIS_grubbs_overtime()`"
            ),
            x = paste(
                "Time point column ('", timepoint_col,
                "') and/or group column ('", group_col,
                "') are missing"
            )
        )
        rlang::abort(err_df)
    }
    p_value_correction <- rlang::arg_match(p_value_correction)
    plot_values <- rlang::arg_match(plot_values)
    values_to_plot <- if (plot_values == "p") {
        paste0("tdist_", p_value_correction)
    } else {
        paste0("minus_log_p_", p_value_correction)
    }
    gene_selection_param <- rlang::arg_match(gene_selection_param)
    stopifnot(is.numeric(prune_tp_treshold))
    stopifnot(is.logical(fill_0_selection))
    fill_0_selection <- fill_0_selection[1]
    stopifnot(is.logical(fill_NA_in_heatmap))
    fill_NA_in_heatmap <- fill_NA_in_heatmap[1]
    stopifnot(is.logical(save_as_files))
    save_as_files <- save_as_files[1]
    stopifnot(is.null(folder_path) || is.character(folder_path))
    if (is.character(folder_path[1])) {
        fs::dir_create(folder_path[1])
    }
    stopifnot(
        (is.character(heatmap_color_palette) &
            heatmap_color_palette == "default") ||
            is.function(heatmap_color_palette)
    )
    files_format <- rlang::arg_match(files_format)
    if (save_as_files == TRUE & is.null(folder_path) &
        getOption("ISAnalytics.verbose", TRUE) == TRUE) {
        warn_msg <- c("Warning: you did not set a folder for saving files",
            i = "Returning heatmaps as a list in R env"
        )
        rlang::inform(warn_msg, class = "no_fold_files")
    }
    stopifnot(is.null(title_generator) || is.function(title_generator))
    dots <- rlang::dots_list(..., .named = TRUE)
    # --- If input is in list form, convert in single df
    if (is.list(x) & !is.data.frame(x)) {
        x <- purrr::map2(x, names(x), ~ {
            tmp <- .x
            if (is.list(tmp) & !is.data.frame(tmp)) {
                purrr::map2(
                    tmp, names(tmp),
                    ~ .x |> dplyr::mutate(!!timepoint_col := .y)
                ) |>
                    purrr::list_rbind() |>
                    dplyr::mutate(!!group_col := .y)
            } else if (is.data.frame(tmp)) {
                tmp |> dplyr::mutate(!!group_col := .y)
            } else {
                non_list_el <- c("Element is not list or df",
                    x = paste(
                        "Element", .y, "in x is not a list or a",
                        "data frame"
                    )
                )
                rlang::abort(non_list_el)
            }
        }) |>
            purrr::list_rbind()
    }
    gene_sym_col <- annotation_IS_vars(TRUE) |>
        dplyr::filter(.data$tag == "gene_symbol") |>
        dplyr::pull(.data$names)
    # --- Expand the df with gene info
    expanded <- .expand_cis_df(
        x, gene_sym_col,
        onco_db_file, tumor_suppressors_db_file,
        species, known_onco, suspicious_genes
    )
    # --- Add log conversions if needed
    if (plot_values == "minus_log_p") {
        expanded <- expanded |>
            dplyr::mutate(
                minus_log_p_bonferroni = -log(.data$tdist_bonferroni,
                    base = 10
                ),
                minus_log_p_fdr = -log(.data$tdist_fdr, base = 10)
            )
    }
    # --- For each combo (group, tp) arrange and slice the top n
    arrange_slice_top <- function(group_df, ...) {
        if (nrow(group_df) < prune_tp_treshold) {
            return(NULL)
        }
        if (plot_values == "p") {
            return(
                group_df |>
                    dplyr::arrange(.data[[values_to_plot]]) |>
                    dplyr::slice_head(n = n_genes)
            )
        } else {
            return(
                group_df |>
                    dplyr::arrange(dplyr::desc(.data[[values_to_plot]])) |>
                    dplyr::slice_head(n = n_genes)
            )
        }
    }
    slice_groups_tps <- expanded |>
        dplyr::select(
            dplyr::all_of(c(
                gene_sym_col, timepoint_col,
                group_col, values_to_plot
            ))
        ) |>
        dplyr::group_by(dplyr::across(dplyr::all_of(
            c(group_col, timepoint_col)
        ))) |>
        dplyr::group_map(.f = arrange_slice_top, .keep = TRUE) |>
        purrr::reduce(dplyr::bind_rows)

    groups <- unique(slice_groups_tps[[group_col]])

    # --- Evaluate statistics for genes in top n slices
    eval_candidates <- function(group_id) {
        df <- slice_groups_tps |>
            dplyr::filter(.data[[group_col]] == group_id) |>
            dplyr::select(-.data[[group_col]])
        if (fill_0_selection == TRUE) {
            value_fill <- list(0)
            names(value_fill) <- values_to_plot
            df <- df |>
                tidyr::complete(.data[[timepoint_col]], .data[[gene_sym_col]],
                    fill = value_fill
                )
        }
        df |>
            dplyr::group_by(dplyr::across(dplyr::all_of(gene_sym_col))) |>
            dplyr::summarise(
                n = dplyr::n(),
                mean = mean(.data[[values_to_plot]], na.rm = TRUE),
                sd = stats::sd(.data[[values_to_plot]], na.rm = TRUE),
                median = stats::median(.data[[values_to_plot]], na.rm = TRUE),
                trimmed = mean(.data[[values_to_plot]],
                    na.rm = TRUE,
                    trim = 0.1
                ),
                mad = stats::mad(.data[[values_to_plot]], na.rm = TRUE),
                min = min(.data[[values_to_plot]], na.rm = TRUE),
                max = max(.data[[values_to_plot]], na.rm = TRUE),
                .groups = "drop"
            )
    }
    candidates <- purrr::map(groups, eval_candidates) |>
        purrr::set_names(groups)

    # --- Select from candidates according to gene_selection_param
    select_from_candidates <- function(group_df) {
        if (gene_selection_param == "n") {
            return(
                group_df |>
                    dplyr::arrange(dplyr::desc(.data$n)) |>
                    dplyr::slice_head(n = n_genes) |>
                    dplyr::pull(.data[[gene_sym_col]])
            )
        }
        if (plot_values == "p") {
            ## Order ascending
            return(
                group_df |>
                    dplyr::arrange(.data[[gene_selection_param]]) |>
                    dplyr::slice_head(n = n_genes) |>
                    dplyr::pull(.data[[gene_sym_col]])
            )
        } else {
            ## Order descending
            return(
                group_df |>
                    dplyr::arrange(dplyr::desc(
                        .data[[gene_selection_param]]
                    )) |>
                    dplyr::slice_head(n = n_genes) |>
                    dplyr::pull(.data[[gene_sym_col]])
            )
        }
    }
    gene_selection <- purrr::map(candidates, select_from_candidates)

    # --- Extract only relevant genes from input
    genes_to_map <- purrr::map(groups, ~ expanded |>
        dplyr::filter(group == .x) |>
        dplyr::filter(.data[[timepoint_col]] %in%
            unique((slice_groups_tps |>
                dplyr::filter(group == .x))[[timepoint_col]])) |>
        dplyr::filter(.data[[gene_sym_col]] %in%
            gene_selection[[.x]]) |>
        dplyr::select(
            dplyr::all_of(c(
                gene_sym_col, timepoint_col,
                group_col, values_to_plot
            )),
            .data$CriticalForInsMut, .data$KnownGeneClass,
            .data$ClinicalRelevance
        )) |>
        purrr::set_names(groups)

    # --- Obtain heatmaps
    ## --- Define global defaults
    if (is.character(heatmap_color_palette)) {
        heatmap_color_palette <- if (plot_values == "minus_log_p") {
            grDevices::colorRampPalette(
                c(
                    "steelblue", "white", "red", "firebrick", "firebrick",
                    "darkred",
                    "darkred", "violet", "violet"
                )
            )
        } else {
            grDevices::colorRampPalette(
                rev(c(
                    "steelblue", "white", "red", "firebrick", "darkred",
                    "violet"
                )),
                bias = 10
            )
        }
    }
    annotation_palette <- if ("annotation_colors" %in% names(dots)) {
        dots$annotation_colors
    } else {
        list(
            KnownGeneClass = c(
                "OncoGene" = "red2",
                "Other" = "palegreen",
                "TumSuppressor" = "dodgerblue4"
            ),
            ClinicalRelevance = c(
                "TRUE" = "gold",
                "FALSE" = "gray90"
            ),
            CriticalForInsMut = c(
                "TRUE" = "red2",
                "FALSE" = "gray90"
            )
        )
    }
    plotting_step <- if (plot_values == "minus_log_p") {
        round((-log(0.05, base = 10) / 3), 3)
    } else {
        round(1 / 100, 3)
    }

    trace_heatmap <- function(data_df) {
        ## --- Fix annotations (fill na, convert to char)
        data_df <- data_df |>
            tidyr::replace_na(list(ClinicalRelevance = FALSE)) |>
            dplyr::mutate(
                CriticalForInsMut = as.character(.data$CriticalForInsMut),
                ClinicalRelevance = as.character(.data$ClinicalRelevance)
            )
        ## --- Obtain data matrix and annotations
        wide <- if (fill_NA_in_heatmap == TRUE) {
            data_df |>
                tidyr::pivot_wider(
                    names_from = timepoint_col,
                    values_from = values_to_plot,
                    values_fill = dplyr::if_else(plot_values == "p",
                        1, 0
                    )
                )
        } else {
            data_df |>
                tidyr::pivot_wider(
                    names_from = timepoint_col,
                    values_from = values_to_plot
                )
        }
        matrix <- wide |>
            dplyr::select(dplyr::all_of(unique(data_df[[timepoint_col]]))) |>
            as.matrix()
        rownames(matrix) <- wide[[gene_sym_col]]
        annotations <- wide |>
            dplyr::select(
                .data$CriticalForInsMut, .data$KnownGeneClass,
                .data$ClinicalRelevance
            ) |>
            as.data.frame()
        rownames(annotations) <- wide[[gene_sym_col]]
        ## --- Obtain other params
        plot_breaks <- if (plot_values == "minus_log_p") {
            seq(
                0, ceiling(max(data_df[[values_to_plot]], na.rm = TRUE)),
                plotting_step
            )
        } else {
            seq(0, 1, plotting_step)
        }
        params_to_pass <- list()
        params_to_pass$color <- heatmap_color_palette(length(plot_breaks) + 1)
        params_to_pass$annotation_colors <- annotation_palette
        params_to_pass$annotation_row <- annotations
        params_to_pass$breaks <- plot_breaks
        if (save_as_files == TRUE & !is.null(folder_path)) {
            file_name <- paste0(
                data_df$group[1], ".", lubridate::today(),
                "_top", n_genes, "-CIS-overtime_using-",
                gene_selection_param, ".", files_format
            )
            params_to_pass$filename <- fs::path(folder_path[1], file_name)
        }
        params_to_pass$cluster_rows <- if ("cluster_rows" %in% names(dots)) {
            dots$cluster_rows
        } else {
            TRUE
        }
        params_to_pass$cluster_cols <- if ("cluster_cols" %in% names(dots)) {
            dots$cluster_cols
        } else {
            FALSE
        }
        params_to_pass$scale <- if ("scale" %in% names(dots)) {
            dots$scale
        } else {
            "none"
        }
        params_to_pass$display_numbers <- if ("display_numbers" %in%
            names(dots)) {
            dots$display_numbers
        } else {
            TRUE
        }
        params_to_pass$number_format <- if ("number_format" %in% names(dots)) {
            dots$number_format
        } else {
            "%.2f"
        }
        params_to_pass$main <- if (!is.null(title_generator)) {
            title_generator(data_df$group[1], n_genes)
        } else {
            method_str <- if (plot_values == "minus_log_p") {
                paste0(
                    "-log(p-value/", p_value_correction, ")\n",
                    "[CIS iif p-value < 0.05; -log(0.05) = ",
                    (round(-log(0.05, base = 10), 3)), "]"
                )
            } else {
                paste(
                    "p-value/", p_value_correction, "\n",
                    "[CIS iif p-value < 0.05]"
                )
            }
            paste0(
                "Patient ", data_df$group[1],
                " - Top ", n_genes,
                " hotspot genes\nAnalysis over time, ",
                method_str
            )
        }
        dots <- dots[!names(dots) %in% c(names(params_to_pass), "mat")]
        params_to_pass <- append(params_to_pass, dots)
        ## --- Plot
        map <- rlang::exec(pheatmap::pheatmap, mat = matrix, !!!params_to_pass)
        if ("filename" %in% names(params_to_pass)) {
            map <- params_to_pass[["filename"]]
        }
        return(map)
    }

    purrr::map(genes_to_map, trace_heatmap)
}

#' Plot results of gene frequency Fisher's exact test.
#'
#' @description
#' `r lifecycle::badge("stable")`
#' Plots results of Fisher's exact test on gene frequency obtained via
#' `gene_frequency_fisher()` as a scatterplot.
#'
#' @details
#' ## Specifying genes to avoid highlighting
#'
#' In some cases, users might want to avoid highlighting certain genes
#' even if their p-value is below the threshold. To do so, use the
#' argument `do_not_highlight`: character vectors are appropriate for specific
#' genes that are to be excluded, expressions or lambdas allow a finer control.
#' For example we can supply:
#' ```{r eval=FALSE}
#' expr <- rlang::expr(!stringr::str_starts(GeneName, "MIR") &
#'                       average_TxLen_1 >= 300)
#' ```
#'
#' with this expression, genes that have a p-value < threshold and start with
#' "MIR" or have an average_TxLen_1 lower than 300 are excluded from the
#' highlighted points.
#' NOTE: keep in mind that expressions are evaluated inside a `dplyr::filter`
#' context.
#'
#' Similarly, lambdas are passed to the filtering function but only operate
#' on the column containing the gene symbol.
#' ```{r eval=FALSE}
#' lambda <- ~ stringr::str_starts(.x, "MIR")
#' ```
#'
#' @param fisher_df Test results obtained via `gene_frequency_fisher()`
#' @param p_value_col Name of the column containing the p-value to consider
#' @param annot_threshold Annotate with a different color if a point is below
#' the significance threshold. Single numerical value.
#' @param annot_color The color in which points below the threshold should be
#' annotated
#' @param gene_sym_col The name of the column containing the gene symbol
#' @param do_not_highlight Either `NULL`, a character vector, an expression
#' or a purrr-style lambda. Tells the function to ignore the highlighting
#' and labeling of these genes even if their p-value is below the threshold.
#' See details.
#' @param keep_not_highlighted If present, how should not highlighted genes
#' be treated? If set to `TRUE` points are plotted and colored with the
#' chosen color scale. If set to `FALSE` the points are removed entirely from
#' the plot.
#'
#' @family Plotting functions
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
#' cis <- CIS_grubbs(aggreg, by = "SubjectID")
#' fisher <- gene_frequency_fisher(cis$cis$PT001, cis$cis$PT002,
#'     min_is_per_gene = 2
#' )
#' fisher_scatterplot(fisher)
fisher_scatterplot <- function(
        fisher_df,
        p_value_col = "Fisher_p_value_fdr",
        annot_threshold = 0.05,
        annot_color = "red",
        gene_sym_col = "GeneName",
        do_not_highlight = NULL,
        keep_not_highlighted = TRUE) {
    stopifnot(is.null(do_not_highlight) || is.character(do_not_highlight) ||
        rlang::is_expression(do_not_highlight))
    if (is.null(do_not_highlight)) {
        below_threshold <- fisher_df |>
            dplyr::filter(.data[[p_value_col]] < annot_threshold)
        above_threshold <- fisher_df |>
            dplyr::filter(.data[[p_value_col]] >= annot_threshold)
    } else if (is.character(do_not_highlight)) {
        below_threshold <- fisher_df |>
            dplyr::filter(.data[[p_value_col]] < annot_threshold &
                !.data[[gene_sym_col]] %in% do_not_highlight)
        if (keep_not_highlighted) {
            above_threshold <- fisher_df |>
                dplyr::anti_join(below_threshold, by = gene_sym_col)
        } else {
            above_threshold <- fisher_df |>
                dplyr::filter(.data[[p_value_col]] >= annot_threshold)
        }
    } else if (rlang::is_formula(do_not_highlight)) {
        below_threshold <- fisher_df |>
            dplyr::filter(.data[[p_value_col]] < annot_threshold &
                !rlang::as_function(do_not_highlight)(
                    .data[[gene_sym_col]]))
        if (keep_not_highlighted) {
            above_threshold <- fisher_df |>
                dplyr::anti_join(below_threshold, by = gene_sym_col)
        } else {
            above_threshold <- fisher_df |>
                dplyr::filter(.data[[p_value_col]] >= annot_threshold)
        }
    } else {
        below_threshold <- fisher_df |>
            dplyr::filter(.data[[p_value_col]] < annot_threshold &
                (rlang::eval_tidy(do_not_highlight)))
        if (keep_not_highlighted) {
            above_threshold <- fisher_df |>
                dplyr::anti_join(below_threshold, by = gene_sym_col)
        } else {
            above_threshold <- fisher_df |>
                dplyr::filter(.data[[p_value_col]] >= annot_threshold)
        }
    }
    plot <- ggplot2::ggplot(
        above_threshold,
        ggplot2::aes(
            x = .data[["IS_per_kbGeneLen_1"]],
            y = .data[["IS_per_kbGeneLen_2"]],
            color = .data[[p_value_col]]
        )
    ) +
        ggplot2::geom_point(alpha = 0.65) +
        ggplot2::geom_point(data = below_threshold, color = annot_color) +
        ggplot2::theme_bw() +
        ggplot2::labs(
            x = "Gene frequency G1", y = "Gene frequency G2",
            color = "Fisher's test p-value"
        ) +
        ggrepel::geom_label_repel(
            data = below_threshold,
            ggplot2::aes(label = .data[[gene_sym_col]]),
            box.padding = ggplot2::unit(0.35, "lines"),
            point.padding = ggplot2::unit(0.3, "lines"),
            max.overlaps = Inf, color = annot_color,
            fill = ggplot2::alpha("white", alpha = 0.6)
        ) +
        ggplot2::scale_x_continuous() +
        ggplot2::scale_y_continuous()
    return(plot)
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
#' p <- HSC_population_plot(estimate$est, "PJ01")
#' p
HSC_population_plot <- function(
        estimates,
        project_name,
        timepoints = "Consecutive",
        models = "Mth Chao (LB)") {
    if (is.null(estimates)) {
        return(NULL)
    }
    ## Pre-filter
    df <- estimates |>
        dplyr::filter(
            .data$Timepoints %in% timepoints,
            .data$Model %in% models
        )
    p <- ggplot2::ggplot(
        data = df,
        ggplot2::aes(
            y = .data[["PopSize"]],
            x = .data[["TimePoint_to"]],
            color = .data[["SubjectID"]]
        ),
        na.rm = TRUE, se = TRUE
    ) +
        ggplot2::geom_point(alpha = .5) +
        ggplot2::geom_line(linewidth = 2, alpha = .7) +
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
#' @description
#' `r lifecycle::badge("stable")`
#' Alluvial plots allow the visualization of integration sites distribution
#' in different points in time in the same group.
#' This functionality requires the suggested package
#' [ggalluvial](https://corybrunson.github.io/ggalluvial/).
#'
#' @details
#' ## Input data frame
#' The input data frame must contain all the columns specified in the
#' arguments `group`, `plot_x`, `plot_y` and `alluvia`. The standard
#' input for this function is the data frame obtained via the
#' \link{compute_abundance} function.
#'
#' ## Plotting threshold on y
#' The plotting threshold on the quantification on the y axis has the
#' function to highlight only relevant information on the plot and reduce
#' computation time. The default value is 1, that acts on the default column
#' plotted on the y axis which contains a percentage value. This translates
#' in natural language roughly as "highlight with colors only those
#' integrations (alluvia) that at least in 1 point in time have an
#' abundance value >= 1 %". The remaining integrations will be plotted
#' as a unique layer in the column, colored as specified by the argument
#' `empty_space_color`.
#'
#' ## Customizing the plot
#' The returned plots are ggplot2 objects and can therefore further modified
#' as any other ggplot2 object. For example, if the user decides to change the
#' fill scale it is sufficient to do
#'
#' ```{r eval=FALSE}
#' plot +
#'   ggplot2::scale_fill_viridis_d(...) + # or any other discrete fill scale
#'   ggplot2::theme(...) # change theme options
#' ```
#' NOTE: if you requested the computation of the top ten abundant tables and
#' you want the colors to match you should re-compute them
#'
#' ## A note on strata ordering
#' Strata in each column are ordered first by time of appearance and secondly
#' in decreasing order of abundance (value of y). It means, for example,
#' that if the plot has 2 or more columns, in the second column, on top,
#' will appear first appear IS that appeared in the previous columns and then
#' all other IS, ordered in decreasing order of abundance.
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
#' @param empty_space_color Color of the empty portion of the bars (IS below
#' the threshold). Can be either a string of known colors, an hex code or
#' `NA_character` to set the space transparent. All color
#' specs accepted in ggplot2
#' are suitable here.
#' @param ... Additional arguments to pass on to \link{top_abund_tableGrob}
#'
#' @family Plotting functions
#' @importFrom rlang abort eval_tidy call2 inform .data fn_fmls_names dots_list
#' @importFrom dplyr group_by across group_keys everything pull group_split
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
integration_alluvial_plot <- function(
        x,
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
        empty_space_color = "grey90",
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
    groups_to_plot <- x |>
        dplyr::group_by(dplyr::across({{ group }}))
    group_names <- groups_to_plot |>
        dplyr::group_keys() |>
        tidyr::unite(col = "id", dplyr::everything()) |>
        dplyr::pull(.data$id)
    groups_to_plot <- groups_to_plot |>
        dplyr::group_split() |>
        purrr::set_names(group_names)

    # # Compute plots in parallel
    FUN <- function(
        group_df,
        plot_x,
        plot_y,
        alluvia,
        alluvia_plot_y_threshold,
        top_abundant_tbl,
        empty_space_color,
        other_params,
        progress) {
        cleaned <- .clean_data(
            group_df, plot_x, plot_y,
            alluvia, alluvia_plot_y_threshold
        )
        alluv_plot <- .alluvial_plot(
            cleaned, plot_x, plot_y,
            empty_space_color
        )
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
                            msg <- paste(
                                "Failed to produce top",
                                "abundance tables, skipping"
                            )
                            rlang::inform(msg)
                        }
                    )
                },
                error = function(cnd) {
                    rlang::inform(conditionMessage(cnd))
                    invokeRestart("missing")
                }
            )
            if (!is.null(progress)) {
                progress()
            }
            return(list(plot = alluv_plot, tables = summary_tbls))
        }
        if (!is.null(progress)) {
            progress()
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
    results <- .execute_map_job(
        data_list = groups_to_plot,
        fun_to_apply = FUN,
        fun_args = list(
            plot_x = plot_x,
            plot_y = plot_y,
            alluvia = alluvia,
            alluvia_plot_y_threshold = alluvia_plot_y_threshold,
            top_abundant_tbl = top_abundant_tbl,
            empty_space_color = empty_space_color,
            other_params = dot_args
        ),
        stop_on_error = FALSE
    )
    if (!all(purrr::map_lgl(results$err, is.null))) {
        errors_msg <- purrr::list_c(results$err)
        names(errors_msg) <- "x"
        errors_msg <- c(
            "Errors occurred during computation",
            errors_msg
        )
        rlang::inform(errors_msg)
    }
    return(results$res)
}

# Internal, used in integration_alluvial_plot to obtain the data frame
# with data to plot
#' @importFrom tidyr unite
#' @importFrom dplyr group_by across summarise n filter distinct pull rename
#' @importFrom rlang .data
.clean_data <- function(
        df,
        plot_x,
        plot_y,
        alluvia,
        alluvia_plot_y_threshold) {
    tbl <- if (length(alluvia) > 1) {
        df |>
            tidyr::unite(col = "alluvia_id", {{ alluvia }})
    } else {
        df |>
            dplyr::rename(alluvia_id = alluvia)
    }
    ## Getting counts by x
    counts_by_x <- tbl |>
        dplyr::group_by(dplyr::across({{ plot_x }})) |>
        dplyr::summarise(count = dplyr::n(), .groups = "drop")
    # Filtering alluvia to plot
    alluvia_to_plot <- tbl |>
        dplyr::filter(.data[[plot_y]] >= alluvia_plot_y_threshold[1]) |>
        dplyr::distinct(.data[["alluvia_id"]]) |>
        dplyr::pull(.data[["alluvia_id"]])
    # Modify ids - identifiers that are below threshold are converted to
    # NA and the quantities in y are aggregated
    tbl <- tbl |>
        dplyr::mutate(
            alluvia_id = dplyr::if_else(
                .data$alluvia_id %in% alluvia_to_plot,
                .data$alluvia_id, NA_character_
            )
        ) |>
        dplyr::group_by(dplyr::across(dplyr::all_of(c(plot_x, "alluvia_id")))) |>
        dplyr::summarise(!!plot_y := sum(.data[[plot_y]], na.rm = TRUE),
            .groups = "drop"
        )
    # Add counts
    tbl <- tbl |>
        dplyr::left_join(counts_by_x, by = plot_x)
    tbl
}

# Internal, used in integration_alluvial_plot to obtain the alluvial plots
# for a single group. NOTE: tbl must contain the column "alluvia_id" and
# "counts"
#' @importFrom ggplot2 ggplot geom_text scale_fill_viridis_d sym
#' @importFrom dplyr group_by summarise pull across all_of
#' @importFrom rlang .data
.alluvial_plot <- function(
        tbl, plot_x, plot_y,
        empty_space_color) {
    sums_y <- tbl |>
        dplyr::group_by(dplyr::across(dplyr::all_of(plot_x))) |>
        dplyr::summarise(sums = sum(.data[[plot_y]])) |>
        dplyr::pull(.data$sums)
    max_y <- max(sums_y)
    labels <- tbl |>
        dplyr::group_by(dplyr::across(dplyr::all_of(plot_x))) |>
        dplyr::summarise(count = .data$count[1], .groups = "drop")
    tbl <- tbl |>
        dplyr::group_by(dplyr::across(dplyr::all_of(plot_x))) |>
        dplyr::arrange(dplyr::desc(dplyr::across(dplyr::all_of(plot_y))),
            .by_group = TRUE
        ) |>
        dplyr::mutate(alluvia_id = forcats::as_factor(.data$alluvia_id)) |>
        dplyr::ungroup()
    alluv <- ggplot2::ggplot(
        tbl,
        ggplot2::aes(
            x = .data[[plot_x]],
            y = .data[[plot_y]],
            alluvium = .data[["alluvia_id"]]
        )
    ) +
        ggalluvial::stat_stratum(ggplot2::aes(stratum = .data[["alluvia_id"]]),
            na.rm = FALSE,
            fill = empty_space_color
        ) +
        ggalluvial::geom_alluvium(ggplot2::aes(fill = .data[["alluvia_id"]]),
            na.rm = FALSE,
            alpha = .75,
            aes.bind = "alluvia"
        )
    alluv <- alluv +
        ggplot2::scale_fill_viridis_d(na.value = NA_character_) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none") +
        ggplot2::geom_text(
            data = labels,
            ggplot2::aes(
                x = .data[[plot_x]],
                y = max_y + 5,
                label = .data[["count"]]
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
#' @param transform_by Either a function or a purrr-style lambda. This
#' function is applied to the column `by` before separating columns. If
#' `NULL` no function is applied. Useful to modify column order in final table.
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
#'
#' # with transform
#' grob <- top_abund_tableGrob(abund, transform_by = ~ as.numeric(.x))
top_abund_tableGrob <- function(
        df,
        id_cols = mandatory_IS_vars(),
        quant_col = "fragmentEstimate_sum_PercAbundance",
        by = "TimePoint",
        alluvial_plot = NULL,
        top_n = 10,
        tbl_cols = "GeneName",
        include_id_cols = FALSE,
        digits = 2,
        perc_symbol = TRUE,
        transform_by = NULL) {
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
        tibble::tibble(g$data[[2]]["alluvium"], g$data[[2]]["fill"]) |>
            dplyr::rename(alluvia_id = "alluvium")
    } else {
        NULL
    }
    if (!is.null(transform_by)) {
        df <- df |>
            dplyr::mutate(
                dplyr::across(
                    .cols = dplyr::all_of(by),
                    .fns = transform_by
                )
            )
    }
    top <- df |>
        tidyr::unite({{ id_cols }}, col = "alluvia_id") |>
        dplyr::select(dplyr::all_of(c(
            "alluvia_id", tbl_cols,
            quant_col, by
        ))) |>
        dplyr::group_by(dplyr::across({{ by }})) |>
        dplyr::group_modify(~ {
            .x |>
                dplyr::arrange(dplyr::desc(.data[[quant_col[1]]])) |>
                dplyr::slice_head(n = top_n)
        }) |>
        dplyr::ungroup() |>
        dplyr::mutate(font_col = "black") |>
        dplyr::rename(Ab = quant_col)
    if (perc_symbol == TRUE) {
        top <- top |>
            dplyr::mutate(Ab = paste(format(
                round(.data$Ab, digits = digits),
                nsmall = digits
            ), "%"))
    } else {
        top <- top |>
            dplyr::mutate(Ab = format(
                round(.data$Ab,
                    digits = digits
                ),
                nsmall = digits
            ))
    }
    if (!is.null(color_map)) {
        top <- top |>
            dplyr::left_join(color_map, by = c("alluvia_id")) |>
            dplyr::distinct()
    }
    if (include_id_cols == FALSE) {
        top <- top |>
            dplyr::select(-.data$alluvia_id)
    }
    distinct_x <- unique(top[[by]])
    sep_x <- function(x) {
        tmp <- if (is.na(x)) {
            top |>
                dplyr::filter(
                    is.na(.data[[by]])
                )
        } else {
            top |>
                dplyr::filter(
                    .data[[by]] == x
                )
        }
        tmp <- tmp |>
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
        tmp <- tmp |>
            dplyr::rename_with(.fn = ~ paste0(.x, " ", x))
    }
    tops_by_x <- purrr::map(distinct_x, sep_x) |> purrr::set_names(distinct_x)

    obtain_grobs <- function(df, x) {
        fill_var <- colnames(df)[stringr::str_detect(colnames(df), "fill")]
        col_var <- colnames(df)[stringr::str_detect(colnames(df), "font_col")]
        fill_vec <- if (!length(fill_var) == 0) {
            df |> dplyr::pull(fill_var)
        } else {
            NULL
        }
        col_vec <- df |> dplyr::pull(col_var)
        df_minus <- df |>
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
#' @description
#' `r lifecycle::badge("stable")`
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
#' @importFrom ggplot2 ggplot geom_raster scale_fill_gradientn geom_text
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
sharing_heatmap <- function(
        sharing_df,
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
        ggplot2::aes(
            x = .data[[show_on_x[1]]],
            y = .data[[show_on_y[1]]],
            fill = .data[[absolute_sharing_col[1]]],
            alpha = .data[[absolute_sharing_col[1]]],
            label = .data[[absolute_sharing_col[1]]]
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
        sharing_df_rounding <- sharing_df |>
            dplyr::mutate(dplyr::across(
                .cols = dplyr::all_of(rel_sharing_col),
                .fns = ~ round(.x, digits = 2)
            ))
        plot_rel_heat <- function(col, df) {
            plot <- ggplot2::ggplot(
                df,
                ggplot2::aes(
                    x = .data[[show_on_x[1]]],
                    y = .data[[show_on_y[1]]],
                    fill = .data[[col]],
                    alpha = .data[[col]]
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
                    ggplot2::geom_text(
                        ggplot2::aes(
                            alpha = 1,
                            label = scales::percent(
                                .data[[col]],
                                scale = 1,
                                accuracy = 0.1
                            )
                        ),
                        size = 3
                    )
            } else {
                plot <- plot +
                    ggplot2::geom_text(
                        ggplot2::aes(
                            alpha = 1,
                            label = .data[[col]]
                        ),
                        size = 3
                    )
            }
            plot
        }
        heatmap_rel <- purrr::map(
            unique(rel_sharing_col),
            ~ plot_rel_heat(.x, df = sharing_df_rounding)
        ) |>
            purrr::set_names(rel_sharing_col)
    }

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
#' @description
#' `r lifecycle::badge("stable")`
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
sharing_venn <- function(
        sharing_df,
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
    filtered_df <- sharing_df[row_range, ]
    if (nrow(filtered_df) == 0) {
        rlang::inform("Empty table, nothing to compute")
        return(NULL)
    }
    plot_fun <- if (euler) {
        eulerr::euler
    } else {
        eulerr::venn
    }
    transform_and_compute <- function(x, plot_fun) {
        as_matrix <- x |>
            tibble::column_to_rownames("int_id") |>
            as.matrix()
        rlang::exec(plot_fun, combinations = as_matrix)
    }
    fixed_tbls <- purrr::map(filtered_df$truth_tbl_venn,
        transform_and_compute,
        plot_fun = plot_fun
    )
    fixed_tbls
}

#' Trace a circos plot of genomic densities.
#'
#' @description
#' `r lifecycle::badge("stable")`
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
#' by_subj <- aggreg |>
#'     dplyr::group_by(.data$SubjectID) |>
#'     dplyr::group_split()
#' circos_genomic_density(by_subj,
#'     track_colors = c("navyblue", "gold"),
#'     grDevice = "default", track.height = 0.1
#' )
#' }
circos_genomic_density <- function(
        data,
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
        df |>
            dplyr::select(dplyr::all_of(mandatory_IS_vars())) |>
            dplyr::distinct() |>
            dplyr::mutate(
                chr = paste0("chr", .data$chr),
                end = .data$integration_locus
            ) |>
            dplyr::rename(start = "integration_locus") |>
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
        if (device == "pdf") {
            device_args <- device_args[!names(device_args) %in% c("file")]
            rlang::exec(.fn = device, file = file_path, !!!device_args)
        } else {
            device_args <- device_args[!names(device_args) %in% c("filename")]
            rlang::exec(.fn = device, filename = file_path, !!!device_args)
        }
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
