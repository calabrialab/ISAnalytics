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
            dplyr::filter(stringr::str_to_lower(.data$GeneName) %in%
                stringr::str_to_lower(highlight_genes))
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
