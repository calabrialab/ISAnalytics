#------------------------------------------------------------------------------#
# Population size estimate
#------------------------------------------------------------------------------#

#' Hematopoietic stem cells population size estimate.
#'
#' @description \lifecycle{experimental}
#' Hematopoietic stem cells population size estimate with capture-recapture
#' models.
#'
#' @details
#' # Input formats
#' Both `x` and `metadata` should be supplied to the function in aggregated
#' format (ideally through the use of \code{\link{aggregate_metadata}}
#' and \code{\link{aggregate_values_by_key}}).
#' Note that the `aggregation_key`, aka the vector of column names used for
#' aggregation, must contain at least the columns SubjectID, CellMarker,
#' Tissue and a time point column (the user can specify the name of the
#' column in the argument `timepoint_column`).
#'
#' # On time points
#' If `stable_timepoints` is a vector with length > 1, the function will look
#' for the first available stable time point and slice the data from that
#' time point onward. If NULL is supplied instead, it means there are no
#' stable time points available. Note that 0 time points are ALWAYS discarded.
#' Also, to be included in the analysis, a group must have at least 2
#' distinct non-zero time points.
#'
#' @param x An aggregated integration matrix. See details.
#' @param metadata An aggregated association file. See details.
#' @param stable_timepoints A numeric vector or NULL if there are no
#' stable time points.
#' @param aggregation_key A character vector indicating the key used for
#' aggregating x and metadata. Note that x and metadata should always be
#' aggregated with the same key.
#' @param blood_lineages A data frame containing information on the blood
#' lineages. Users can supply their own, provided the columns `CellMarker` and
#' `CellType` are present.
#' @param timepoint_column What is the name of the time point column to use?
#' Note that this column must be present in the key.
#' @param seqCount_column What is the name of the column in x holding the
#' values of sequence count quantification?
#' @param seqCount_threshold A single numeric value. After re-aggregating `x`,
#' rows with a value greater or equal will be kept, the others will be
#' discarded.
#' @param nIS_threshold A single numeric value. If a group (row) in the
#' metadata data frame has a count of distinct integration sites strictly
#' greater than this number it will be kept, otherwise discarded.
#' @param cell_type The cell types to include in the models. Note that
#' the matching is case-insensitive.
#' @param tissue_type The tissue types to include in the models. Note that
#' the matching is case-insensitive.
#'
#' @return A data frame with the results of the estimates
#' @family Population estimates
#'
#' @import dplyr
#' @import BiocParallel
#' @importFrom rlang abort inform
#' @importFrom stringr str_to_upper str_detect
#' @importFrom purrr reduce
#'
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
#' estimate <- HSC_population_size_estimate(x = agg,
#' metadata = aggregated_meta,
#' stable_timepoints = NULL)
#' options(op)
HSC_population_size_estimate <- function(x,
    metadata,
    stable_timepoints = NULL,
    aggregation_key = c("SubjectID", "CellMarker", "Tissue", "TimePoint"),
    blood_lineages = blood_lineages_default(),
    timepoint_column = "TimePoint",
    seqCount_column = "seqCount_sum",
    seqCount_threshold = 3,
    nIS_threshold = 5,
    cell_type = "MYELOID",
    tissue_type = "PB") {
    # Param check
    ## Basic checks on types
    stopifnot(is.data.frame(x))
    stopifnot(is.data.frame(metadata))
    stopifnot(is.null(stable_timepoints) | is.numeric(stable_timepoints))
    stopifnot(is.character(aggregation_key))
    stopifnot(is.data.frame(blood_lineages))
    stopifnot(is.character(timepoint_column) & length(timepoint_column) == 1)
    stopifnot(is.character(seqCount_column) & length(seqCount_column) == 1)
    stopifnot(is.numeric(seqCount_threshold))
    seqCount_threshold <- seqCount_threshold[1]
    stopifnot(is.numeric(nIS_threshold))
    nIS_threshold <- nIS_threshold[1]
    stopifnot(is.character(cell_type))
    stopifnot(is.character(tissue_type))
    ## Convert cell_type and tissue_type to uppercase (for case insensitivity)
    cell_type <- stringr::str_to_upper(cell_type)
    tissue_type <- stringr::str_to_upper(tissue_type)
    ## Assumptions on aggregation key
    minimum_key <- c("SubjectID", "CellMarker", "Tissue", timepoint_column)
    if (!all(minimum_key %in% aggregation_key)) {
        rlang::abort(.not_min_key_err(
            minimum_key[!minimum_key %in% aggregation_key]
        ))
    }
    ## Aggregation key must be found in both data and meta
    if (!all(aggregation_key %in% colnames(x))) {
        rlang::abort(.agg_key_not_found_err("x", aggregation_key))
    }

    if (!all(aggregation_key %in% colnames(metadata))) {
        rlang::abort(.agg_key_not_found_err("metadata", aggregation_key))
    }
    ## Check actual aggregation
    distinct_agg_groups <- dplyr::distinct(
        metadata,
        dplyr::across(dplyr::all_of(aggregation_key))
    )
    meta_is_aggregated <- if (nrow(metadata) == nrow(distinct_agg_groups)) {
        TRUE
    } else {
        FALSE
    }
    if (!meta_is_aggregated) {
        rlang::abort(.meta_not_agg_err())
    }
    ## Check seqCount col in x
    if (!seqCount_column %in% colnames(x)) {
        rlang::abort(c("Sequence count column not found in x"))
    }
    ## Reorder stable timepoints
    if (is.null(stable_timepoints)) {
      stable_timepoints <- numeric(0)
    }
    stable_timepoints <- sort(stable_timepoints)
    ## Check presence of NumIS column
    if (!"NumIS" %in% colnames(metadata)) {
        if (getOption("ISAnalytics.verbose") == TRUE) {
            rlang::inform(c("Calculating number of IS for each group..."))
        }
        numIs <- x %>%
          dplyr::left_join(metadata, by = dplyr::all_of(aggregation_key)) %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(aggregation_key))) %>%
          dplyr::distinct(dplyr::across(dplyr::all_of(mandatory_IS_vars()))) %>%
          dplyr::count(name = "NumIS")
        metadata <- metadata %>%
          dplyr::left_join(numIs, by = dplyr::all_of(aggregation_key)) %>%
          dplyr::distinct()
    }
    ## Check blood lineages
    if (!all(c("CellMarker", "CellType") %in% colnames(blood_lineages))) {
        rlang::abort(c(paste(
            "The blood lineages table must contain at least",
            "the columns `CellMarker` and `CellType`"
        )))
    }
    # --- METADATA
    ### Join meta with blood lineages
    metadata <- metadata %>%
        dplyr::filter(.data$NumIS > nIS_threshold) %>%
        dplyr::left_join(blood_lineages, by = "CellMarker")
    # --- SPLIT THE INPUT AGGREGATED MATRIX BY SubjectID
    x_subj_split <- x %>%
        dplyr::group_by(.data$SubjectID) %>%
        dplyr::group_split()
    #### Process in parallel
    annotation_cols <- if (.is_annotated(x)) {
        annotation_IS_vars()
    } else {
        NULL
    }
    if (.Platform$OS.type == "windows") {
        p <- BiocParallel::SnowParam(
            stop.on.error = FALSE,
            tasks = length(x_subj_split),
            progressbar = getOption("ISAnalytics.verbose"),
            exportglobals = FALSE
        )
    } else {
        p <- BiocParallel::MulticoreParam(
            stop.on.error = FALSE,
            tasks = length(x_subj_split),
            progressbar = getOption("ISAnalytics.verbose"),
            exportglobals = FALSE
        )
    }
    FUN <- function(df, metadata) {
        # --- RE-AGGREGATION
        re_agg <- aggregate_values_by_key(
            x = df,
            association_file = metadata,
            value_cols = seqCount_column,
            key = c("CellType", "Tissue", timepoint_column),
            join_af_by = aggregation_key
        )
        re_agg <- re_agg %>%
            dplyr::filter(!is.na(.data$CellType))
        seqCount_column <- colnames(re_agg)[stringr::str_detect(
            colnames(re_agg),
            seqCount_column
        )]
        ### Keep only sequence count value greater or equal to
        ### seqCount_threshold
        ### (for each is)
        per_int_sums_low <- re_agg %>%
            dplyr::group_by(dplyr::across(
              dplyr::all_of(mandatory_IS_vars()))) %>%
            dplyr::summarise(
                sum = sum(.data[[seqCount_column]]),
                .groups = "drop"
            ) %>%
            dplyr::filter(.data$sum >= seqCount_threshold)
        re_agg <- re_agg %>%
            dplyr::filter(
                .data$chr %in% per_int_sums_low$chr,
                .data$integration_locus %in% per_int_sums_low$integration_locus,
                .data$strand %in% per_int_sums_low$strand
            )
        ### Keep values whose cell type is in cell_type, tissue
        ### in tissue_type and
        ### have timepoint greater than zero
        patient_slice <- re_agg %>%
            dplyr::filter(
                stringr::str_to_upper(.data$CellType) %in% cell_type,
                stringr::str_to_upper(.data$Tissue) %in% tissue_type,
                as.numeric(.data[[timepoint_column]]) > 0
            )
        ### If selected patient does not have at least 3 distinct timepoints
        ### simply return (do not consider for population estimate)
        if (length(unique(patient_slice[[timepoint_column]])) <= 2) {
            return(NULL)
        }
        estimate <- .estimate_pop(
            df = patient_slice,
            seqCount_column = seqCount_column,
            timepoint_column = timepoint_column,
            annotation_cols = annotation_cols,
            subject = df$SubjectID[1],
            stable_timepoints = stable_timepoints
        )
    }

    population_size <- BiocParallel::bplapply(x_subj_split,
        FUN = FUN,
        metadata = metadata,
        BPPARAM = p
    )
    BiocParallel::bpstop(p)
    population_size <- purrr::reduce(population_size, function(x, y) {
      if (is.null(x) & is.null(y)) {
        return(NULL)
      } else {
        return(dplyr::bind_rows(x, y))
      }
    })
    if (!is.null(population_size)) {
      population_size <- population_size %>%
        dplyr::mutate(PopSize = round(.data$abundance - .data$stderr))
    }
    population_size
}


#---- Default blood lineages --------------------------------------------------#
#' Default blood lineages info
#'
#' A default table with info relative to different blood lineages associated
#' with cell markers that can be supplied as a parameter to
#'  \code{\link{HSC_population_size_estimate}}
#'
#' @return A data frame
#' @importFrom tibble tribble
#' @export
#'
#' @examples
#' blood_lineages_default()
blood_lineages_default <- function() {
    tibble::tribble(
        ~CellMarker, ~Keywords, ~CellType, ~HematoLineage, ~SuperGroup,
        ~LineageByPurity,
        "CD13", "MYELO", "Myeloid", "Myeloid", "CD13", "Myeloid",
        "CD14", "MYELO", "Myeloid", "Myeloid", "CD14", "Myeloid",
        "CD15", "MYELO", "Myeloid", "Myeloid", "CD15", "Myeloid",
        "CD19", "B", "B", "Lymphoid", "CD19", "Lymphoid",
        "CD3", "T", "T", "Lymphoid", "CD3", "Lymphoid",
        "CD34", "CD34", "CD34", "CD34", "CD34", "CD34",
        "CD34Molmed", "CD34", "Other", "Other", "Other", "Other",
        "CD34NEG", "CD34", "Other", "Other", "Other", "Other",
        "CD36", "TE", "Erythroid", "Erythroid", "CD36", "Erythroid",
        "CD4", "T", "T", "Lymphoid", "CD3", "Lymphoid",
        "CD56", "NK", "NK", "Lymphoid", "CD56", "Other",
        "CD61", "MEGACARIO", "Megacario", "Megacario", "CD61", "Megacario",
        "CD8", "T", "T", "Lymphoid", "CD3", "Lymphoid",
        "CellLine", "CellLine", "Other", "Other", "Other", "Other",
        "CFC", "CFC", "Other", "Other", "Other", "Other",
        "CFC-BFUE", "CFC", "Other", "Other", "Other", "Other",
        "CFC-BULK", "CFC", "Other", "Other", "Other", "Other",
        "CFC-CFUGM", "CFC", "Other", "Other", "Other", "Other",
        "CFCPOOL", "CFC", "Other", "Other", "Other", "Other",
        "CMP", "CD34", "CD34", "CD34", "CD34", "CD34",
        "GLY", "TE", "Erythroid", "Erythroid", "GLY", "Erythroid",
        "GLYA", "TE", "Erythroid", "Erythroid", "GLY", "Erythroid",
        "GMP", "CD34", "CD34", "CD34", "CD34", "CD34",
        "Granulo", "Granulo", "Granulo", "Other", "Other", "Other",
        "H2O", "H2O", "Other", "Other", "Other", "Other",
        "HSC", "CD34", "CD34", "CD34", "CD34", "CD34",
        "LCE", "LCE", "Other", "Other", "Other", "Other",
        "MEP", "CD34", "CD34", "CD34", "CD34", "CD34",
        "MLP", "CD34", "CD34", "CD34", "CD34", "CD34",
        "MNC", "MNC", "Other", "Other", "Other", "Other",
        "Molmed-CD34", "CD34", "Other", "Other", "Other", "Other",
        "MPP", "CD34", "CD34", "CD34", "CD34", "CD34",
        "NEGCD34", "CD34", "Other", "Other", "Other", "Other",
        "NONE", "NONE", "Other", "Other", "Other", "Other",
        "PBMC", "Whole", "Other", "Other", "Other", "Other",
        "PHA", "CellLine", "Other", "Other", "Other", "Other",
        "Plasma", "PLASMA", "Plasma", "Plasma", "Plasma", "Plasma",
        "PreB-NK", "NK", "NK", "Lymphoid", "CD56", "Lymphoid",
        "PreBNK", "NK", "NK", "Lymphoid", "CD56", "Lymphoid",
        "UTR", "UTR", "Other", "Other", "Other", "Other",
        "WBM", "BM", "Other", "Other", "Other", "Other",
        "Whole", "Whole", "Other", "Other", "Other", "Other",
        "WPB", "PB", "Other", "Other", "Other", "Other",
        "CFCBFUE", "CFC", "Other", "Other", "Other", "Other",
        "CFCBULK", "CFC", "Other", "Other", "Other", "Other",
        "CFCCFUGM", "CFC", "Other", "Other", "Other", "Other",
        "Ctlin", "Ctlin", "Other", "Other", "Other", "Other"
    )
}
