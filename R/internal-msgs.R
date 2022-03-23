### Convenience functions for errors and warnings ###

.missing_req_cols <- function(requested, missing) {
  err_msg <- c("Some of the required columns were not found",
               i = paste0("Requested: ",
                          paste0(requested, collapse = ", ")),
               x = paste0("Not found: ",
                          paste0(missing, collapse = ", ")))
  return(err_msg)
}

# Signals the user that file system alignment can't be performed because
# column containing the path to project folder is missing.
# USED IN:
# - .check_file_system_alignment
.af_missing_pathfolder_error <- function(folder_col) {
    c(paste("Column", folder_col, "is missing"),
        x = "Can't proceed with file system alignment",
        i = paste(
            "To import metadata without performing file system aligment",
            "set the parameter 'root' to NULL"
        )
    )
}

# Signals the user that the column
# 'Path' is missing (needed for matrix import).
# USED IN:
# - import_parallel_Vispa2Matrices
.af_missing_path_error <- function(col) {
    c(paste0("Column '", col, "' not found in the association file"),
        x = paste(
            "File system alignment is necessary for this step. Please",
            "re-import the association file with the alignment feature"
        )
    )
}

# Produces a mini-report to print after file reading only if verbose is active
# USED IN:
# - import_associaton_file
.summary_af_import_msg <- function(pars_prob, dates_prob, cols_prob, crit_na,
    checks) {
    c(
        "*** Association file import summary ***",
        i = paste(
            "For detailed report please set option",
            "'ISAnalytics.reports' to TRUE"
        ),
        paste0("Parsing problems detected: ", !is.null(pars_prob)),
        paste0("Date parsing problems: ", !is.null(dates_prob)),
        paste0("Column problems detected: ", !is.null(cols_prob)),
        paste0("NAs found in important columns: ", !is.null(crit_na)),
        paste0("File system alignment: ", checks)
    )
}

# Warns the user that the input file is in excel format and it is not
# reccommended for potential parsing issues.
# USED IN:
# - .read_af
.xls_file_warning <- function() {
    c("Warning: file in xls/xlsx format",
        i = paste(
            "The use of xls/xlsx is discouraged as it can lead to",
            "potential problems in parsing data. Use of *.tsv or *.csv",
            "is recommended."
        ),
        "Carefully review the data after importing!"
    )
}

# Signals that the association file must be imported and aligned with fs.
# USED IN:
# - import_Vispa2_stats
.af_not_imported_err <- function() {
    c("The association file must be a data frame",
        paste(
            "Import the association file via",
            "`import_association_file()` with file system alignment"
        ),
        i = "See `?import_association_file` and `?import_Vispa2_stats`"
    )
}

# Signals that the association file must be aligned with fs.
# USED IN:
# - import_Vispa2_stats
.af_not_aligned_err <- function() {
    c("The association file has been imported without file system alignment",
        paste(
            "Import the association file via",
            "`import_association_file()` with file system alignment"
        ),
        i = "See `?import_association_file` and `?import_Vispa2_stats`"
    )
}

# USED IN:
# - import_Vispa2_stats
.missing_needed_cols <- function(missing) {
    c("Some required columns are missing",
        i = paste("Missing columns:", paste0(missing, collapse = ", "))
    )
}

# Produces a mini-report to print after file reading only if verbose is active
# USED IN:
# - import_single_Vispa2Matrix
.summary_ism_import_msg <- function(annotated,
                                    dims,
                                    mode,
                                    sample_count) {
    c(
        "*** File info *** ",
        paste("--- Annotated:", annotated),
        paste("--- Dimensions:", paste0(dims, collapse = " x ")),
        paste("--- Read mode:", mode),
        paste("--- Sample count:", sample_count)
    )
}


# Warns the user that the file in input is compressed but the compression
# format is not supported by fread.
# USED IN:
# - import_single_Vispa2Matrix
# - .read_af
.unsupported_comp_format_inf <- function() {
    c(paste(
        "Warning: compression format not",
        "supported by fread"
    ),
    i = "File will be read using readr"
    )
}

# Error that signals the integration matrix to import is malformed, aka
# it does not contain mandatory variables or does not have a standard structure
# USED IN:
# - import_single_Vispa2Matrix
.malformed_ISmatrix_error <- function() {
    c("The input integration matrix seems to be malformed",
        x = "Non standard column structure detected",
        i = paste(
            "Matrix should contain either these columns: [",
            paste0(mandatory_IS_vars(), collapse = ", "),
            "], or",
            "'IS_genomicID'",
            "but not both or a combination of the two."
        )
    )
}

# General error, used in multiple functions: signals that the input is
# formally not considered an integration matrix
.non_ISM_error <- function() {
    paste(
        "One or more elements in x are not integration matrices.",
        "Aborting."
    )
}

# @keywords internal
.missing_complAmpID_error <- function() {
    paste(
        "The `CompleteAmplificationID` column is missing.",
        "The column is needed for this operation.",
        "Aborting."
    )
}

# USED IN:
# - remove_collisions
# - compute_near_integrations
.missing_cAmp_sub_msg <- function() {
    paste0(
        "The `CompleteAmplificationID` column is missing",
        " and it is needed for this functionality"
    )
}

# USED IN:
# - remove_collisions
# - compute_near_integrations
# - circos_genomic_density
.missing_mand_vars <- function() {
    paste(
        "Missing mandatory vars: ",
        paste0(mandatory_IS_vars(), collapse = ", ")
    )
}

# @keywords internal
.missing_num_cols_error <- function() {
    paste("No numeric columns found")
}

# @keywords internal
.non_quant_cols_msg <- function(x) {
    c(paste(
        "Found numeric columns that are not quantification values -",
        "these columns will be copied in all resulting matrices."
    ),
    i = paste0("Found: ", paste0(x, collapse = ", "))
    )
}

# @keywords internal
.non_quant_cols_error <- function() {
    paste(
        "No quantification values columns found. Did you set the function",
        "parameters correctly?"
    )
}

# @keywords internal
.max_val_col_warning <- function(x) {
    paste0("Column for max value `", x, "` not found in numeric columns.")
}

# @keywords internal
.using_val_col_warning <- function(x) {
    paste(c(
        .max_val_col_warning(x),
        "Using `Value` column as reference instead."
    ), collapse = "\n")
}

# @keywords internal
.max_val_stop_error <- function(x) {
    paste(c(
        .max_val_col_warning(x),
        "Did you set `max_value_column` parameter correctly?"
    ),
    collapse = "\n"
    )
}

# @keywords internal
.nas_introduced_msg <- function() {
    c("NAs were introduced while producing the data frame.",
        i = paste(
            "The possible cause for this is:",
            "some quantification matrices were not imported for all pools"
        )
    )
}

.no_collisions_msg <- function() {
    paste("No collisions found for the given matrix, nothing to do")
}


# @keywords internal
.threshold_err <- function() {
    paste("The parameter `threshold` must be numeric or integer")
}

# @keywords internal
.comparators_err <- function() {
    paste(c(
        "The parameter `comparators` must be a character vector",
        "and in the allowed range: ", c("<", ">", "==", "!=", ">=", "<=")
    ),
    collapse = " "
    )
}

# @keywords internal
.names_list_param_err <- function() {
    paste(
        "If provided as lists, parameters must have",
        "names equal to names in x"
    )
}

# @keywords internal
.diff_leng_args <- function() {
    paste0("Parameters `threshold`, `cols_to_compare` and `comparators`",
        "have different lengths",
        collapse = " "
    )
}

# Error message: notifies the user that one or more column names specified
# as parameters were not found in the data frame column names
# USED IN:
# - compute_abundance
# - top_integrations
# - outliers_by_pool_fragments
# - sample_statistics
# - compute_near_integrations
# - cumulative_is
.missing_user_cols_error <- function(missing_cols) {
    c(paste(
        "Some or all of the input column names were not found",
        "in the data frame"
    ),
    i = paste(
        "Columns missing:",
        paste0(missing_cols, collapse = ", ")
    )
    )
}

# Error message: notifies the user that one or more column names specified
# as parameters were not found in metadata column names
# USED IN:
# - sample_statistics
.missing_user_cols_meta_error <- function(missing_cols) {
    c(paste(
        "Some or all of the input column names were not found",
        "in metadata"
    ),
    i = paste(
        "Columns missing:",
        paste0(missing_cols, collapse = ", ")
    )
    )
}

# @keywords internal
.missing_user_cols_list_error <- function(set, name) {
    paste0(c(
        "Columns ",
        paste0(set, collapse = ", "),
        " not found in x$", name
    ), collapse = "")
}

# Error message: notifies the user that one or more column names specified
# as parameters are not numeric
# USED IN:
# - compute_abundance
# - sample_statistics
.non_num_user_cols_error <- function(non_num) {
    paste(
        paste(
            "Some or all of the input column names are not numeric",
            "or integer in the data frame"
        ),
        "Non-numeric columns:",
        paste0(non_num, collapse = ", "),
        sep = "\t"
    )
}

# Error: the function expects all element in the list to be either
# functions or formulas
# USED IN:
# - sample_statistics
.non_function_elem_error <- function() {
    c("Wrong element type",
        x = paste(
            "The function parameter should contain a list",
            "of either functions or formulas."
        ),
        i = "See ?sample_statistics for details"
    )
}

# USED IN:
# - sample_statistics
.inform_skip_count_is <- function() {
    c("Mandatory IS vars are missing, skipping count IS")
}

# @keywords internal
.list_params_err <- function() {
    paste(
        "Provided parameters `threshold`, `comparators`",
        "or `cols_to_compare` are lists when input is a",
        "data frame. Please see ?threshold_filter for more",
        "info."
    )
}

# @keywords internal
.unnamed_list_err <- function() {
    paste(
        "Some parameters were provided as named lists",
        "but x is unnamed list. Give names to x or",
        "change the parameters, see ?threshold_filter",
        "for details"
    )
}

# @keywords internal
.param_names_not_in_list_err <- function() {
    paste(
        "Named lists in parameters must have the same",
        "names as x, see ?threshold_filter",
        "for details"
    )
}

# @keywords internal
.param_names_not_equal_err <- function() {
    paste(
        "Some list parameters between `threshold`,",
        "`cols_to_compare` or `comparators` miss elements.",
        "See ?threshold_filter for details"
    )
}

#---- USED IN : CIS_grubbs ----
.missing_annot <- function() {
    paste(
        "Annotation columns are missing but are required for",
        "the function correct execution"
    )
}

.non_standard_annotation_structure <- function() {
    paste(
        "The genomic annotation file must have the standard UCSC format,",
        "see ?CIS_grubbs for details"
    )
}

#---- USED IN : cumulative_count_union ----
.agg_with_null_meta_err <- function() {
    paste(
        "Matrix aggregation can't be performed without the ",
        "association file, please specify one in the `metadata` parameter"
    )
}

.key_without_tp_err <- function() {
    paste("The sample key must contain the time point column")
}

.key_not_found <- function() {
    paste("One or more columns in the sample keys were not found in x")
}

.not_min_key_err <- function(missing) {
    c("The aggregation key must contain the minimal required key",
        x = paste(
            "Missing columns:",
            paste0(missing, collapse = ", ")
        )
    )
}

.agg_key_not_found_err <- function(df, key) {
    c(paste("The aggregation key was not found in", df),
        i = paste("Aggregation key used:", paste0(key, collapse = ", "))
    )
}

.meta_not_agg_err <- function() {
    c("Metadata appears to not be aggregated by the provided aggregation key",
        i = paste("See `?aggregate_metadata`")
    )
}

# Error message displayed for suggestion packages that are not installed
# but required by the called function
.missing_pkg_error <- function(pkg, lib = "CRAN") {
    inst_sugg <- if (lib == "CRAN") {
        paste0('To install: `install.packages("', pkg, '")`')
    } else if (lib == "BIOC") {
        paste0('To install: `BiocManager::install("', pkg, '")`')
    }
    c("Missing package",
        x = paste0("Package '", pkg, "' is required for this functionality."),
        i = inst_sugg
    )
}

# Error message displayed for functions that require names of a list to be
# quantification types.
# USED IN:
# - remove_collisions
.quantifications_names_err <- function(non_quant) {
    c("Unknown quantifications",
        x = paste0(
            "Names of the input list must be quantification ",
            "types in `quantification_types()`"
        ),
        i = paste0("Unknown quantifications: ", paste0(non_quant, collapse = ", "))
    )
}

# Error message displayed when sequence count matrix is needed but not detected
# IN A LIST.
# USED IN:
# - remove_collisions
.seqCount_df_err <- function() {
    c(x = paste(
        "Sequence count data frame is required for collision removal",
        "but none was detected in x"
    ))
}

# Error message displayed when sequence count column is not detected in the
# input matrix
# USED IN:
# - remove_collisions
.seqCount_col_err <- function() {
    c("Sequence count column not found in the input data frame",
        i = "Check the `quant_cols` parameter. Did you set it correctly?"
    )
}

# Error message for missing needed AF columns
# USED IN:
# - remove_collisions
.missing_af_needed_cols <- function(missing) {
    c(
        "Missing needed columns in the association file: ",
        x = paste(
            "Columns missing:",
            paste0(missing, collapse = ", ")
        )
    )
}

# USED IN:
# - remove_collisions
.na_in_date_col <- function() {
    c(paste(
        "Selected date column contains NA values"
    ))
}

# Warns the user that some samples contained in the matrices are not found
# in metadata and are going to be excluded from computation.
# USED IN:
# - remove_collisions
.missing_af_samples_msg <- function(n_missing) {
    c("Missing samples from association file",
        i = paste(
            "Could not find some samples contained in",
            "matrices in the association file. For more detail",
            "please consult the report. Missing samples",
            "will be removed from the matrix."
        ),
        paste0("-- Number of samples missing: ", n_missing)
    )
}

# Warns the user that some samples are present in the association file but
# not in matrix.
# USED IN:
# - remove_collisions
.additional_ad_samples_msg <- function() {
    c("Additional samples found in the association file",
        i = paste(
            "Some samples were found in the association file but not",
            "in the input matrix. For more detail consult the report."
        )
    )
}

# Error displayed if no root argument was set but import of ISS was requested
# USED IN:
# - import_associaton_file
.no_stats_import_err <- function() {
    c(paste(
        "Can't import Vispa2 stats files without",
        "file system alignment. Provide the appropriate",
        "root."
    ))
}

.no_stat_files_imported <- function() {
    c("Warning: no VISPA2 stats files imported")
}

# Message that warns the user no matrices can be imported because all paths
# are NA
# USED IN:
# - import_parallel_Vispa2Matrices_interactive/auto
.af_empty_msg <- function() {
    c("Nothing to import",
        i = paste(
            "All paths are NA, no matrices to import.",
            "Check your association file."
        )
    )
}

# USED IN:
# - is_sharing
.non_df_input_err <- function() {
    c("Non-data frame found",
        x = paste("Input list contains objects that are not data.frame")
    )
}

# USED IN:
# - is_sharing
.unnamed_keys_warn <- function() {
    c("Unnamed keys",
        i = paste("Group keys should be named, using default names")
    )
}

# USED IN:
# - is_sharing
.keys_not_char_err <- function() {
    c("Keys need to be character vectors")
}

# USED IN:
# - is_sharing
.no_data_supp <- function() {
    c("No data supplied",
        x = paste("You must provide at least 1 data frame to the function")
    )
}

# USED IN:
# - cumulative_is
.only_zero_tp <- function() {
    c("Data frame empty after filtering",
        i = paste(
            "Data frame has 0 rows after filtering",
            "out zero-timepoints. Nothing to do."
        )
    )
}
