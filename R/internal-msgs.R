### Convenience functions for errors and warnings ###

.af_correctness_warning <- function() {
    paste(
        "Columns found in the imported association file differ from",
        "the default values in 'association_file_columns()' ",
        "problems may arise in some functions"
    )
}

.af_missing_path_warning <- function(alignment) {
    if (alignment) {
        paste(
            "Column 'PathToFolderProjectID' is missing, can't proceed",
            "with file system alignment"
        )
    } else {
        paste("Column 'PathToFolderProjectID' is missing")
    }
}

.af_missing_path_error <- function() {
    paste(
        "Column 'Path' not found in the association file,",
        "file system alignment is necessary for this step. Please",
        "re-import the association file with the alignment feature"
    )
}

.widgets_error <- function() {
    paste("Unable to produce widget report, skipping this step")
}

.widgets_save_error <- function() {
    paste("Unable to save widget to file, skipping this step")
}

# @keywords internal
.malformed_ISmatrix_warning <- function() {
    paste(c(
        "Mandatory integration matrix variables, ", mandatory_IS_vars(),
        ", were not detected"
    ), collapse = " ")
}
# @keywords internal
.non_ISM_error <- function() {
    paste(
        "One or more elements in x are not integration matrices.",
        "Aborting."
    )
}

# @keywords internal
.missing_value_col_error <- function() {
    paste(
        "The value column is missing or it contains non-numeric data.",
        "The column is needed for this operation.",
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

.warning_update_after_alignment <- function(root) {
    paste0("One or more projects were not found in the file ",
        "system starting from ", root, ", please check your ",
        "association file for errors and/or your file system.",
        "Until you re-import the association file these ",
        "missing files will be ignored.",
        collapse = ""
    )
}

# @keywords internal
.quant_types_error <- function() {
    paste(
        "The list names must be quantification types",
        ", see quantification_types() for reference"
    )
}

# @keywords internal
.missing_num_cols_error <- function() {
    paste("No numeric columns found")
}

# @keywords internal
.non_quant_cols_msg <- function(x) {
    paste(c(
        "Found numeric columns that are not quantification values:",
        "these columns will be copied in all resulting matrices.",
        "Found: ", x
    ), collapse = "\n")
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
    paste("NAs were introduced while producing the data frame.",
        "The possible cause for this is:",
        "some quantification matrices were not imported for all pools",
        sep = "\n"
    )
}

#---- USED IN : remove_collisions ----
.no_collisions_msg <- function() {
    paste("No collisions found for the given matrix, nothing to do")
}

#---- USED IN : threshold_filter ----

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

# @keywords internal
.missing_user_cols_error <- function() {
    paste(
        "Some or all of the column names specified were not found",
        "in the data frame"
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

# @keywords internal
.non_num_user_cols_error <- function() {
    paste(
        "Some or all of the column names specified are not numeric",
        "or integer columns in the data frame"
    )
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
