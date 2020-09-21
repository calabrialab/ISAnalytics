### Convenience functions for errors and warnings ###

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
        "The `Value` column is missing or it contains non-numeric data.",
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
