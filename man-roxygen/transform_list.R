#' @param transformations Either `NULL` or a named list of purrr-style lambdas
#' where names
#' are column names the function should be applied to.
#' @details
#' ## Transformations
#' Lambdas provided in input in the `transformations` argument,
#' must be transformations, aka functions that take
#' in input a vector and return a vector of the same length as the input.
#'
#' If the transformation list contains column names that are not present
#' in the data frame, they are simply ignored.
#'
#' @seealso \code{\link{transform_columns}}
