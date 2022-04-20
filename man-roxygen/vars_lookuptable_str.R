#' @details
#' The user can supply specifications in the form of a named vector or a
#' data frame.
#'
#' ## Named vector
#' When using a named vector, names should be the names of the columns,
#' values should be the type associated with each column in the form
#' of a string. The vector gets automatically converted into a data frame
#' with the right format (default values for the columns `transform` and
#' `flag` are `NULL` and `required` respectively). Use of this method is
#' however discouraged: data frame inputs are preferred since they offer more
#' control.
#'
#' ## Look-up table structure
#' The look-up table for dynamic vars should always follow this structure:
#' names                | types  | transform          | flag   | tag
#' ---------------------|--------|--------------------|--------|-----
#' `<name of the column>` | `<type>` | `<a lambda or NULL>` | `<flag>` | `<tag>`
#'
#' where
#' * `names` contains the name of the column as a character
#' * `types` contains the type of the column. Type should be expressed as a
#' string and should be in one of the allowed types
#' * `char` for character (strings)
#' * `int` for integers
#' * `logi` for logical values (TRUE / FALSE)
#' * `numeric` for numeric values
#' * `factor` for factors
#' * `date` for generic date format - note that functions that
#' need to read and parse files will try to guess the format and parsing
#' may fail
#' * One of the accepted date/datetime formats by `lubridate`,
#' you can use `ISAnalytics::date_formats()` to view the accepted formats
#' * `transform`: a purrr-style lambda that is applied immediately after
#' importing.
#' This is useful to operate simple transformations like removing unwanted
#' characters or rounding to a certain precision. Please note that these lambdas
#' need to be functions that accept a vector as input and only operate a
#' transformation, aka they output a vector of the same length as the
#' input. For more complicated applications that may require the value of other
#' columns, appropriate functions should be manually applied post-import.
#' * `flag`: as of now, it should be set either to `required` or `optional` -
#'  some functions internally check for only required tags presence and if those
#' are missing from inputs they fail, signaling failure to the user
#' * `tag`: a specific tag expressed as a string
#'
#' ## Column types:
#' Type should be expressed as a
#' string and should be in one of the allowed types
#' * `char` for character (strings)
#' * `int` for integers
#' * `logi` for logical values (TRUE / FALSE)
#' * `numeric` for numeric values
#' * `factor` for factors
#' * `date` for generic date format - note that functions that
#' need to read and parse files will try to guess the format and parsing
#' may fail
#' * One of the accepted date/datetime formats by `lubridate`,
#' you can use `ISAnalytics::date_formats()` to view the accepted formats
