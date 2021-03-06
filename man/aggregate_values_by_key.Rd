% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregate-functions.R
\name{aggregate_values_by_key}
\alias{aggregate_values_by_key}
\title{Aggregates matrices values based on specified key.}
\usage{
aggregate_values_by_key(
  x,
  association_file,
  value_cols = "Value",
  key = c("SubjectID", "CellMarker", "Tissue", "TimePoint"),
  lambda = list(sum = ~sum(.x, na.rm = TRUE)),
  group = c(mandatory_IS_vars(), annotation_IS_vars()),
  join_af_by = "CompleteAmplificationID"
)
}
\arguments{
\item{x}{A single integration matrix (tibble) or a list of imported
integration matrices (tibble)}

\item{association_file}{The imported association file}

\item{value_cols}{A character vector containing the names of the
columns to apply the given lambdas. Must be numeric or integer
columns.}

\item{key}{A string or a character vector with column names of the
association file to take as key}

\item{lambda}{A named list of functions or purrr-style lambdas.
See details section.}

\item{group}{Other variables to include in the grouping besides \code{key},
can be set to NULL}

\item{join_af_by}{A character vector representing the joining key
between the matrix and the metadata. Useful to re-aggregate already
aggregated matrices.}
}
\value{
A list of tibbles or a single tibble aggregated according to
the specified arguments
}
\description{
\lifecycle{maturing}
Performs aggregation on values contained in the integration matrices based
on the key and the specified lambda. For more details on how to use this
function:
\code{vignette("Working with aggregate functions", package = "ISAnalytics")}
}
\details{
\subsection{Setting the lambda parameter}{

The lambda parameter should always contain a named list of either
functions or purrr-style lambdas.
It is also possible to specify the namespace of the function in both
ways, for example:\if{html}{\out{<div class="r">}}\preformatted{lambda = list(sum = sum, desc = psych::describe)
}\if{html}{\out{</div>}}

Using purrr-style lambdas allows to specify arguments for the functions,
keeping in mind that the first parameter should always be \code{.x}:\if{html}{\out{<div class="r">}}\preformatted{lambda = list(sum = ~sum(.x, na.rm = TRUE))
}\if{html}{\out{</div>}}

It is also possible to use custom user-defined functions, keeping in
mind that the symbol will be evaluated in the calling environment,
for example if the function is called in the global environment
and lambda contains "foo" as a function, "foo" will be evaluated in
the global environment.\if{html}{\out{<div class="r">}}\preformatted{foo <- function(x) \{
  sum(x)
\}

lambda = list(sum = ~sum(.x, na.rm = TRUE), foo = foo)

# Or with lambda notation
lambda = list(sum = ~sum(.x, na.rm = TRUE), foo = ~foo(.x))
}\if{html}{\out{</div>}}
}

\subsection{Constraints on aggregation functions}{

Functions passed in the lambda parameters must respect a few constraints
to properly work and it's the user responsibility to ensure this.
\itemize{
\item Functions have to accept as input a numeric or integer vector
\item Function should return a single value or a list/data frame:
if a list or a data frame is returned as a result, all the columns
will be added to the final data frame.
}
}
}
\examples{
op <- options("ISAnalytics.widgets" = FALSE, "ISAnalytics.verbose" = FALSE)
path_AF <- system.file("extdata", "ex_association_file.tsv",
    package = "ISAnalytics"
)
root_correct <- system.file("extdata", "fs.zip", package = "ISAnalytics")
root_correct <- unzip_file_system(root_correct, "fs")
association_file <- import_association_file(path_AF, root_correct,
    dates_format = "dmy"
)
matrices <- import_parallel_Vispa2Matrices_auto(
    association_file = association_file, root = NULL,
    quantification_type = c("fragmentEstimate", "seqCount"),
    matrix_type = "annotated", workers = 2, matching_opt = "ANY"
)
agg <- aggregate_values_by_key(
    x = matrices,
    association_file = association_file,
    value_cols = c("fragmentEstimate", "seqCount")
)
options(op)
}
\seealso{
Other Aggregate functions: 
\code{\link{aggregate_metadata}()},
\code{\link{default_meta_agg}()}
}
\concept{Aggregate functions}
