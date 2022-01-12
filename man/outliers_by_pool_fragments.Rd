% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/raw-reads-filtering.R
\name{outliers_by_pool_fragments}
\alias{outliers_by_pool_fragments}
\title{Identify and flag outliers based on pool fragments.}
\usage{
outliers_by_pool_fragments(
  metadata,
  key = "BARCODE_MUX",
  outlier_p_value_threshold = 0.01,
  normality_test = FALSE,
  normality_p_value_threshold = 0.05,
  transform_log2 = TRUE,
  per_pool_test = TRUE,
  pool_col = "PoolID",
  min_samples_per_pool = 5,
  flag_logic = "AND",
  keep_calc_cols = TRUE,
  report_path = default_report_path()
)
}
\arguments{
\item{metadata}{The metadata data frame}

\item{key}{A character vector of numeric column names}

\item{outlier_p_value_threshold}{The p value threshold for a read to be
considered an outlier}

\item{normality_test}{Perform normality test? Normality is assessed for
each column in the key using Shapiro-Wilk test and if the values do not
follow a normal distribution, other calculations are skipped}

\item{normality_p_value_threshold}{Normality threshold}

\item{transform_log2}{Perform a log2 trasformation on values prior the
actual calculations?}

\item{per_pool_test}{Perform the test for each pool?}

\item{pool_col}{A character vector of the names of the columns that
uniquely identify a pool}

\item{min_samples_per_pool}{The minimum number of samples that a pool
needs to contain in order to be processed - relevant only if
\code{per_pool_test = TRUE}}

\item{flag_logic}{A character vector of logic operators to obtain a
global flag formula - only relevant if the key is longer than one.
All operators must be chosen between:
AND, OR, XOR, NAND, NOR, XNOR}

\item{keep_calc_cols}{Keep the calculation columns in the output data frame?}

\item{report_path}{The path where the report file should be saved.
Can be a folder, a file or NULL if no report should be produced.
Defaults to \code{{user_home}/ISAnalytics_reports}.}
}
\value{
A data frame of metadata with the column \code{to_remove}
}
\description{
\lifecycle{experimental}
Identify and flag outliers
}
\details{
This particular test calculates for each column in the key
\itemize{
\item The zscore of the values
\item The tstudent of the values
\item The the distribution of the tstudent values
}

Optionally the test can be performed for each pool and a normality test
can be run prior the actual calculations.
Samples are flagged if this condition is respected:
\itemize{
\item tdist < outlier_p_value_threshold & zscore < 0
}

If the key contains more than one column an additional flag logic can be
specified for combining the results.
Example:
let's suppose the key contains the names of two columns, X and Y
\code{key = c("X", "Y")}
if we specify the the argument \code{flag_logic = "AND"} then the reads will
be flagged based on this global condition:
(tdist_X < outlier_p_value_threshold & zscore_X < 0) AND
(tdist_Y < outlier_p_value_threshold & zscore_Y < 0)

The user can specify one or more logical operators that will be applied
in sequence.
}
\examples{
data("association_file", package = "ISAnalytics")
flagged <- outliers_by_pool_fragments(association_file,
    report_path = NULL
)
head(flagged)
}
\seealso{
Other Outlier tests: 
\code{\link{available_outlier_tests}()}
}
\concept{Outlier tests}