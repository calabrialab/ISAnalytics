% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis-functions.R
\name{comparison_matrix}
\alias{comparison_matrix}
\title{obtain a single integration matrix from individual quantification
matrices.}
\usage{
comparison_matrix(
  x,
  fragmentEstimate = "fragmentEstimate",
  seqCount = "seqCount",
  barcodeCount = "barcodeCount",
  cellCount = "cellCount",
  ShsCount = "ShsCount"
)
}
\arguments{
\item{x}{A named list of integration matrices, ideally obtained via
\link{import_parallel_Vispa2Matrices_interactive} or
\link{import_parallel_Vispa2Matrices_auto}. Names must be
quantification types.}

\item{fragmentEstimate}{The name of the output column for fragment
estimate values}

\item{seqCount}{The name of the output column for sequence
count values}

\item{barcodeCount}{The name of the output column for barcode count
values}

\item{cellCount}{The name of the output column for cell count values}

\item{ShsCount}{The name of the output column for Shs count values}
}
\value{
A tibble
}
\description{
\lifecycle{maturing}
Takes a list of integration matrices referring to different qunatification
types and merges them in a single data frame that has multiple
value columns, each renamed according to their quantification type
of reference.
}
\examples{
op <- options("ISAnalytics.widgets" = FALSE)
path <- system.file("extdata", "ex_association_file.tsv",
    package = "ISAnalytics"
)
root_pth <- system.file("extdata", "fs.zip", package = "ISAnalytics")
root <- unzip_file_system(root_pth, "fs")
matrices <- import_parallel_Vispa2Matrices_auto(
    association_file = path, root = root,
    quantification_type = c("fragmentEstimate", "seqCount"),
    matrix_type = "annotated", workers = 2, patterns = NULL,
    matching_opt = "ANY",
    dates_format = "dmy", multi_quant_matrix = FALSE
)
total_matrix <- comparison_matrix(matrices)
options(op)
}
\seealso{
\link{quantification_types}

Other Analysis functions: 
\code{\link{CIS_grubbs}()},
\code{\link{compute_abundance}()},
\code{\link{cumulative_count_union}()},
\code{\link{is_sharing}()},
\code{\link{sample_statistics}()},
\code{\link{separate_quant_matrices}()},
\code{\link{threshold_filter}()},
\code{\link{top_integrations}()}
}
\concept{Analysis functions}
