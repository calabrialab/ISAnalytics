#' Example file of annotated integration matrix.
#'
#' @description This file is a very simple and compressed example to showcase
#' some of the functionalities of ISAnalytics package. The general structure of
#' the matrix is obtained as a product of the Vispa2 pipeline + create_matrix +
#' annotate_matrix functions.
#' For more information regarding Vispa2 please read this article:
#'
#' \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5702242/}{VISPA2:
#' A Scalable Pipeline for High-Throughput Identification
#' and Annotation of Vector Integration Sites}.\cr
#'
#' The headers of the dataset are standard except for experimental data:
#' * chr : the chromosome number
#' * integration_locus: the site of integration
#' * strand: the DNA strand on which the integration took place
#' * GeneName : name of the gene
#' * GeneStrand : strand of the gene
#' * exp_... : names of the experimental data (not standard, can be anything)
#'
"ex_annotated_ISMatrix"

#' Example file of old style, not annotated integration matrix.
#'
#' @description This file is a very simple and compressed example to showcase
#' some of the functionalities of ISAnalytics package. The general structure of
#' the matrix is obtained as a product of the Vispa2 pipeline + create_matrix +
#' annotate_matrix functions.
#' For more information regarding Vispa2 please read this article:
#'
#' \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5702242/}{VISPA2:
#' A Scalable Pipeline for High-Throughput Identification
#' and Annotation of Vector Integration Sites}.\cr
#'
#' The difference with standard Vispa2 annotated matrices lies in the fact that
#' the genomic coordinates are stored in a single column called IS_genomicID.
#'
#' Please note that this example and associated functionality is for
#' compatibility reasons only: if you're using the current version of Vispa2
#' you should always obtain separated columns for genomic coordinates and
#' the corresponding annotation columns.
#'
"ex_old_style_ISMatrix"

#' Example file of not annotate integration matrix.
#'
#' @description This file is a very simple and compressed example to showcase
#' some of the functionalities of ISAnalytics package. The general structure of
#' the matrix is obtained as a product of the Vispa2 pipeline + create_matrix +
#' annotate_matrix functions.
#' For more information regarding Vispa2 please read this article:
#'
#' \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5702242/}{VISPA2:
#' A Scalable Pipeline for High-Throughput Identification
#' and Annotation of Vector Integration Sites}.\cr
#'
#' The difference with standard Vispa2 annotated matrices lies in the fact that
#' there are no annotation columns, namely "GeneName" and "GeneStrand".
#'
#'Please note that this example and associated functionality is for
#' compatibility reasons only: if you're using the current version of Vispa2
#' you should always obtain separated columns for genomic coordinates and
#' the corresponding annotation columns.
#'
"ex_notann_ISMatrix"

#' Example file of malformed integration matrix.
#'
#' @description This file is a very simple and compressed example to showcase
#' some of the functionalities of ISAnalytics package. The general structure of
#' the matrix is obtained as a product of the Vispa2 pipeline + create_matrix +
#' annotate_matrix functions.
#' For more information regarding Vispa2 please read this article:
#'
#' \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5702242/}{VISPA2:
#' A Scalable Pipeline for High-Throughput Identification
#' and Annotation of Vector Integration Sites}.\cr
#'
#' This matrix is malformed on purpose, missing one of the fundamental columns.
#' Trying to import this matrix or convert it to an ISADataFrame object should
#' always result in some kind of error.
#'
"ex_malformed_ISMatrix"
