#' Example of imported multi-quantification integration matrices.
#'
#' The data was obtained manually by simulating real research
#' data.
#'
#' @usage data("integration_matrices")
#'
#' @format Data frame with 1689 rows and 8 columns
#' \describe{
#'  \item{chr}{The chromosome number (as character)}
#'  \item{integration_locus}{Number of the base at
#'  which the viral insertion occurred}
#'  \item{strand}{Strand of the integration}
#'  \item{GeneName}{Symbol of the closest gene}
#'  \item{GeneStrand}{Strand of the closest gene}
#'  \item{CompleteAmplificationID}{Unique sample identifier}
#'  \item{seqCount}{Value of the sequence count quantification}
#'  \item{fragmentEstimate}{Value of the fragment estimate quantification}
#' }
"integration_matrices"

#' Example of association file.
#'
#' The data was obtained manually by simulating real research
#' data.
#'
#' @usage data("association_file")
#'
#' @description This file is a simple example of association file. Use it as
#' reference to properly fill out yours.
#' To generate an empty association file to fill see the
#' `generate_blank_association_file()` function.
#' @seealso \code{\link{generate_blank_association_file}}
"association_file"

#' Gene annotation files for hg19, mm9 and mm10.
#'
#' @description
#' This file was obtained following this steps:
#'
#' 1. Download from {http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/}
#' the refGene.sql, knownGene.sql, knownToRefSeq.sql, kgXref.sql tables
#' 2. Import everything it in mysql
#' 3. Generate views for annotation:
#'
#' ```
#' SELECT kg.`chrom`, min(kg.cdsStart) as CDS_minStart,
#' max(kg.`cdsEnd`) as CDS_maxEnd, k2a.geneSymbol,
#' kg.`strand` as GeneStrand, min(kg.txStart) as TSS_minStart,
#' max(kg.txEnd) as TSS_maxStart,
#' kg.proteinID as ProteinID, k2a.protAcc as ProteinAcc, k2a.spDisplayID
#' FROM `knownGene` AS kg JOIN kgXref AS k2a
#' ON BINARY kg.name = k2a.kgID COLLATE latin1_bin
#' -- latin1_swedish_ci
#' -- WHERE k2a.spDisplayID IS NOT NULL and (k2a.`geneSymbol` LIKE 'Tcra%' or
#' k2a.`geneSymbol` LIKE 'TCRA%')
#' WHERE (k2a.spDisplayID IS NOT NULL or k2a.spDisplayID NOT LIKE '')
#' and k2a.`geneSymbol` LIKE 'Tcra%'
#' group by kg.`chrom`, k2a.geneSymbol
#' ORDER BY kg.chrom ASC , kg.txStart ASC
#' ```
#' @usage data("refGenes_hg19")
"refGenes_hg19"
#' @describeIn refGenes_hg19 Data frame for murine mm9 genome
#' @usage data("refGenes_mm9")
"refGenes_mm9"

#' Data frames for proto-oncogenes (human and mouse)
#' amd tumor-suppressor genes from UniProt.
#'
#' @description
#' The file is simply a result of a research with the keywords
#' "proto-oncogenes" and "tumor suppressor" for the target genomes
#' on UniProt database.
#' @usage data("proto_oncogenes")
"proto_oncogenes"
#' @describeIn proto_oncogenes Data frame for tumor suppressor genes
#' @usage data("tumor_suppressors")
"tumor_suppressors"
