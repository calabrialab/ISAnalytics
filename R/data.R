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

#' Gene annotation files for hg19, mm9.
#'
#' @name refGenes_hg19
#' @description
#' This file was obtained following this steps:
#'
#' 1. Download from \url{http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/}
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

#' Reference gene annotation for hg38 or mm10.
#'
#' @name refGenes_hg38
#' @description
#' A gene-level annotation dataset derived from the UCSC knownGene and kgXref tables
#' for the hg38 or mm10 genome assembly. This data aggregates transcript-level information into
#' gene-level summary statistics, including transcript span, CDS length, and average values
#' across isoforms. It is the hg38 equivalent of `refGenes_hg19`,
#' or mm10 equivalent of `refGenes_mm9`, updated using Ensembl-based
#' transcript IDs from GENCODE.
#'
#' These objects are tibbles (`tbl_df`) and inherit from `data.frame`.
#'
#' @format A tibble with one row per gene and the following columns:
#' \describe{
#'   \item{name2}{Gene symbol (e.g., A1CF)}
#'   \item{chrom}{Chromosome (e.g., chr10)}
#'   \item{strand}{Strand direction, "+" or "-"}
#'   \item{min_txStart}{Minimum transcript start position across all isoforms}
#'   \item{max_txEnd}{Maximum transcript end position across all isoforms}
#'   \item{minmax_TxLen}{Gene length computed as max_txEnd - min_txStart}
#'   \item{average_TxLen}{Average transcript length across isoforms}
#'   \item{name}{Transcript ID (typically Ensembl ID in hg38, e.g., ENST00000...)}
#'   \item{min_cdsStart}{Minimum CDS start position}
#'   \item{max_cdsEnd}{Maximum CDS end position}
#'   \item{minmax_CdsLen}{CDS length computed as max_cdsEnd - min_cdsStart}
#'   \item{average_CdsLen}{Average CDS length across isoforms}
#' }
#' @source UCSC Genome Browser: \url{https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/}
#' @source UCSC Genome Browser: \url{https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/}
#' 
#' @usage data("refGenes_hg38")
"refGenes_hg38"
#' @describeIn refGenes_hg38 Data frame for murine mm10 genome
#' @usage data("refGenes_mm10")
"refGenes_mm10"

#' Data frames for proto-oncogenes (human and mouse)
#' and tumor-suppressor genes from UniProt.
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
