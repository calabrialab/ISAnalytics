# Reference gene annotation for hg38 or mm10.

A gene-level annotation dataset derived from the UCSC knownGene and
kgXref tables for the hg38 or mm10 genome assembly. This data aggregates
transcript-level information into gene-level summary statistics,
including transcript span, CDS length, and average values across
isoforms. It is the hg38 equivalent of `refGenes_hg19`, or mm10
equivalent of `refGenes_mm9`, updated using Ensembl-based transcript IDs
from GENCODE.

These objects are tibbles (`tbl_df`) and inherit from `data.frame`.

## Usage

``` r
data("refGenes_hg38")

data("refGenes_mm10")
```

## Format

A tibble with one row per gene and the following columns:

- name2:

  Gene symbol (e.g., A1CF)

- chrom:

  Chromosome (e.g., chr10)

- strand:

  Strand direction, "+" or "-"

- min_txStart:

  Minimum transcript start position across all isoforms

- max_txEnd:

  Maximum transcript end position across all isoforms

- minmax_TxLen:

  Gene length computed as max_txEnd - min_txStart

- average_TxLen:

  Average transcript length across isoforms

- name:

  Transcript ID (typically Ensembl ID in hg38, e.g., ENST00000...)

- min_cdsStart:

  Minimum CDS start position

- max_cdsEnd:

  Maximum CDS end position

- minmax_CdsLen:

  CDS length computed as max_cdsEnd - min_cdsStart

- average_CdsLen:

  Average CDS length across isoforms

An object of class `tbl_df` (inherits from `tbl`, `data.frame`) with
55316 rows and 12 columns.

## Source

UCSC Genome Browser:
<https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/>

UCSC Genome Browser:
<https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/>

## Functions

- `refGenes_mm10`: Data frame for murine mm10 genome
