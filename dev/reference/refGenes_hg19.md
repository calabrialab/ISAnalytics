# Gene annotation files for hg19, mm9.

This file was obtained following this steps:

1.  Download from
    http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ the
    refGene.sql, knownGene.sql, knownToRefSeq.sql, kgXref.sql tables

2.  Import everything it in mysql

3.  Generate views for annotation:

    SELECT kg.`chrom`, min(kg.cdsStart) as CDS_minStart,
    max(kg.`cdsEnd`) as CDS_maxEnd, k2a.geneSymbol,
    kg.`strand` as GeneStrand, min(kg.txStart) as TSS_minStart,
    max(kg.txEnd) as TSS_maxStart,
    kg.proteinID as ProteinID, k2a.protAcc as ProteinAcc, k2a.spDisplayID
    FROM `knownGene` AS kg JOIN kgXref AS k2a
    ON BINARY kg.name = k2a.kgID COLLATE latin1_bin
    -- latin1_swedish_ci
    -- WHERE k2a.spDisplayID IS NOT NULL and (k2a.`geneSymbol` LIKE 'Tcra%' or
    k2a.`geneSymbol` LIKE 'TCRA%')
    WHERE (k2a.spDisplayID IS NOT NULL or k2a.spDisplayID NOT LIKE '')
    and k2a.`geneSymbol` LIKE 'Tcra%'
    group by kg.`chrom`, k2a.geneSymbol
    ORDER BY kg.chrom ASC , kg.txStart ASC

## Usage

``` r
data("refGenes_hg19")

data("refGenes_mm9")
```

## Format

An object of class `tbl_df` (inherits from `tbl`, `data.frame`) with
27275 rows and 12 columns.

An object of class `tbl_df` (inherits from `tbl`, `data.frame`) with
24487 rows and 12 columns.

## Functions

- `refGenes_mm9`: Data frame for murine mm9 genome
