# Required packages
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(dplyr)

# Download the files if needed
download_if_missing <- function(url, dest) {
  if (!file.exists(dest)) download.file(url, dest)
}

base_url <- "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/"
download_if_missing(paste0(base_url, "knownGene.txt.gz"), "knownGene.txt.gz")
download_if_missing(paste0(base_url, "kgXref.txt.gz"), "kgXref.txt.gz")

# Read knownGene table
knownGene <- read.delim(gzfile("knownGene.txt.gz"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(knownGene) <- c(
  "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
  "exonCount", "exonStarts", "exonEnds", "proteinID", "alignID"
)

# Read kgXref table
kgXref <- read.delim(gzfile("kgXref.txt.gz"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(kgXref) <- c(
  "kgID", "mRNA", "spID", "spDisplayID", "geneSymbol",
  "refseq", "protAcc", "description", "alias1", "alias2"
)

# Join knownGene and kgXref
joined <- inner_join(knownGene, kgXref, by = c("name" = "kgID"))

# Compute summary per geneSymbol
refGenes_hg38 <- joined %>%
  filter(!is.na(geneSymbol) & geneSymbol != "") %>%
  group_by(geneSymbol, chrom, strand) %>%
  summarise(
    min_txStart = min(txStart, na.rm = TRUE),
    max_txEnd = max(txEnd, na.rm = TRUE),
    minmax_TxLen = max(txEnd, na.rm = TRUE) - min(txStart, na.rm = TRUE),
    average_TxLen = mean(txEnd - txStart, na.rm = TRUE),
    name = first(name),  # First transcript ID
    min_cdsStart = min(cdsStart, na.rm = TRUE),
    max_cdsEnd = max(cdsEnd, na.rm = TRUE),
    minmax_CdsLen = max(cdsEnd, na.rm = TRUE) - min(cdsStart, na.rm = TRUE),
    average_CdsLen = mean(cdsEnd - cdsStart, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(
    name2 = geneSymbol,
    strand = strand,
    chrom = chrom
  ) %>%
  select(
    name2, chrom, strand, min_txStart, max_txEnd, minmax_TxLen, average_TxLen,
    name, min_cdsStart, max_cdsEnd, minmax_CdsLen, average_CdsLen
  ) %>%
  arrange(chrom, min_txStart)

refGenes_hg38$chrom <- sub("^chr", "", refGenes_hg38$chrom)

refGenes_hg38_accepted <- refGenes_hg38 %>%
  filter(chrom %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
                      "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y"))

# Save as .RData
save(refGenes_hg38, file = "data/refGenes_mm10.RData")

message("✅ Saved refGenes_hg38.RData and .tsv with matched structure.")
