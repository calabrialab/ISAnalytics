library(ISAnalytics)

#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
op <- withr::local_options(
    ISAnalytics.widgets = FALSE,
    ISAnalytics.verbose = FALSE
)
# Path to example association file
path_af <- system.file("extdata", "ex_association_file.tsv",
    package = "ISAnalytics"
)
# Path to correct file system example
path_root_correct <- system.file("extdata", "fs.zip",
    package = "ISAnalytics"
)
root_correct <- unzip_file_system(path_root_correct, "fs")
as_file <- import_association_file(path_af, root_correct, dates_format = "dmy")

#------------------------------------------------------------------------------#
# Tests generate_blank_association_file
#------------------------------------------------------------------------------#
test_that("generate_blank_association_file stops if path is not char", {
    expect_error({
        generate_blank_association_file(1)
    })
})

test_that("generate_blank_association_file works correctly", {
    temp <- tempfile()
    generate_blank_association_file(temp)
    af <- read.csv(temp, sep = "\t", check.names = FALSE, header = TRUE)
    expect_true(all(colnames(af) %in% association_file_columns()))
})

#------------------------------------------------------------------------------#
# Tests generate_Vispa2_launch_AF
#------------------------------------------------------------------------------#
## Testing input
test_that("generate_Vispa2_launch_AF stops if association_file is not df", {
    expect_error({
        generate_Vispa2_launch_AF(path_af, "x", "y", "z")
    })
})

test_that("generate_Vispa2_launch_AF stops if project is not char", {
    expect_error({
        generate_Vispa2_launch_AF(as_file, 1, "y", "z")
    })
})

test_that("generate_Vispa2_launch_AF stops if pool is not char", {
    expect_error({
        generate_Vispa2_launch_AF(as_file, "CLOEXP", 1, "z")
    })
})

test_that("generate_Vispa2_launch_AF stops if lengths of projects and pool is
          not the same", {
    expect_error({
        generate_Vispa2_launch_AF(
            as_file, c("CLOEXP", "PROJECT1100"),
            c("POOL6"), "z"
        )
    })
})

test_that("generate_Vispa2_launch_AF stops if path is incorrect", {
    expect_error({
        generate_Vispa2_launch_AF(as_file, c("CLOEXP"), c("POOL6"), 1)
    })
    expect_error({
        generate_Vispa2_launch_AF(as_file, c("CLOEXP"), c("POOL6"), c("x", "y"))
    })
})

test_that("generate_Vispa2_launch_AF stops if af is malformed", {
    af <- as_file %>% dplyr::select(-c(.data$ProjectID))
    expect_error({
        generate_Vispa2_launch_AF(af, c("CLOEXP"), c("POOL6"), 1)
    })
})

## Testing output
test_that("generate_Vispa2_launch_AF works for single pair", {
    temp <- withr::local_tempdir()
    project <- c("CLOEXP")
    pool <- c("POOL6")
    name <- paste0(project, "-", pool, "_AF.tsv")
    complete_path <- file.path(temp, name)
    complete_path <- gsub('"', "", gsub("\\\\", "/", complete_path))
    generate_Vispa2_launch_AF(as_file, project, pool, temp)
    df <- read.csv(complete_path, sep = "\t", header = FALSE)
    expect_true(ncol(df) == 11)
    expect_equal(df[, 1], df[, 2])
    expect_message({
        generate_Vispa2_launch_AF(as_file, c("CLOEXP"), c("x"), temp)
    })
})

test_that("generate_Vispa2_launch_AF works for multiple pair", {
    temp <- withr::local_tempdir()
    project <- c("CLOEXP", "PROJECT1100")
    pool <- c("POOL6", "ABX-LR-PL5-POOL14")
    name <- paste0(project, "-", pool, "_AF.tsv")
    complete_path <- file.path(temp, name)
    complete_path <- gsub('"', "", gsub("\\\\", "/", complete_path))
    generate_Vispa2_launch_AF(as_file, project, pool, temp)
    df <- read.csv(complete_path[1], sep = "\t", header = FALSE)
    expect_true(ncol(df) == 11)
    expect_equal(df[, 1], df[, 2])
    df <- read.csv(complete_path[2], sep = "\t", header = FALSE)
    expect_true(ncol(df) == 11)
    expect_equal(df[, 1], df[, 2])
})

#------------------------------------------------------------------------------#
# Tests as_sparse_matrix
#------------------------------------------------------------------------------#
smpl <- tibble::tibble(
    chr = c(1, 2, 3), integration_locus = c(1354, 5634, 4765),
    strand = c("+", "+", "+"),
    GeneName = c("GENE1", "GENE2", "GENE3"),
    GeneStrand = c("+", "+", "+"),
    CompleteAmplificationID = c("ID1", "ID2", "ID3"),
    Value = c(46, 546, 587)
)

smpl_multi <- tibble::tibble(
    chr = c(1, 2, 3),
    integration_locus = c(1354, 5634, 4765),
    strand = c("+", "+", "+"),
    GeneName = c("GENE1", "GENE2", "GENE3"),
    GeneStrand = c("+", "+", "+"),
    CompleteAmplificationID = c("ID1", "ID2", "ID3"),
    seqCount = c(46, 546, 587),
    fragmentEstimate = c(4234.5, 533.45, 5431.43),
    barcodeCount = c(46, 6456, 456)
)
test_that("as_sparse_matrix works with single matrix", {
    sparse <- as_sparse_matrix(smpl)
    expect_true(tibble::is_tibble(sparse))
    expect_true(all(c(
        "chr", "integration_locus", "strand",
        "GeneName", "GeneStrand", "ID1", "ID2", "ID3"
    ) %in%
        colnames(sparse)))
    expect_equal(nrow(sparse), 3)
})

test_that("as_sparse_matrix works with multi-quant matrix", {
    sparse <- as_sparse_matrix(smpl_multi)
    expect_true(!tibble::is_tibble(sparse) & is.list(sparse))
    expect_equal(names(sparse), c(
        "fragmentEstimate",
        "seqCount",
        "barcodeCount"
    ))
    expect_true(all(c(
        "chr", "integration_locus", "strand",
        "GeneName", "GeneStrand", "ID1", "ID2", "ID3"
    ) %in%
        colnames(sparse$seqCount)))
    expect_equal(nrow(sparse$seqCount), 3)
    expect_true(all(c(
        "chr", "integration_locus", "strand",
        "GeneName", "GeneStrand", "ID1", "ID2", "ID3"
    ) %in%
        colnames(sparse$fragmentEstimate)))
    expect_equal(nrow(sparse$fragmentEstimate), 3)
})

test_that("as_sparse_matrix works with list of matrices", {
    sparse <- as_sparse_matrix(list(smpl, smpl))
    expect_true(!tibble::is_tibble(sparse) & is.list(sparse))
    expect_true(all(c(
        "chr", "integration_locus", "strand",
        "GeneName", "GeneStrand", "ID1", "ID2", "ID3"
    ) %in%
        colnames(sparse[[1]])))
    expect_equal(nrow(sparse[[1]]), 3)
    expect_true(all(c(
        "chr", "integration_locus", "strand",
        "GeneName", "GeneStrand", "ID1", "ID2", "ID3"
    ) %in%
        colnames(sparse[[2]])))
    expect_equal(nrow(sparse[[2]]), 3)
})
