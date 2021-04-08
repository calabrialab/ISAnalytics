library(ISAnalytics)

#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
sample_xz_file <- system.file("extdata", "ex_annotated_ISMatrix.tsv.xz",
    package = "ISAnalytics"
)
sample_zip_file <- system.file("extdata", "fs.zip",
    package = "ISAnalytics"
)
sample_tsv_file <- system.file("extdata", "ex_association_file.tsv",
    package = "ISAnalytics"
)

#------------------------------------------------------------------------------#
# Test .check_file_extension
#------------------------------------------------------------------------------#
test_that(".check_file_extension works with compressed file", {
    expect_equal(.check_file_extension(sample_xz_file), "tsv")
})

test_that(".check_file_extension works with compressed folder", {
    expect_equal(.check_file_extension(sample_zip_file), "")
})

test_that(".check_file_extension works with non comp file", {
    expect_equal(.check_file_extension(sample_tsv_file), "tsv")
})

test_that(".check_file_extension works with multiple input", {
    checks <- .check_file_extension(
        c(
            sample_xz_file,
            sample_zip_file,
            sample_tsv_file
        )
    )
    expected <- c("tsv", "", "tsv")
    expect_equal(checks, expected)
})
