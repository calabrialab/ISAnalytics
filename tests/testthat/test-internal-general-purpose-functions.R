
#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
# sample_zip_file <- system.file("testdata", "fs.zip",
#     package = "ISAnalytics"
# )
# sample_tsv_file <- system.file("testdata", "ex_association_file.tsv.gz",
#     package = "ISAnalytics"
# )

#------------------------------------------------------------------------------#
# Test .check_file_extension
#------------------------------------------------------------------------------#
# test_that(".check_file_extension works with compressed folder", {
#     expect_equal(.check_file_extension(sample_zip_file), "")
# })
#
# test_that(".check_file_extension works with non comp file", {
#     expect_equal(.check_file_extension(sample_tsv_file), "tsv")
# })
#
# test_that(".check_file_extension works with multiple input", {
#     checks <- .check_file_extension(
#         c(
#             sample_zip_file,
#             sample_tsv_file
#         )
#     )
#     expected <- c("", "tsv")
#     expect_equal(checks, expected)
# })
