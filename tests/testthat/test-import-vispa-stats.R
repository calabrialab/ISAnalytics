library(ISAnalytics)
func_name <- c(
    "import_Vispa2_stats"
)

#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
withr::local_options(
    list(
        ISAnalytics.reports = FALSE,
        ISAnalytics.verbose = FALSE
    )
)

test_af_path <- system.file("extdata", "asso.file.tsv.gz",
    package = "ISAnalytics"
)
fs_path <- system.file("extdata", "fs.zip", package = "ISAnalytics")
root <- unzip_file_system(fs_path, "fs")

af <- import_association_file(test_af_path, root = root)

#------------------------------------------------------------------------------#
# Test import_Vispa2_stats
#------------------------------------------------------------------------------#
test_that(paste(func_name[1], "correct output no issues"), {
    imported <- import_Vispa2_stats(af, join_with_af = FALSE)
    expected_pools_tags <- af %>%
        dplyr::distinct(.data$concatenatePoolIDSeqRun, .data$TagSequence)
    expect_true(all(expected_pools_tags$concatenatePoolIDSeqRun %in%
        unique(imported$POOL)))
    expect_true(all(expected_pools_tags$TagSequence %in%
        unique(imported$TAG)))
})
