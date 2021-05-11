library(ISAnalytics)
func_name <- c(
  "import_Vispa2_stats"
)

#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
withr::local_options(
  list(ISAnalytics.widgets = FALSE,
       ISAnalytics.verbose = FALSE)
)

test_af_path <- system.file("extdata", "ex_association_file.tsv",
                            package = "ISAnalytics"
)

# Path to correct file system example
path_root_correct <- system.file("extdata", "fs.zip",
                                 package = "ISAnalytics"
)
root_correct <- unzip_file_system(path_root_correct, "fs")

# Path to incorrect file system example
path_root_err <- system.file("extdata", "fserr.zip",
                             package = "ISAnalytics"
)
root_err <- unzip_file_system(path_root_err, "fserr")

af_corr <- import_association_file(test_af_path, root = root_correct)
af_err <- import_association_file(test_af_path, root = root_err)

#------------------------------------------------------------------------------#
# Test import_Vispa2_stats
#------------------------------------------------------------------------------#
test_that(paste(func_name[1], "correct output no issues"), {
  imported <- import_Vispa2_stats(af_corr, join_with_af = FALSE)
  expected_pools_tags <- af_corr %>%
    dplyr::distinct(.data$concatenatePoolIDSeqRun, .data$TagSequence)
  expect_true(all(expected_pools_tags$concatenatePoolIDSeqRun %in%
                    unique(imported$POOL)))
  expect_true(all(expected_pools_tags$TagSequence %in%
                    unique(imported$TAG)))
})

test_that(paste(func_name[1], "correct output issues"), {
  imported <- import_Vispa2_stats(af_err, join_with_af = FALSE)
  expected_pools_tags <- af_err %>%
    dplyr::distinct(.data$concatenatePoolIDSeqRun, .data$TagSequence) %>%
    dplyr::filter(.data$concatenatePoolIDSeqRun == "POOL6-1")
  expect_true(all(expected_pools_tags$concatenatePoolIDSeqRun %in%
                    unique(imported$POOL)))
  expect_true(all(expected_pools_tags$TagSequence %in%
                    unique(imported$TAG)))
  expect_true(all(unique(imported$POOL) == "POOL6-1"))
})
