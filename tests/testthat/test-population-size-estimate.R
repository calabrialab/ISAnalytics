library(ISAnalytics)
func_name <- c("HSC_population_size_estimate")
#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
test_meta <- readRDS(system.file("testdata", "test_population_meta.Rds",
                                 package = "ISAnalytics"))

test_data <- readRDS(system.file("testdata", "test_population_iss.Rds",
                                 package = "ISAnalytics"))

test_expected <- readRDS(system.file("testdata", "test_population_expected.Rds",
                                 package = "ISAnalytics"))

#------------------------------------------------------------------------------#
# Test HSC_population_size_estimate
#------------------------------------------------------------------------------#
test_that(paste(func_name[1], "produces expected output"), {
  popul <- HSC_population_size_estimate(x = test_data,
                                        metadata = test_meta,
                                        stable_timepoints = c(9,12,13,18,30))
  expect_equal(popul, test_expected)
})


test_that(paste(func_name[1], "produces output missing NumIS"), {
  mod_meta <- test_meta %>% dplyr::select(-.data$NumIS)
  popul <- HSC_population_size_estimate(x = test_data,
                                        metadata = mod_meta,
                                        stable_timepoints = c(9,12,13,18,30))
  expect_equal(popul, test_expected)
})
