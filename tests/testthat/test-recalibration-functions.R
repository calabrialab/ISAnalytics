context("Recalibration functions")

library(ISAnalytics)

#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
# Annotated new matrix
example_path <- system.file("extdata", "ex_annotated_ISMatrix.tsv.xz",
                            package = "ISAnalytics"
)

isa_df <- import_single_Vispa2Matrix(example_path)

# Samples
sample_group1 <- tibble::tibble(chr = c(1,1),
                                integration_locus = c(14568, 14565),
                                strand = c("+", "+"),
                                Value = c(126, 20))

sample_group2 <- tibble::tibble(chr = c(rep_len(1, 6)),
                                integration_locus = c(14572, 15172,
                                                      14575, 14571,
                                                      14577, 14581),
                                strand = c(rep_len("+", 6)),
                                Value = c(70, 150, 120, 1400, 36, 15))

sample_group3 <- tibble::tibble(chr = rep_len(3, 7),
                                integration_locus = c(16380,
                                                      16396,
                                                      16402,
                                                      16395,
                                                      16378,
                                                      16399,
                                                      16387),
                                strand = rep_len("+", 7),
                                Value = c(1846, 64, 543, 89, 123, 886,
                                          48))

#------------------------------------------------------------------------------#
# Tests .locus_distance
#------------------------------------------------------------------------------#
test_that(".locus_distance calculates distance correctly", {
  locus1 <- 14598
  locus2 <- 14595
  dist <- .locus_distance(locus1, locus2)
  expect_equal(dist, 3)
  dist <- .locus_distance(locus2, locus1)
  expect_equal(dist, 3)
})

#------------------------------------------------------------------------------#
# Tests .keep_first
#------------------------------------------------------------------------------#
test_that(".keep_first returns tibble for tibble", {
  res <- .keep_first(sample_group1, inverted = FALSE)
  expect_true(tibble::is_tibble(res))
  expect_equal(res$integration_locus, c(14568))
  expect_equal(res$Value, 146)
})

test_that(".keep_first returns vector for vector", {
  res <- .keep_first(sample_group1$Value, inverted = FALSE)
  expect_true(is.numeric(res))
  expect_equal(res, 126)
  res <- .keep_first(sample_group1$Value, inverted = TRUE)
  expect_true(is.numeric(res))
  expect_equal(res, 20)
})

#------------------------------------------------------------------------------#
# Tests .keep_central
#------------------------------------------------------------------------------#
test_that(".keep_central returns tibble for tibble", {
  res <- .keep_central(sample_group1, inverted = FALSE)
  expect_true(tibble::is_tibble(res))
  expect_equal(res$integration_locus, c(14565))
  expect_equal(res$Value, 146)
})

test_that(".keep_central returns vector for vector", {
  res <- .keep_central(sample_group1$Value, inverted = FALSE)
  expect_true(is.numeric(res))
  expect_equal(res, 20)
  res <- .keep_central(sample_group1$Value, inverted = TRUE)
  expect_true(is.numeric(res))
  expect_equal(res, 126)
})

#------------------------------------------------------------------------------#
# Tests .keep_max_value
#------------------------------------------------------------------------------#
test_that(".keep_max_value returns tibble for tibble", {
  res <- .keep_max_value(sample_group1, second_choice = "keep_first",
                         inverted = FALSE)
  expect_true(tibble::is_tibble(res))
  expect_equal(res$integration_locus, c(14568))
  expect_equal(res$Value, c(146))
})

test_that(".keep_max_value returns vector for vector", {
  res <- .keep_max_value(sample_group1$Value, second_choice = "keep_first",
                         inverted = FALSE)
  expect_true(is.numeric(res))
  expect_equal(res, 1)
  res <- .keep_max_value(sample_group1$Value, second_choice = "keep_first",
                         inverted = TRUE)
  expect_true(is.numeric(res))
  expect_equal(res, 2)
})

test_that(".keep_max_value uses second criteria if value not comparable", {
  vec <- c(first = 126, second = 126)
  res <- .keep_max_value(vec, second_choice = "keep_first",
                         inverted = FALSE)
  expect_true(names(res) == c("first"))
  res <- .keep_max_value(vec, second_choice = "keep_first",
                         inverted = TRUE)
  expect_true(names(res) == c("second"))

  res <- .keep_max_value(vec, second_choice = "keep_central",
                         inverted = FALSE)
  expect_true(names(res) == c("second"))
  res <- .keep_max_value(vec, second_choice = "keep_central",
                         inverted = TRUE)
  expect_true(names(res) == c("first"))
})

#------------------------------------------------------------------------------#
# Tests .check_window_criteria
#------------------------------------------------------------------------------#
test_that(".check_window_criteria runs correctly", {
  mod <- sample_group2 %>% dplyr::arrange(.data$integration_locus) %>%
    dplyr::slice(1:3)
  values <- c(first = mod$Value[1], center = mod$Value[2],
              last =  mod$Value[3])
  check <- .check_window_criteria(c("max_value", "keep_first"), values,
                                  subset = NULL)
  expect_equal(names(check), c("center", "last"))
  check <- .check_window_criteria(c("keep_first", "keep_central"), values,
                                  subset = NULL)
  expect_equal(names(check), c("center", "last"))
  check <- .check_window_criteria(c("keep_central", "keep_first"), values,
                                  subset = NULL)
  expect_equal(names(check), c("first", "last"))
  check <- .check_window_criteria(c("max_value", "keep_first"), values,
                                  subset = c(1,2))
  expect_equal(names(check), c("center"))
  check <- .check_window_criteria(c("keep_first", "keep_central"), values,
                                  subset = c(1,2))
  expect_equal(names(check), c("center"))
  check <- .check_window_criteria(c("keep_central", "keep_first"), values,
                                  subset = c(1,2))
  expect_equal(names(check), c("first"))
})

#------------------------------------------------------------------------------#
# Tests .window
#------------------------------------------------------------------------------#
test_that(".window returns correctly NULL if no row to drop", {
  mod <- sample_group2 %>% dplyr::arrange(.data$integration_locus)
  wind_res <- .window(indexes = c(first = 1, center = 2, last = 3),
                      loci = c(first = mod$integration_locus[1],
                                  center = mod$integration_locus[2],
                                  last = mod$integration_locus[3]),
                      values = c(first = mod$Value[1],
                                    center = mod$Value[2],
                                    last = mod$Value[3]),
                      criterias = c("max_value", "keep_first"),
                      threshold = 0)
  expect_true(is.null(wind_res))
})

test_that(".window returns correctly 1/3 rows", {
  mod <- sample_group3 %>% dplyr::arrange(.data$integration_locus)
  wind_res <- .window(indexes = c(first = 3, center = 4, last = 5),
                      loci = c(first = mod$integration_locus[3],
                               center = mod$integration_locus[4],
                               last = mod$integration_locus[5]),
                      values = c(first = mod$Value[3],
                                 center = mod$Value[4],
                                 last = mod$Value[5]),
                      criterias = c("max_value", "keep_first"),
                      threshold = 4)
  expected <- list(drop = c(last = 5), collapse_on = c(center = 4), value = 153)
  expect_equal(wind_res, expected)
})


test_that(".window returns correctly indexes if 2/3 row to drop", {
  mod <- sample_group2 %>% dplyr::arrange(.data$integration_locus)
  wind_res <- .window(indexes = c(first = 1, center = 2, last = 3),
                      loci = c(first = mod$integration_locus[1],
                               center = mod$integration_locus[2],
                               last = mod$integration_locus[3]),
                      values = c(first = mod$Value[1],
                                 center = mod$Value[2],
                                 last = mod$Value[3]),
                      criterias = c("max_value", "keep_first"),
                      threshold = 4)
  expected <- list(drop = c(center = 2, last = 3),
                collapse_on = c(first = 1), value = 1590)
  expect_equal(wind_res, expected)
})

#------------------------------------------------------------------------------#
# Tests .window_slide
#------------------------------------------------------------------------------#
test_that(".window_slide produces expected output", {
  mod <- sample_group2 %>% dplyr::arrange(.data$integration_locus)
  res <- .window_slide(mod, start = 1, criterias = c("max_value", "keep_first"),
                       threshold = 4)
  expected <- tibble::tibble(chr = rep_len(1, 4),
                           integration_locus = c(14571,
                                                 14577,
                                                 14581,
                                                 15172),
                           strand = rep_len("+", 4),
                           Value = c(1590, 36, 15, 150))
  expect_identical(res, expected)
  mod <- sample_group3 %>% dplyr::arrange(.data$integration_locus)
  res <- .window_slide(mod, start = 1, criterias = c("max_value", "keep_first"),
                       threshold = 4)
  expected <- tibble::tibble(chr = rep_len(3, 4),
                             integration_locus = c(16380,
                                                   16387,
                                                   16395,
                                                   16399),
                             strand = rep_len("+", 4),
                             Value = c(1969, 48, 153, 1429))
  expect_identical(res, expected)
})

#------------------------------------------------------------------------------#
# Tests .sliding_window
#------------------------------------------------------------------------------#
test_that(".sliding_window returns correct output for group of 2", {
  res <- .sliding_window(x = sample_group1, threshold = 4,
                         keep_criteria = c("max_value", "keep_first"))
  expected <- tibble::tibble(chr = c(1),
                             integration_locus = c(14568),
                             strand = c("+"),
                             Value = c(146))
  expect_identical(res, expected)
})

test_that(".sliding_window returns correct output for group of 3+", {
  res <- .sliding_window(x = sample_group2, threshold = 4,
                         keep_criteria = c("max_value", "keep_first"))
  expected <- tibble::tibble(chr = rep_len(1, 4),
                             integration_locus = c(14571,
                                                   14577,
                                                   14581,
                                                   15172),
                             strand = rep_len("+", 4),
                             Value = c(1590, 36, 15, 150))
  expect_identical(res, expected)
  res <- .sliding_window(x = sample_group3, threshold = 4,
                         keep_criteria = c("max_value", "keep_first"))
  expected <- tibble::tibble(chr = rep_len(3, 4),
                             integration_locus = c(16380,
                                                   16387,
                                                   16395,
                                                   16399),
                             strand = rep_len("+", 4),
                             Value = c(1969, 48, 153, 1429))
  expect_identical(res, expected)
})

#------------------------------------------------------------------------------#
# Tests compute_near_integrations
#------------------------------------------------------------------------------#
## Test input
test_that("compute_near_integrations stops if x is not single tibble", {
  expect_error({
    res <- compute_near_integrations(x = 1)
  })
  expect_error({
    res <- compute_near_integrations(x = c(isa_df, isa_df))
  })
})
test_that("compute_near_integrations stops if x lacks columns", {
  isa_df_mod <- isa_df %>% dplyr::select(-c("chr"))
  expect_error({
    res <- compute_near_integrations(x = isa_df_mod)
  })
  isa_df_mod <- isa_df %>% dplyr::select(-c("strand"))
  expect_error({
    res <- compute_near_integrations(x = isa_df_mod)
  })
  isa_df_mod <- isa_df %>% dplyr::select(-c("integration_locus"))
  expect_error({
    res <- compute_near_integrations(x = isa_df_mod)
  })
  isa_df_mod <- isa_df %>% dplyr::select(-c("Value"))
  expect_error({
    res <- compute_near_integrations(x = isa_df_mod)
  })
})
test_that("compute_near_integrations stops if treshold is incorrect", {
  expect_error({
    res <- compute_near_integrations(x = isa_df, threshold = "5")
  })
  expect_error({
    res <- compute_near_integrations(x = isa_df, threshold = c(1,2))
  })
})
test_that("compute_near_integrations stops if keep_criteria is incorrect", {
  expect_error({
    res <- compute_near_integrations(x = isa_df, threshold = 4,
                                     keep_criteria = c(1,2))
  })
  expect_error({
    res <- compute_near_integrations(x = isa_df, threshold = 4,
                                     keep_criteria = c("a", "b"))
  })
  expect_error({
    res <- compute_near_integrations(x = isa_df, threshold = 4,
                                     keep_criteria = c("max_value", "first"))
  })
})
test_that("compute_near_integrations stops if strand_specific is incorrect", {
  expect_error({
    res <- compute_near_integrations(x = isa_df, strand_specific = "TRUE")
  })
  expect_error({
    res <- compute_near_integrations(x = isa_df,
                                     strand_specific = c(TRUE,FALSE))
  })
})
## TEST VALUES
test_that("compute_near_integrations produces correct output for sample", {
  matrix <- sample_group2 %>% dplyr::bind_rows(sample_group3)
  near <- compute_near_integrations(matrix)
  expected <- tibble::tibble(chr = c(rep_len(1, 4), rep_len(3, 4)),
                             integration_locus = c(14571,
                                                   14577,
                                                   14581,
                                                   15172,
                                                   16380,
                                                   16387,
                                                   16395,
                                                   16399),
                             strand = rep_len("+", 8),
                             Value = c(1590, 36, 15, 150,
                                       1969, 48, 153, 1429))
  expect_equal(near, expected)
})
