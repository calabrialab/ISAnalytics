library(ISAnalytics)

#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
# op <- withr::local_options(
#     ISAnalytics.reports = FALSE,
#     ISAnalytics.verbose = FALSE
# )
# # Path to example association file
# path_af <- system.file("testdata", "ex_association_file.tsv.gz",
#     package = "ISAnalytics"
# )
# # Path to correct file system example
# path_root_correct <- system.file("testdata", "fs.zip",
#     package = "ISAnalytics"
# )
# root_correct <- unzip_file_system(path_root_correct, "fs")
# as_file <- import_association_file(path_af, root_correct, dates_format = "dmy")

#------------------------------------------------------------------------------#
# Tests .process_af_for_gen
#------------------------------------------------------------------------------#
test_that(".process_af_for_gen works for default af", {
  res <- .process_af_for_gen(af = "default")
  af_out <- res$af
  tag_out <- res$tag_list
  expected_tags <- list("pcr_repl_id" = "CompleteAmplificationID",
                        "project_id" = "ProjectID",
                        "tag_seq" = "TagSequence",
                        "vispa_concatenate" = "concatenatePoolIDSeqRun")
  af_sym <- "association_file"
  utils::data(list = af_sym, envir = rlang::current_env())
  association_file <- rlang::eval_tidy(rlang::sym(af_sym))
  expect_equal(tag_out, expected_tags)
  expect_equal(af_out,  association_file)
})

test_that(".process_af_for_gen works for custom af", {
  af_custom <- tibble::tribble(
    ~ ProjectID, ~ TagSequence, ~ CompleteAmplificationID,
    ~ concatenatePoolIDSeqRun, ~ X,
    "PJ01", "LTR75LC38", "ID1", "POOL1", "a",
    "PJ01", "LTR53LC32", "ID2", "POOL1", "b",
    "PJ01", "LTR27LC94", "ID3", "POOL2", "c",
    "PJ01", "LTR69LC52", "ID4", "POOL3", "d"
  )
  res <- .process_af_for_gen(af = af_custom)
  af_out <- res$af
  tag_out <- res$tag_list
  expected_tags <- list("pcr_repl_id" = "CompleteAmplificationID",
                        "project_id" = "ProjectID",
                        "tag_seq" = "TagSequence",
                        "vispa_concatenate" = "concatenatePoolIDSeqRun")
  expect_equal(tag_out, expected_tags)
  expect_equal(af_out,  af_custom)

  af_custom_err <- af_custom %>%
    dplyr::select(-.data$ProjectID)
  expect_error({
    res <- .process_af_for_gen(af = af_custom_err)
  }, class = "missing_req_col_err")
})

#------------------------------------------------------------------------------#
# Tests .process_m_for_gen
#------------------------------------------------------------------------------#
test_that(".process_m_for_gen works with default matrices", {
  af <- .process_af_for_gen(af = "default")
  tags <- af$tag_list
  af <- af$af
  res <- .process_m_for_gen(matrices = "default", af, tags)
  expect_equal(names(res), "PJ01")
  expect_equal(names(res$PJ01), c("POOL01-1", "POOL02-1",
                                  "POOL03-1", "POOL04-1"))
})

#------------------------------------------------------------------------------#
# Tests .process_iss_for_gen
#------------------------------------------------------------------------------#
test_that(".process_iss_for_gen works with default stats", {
  af <- .process_af_for_gen(af = "default")
  tags <- af$tag_list
  af <- af$af
  res <- .process_iss_for_gen(af, tags)
  expect_equal(names(res), "PJ01")
  expect_equal(names(res$PJ01), c("POOL01-1", "POOL02-1",
                                  "POOL03-1", "POOL04-1"))
})

#------------------------------------------------------------------------------#
# Tests .generate_correct
#------------------------------------------------------------------------------#
test_that(".generate_correct works correctly with defaults", {
  af <- .process_af_for_gen(af = "default")
  tags <- af$tag_list
  af <- af$af
  mat <- .process_m_for_gen(matrices = "default", af, tags)
  iss <- .process_iss_for_gen(af, tags)
  temp_dir <- tempdir()
  expect_error({
    path_out <- .generate_correct(temp_dir, sep_matrices = mat,
                                  sep_stats = iss)
  }, regexp = NA)
  proj_fold <- fs::path(path_out, "PJ01")
  quant_fold <- fs::path(proj_fold, "quantification")
  iss_fold <- fs::path(proj_fold, "iss")
  pools <- paste0(paste0("POOL0", 1:4), "-1")
  pools_quant <- fs::path(quant_fold, pools)
  pools_iss <- fs::path(iss_fold, pools)
  expect_true(all(pools_quant %in% fs::dir_ls(path_out, recurse = TRUE)))
  expect_true(all(pools_iss %in% fs::dir_ls(path_out, recurse = TRUE)))
  expect_true(all(purrr::map_lgl(pools_quant,
                                    ~ length(fs::dir_ls(.x)) == 2)))
  expect_true(all(purrr::map_lgl(pools_iss,
                                    ~ length(fs::dir_ls(.x)) == 1)))
})

#------------------------------------------------------------------------------#
# Tests .generate_incorrect
#------------------------------------------------------------------------------#
test_that(".generate_incorrect works correctly with defaults", {
  af <- .process_af_for_gen(af = "default")
  tags <- af$tag_list
  af <- af$af
  mat <- .process_m_for_gen(matrices = "default", af, tags)
  iss <- .process_iss_for_gen(af, tags)
  temp_dir <- tempdir()
  expect_error({
    path_out <- .generate_incorrect(temp_dir, sep_matrices = mat,
                                    sep_stats = iss)
  }, regexp = NA)
  proj_fold <- fs::path(path_out, "PJ01")
  quant_fold <- fs::path(proj_fold, "quantification")
  iss_fold <- fs::path(proj_fold, "iss")
  pools <- c("POOL01-1-err", paste0(paste0("POOL0", 2:4), "-1"))
  pools_quant <- fs::path(quant_fold, pools)
  pools_iss <- fs::path(iss_fold, pools)
  expect_true(all(pools_quant %in% fs::dir_ls(path_out, recurse = TRUE)))
  expect_true(all(pools_iss %in% fs::dir_ls(path_out, recurse = TRUE)))
  expect_true(all(purrr::map_lgl(pools_quant, ~ {
    split <- fs::path_split(.x)
    pool_fold <- split[[1]][length(split[[1]])]
    if (pool_fold == "POOL02-1" & length(fs::dir_ls(.x)) == 1) {
      return(TRUE)
    }
    return(length(fs::dir_ls(.x)) == 2)
  })))
  expect_true(all(purrr::map_lgl(pools_iss, ~ {
    split <- fs::path_split(.x)
    pool_fold <- split[[1]][length(split[[1]])]
    if (pool_fold == "POOL02-1" & length(fs::dir_ls(.x)) == 0) {
      return(TRUE)
    }
    return(length(fs::dir_ls(.x)) == 1)
  })))
})

#------------------------------------------------------------------------------#
# Tests transform_columns
#------------------------------------------------------------------------------#
test_that("transform_columns applies transform as expected", {
  test_df <- tibble::tribble(
    ~ A, ~B, ~C,
    1, 2, 3,
    4, 5, 6,
    7, 8, 9
  )
  to_apply <- list(A = ~stringr::str_pad(
    as.character(.x), pad = "0", side = "left", width = 2), B = ~ .x + 2)
  result <- transform_columns(test_df, to_apply)
  expect_equal(result$A, c("01", "04", "07"))
  expect_equal(result$B, c(4, 7, 10))
})

#------------------------------------------------------------------------------#
# Tests pcr_id_column
#------------------------------------------------------------------------------#
test_that("pcr_id_column errors if no pcr_repl_id in vars", {
  withr::local_options(list(ISAnalytics.af_specs = "default"))
  af_specs <- association_file_columns(TRUE)
  af_specs_mod <- af_specs %>%
    dplyr::filter(!.data$tag %in% c("pcr_repl_id"))
  withr::with_options(list(ISAnalytics.af_specs = af_specs_mod), {
    expect_error({
      col <- pcr_id_column()
    })
  })
})

#------------------------------------------------------------------------------#
# Tests generate_blank_association_file
#------------------------------------------------------------------------------#
# test_that("generate_blank_association_file stops if path is not char", {
#     expect_error({
#         generate_blank_association_file(1)
#     })
# })
#
# test_that("generate_blank_association_file works correctly", {
#     temp <- tempfile()
#     generate_blank_association_file(temp)
#     af <- read.csv(temp, sep = "\t", check.names = FALSE, header = TRUE)
#     expect_true(all(colnames(af) %in% association_file_columns()))
# })

#------------------------------------------------------------------------------#
# Tests generate_Vispa2_launch_AF
#------------------------------------------------------------------------------#
# ## Testing input
# test_that("generate_Vispa2_launch_AF stops if association_file is not df", {
#     expect_error({
#         generate_Vispa2_launch_AF(path_af, "x", "y", "z")
#     })
# })
#
# test_that("generate_Vispa2_launch_AF stops if project is not char", {
#     expect_error({
#         generate_Vispa2_launch_AF(as_file, 1, "y", "z")
#     })
# })
#
# test_that("generate_Vispa2_launch_AF stops if pool is not char", {
#     expect_error({
#         generate_Vispa2_launch_AF(as_file, "CLOEXP", 1, "z")
#     })
# })
#
# test_that("generate_Vispa2_launch_AF stops if lengths of projects and pool is
#           not the same", {
#     expect_error({
#         generate_Vispa2_launch_AF(
#             as_file, c("CLOEXP", "PROJECT1100"),
#             c("POOL6"), "z"
#         )
#     })
# })
#
# test_that("generate_Vispa2_launch_AF stops if path is incorrect", {
#     expect_error({
#         generate_Vispa2_launch_AF(as_file, c("CLOEXP"), c("POOL6"), 1)
#     })
#     expect_error({
#         generate_Vispa2_launch_AF(as_file, c("CLOEXP"), c("POOL6"), c("x", "y"))
#     })
# })
#
# test_that("generate_Vispa2_launch_AF stops if af is malformed", {
#     af <- as_file %>% dplyr::select(-c(.data$ProjectID))
#     expect_error({
#         generate_Vispa2_launch_AF(af, c("CLOEXP"), c("POOL6"), 1)
#     })
# })
#
# ## Testing output
# test_that("generate_Vispa2_launch_AF works for single pair", {
#     temp <- withr::local_tempdir()
#     project <- c("CLOEXP")
#     pool <- c("POOL6")
#     name <- paste0(project, "-", pool, "_AF.tsv")
#     complete_path <- file.path(temp, name)
#     complete_path <- gsub('"', "", gsub("\\\\", "/", complete_path))
#     generate_Vispa2_launch_AF(as_file, project, pool, temp)
#     df <- read.csv(complete_path, sep = "\t", header = FALSE)
#     expect_true(ncol(df) == 11)
#     expect_equal(df[, 1], df[, 2])
#     expect_message({
#         generate_Vispa2_launch_AF(as_file, c("CLOEXP"), c("x"), temp)
#     })
# })
#
# test_that("generate_Vispa2_launch_AF works for multiple pair", {
#     temp <- withr::local_tempdir()
#     project <- c("CLOEXP", "PROJECT1100")
#     pool <- c("POOL6", "ABX-LR-PL5-POOL14")
#     name <- paste0(project, "-", pool, "_AF.tsv")
#     complete_path <- file.path(temp, name)
#     complete_path <- gsub('"', "", gsub("\\\\", "/", complete_path))
#     generate_Vispa2_launch_AF(as_file, project, pool, temp)
#     df <- read.csv(complete_path[1], sep = "\t", header = FALSE)
#     expect_true(ncol(df) == 11)
#     expect_equal(df[, 1], df[, 2])
#     df <- read.csv(complete_path[2], sep = "\t", header = FALSE)
#     expect_true(ncol(df) == 11)
#     expect_equal(df[, 1], df[, 2])
# })

#------------------------------------------------------------------------------#
# Tests as_sparse_matrix
#------------------------------------------------------------------------------#
# smpl <- tibble::tibble(
#     chr = c(1, 2, 3), integration_locus = c(1354, 5634, 4765),
#     strand = c("+", "+", "+"),
#     GeneName = c("GENE1", "GENE2", "GENE3"),
#     GeneStrand = c("+", "+", "+"),
#     CompleteAmplificationID = c("ID1", "ID2", "ID3"),
#     Value = c(46, 546, 587)
# )
#
# smpl_multi <- tibble::tibble(
#     chr = c(1, 2, 3),
#     integration_locus = c(1354, 5634, 4765),
#     strand = c("+", "+", "+"),
#     GeneName = c("GENE1", "GENE2", "GENE3"),
#     GeneStrand = c("+", "+", "+"),
#     CompleteAmplificationID = c("ID1", "ID2", "ID3"),
#     seqCount = c(46, 546, 587),
#     fragmentEstimate = c(4234.5, 533.45, 5431.43),
#     barcodeCount = c(46, 6456, 456)
# )
# test_that("as_sparse_matrix works with single matrix", {
#     sparse <- as_sparse_matrix(smpl)
#     expect_true(tibble::is_tibble(sparse))
#     expect_true(all(c(
#         "chr", "integration_locus", "strand",
#         "GeneName", "GeneStrand", "ID1", "ID2", "ID3"
#     ) %in%
#         colnames(sparse)))
#     expect_equal(nrow(sparse), 3)
# })
#
# test_that("as_sparse_matrix works with multi-quant matrix", {
#     sparse <- as_sparse_matrix(smpl_multi)
#     expect_true(!tibble::is_tibble(sparse) & is.list(sparse))
#     expect_equal(names(sparse), c(
#         "fragmentEstimate",
#         "seqCount",
#         "barcodeCount"
#     ))
#     expect_true(all(c(
#         "chr", "integration_locus", "strand",
#         "GeneName", "GeneStrand", "ID1", "ID2", "ID3"
#     ) %in%
#         colnames(sparse$seqCount)))
#     expect_equal(nrow(sparse$seqCount), 3)
#     expect_true(all(c(
#         "chr", "integration_locus", "strand",
#         "GeneName", "GeneStrand", "ID1", "ID2", "ID3"
#     ) %in%
#         colnames(sparse$fragmentEstimate)))
#     expect_equal(nrow(sparse$fragmentEstimate), 3)
# })
#
# test_that("as_sparse_matrix works with list of matrices", {
#     sparse <- as_sparse_matrix(list(smpl, smpl))
#     expect_true(!tibble::is_tibble(sparse) & is.list(sparse))
#     expect_true(all(c(
#         "chr", "integration_locus", "strand",
#         "GeneName", "GeneStrand", "ID1", "ID2", "ID3"
#     ) %in%
#         colnames(sparse[[1]])))
#     expect_equal(nrow(sparse[[1]]), 3)
#     expect_true(all(c(
#         "chr", "integration_locus", "strand",
#         "GeneName", "GeneStrand", "ID1", "ID2", "ID3"
#     ) %in%
#         colnames(sparse[[2]])))
#     expect_equal(nrow(sparse[[2]]), 3)
# })
#
# #------------------------------------------------------------------------------#
# # Tests annotation_issues
# #------------------------------------------------------------------------------#
# test_df_issues <- tibble::tribble(
#     ~chr, ~integration_locus, ~strand, ~GeneName, ~GeneStrand,
#     ~CompleteAmplificationID, ~Value,
#     "1", 123456, "+", "ABCDE", "+", "ID1", 56,
#     "1", 123456, "+", "ABCDE", "-", "ID2", 675,
#     "1", 123456, "+", "FGHIL", "-", "ID3", 67,
#     "2", 5674653, "-", "FGHIL", "-", "ID2", 873,
#     "1", 4578768, "-", "RSPQX", "-", "ID3", 983,
# )
#
# test_df_no_issues <- tibble::tribble(
#     ~chr, ~integration_locus, ~strand, ~GeneName, ~GeneStrand,
#     ~CompleteAmplificationID, ~Value,
#     "1", 123456, "+", "ABCDE", "+", "ID1", 56,
#     "1", 123456, "+", "ABCDE", "+", "ID2", 675,
#     "1", 123456, "+", "ABCDE", "+", "ID3", 67,
#     "2", 5674653, "-", "FGHIL", "-", "ID2", 873,
#     "1", 4578768, "-", "RSPQX", "-", "ID3", 983,
# )
#
# test_that("annotation_issues returns df if issues", {
#     res <- annotation_issues(test_df_issues)
#     expect_true(!is.null(res))
#     expect_true(nrow(res) == 1 & res$chr[1] == "1" &
#         res$integration_locus[1] == 123456 & res$strand == "+" &
#         res$distinct_genes == 3)
# })
#
# test_that("annotation_issues returns null if no issues", {
#   expect_message({
#     res <- annotation_issues(test_df_no_issues)
#   })
#   expect_null(res)
# })
#
# test_that("annotation_issues works with lists", {
#     res <- annotation_issues(list(a = test_df_issues, b = test_df_no_issues))
#     expect_true(!is.null(res))
#     expect_true(is.null(res[[2]]))
#     expect_true(nrow(res[[1]]) == 1)
# })
