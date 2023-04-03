#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
withr::local_options(
    ISAnalytics.reports = FALSE,
    ISAnalytics.verbose = FALSE
)

set.seed(125)
normal_vect <- extraDistr::rtnorm(
    n = 5000, mean = 2500, sd = 1,
    a = 1, b = 5000
)

min_example <- tibble::tibble(
    CompleteAmplificationID = paste0(
        "ID",
        seq_len(1000)
    ),
    PoolID = paste0(
        "POOL",
        sample(c(1, 2, 3, 4), 1000,
            replace = TRUE,
            prob = c(
                0.1, 0.3,
                0.3, 0.3
            )
        )
    ),
    X = extraDistr::rtnorm(
        n = 1000, mean = 10000,
        sd = 1, a = 1, b = 20000
    ),
    Y = runif(1000, min = -2, max = 1000),
    Z = sample(c(NA, seq_len(4000)), 1000,
        replace = TRUE
    )
)

#------------------------------------------------------------------------------#
# .outlier_test_verify_logiop
#------------------------------------------------------------------------------#
test_that(".outlier_test_verify_logiop warns when flag_logic too long", {
    withr::local_options(list(ISAnalytics.verbose = TRUE))
    key <- c("A", "B")
    flag_logic <- c("AND", "OR")
    expect_message(
        {
            .outlier_test_verify_logiop(key, flag_logic, "flag_logic")
        },
        class = "flag_logic_long"
    )
    expect_true(length(flag_logic) == 1)
})

test_that(".outlier_test_verify_logiop warns flags short but not 1", {
    withr::local_options(list(ISAnalytics.verbose = TRUE))
    key <- c("A", "B", "C", "D")
    flag_logic <- c("AND", "OR")
    expect_message(
        {
            .outlier_test_verify_logiop(key, flag_logic, "flag_logic")
        },
        class = "flag_logic_short"
    )
    expect_true(length(flag_logic) == 1)
})

test_that(".outlier_test_verify_logiop errors if unsupported ops", {
    withr::local_options(list(ISAnalytics.verbose = TRUE))
    key <- c("A", "B", "C", "D")
    flag_logic <- c("AND", "OR", "X")
    expect_error(
        {
            .outlier_test_verify_logiop(key, flag_logic, "flag_logic")
        },
        class = "unsupp_logi_op"
    )
})

#------------------------------------------------------------------------------#
# .tests_pool_frag
#------------------------------------------------------------------------------#
test_that(".tests_pool_frag - normal dist with normality test", {
    calc <- .tests_pool_frag(normal_vect,
        suffix = "test",
        normality_test = TRUE, normality_threshold = 0.05,
        log2 = FALSE
    )
    expect_true("normality_test" %in% colnames(calc))
    expect_true(all(!is.na(calc$normality_test)) &
        all(calc$normality_test == TRUE))
    expect_true(all(!is.na(calc$zscore_test)))
    expect_true(all(!is.na(calc$tstudent_test)))
    expect_true(all(!is.na(calc$tdist_test)))
})

test_that(".tests_pool_frag - skips normality test with no error", {
    set.seed(125)
    x <- extraDistr::rtnorm(
        n = 6000, mean = 3000, sd = 1,
        a = 1, b = 6000
    )
    expect_message({
        expect_message({
            calc <- .tests_pool_frag(x,
                suffix = "test",
                normality_test = TRUE, normality_threshold = 0.05,
                log2 = FALSE
            )
        })
    })
    expect_true(all(is.na(calc)))
})

test_that(".tests_pool_frag - performs log2 no norm test", {
    set.seed(125)
    x <- extraDistr::rtnorm(
        n = 6000, mean = 3000, sd = 1,
        a = 1, b = 6000
    )
    calc <- .tests_pool_frag(x,
        suffix = "test",
        normality_test = FALSE,
        normality_threshold = 0.05,
        log2 = TRUE
    )
    expect_true(all(!is.na(calc$log2_test)))
})

test_that(".tests_pool_frag - performs log2 with normality test", {
    calc <- .tests_pool_frag(normal_vect,
        suffix = "test",
        normality_test = TRUE,
        normality_threshold = 0.05,
        log2 = TRUE
    )
    expect_true("normality_test" %in% colnames(calc))
    expect_true(all(!is.na(calc$normality_test)) &
        all(calc$normality_test == TRUE))
    expect_true(all(!is.na(calc$zscore_test)))
    expect_true(all(!is.na(calc$tstudent_test)))
    expect_true(all(!is.na(calc$tdist_test)))
    expect_true(all(!is.na(calc$log2_test)))
})

#------------------------------------------------------------------------------#
# .flag_cond_outliers_pool_frag
#------------------------------------------------------------------------------#
test_that(".flag_cond_outliers_pool_frag flags false ok", {
    ## If processed is false flag must be false in any case
    flag <- .flag_cond_outliers_pool_frag(FALSE, 0, 0, 0)
    expect_false(flag)
    ## tdist not less than threshold
    flag <- .flag_cond_outliers_pool_frag(
        proc = TRUE,
        tdist = 1.5,
        zscore = -2,
        outlier_threshold = 0.05
    )
    expect_false(flag)
    flag <- .flag_cond_outliers_pool_frag(
        proc = TRUE,
        tdist = 0.05,
        zscore = -2,
        outlier_threshold = 0.05
    )
    expect_false(flag)
    ## zscore not less than 0
    flag <- .flag_cond_outliers_pool_frag(
        proc = TRUE,
        tdist = 0.01,
        zscore = 1,
        outlier_threshold = 0.05
    )
    expect_false(flag)
    flag <- .flag_cond_outliers_pool_frag(
        proc = TRUE,
        tdist = 0.01,
        zscore = 0,
        outlier_threshold = 0.05
    )
    expect_false(flag)
})

test_that(".flag_cond_outliers_pool_frag flags true ok", {
    flag <- .flag_cond_outliers_pool_frag(
        proc = TRUE,
        tdist = 0.01,
        zscore = -1,
        outlier_threshold = 0.05
    )
    expect_true(flag)
})

#------------------------------------------------------------------------------#
# .apply_flag_logic
#------------------------------------------------------------------------------#
test_that(".apply_flag_logic correct AND logic", {
    flag <- .apply_flag_logic(x = TRUE, y = FALSE, logic = "AND")
    expect_false(flag)
    flag <- .apply_flag_logic(x = FALSE, y = FALSE, logic = "AND")
    expect_false(flag)
    flag <- .apply_flag_logic(x = TRUE, y = TRUE, logic = "AND")
    expect_true(flag)
})

test_that(".apply_flag_logic correct OR logic", {
    flag <- .apply_flag_logic(x = TRUE, y = FALSE, logic = "OR")
    expect_true(flag)
    flag <- .apply_flag_logic(x = FALSE, y = FALSE, logic = "OR")
    expect_false(flag)
    flag <- .apply_flag_logic(x = TRUE, y = TRUE, logic = "OR")
    expect_true(flag)
})

test_that(".apply_flag_logic correct NOR logic", {
    flag <- .apply_flag_logic(x = TRUE, y = FALSE, logic = "NOR")
    expect_false(flag)
    flag <- .apply_flag_logic(x = FALSE, y = FALSE, logic = "NOR")
    expect_true(flag)
    flag <- .apply_flag_logic(x = TRUE, y = TRUE, logic = "NOR")
    expect_false(flag)
})

test_that(".apply_flag_logic correct XOR logic", {
    flag <- .apply_flag_logic(x = TRUE, y = FALSE, logic = "XOR")
    expect_true(flag)
    flag <- .apply_flag_logic(x = FALSE, y = FALSE, logic = "XOR")
    expect_false(flag)
    flag <- .apply_flag_logic(x = TRUE, y = TRUE, logic = "XOR")
    expect_false(flag)
})

test_that(".apply_flag_logic correct XNOR logic", {
    flag <- .apply_flag_logic(x = TRUE, y = FALSE, logic = "XNOR")
    expect_false(flag)
    flag <- .apply_flag_logic(x = FALSE, y = FALSE, logic = "XNOR")
    expect_true(flag)
    flag <- .apply_flag_logic(x = TRUE, y = TRUE, logic = "XNOR")
    expect_true(flag)
})

test_that(".apply_flag_logic correct NAND logic", {
    flag <- .apply_flag_logic(x = TRUE, y = FALSE, logic = "NAND")
    expect_true(flag)
    flag <- .apply_flag_logic(x = FALSE, y = FALSE, logic = "NAND")
    expect_true(flag)
    flag <- .apply_flag_logic(x = TRUE, y = TRUE, logic = "NAND")
    expect_false(flag)
})

test_that(".apply_flag_logic correct combined logic AND cross", {
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("AND", "OR")
    )
    expect_false(flag)
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("AND", "XOR")
    )
    expect_false(flag)
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("AND", "NOR")
    )
    expect_true(flag)
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("AND", "NAND")
    )
    expect_true(flag)
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("AND", "XNOR")
    )
    expect_true(flag)
})

test_that(".apply_flag_logic correct combined logic OR cross", {
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("OR", "AND")
    )
    expect_false(flag)
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("OR", "XOR")
    )
    expect_true(flag)
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("OR", "NOR")
    )
    expect_false(flag)
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("OR", "NAND")
    )
    expect_true(flag)
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("OR", "XNOR")
    )
    expect_false(flag)
})

test_that(".apply_flag_logic correct combined logic XOR cross", {
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("XOR", "AND")
    )
    expect_false(flag)
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("XOR", "OR")
    )
    expect_true(flag)
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("XOR", "NOR")
    )
    expect_false(flag)
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("XOR", "NAND")
    )
    expect_true(flag)
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("XOR", "XNOR")
    )
    expect_false(flag)
})

test_that(".apply_flag_logic correct combined logic NOR cross", {
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("NOR", "AND")
    )
    expect_false(flag)
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("NOR", "OR")
    )
    expect_false(flag)
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("NOR", "XOR")
    )
    expect_false(flag)
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("NOR", "NAND")
    )
    expect_true(flag)
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("NOR", "XNOR")
    )
    expect_true(flag)
})

test_that(".apply_flag_logic correct combined logic NAND cross", {
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("NAND", "AND")
    )
    expect_false(flag)
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("NAND", "OR")
    )
    expect_true(flag)
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("NAND", "XOR")
    )
    expect_true(flag)
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("NAND", "NOR")
    )
    expect_false(flag)
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("NAND", "XNOR")
    )
    expect_false(flag)
})

test_that(".apply_flag_logic correct combined logic XNOR cross", {
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("XNOR", "AND")
    )
    expect_false(flag)
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("XNOR", "OR")
    )
    expect_false(flag)
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("XNOR", "XOR")
    )
    expect_false(flag)
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("XNOR", "NOR")
    )
    expect_true(flag)
    flag <- .apply_flag_logic(
        x = TRUE, y = FALSE,
        k = FALSE, logic = c("XNOR", "NAND")
    )
    expect_true(flag)
})

#------------------------------------------------------------------------------#
# .process_pool_frag
#------------------------------------------------------------------------------#
test_that(".process_pool_frag - single key norm test", {
    res <- .process_pool_frag(min_example,
        key = "X", normality_test = TRUE,
        normality_threshold = 0.05, log2 = FALSE
    )
    expect_true(all(res$normality_X))
    expect_true(all(!is.na(res$zscore_X)))
    expect_true(all(!is.na(res$tdist_X)))
    res <- .process_pool_frag(min_example,
        key = "X", normality_test = TRUE,
        normality_threshold = 0.05, log2 = TRUE
    )
    expect_true(all(res$normality_X))
    expect_true(all(!is.na(res$zscore_X)))
    expect_true(all(!is.na(res$tdist_X)))

    res <- .process_pool_frag(min_example,
        key = "Y", normality_test = TRUE,
        normality_threshold = 0.05, log2 = FALSE
    )
    expect_true(all(!res$normality_Y))
    expect_true(all(is.na(res$zscore_Y)))
    expect_true(all(is.na(res$tdist_Y)))
})

#------------------------------------------------------------------------------#
# .pool_frag_calc
#------------------------------------------------------------------------------#
test_that(".pool_frag_calc - single key", {
    samples_per_pool <- 5
    samples_in_pools <- min_example |>
        dplyr::group_by(.data$PoolID) |>
        dplyr::summarise(count = dplyr::n())
    ## ON X
    res <- .pool_frag_calc(min_example,
        key = "X", by_pool = TRUE,
        min_samples_per_pool = samples_per_pool,
        pool_col = "PoolID",
        normality_test = FALSE,
        normality_threshold = 0.05,
        log2 = FALSE,
        pcr_id_col = "CompleteAmplificationID"
    )
    processed_pools <- purrr::pmap_lgl(samples_in_pools, function(...) {
        pool <- tibble::tibble(...)
        should_be <- if (pool$count < samples_per_pool) {
            FALSE
        } else {
            TRUE
        }
        all(res$metadata |> dplyr::filter(.data$PoolID == pool$PoolID) |>
            dplyr::pull(.data$processed) == should_be)
    })
    non_na_if_proc <- purrr::pmap_lgl(res$metadata, function(...) {
        row <- tibble::tibble(...)
        if (row$processed == FALSE) {
            return(is.na(row$tstudent_X) & is.na(row$tdist_X))
        } else {
            return(!is.na(row$tstudent_X) & !is.na(row$tdist_X))
        }
    })
    expect_true(all(processed_pools))
    expect_true(all(non_na_if_proc))
    expect_true(all(is.na(res$normality_X)))
    res <- .pool_frag_calc(min_example,
        key = "X", by_pool = FALSE,
        min_samples_per_pool = samples_per_pool,
        pool_col = "PoolID",
        normality_test = FALSE,
        normality_threshold = 0.05,
        log2 = FALSE, pcr_id_col = "CompleteAmplificationID"
    )
    non_na_if_proc <- purrr::pmap_lgl(res$metadata, function(...) {
        row <- tibble::tibble(...)
        if (row$processed == FALSE) {
            return(is.na(row$tstudent_X) & is.na(row$tdist_X))
        } else {
            return(!is.na(row$tstudent_X) & !is.na(row$tdist_X))
        }
    })
    expect_true(all(non_na_if_proc))
    expect_true(all(is.na(res$normality_X)))
    ## ON Y
    identified_zeros <- min_example |>
        dplyr::filter(.data$Y <= 0)
    res <- .pool_frag_calc(min_example,
        key = "Y", by_pool = TRUE,
        min_samples_per_pool = samples_per_pool,
        pool_col = "PoolID",
        normality_test = FALSE,
        normality_threshold = 0.05,
        log2 = TRUE, pcr_id_col = "CompleteAmplificationID"
    )
    samples_in_pools <- min_example |>
        dplyr::filter(!.data$CompleteAmplificationID %in%
            identified_zeros$CompleteAmplificationID) |>
        dplyr::group_by(.data$PoolID) |>
        dplyr::summarise(count = dplyr::n())
    processed_pools <- purrr::pmap_lgl(samples_in_pools, function(...) {
        pool <- tibble::tibble(...)
        should_be <- if (pool$count < samples_per_pool) {
            FALSE
        } else {
            TRUE
        }
        all(res$metadata |>
            dplyr::filter(
                !.data$CompleteAmplificationID %in%
                    identified_zeros$CompleteAmplificationID,
                .data$PoolID == pool$PoolID
            ) |>
            dplyr::pull(.data$processed) == should_be)
    })
    non_na_if_proc <- purrr::pmap_lgl(res$metadata, function(...) {
        row <- tibble::tibble(...)
        if (row$processed == FALSE) {
            return(is.na(row$tstudent_Y) & is.na(row$tdist_Y))
        } else {
            return(!is.na(row$tstudent_Y) & !is.na(row$tdist_Y))
        }
    })
    expect_true(all(processed_pools))
    expect_true(all(non_na_if_proc))
    expect_true(all(is.na(res$normality_Y)))
    expect_true(all((res$metadata |>
        dplyr::filter(.data$CompleteAmplificationID %in%
            identified_zeros$CompleteAmplificationID) |>
        dplyr::pull(.data$processed)
    ) == FALSE))
})

#------------------------------------------------------------------------------#
# outliers_by_pool_fragments
#------------------------------------------------------------------------------#
test_that("outliers_by_pool_fragments - processes with nas", {
    nas <- min_example |>
        dplyr::filter(is.na(.data$Z))
    res <- outliers_by_pool_fragments(min_example,
        key = "Z",
        keep_calc_cols = TRUE
    )
    expect_true(all(res |>
        dplyr::filter(.data$CompleteAmplificationID %in%
            nas$CompleteAmplificationID) |>
        dplyr::pull(.data$processed) == FALSE))
    expect_true(all(res |>
        dplyr::filter(.data$CompleteAmplificationID %in%
            nas$CompleteAmplificationID) |>
        dplyr::pull(.data$to_remove) == FALSE))
    res <- outliers_by_pool_fragments(min_example,
        key = c("X", "Z"),
        keep_calc_cols = TRUE
    )
    expect_true(all(res |>
        dplyr::filter(.data$CompleteAmplificationID %in%
            nas$CompleteAmplificationID) |>
        dplyr::pull(.data$processed) == FALSE))
    expect_true(all(res |>
        dplyr::filter(.data$CompleteAmplificationID %in%
            nas$CompleteAmplificationID) |>
        dplyr::pull(.data$to_remove) == FALSE))
})

test_that("outliers_by_pool_fragments - produces report", {
    withr::local_options(list(ISAnalytics.reports = TRUE))
    withr::local_options(list(ISAnalytics.reports = TRUE))
    tmp_dir <- withr::local_tempdir()
    res <- outliers_by_pool_fragments(min_example,
        key = c("X", "Z"),
        keep_calc_cols = TRUE,
        report_path = tmp_dir
    )
    path_to_file <- fs::path(tmp_dir, .generate_report_filename("outlier_flag"))
    expect_true(fs::file_exists(path_to_file))
})

#------------------------------------------------------------------------------#
# outlier_filter
#------------------------------------------------------------------------------#
test_that("outlier_filter - works with calls to functions", {
    example <- tibble::tribble(
        ~CompleteAmplificationID,
        ~PoolID, ~X, ~Y, ~Z,
        "ID1", "POOL2", 10000.143, 160.42487, NA,
        "ID2", "POOL1", 10000.257, 214.61753, 363,
        "ID3", "POOL3", 10000.716, 931.20962, 956,
        "ID4", "POOL2", 10001.609, 211.75091, NA,
        "ID5", "POOL1", 10000.667, 413.37500, 541,
        "ID6", "POOL1", 9999.257, 33.85144, 398,
        "ID7", "POOL1", 10000.065, 247.01316, 1549,
        "ID8", "POOL2", 9999.991, -610.04281, 578,
        "ID9", "POOL3", 9999.359, 385.85765, 2520,
        "ID10", "POOL2", 9999.490, 955.12338, 2032,
        "ID11", "POOL3", 10001.566, 888.10511, 2692,
        "ID12", "POOL3", 10000.270, 451.15594, 3542,
        "ID13", "POOL1", 9999.588, 683.70609, 3477,
        "ID14", "POOL3", 10001.123, 262.19671, 1891,
        "ID15", "POOL3", 10000.218, 698.22728, 3133,
        "ID16", "POOL3", 9997.922, 185.25041, 2699,
        "ID17", "POOL2", 10001.322, 747.29494, 2316,
        "ID18", "POOL3", 9999.698, 889.56291, 3751,
        "ID19", "POOL1", 10001.496, 484.01714, 542,
        "ID20", "POOL3", 9998.491, 277.64779, 641
    )
    res <- outlier_filter(example,
        key = "Z",
        keep_calc_cols = TRUE,
        outlier_p_value_threshold = 0.05
    )
    expected_flagged <- c("ID20")
    expect_true(all(!expected_flagged %in% res$CompleteAmplificationID))

    ## With custom function
    foo <- function(metadata) {
        processed <- metadata
        processed <- processed |>
            dplyr::mutate(
                to_remove = dplyr::if_else(
                    .data$X < 10000, TRUE, FALSE
                )
            )
        return(processed[, c("CompleteAmplificationID", "to_remove")])
    }
    res <- outlier_filter(
        metadata = example,
        outlier_test = c(outliers_by_pool_fragments, foo),
        key = "Z",
        outlier_p_value_threshold = 0.05
    )
    expected_flagged <- c("ID20")
    expect_true(all(!expected_flagged %in% res$CompleteAmplificationID))
    res <- outlier_filter(
        metadata = example,
        outlier_test = c(outliers_by_pool_fragments, foo),
        key = "Z",
        outlier_p_value_threshold = 0.05,
        combination_logic = c("OR")
    )
    expected_flagged <- c(
        "ID20", "ID6", "ID13", "ID9", "ID16", "ID18", "ID8",
        "ID10"
    )
    expect_true(all(!expected_flagged %in% res$CompleteAmplificationID))
})

test_that("outlier_filter - works with tests outputs", {
    example <- tibble::tribble(
        ~CompleteAmplificationID,
        ~PoolID, ~X, ~Y, ~Z,
        "ID1", "POOL2", 10000.143, 160.42487, NA,
        "ID2", "POOL1", 10000.257, 214.61753, 363,
        "ID3", "POOL3", 10000.716, 931.20962, 956,
        "ID4", "POOL2", 10001.609, 211.75091, NA,
        "ID5", "POOL1", 10000.667, 413.37500, 541,
        "ID6", "POOL1", 9999.257, 33.85144, 398,
        "ID7", "POOL1", 10000.065, 247.01316, 1549,
        "ID8", "POOL2", 9999.991, -610.04281, 578,
        "ID9", "POOL3", 9999.359, 385.85765, 2520,
        "ID10", "POOL2", 9999.490, 955.12338, 2032
    )
    test_result_1 <- tibble::tribble(
        ~CompleteAmplificationID, ~to_remove,
        "ID1", TRUE,
        "ID2", FALSE,
        "ID3", FALSE,
        "ID4", FALSE,
        "ID5", TRUE,
        "ID6", TRUE,
        "ID7", FALSE,
        "ID8", FALSE,
        "ID9", FALSE,
        "ID10", FALSE
    )
    test_result_2 <- tibble::tribble(
        ~CompleteAmplificationID, ~to_remove,
        "ID1", TRUE,
        "ID2", FALSE,
        "ID3", TRUE,
        "ID4", FALSE,
        "ID5", FALSE,
        "ID6", TRUE,
        "ID7", FALSE,
        "ID8", FALSE,
        "ID9", FALSE,
        "ID10", TRUE
    )
    res <- outlier_filter(
        metadata = example,
        outlier_test_outputs = test_result_1
    )
    expected_flagged <- c("ID1", "ID5", "ID6")
    expect_true(all(!expected_flagged %in% res$CompleteAmplificationID))

    res <- outlier_filter(
        metadata = example,
        outlier_test_outputs = list(
            test_result_1,
            test_result_2
        ),
        combination_logic = c("AND")
    )
    expected_flagged <- c("ID1", "ID6")
    expect_true(all(!expected_flagged %in% res$CompleteAmplificationID))
})

test_that("outlier_filter - fails with wrong output format", {
    example <- tibble::tribble(
        ~CompleteAmplificationID,
        ~PoolID, ~X, ~Y, ~Z,
        "ID1", "POOL2", 10000.143, 160.42487, NA,
        "ID2", "POOL1", 10000.257, 214.61753, 363,
        "ID3", "POOL3", 10000.716, 931.20962, 956,
        "ID4", "POOL2", 10001.609, 211.75091, NA
    )
    ## Missing IDs
    test_result_1 <- tibble::tribble(
        ~CompleteAmplificationID, ~to_remove,
        "ID1", TRUE,
        "ID3", FALSE,
        "ID4", FALSE
    )
    ## Missing columns
    test_result_2 <- tibble::tribble(
        ~SampleID, ~remove,
        "ID1", TRUE,
        "ID2", FALSE,
        "ID3", FALSE,
        "ID4", FALSE
    )
    expect_error(
        {
            res <- outlier_filter(
                metadata = example,
                outlier_test_outputs = test_result_1
            )
        },
        class = "outlier_format_err"
    )
    expect_error(
        {
            res <- outlier_filter(
                metadata = example,
                outlier_test_outputs = test_result_2
            )
        },
        class = "outlier_format_err"
    )
    ## Additional IDs should not raise errors but should not appear in final
    ## output
    test_result_3 <- tibble::tribble(
        ~CompleteAmplificationID, ~to_remove,
        "ID1", TRUE,
        "ID2", FALSE,
        "ID3", FALSE,
        "ID4", FALSE,
        "ID5", FALSE
    )
    res <- outlier_filter(
        metadata = example,
        outlier_test_outputs = test_result_3
    )
    expect_true(all(c("ID2", "ID3", "ID4") %in%
        res$CompleteAmplificationID))
    expect_true(all(!c("ID1", "ID5") %in% res$CompleteAmplificationID))
})
