
func_name <- c("HSC_population_size_estimate", ".re_agg_and_estimate")
withr::local_options(list(ISAnalytics.verbose = FALSE))
#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
test_meta <- readRDS(system.file("testdata", "test_population_meta.Rds",
    package = "ISAnalytics"
))

test_data <- readRDS(system.file("testdata", "test_population_iss.Rds",
    package = "ISAnalytics"
))

test_expected <- readRDS(system.file("testdata", "test_population_expected.Rds",
    package = "ISAnalytics"
))

#------------------------------------------------------------------------------#
# Test .re_agg_and_estimate
#------------------------------------------------------------------------------#
test_that(paste(func_name[2], "produces expected output"), {
    meta_with_lineage <- test_meta %>%
        dplyr::filter(.data$NumIS > 5) %>%
        dplyr::left_join(blood_lineages_default(), by = "CellMarker")
    est_slice <- .re_agg_and_estimate(test_data,
        metadata = meta_with_lineage,
        fragmentEstimate_column = NULL,
        seqCount_column = "seqCount_sum",
        tissue_col = "Tissue",
        timepoint_column = "TimePoint",
        aggregation_key = c(
            "SubjectID",
            "CellMarker",
            "Tissue",
            "TimePoint"
        ),
        seqCount_threshold = 3,
        fragmentEstimate_threshold = 3,
        cell_type = "MYELOID",
        tissue_type = "PB",
        annotation_cols = c(
            "GeneName",
            "GeneStrand"
        ),
        subj_col = "SubjectID",
        stable_timepoints = c(9, 12, 13, 18, 30)
    )
    expected <- tibble::tibble(
        Model = c(
            "M0", "Mh Chao (LB)", "Mh Poisson2",
            "Mh Darroch", "Mh Gamma3.5", "M0",
            "Mh Chao (LB)", "Mh Poisson2",
            "Mh Darroch", "Mh Gamma3.5", "M0",
            "M0", "M0", "M0", "M0", "Mth Chao (LB)",
            "Mth Chao (LB)"
        ),
        abundance = c(
            5121.04002304384, 5561.8125,
            5530.54501895936, 6044.43297754501,
            6787.4831931083, 4925.2909486433,
            5343.28561501042, 6005.49169587238,
            7187.70770661454, 9334.62149492534,
            5120.93487419653, 4925.14379335418,
            5493.83500573394, 4061.20124501992,
            4831.28292631282, 4403.02434948464,
            5119.05887032526
        ),
        stderr = c(
            29.6524454059072, 59.5674353240835,
            58.9086290049903, 115.962545291672,
            222.211267809243, 26.925642017495,
            53.181852179626, 125.385406646309,
            296.997069135912, 680.47720113112,
            29.6477624987035, 26.9188940699876,
            189.007141446797, 50.921704530876,
            45.3038415111286, 60.5582122530621,
            42.8131170229558
        ),
        SubjectID = c(
            "PATIENT01", "PATIENT01",
            "PATIENT01", "PATIENT01",
            "PATIENT01", "PATIENT01",
            "PATIENT01", "PATIENT01",
            "PATIENT01", "PATIENT01",
            "PATIENT01", "PATIENT01",
            "PATIENT01", "PATIENT01",
            "PATIENT01", "PATIENT01",
            "PATIENT01"
        ),
        Timepoints = c(
            "All", "All", "All", "All",
            "All", "Stable", "Stable",
            "Stable", "Stable", "Stable",
            "All", "Stable", "Consecutive",
            "Consecutive", "Consecutive",
            "Consecutive", "Consecutive"
        ),
        CellType = c(
            "Myeloid", "Myeloid", "Myeloid",
            "Myeloid", "Myeloid", "Myeloid",
            "Myeloid", "Myeloid", "Myeloid",
            "Myeloid", "Myeloid", "Myeloid",
            "Myeloid", "Myeloid", "Myeloid",
            "Myeloid", "Myeloid"
        ),
        Tissue = c(
            "PB", "PB", "PB", "PB", "PB", "PB",
            "PB", "PB", "PB", "PB", "PB", "PB",
            "PB", "PB", "PB", "PB", "PB"
        ),
        TimePoint_from = c(
            6, 6, 6, 6, 6, 12, 12, 12,
            12, 12, 6, 12, 6, 12, 18, 6,
            12
        ),
        TimePoint_to = c(
            24, 24, 24, 24, 24, 24, 24, 24,
            24, 24, 24, 24, 12, 18, 24, 18,
            24
        ),
        ModelType = c(
            "ClosedPopulation",
            "ClosedPopulation",
            "ClosedPopulation",
            "ClosedPopulation",
            "ClosedPopulation",
            "ClosedPopulation",
            "ClosedPopulation",
            "ClosedPopulation",
            "ClosedPopulation",
            "ClosedPopulation",
            "ClosedPopulation",
            "ClosedPopulation",
            "ClosedPopulation",
            "ClosedPopulation",
            "ClosedPopulation",
            "ClosedPopulation",
            "ClosedPopulation"
        ),
        ModelSetUp = c(
            "models0", "models0",
            "models0", "models0",
            "models0", "models0",
            "models0", "models0",
            "models0", "models0",
            "mthchaobc", "mthchaobc",
            "models0BC", "models0BC",
            "models0BC", "modelMTHBC",
            "modelMTHBC"
        )
    )
    expect_equal(est_slice$est, expected, ignore_attr = TRUE)
})

#------------------------------------------------------------------------------#
# Test HSC_population_size_estimate
#------------------------------------------------------------------------------#
test_that(paste(func_name[1], "produces expected output"), {
    popul <- HSC_population_size_estimate(
        x = test_data,
        metadata = test_meta,
        fragmentEstimate_column = NULL,
        stable_timepoints = c(9, 12, 13, 18, 30)
    )
    expect_equal(popul$est, test_expected)
})


test_that(paste(func_name[1], "produces output missing NumIS"), {
    mod_meta <- test_meta %>% dplyr::select(-dplyr::all_of("NumIS"))
    popul <- HSC_population_size_estimate(
        x = test_data,
        metadata = mod_meta,
        fragmentEstimate_column = NULL,
        stable_timepoints = c(9, 12, 13, 18, 30)
    )
    expect_equal(popul$est, test_expected)
})
