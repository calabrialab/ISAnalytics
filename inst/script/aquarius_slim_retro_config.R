# ISAnalytics pipeline alpha config for aquarius.
#
# Usage from the package root or from a machine where ISAnalytics is installed:
# Rscript inst/script/isa_pipeline_alpha.R inst/script/aquarius_slim_retro_config.R

config <- list(
    project = "SLiM-Retro",
    subproject = "AllPools",

    run_dir = "/storage9/workspace/lomagnoa/ISAnalytics_tests/SLiM-Retro_AllPools_run",
    vispa_results_dir = "/storage9/workspace/lomagnoa/Vispa3/full_runs/SLiM_Retro_PoolSA111/results/SLiM-Retro/AllPools",

    association_file_path = "/storage9/workspace/lomagnoa/Vispa3/full_runs/SLiM_Retro_PoolSA111/input/asso.SLiM-Retro.PoolSA111.tsv",
    vispa_root = NULL,
    settings_path = NULL,

    package = list(
        load_from_source = FALSE,
        source_dir = NULL
    ),

    resume = TRUE,
    start_from = NULL,
    stop_at = NULL,

    steps = list(
        probe = TRUE,
        import_matrices = TRUE,
        annotation_qc = TRUE,
        near_integrations = TRUE,

        association = TRUE,
        remove_collisions = TRUE,
        aggregate = TRUE,
        abundance = TRUE,
        sharing = TRUE,
        cis = TRUE,
        plots = TRUE
    ),

    # If discovery finds more than one file for a quantification, set the
    # corresponding path explicitly here.
    direct_matrix_files = list(
        seqCount = NULL,
        fragmentEstimate = NULL
    ),

    matrix_discovery = list(
        seqCount = "seqCount.*matrix.*annotated.*\\.tsv(\\.gz)?$",
        fragmentEstimate = "fragmentEstimate.*matrix.*annotated.*\\.tsv(\\.gz)?$"
    ),

    matrix_import = list(
        separator = "\t",
        sample_names_to = "CompleteAmplificationID",
        additional_cols = NULL
    ),

    near_integrations = list(
        threshold = 4,
        value_columns = c("seqCount", "fragmentEstimate"),
        max_value_column = "seqCount",
        map_as_file = TRUE
    ),

    collision_removal = list(
        report_path = "reports"
    ),

    aggregation = list(
        key = c("SubjectID", "CellMarker", "Tissue", "TimePoint"),
        value_cols = c("seqCount", "fragmentEstimate")
    ),

    abundance = list(
        columns = c("seqCount_sum", "fragmentEstimate_sum"),
        key = c("SubjectID", "CellMarker", "Tissue", "TimePoint")
    ),

    sharing = list(
        group_key = c("SubjectID", "CellMarker", "Tissue", "TimePoint")
    ),

    cis = list(
        genomic_annotation_file = "hg19"
    ),

    plots = list(
        cis_volcano = TRUE,
        sharing_heatmap = TRUE,
        alluvial = TRUE,
        diversity = TRUE,
        cis_title_prefix = "SLiM-Retro AllPools",
        alluvial_group = c("SubjectID", "CellMarker", "Tissue"),
        alluvial_plot_x = "TimePoint",
        alluvial_plot_y = "fragmentEstimate_sum_PercAbundance",
        alluvia_plot_y_threshold = 1,
        alluvial_max_png = 8,
        diversity_value_col = "fragmentEstimate_sum_PercAbundance",
        timepoint_col = "TimePoint"
    )
)
