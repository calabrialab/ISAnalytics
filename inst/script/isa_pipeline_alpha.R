#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
    if (is.null(x)) {
        y
    } else {
        x
    }
}

parse_args <- function(args) {
    if (length(args) < 1) {
        stop(
            "Usage: Rscript inst/script/isa_pipeline_alpha.R CONFIG.R ",
            "[--from=stage] [--to=stage] [--no-resume]",
            call. = FALSE
        )
    }
    opts <- list(config_path = args[[1]], from = NULL, to = NULL, resume = NULL)
    rest <- args[-1]
    for (arg in rest) {
        if (grepl("^--from=", arg)) {
            opts$from <- sub("^--from=", "", arg)
        } else if (grepl("^--to=", arg)) {
            opts$to <- sub("^--to=", "", arg)
        } else if (identical(arg, "--no-resume")) {
            opts$resume <- FALSE
        } else {
            stop("Unknown argument: ", arg, call. = FALSE)
        }
    }
    opts
}

load_config <- function(path) {
    env <- new.env(parent = baseenv())
    sys.source(path, envir = env)
    if (!exists("config", envir = env, inherits = FALSE)) {
        stop("Config file must define an object named `config`.", call. = FALSE)
    }
    get("config", envir = env)
}

load_isanalytics <- function(config) {
    pkg_cfg <- config$package %||% list()
    if (isTRUE(pkg_cfg$load_from_source)) {
        source_dir <- pkg_cfg$source_dir %||% getwd()
        if (!requireNamespace("pkgload", quietly = TRUE)) {
            stop("Package `pkgload` is required to load ISAnalytics from source.")
        }
        pkgload::load_all(source_dir, quiet = TRUE)
    } else {
        suppressPackageStartupMessages(library(ISAnalytics))
    }
}

make_dirs <- function(config) {
    dirs <- list(
        run = config$run_dir,
        artifacts = file.path(config$run_dir, "artifacts"),
        tables = file.path(config$run_dir, "tables"),
        reports = file.path(config$run_dir, "reports"),
        logs = file.path(config$run_dir, "logs")
    )
    lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
    dirs
}

log_msg <- function(..., dirs = NULL) {
    msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""))
    message(msg)
    if (!is.null(dirs)) {
        cat(msg, "\n", file = file.path(dirs$logs, "pipeline.log"), append = TRUE)
    }
}

write_table <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    if (requireNamespace("readr", quietly = TRUE)) {
        readr::write_tsv(x, path)
    } else {
        utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
    }
}

artifact_path <- function(dirs, name) {
    file.path(dirs$artifacts, paste0(name, ".rds"))
}

save_artifact <- function(object, dirs, name) {
    saveRDS(object, artifact_path(dirs, name))
}

load_artifact <- function(dirs, name) {
    readRDS(artifact_path(dirs, name))
}

artifact_exists <- function(dirs, names) {
    all(file.exists(vapply(names, artifact_path, character(1), dirs = dirs)))
}

stage_defs <- list(
    probe = "file_manifest",
    import_matrices = "matrices_raw",
    annotation_qc = "annotation_issues",
    near_integrations = "matrices_recalibrated",
    association = "association_file",
    remove_collisions = "matrices_no_collisions",
    aggregate = c("data_aggreg", "meta_aggreg"),
    abundance = "abundance",
    sharing = "sharing",
    cis = "cis"
)

stage_names <- names(stage_defs)

stage_index <- function(stage) {
    match(stage, stage_names)
}

validate_stage <- function(stage, arg_name) {
    if (is.null(stage)) {
        return(NULL)
    }
    if (!stage %in% stage_names) {
        stop("Unknown ", arg_name, " stage: ", stage, call. = FALSE)
    }
    stage
}

stage_enabled <- function(config, stage) {
    steps <- config$steps %||% list()
    value <- steps[[stage]]
    if (is.null(value)) {
        TRUE
    } else {
        isTRUE(value)
    }
}

load_existing_stage <- function(stage, state, dirs) {
    names <- stage_defs[[stage]]
    if (!artifact_exists(dirs, names)) {
        return(state)
    }
    for (name in names) {
        state[[name]] <- load_artifact(dirs, name)
    }
    state
}

discover_files <- function(config) {
    root <- config$vispa_results_dir
    if (is.null(root) || !dir.exists(root)) {
        stop("VISPA results directory not found: ", root, call. = FALSE)
    }
    paths <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
    info <- file.info(paths)
    data.frame(
        path = paths,
        directory = dirname(paths),
        file = basename(paths),
        size = info$size,
        modified = as.character(info$mtime),
        stringsAsFactors = FALSE
    )
}

matrix_candidates <- function(manifest, config) {
    patterns <- config$matrix_discovery %||% list()
    out <- lapply(names(patterns), function(q) {
        pat <- patterns[[q]]
        hit <- grepl(pat, manifest$file, ignore.case = TRUE) |
            grepl(pat, manifest$path, ignore.case = TRUE)
        if (!any(hit)) {
            return(data.frame())
        }
        cbind(quantification = q, manifest[hit, , drop = FALSE])
    })
    do.call(rbind, out)
}

resolve_matrix_files <- function(config, manifest, dirs) {
    direct <- config$direct_matrix_files %||% list()
    quantifications <- names(config$matrix_discovery %||% direct)
    resolved <- list()
    candidates <- matrix_candidates(manifest, config)
    write_table(candidates, file.path(dirs$tables, "matrix_candidates.tsv"))
    for (q in quantifications) {
        explicit <- direct[[q]]
        if (!is.null(explicit) && nzchar(explicit)) {
            if (!file.exists(explicit)) {
                stop("Configured matrix file for ", q, " does not exist: ", explicit, call. = FALSE)
            }
            resolved[[q]] <- explicit
            next
        }
        q_candidates <- candidates[candidates$quantification == q, , drop = FALSE]
        if (nrow(q_candidates) != 1) {
            stop(
                "Expected exactly one candidate matrix for ", q, ", found ",
                nrow(q_candidates), ". Inspect ",
                file.path(dirs$tables, "matrix_candidates.tsv"),
                " and set `direct_matrix_files$", q, "` in the config.",
                call. = FALSE
            )
        }
        resolved[[q]] <- q_candidates$path[[1]]
    }
    resolved
}

best_matrix <- function(state) {
    state$matrices_no_collisions %||%
        state$matrices_recalibrated %||%
        state$matrices_raw
}

run_stage <- function(stage, state, config, dirs) {
    log_msg("Running stage: ", stage, dirs = dirs)
    if (identical(stage, "probe")) {
        manifest <- discover_files(config)
        save_artifact(manifest, dirs, "file_manifest")
        write_table(manifest, file.path(dirs$tables, "file_manifest.tsv"))
        candidates <- matrix_candidates(manifest, config)
        write_table(candidates, file.path(dirs$tables, "matrix_candidates.tsv"))
        state$file_manifest <- manifest
    } else if (identical(stage, "import_matrices")) {
        manifest <- state$file_manifest %||% load_artifact(dirs, "file_manifest")
        files <- resolve_matrix_files(config, manifest, dirs)
        import_cfg <- config$matrix_import %||% list()
        matrices <- lapply(files, function(path) {
            ISAnalytics::import_single_Vispa2Matrix(
                path = path,
                separator = import_cfg$separator %||% "\t",
                additional_cols = import_cfg$additional_cols,
                sample_names_to = import_cfg$sample_names_to %||% "CompleteAmplificationID",
                values_to = "Value"
            )
        })
        matrices_raw <- ISAnalytics::comparison_matrix(matrices)
        save_artifact(matrices_raw, dirs, "matrices_raw")
        write_table(utils::head(matrices_raw, 100), file.path(dirs$tables, "matrices_raw_head.tsv"))
        state$matrices_raw <- matrices_raw
    } else if (identical(stage, "annotation_qc")) {
        mat <- best_matrix(state)
        issues <- ISAnalytics::annotation_issues(mat)
        save_artifact(issues, dirs, "annotation_issues")
        if (!is.null(issues) && is.data.frame(issues)) {
            write_table(issues, file.path(dirs$tables, "annotation_issues.tsv"))
        }
        state$annotation_issues <- issues
    } else if (identical(stage, "near_integrations")) {
        mat <- best_matrix(state)
        rec_cfg <- config$near_integrations %||% list()
        matrices_recalibrated <- ISAnalytics::compute_near_integrations(
            x = mat,
            threshold = rec_cfg$threshold %||% 4,
            value_columns = rec_cfg$value_columns %||% c("seqCount", "fragmentEstimate"),
            max_value_column = rec_cfg$max_value_column %||% "seqCount",
            map_as_file = rec_cfg$map_as_file %||% TRUE,
            file_path = dirs$reports
        )
        save_artifact(matrices_recalibrated, dirs, "matrices_recalibrated")
        state$matrices_recalibrated <- matrices_recalibrated
    } else if (identical(stage, "association")) {
        if (is.null(config$association_file_path)) {
            stop("`association_file_path` is NULL; cannot run association stage.", call. = FALSE)
        }
        af <- ISAnalytics::import_association_file(
            path = config$association_file_path,
            root = config$vispa_root,
            report_path = dirs$reports
        )
        save_artifact(af, dirs, "association_file")
        state$association_file <- af
    } else if (identical(stage, "remove_collisions")) {
        mat <- best_matrix(state)
        af <- state$association_file %||% load_artifact(dirs, "association_file")
        no_coll <- ISAnalytics::remove_collisions(
            x = mat,
            association_file = af,
            report_path = dirs$reports
        )
        save_artifact(no_coll, dirs, "matrices_no_collisions")
        state$matrices_no_collisions <- no_coll
    } else if (identical(stage, "aggregate")) {
        mat <- best_matrix(state)
        af <- state$association_file %||% load_artifact(dirs, "association_file")
        agg_cfg <- config$aggregation %||% list()
        data_aggreg <- ISAnalytics::aggregate_values_by_key(
            x = mat,
            association_file = af,
            key = agg_cfg$key %||% c("SubjectID", "CellMarker", "Tissue", "TimePoint"),
            value_cols = agg_cfg$value_cols %||% c("seqCount", "fragmentEstimate")
        )
        meta_aggreg <- ISAnalytics::aggregate_metadata(
            association_file = af,
            grouping_keys = agg_cfg$key %||% c("SubjectID", "CellMarker", "Tissue", "TimePoint")
        )
        save_artifact(data_aggreg, dirs, "data_aggreg")
        save_artifact(meta_aggreg, dirs, "meta_aggreg")
        state$data_aggreg <- data_aggreg
        state$meta_aggreg <- meta_aggreg
    } else if (identical(stage, "abundance")) {
        data_aggreg <- state$data_aggreg %||% load_artifact(dirs, "data_aggreg")
        ab_cfg <- config$abundance %||% list()
        abundance <- ISAnalytics::compute_abundance(
            x = data_aggreg,
            columns = ab_cfg$columns %||% c("seqCount_sum", "fragmentEstimate_sum"),
            key = ab_cfg$key %||% c("SubjectID", "CellMarker", "Tissue", "TimePoint")
        )
        save_artifact(abundance, dirs, "abundance")
        state$abundance <- abundance
    } else if (identical(stage, "sharing")) {
        data_aggreg <- state$data_aggreg %||% load_artifact(dirs, "data_aggreg")
        sh_cfg <- config$sharing %||% list()
        sharing <- ISAnalytics::is_sharing(
            x = data_aggreg,
            group_key = sh_cfg$group_key %||% c("SubjectID", "CellMarker", "Tissue", "TimePoint")
        )
        save_artifact(sharing, dirs, "sharing")
        state$sharing <- sharing
    } else if (identical(stage, "cis")) {
        data_aggreg <- state$data_aggreg %||% load_artifact(dirs, "data_aggreg")
        cis_cfg <- config$cis %||% list()
        cis <- ISAnalytics::CIS_grubbs(
            x = data_aggreg,
            genomic_annotation_file = cis_cfg$genomic_annotation_file %||% "hg19"
        )
        save_artifact(cis, dirs, "cis")
        state$cis <- cis
    } else {
        stop("No runner implemented for stage: ", stage, call. = FALSE)
    }
    log_msg("Completed stage: ", stage, dirs = dirs)
    state
}

main <- function() {
    opts <- parse_args(commandArgs(trailingOnly = TRUE))
    config <- load_config(opts$config_path)
    load_isanalytics(config)
    if (!is.null(config$settings_path)) {
        if (!file.exists(config$settings_path)) {
            stop("Configured settings file does not exist: ", config$settings_path, call. = FALSE)
        }
        ISAnalytics::import_ISA_settings(config$settings_path)
    }

    dirs <- make_dirs(config)
    from <- validate_stage(opts$from %||% config$start_from, "from")
    to <- validate_stage(opts$to %||% config$stop_at, "to")
    resume <- opts$resume
    if (is.null(resume)) {
        resume <- isTRUE(config$resume)
    }

    from_idx <- stage_index(from %||% stage_names[[1]])
    to_idx <- stage_index(to %||% stage_names[[length(stage_names)]])
    if (from_idx > to_idx) {
        stop("`from` stage must come before or equal to `to` stage.", call. = FALSE)
    }

    log_msg("Starting ISAnalytics alpha pipeline", dirs = dirs)
    log_msg("Project: ", config$project, " / ", config$subproject, dirs = dirs)
    log_msg("Run dir: ", config$run_dir, dirs = dirs)

    state <- list()
    completed <- character(0)
    skipped <- character(0)

    for (stage in stage_names) {
        idx <- stage_index(stage)
        if (!stage_enabled(config, stage)) {
            skipped <- c(skipped, stage)
            next
        }

        artifacts <- stage_defs[[stage]]
        before_from <- idx < from_idx
        in_range <- idx >= from_idx && idx <= to_idx

        if (before_from) {
            if (artifact_exists(dirs, artifacts)) {
                state <- load_existing_stage(stage, state, dirs)
                log_msg("Loaded prerequisite stage: ", stage, dirs = dirs)
            } else {
                stop(
                    "Stage ", stage, " is before `from` but its artifact is missing. ",
                    "Run from an earlier stage.",
                    call. = FALSE
                )
            }
            next
        }

        if (!in_range) {
            next
        }

        if (resume && is.null(from) && artifact_exists(dirs, artifacts)) {
            state <- load_existing_stage(stage, state, dirs)
            skipped <- c(skipped, stage)
            log_msg("Resumed existing stage: ", stage, dirs = dirs)
            next
        }

        state <- run_stage(stage, state, config, dirs)
        completed <- c(completed, stage)
        saveRDS(
            list(
                project = config$project,
                subproject = config$subproject,
                config_path = normalizePath(opts$config_path, mustWork = FALSE),
                completed = completed,
                skipped = skipped,
                timestamp = Sys.time(),
                artifacts = stage_defs
            ),
            file.path(config$run_dir, "manifest.rds")
        )
    }

    log_msg("Pipeline finished", dirs = dirs)
    invisible(state)
}

main()
