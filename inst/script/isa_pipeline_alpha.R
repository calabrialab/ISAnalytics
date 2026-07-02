#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
    if (is.null(x)) {
        y
    } else {
        x
    }
}

read_isa_pipeline_config <- function(path) {
    env <- new.env(parent = baseenv())
    sys.source(path, envir = env)
    if (!exists("config", envir = env, inherits = FALSE)) {
        stop("Config file must define an object named `config`.", call. = FALSE)
    }
    get("config", envir = env)
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
        plots = file.path(config$run_dir, "plots"),
        logs = file.path(config$run_dir, "logs")
    )
    lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
    dirs
}

log_msg <- function(..., dirs = NULL) {
    msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""))
    message(msg)
    flush.console()
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

save_ggplot <- function(plot, path, width = 9, height = 6) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(
        filename = path,
        plot = plot,
        width = width,
        height = height,
        dpi = 150,
        limitsize = FALSE
    )
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
    cis = "cis",
    plots = "plots"
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
        !identical(stage, "plots")
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

artifact_summary_row <- function(name, object) {
    cls <- paste(class(object), collapse = "/")
    if (is.data.frame(object)) {
        rows <- nrow(object)
        cols <- ncol(object)
        detail <- ""
    } else if (is.null(object)) {
        rows <- 0L
        cols <- 0L
        detail <- "NULL"
    } else if (is.list(object) && name == "cis" && is.data.frame(object$cis)) {
        rows <- nrow(object$cis)
        cols <- ncol(object$cis)
        detail <- "cis$cis"
    } else if (is.list(object) && name == "plots") {
        rows <- length(object)
        cols <- NA_integer_
        detail <- paste(names(object), collapse = ";")
    } else if (is.list(object)) {
        rows <- length(object)
        cols <- NA_integer_
        detail <- paste(names(object), collapse = ";")
    } else {
        rows <- NA_integer_
        cols <- NA_integer_
        detail <- ""
    }
    data.frame(
        artifact = name,
        class = cls,
        rows = rows,
        cols = cols,
        size_mb = round(as.numeric(utils::object.size(object)) / 1024^2, 3),
        detail = detail,
        stringsAsFactors = FALSE
    )
}

dimension_report <- function(dirs) {
    artifact_names <- unique(unlist(stage_defs, use.names = FALSE))
    existing <- artifact_names[file.exists(vapply(
        artifact_names,
        artifact_path,
        character(1),
        dirs = dirs
    ))]
    rows <- lapply(existing, function(name) {
        artifact_summary_row(name, load_artifact(dirs, name))
    })
    summary <- do.call(rbind, rows)
    if (is.null(summary) || nrow(summary) == 0) {
        return(data.frame())
    }

    lineage <- c(
        "matrices_raw",
        "matrices_recalibrated",
        "matrices_no_collisions",
        "data_aggreg",
        "abundance"
    )
    summary$previous_artifact <- NA_character_
    summary$row_delta_from_previous <- NA_integer_
    prev <- NULL
    for (artifact in lineage) {
        idx <- which(summary$artifact == artifact)
        if (length(idx) == 0) {
            next
        }
        if (!is.null(prev)) {
            prev_idx <- which(summary$artifact == prev)
            summary$previous_artifact[[idx]] <- prev
            summary$row_delta_from_previous[[idx]] <-
                summary$rows[[idx]] - summary$rows[[prev_idx]]
        }
        prev <- artifact
    }
    summary
}

write_dimension_report <- function(dirs) {
    summary <- dimension_report(dirs)
    if (nrow(summary) > 0) {
        write_table(summary, file.path(dirs$tables, "pipeline_dimensions.tsv"))
    }
    summary
}

log_dimension_report <- function(summary, dirs) {
    if (nrow(summary) == 0) {
        log_msg("No dimension summary available", dirs = dirs)
        return(invisible(NULL))
    }
    log_msg("Dimension summary:", dirs = dirs)
    for (i in seq_len(nrow(summary))) {
        delta <- summary$row_delta_from_previous[[i]]
        delta_txt <- if (is.na(delta)) {
            ""
        } else {
            paste0(", delta rows vs ", summary$previous_artifact[[i]], ": ", delta)
        }
        log_msg(
            "  - ", summary$artifact[[i]],
            ": ", summary$rows[[i]], " rows x ", summary$cols[[i]], " cols",
            delta_txt,
            dirs = dirs
        )
    }
    invisible(NULL)
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

sanitize_filename <- function(x) {
    x <- gsub("[^A-Za-z0-9._-]+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    ifelse(nzchar(x), x, "plot")
}

compute_diversity_metrics <- function(abundance, config) {
    plot_cfg <- config$plots %||% list()
    value_col <- plot_cfg$diversity_value_col %||%
        "fragmentEstimate_sum_PercAbundance"
    group_key <- plot_cfg$diversity_key %||%
        config$abundance$key %||%
        config$aggregation$key %||%
        c("SubjectID", "CellMarker", "Tissue", "TimePoint")
    group_key <- group_key[group_key %in% colnames(abundance)]

    if (!value_col %in% colnames(abundance) || length(group_key) == 0) {
        return(NULL)
    }

    abundance |>
        dplyr::filter(!is.na(.data[[value_col]]), .data[[value_col]] > 0) |>
        dplyr::group_by(dplyr::across(dplyr::all_of(group_key))) |>
        dplyr::mutate(
            .isa_diversity_p = .data[[value_col]] /
                sum(.data[[value_col]], na.rm = TRUE)
        ) |>
        dplyr::summarise(
            nIS = dplyr::n(),
            shannon = -sum(
                dplyr::if_else(
                    .data$.isa_diversity_p > 0,
                    .data$.isa_diversity_p * log(.data$.isa_diversity_p),
                    0
                ),
                na.rm = TRUE
            ),
            simpson = 1 - sum(.data$.isa_diversity_p^2, na.rm = TRUE),
            inverse_simpson = dplyr::if_else(
                sum(.data$.isa_diversity_p^2, na.rm = TRUE) > 0,
                1 / sum(.data$.isa_diversity_p^2, na.rm = TRUE),
                NA_real_
            ),
            max_abundance = max(.data[[value_col]], na.rm = TRUE),
            .groups = "drop"
        )
}

plot_diversity_metrics <- function(diversity, config) {
    plot_cfg <- config$plots %||% list()
    timepoint_col <- plot_cfg$timepoint_col %||% "TimePoint"
    group_cols <- intersect(
        c("SubjectID", "Tissue", "CellMarker"),
        colnames(diversity)
    )
    diversity$.isa_group_label <- if (length(group_cols) == 0) {
        "all"
    } else {
        do.call(paste, c(diversity[group_cols], sep = " / "))
    }

    if (timepoint_col %in% colnames(diversity)) {
        plot <- ggplot2::ggplot(
            diversity,
            ggplot2::aes(
                x = .data[[timepoint_col]],
                y = .data$shannon,
                color = .data$.isa_group_label,
                group = .data$.isa_group_label
            )
        ) +
            ggplot2::geom_line(linewidth = 0.6, alpha = 0.8) +
            ggplot2::geom_point(size = 2) +
            ggplot2::labs(
                title = "Diversity over time",
                x = timepoint_col,
                y = "Shannon diversity",
                color = NULL
            )
    } else {
        plot <- ggplot2::ggplot(
            diversity,
            ggplot2::aes(
                x = .data$.isa_group_label,
                y = .data$shannon,
                fill = .data$.isa_group_label
            )
        ) +
            ggplot2::geom_col() +
            ggplot2::labs(
                title = "Diversity by group",
                x = NULL,
                y = "Shannon diversity",
                fill = NULL
            )
    }

    if (all(c("Tissue", "CellMarker") %in% colnames(diversity))) {
        plot <- plot +
            ggplot2::facet_grid(stats::as.formula("Tissue ~ CellMarker"))
    }

    plot +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            legend.position = "bottom"
        )
}

save_plot_list <- function(plots, dirs, prefix, max_png = Inf) {
    if (length(plots) == 0) {
        return(invisible(character(0)))
    }
    saved <- character(0)
    plot_names <- names(plots)
    if (is.null(plot_names)) {
        plot_names <- paste0("plot_", seq_along(plots))
    }
    n <- min(length(plots), max_png)
    for (i in seq_len(n)) {
        plot <- plots[[i]]
        if (is.list(plot) && inherits(plot$plot, "ggplot")) {
            plot <- plot$plot
        }
        if (!inherits(plot, "ggplot")) {
            next
        }
        path <- file.path(
            dirs$plots,
            paste0(prefix, "_", sanitize_filename(plot_names[[i]]), ".png")
        )
        save_ggplot(plot, path)
        saved <- c(saved, path)
    }
    invisible(saved)
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
        log_msg(
            "Matrix files: ",
            paste(paste(names(files), basename(unlist(files)), sep = "="), collapse = ", "),
            dirs = dirs
        )
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
        log_msg(
            "Near integration threshold: ", rec_cfg$threshold %||% 4,
            dirs = dirs
        )
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
    } else if (identical(stage, "plots")) {
        plot_cfg <- config$plots %||% list()
        plots <- list()

        if (isTRUE(plot_cfg$cis_volcano %||% TRUE)) {
            cis <- state$cis
            if (is.null(cis) && artifact_exists(dirs, "cis")) {
                cis <- load_artifact(dirs, "cis")
            }
            cis_df <- if (is.list(cis) && is.data.frame(cis$cis)) {
                cis$cis
            } else if (is.data.frame(cis)) {
                cis
            } else {
                NULL
            }
            if (!is.null(cis_df) && nrow(cis_df) > 0) {
                log_msg("Producing CIS volcano plot", dirs = dirs)
                plots$cis_volcano <- ISAnalytics::CIS_volcano_plot(
                    x = cis_df,
                    title_prefix = plot_cfg$cis_title_prefix %||%
                        paste(config$project, config$subproject)
                )
                save_ggplot(
                    plots$cis_volcano,
                    file.path(dirs$plots, "cis_volcano.png")
                )
            } else {
                log_msg("Skipping CIS volcano plot: no CIS data frame available", dirs = dirs)
            }
        }

        if (isTRUE(plot_cfg$sharing_heatmap %||% TRUE)) {
            sharing <- state$sharing
            if (is.null(sharing) && artifact_exists(dirs, "sharing")) {
                sharing <- load_artifact(dirs, "sharing")
            }
            if (is.data.frame(sharing) && nrow(sharing) > 0) {
                log_msg("Producing sharing heatmaps", dirs = dirs)
                plots$sharing_heatmaps <- ISAnalytics::sharing_heatmap(
                    sharing_df = sharing,
                    interactive = FALSE
                )
                save_plot_list(plots$sharing_heatmaps, dirs, "sharing_heatmap")
            } else {
                log_msg("Skipping sharing heatmaps: no sharing data frame available", dirs = dirs)
            }
        }

        if (isTRUE(plot_cfg$alluvial %||% TRUE)) {
            if (!requireNamespace("ggalluvial", quietly = TRUE)) {
                log_msg("Skipping alluvial plots: package `ggalluvial` is not installed", dirs = dirs)
            } else {
                abundance <- state$abundance
                if (is.null(abundance) && artifact_exists(dirs, "abundance")) {
                    abundance <- load_artifact(dirs, "abundance")
                }
                alluvial_plot_y <- plot_cfg$alluvial_plot_y %||%
                    "fragmentEstimate_sum_PercAbundance"
                if (is.data.frame(abundance) &&
                    alluvial_plot_y %in% colnames(abundance)) {
                    log_msg("Producing alluvial plots", dirs = dirs)
                    plots$alluvial <- ISAnalytics::integration_alluvial_plot(
                        x = abundance,
                        group = plot_cfg$alluvial_group %||%
                            c("SubjectID", "CellMarker", "Tissue"),
                        plot_x = plot_cfg$alluvial_plot_x %||% "TimePoint",
                        plot_y = alluvial_plot_y,
                        alluvia_plot_y_threshold =
                            plot_cfg$alluvia_plot_y_threshold %||% 1,
                        top_abundant_tbl = FALSE
                    )
                    save_plot_list(
                        plots$alluvial,
                        dirs,
                        "alluvial",
                        max_png = plot_cfg$alluvial_max_png %||% 8
                    )
                } else {
                    log_msg("Skipping alluvial plots: abundance data or y column is missing", dirs = dirs)
                }
            }
        }

        if (isTRUE(plot_cfg$diversity %||% TRUE)) {
            abundance <- state$abundance
            if (is.null(abundance) && artifact_exists(dirs, "abundance")) {
                abundance <- load_artifact(dirs, "abundance")
            }
            diversity <- if (is.data.frame(abundance)) {
                compute_diversity_metrics(abundance, config)
            } else {
                NULL
            }
            if (!is.null(diversity) && nrow(diversity) > 0) {
                log_msg("Producing diversity summary and plot", dirs = dirs)
                write_table(diversity, file.path(dirs$tables, "diversity_metrics.tsv"))
                plots$diversity_metrics <- diversity
                plots$diversity_plot <- plot_diversity_metrics(diversity, config)
                save_ggplot(
                    plots$diversity_plot,
                    file.path(dirs$plots, "diversity_shannon.png")
                )
            } else {
                log_msg("Skipping diversity plot: required abundance column is missing", dirs = dirs)
            }
        }

        save_artifact(plots, dirs, "plots")
        state$plots <- plots
    } else {
        stop("No runner implemented for stage: ", stage, call. = FALSE)
    }
    log_msg("Completed stage: ", stage, dirs = dirs)
    state
}

run_isa_pipeline_alpha <- function(
        config,
        from = NULL,
        to = NULL,
        resume = NULL,
        config_path = NULL) {
    load_isanalytics(config)
    if (!is.null(config$settings_path)) {
        if (!file.exists(config$settings_path)) {
            stop("Configured settings file does not exist: ", config$settings_path, call. = FALSE)
        }
        ISAnalytics::import_ISA_settings(config$settings_path)
    }

    dirs <- make_dirs(config)
    from <- validate_stage(from %||% config$start_from, "from")
    to <- validate_stage(to %||% config$stop_at, "to")
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
    config_path_saved <- if (is.null(config_path)) {
        NA_character_
    } else {
        normalizePath(config_path, mustWork = FALSE)
    }

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
                config_path = config_path_saved,
                completed = completed,
                skipped = skipped,
                timestamp = Sys.time(),
                artifacts = stage_defs
            ),
            file.path(config$run_dir, "manifest.rds")
        )
    }

    dim_summary <- write_dimension_report(dirs)
    log_dimension_report(dim_summary, dirs)
    log_msg("Pipeline finished", dirs = dirs)
    invisible(list(
        state = state,
        dimensions = dim_summary,
        completed = completed,
        skipped = skipped,
        run_dir = config$run_dir
    ))
}

main <- function() {
    opts <- parse_args(commandArgs(trailingOnly = TRUE))
    config <- read_isa_pipeline_config(opts$config_path)
    run_isa_pipeline_alpha(
        config = config,
        from = opts$from,
        to = opts$to,
        resume = opts$resume,
        config_path = opts$config_path
    )
}

if (sys.nframe() == 0) {
    main()
}
