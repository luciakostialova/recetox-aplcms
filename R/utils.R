#' @import dplyr tidyr tibble stringr arrow
NULL
#> NULL

register_functions_to_cluster <- function(cluster) {
    snow::clusterExport(cluster, list(
        'adaptive.bin',
        'add_feature_ids',
        'adjust.time',
        'aggregate_by_rt',
        'as_feature_crosstab',
        'as_feature_sample_table',
        'as_wide_aligned_table',
        'augment_known_table',
        'bigauss.esti',
        'bigauss.esti.EM',
        'bigauss.mix',
        'bind_batch_label_column',
        'characterize',
        'check_files',
        'clean_data_matrix',
        'comb',
        'compute_boundaries',
        'compute_bounds',
        'compute_breaks',
        'compute_breaks_3',
        'compute_chromatographic_profile',
        'compute_clusters',
        'compute_clusters_simple',
        'compute_comb',
        'compute_corrected_features',
        'compute_corrected_features_v2',
        'compute_curr_rec_with_enough_peaks',
        'compute_delta_rt',
        'compute_densities',
        'compute_dx',
        'compute_e_step',
        'compute_end_bound',
        'compute_gaussian_peak_shape',
        'compute_initiation_params',
        'compute_intensity_medians',
        'compute_mass_density',
        'compute_mass_values',
        'compute_min_mz_tolerance',
        'compute_mu_sc_std',
        'compute_mz_sd',
        'compute_peaks_and_valleys',
        'compute_pks_vlys_rt',
        'compute_rectangle',
        'compute_rt_intervals_indices',
        'compute_scale',
        'compute_sel',
        'compute_start_bound',
        'compute_target_times',
        'compute_template',
        'compute_template_adjusted_rt',
        'compute_uniq_grp',
        'concatenate_feature_tables',
        'correct_time',
        'correct_time_v2',
        'count_peaks',
        'create_aligned_feature_table',
        'create_features_from_cluster',
        'create_intensity_row',
        'create_metadata',
        'create_output',
        'create_rt_row',
        'draw_plot',
        'draw_rt_correction_plot',
        'draw_rt_normal_peaks',
        'duplicate.row.remove',
        'enrich_table_by_known_features',
        'extract_pattern_colnames',
        'feature_recovery',
        'fill_missing_values',
        'filter_based_on_density',
        'filter_features_by_presence',
        'find.match',
        'find.tol.time',
        'find.turn.point',
        'find_local_maxima',
        'find_min_position',
        'find_mz_match',
        'find_mz_tolerance',
        'find_optima',
        'get_custom_rt_tol',
        'get_features_in_rt_range',
        'get_mzrange_bound_indices',
        'get_num_workers',
        'get_rt_region_indices',
        'get_sample_name',
        'get_single_occurrence_mask',
        'get_times_to_use',
        'hybrid',
        'increment_counter',
        'interpol.area',
        'l2normalize',
        'label_val_to_keep',
        'load.lcms',
        'load.lcms.raw',
        'load_aligned_features',
        'load_data',
        'load_file',
        'long_to_wide_feature_table',
        'match_peaks',
        'merge_features_and_known_table',
        'merge_known_tables',
        'msExtrema',
        'normix',
        'normix.bic',
        'pivot_feature_values',
        'plot_normix_bic',
        'plot_peak_summary',
        'plot_raw_profile_histogram',
        'plot_rt_profile',
        'predict_mz_break_indices',
        'predict_smoothed_rt',
        'prep_uv',
        'preprocess_bandwidth',
        'preprocess_profile',
        'process_chunk',
        'prof.to.features',
        'read_parquet_files',
        'readjust_times',
        'recover.weaker',
        'recover_weaker_signals',
        'refine_selection',
        'remove_noise',
        'rev_cum_sum',
        'rm.ridge',
        'run_filter',
        'semi.sup',
        'semisup_to_hybrid_adapter',
        'solve_a',
        'solve_sigma',
        'sort_data',
        'span',
        'tolerance_plot',
        'two.step.hybrid',
        'unsupervised',
        'validate_contents',
        'validate_model_method_input',
        'wide_to_long_feature_table'
    ))
    snow::clusterEvalQ(cluster, library("dplyr"))
    snow::clusterEvalQ(cluster, library("stringr"))
    snow::clusterEvalQ(cluster, library("tidyselect"))
}

#' Concatenate multiple feature lists and add the sample id (origin of feature) as additional column.
#' 
#' @param features list List of tibbles containing extracted feature tables.
#' @export
concatenate_feature_tables <- function(features, sample_names) {
    for (i in seq_along(features)) {
        if(!("sample_id" %in% colnames(features[[i]]))) {
            features[[i]] <- tibble::add_column(features[[i]], sample_id = sample_names[i])
        }
    }
    
    merged <- dplyr::bind_rows(features)
    return(merged)
}

#' @export
load_aligned_features <- function(metadata_file, intensities_file, rt_file) {
    metadata <- arrow::read_parquet(metadata_file)
    intensities <- arrow::read_parquet(intensities_file)
    rt <- arrow::read_parquet(rt_file)
    
    result <- list()
    result$metadata <- as_tibble(metadata)
    result$intensity <- as_tibble(intensities)
    result$rt <- as_tibble(rt)
    return(result)
}

#' Calculate the span of a numeric vector.
#' @description
#' This function calculates the span (range) of a numeric vector, ignoring NA values.
#' @param x A numeric vector.
#' @return A numeric value representing the span of the vector.
#' @export
span <- function(x) {
    diff(range(x, na.rm = TRUE))
}

#' Compute standard deviation of m/z values for feature groups.
#' @description
#' This function computes the standard deviation of m/z values for each group of features.
#' @param feature_groups A list of data frames, where each data frame represents a group of features with m/z values.
#' @return A numeric vector of standard deviations of m/z values for each feature group.
#' @export
compute_mz_sd <- function(feature_groups) {
    mz_sd <- c()
    for (i in seq_along(feature_groups)) {
        group <- feature_groups[[i]]
        
        if (nrow(group > 1)) {
            group_sd <- sd(group[, "mz"])
            mz_sd <- append(mz_sd, group_sd)
        }
    }
    return(mz_sd)
}

#' Get the number of available worker cores.
#' @description
#' This function determines the number of available worker cores, taking into account CRAN's limit on the number of cores.
#' @return An integer representing the number of available worker cores.
#' @export
get_num_workers <- function() {
    # CRAN limits the number of cores available to packages to 2
    # source https://stackoverflow.com/questions/50571325/r-cran-check-fail-when-using-parallel-functions#50571533
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    
    if (nzchar(chk) && chk == "TRUE") {
        # use 2 cores in CRAN/Travis/AppVeyor
        num_workers <- 2L
    } else {
        # use all cores in devtools::test()
        num_workers <- parallel::detectCores()
    }
    return(num_workers)
}

#' @export
read_parquet_files <- function(filename, folder, pattern) {
  testdata <- file.path("..", "testdata")
  
  input <- lapply(filename, function(x) {
    tibble::as_tibble(arrow::read_parquet(file.path(testdata, folder, paste0(x, pattern))))
  })

  return(input)
}
