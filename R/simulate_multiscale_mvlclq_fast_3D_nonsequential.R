#' function to perform single-image, multi-scale 3D permutation
#'
#' @param sce_obj SpatialExperiment or SingleCellExperiment object with one image
#' @param spatial_weights_multi spatial weights for the image
#' @param type_column column name for cell type
#' @param point_type reference cell type
#' @param target_types query cell types
#' @param observed_mvlclq_multi observed mlCLQ scores calculated with calculate_mvlclq_multi_image
#' @param occurrence_probs_multi occurrence probability calculated with calculate_occurrence_probability_multi_image
#' @param image_column Column name for image index. NULL for direct inputs, or the column containing image ID "1" when the upstream inputs are one-image lists filled by "1"
#' @param scales cube side lengths; scoring neighborhoods do not change
#' @param coord_names column names for spatial coordinates
#' @param simulation_times number of independent permutations per scale
#' @param min_ref_number minimum number of reference cells. if reference cell count is lower than min_ref_number, then permutation will not be performed
#' @param min_query_number minimum number of query cells. if query cell count is lower than min_query_number, then permutation will not be performed
#' @param downsample_fraction fraction of reference cells permuted per iteration
#' @param significant_pval threshold before correction across interactions
#' @param seed seed for random permutation
#' @param record_timing include computation timing in the result
#' @return A list containing `results_by_image`, `processed_images`,
#'   `scales`, `simulation_times`, and `point_type`; `timing_info` is also
#'   included when `record_timing = TRUE`
#' @importFrom methods is
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @export   

# parent function here, 9 helper functions documented below
simulate_multiscale_mvlclq_fast_3D_nonsequential <- function(
    sce_obj,
    spatial_weights_multi,
    type_column,
    point_type,
    target_types,
    observed_mvlclq_multi,
    occurrence_probs_multi,
    image_column = NULL,
    scales = c(5, 10, 50, 100, 150, 200),
    coord_names = c("x", "y", "z"),
    simulation_times = 999,
    min_ref_number = 25,
    min_query_number = 25,
    downsample_fraction = 1.0,
    significant_pval = 0.05,
    seed = 123,
    record_timing = TRUE) {

  timing_info <- list()
  if (record_timing) {
    start_time_total <- Sys.time()
    timing_info$phase_times <- list()
  }
  
  set.seed(seed)
  
  if (!is(sce_obj, "SpatialExperiment") &&
      !is(sce_obj, "SingleCellExperiment")) {
    stop("Input must be a SpatialExperiment or SingleCellExperiment object")
  }
  
  cell_data <- colData(sce_obj)

  # one-image implementation instead of silently treating several slides as one 
  # note that the deepSTARmap (mouse brain thick slice) input uses image_column = "image" and image = "1".
  if (is.null(image_column) && "image" %in% colnames(cell_data)) {
    image_values <- unique(as.character(cell_data[["image"]]))
  } else if (!is.null(image_column)) {
    image_values <- unique(as.character(cell_data[[image_column]]))
  } else {image_values <- "1"}
  
  #if (length(image_values) != 1L || is.na(image_values) ||
  #    image_values != "1") {
  #  stop("this version requires exactly one image with image ID '1', multi-image version to be updated")
  #}
  
  coords <- .mango3d_extract_spatial_coord(sce_obj)
  
  all_images_results <- list()

  # one-element loop (replaced the original multi-image. multi-image version to be updated)
  for (current_image in "1") {
    if (record_timing) image_start_time <- Sys.time()

    if (is.null(image_column)) {
      current_spatial_weights <- spatial_weights_multi
      current_observed_mvlclq <- observed_mvlclq_multi
      current_occurrence_probs <- occurrence_probs_multi
    } else {
      current_spatial_weights <- spatial_weights_multi[[current_image]]
      current_observed_mvlclq <- observed_mvlclq_multi[[current_image]]$results
      current_occurrence_probs <- occurrence_probs_multi[[current_image]]
    }
    
    # skip if no data for this image
    # commented this out, since we only have one single image now, there should be data in this image
    #if (is.null(current_spatial_weights) || is.null(current_observed_mvlclq) ||
    #    is.null(current_occurrence_probs)) {
    #  warning(paste("Missing data for image", current_image, "- skipping"))
    #  next
    #}
    
    # commented this out, since we only have one single image now, please select cell types that have nonzero counts
    # in the 2D multi-image version, this was included as some images might not contain the target / query cell type
    #if (length(current_spatial_weights$point_indices) == 0L) {
    #  warning(paste("No cells of type", point_type, "found in image", current_image, "- skipping"))
    #  next
    #}
    
    image_types <- as.character(cell_data[[type_column]])
    
    point_indices <- current_spatial_weights$point_indices
    other_indices <- current_spatial_weights$other_indices
    n_ref <- length(point_indices)
    if (anyNA(c(point_indices, other_indices)) ||
        any(c(point_indices, other_indices) < 1L) ||
        any(c(point_indices, other_indices) > length(image_types))) {
      stop("Spatial-weight indices do not match the single-image sce_obj")
    }
    
    if (n_ref < min_ref_number) {
      warning(paste("Less than", min_ref_number, "reference cells in image", current_image, "- skipping"))
      next}
    
    interaction_names <- names(current_observed_mvlclq)
    interaction_names <- interaction_names[
      interaction_names != "point_indices"
    ]
    
    interactions_to_process <- character(0)
    valid_interactions_info <- list()
    # retrieves the precomputed count for each cell type
    type_counts <- current_occurrence_probs$type_counts
    
    for (interaction_name in interaction_names) {
      current_interaction <- current_observed_mvlclq[[interaction_name]]
      current_target_types <- if (
        "target_types" %in% names(current_interaction)
      ) {
        current_interaction$target_types
      } else {
        current_interaction$target_type
      }
      
      if (length(current_target_types) == 1L) {
        target_type <- current_target_types[1L]
        
        # commented out for the same reason - remember to check cell count before testing
        #if (!(target_type %in% names(type_counts))) {
        #  warning(paste(
        #    "Skipping interaction", interaction_name, "for image",
        #    current_image, "- query type", target_type, "has 0 cells"
        #  ))
        #  next
        #}
        
        # kept these for minimum cell count
        count <- as.numeric(type_counts[[target_type]])
        if (is.na(count) || count < min_query_number) {
          warning(paste(
            "Skipping interaction", interaction_name, "for image",
            current_image, "- insufficient query cells for type:",
            target_type, "(count:", count, ", minimum required:",
            min_query_number, ")"
          ))
          next
        }
      } else {
        all_insufficient <- all(vapply(
          current_target_types,
          function(target_type) {
            if (!(target_type %in% names(type_counts))) {
              return(TRUE)
            }
            count <- as.numeric(type_counts[[target_type]])
            is.na(count) || count < min_query_number
          },
          logical(1)
        ))
        
        if (all_insufficient) {
          warning(paste(
            "Skipping combined interaction", interaction_name,
            "for image", current_image,
            "- all query types have insufficient cells (minimum required:",
            min_query_number, ")"
          ))
          next
        }
      }
      
      # adds the current interaction to the list of ones that passed the checks
      interactions_to_process <- c(
        interactions_to_process,
        interaction_name
      )
      valid_interactions_info[[interaction_name]] <- list(
        target_types = current_target_types,
        is_pairwise = length(current_target_types) == 1L
      )
    }
    
    interaction_names <- interactions_to_process
    
    # If no interactions can be processed, skip this image entirely
    if (length(interaction_names) == 0) {
      warning(paste("No interactions can be processed for image", current_image,
                    "due to insufficient query cells. Skipping image."))
      next
    }
    
    ## Step 1 Precalculation
    
    if (record_timing) phase1_start <- Sys.time()
    
    N_other <- length(other_indices)
    all_types <- unique(image_types)
    other_types <- image_types[other_indices]
    
    A_normalized <- .mango3d_build_A_normalized(current_spatial_weights,
                                                n_ref, other_indices)
    
    if (record_timing) {
      timing_info$phase_times[[as.character(current_image)]]$phase1 <-
        as.numeric(difftime(Sys.time(), phase1_start, units = "secs"))
    }
    
    ## Step 2 Extract observed scores from input
    
    if (record_timing) phase2_start <- Sys.time()
    
    observed_scores <- list()
    for (interaction_name in interaction_names) {
      target_types_vec <-
        valid_interactions_info[[interaction_name]]$target_types
      is_pairwise <-
        valid_interactions_info[[interaction_name]]$is_pairwise
      v_obs <- current_observed_mvlclq[[interaction_name]]$local_mvlclq
      
      if (length(v_obs) != n_ref) {
        stop(paste(
          "Length mismatch for interaction", interaction_name,
          ": expected", n_ref, "values but got", length(v_obs)
        ))
      }
      
      observed_scores[[interaction_name]] <- list(
        observed = v_obs,
        target_types = target_types_vec,
        is_pairwise = is_pairwise
      )
    }
    
    total_count <- current_occurrence_probs$total_count
    interaction_info <- .mango3d_prepare_interactions(
      interaction_names = interaction_names,
      observed_scores = observed_scores,
      type_counts = type_counts,
      total_count = total_count,
      all_types = all_types
    )
    
    if (record_timing) {
      timing_info$phase_times[[as.character(current_image)]]$phase2 <-
        as.numeric(difftime(Sys.time(), phase2_start, units = "secs"))
    }
    
    # Step 3: scale-specific simulation (vectorized)
    
    ref_coords <- coords[point_indices, , drop = FALSE]
    other_coords <- coords[other_indices, , drop = FALSE]
    
    results_by_scale <- list()
    
    for (scale in scales) {
      
      if (record_timing) scale_start <- Sys.time()
      
      # Step 3a Build neighborhoods
      
      if (record_timing) neighbor_build_start <- Sys.time()
      cube_info <- .mango3d_build_cube_matrix(ref_coords = ref_coords, other_coords = other_coords,
        n_ref = n_ref, scale = scale)
      
      if (record_timing) {
        neighbor_build_time <- as.numeric(difftime(Sys.time(), neighbor_build_start, units = "secs"))
      }
      
      # Step 3b: Run simulations
      
      if (record_timing) simulation_start <- Sys.time()
      
      simulated_scores <- .mango3d_run_nonsequential(
        simulation_times = simulation_times,
        other_types = other_types,
        n_ref = n_ref,
        downsample_fraction = downsample_fraction,
        neighbor_matrix = cube_info$neighbor_matrix,
        neighbor_counts = cube_info$neighbor_counts,
        N_other = N_other,
        all_types = all_types,
        A_normalized = A_normalized,
        interaction_names = interaction_names,
        interaction_info = interaction_info
      )
      
      if (record_timing) {
        simulation_time <- as.numeric(difftime(Sys.time(), simulation_start, units = "secs"))
      }
      
      scale_results <- list(
        scale = scale,
        image_id = current_image,
        interactions = .mango3d_calculate_significance(
          observed_scores = observed_scores,
          simulated_scores = simulated_scores,
          interaction_names = interaction_names,
          n_ref = n_ref,
          simulation_times = simulation_times,
          significant_pval = significant_pval,
          scale = scale
        )
      )
      
      if (record_timing) {
        stats_time <- as.numeric(difftime(Sys.time(), stats_start, units = "secs"))
        total_scale_time <- as.numeric(difftime(Sys.time(), scale_start, units = "secs"))
        
        if (is.null(
          timing_info$phase_times[[as.character(current_image)]]$scales
        )) {
          timing_info$phase_times[[as.character(current_image)]]$scales <- list()
        }
        
        timing_info$phase_times[[as.character(current_image)]]$scales[[
          as.character(scale)
        ]] <- list(
          neighbor_building = neighbor_build_time, simulation = simulation_time,
          statistics = stats_time, total = total_scale_time
        )

      }
      
      results_by_scale[[as.character(scale)]] <- scale_results
    }
    
    all_images_results[[as.character(current_image)]] <- list(
      image_id = current_image,
      results_by_scale = results_by_scale,
      scales = scales,
      point_indices = point_indices,
      interaction_names = interaction_names
    )
    
    if (record_timing) {
      image_total_time <- as.numeric(difftime(Sys.time(), image_start_time, units = "secs"))
      timing_info$phase_times[[as.character(current_image)]]$total_image_time <- image_total_time
      
    }
    
  }
  
  if (record_timing) {
    timing_info$total_runtime <- as.numeric(difftime(Sys.time(), start_time_total, units = "secs"))
    
  }
  
  result <- list(
    results_by_image = all_images_results,
    processed_images = names(all_images_results),
    scales = scales,
    simulation_times = simulation_times,
    point_type = point_type
  )
  
  if (record_timing) result$timing_info <- timing_info
  
  result
}

#' function extract_spatial_coord (only called in simulate_multiscale_mvlclq_fast_downsample)
#' extract spatial coordinates from single cell and spatial experiment objects
#' @param sce_obj SpatialExperiment or SingleCellExperiment object
#' @noRd

.mango3d_extract_spatial_coord <- function(sce_obj) {
  if (is(sce_obj, "SpatialExperiment")) {
    SpatialExperiment::spatialCoords(sce_obj)
  } else if ("spatial" %in% reducedDimNames(sce_obj)) {
    reducedDim(sce_obj, "spatial")
  } else {
    stop("Cannot find spatial coordinates in the SingleCellExperiment object")
  }
}

#' function .mango3d_build_A_normalized (only called in simulate_multiscale_mvlclq_fast_downsample)
#' build normalized weight matrix
#' @param current_spatial_weights spatial weights being taken in
#' @param n_ref number of reference cells
#' @param other_indices neighboring cell indices
#' @noRd

.mango3d_build_A_normalized <- function(current_spatial_weights,
                                        n_ref,
                                        other_indices) {
  N_other <- length(other_indices)
  A_normalized <- matrix(0, nrow = n_ref, ncol = N_other)
  
  for (i in seq_len(n_ref)) {
    neighbors <- current_spatial_weights$neighbors[[i]]
    weights <- current_spatial_weights$weights[[i]]
    
    if (length(weights) > 0L) {
      neighbor_positions <- match(neighbors, other_indices)
      sum_weights <- sum(weights)
      
      normalized_weights <- if (sum_weights > 0) {
        weights / sum_weights
      } else {
        rep(1 / length(weights), length(weights))
      }
      
      A_normalized[i, neighbor_positions] <- normalized_weights
    }
  }
  
  A_normalized
}

#' function: .mango3d_build_cube_matrix 
#' this function builds permutation neighborhood matrix, not the constructed cellular neighborhoods
#' @param ref_coords reference cell coordinates
#' @param other_coords neighboring cell coordinates
#' @param n_ref number of reference cells
#' @param scale permutation scale
#' @noRd

.mango3d_build_cube_matrix <- function(ref_coords,
                                       other_coords,
                                       n_ref,
                                       scale) {

  
  cube_list <- vector("list", n_ref)
  max_neighbors <- 0L
  half_scale <- scale / 2
  
  use_grid <- length(half_scale) == 1L &&
    is.finite(half_scale) &&
    half_scale > 0 &&
    all(is.finite(ref_coords[, 1:3, drop = FALSE])) &&
    all(is.finite(other_coords[, 1:3, drop = FALSE]))
  
  # bins are built to speed up the permutation square constructions
  # each eligible cell is assigned to a grid bin of width h:
  # floor(x / h), floor(y / h), floor(z / h)
  # for each reference cell, we retrieve the cells from the reference bin
  # and its 26 adjacent bins
  if (use_grid) {
    other_bin_x <- floor(other_coords[, 1L] / half_scale)
    other_bin_y <- floor(other_coords[, 2L] / half_scale)
    other_bin_z <- floor(other_coords[, 3L] / half_scale)
    other_keys <- paste(other_bin_x, other_bin_y, other_bin_z, sep = "\r")
    buckets <- split(seq_len(nrow(other_coords)), other_keys)
    
    for (i in seq_len(n_ref)) {
      ref_bin_x <- floor(ref_coords[i, 1L] / half_scale)
      ref_bin_y <- floor(ref_coords[i, 2L] / half_scale)
      ref_bin_z <- floor(ref_coords[i, 3L] / half_scale)
      
      x_bins <- ref_bin_x + c(-1, 0, 1)
      y_bins <- ref_bin_y + c(-1, 0, 1)
      z_bins <- ref_bin_z + c(-1, 0, 1)
      
      candidate_keys <- paste(
        rep(x_bins, each = 9L),
        rep(rep(y_bins, each = 3L), times = 3L),
        rep(z_bins, times = 9L),
        sep = "\r"
      )
      
      candidates <- sort(
        unlist(buckets[candidate_keys], use.names = FALSE)
      )
      
      if (length(candidates) > 0L) {
        keep <- abs(other_coords[candidates, 1L] - ref_coords[i, 1L]) <=
          half_scale &
          abs(other_coords[candidates, 2L] - ref_coords[i, 2L]) <=
          half_scale &
          abs(other_coords[candidates, 3L] - ref_coords[i, 3L]) <=
          half_scale
        in_cube <- candidates[keep]
      } else {
        in_cube <- integer(0)
      }
      
      cube_list[[i]] <- in_cube
      max_neighbors <- max(max_neighbors, length(in_cube))
    }
  } else {
    for (i in seq_len(n_ref)) {
      in_cube <- which(
        abs(other_coords[, 1L] - ref_coords[i, 1L]) <= half_scale &
          abs(other_coords[, 2L] - ref_coords[i, 2L]) <= half_scale &
          abs(other_coords[, 3L] - ref_coords[i, 3L]) <= half_scale
      )
      
      cube_list[[i]] <- in_cube
      max_neighbors <- max(max_neighbors, length(in_cube))
    }
  }
  
  neighbor_matrix <- matrix(
    NA_integer_,
    nrow = n_ref,
    ncol = max_neighbors
  )
  neighbor_counts <- lengths(cube_list)
  
  for (i in seq_len(n_ref)) {
    n_neighbors <- neighbor_counts[i]
    if (n_neighbors > 0L) {
      neighbor_matrix[i, seq_len(n_neighbors)] <- cube_list[[i]]
    }
  }
  
  list(
    neighbor_matrix = neighbor_matrix,
    neighbor_counts = neighbor_counts
  )
}

#' function: .mango3d_neighbor_rows
#' unpack cells in the permutation squares
#' @param neighbor_matrix one row per reference cell; entries are positions in the set of eligible cells inside the reference's permutation square. shorter rows contain NAs
#' @param neighbor_counts how many entries in each row are real
#' @param n_ref number of reference cells.
#' @noRd
.mango3d_neighbor_rows <- function(neighbor_matrix,
                                   neighbor_counts,
                                   n_ref) {
  lapply(seq_len(n_ref), function(i) {
    n_neighbors <- neighbor_counts[i]
    if (n_neighbors > 0L) {
      as.integer(neighbor_matrix[i, seq_len(n_neighbors)])
    } else {
      integer(0)
    }
  })
}

#' function: .mango3d_prepare_interactions
#' this function builds permutation neighborhood matrix, not the constructed cellular neighborhoods
#' @param ref_coords reference cell coordinates
#' @param other_coords neighboring cell coordinates
#' @param n_ref number of reference cells
#' @param scale permutation scale
#' @noRd
.mango3d_prepare_interactions <- function(interaction_names,
                                          observed_scores,
                                          type_counts,
                                          total_count,
                                          all_types) {
  requested_types <- unique(unlist(
    lapply(
      interaction_names,
      function(interaction_name) {
        as.character(observed_scores[[interaction_name]]$target_types)
      }
    ),
    use.names = FALSE
  ))
  
  missing_types <- setdiff(requested_types, all_types)
  if (length(missing_types) > 0L) {
    stop(
      "target cell type(s) absent from the image: ",
      paste(missing_types, collapse = ", ")
    )
  }
  
  # Keep the same ordering inherited from unique(image_types)
  scored_types <- all_types[all_types %in% requested_types]
  plans <- setNames(vector("list", length(interaction_names)), interaction_names)
  
  for (interaction_name in interaction_names) {
    target_types_vec <- as.character(
      observed_scores[[interaction_name]]$target_types
    )
    is_pairwise <- observed_scores[[interaction_name]]$is_pairwise
    
    denominator <- if (is_pairwise) {
      as.numeric(type_counts[target_types_vec[1L]]) / (total_count - 1)
    } else {
      prod(as.numeric(type_counts[target_types_vec])) /
        (total_count - 1)^length(target_types_vec)
    }
    
    plans[[interaction_name]] <- list(
      target_columns = match(target_types_vec, scored_types),
      denominator = denominator,
      is_pairwise = is_pairwise
    )
  }
  
  list(scored_types = scored_types, plans = plans)
}

#' function: .mango3d_store_interaction_scores
#' stores one simulations's mLCLQ scores
#' @param P_null matrix with one row per refernece cell, one column per query type each entry is the sum of normalized weights assigned to that type unter the current permutation. i.e. local weighted proportion
#' @param simulated_scores list of matrices: reference cells * simulation iterations
#' @param interaction_names cellular interaction names
#' @param interactoin_pans include target_columns, denominator, and if it is pairwise
#' @param simulation_index iteration column
#' @noRd

.mango3d_store_interaction_scores <- function(P_null,
                                              simulated_scores,
                                              interaction_names,
                                              interaction_plans,
                                              simulation_index) {
  for (interaction_name in interaction_names) {
    plan <- interaction_plans[[interaction_name]]
    target_columns <- plan$target_columns
    
    if (plan$is_pairwise) {
      v_null <- P_null[, target_columns[1L]] / plan$denominator
    } else {
      v_null <- P_null[, target_columns[1L]]
      for (j in 2:length(target_columns)) {
        v_null <- v_null * P_null[, target_columns[j]]
      }
      v_null <- v_null / plan$denominator
    }
    
    simulated_scores[[interaction_name]][, simulation_index] <- v_null
  }
  
  simulated_scores
}

#' function: .mango3d_run_nonsequential
#' @param simulation_times number of permutations
#' @param other_types cell type labels of the cells that could be neighbors in the current image
#' @param n_ref number of reference cells
#' @param downsample_fraction downsample reference cells if needed. added to reduce computing time. input as decimal, [0, 1] always kept as 1 for all analyses performed in manuscript v9.1
#' @param perm_neighbor_matrix store neighbor information for all reference cells one row per reference cell, each row contains all cells inside the reference cell's permutation square. note: here, the neighbors are not the knn neighbors, but the cells in the permutation squares
#' @param perm_neighbor_counts number of cells in each reference's permutation square
#' @param N_other number of neighboring cells in reference cell's neighborhoods. note: here, the neighbors are neighbor cells in the constructed neighborhoods. neighborhoods constructed through either knn, fixed bandwidth, or adaptive bandwidth; not the neighbors in permutation squares
#' @param C number of cell types in the current image
#' @param A_normalized normalized weights matrix (w_ij / sum(w_ik))
#' @param type_counts cell count by cell type
#' @param total_count sum of type_counts, global cell count
#' @noRd

.mango3d_run_nonsequential <- function(simulation_times,
                                       other_types,
                                       n_ref,
                                       downsample_fraction,
                                       neighbor_matrix,
                                       neighbor_counts,
                                       N_other,
                                       all_types,
                                       A_normalized,
                                       interaction_names,
                                       interaction_info) {
  N_other <- as.numeric(N_other)
  stopifnot(length(other_types) == N_other)

  n_to_shuffle <- if (downsample_fraction < 1.0) {
    max(1L, floor(n_ref * downsample_fraction))
  } else {
    n_ref
  }
  
  scored_types <- interaction_info$scored_types
  interaction_plans <- interaction_info$plans
  n_scored_types <- length(scored_types)
  scored_type_codes <- seq_len(n_scored_types)
  
  # in 2D version, interaction_info <- .mango_prepare_interactions is extracted in this function, here.
  # we moved this to the main parent function for 3D
  
  # types that are not requested by a retained interaction can be represented as NA.
  # They still occupy positions in each permutation and are sampled exactly as before; they just simply do not require a temporary score column
  original_codes <- match(other_types, scored_types)
  neighbor_rows <- .mango3d_neighbor_rows(
    neighbor_matrix,
    neighbor_counts,
    n_ref
  )
  
  update_info <- vector("list", n_ref)
  P_original <- matrix(0, nrow = n_ref, ncol = n_scored_types)
  
  for (i in seq_len(n_ref)) {
    cube_positions <- neighbor_rows[[i]]
    weighted_positions <- which(A_normalized[i, ] != 0)
    score_weights <- A_normalized[i, weighted_positions]
    original_weighted_codes <- original_codes[weighted_positions]
    
    # map nonzero scoring positions to their locations within cube_positions.
    # A zero indicates a weighted position outside the permutation cube.
    weighted_to_cube <- match(
      weighted_positions,
      cube_positions,
      nomatch = 0L
    )
    weighted_positions_in_cube <- which(weighted_to_cube > 0L)
    cube_positions_for_weighted <- weighted_to_cube[
      weighted_positions_in_cube
    ]
    
    original_score <- vapply(
      scored_type_codes,
      function(type_code) {
        sum(score_weights[which(original_weighted_codes == type_code)])
      },
      numeric(1)
    )
    
    P_original[i, ] <- original_score
    update_info[[i]] <- list(
      cube_positions = cube_positions,
      original_weighted_codes = original_weighted_codes,
      score_weights = score_weights,
      weighted_positions_in_cube = weighted_positions_in_cube,
      cube_positions_for_weighted = cube_positions_for_weighted
    )
  }
  
  simulated_scores <- setNames(
    lapply(
      interaction_names,
      function(x) matrix(0, nrow = n_ref, ncol = simulation_times)
    ),
    interaction_names
  )

  
  all_references <- seq_len(n_ref)

  for (k in seq_len(simulation_times)) {
    cells_to_shuffle <- if (downsample_fraction < 1.0) {
      sample(all_references, size = n_to_shuffle, replace = FALSE)
    } else {
      all_references
    }
    shuffle_reference <- logical(n_ref)
    shuffle_reference[cells_to_shuffle] <- TRUE
    active_references <- which(shuffle_reference & neighbor_counts > 1L)
    P_null <- P_original
    
    # ascending reference order matches the 2D implementation's order
    for (i in active_references) {
      info <- update_info[[i]]
      
      # always sample from the original labels for this cube. This is the key
      # difference from the shared, sequential shuffled_types vector.
      permuted_codes <- sample(
        original_codes[info$cube_positions],
        replace = FALSE
      )
      
      # sampling is still performed when the scoring neighborhood and cube do not intersect, preserving RNG progression and the original loop logic
      # an extreme case: cells within the knn neighborhood might acually locate outside of the shuffling cube. 
      # for example, a cell at (1.1, 0, 0) is outside a cube extending from −1 to 1 on each axis, but is closer to the center than a cell at (0.9, 0.9, 0.9) inside the cube
      # therefore, make sure to include larger permutation scales
      if (length(info$weighted_positions_in_cube) > 0L) {
        current_weighted_codes <- info$original_weighted_codes
        current_weighted_codes[info$weighted_positions_in_cube] <- permuted_codes[info$cube_positions_for_weighted]
        
        P_null[i, ] <- vapply(
          scored_type_codes,
          function(type_code) {
            sum(
              info$score_weights[
                which(current_weighted_codes == type_code)
              ]
            )
          }, numeric(1))
      }
    }
    
    simulated_scores <- .mango3d_store_interaction_scores(
      P_null = P_null,
      simulated_scores = simulated_scores,
      interaction_names = interaction_names,
      interaction_plans = interaction_plans,
      simulation_index = k
    )
  }
  
  simulated_scores
}

#' function: .mango3d_calculate_significance
#' calculate single-cell-level ses scores and p values.
#' single-cell-level p values are not used to determine significance,
#' full-tissue permutation p values are
#' @param observed_scores observed mLCLQ scores
#' @param simulated_scores simulated mlCLQ scores (mLCLQ scores calculated after permutations)
#' @param interaction_names cellular interaction names
#' @param n_ref number of reference cells
#' @param simulation_times number of permutations
#' @param significant_pval significant p value threshold before correction for multiple testing
#' @param scale permutation scale; permutation squares constructed with side lengths of scales, neighborhood sizes would not be changed
#' @noRd

.mango3d_calculate_significance <- function(observed_scores,
                                            simulated_scores,
                                            interaction_names,
                                            n_ref,
                                            simulation_times,
                                            significant_pval,
                                            scale) {
  interactions_out <- list()
  has_matrix_stats <- requireNamespace("matrixStats", quietly = TRUE)
  
  for (interaction_name in interaction_names) {
    v_obs <- observed_scores[[interaction_name]]$observed
    S_sim <- simulated_scores[[interaction_name]]
    target_types_vec <- observed_scores[[interaction_name]]$target_types
    is_pairwise <- observed_scores[[interaction_name]]$is_pairwise
    
    mu_sim <- rowMeans(S_sim)
    sigma_sim <- if (has_matrix_stats) {
      matrixStats::rowSds(S_sim)
    } else {
      apply(S_sim, 1L, sd)
    }
    null_medians <- if (has_matrix_stats) {
      matrixStats::rowMedians(S_sim)
    } else {
      apply(S_sim, 1L, median)
    }
    
    sigma_sim[sigma_sim == 0] <- 1e-10
    ses_scores <- (v_obs - mu_sim) / sigma_sim
    ses_scores[is.na(ses_scores)] <- 0.0
    sigma_sim[is.na(sigma_sim)] <- 0.0
    
    # Match the 2D cell-level tail choice around the null median.
    upper_tail_vec <- v_obs >= null_medians
    upper_counts <- rowSums(S_sim >= v_obs)
    lower_counts <- rowSums(S_sim <= v_obs)
    n_greater_equal <- ifelse(upper_tail_vec, upper_counts, lower_counts)
    p_values <- (n_greater_equal + 1) / (simulation_times + 1)
    p_values[is.na(p_values)] <- 1.0
    
    # Bonferroni correction across the interactions tested in this image.
    n_tests <- length(interaction_names)
    corrected_pval_threshold <- significant_pval / n_tests
    
    classified <- v_obs
    not_significant <- is.na(p_values) |
      p_values > corrected_pval_threshold
    classified[not_significant] <- ifelse(
      classified[not_significant] > 1,
      -1,
      -2
    )
    
    interactions_out[[interaction_name]] <- list(
      simulated_mvlclq = t(S_sim),
      observed = v_obs,
      ses_scores = ses_scores,
      p_values = p_values,
      std_deviations = sigma_sim,
      classified_mvlclq = classified,
      target_types = target_types_vec,
      is_pairwise = is_pairwise,
      categories = c(
        "Not Significant" = sum(classified < 0),
        "0.0 <= mvlclq < 1.0" =
          sum(classified >= 0 & classified < 1),
        "1.0 <= mvlclq < 1.5" =
          sum(classified >= 1.0 & classified < 1.5),
        "1.5 <= mvlclq < 2.0" =
          sum(classified >= 1.5 & classified < 2.0),
        ">= 2.0" = sum(classified >= 2.0)
      )
    )
    
  }
  
  interactions_out
}
