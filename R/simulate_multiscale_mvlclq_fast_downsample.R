
#' function to perform multi-scale permutation

#' @param sce_obj SpatialExperiment or SingleCellExperiment object
#' @param spatial_weights_multi spatial weights calculated using compute_spatial_weights_multi_image
#' @param type_column Column name for cell type
#' @param point_type Reference cell type, input as character object
#' @param target_types Query cell types, input as a vector of character objects
#' @param observed_mvlclq_multi observed mlCLQ scores calculated with calculate_mvlclq_multi_image
#' @param occurrence_probs_multi occurrence probability calcualted with calculate_occurrence_probability_multi_image
#' @param image_column Column name for image index
#' @param scales permutation scales; permutation squares constructed with side lengths of scales, neighborhood sizes would not be changed
#' @param coord_names column names for spatial coordinates
#' @param simulation_times number of permutations
#' @param min_ref_number minimum number of reference cells. if reference cell count is lower than min_ref_number, then permutation will not be performed
#' @param min_query_number minimum number of query cells. if query cell count is lower than min_query_number, then permutation will not be performed
#' @param downsample_fraction downsample reference cells if needed. added to reduce computing time. input as decimal, a fraction between 0 and 1, inclusive; always kept as 1 for all analyses performed in manuscript v9.1
#' @param significant_pval significant p value threshold before correction for multiple testing this threshold is then corrected for testing multiple interactions.
#' @param seed seed for random permutation
#' @param record_timing record computation timing
#' @returns A list containing `results_by_image`, `processed_images`,
#'   `scales`, `simulation_times`, and `point_type`; `timing_info` is also
#'   included when `record_timing = TRUE`.
#' @importFrom methods is
#' @importFrom SummarizedExperiment colData
#' @export

# parent function, 9 helper functions documented below
simulate_multiscale_mvlclq_fast_downsample <- function(sce_obj,
                                                       spatial_weights_multi,
                                                       type_column,
                                                       point_type,
                                                       target_types,
                                                       observed_mvlclq_multi,
                                                       occurrence_probs_multi,
                                                       image_column = NULL,
                                                       scales = c(5, 10, 50, 100, 150, 200),
                                                       coord_names = c("x", "y"),
                                                       simulation_times = 999,
                                                       min_ref_number = 25,
                                                       min_query_number = 25,
                                                       downsample_fraction = 1.0,
                                                       significant_pval = 0.05,
                                                       seed = 123,
                                                       record_timing = TRUE) {

  # timing record
  timing_info <- list()
  if (record_timing) {start_time_total <- Sys.time()
  timing_info$phase_times <- list()}

  set.seed(seed)

  # cell type info
  if (is(sce_obj, "SpatialExperiment") || is(sce_obj, "SingleCellExperiment")) {
    cell_data <- colData(sce_obj)
  } else {
    stop("Input must be a SpatialExperiment or SingleCellExperiment object")
  }

  # single image case (image_column = NULL)
  if (is.null(image_column)) {
    unique_images <- "single_image"
    cell_data$single_image <- "single_image"
    image_column <- "single_image"} else {unique_images <- unique(cell_data[[image_column]])}

  # spatial coordinates
  coords <- extract_spatial_coord(sce_obj)

  # results list for all images
  all_images_results <- list()

  all_labels_global <- as.character(cell_data[[type_column]])

  # process each image
  for (current_image in unique_images) {

    if (record_timing) {
      image_start_time <- Sys.time()
    }

    # Get image-specific data
    if (is.null(image_column) || image_column == "single_image") {
      current_spatial_weights <- spatial_weights_multi
      current_observed_mvlclq <- observed_mvlclq_multi
      current_occurrence_probs <- occurrence_probs_multi
      image_mask <- rep(TRUE, nrow(cell_data))
    } else {
      image_key <- as.character(current_image)
      current_spatial_weights <- spatial_weights_multi[[image_key]]
      current_observed_mvlclq <- observed_mvlclq_multi[[image_key]]$results
      current_occurrence_probs <- occurrence_probs_multi[[image_key]]
      image_mask <- cell_data[[image_column]] == current_image
    }

    # skip if no data for this image
    if (is.null(current_spatial_weights) || is.null(current_observed_mvlclq) ||
        is.null(current_occurrence_probs)) {
      warning(paste("Missing data for image", current_image, "- skipping"))
      next
    }

    # extract image-specific data
    # fix for multi-image
    image_indices <- which(image_mask)
    image_types <- all_labels_global[image_indices]
    #image_types <- as.character(cell_data[[type_column]][image_mask])
    #image_coords <- coords[image_mask, , drop = FALSE]

    # Get reference cell information
    point_indices <- current_spatial_weights$point_indices
    other_indices <- current_spatial_weights$other_indices
    n_ref <- length(point_indices)

    # Check minimum reference cells
    if (n_ref < min_ref_number) {
      warning(paste("Less than", min_ref_number, "reference cells in image",
                    current_image, "- skipping"))
      next
    }


    # extract interaction names (there is a "point_indices" in observed_mvlclq_multi[[image_key]]$results,
    # recording all refernce cell indices, remove it from the interaction list)
    interaction_names <- names(current_observed_mvlclq)
    interaction_names <- interaction_names[
      interaction_names != "point_indices"
    ]

    # filter interactions based on minimum query cell requirements
    interactions_to_process <- c()
    valid_interactions_info <- list()

    for (interaction_name in interaction_names) {
      # extract target types for this interaction
      if("target_types" %in% names(current_observed_mvlclq[[interaction_name]])) {
        current_target_types <- current_observed_mvlclq[[interaction_name]]$target_types
      } else {
        current_target_types <- current_observed_mvlclq[[interaction_name]]$target_type
      }

      # Get type counts from occurrence probabilities
      type_counts <- current_occurrence_probs$type_counts

      # For pairwise interactions: skip if query type < 25
      if (length(current_target_types) == 1) {
        target_type <- current_target_types[1]
        # Check if type exists in type_counts (handles 0 cell case)
        # this is added for processing multiple images. when processing across
        # multiple images, some images may not contain the query cell types at all
        if (!(target_type %in% names(type_counts))) {
          warning(paste("Skipping interaction", interaction_name, "for image", current_image,
                        "- query type", target_type, "has 0 cells"))
          next
        }
        count <- as.numeric(type_counts[[target_type]])
        if (is.na(count) || count < min_query_number) {
          warning(paste("Skipping interaction", interaction_name, "for image", current_image,
                        "- insufficient query cells for type:", target_type,
                        "(count:", count, ", minimum required:", min_query_number, ")"))
          next
        }
      }
      # For combined interactions: skip only if ALL query types < 25
      else {
        # Check if ALL types have < 25 cells
        all_insufficient <- all(sapply(current_target_types, function(target_type) {
          # If type doesn't exist in type_counts, treat as 0 cells
          if (!(target_type %in% names(type_counts))) {
            return(TRUE)
          }
          count <- as.numeric(type_counts[[target_type]])
          return(is.na(count) || count < min_query_number)
        }))

        if (all_insufficient) {
          warning(paste("Skipping combined interaction", interaction_name, "for image", current_image,
                        "- all query types have insufficient cells (minimum required:", min_query_number, ")"))
          next
        }
      }

      interactions_to_process <- c(interactions_to_process, interaction_name)
      valid_interactions_info[[interaction_name]] <- list(
        target_types = current_target_types,
        is_pairwise = length(current_target_types) == 1
      )
    }

    # update interaction_names to only include interactions with sufficient cells
    interaction_names <- interactions_to_process

    # If no interactions can be processed, skip this image entirely
    if (length(interaction_names) == 0) {
      warning(paste("No interactions can be processed for image", current_image,
                    "due to insufficient query cells. Skipping image."))
      next
    }


    ## Step 1 Precalculation

    if (record_timing) phase1_start <- Sys.time()

    # function extracted to an individual one
    precalculation_list <- simulation_precalculation(current_spatial_weights,
                                                     n_ref,
                                                     image_types,
                                                     other_indices,
                                                     all_labels_global)

    A_normalized <- precalculation_list$A_normalized
    other_types <- precalculation_list$other_types
    N_other <- as.numeric(precalculation_list$N_other)
    C <- as.numeric(precalculation_list$C)
    all_types <- precalculation_list$all_types

    if (record_timing) {
      timing_info$phase_times[[as.character(current_image)]]$phase1 <-
        as.numeric(difftime(Sys.time(), phase1_start, units = "secs"))
    }

    ## Step 2 Extract observed scores from input

    if (record_timing) phase2_start <- Sys.time()

    # function extracted to an individual one
    observed_scores <- extract_observed_mvlclq(current_observed_mvlclq,
                                               valid_interactions_info,
                                               interaction_names,
                                               n_ref = n_ref
                                               )

    if (record_timing) {
      timing_info$phase_times[[as.character(current_image)]]$phase2 <-
        as.numeric(difftime(Sys.time(), phase2_start, units = "secs"))
    }

    # Step 3: scale-specific simulation (vectorized)

    # Get type counts and total count for MVLCLQ denominators
    type_counts <- current_occurrence_probs$type_counts
    total_count <- current_occurrence_probs$total_count

    ref_coords <- coords[point_indices, , drop = FALSE]
    other_coords <- coords[other_indices, , drop = FALSE]
    # ref_coords <- image_coords[point_indices, , drop = FALSE]
    # other_coords <- image_coords[other_indices, , drop = FALSE]

    results_by_scale <- list()

    for (scale in scales) {

      if (record_timing) scale_start <- Sys.time()

      # Step 3a Build neighborhoods

      if (record_timing) neighbor_build_start <- Sys.time()

      nbr <- build_neighborhood_matrix(ref_coords, other_coords, n_ref, scale)
      perm_neighbor_matrix <- nbr$perm_neighbor_matrix
      perm_neighbor_counts <- nbr$perm_neighbor_counts

      if (record_timing) {
        neighbor_build_time <- as.numeric(difftime(Sys.time(), neighbor_build_start, units = "secs"))
      }

      # Step 3b: Run simulations

      if (record_timing) simulation_start <- Sys.time()

      simulated_scores <- run_simulations_nonsequential(
        simulation_times = simulation_times,
        other_types = other_types,
        n_ref = n_ref,
        downsample_fraction = downsample_fraction,
        perm_neighbor_matrix = perm_neighbor_matrix,
        perm_neighbor_counts = perm_neighbor_counts,
        N_other = N_other,
        all_types = all_types,
        C = C,
        A_normalized = A_normalized,
        interaction_names = interaction_names,
        observed_scores = observed_scores,
        type_counts = type_counts,
        total_count = total_count
      )

      if (record_timing) {
        simulation_time <- as.numeric(difftime(Sys.time(), simulation_start, units = "secs"))
      }

      # Step 4: statistical significance (VECTORIZED)

      if (record_timing) stats_start <- Sys.time()

      scale_results <- list(
        scale = scale,
        image_id = current_image,
        interactions = calculate_significance(
          observed_scores = observed_scores,
          simulated_scores  = simulated_scores,
          interaction_names = interaction_names,
          n_ref = n_ref,
          simulation_times = simulation_times,
          significant_pval = significant_pval,
          scale = scale,
          current_image = current_image)
      )

      if (record_timing) {
        stats_time <- as.numeric(difftime(Sys.time(), stats_start, units = "secs"))
        total_scale_time <- as.numeric(difftime(Sys.time(), scale_start, units = "secs"))

        scale_timing <- list(
          neighbor_building = neighbor_build_time,
          simulation = simulation_time,
          statistics = stats_time,
          total = total_scale_time
        )

        if (is.null(timing_info$phase_times[[as.character(current_image)]]$scales)) {
          timing_info$phase_times[[as.character(current_image)]]$scales <- list()
        }
        timing_info$phase_times[[as.character(current_image)]]$scales[[as.character(scale)]] <- scale_timing

      }

      results_by_scale[[as.character(scale)]] <- scale_results
    }

    # Store results for this image
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

  # Calculate total runtime
  if (record_timing) {
    timing_info$total_runtime <- as.numeric(difftime(Sys.time(), start_time_total, units = "secs"))

  }

  # Return multi-image, multi-scale, multi-interaction results
  result <- list(
    results_by_image = all_images_results,
    processed_images = names(all_images_results),
    scales = scales,
    simulation_times = simulation_times,
    point_type = point_type
  )

  # Add timing information if requested
  if (record_timing) { result$timing_info <- timing_info }

  return(result)
}

#' function extract_spatial_coord (only called in simulate_multiscale_mvlclq_fast_downsample)
#' extract spatial coordinates from single cell and spatial exepriment objects
#' @param sce_obj SpatialExperiment or SingleCellExperiment object
#' @param coords spatial coordinate
#' @noRd

extract_spatial_coord <- function(sce_obj) {
  if (methods::is(sce_obj, "SpatialExperiment")) {
    coords <- SpatialExperiment::spatialCoords(sce_obj)
  } else if (methods::is(sce_obj, "SingleCellExperiment")) {
    if (!("spatial" %in% SingleCellExperiment::reducedDimNames(sce_obj))) {
      stop("Cannot find spatial coordinates in the SingleCellExperiment object")
    }
    coords <- SingleCellExperiment::reducedDim(sce_obj, "spatial")
  } else {
    stop("Input must be a SpatialExperiment or SingleCellExperiment object")
  }

  coords
}

#' function simulation_precalculation
#' prepare label metadata and normalized weight matrix
#' @param current_spatial_weights spatial weight object
#' @param n_ref number of reference cells
#' @param image_types cell type label of every cell in the current image
#' @param other_indices cells that could be neighbors in the current image
#' @noRd

simulation_precalculation <- function(current_spatial_weights,
                                      n_ref,
                                      image_types,
                                      other_indices,
                                      all_labels_global) {

  # The unused initial one-hot matrix L is removed. -July 2026
  N_other <- length(other_indices)
  all_types <- unique(image_types)
  C <- length(all_types)
  other_types <- all_labels_global[other_indices]
  #other_types <- image_types[other_indices]

  # A_normalized = w_ij / sum(w_ik)
  # A_noramlized match to eqn 3 (the weights part, without the cell type proportions)
  A_normalized <- matrix(0, nrow = n_ref, ncol = N_other)

  for (i in seq_len(n_ref)) {
    neighbors <- current_spatial_weights$neighbors[[i]]
    weights <- current_spatial_weights$weights[[i]]

    if (length(weights) > 0) {
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



  return(list(C = C,
              A_normalized = A_normalized,
              other_types = other_types,
              N_other = N_other,
              all_types = all_types))
}

#' function extract_observed_mvlclq
#' extract observed scores directly from observed_mvlclq_multi
#' @param current_observed_mvlclq source of the already calcualted scores for the current image
#' @param valid_interactions_info target types & whether if it is pairwise
#' @param interaction_names cellular interaction names
#' @param n_ref number of reference cells
#' @noRd

extract_observed_mvlclq <- function(current_observed_mvlclq,
                                    valid_interactions_info,
                                    interaction_names,
                                    n_ref = NULL) {
  observed_scores <- list()

  for (interaction_name in interaction_names) {
    target_types_vec <- valid_interactions_info[[interaction_name]]$target_types
    is_pairwise <- valid_interactions_info[[interaction_name]]$is_pairwise

    # Extract local_mvlclq directly from input
    v_obs <- current_observed_mvlclq[[interaction_name]]$local_mvlclq

    # Verify length matches number of reference cells
    if (length(v_obs) != n_ref) {
      stop(paste("Length mismatch for interaction", interaction_name,
                 ": expected", n_ref, "values but got", length(v_obs)))
    }

    observed_scores[[interaction_name]] <- list(
      observed = v_obs,
      target_types = target_types_vec,
      is_pairwise = is_pairwise
    )
  }

  return(observed_scores)

}

#' function: build_neighborhood_matrix
#' this function builds permutation neighborhood matrix, not the constructed cellular neighborhoods
#' @param ref_coords reference cell coordinates
#' @param other_coords neighboring cell coordinates
#' @param n_ref number of reference cells
#' @param scale permutation scale
#' @noRd

build_neighborhood_matrix <- function(ref_coords,
                                      other_coords,
                                      n_ref,
                                      scale) {

  perm_neighbor_list <- vector("list", n_ref)
  max_neighbors <- 0
  half_scale <- scale / 2

  use_grid <- length(half_scale) == 1L &&
    is.finite(half_scale) &&
    half_scale > 0 &&
    ncol(ref_coords) >= 2L &&
    ncol(other_coords) >= 2L &&
    all(is.finite(ref_coords[, 1:2, drop = FALSE])) &&
    all(is.finite(other_coords[, 1:2, drop = FALSE]))

  # bins are built to speed up the permutation square constructions
  # each eligible cell is assigned to a grid bin of width h:
  # floor(x / h) and floor(y / h)
  # for each reference cell, we retrieve the cells from the reference bin
  # and its 8 adjacent bins
  if (use_grid) {
    other_bin_x <- floor(other_coords[, 1L] / half_scale)
    other_bin_y <- floor(other_coords[, 2L] / half_scale)
    other_keys <- paste(other_bin_x, other_bin_y, sep = "\r")
    buckets <- split(seq_len(nrow(other_coords)), other_keys)

    for (i in seq_len(n_ref)) {
      ref_bin_x <- floor(ref_coords[i, 1L] / half_scale)
      ref_bin_y <- floor(ref_coords[i, 2L] / half_scale)

      x_bins <- ref_bin_x + c(-1, 0, 1)
      y_bins <- ref_bin_y + c(-1, 0, 1)
      candidate_keys <- paste(
        rep(x_bins, each = 3L),
        rep(y_bins, times = 3L),
        sep = "\r"
      )

      candidates <- sort(
        unlist(buckets[candidate_keys], use.names = FALSE)
      )

      if (length(candidates) > 0L) {
        keep <- abs(other_coords[candidates, 1L] - ref_coords[i, 1L]) <=
          half_scale &
          abs(other_coords[candidates, 2L] - ref_coords[i, 2L]) <=
          half_scale
        in_square <- candidates[keep]
      } else {
        in_square <- integer(0)
      }

      perm_neighbor_list[[i]] <- in_square
      max_neighbors <- max(max_neighbors, length(in_square))
    }
  } else {
    for (i in seq_len(n_ref)) {
      in_square <- which(
        abs(other_coords[, 1L] - ref_coords[i, 1L]) <= half_scale &
          abs(other_coords[, 2L] - ref_coords[i, 2L]) <= half_scale
      )
      perm_neighbor_list[[i]] <- in_square
      max_neighbors <- max(max_neighbors, length(in_square))
    }
  }

  perm_neighbor_matrix <- matrix(
    NA_integer_,
    nrow = n_ref,
    ncol = max_neighbors
  )
  perm_neighbor_counts <- lengths(perm_neighbor_list)

  for (i in seq_len(n_ref)) {
    n_neighbors <- perm_neighbor_counts[i]
    if (n_neighbors > 0L) {
      perm_neighbor_matrix[i, seq_len(n_neighbors)] <- perm_neighbor_list[[i]]
    }
  }


  list(
    perm_neighbor_matrix = perm_neighbor_matrix,
    perm_neighbor_counts = perm_neighbor_counts
  )
}

#' function: .mango_prepare_interactions
#' this function builds permutation neighborhood matrix, not the constructed cellular neighborhoods
#' @param ref_coords reference cell coordinates
#' @param other_coords neighboring cell coordinates
#' @param n_ref number of reference cells
#' @param scale permutation scale
#' @noRd

.mango_prepare_interactions <- function(interaction_names,
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
      "Target cell type(s) absent from the image: ",
      paste(missing_types, collapse = ", ")
    )
  }

  # keep the original all_types ordering for deterministic column lookup
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

#' function: .mango_store_interaction_scores
#' stores one simulations's mLCLQ scores
#' @param P_null matrix with one row per refernece cell, one column per query type each entry is the sum of normalized weights assigned to that type unter the current permutation. i.e. local weighted proportion
#' @param simulated_scores list of matrices: reference cells * simulation iterations
#' @param interaction_names cellular interaction names
#' @param interactoin_pans include target_columns, denominator, and if it is pairwise
#' @param simulation_index iteration column
#' @noRd

.mango_store_interaction_scores <- function(P_null,
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

#' function: .mango_perm_neighbor_rows
#' unpack cells in the permutation squares
#' @param perm_neighbor_matrix one row per reference cell; entries are positions in the set of eligible cells inside the reference's permutation square. shorter rows contain NAs
#' @param perm_neighbor_counts how many entries in each row are real
#' @param n_ref number of reference cells.
#' @noRd

.mango_perm_neighbor_rows <- function(perm_neighbor_matrix, perm_neighbor_counts, n_ref) {
  lapply(seq_len(n_ref), function(i) {
    n_neighbors <- perm_neighbor_counts[i]
    if (n_neighbors > 0L) {
      as.integer(perm_neighbor_matrix[i, seq_len(n_neighbors)])
    } else {
      integer(0)
    }
  })
}

#' function: run_simulations_nonsequential
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

run_simulations_nonsequential <- function(simulation_times,
                                          other_types,
                                          n_ref,
                                          downsample_fraction,
                                          perm_neighbor_matrix,
                                          perm_neighbor_counts,
                                          N_other,
                                          all_types,
                                          C,
                                          A_normalized,
                                          interaction_names,
                                          observed_scores,
                                          type_counts,
                                          total_count) {
  N_other <- as.numeric(N_other)
  stopifnot(length(other_types) == N_other, length(all_types) == C)

  n_to_shuffle <- if (downsample_fraction < 1.0) {
    max(1, floor(n_ref * downsample_fraction))
  } else {
    n_ref
  }

  interaction_info <- .mango_prepare_interactions(
    interaction_names = interaction_names,
    observed_scores = observed_scores,
    type_counts = type_counts,
    total_count = total_count,
    all_types = all_types
  )
  scored_types <- interaction_info$scored_types
  interaction_plans <- interaction_info$plans
  n_scored_types <- length(scored_types)

  original_codes <- match(other_types, scored_types)
  perm_neighbor_rows <- .mango_perm_neighbor_rows(
    perm_neighbor_matrix,
    perm_neighbor_counts,
    n_ref
  )
  scored_type_codes <- seq_len(n_scored_types)

  # precompute everything that is constant for a reference cell at this scale
  update_info <- lapply(seq_len(n_ref), function(i) {
    perm_neighbors_idx <- perm_neighbor_rows[[i]]
    weighted_positions <- which(A_normalized[i, ] != 0)
    score_weights <- A_normalized[i, weighted_positions]
    original_weighted_codes <- original_codes[weighted_positions]

    # map each nonzero scoring position to its location within this reference
    # cell's permutation square. A zero means that position is not shuffled
    weighted_to_square <- match(
      weighted_positions,
      perm_neighbors_idx,
      nomatch = 0L
    )
    weighted_positions_in_square <- which(weighted_to_square > 0L)
    square_positions_for_weighted <- weighted_to_square[
      weighted_positions_in_square
    ]

    original_score <- vapply(
      scored_type_codes,
      function(type_code) {
        sum(score_weights[which(original_weighted_codes == type_code)])
      },
      numeric(1)
    )

    list(
      neighbors = perm_neighbors_idx,
      original_weighted_codes = original_weighted_codes,
      score_weights = score_weights,
      weighted_positions_in_square = weighted_positions_in_square,
      square_positions_for_weighted = square_positions_for_weighted,
      original_score = original_score
    )
  })

  P_original <- do.call(
    rbind,
    lapply(update_info, function(info) info$original_score)
  )

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
      sample(
        all_references,
        size = n_to_shuffle,
        replace = FALSE
      )
    } else {
      all_references
    }

    # Logical indexed assignment is faster than seq_len(n_ref) %in% cells.
    shuffle_reference <- logical(n_ref)
    shuffle_reference[cells_to_shuffle] <- TRUE

    # which() restores ascending reference order, matching the original lapply
    # order and therefore preserving the assignment of RNG draws to cells.
    active_references <- which(
      shuffle_reference & perm_neighbor_counts > 1L
    )

    P_null <- P_original

    for (i in active_references) {
      info <- update_info[[i]]

      # The original implementation calls sample() for every shuffled square
      # with >1 cells, even if none of those cells has a nonzero scoring weight.
      # Keep that call to preserve both logic and RNG progression.
      permuted_codes <- sample(
        original_codes[info$neighbors],
        replace = FALSE
      )

      if (length(info$weighted_positions_in_square) > 0L) {
        current_weighted_codes <- info$original_weighted_codes
        current_weighted_codes[
          info$weighted_positions_in_square
        ] <- permuted_codes[info$square_positions_for_weighted]

        # This is the original
        # sum(A_normalized[i, which(local_types == current_type)])
        # calculation restricted to the nonzero entries of A_normalized row i.
        P_null[i, ] <- vapply(
          scored_type_codes,
          function(type_code) {
            sum(
              info$score_weights[
                which(current_weighted_codes == type_code)
              ]
            )
          },
          numeric(1)
        )
      }
    }

    simulated_scores <- .mango_store_interaction_scores(
      P_null = P_null,
      simulated_scores = simulated_scores,
      interaction_names = interaction_names,
      interaction_plans = interaction_plans,
      simulation_index = k
    )

  }

  simulated_scores
}

#' function: calculate_significance
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
#' @param current_image current image being processed
#' @noRd

calculate_significance <- function(observed_scores,
                                   simulated_scores,
                                   interaction_names,
                                   n_ref,
                                   simulation_times,
                                   significant_pval,
                                   scale,
                                   current_image) {
  interactions_out <- list()

  #has_matrix_stats <- requireNamespace("matrixStats", quietly = TRUE)

  for (interaction_name in interaction_names) {
    v_obs <- observed_scores[[interaction_name]]$observed
    S_sim <- simulated_scores[[interaction_name]]
    target_types_vec <- observed_scores[[interaction_name]]$target_types
    is_pairwise <- observed_scores[[interaction_name]]$is_pairwise

    mu_sim <- rowMeans(S_sim)

    #if (has_matrix_stats) {
    #  sigma_sim <- matrixStats::rowSds(S_sim)
    #  null_medians <- matrixStats::rowMedians(S_sim)
    #} else {
    #  sigma_sim <- apply(S_sim, 1L, sd)
    #  null_medians <- apply(S_sim, 1L, median)
    #}

    sigma_sim <- matrixStats::rowSds(S_sim)
    null_medians <- matrixStats::rowMedians(S_sim)

    sigma_sim[sigma_sim == 0] <- 1e-10

    ses_scores <- (v_obs-mu_sim) / sigma_sim
    ses_scores[is.na(ses_scores)] <- 0.0
    sigma_sim[is.na(sigma_sim)] <- 0.0

    upper_tail_vec <- v_obs >= null_medians

    # Matrix-vector comparison recycles v_obs once per column, producing the
    # same per-row counts as the original cell loop.
    upper_counts <- rowSums(S_sim >= v_obs)
    lower_counts <- rowSums(S_sim <= v_obs)
    n_greater_equal <- ifelse(upper_tail_vec, upper_counts, lower_counts)

    # these are single cell level p values, currently not used in manuscript v 9.1
    # all p values used in manuscript v9.1 come from full-tissue permutation, permuted in non-overlapping squares
    # see function simulate_mvlclq_global_core for full-tissue permutation p values
    p_values <- (1 + n_greater_equal) / (1 + simulation_times)

    p_values[is.na(p_values)] <- 1.0

    # multiple test correction
    n_tests <- length(interaction_names)
    corrected_pval_threshold <- significant_pval / n_tests

    # we can also classify the observed mLCLQ scores.
    # enrichment correspond to mLCLQ scores above 1,
    # depletion correspond to mLCLQ scores below 1
    # this categorization is used in CAF cell classification. Aug, 2025.
    # insignificant mLCLQ scores are assigned negative values to differentiate from significant ones,
    # this negative score assignment is only used in classification step, would not affect SES scores.

    # overall, raw observed mLCLQ score classification is not used anymore. March, 2026.
    # CAF cell classification is now based on single-cell-level SES Figure 4g,
    # not the raw mLCLQ scores any more. Oct, 2026.
    classified <- v_obs
    not_significant <- is.na(p_values) |
      p_values > corrected_pval_threshold
    classified[not_significant] <- ifelse(classified[not_significant] > 1, -1, -2)

    interactions_out[[interaction_name]] <- list(
      simulated_mvlclq = t(S_sim),
      observed = v_obs,
      ses_scores = ses_scores,
      p_values = p_values,
      std_deviations = sigma_sim,
      classified_mvlclq = classified,
      target_types = target_types_vec,
      is_pairwise = is_pairwise,
      categories = c("Not Significant" = sum(classified < 0),
                     "0.0 <= mvlclq < 1.0" = sum(classified >= 0 & classified < 1),
                     "1.0 <= mvlclq < 1.5" = sum(classified >= 1.0 & classified < 1.5),
                     "1.5 <= mvlclq < 2.0" = sum(classified >= 1.5 & classified < 2.0),
                     ">= 2.0" = sum(classified >= 2.0))
    )
  }


  interactions_out
}
