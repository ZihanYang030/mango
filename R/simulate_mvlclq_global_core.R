
# Full-tissue permutation
# ==============================================================================
### Note: The permutation here is a full-tissue level cell-type label shuffle ###
# scale-specific permutation (simulate_multiscale_mvlclq_fast_downsample) constructs shuffling squares around each reference cell
#   -> calculate single-cell-level SES score
# full-tissue permutation (simulate_mvlclq_global_core) shuffle all cell types across the tissue section
#   -> calculate tissue-level p value
#
# how it differs from scale-specific permutation
# 1. null model: full tissue permutation; scale windows are no longer used;
#                build_neighborhood_matrix function is no longer used
# 2. scale: permutation scale loop removed (neighborhood size still comes from compute_spatial_weights_multi_image)
# 3. permute_ref: choose to either keep the reference cells fixed at their original locations or not
#                 permute_ref = FALSE -> shuffle only non-reference cells
#                 permute_ref = TRUE -> shuffle both non-reference and reference cells
#                 the original reference cell locations always stay as the fixed centers of the constructed neighborhoods either way
# 4. downsample fraction: select a random fraction of the all the shuffled cells
# ==============================================================================

#' function: simulate_mvlclq_global_core
#' @details Each simulation starts from the observed cell-type labels and independently shuffles labels among eligible cells while keeping cell locations and the constructed spatial neighborhoods fixed. No scale-specific permutation squares are built. `permute_ref` controls whether reference-type cells can participate, and `downsample_fraction` controls how many eligible cells are sampled for each shuffle.
#' @param sce_obj SpatialExperiment or SingleCellExperiment object
#' @param spatial_weights_multi spatial weights calculated using compute_spatial_weights_multi_image
#' @param type_column Column name for cell type
#' @param point_type Reference cell type, input as character object
#' @param target_types Query cell types, input as a vector of character objects
#' @param observed_mvlclq_multi observed mlCLQ scores calculated with calculate_mvlclq_multi_image
#' @param occurrence_probs_multi occurrence probability calcualted with calculate_occurrence_probability_multi_image
#' @param image_column Column name for image index. If no image_column input, then all cells are considered on the same image
#' @param coord_names column names for spatial coordinates
#' @param permute_ref choice to either permute reference cells or not
#' @param two_sided whether to use two-sided p values for the significance calculations
#' @param simulation_times number of permutations
#' @param min_ref_number minimum number of reference cells. if reference cell count is lower than min_ref_number, then permutation will not be performed
#' @param min_query_number minimum number of query cells. if query cell count is lower than min_query_number, then permutation will not be performed
#' @param downsample_fraction downsample all the cells being permuted
#' @param log_transform choice to either or not perform log transform on the observed and permuted mLCLQ scores
#' @param significant_pval significant p value threshold before correction for multiple testing this threshold is then corrected for testing multiple interactions.
#' @param seed seed for random permutation
#' @returns A list containing image-keyed full-tissue permutation results, processed image IDs, the `"global"` scale label, permutation settings, and the reference cell type.
#' @importFrom SummarizedExperiment colData
#' @importFrom stats median sd setNames
#' @export

simulate_mvlclq_global_core <- function(sce_obj,
                                        spatial_weights_multi,
                                        type_column,
                                        point_type,
                                        target_types,
                                        observed_mvlclq_multi,
                                        occurrence_probs_multi,
                                        image_column = NULL,
                                        coord_names = c("x", "y"),
                                        permute_ref = FALSE,
                                        two_sided = TRUE,
                                        simulation_times = 999,
                                        min_ref_number = 25,
                                        min_query_number = 25,
                                        downsample_fraction = 1.0,
                                        log_transform = TRUE,
                                        significant_pval = 0.05,
                                        seed = 123) {
  set.seed(seed)

  cell_data <- colData(sce_obj)

  # single image case
  # assign a single image label
  # when there is no image column input, we consider all cells as in the same img
  if (is.null(image_column)) {
    unique_images <- "single_image"
    cell_data$single_image <- "single_image"
    image_column <- "single_image"} else {
    if (!(image_column %in% colnames(cell_data))) {
      stop(paste("image column not found"))}
    unique_images <- unique(cell_data[[image_column]])
  }

  coords <- extract_spatial_coord(sce_obj)
  all_labels_global <- as.character(cell_data[[type_column]])   # gloal label vector (corrected for multiple images)
  all_images_results <- list()

  for (current_image in unique_images) {

    # image-specific inputs
    if (image_column == "single_image") {
      current_spatial_weights <- spatial_weights_multi
      current_observed_mvlclq <- observed_mvlclq_multi
      current_occurrence_probs <- occurrence_probs_multi
      image_mask <- rep(TRUE, nrow(cell_data))
    } else {
      current_spatial_weights <- spatial_weights_multi[[as.character(current_image)]]
      current_observed_mvlclq <- observed_mvlclq_multi[[as.character(current_image)]]$results
      current_occurrence_probs <- occurrence_probs_multi[[as.character(current_image)]]
      image_mask <- cell_data[[image_column]] == current_image
    }

    if (is.null(current_spatial_weights) || is.null(current_observed_mvlclq) || is.null(current_occurrence_probs)) {
      warning(paste("missing data for image")); next}

    image_indices <- which(image_mask)
    image_types_local <- as.character(cell_data[[type_column]][image_mask])

    point_indices <- current_spatial_weights$point_indices   # global from compute_spatial_weights
    other_indices <- current_spatial_weights$other_indices   # global
    n_ref <- length(point_indices)

    if (n_ref < min_ref_number) {warning(paste("not enough ref cels")); next}

    # fixed for global index
    all_labels <- all_labels_global

    #C <- length(sort(unique(as.character(cell_data[[type_column]][image_mask]))))
    all_types <- sort(unique(image_types_local))
    C <- length(all_types)

    # min query count
    interaction_names <- names(current_observed_mvlclq)
    interaction_names <- interaction_names[interaction_names != "point_indices"]

    type_counts_all <- current_occurrence_probs$type_counts
    interactions_to_process <- c()
    valid_interactions_info <- list()

    for (interaction_name in interaction_names) {
      if ("target_types" %in% names(current_observed_mvlclq[[interaction_name]])) {
        current_target_types <- current_observed_mvlclq[[interaction_name]]$target_types
      } else {current_target_types <- current_observed_mvlclq[[interaction_name]]$target_type}

      if (length(current_target_types) == 1) {
        # for single query type, skip if the query type has count lower than min query count
        target_type <- current_target_types[1]
        if (!(target_type %in% names(type_counts_all))) { warning(paste("skip", interaction_name, "img", current_image)); next}
        target_count <- as.numeric(type_counts_all[[target_type]])
        if (is.na(target_count) || target_count < min_query_number) {
          warning(paste("Skip", interaction_name, "img", current_image, "- not enough query cells for", target_type)); next
        }
      } else {
        # for multiple query types, skip if ALL the query types have count lower than min query count
        all_insufficient <- all(sapply(current_target_types, function(target_type) {
          if (!(target_type %in% names(type_counts_all))) return(TRUE)
          target_count <- as.numeric(type_counts_all[[target_type]])
          is.na(target_count) || target_count < min_query_number
        }))
        if (all_insufficient) { warning(paste("img skipped, all query types insufficient")); next }
      }

      interactions_to_process <- c(interactions_to_process, interaction_name)
      valid_interactions_info[[interaction_name]] <- list(
        target_types = current_target_types, is_pairwise = length(current_target_types) == 1)
    }

    interaction_names <- interactions_to_process

    # adjacency + observed scores
    A_normalized <- build_A_normalized(current_spatial_weights, n_ref, other_indices)

    observed_scores <- extract_observed_mvlclq(current_observed_mvlclq, valid_interactions_info,
                                               interaction_names, n_ref = n_ref)

    type_counts <- current_occurrence_probs$type_counts
    total_count <- current_occurrence_probs$total_count

    # single global-null simulation
    simulated_scores <- run_simulations_global(simulation_times = simulation_times, all_labels = all_labels,
      point_indices = point_indices, other_indices = other_indices,
      permute_ref = permute_ref, downsample_fraction = downsample_fraction,
      all_types = all_types, C = C, A_normalized = A_normalized,
      interaction_names = interaction_names,
      observed_scores = observed_scores,
      type_counts = type_counts, total_count  = total_count
    )

    interactions <- calculate_significance_fullt(
      observed_scores = observed_scores, simulated_scores = simulated_scores,
      interaction_names = interaction_names, n_ref = n_ref,
      simulation_times = simulation_times, log_transform = log_transform,
      significant_pval, scale = "global",
      current_image = current_image, two_sided = two_sided)

    interactions <- full_tissue_significance(interactions,
                                             interaction_names = interaction_names,
                                             significant_pval,
                                             simulation_times = simulation_times,
                                             two_sided = two_sided)

    # keep a scale-keyed (single "global" entry) so downstream code that
    # iterates results_by_scale still works
    all_images_results[[as.character(current_image)]] <- list(
      image_id = current_image,
      results_by_scale = list(global = list(scale = "global", image_id = current_image, interactions = interactions)),
      scales = "global", point_indices = point_indices, interaction_names = interaction_names
    )

  }

  list(
    results_by_image = all_images_results, processed_images = names(all_images_results),
    scales = "global", null_model = "tissue_global_permutation",
    permute_ref = permute_ref, two_sided = two_sided, simulation_times = simulation_times,
    point_type = point_type
  )
}

#' function: build_A_normalized
#' build normalized weight matrix
#' @param current_spatial_weights spatial weights results selected for the current image
#' @param n_ref number of reference cells
#' @param other_indices global positions of the image's neighboring cells
#' @noRd

build_A_normalized <- function(current_spatial_weights, n_ref, other_indices) {
  N_other <- length(other_indices)
  A_normalized <- matrix(0, nrow = n_ref, ncol = N_other)

  for (i in seq_len(n_ref)) {
    neighbors <- current_spatial_weights$neighbors[[i]]
    weights <- current_spatial_weights$weights[[i]]

    if (length(weights) > 0) {
      neighbor_positions <- match(neighbors, other_indices) # global vs global
      sum_w <- sum(weights)
      norm_w <- if (sum_w > 0) weights/sum_w else rep(1/length(weights), length(weights))

      keep <- !is.na(neighbor_positions) # drop unmatched safely
      A_normalized[i, neighbor_positions[keep]] <- norm_w[keep]
    }
  }
  A_normalized
}

#' function: run_simulations_global
#' @param simulation_times number of permutations
#' @param all_labels cell type label for every cell in the sce object
#' @param point_indices global positions of the image's reference cells
#' @param other_indices global positions of the image's neighboring cellss
#' @param permute_ref choice to either permute reference cells or not
#' @param downsample_fraction downsample all the cells being permuted
#' @param all_types all types
#' @param C number of cell types in the current image
#' @param A_normalized normalized weights matrix (w_ij / sum(w_ik))
#' @param interaction_names cellular interaction names
#' @param observed_scores extracted observed mLCLQ scores
#' @param type_counts cell count by cell type
#' @param total_count sum of type_counts, global cell count
#' @noRd

run_simulations_global <- function(simulation_times, all_labels, point_indices,
                                   other_indices, permute_ref, downsample_fraction,
                                   all_types, C, A_normalized, interaction_names,
                                   observed_scores, type_counts, total_count) {

  n_ref <- nrow(A_normalized)
  N_other <- ncol(A_normalized)
  stopifnot(length(other_indices) == N_other)

  # shuffle all cells (positions permuted each iteration)
  universe <- if (permute_ref) {
    sort(unique(c(point_indices, other_indices))) # every cell in the tissue
  } else {setdiff(other_indices, point_indices)} # non-reference cells only

  n_universe <- length(universe)

  n_to_shuffle <- if (downsample_fraction < 1.0) {
    max(2L, as.integer(floor(n_universe * downsample_fraction)))
  } else {n_universe}

  simulated_scores <- setNames(
    lapply(interaction_names, function(x) matrix(0, nrow = n_ref, ncol = simulation_times)),
    interaction_names
  )

  for (k in seq_len(simulation_times)) {

    labels_k <- all_labels

    # pick participating cells from the universe
    pool <- if (downsample_fraction < 1.0 && n_to_shuffle < n_universe) {
      universe[sample.int(n_universe, n_to_shuffle)]} else {universe}

    # global permutation of labels among the pool (do not change the locations)
    labels_k[pool] <- labels_k[pool][sample.int(length(pool))]

    # labels seen at the scoring-neighbor positions (order defines A_normalized cols)
    other_types_k <- labels_k[other_indices]

    # NA-safe one-hot: an out-of-range / unknown label -> all-zero row
    # This reproduces the original loop's behavior (which(all_types==NA) is empty)
    # and is what makes the non-fixed multi-image version degrade rather than error
    # turns the permuted cell-type labels into a one-hot matrix so the code can sum neighborhood weights by cell type
    type_idx <- match(other_types_k, all_types)
    valid <- !is.na(type_idx)
    L_perm <- matrix(0, nrow = N_other, ncol = C)
    colnames(L_perm) <- all_types
    if (any(valid)) {L_perm[cbind(which(valid), type_idx[valid])] <- 1 }

    P_null <- A_normalized %*% L_perm # n_ref x C, columns named by all_types

    for (interaction_name in interaction_names) {
      target_types_vec <- observed_scores[[interaction_name]]$target_types
      is_pairwise <- observed_scores[[interaction_name]]$is_pairwise

      if (is_pairwise) {
        current_target <- target_types_vec[1]
        type_count <- as.numeric(type_counts[current_target])
        mvlclq_denominator <- type_count / (total_count - 1)
        v_null <- P_null[, current_target] / mvlclq_denominator
      } else {
        product_counts <- prod(as.numeric(type_counts[target_types_vec]))
        mvlclq_denominator <- product_counts / (total_count - 1)^length(target_types_vec)
        v_null <- P_null[, target_types_vec[1]]
        for (j in 2:length(target_types_vec)) {
          v_null <- v_null * P_null[, target_types_vec[j]]
        }
        v_null <- v_null / mvlclq_denominator
      }

      simulated_scores[[interaction_name]][, k] <- v_null
    }

  }
  simulated_scores
}

#' function: calculate_sifnificance
#' @param observed_scores extracted observed mLCLQ scores
#' @param simulated_scores mLCLQ scores calculated after permutation
#' @param interaction_names cellular interaction names
#' @param n_ref number of references
#' @param simulation_times permutation times
#' @param log_transform select whether or not to perform log transform on observed and simulated mLCLQ scores
#' @param significant_pval p value significant threshold before correction for multiple testing
#' @param scale set as global for full tissue permutation
#' @param current_image current image name
#' @param two_sided whether the test is two sided or not
#' @noRd

calculate_significance_fullt <- function(observed_scores,
                                         simulated_scores,
                                         interaction_names,
                                         n_ref,
                                         simulation_times,
                                         log_transform,
                                         significant_pval,
                                         scale,
                                         current_image,
                                         two_sided = TRUE) {

  tail_mult    <- if (two_sided) 2 else 1

  interactions_out <- list()

  for (interaction_name in interaction_names) {
    if (log_transform) {
      v_obs <- log(observed_scores[[interaction_name]]$observed + 1)
      S_sim <- log(simulated_scores[[interaction_name]] + 1)
    } else {
      v_obs <- observed_scores[[interaction_name]]$observed
      S_sim <- simulated_scores[[interaction_name]]
    }

    target_types_vec <- observed_scores[[interaction_name]]$target_types
    is_pairwise <- observed_scores[[interaction_name]]$is_pairwise

    mu_sim <- rowMeans(S_sim)
    sigma_sim <- apply(S_sim, 1, sd)
    sigma_sim[sigma_sim == 0] <- 1e-10

    ses_scores <- (v_obs - mu_sim) / sigma_sim
    ses_scores[is.na(ses_scores)]   <- 0.0
    sigma_sim[is.na(sigma_sim)] <- 0.0

    null_medians <- apply(S_sim, 1, median)
    upper_tail_vec <- v_obs >= null_medians

    # tail-specific exceedance counts (one-tailed, in the selected tail)
    n_greater_equal <- numeric(n_ref)
    for (ref in seq_len(n_ref)) {
      if (upper_tail_vec[ref]) {
        n_greater_equal[ref] <- sum(S_sim[ref, ] >= v_obs[ref])
      } else {
        n_greater_equal[ref] <- sum(S_sim[ref, ] <= v_obs[ref])
      }
    }

    p_values <- numeric(n_ref)

    for (ref in seq_len(n_ref)) {

      upper_tail <- upper_tail_vec[ref]
      p_one_sided   <- (1 + n_greater_equal[ref]) / (1 + simulation_times)
      p_values[ref] <- min(1, tail_mult * p_one_sided)
    }

    p_values[is.na(p_values)] <- 1.0

    # Bonferroni across interactions
    n_tests <- length(target_types_vec) + 1
    corrected_pval_threshold <- significant_pval/n_tests

    interactions_out[[interaction_name]] <- list(
      simulated_mvlclq   = t(S_sim),
      observed = v_obs,
      ses_scores = ses_scores,
      p_values = p_values,
      std_deviations = sigma_sim,
      target_types = target_types_vec,
      is_pairwise = is_pairwise,
      two_sided = two_sided
    )

  }

  interactions_out
}

#' function: full_tissue_significance
#' calculate the empirical p value of the interaction at a full tissue permutation level
#' @param interactions list storing the singel-cell-level mLCLQ scores for each interaction
#' @param interaction_names list of cellular interaction names
#' @param significant_pval significant p value threshold before correction for multiple testing
#' @param simulation_times permutation times
#' @param two_sided Whether to classify significance using the two-sided p
#'   value (`TRUE`) or the selected-tail p value (`FALSE`).
#' @noRd

full_tissue_significance <- function(interactions, interaction_names,
                                     significant_pval,
                                     simulation_times = 999,
                                     two_sided = TRUE) {
  full_tissue_p <- numeric(length(interaction_names))
  for(interaction_name in interaction_names) {
    # rowmeans for simulated scores
    # simulated scores are organized as simulation_times * n_ref
    rowmeans <- rowMeans(interactions[[interaction_name]][["simulated_mvlclq"]])
    mean_observed <- mean(interactions[[interaction_name]][["observed"]])

    null_medians <- median(rowmeans)
    upper_tail <- mean_observed >= null_medians

    n_greater_equal <- numeric()
    if (upper_tail) {
      # eqn 7 as of manuscript v 9.1
      n_greater_equal <- sum(rowmeans >= mean_observed)
    } else {
      # eqn 8 as of manuscript v 9.1
      n_greater_equal <- sum(rowmeans <= mean_observed)
    }

    # empirical p
    interaction_p <- (1+n_greater_equal) / (1+simulation_times)

    interactions[[interaction_name]][["onetail_fulltissue_pval"]] <- interaction_p
    interactions[[interaction_name]][["tail"]] <- ifelse(upper_tail, "upper_tail", "lower_tail")
    # eqn 9 as of manuscript v 9.1
    two_tail_p <- min(2 * interaction_p, 1)
    interactions[[interaction_name]][["twotail_fulltissue_pval"]] <- two_tail_p

    p_two_thres <- significant_pval / length(interaction_names)

    interactions[[interaction_name]][["n_tests"]] <- length(interaction_names)
    interactions[[interaction_name]][["pval_significance_threshold_corrected"]] <- p_two_thres
    if (two_sided) {
      interactions[[interaction_name]][["fulltissue_pval_significance"]] <- two_tail_p <= p_two_thres
    } else {
      interactions[[interaction_name]][["fulltissue_pval_significance"]] <- interaction_p <= p_two_thres
    }

  }
  interactions
}

