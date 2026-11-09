
#' compute multivariate colocation quotient
#'
#' @param sce_obj SpatialExperiment or SingleCellExperiment object
#' @param spatial_weights_multi spatial weights calculated using compute_spatial_weights_multi_image
#' @param type_column Column name for cell type
#' @param point_type Reference cell type, input as character object
#' @param target_types Query cell types, input as a vector of character objects
#' @param capped when multiple query cell types are input, user can select capped = TRUE, to limit the number of query cell types being tested in interactions. For example, when four query cell types, B, C, D, E are input, set capped = TRUE and capped_query_num = 2 will restrict it to triple-cell interaction testing (A-B-C, A-B-D, A-B-E, etc.)
#' @param capped_query_num number of query cell types being tested
#' @param occurrence_probs occurrence probability calculated with calculate_occurrence_probability_multi_image
#' @param image_column Column name for image index
#' @param n_cores number of CPUs used to perform parallel computation
#' @returns A list keyed by image ID. Each element contains `image_id`, `results` (interaction scores and reference-cell indices), and `image_stats` (cell counts).
#' @importFrom utils combn
#' @export

calculate_mvlclq_multi_image <- function(
    sce_obj, spatial_weights_multi, type_column, point_type,
    target_types, capped = TRUE, capped_query_num = 2, occurrence_probs = NULL,
    image_column = "slide", n_cores = NULL) {

  #n_cores_use <- max(1, n_cores_use)
  available_cores <- max(1L, parallel::detectCores(logical = FALSE) - 1L, na.rm = TRUE)
  n_cores_use <- if (is.null(n_cores)) available_cores else as.integer(n_cores)
  n_cores_use <- max(1L, min(n_cores_use, available_cores))

  cl <- parallel::makeCluster(n_cores_use)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)   # clean up cluster

  # when capped is set TRUE:
  # iter 1: every single query type (pariwise interaction) Ref-Q1
  # iter 2: every paired query types (triple cell interaction) Ref-Q1-Q2
  # iter x: every x query types (x+1 cell interaction) Ref-Q1-Q2-...-Qx
  # the order of the query cell types do not matter Ref-Q1-Q2 = Ref-Q2-Q1

  max_combo_size <- if (capped) as.integer(capped_query_num) else length(target_types)
  if (max_combo_size < 1) stop("`capped_query_num` must be >= 1.")
  if (max_combo_size > length(target_types))
    max_combo_size <- length(target_types)

  # pre-build the full combination list (shared across all images)
  # each element is a list: list(iter = k, targets = character vector of query types)
  combo_schedule <- list()
  for (iter in 1:max_combo_size) {
    combos <- combn(target_types, iter, simplify = FALSE)
    for (combo in combos)
      combo_schedule <- c(combo_schedule, list(list(iter = iter, targets = combo)))
  }

  # occur prob
  if (is.null(occurrence_probs)) {
    occurrence_data <- calculate_occurrence_probability_multi_image(
      sce_obj, type_column, image_column)} else { occurrence_data <- occurrence_probs }

  image_names <- names(spatial_weights_multi)

  # compute MVLCLQ values for one combination in one image
  #
  # smallest calculation loop, extracted so it can be organized in parallel
  compute_one_combo <- function(combo_entry, spatial_weights, type_counts, total_count) {

    iter <- combo_entry$iter
    current_targets <- combo_entry$targets

    n_points <- length(spatial_weights$point_indices)

    # denominator
    # unweighted global cell abundance
    # eqn 4 as of manuscript v 9.1
    if (iter == 1) {
      # singel query type
      type_count <- type_counts[current_targets]
      mvlclq_denominator <- type_count / ((total_count-1)^iter)
    } else {
      # multiple query types
      product_counts <- prod(as.numeric(type_counts[current_targets]))
      mvlclq_denominator <- product_counts / (total_count-1)^iter
    }

    # per-point loop
    mvlclq_values <- numeric(n_points)

    for (i in seq_len(n_points)) {
      neighbors_types <- spatial_weights$types[[i]]
      weights <- spatial_weights$weights[[i]]

      if (length(weights) == 0) next # if no neighbior, leaves mvlclq_values[i] = 0

      sum_weights <- sum(weights)

      # eqn 2, 3 as of manuscript v 9.1
      # mLCLQ calculation
      if (iter == 1) {
        # single query type
        type_matches <- neighbors_types == current_targets
        type_weight <- if (any(type_matches)) sum(weights[type_matches]) else 0

        if (type_weight == 0) {
          mvlclq_value <- 0
        } else if (sum_weights > 0) {
          mvlclq_value <- (type_weight/(sum_weights^iter)) / mvlclq_denominator
        } else {
          # if sum_weights <= 0, divides by number of neighbor cell types
          mvlclq_value <- if (length(neighbors_types) > 0) (type_weight/length(neighbors_types)) / mvlclq_denominator else 0
        }

      } else {
        # multi-type combination
        type_weights <- setNames(numeric(iter), current_targets)
        for (j in seq_len(iter)) {
          tm <- neighbors_types == current_targets[j]
          type_weights[j] <- if (any(tm)) sum(weights[tm]) else 0
        }

        if (prod(type_weights) == 0) {
          mvlclq_value <- 0
        } else if (sum_weights > 0) {
          mvlclq_value <- (prod(type_weights)/(sum_weights^iter)) / mvlclq_denominator
        } else {
          mvlclq_value <- if (length(neighbors_types) > 0) (prod(type_weights)/length(neighbors_types)) / mvlclq_denominator else 0
        }
      }

      mvlclq_values[i] <- mvlclq_value
    }

    list( local_mvlclq = mvlclq_values, global_mvlclq = mean(mvlclq_values), target_types = current_targets, iter = iter )
  }

  # main image loop
  all_image_results <- list()

  for (current_image in image_names) {


    image_key <- as.character(current_image)

    #spatial_weights <- spatial_weights_multi[[current_image]]
    #current_image_occur <- occurrence_data[[current_image]]
    spatial_weights <- spatial_weights_multi[[image_key]]
    current_image_occur <- occurrence_data[[image_key]]

    # missing occurrence data, skip current image
    if (is.null(current_image_occur)) {
      all_image_results[[image_key]] <- list(image_id = image_key, results = list(), note = "no occurrence data")
      next
    }

    type_counts <- current_image_occur$type_counts
    total_count <- current_image_occur$total_count

    # no reference cells in this image
    if (length(spatial_weights$point_indices) == 0) {
      all_image_results[[image_key]] <- list(
        image_id = image_key, results = list(), note = "No points of interest")
      next }

    # filter combo_schedule to types present in this image
    # here, only cell types that do not exist are filtered out
    # later, in the permutation step, there is a check for cell types with count < 25
    # Combinations whose types are all present - compute normally
    # Combinations with any missing type - zero results
    present_types <- names(type_counts)

    valid_combos <- Filter(function(ce) all(ce$targets %in% present_types), combo_schedule)
    missing_combos <- Filter(function(ce) !all(ce$targets %in% present_types), combo_schedule)

    # parallel dispatch over valid combinations
    # foreach iterates over valid_combos; each worker calls compute_one_combo
    # Variables captured in the closure (spatial_weights, type_counts,
    # total_count, compute_one_combo) are exported explicitly.
    `%dopar%` <- foreach::`%dopar%`

    computed_results <- foreach::foreach(
      ce = valid_combos,
      .packages = character(0),
      .errorhandling = "stop"
    ) %dopar% {
      compute_one_combo(ce, spatial_weights, type_counts, total_count)
    }

    result_list <- list()

    # zero-results for combinations with missing types
    for (ce in missing_combos) {
      combo_name <- if (ce$iter == 1) ce$targets else paste(sort(ce$targets), collapse = "_")
      missing <- ce$targets[!ce$targets %in% present_types]
      result_list[[combo_name]] <- list(
        local_mvlclq = rep(0, length(spatial_weights$point_indices)),
        global_mvlclq = 0,
        target_types = ce$targets,
        note = paste("missing:", paste(missing, collapse = ",")))}

    # computed results
    for (idx in seq_along(computed_results)) {
      res <- computed_results[[idx]]
      combo_name <- if (res$iter == 1) res$target_types else paste(sort(res$target_types), collapse = "_")

      if (res$iter == 1) {
        result_list[[combo_name]] <- list( local_mvlclq = res$local_mvlclq, global_mvlclq = res$global_mvlclq,
          target_type = res$target_types  # this is a scalar for single-type entries
        )
      } else {
        result_list[[combo_name]] <- list(local_mvlclq = res$local_mvlclq, global_mvlclq = res$global_mvlclq, target_types = res$target_types)
      }
    }

    result_list$point_indices <- spatial_weights$point_indices

    all_image_results[[image_key]] <- list(
      image_id = image_key, results = result_list,
      image_stats = list(total_cells = total_count, type_counts = type_counts))


  }

  return(all_image_results)
}
