
#' function: compute_spatial_weights_multi_image
#' @param sce_obj SpatialExperiment or SingleCellExperiment object
#' @param type_column Column name for cell type
#' @param point_type Reference cell type, input as character
#' @param image_column Column name for image index
#' @param method neighborhood construction method, 
#'        select one of c("k_nearest", "fixed_distance", "adaptive_bandwidth)
#' @param k number of k nearest neighborhoods 
#' @param max_distance maximum distance for adaptive_bandwidth neighborhood construction
#' @param bandwidth fixed distance for fixed_distance neighborhood construction
#' @param coords column names for spatial coordinates
#' @param exclude_ref_as_neighbor whether or not to consider other reference 
#'        cells as neighbors of the central reference cell for each neighborhood
#'        for example, for each A cell neighborhood, should other nearby A cells 
#'        be considered as neighbors of the central A reference cell.
#' @param kernel weight kernel, select one from c("gaussian", "matern", "exponential")
#' @returns A list keyed by image ID. Each element contains global `point_indices` and `other_indices`, plus per-reference-cell `neighbors`, `distances`, `types`, and `weights`.
#' @importFrom methods is
#' @importFrom SummarizedExperiment colData
#' @importFrom SpatialExperiment spatialCoords
#' @export

compute_spatial_weights_multi_image <- function(sce_obj, 
                                                type_column, 
                                                point_type, 
                                                image_column = "slide",
                                                method = "k_nearest", 
                                                k = 15, 
                                                max_distance = NULL, 
                                                bandwidth = NULL,
                                                coords = c("x", "y"),
                                                exclude_ref_as_neighbor = FALSE,
                                                kernel = "gaussian") {
  kernel <- match.arg(kernel, c("gaussian", "matern", "exponential"))
  
  cell_data <- colData(sce_obj)
  
  # extract coordinates and cell type information
  if (is(sce_obj, "SpatialExperiment")) {spatial_coords <- spatialCoords(sce_obj)} 
  else if (is(sce_obj, "SingleCellExperiment")) {
    if (!all(coords %in% colnames(colData(sce_obj)))) {
      stop("coord columns not found ")}
    spatial_coords <- as.matrix(colData(sce_obj)[, coords])}
  
  unique_images <- unique(cell_data[[image_column]])
  
  all_results <- list()
  
  # loop through each image
  for (current_image in unique_images) {
    
    # to avoid missing img labels
    #image_coords <- spatial_coords[cell_data[[image_column]] == current_image, ]
    #image_cell_data <- cell_data[cell_data[[image_column]] == current_image, ]
    image_mask <- cell_data[[image_column]] == current_image
    image_indices <- which(image_mask)
    image_coords <- spatial_coords[image_indices, , drop = FALSE]
    image_cell_data <- cell_data[image_indices, , drop = FALSE]
    
    point_indices_local <- which(image_cell_data[[type_column]] == point_type)
    
    if (exclude_ref_as_neighbor) {
      # index for query cells (point_type) and others within this image
      other_indices_local <- which(image_cell_data[[type_column]] != point_type)
    } else {
      # index for query cells (point_type) and all cells within this image
      other_indices_local <- seq_len(nrow(image_cell_data))
      # Include all cells as potential neighbors
    }
    
    # local indices to global indices
    # this step is needed for multi-image data to convert the cell indices 
    # within each single image back to the global index
    # use image_indices as a look-up table
    point_indices_global <- image_indices[point_indices_local]
    other_indices_global <- image_indices[other_indices_local]
    
    if (length(point_indices_local) == 0) {
      
      # empty result for this image
      all_results[[as.character(current_image)]] <- list(
        point_indices = integer(0), other_indices = other_indices_global,
        weights = list(), neighbors = list(), distances = list(),
        types = list(), image_id = current_image)
      next
    }
    
    # initialize result structure for current image
    image_result <- list()
    image_result$point_indices <- point_indices_global
    image_result$other_indices <- other_indices_global
    image_result$weights <- list()
    image_result$neighbors <- list()
    image_result$distances <- list()
    image_result$types <- list()
    image_result$image_id <- current_image
    
    for (i in seq_along(point_indices_local)) {
      current_idx_local <- point_indices_local[i]
      current_point <- image_coords[current_idx_local, ]
      
      if (exclude_ref_as_neighbor) {
        neighbor_candidates <- other_indices_local
      } else {
        neighbor_candidates <- other_indices_local[other_indices_local != current_idx_local]
      }
      
      # distances to all other points within the same image
      if (length(neighbor_candidates) == 0) {
        image_result$weights[[i]] <- numeric(0)
        image_result$neighbors[[i]] <- integer(0)
        image_result$distances[[i]] <- numeric(0)
        image_result$types[[i]] <- character(0)
        next
      }
      
      distances <- sqrt(rowSums((image_coords[neighbor_candidates, , drop = FALSE] - 
                                   matrix(current_point, 
                                          nrow = length(neighbor_candidates), 
                                          ncol = length(current_point), 
                                          byrow = TRUE))^2))
      
      if (method == "k_nearest") {
        # k nearest neighbors
        if (length(distances) <= k) { selected_indices <- seq_along(distances) } else {
          selected_indices <- order(distances)[1:k] }
        
        selected_distances <- distances[selected_indices]
        image_result$neighbors[[i]] <- image_indices[neighbor_candidates[selected_indices]]
        image_result$distances[[i]] <- selected_distances
        image_result$types[[i]] <- image_cell_data[[type_column]][neighbor_candidates[selected_indices]]
        max_dist <- max(selected_distances)
        image_result$max_distances <- max_dist
        
        weights <- .mango_kernel_weights(selected_distances, max_dist, kernel)
        image_result$weights[[i]] <- weights
        
      } else if (method == "fixed_distance") {
        #if (is.null(bandwidth)) {stop("bndwidth missing")}
        if (is.null(bandwidth) || length(bandwidth) != 1L || !is.numeric(bandwidth) || !is.finite(bandwidth) || bandwidth <= 0) {
          stop("bandwidth must be a finite and positive ")}
        
        # neighbors within fixed distance
        selected_indices <- which(distances <= bandwidth)
        
        if (length(selected_indices) == 0) {
          
          image_result$weights[[i]] <- numeric(0)
          image_result$neighbors[[i]] <- integer(0)
          image_result$distances[[i]] <- numeric(0)
          image_result$types[[i]] <- character(0)
        } else {
          selected_distances <- distances[selected_indices]
          image_result$neighbors[[i]] <- image_indices[neighbor_candidates[selected_indices]]
          image_result$distances[[i]] <- selected_distances
          
          #image_result$weights[[i]] <- .mango_kernel_weights(selected_distances, bandwidth, kernel)
          
          weights <- .mango_kernel_weights(selected_distances, bandwidth, kernel)
          image_result$weights[[i]] <- weights
          
          image_result$types[[i]] <- image_cell_data[[type_column]][neighbor_candidates[selected_indices]]
        }
        
      } else if (method == "adaptive_bandwidth") {
        if (is.null(max_distance) || is.null(k)) { stop("both max_distance and k needed") }
        
        # Select neighbors within max distance
        selected_indices <- which(distances <= max_distance)
        
        if (length(selected_indices) <= k) {
          if (length(selected_indices) == 0) {
            image_result$weights[[i]] <- numeric(0)
            image_result$neighbors[[i]] <- integer(0)
            image_result$distances[[i]] <- numeric(0)
            image_result$types[[i]] <- character(0)
            next
          }
        } else {
          # k nearest neighbors among those within max distance
          selected_indices <- selected_indices[order(distances[selected_indices])[1:k]]
        }
        image_result$neighbors[[i]] <- image_indices[neighbor_candidates[selected_indices]]
        image_result$types[[i]] <- image_cell_data[[type_column]][neighbor_candidates[selected_indices]]
        
        selected_distances <- distances[selected_indices]
        image_result$distances[[i]] <- selected_distances
        max_dist <- max(selected_distances)
        image_result$max_distances <- max_dist
        
        weights <- .mango_kernel_weights(selected_distances, max_dist, kernel)
        image_result$weights[[i]] <- weights
        
      }
    }
    

    all_results[[as.character(current_image)]] <- image_result
    
  }
  
  all_results
}

