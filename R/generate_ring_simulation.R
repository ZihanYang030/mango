
#' Generate a simulated tissue with a core and concentric rings
#' @param num_cells total cells in the simulated tissue
#' @param tissue_size side length of tissue square
#' @param center center of the co-centric ring structure, use tissue center if NULL
#' @param core_radius radius of the central core domain
#' @param core_types cell types filling the core
#' @param ring_radii inner and outer radii of rings, one set c(inner_radius, outer_radius) for each ring
#' @param ring_types cell types filling the ring, one entry per ring; list length should be same as ring_radii
#' @param bg_type primary background cell type
#' @param bg_prop proportion for the primary background cell type
#' @param seed random seed
#' @returns A data frame with `x`, `y`, `cell_id`, and `cell_type` columns.
#' @importFrom stats runif
#' @export

generate_ring_simulation <- function(num_cells = 3000, tissue_size = 2000,
    center = NULL, core_radius = 60, core_types = c("A"),
    ring_radii = list(c(60, 120), c(120, 180)), ring_types = list(c("B"), c("C")),
    bg_type = "D", bg_prop = 0.997, seed = 123) {

  n_ring <- length(ring_radii)

  set.seed(seed)

  if(is.null(center)) { center <- c(tissue_size/2, tissue_size/2) }

  cx0 <- center[1L]
  cy0 <- center[2L]

  #x <- runif(num_cells, 0, tissue_size)
  #y <- runif(num_cells, 0, tissue_size)
  all_positions <- data.frame(
    x = runif(num_cells, 0, tissue_size), y = runif(num_cells, 0, tissue_size),
    cell_id = seq_len(num_cells), cell_type = NA_character_, stringsAsFactors = FALSE)

  # define core and ring regions
  region_defs <- vector("list", n_ring + 1L)
  region_defs[[1L]] <- list(types = core_types,
                            inner = 0,
                            outer = core_radius)

  for(i in seq_len(n_ring)) {
    region_defs[[i+1L]] <- list(types = ring_types[[i]],
      inner = ring_radii[[i]][1L],
      outer = ring_radii[[i]][2L])}

  # assign cell types
  dist_to_center <- sqrt((all_positions$x - cx0)^2 + (all_positions$y - cy0)^2)

  # membership[cell, region] = TRUE if that cell sits in that region
  membership <- do.call(cbind, lapply(region_defs, function(rg) {
    dist_to_center >= rg$inner & dist_to_center <= rg$outer
  }))

  for (i in which(rowSums(membership) > 0)) {
    regs <- which(membership[i, ])
    pool <- unique(unlist(lapply(regs, function(k) region_defs[[k]]$types)))
    all_positions$cell_type[i] <- pool[sample.int(length(pool), 1)]}

  # bg cells
  background <- which(is.na(all_positions$cell_type))
  if (length(background)) {
    all_circle_types <- unique(unlist(c(core_types, ring_types)))
    # background query cell types divided equally
    query_prop <- (1-bg_prop) / length(all_circle_types)
    bg_props <- setNames(c(rep(query_prop, length(all_circle_types)), bg_prop), c(all_circle_types, bg_type))
    all_positions$cell_type[background] <- sample(names(bg_props), length(background), replace = TRUE, prob = bg_props)
  }

  all_positions
}
