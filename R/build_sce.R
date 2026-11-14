#' helper function for converting the generated simulation data to a sce object
#'
#' @param sim_dat simulation data generated using `generate_null_simulation_v2`, `generate_alternative_circle_simulation`, or `generate_ring_simulation` functions
#' @returns sce object for mango analysis
#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export

build_sce <- function(sim_dat) {
  counts <- matrix(0L, nrow = 0L, ncol = nrow(sim_dat), dimnames = list(NULL, paste0("Cell_", sim_dat$cell_id)))

  col_data <- DataFrame(cell_id = sim_dat$cell_id, cell_type = sim_dat$cell_type,
    x_coord = sim_dat$x, y_coord = sim_dat$y, image_number = 1L)

  #spatial_coords <- sim_dat[, c("x", "y")]
  spatial_coords <- as.matrix(sim_dat[, c("x", "y")])
  rownames(spatial_coords) <- paste0("Cell_", sim_dat$cell_id)

  sce <- SingleCellExperiment(
    assays = list(counts = counts), colData = col_data,
    reducedDims = list(spatial = spatial_coords))

  sce
}
