
#' function to compute spatial weights
#' three weight kernels to choose: 
#' gaussian (default), matern, exponential
#' three neighborhood construction methods to choose: 
#' k_nearest (default), fixed_distance, adaptive_bandwidth

#' function: .mango_kernel_weights
#' @param distances distances between neighbors and reference
#' @param max_dist maximum distance between neighbor to reference
#' @param kernel weight kernel, select one from c("gaussian", "matern", "exponential")
#' @returns A numeric vector of kernel weights, one per distance.
#' @noRd

# kernel calculation for gaussian, matern, and exponential
.mango_kernel_weights <- function(distances, max_dist, kernel) {
  if (length(distances) == 0L) return(numeric(0))
  
  if (max_dist == 0) {
    if (all(distances == 0)) return(rep(1, length(distances)))
    stop("Kernel max_dist must be positive for nonzero distances")
  }
  
  u <- distances / max_dist
  switch(kernel,
         gaussian = exp(-0.5 * u^2), # eqn 1, as of manuscript v9.1
         matern = (1 + sqrt(3) * u) * exp(-sqrt(3) * u),
         exponential = exp(-u))
}
