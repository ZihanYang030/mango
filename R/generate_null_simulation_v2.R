
#' Generate a simulated tissue with pairwise circular domains
#' @param num_cells total cells in the simulated tissue
#' @param tissue_size side length of tissue square
#' @param pattern pattern for pairwise cellular domains
#' @param n_per_pair circles per pair type  (total = 3 x n_per_pair)
#' @param circle_radius radius of each circle
#' @param min_circle_gap minimum gap between any two circles (from edge to edge)
#' @param query_types cell types included in domains
#' @param bg_type primary background cell type
#' @param base_props background cell type proportions
#' @param seed random seed
#' @returns A data frame with `x`, `y`, `cell_id`, and `cell_type` columns.
#' @importFrom stats runif
#' @export

generate_null_simulation_v2 <- function(num_cells = 5000, tissue_size = 2000,
    pattern = "MultiCircle", n_per_pair = 5, circle_radius = 100,
    min_circle_gap = 40, query_types = c("A", "B", "C"), bg_type = "D",
    base_props = c(A = 0.001, B = 0.001, C = 0.001, D = 0.997), seed = 123) {

  set.seed(seed)

  pairs <- list(AB = query_types[c(1, 2)], BC = query_types[c(2, 3)], AC = query_types[c(1, 3)])


  pos <- data.frame(x = runif(num_cells, 0, tissue_size), y = runif(num_cells, 0, tissue_size),
                    cell_id = seq_len(num_cells), cell_type = NA_character_)

  ap <- function(idx, pk) {
    pos$cell_type[idx] <<- sample(pairs[[pk]], length(idx), replace = TRUE)}

  n_total <- 3L * n_per_pair
  min_cd  <- 2 * circle_radius + min_circle_gap
  lo <- circle_radius + 2
  hi <- tissue_size - circle_radius - 2

  if (lo >= hi) stop(sprintf("circle_radius too large"))

  # one pair type per circle
  porder <- rep(names(pairs), each = n_per_pair)
  cx <- cy <- numeric(0)
  cp <- character(0)

  for (k in seq_len(n_total)) {
    placed <- FALSE
    for (trial in seq_len(100000)) {
      tx <- runif(1, lo, hi)
      ty <- runif(1, lo, hi)

      if (!length(cx) ||
          min(sqrt((cx - tx)^2 + (cy - ty)^2)) >= min_cd) {
        cx <- c(cx, tx)
        cy <- c(cy, ty)
        cp <- c(cp, porder[k])
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      warning(sprintf("placed only %d circles\n", k - 1L))
      break
    }
  }

  for (i in seq_along(cx)) {
    d <- sqrt((pos$x - cx[i])^2 + (pos$y - cy[i])^2)
    idx <- which(d <= circle_radius & is.na(pos$cell_type))

    if (length(idx)) {pos$cell_type[idx] <- sample(pairs[[cp[i]]], length(idx), replace = TRUE)}}

  # background
  bg_idx <- which(is.na(pos$cell_type))

  bg_labels <- if (length(bg_type) > 1L) bg_type else names(base_props)
  w <- base_props[bg_labels]

  #pos$cell_type[bg_idx] <- sample(bg_labels, length(bg_idx), replace = FALSE, prob = w)
  pos$cell_type[bg_idx] <- sample(bg_labels, length(bg_idx), replace = TRUE, prob = w)

  pos
}
