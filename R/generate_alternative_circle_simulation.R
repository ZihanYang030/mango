#' Generate a simulated tissue with circular domains
#' @param num_cells total cells in the simulated tissue
#' @param tissue_size side length of tissue square
#' @param n_circle number of circle domains
#' @param circle_radius radius of each circular domain
#' @param circle_types cell types in each circle, 50:50 split within each circle
#' @param circle_distances distances between pairs of circles. can input as a single scalar (same distance for all pairs) or a matrix
#' @param bg_type primary background cell type
#' @param bg_prop proportion for the primary background cell type
#' @param seed random seed
#' @returns A data frame with `x`, `y`, `cell_id`, and `cell_type` columns.
#' @importFrom stats runif
#' @export

generate_alternative_circle_simulation <- function(
    num_cells = 3000, tissue_size = 2000, n_circle = 3, circle_radius = 50,
    circle_types = list(c("A", "B"), c("B", "C"), c("A", "C")),
    circle_distances = 500, bg_type = "D", bg_prop = 0.997, seed = 123) {

  set.seed(seed)

  # distance matrix
  if (is.matrix(circle_distances)) {dist_mat <- circle_distances} else {
    dist_mat <- matrix(circle_distances, nrow = n_circle, ncol = n_circle)
    diag(dist_mat) <- 0
  }

  # cell positions
  all_positions <- data.frame(
    x = runif(num_cells, 0, tissue_size),
    y = runif(num_cells, 0, tissue_size),
    cell_id = seq_len(num_cells),
    cell_type = NA_character_,
    stringsAsFactors = FALSE
  )

  # circle centers
  range <- c(circle_radius+1, tissue_size - (circle_radius+1))

  in_bounds <- function(cx, cy) {cx >= range[1] && cx <= range[2] && cy >= range[1] && cy <= range[2]}

  # given two placed centers P1, P2, find a point P3 at distance r1 from P1 and r2 from P2
  find_third_center <- function(p1x, p1y, p2x, p2y, r1, r2) {
    d <- sqrt((p2x - p1x)^2 + (p2y - p1y)^2)
    empty <- data.frame(x = numeric(0), y = numeric(0))

    if (d == 0) {
      stop("The first two circle centers coincide; P3 is not uniquely determined")
    }
    if (d > r1 + r2 || d < abs(r1 - r2)) {
      return(empty)
    }

    a <- (d^2 + r1^2 - r2^2) / (2*d)
    h <- sqrt(max(0, r1^2 - a^2))

    ux <- (p2x - p1x) / d
    uy <- (p2y - p1y) / d
    px <- -uy
    py <- ux
    cx <- p1x + a*ux
    cy <- p1y + a*uy
    candidates <- data.frame(x = c(cx + h * px, cx - h * px),
                             y = c(cy + h * py, cy - h * py))
    candidates[apply(candidates, 1, function(r) in_bounds(r[1], r[2])), ]
  }

  center_x <- numeric(n_circle)
  center_y <- numeric(n_circle)

  max_global_tries <- 10000
  global_placed <- FALSE

  for (global_attempt in seq_len(max_global_tries)) {

    # circle 1: random
    center_x[1] <- runif(1, range[1], range[2])
    center_y[1] <- runif(1, range[1], range[2])

    # circle 2: random angle, at dist_mat[1,2] from circle 1
    placed2 <- FALSE
    for (attempt2 in seq_len(1000)) {
      angle <- runif(1, 0, 2*pi)
      cx2 <- center_x[1] + dist_mat[1, 2]*cos(angle)
      cy2 <- center_y[1] + dist_mat[1, 2]*sin(angle)

      if (!in_bounds(cx2, cy2)) next

      center_x[2] <- cx2
      center_y[2] <- cy2
      placed2 <- TRUE
      break
    }
    if (!placed2) next

    # circles 3+: each must satisfy distances to all previous circles
    all_placed <- TRUE
    for (i in seq_len(n_circle)[-c(1, 2)]) {
      cands <- find_third_center(center_x[1], center_y[1], center_x[2], center_y[2], dist_mat[1, i], dist_mat[2, i])
      if (nrow(cands) == 0) { all_placed <- FALSE; break }

      # filter candidates by remaining distance constraints
      if (i > 3) {
        ok <- apply(cands, 1, function(row) {
          all(sapply(3:(i - 1), function(j) {
            abs(sqrt((row[1] - center_x[j])^2 + (row[2] - center_y[j])^2) -
                  dist_mat[j, i]) < 1e-6
          }))
        })
        cands <- cands[ok, ]
        if (nrow(cands) == 0) { all_placed <- FALSE; break }
      }

      chosen  <- cands[sample(nrow(cands), 1), ]
      center_x[i] <- chosen$x
      center_y[i]  <- chosen$y
    }

    if (all_placed) { global_placed <- TRUE; break }
  }

  if (!global_placed) stop(
    "cannot place all circles"
  )

  circle_centers <- data.frame(x = center_x, y = center_y)

  # assign cell types

  for (cell_idx in seq_len(num_cells)) {
    cx <- all_positions$x[cell_idx]
    cy <- all_positions$y[cell_idx]

    dists_to_centers <- sqrt((cx - circle_centers$x)^2 + (cy - circle_centers$y)^2)
    containing <- which(dists_to_centers <= circle_radius)

    # bg
    if (length(containing) == 0) next

    #union of cell types across all containing circles
    pool <- unique(unlist(circle_types[containing]))
    all_positions$cell_type[cell_idx] <- sample(pool, 1)
  }

  # bg
  background <- which(is.na(all_positions$cell_type))
  if (length(background) > 0) {
    all_circle_types <- unique(unlist(circle_types))
    query_prop <- (1-bg_prop) / length(all_circle_types)
    bg_props <- setNames(c(rep(query_prop, length(all_circle_types)), bg_prop), c(all_circle_types, bg_type))
    all_positions$cell_type[background] <- sample(names(bg_props), length(background), replace = TRUE, prob = bg_props)
  }

  all_positions
}
