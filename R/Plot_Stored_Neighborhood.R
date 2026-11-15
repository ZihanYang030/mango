#' function to plot the stored neighborhood and weights
#'
#' @param sce_obj sce_obj SpatialExperiment or SingleCellExperiment object description
#' @param weights_data spatial weights calculated using compute_spatial_weights_multi_image
#' @param type_column Column name for cell type
#' @param reference_cell the reference cell being plotted
#' @param image_id image id for plotting
#' @param reference_by either use sce_index or cell position as the reference
#' @param coords sptial coordinates
#' @param distance_coords coordinate columns used to measure neighbor distances
#' @param neighbor_color_by either color by cell type or weight: c("cell_type", "weight")
#' @param weight_name legend title for weight
#' @param weight_limits limit of the weight
#' @param zoom crop the panel to the reference cell and its neighbors
#' @param show_links draw a line from the reference cell to each stored neighbor
#' @return A ggplot object.
#' @export

Plot_Stored_Neighborhood <- function(sce_obj, weights_data, type_column, reference_cell,
                                     image_id = NULL, reference_by = c("sce_index", "position"), coords = c("x", "y"),
                                     distance_coords = NULL, neighbor_color_by = c("cell_type", "weight"),
                                     weight_name = NULL, weight_limits = NULL, zoom = FALSE,
                                     show_links = FALSE, title = "Neighborhood Plot") {

  reference_by <- match.arg(reference_by)
  neighbor_color_by <- match.arg(neighbor_color_by)

  background_color <- "#E8E8E8"; background_size <- 3; background_alpha <- 0.8
  neighbor_size <- 4
  reference_color <- "#D73027"
  reference_size <- 4
  link_color <- "#4D4D4D"; link_alpha <- 0.90; link_width <- 1
  cell_type_palette <- c(
    "#F9C74F", "#67A8CD", "#FF9D9F", "#D1EDC2", "#DEB6F0",
    "#A2D2E7", "#CF9F88", "#6FB3A8", "#B3E19B", "#EBD7A5",
    "#F36569", "#FF8831", "#FFEC99", "#A6CA83", "#7D67C7",
    "#FFC17F", "#466399", "#ED729B", "#947876", "#FBABFF",
    "#AA6DEA", "#B0A079", "#92F1E8", "#DFF7A5", "#B7CBB5",
    "#E07A5F", "#BD9E39", "#006073", "#099396", "#33A02C",
    "#D28F38", "#52F02B", "#FAF219", "#DC19FA"
  )

  cell_data <- as.data.frame(SummarizedExperiment::colData(sce_obj))
  #n_cells <- nrow(cell_data)
  cell_types <- as.character(cell_data[[type_column]])

  get_coordinates <- function(column_names) {
    source_data <- if (all(column_names %in% colnames(cell_data))) {
      cell_data
    } else {as.data.frame(SpatialExperiment::spatialCoords(sce_obj))}
    result <- as.matrix(source_data[, column_names, drop = FALSE])
    storage.mode(result) <- "double"
    result
  }

  plot_coordinates <- get_coordinates(coords)

  if (is.null(distance_coords)) {distance_coords <- if (all(c("x", "y", "z") %in% colnames(cell_data))) c("x", "y", "z") else coords}
  distance_coordinates <- get_coordinates(distance_coords)

  # this is for multi-image case: select the image-specific weight result
  if ("point_indices" %in% names(weights_data)) {
    spatial_weights <- weights_data
  } else if (is.null(image_id)) { spatial_weights <- weights_data[[1L]] } else {
    image_key <- as.character(image_id)
    spatial_weights <- if (image_key %in% names(weights_data)) {
      weights_data[[image_key]]
    } else {weights_data[[which(vapply(
        weights_data,
        function(x) identical(as.character(x$image_id)[1L], image_key), logical(1)
      ))[1L]]]
    }
  }

  # find the reference cell
  point_indices <- as.integer(spatial_weights$point_indices)
  reference_cell <- as.integer(reference_cell)

  if (reference_by == "sce_index") {
    reference_index <- reference_cell
    reference_position <- match(reference_index, point_indices)
  } else {
    reference_position <- reference_cell
    reference_index <- point_indices[reference_position]}

  # extract the stored neighborhood
  neighbor_indices <- as.integer(spatial_weights$neighbors[[reference_position]])
  neighbor_weights <- as.numeric(spatial_weights$weights[[reference_position]])
  n_neighbors <- length(neighbor_indices)

  stored_distances <- if (is.null(spatial_weights$distances)) rep(NA_real_, n_neighbors) else as.numeric(spatial_weights$distances[[reference_position]])

  reference_point <- as.numeric(distance_coordinates[reference_index, ])
  distance_to_reference <- function(indices) {
    sqrt(rowSums(sweep(distance_coordinates[indices, , drop = FALSE], 2L, reference_point, "-")^2))
  }

  # plot data
  background_indices <- setdiff(seq_len(nrow(cell_data)), c(reference_index, neighbor_indices))
  background_data <- data.frame(x = plot_coordinates[background_indices, 1L], y = plot_coordinates[background_indices, 2L])

  neighbor_data <- data.frame(
    sce_index = neighbor_indices,
    x = plot_coordinates[neighbor_indices, 1L], y = plot_coordinates[neighbor_indices, 2L],
    cell_type = cell_types[neighbor_indices], weight = neighbor_weights,
    stored_distance = stored_distances,
    recomputed_distance = distance_to_reference(neighbor_indices))

  reference_x <- plot_coordinates[reference_index, 1L]
  reference_y <- plot_coordinates[reference_index, 2L]

  # plot
  plot <- ggplot2::ggplot() +
    ggplot2::geom_point(data = background_data, ggplot2::aes(x, y),
                        color = background_color, size = background_size, alpha = background_alpha)

  if (show_links) {
    plot <- plot + ggplot2::geom_segment(
      data = transform(neighbor_data, xend = reference_x, yend = reference_y),
      ggplot2::aes(x, y, xend = xend, yend = yend),
      color = link_color, alpha = link_alpha, linewidth = link_width)}

  if (neighbor_color_by == "weight") {
    # Draw low-weight points first so high-weight points remain visible.
    neighbor_data <- neighbor_data[order(neighbor_data$weight), ]
    plot <- plot +
      ggplot2::geom_point(data = neighbor_data, ggplot2::aes(x, y, color = weight), size = neighbor_size, alpha = 0.95) +
      ggplot2::scale_color_viridis_c(
        name = weight_name, limits = weight_limits, oob = scales::squish,
        guide = ggplot2::guide_colourbar(
          position = "bottom", direction = "horizontal",
          theme = ggplot2::theme(
            legend.key.height = grid::unit(0.4, "cm"), legend.key.width = grid::unit(1, "null"),
            legend.title = ggplot2::element_text(hjust = 0.5), legend.title.position = "top"
          )
        )
      )
  } else {
    all_cell_types <- sort(unique(cell_types[!is.na(cell_types)]))
    color_palette <- stats::setNames(rep(cell_type_palette, length.out = length(all_cell_types)), all_cell_types)
    plot <- plot +
      ggplot2::geom_point(data = neighbor_data, ggplot2::aes(x, y, color = cell_type),
                          size = neighbor_size, alpha = 0.95) +
      ggplot2::scale_color_manual(values = color_palette, name = "Neighbor cell type") +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2.5, alpha = 1)))
  }

  plot <- plot + ggplot2::geom_point(
    data = data.frame(x = reference_x, y = reference_y), ggplot2::aes(x, y),
    shape = 21, fill = reference_color, color = "black", stroke = 0.6, size = reference_size
  )


  plot <- plot +
    ggplot2::labs(title = title, x = coords[1L], y = coords[2L]) +
    #ggplot2::labs(title = title, x = coords[1L], y = coords[2L]) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 8),
      axis.title = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 8))

  if (zoom) {
    zoom_x <- range(c(reference_x + neighbor_data$x*0.3, reference_x - neighbor_data$x*0.3))
    zoom_y <- range(c(reference_y + neighbor_data$y*0.3, reference_y - neighbor_data$y*0.3))
    padding <- 0.10 * max(diff(zoom_x), diff(zoom_y))
    plot <- plot + ggplot2::coord_fixed(xlim = zoom_x + c(-1, 1) * padding,
                                        ylim = zoom_y + c(-1, 1) * padding, expand = FALSE)
  } else { plot <- plot + ggplot2::coord_fixed() }



  plot
}
