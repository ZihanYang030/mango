
#' function: Plot_Basic_3D, spatial visualization of cell types and locations (3D version)
#' @param sce_obj SpatialExperiment or SingleCellExperiment object
#' @param type_column Column name for cell type
#' @param highlight_celltypes takes a list of cell types, only color these cells to highlight
#' @param coords spatial coordinates column
#' @param cell_size dot size could vary for different cells, enter column name for cell size
#' @param point_color plotting parameter, dot color for cells (if NULL, use the default color palette)
#' @param title plot title, default: "Spatial Plot"
#' @param point_size uniform point size for every cell
#' @returns A Plotly object.
#' @importFrom SummarizedExperiment colData
#' @importFrom stats complete.cases
#' @export

Plot_Basic_3D <- function(sce_obj,
                          type_column,
                          highlight_celltypes = NULL,
                          coords = c("x", "y", "z"),
                          cell_size = NULL,
                          point_color = NULL,
                          title = NULL,
                          point_size = 0.5) {

  # extract data
  if (methods::is(sce_obj, "SpatialExperiment")) {
    x_coord <- SpatialExperiment::spatialCoords(sce_obj)[, 1]
    y_coord <- SpatialExperiment::spatialCoords(sce_obj)[, 2]
    z_coord <- SpatialExperiment::spatialCoords(sce_obj)[, 3]
  } else if (methods::is(sce_obj, "SingleCellExperiment")) {
    x_coord <- colData(sce_obj)[[coords[1]]]
    y_coord <- colData(sce_obj)[[coords[2]]]
    z_coord <- colData(sce_obj)[[coords[3]]]
    }

  cell_types <- colData(sce_obj)[[type_column]]
  cell_data <- colData(sce_obj)

  if (!is.null(cell_size)) {point_size <- cell_data@listData[[cell_size]]}

  # data frame for plotting
  plot_data <- data.frame(x = x_coord, y = y_coord, z = z_coord, cell_type = as.factor(cell_types))
  plot_data <- plot_data[complete.cases(plot_data), ]

  all_cell_types <- levels(plot_data$cell_type)

  # Set up colors
  if (!is.null(point_color)) { color_palette <- point_color } else {
    color_palette <- c("#F9C74F", "#67a8cd", "#ff9d9f", "#d1edc2", "#deb6f0", "#a2d2e7", "#cf9f88",
                       "#6fb3a8", "#b3e19b", "#EBD7A5", "#f36569", "#ff8831", "#ffec99",
                       "#A6CA83", "#7d67c7", "#ffc17f", "#466399", "#ED729B", "#947876",
                       "#fbabff", "#aa6dea", "#b0a079", "#92f1e8", "#dff7a5", "#b7cbb5",
                       "#e07a5f", "#BD9E39", "#006073", "#099396", "#33A02C", "#d28f38", "#52f02b",
                       "#faf219", "#dc19fa")}

  if (length(color_palette) < length(all_cell_types)) {
    color_palette <- rep(color_palette, length.out = length(all_cell_types))
  }

  # assign colors to cell types
  names(color_palette) <- all_cell_types[1:length(all_cell_types)]

  # highlight cell types
  if (!is.null(highlight_celltypes)) {

    final_colors <- rep("#e8e8e8", length(all_cell_types))  # Light gray for other types
    names(final_colors) <- all_cell_types
    # original colors for highlighted cell types
    final_colors[highlight_celltypes] <- color_palette[highlight_celltypes]

  } else {
    final_colors <- color_palette
  }

  if (is.null(title)) { title <- "Spatial Plot" }

  # plot
  plot <- plotly::plot_ly(
    data = plot_data, x = ~x, y = ~y, z = ~z,
    color = ~cell_type, colors = final_colors,
    type = "scatter3d", mode = "markers",
    marker = list(size = point_size)
  )
  plot <- plotly::layout(
    plot,
    title = title,
    scene = list(xaxis = list(title = coords[1]), yaxis = list(title = coords[2]), zaxis = list(title = coords[3]))
  )

  return(plot)
}
