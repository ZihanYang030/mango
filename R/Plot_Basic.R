
#' function: Plot_Basic, spatial visualization of cell types and locations (2D version)
#' @param sce_obj SpatialExperiment or SingleCellExperiment object
#' @param type_column Column name for cell type
#' @param highlight_celltypes takes a list of cell types, only color these cells to highlight
#' @param coords spatial coordinates column
#' @param cell_size dot size could vary for different cells, enter column name for cell size
#' @param point_color plotting parameter, dot color for cells (if NULL, use the default color palette)
#' @param title plot title, default: "Spatial Plot"
#' @param point_size uniform point size for every cell
#' @returns A ggplot object.
#' @export

Plot_Basic <- function(sce_obj, type_column, highlight_celltypes = NULL,
                       coords = c("x", "y"), cell_size = NULL,
                       point_color = NULL, title = NULL, point_size = 0.5) {
  
  # extract data
  if (methods::is(sce_obj, "SpatialExperiment")) {
    x_coord <- SpatialExperiment::spatialCoords(sce_obj)[, 1]
    y_coord <- SpatialExperiment::spatialCoords(sce_obj)[, 2]
  } else if (methods::is(sce_obj, "SingleCellExperiment")) {
    x_coord <- SummarizedExperiment::colData(sce_obj)[[coords[1]]]
    y_coord <- SummarizedExperiment::colData(sce_obj)[[coords[2]]]}
  
  cell_types <- SummarizedExperiment::colData(sce_obj)[[type_column]]
  cell_data <- SummarizedExperiment::colData(sce_obj)
  
  if (!is.null(cell_size)) {point_size <- cell_data@listData[[cell_size]]}
  
  # data frame for plotting
  plot_data <- data.frame(x = x_coord, y = y_coord, cell_type = as.factor(cell_types))
  #plot_data <- plot_data[complete.cases(plot_data), ]
  plot_data <- plot_data[stats::complete.cases(plot_data), ]

  all_cell_types <- levels(plot_data$cell_type)
  
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
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = cell_type)) +
    ggplot2::geom_point(size = point_size, alpha = 0.8) +
    ggplot2::scale_color_manual(values = final_colors, name = "Cell Type") +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1))) + 
    ggplot2::labs(title=title, x = coords[1], y = coords[2]) + 
    ggplot2::theme_minimal() +
    ggplot2::theme(
      aspect.ratio = NULL,
      axis.text = ggplot2::element_text(size = 16),
      axis.title = ggplot2::element_text(size = 20),
      plot.title = ggplot2::element_text(size = 22, hjust = 0.5),
      legend.title = ggplot2::element_text(size = 18),
      legend.text = ggplot2::element_text(size = 14)
    ) + ggplot2::coord_fixed(ratio = 1)
  
  return(plot)
}
