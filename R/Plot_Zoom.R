#' function: Plot_Zoom, spatial visualization of a zoom-in square region (2D)
#' 
#' @param sce_obj SpatialExperiment or SingleCellExperiment object
#' @param type_column Column name for cell type
#' @param center_cell the cell at the center of the zoom-in square, either a column index or a cell name (colnames of sce_obj)
#' @param side_length side length of the zoom-in square, same unit as the spatial coordinates
#' @param highlight_celltypes takes a list of cell types, only color these cells to highlight
#' @param coords spatial coordinates column
#' @param cell_size dot size could vary for different cells, enter column name for cell size
#' @param point_color plotting parameter, dot color for cells (if NULL, use the default color palette)
#' @param title plot title, default: "Zoom-in Spatial Plot"
#' @param point_size uniform point size for every cell
#' @param mark_center mark the center cell with a black cross, default: FALSE
#' @returns A ggplot object
#' @export

Plot_Zoom <- function(sce_obj, type_column, center_cell, side_length,
                      highlight_celltypes = NULL,
                      coords = c("x", "y"), cell_size = NULL,
                      point_color = NULL, title = NULL, point_size = 0.5,
                      mark_center = FALSE) {
  
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
  
  # locate the center cell, accepts either a column index or a cell name
  if (is.character(center_cell)) {center_cell <- match(center_cell, colnames(sce_obj))}
  if (is.na(center_cell) | center_cell < 1 | center_cell > ncol(sce_obj)) {
    stop("center_cell is not found in sce_obj")}
  
  x_center <- x_coord[center_cell]
  y_center <- y_coord[center_cell]
  half_side <- side_length / 2
  
  # data frame for plotting, point size carried along so it stays matched to the cells kept
  plot_data <- data.frame(x = x_coord, y = y_coord, cell_type = as.factor(cell_types),
                          point_size = point_size)
  #plot_data <- plot_data[complete.cases(plot_data), ]
  plot_data <- plot_data[stats::complete.cases(plot_data), ]
  
  # levels come from the whole object, so colors stay consistent with the full tissue plot
  all_cell_types <- levels(plot_data$cell_type)
  
  # keep only the cells inside the square centered on the center cell
  plot_data <- plot_data[abs(plot_data$x - x_center) <= half_side &
                           abs(plot_data$y - y_center) <= half_side, ]
  
  if (nrow(plot_data) == 0) {
    stop("no cell falls inside the zoom-in square, try a larger side_length")}
  
  point_size <- plot_data$point_size
  
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
  
  if (is.null(title)) { title <- "Zoom-in Spatial Plot" }
  
  # plot
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = cell_type)) +
    ggplot2::geom_point(size = point_size, alpha = 0.8) +
    ggplot2::scale_color_manual(values = final_colors, name = "Cell Type") +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1))) +
    ggplot2::labs(title=title, x = coords[1], y = coords[2]) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      aspect.ratio = NULL,
      axis.text = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 8)
    ) +
    ggplot2::coord_fixed(ratio = 1,
                         xlim = c(x_center - half_side, x_center + half_side),
                         ylim = c(y_center - half_side, y_center + half_side))
  
  # mark the center cell
  if (mark_center) {
    plot <- plot + ggplot2::annotate("point", x = x_center, y = y_center,
                                     shape = 4, size = 2, color = "black")}
  
  return(plot)
}