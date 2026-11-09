#' plotting function for single-cell-level SES scores and expression levels
#' 
#' @param sce_obj SpatialExperiment or SingleCellExperiment object in the original simulation cell order
#' @param type_column Column name for cell type
#' @param simulation_results multi-scale single-cell level simulation results
#' @param coords column names for spatial coordinates
#' @param cell_name colData column with cell IDs matching expression_matrix columns; required for expression plotting
#' @param scale specific permutation scale that the SES results will be plotted
#' @param image_name image ID to plot
#' @param query_celltype query cell type name, put "B_C" for A-B-C interaction
#' @param expression_celltype cell type where expression level will be plotted, this may or may not be the query cell type
#' @param target_gene gene which expression level is plotted
#' @param expression_matrix numeric matrix with genes as rows and cell IDs as columns
#' @param reference_celltype reference cell type name
#' @param point_color optional colors for the SES gradient
#' @param title plot title
#' @param point_size plotting point size
#' @param image_column `colData` column containing image IDs; supply it when `simulation_results` contains multiple images.
#' @returns ggplot object
#' @importFrom methods is
#' @importFrom SummarizedExperiment colData
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom stats complete.cases
#' @importFrom grid unit
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradientn guide_colorbar labs theme_minimal theme element_text coord_fixed
#' @export

Scatter_Heatmap <- function(sce_obj, type_column, simulation_results, coords = c("x", "y"), cell_name = NULL,
                            scale = 1000, image_name = "image1", query_celltype = "B", expression_celltype = NULL,
                            target_gene = NULL, expression_matrix = NULL,
                            reference_celltype = NULL, point_color = NULL, title = NULL, point_size = 0.5,
                            image_column = NULL) {

  if (length(coords) != 2L) stop("coords must contain exactly two coordinate names")

  # extract data from sce or spe object
  if (is(sce_obj, "SpatialExperiment")) {
    cell_data <- colData(sce_obj)
    coords_data <- spatialCoords(sce_obj)
    if (ncol(coords_data) < 2L) stop("At least two spatial coordinates are required")
    
    x_coord <- coords_data[, 1]
    y_coord <- coords_data[, 2]
    
  } else if (is(sce_obj, "SingleCellExperiment")) {

    cell_data <- colData(sce_obj)
    if (!all(coords %in% colnames(cell_data))) {stop("Specified coordinates not found in colData")}
    
    x_coord <- cell_data[[coords[1L]]]
    y_coord <- cell_data[[coords[2L]]]
    
  } else {
    stop("please input a SpatialExperiment or SingleCellExperiment")
  }
  
  expression_requested <- any(!is.null(expression_celltype), !is.null(target_gene), !is.null(expression_matrix))
  
  if (expression_requested &&
      (is.null(cell_name) ||
       is.null(expression_celltype) ||
       is.null(target_gene) ||
       is.null(expression_matrix))) {
    stop("pls provide cell_name, expression_celltype, target_gene, and expression_matrix together")
  }
  
  # cell_types <- colData(sce_obj)[[type_column]]
  if (!type_column %in% colnames(cell_data)) stop("type_column not found in colData")
  cell_types <- cell_data[[type_column]]
  
  #ses_scores from simulation results
  image_name <- as.character(image_name)
  scale <- as.character(scale)
  query_celltype <- as.character(query_celltype)
  #ses_scores <- simulation_results[["results_by_image"]][[image_name]][["results_by_scale"]][[scale]][["interactions"]][[query_celltype]][["ses_scores"]]
  image_result <- simulation_results[["results_by_image"]][[image_name]]
  if (is.null(image_result)) stop("image_name not found")
  ses_scores <- image_result[["results_by_scale"]][[scale]][["interactions"]][[query_celltype]][["ses_scores"]]
  
  if (is.null(ses_scores)) {stop("no ses-scores detected")}
  
  #plot_data <- data.frame( x = x_coord, y = y_coord, cell_type = as.factor(cell_types))
  plot_data <- data.frame(x = x_coord, y = y_coord,
                          cell_type = as.factor(cell_types), ses_score = NA_real_)
  
  if (is.null(reference_celltype)) { stop("no reference type provided") }
  
  #ref_indices <- which(plot_data$cell_type == reference_celltype)
  ref_indices <- image_result[["point_indices"]]

  #if(length(ref_indices) != length(ses_scores)) {stop("number of references does not equal number of ses scores, check inputs")}
  if (is.null(ref_indices) || length(ref_indices) == 0L ||
      length(ref_indices) != length(ses_scores) || anyNA(ref_indices) ||
      any(ref_indices < 1L | ref_indices > ncol(sce_obj))) {
    stop("Reference indices and SES scores do not match the original sce_obj")
  }
  if (!all(as.character(plot_data$cell_type[ref_indices]) %in% reference_celltype)) {
    stop("reference_celltype does not match the stored reference indices")
  }

  plot_data$ses_score[ref_indices] <- ses_scores

  #plot_data <- plot_data[complete.cases(plot_data), ]
  if (expression_requested) {
    if (!cell_name %in% colnames(cell_data)) stop("cell_name not found in colData")
    plot_data$cell_id <- as.character(cell_data[[cell_name]])}
  
  if (is.null(image_column) && length(simulation_results[["results_by_image"]]) > 1L) {
    stop("Provide image_column when simulation_results contains multiple images")}
  
  if (!is.null(image_column)) {
    if (!image_column %in% colnames(cell_data)) stop("image_column not found in colData")
    image_mask <- !is.na(cell_data[[image_column]]) &
      as.character(cell_data[[image_column]]) == image_name
    if (!all(image_mask[ref_indices])) stop("Stored reference indices are outside image_name")
    plot_data <- plot_data[image_mask, , drop = FALSE]
  }
  plot_data <- plot_data[complete.cases(plot_data[c("x", "y", "cell_type")]),
                         , drop = FALSE]
  

  if (!is.null(point_color) && length(point_color) < 2L) { stop("point_color must contain at least two colors") }
  color_palette <- if (is.null(point_color)) {
    viridisLite::viridis(256, begin = 0, end = 1)
  } else {
    grDevices::colorRampPalette(point_color)(256)
  }
  

  if (is.null(title)) { title <- paste("Spatial Plot -", query_celltype, "SES-scores") }
  
  #ses_scores_actual <- plot_data$ses_score[!is.na(plot_data$ses_score)]
  if (any(!is.finite(plot_data$ses_score[!is.na(plot_data$ses_score)]))) {
    stop("SES scores must be finite")
  }
  min_ses <- min(plot_data$ses_score, na.rm = TRUE)
  max_ses <- max(plot_data$ses_score, na.rm = TRUE)
  
  if (expression_requested) {
    expression_plotdata <- plot_data[plot_data$cell_type %in% expression_celltype, , drop = FALSE]
    if (!target_gene %in% rownames(expression_matrix)) { stop("wrong gene name, change one")}
    expression_ids <- expression_plotdata$cell_id
    if (anyNA(expression_ids) || anyDuplicated(expression_ids) || !all(expression_ids %in% colnames(expression_matrix))) {
      stop("expression cell IDs must uniquely match expression_matrix columns")}
    expression_plotdata$expression_level <- as.numeric(as.matrix(
      expression_matrix[target_gene, expression_ids, drop = FALSE]
    ))
    if (all(is.na(expression_plotdata$expression_level))) { stop("no expression to plot") }
  }
  
  plot <- ggplot(plot_data, aes(x = x, y = y)) +
    geom_point(data = plot_data[is.na(plot_data$ses_score), ], aes(x = x, y = y), 
               color = "#e8e8e8", size = point_size, alpha = 0.6) +
    geom_point(data = plot_data[!is.na(plot_data$ses_score), ], aes(x = x, y = y, color = ses_score),
               # size = point_size, alpha = 0.8) +
               size = if (expression_requested) point_size + 2.5 else point_size,
               alpha = 0.8) +
    scale_color_gradientn(
      name = "SES",
      # colors = viridisLite::viridis(256, begin = 0, end = 1),
      colors = color_palette,
      values = seq(0, 1, length.out = 256)^1.5,
      # limits = c(-2, max_ses),
      limits = c(min(-2, min_ses), max(0, max_ses)),
      oob = scales::squish,
      na.value = "#e8e8e8",
      guide = guide_colorbar(direction = "horizontal", title.position = "top",
        title.hjust = 0.5, barwidth = unit(0.6, "npc"), barheight = unit(0.35, "cm"))
    ) +
    labs(title = title, x = coords[1], y = coords[2]) + 
    theme_minimal() +
    theme(
      aspect.ratio = NULL,
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 20),
      plot.title = element_text(size = 22, hjust = 0.5),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 14),
      legend.position = "bottom"
    ) +
    coord_fixed(ratio = 1) 
  
  if(expression_requested) {
    plot <- plot + ggnewscale::new_scale_color() + 
      geom_point(data = expression_plotdata,
                 aes(x = x, y = y, color = expression_level),
                 size = point_size, alpha = 0.8) +
      scale_color_gradientn(
        name = "Expression Level",
        colors = viridisLite::viridis(option = "F", 256, begin = 0, end = 1, direction = -1),
        values = seq(0, 1, length.out = 256)^2,
        #limits = c(min(-0.1 * max_expr_level, min_expr_level), max_expr_level),
        oob = scales::squish,
        na.value = "#e8e8e8",
        guide = guide_colorbar(
          direction = "vertical",
          title.position = "top",
          title.hjust = 0.5,
          barwidth = unit(0.35, "cm"),
          barheight = unit(4, "cm")
        )
      ) 
  }
  
  return(plot)
}
