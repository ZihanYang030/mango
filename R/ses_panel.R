#' function to organize plots into figure
#' 
#' @param sce SpatialExperiment or SingleCellExperiment object
#' @param results multi-scale permutation results
#' @param reference reference cell type
#' @param query query cell types
#' @param type_column Column name for cell type
#' @param coords spatial coordinates
#' @param scale permutation scale to be plotted (for spatial heatmap)
#' @param image_name single image to be plotted
#' @param point_size plotting parameter
#' @param legend_ncol number of column - specify this when there is a large number of cell types
#' @returns A \code{patchwork} object, which inherits from \code{ggplot}.
#' @importFrom patchwork wrap_plots
#' @importFrom ggplot2 guides guide_legend theme element_blank element_text
#' @importFrom grid unit
#' @export

ses_panel <- function(sce, results, reference, queries,
                      type_column, coords, scale, image_name,
                      point_size = 1, legend_ncol = 2) {
  
  # newly added - for organizing ses score + heatmap panels
  # just an organizing function, not an actual plotting function
  heatmap_one <- function(q) {
    Scatter_Heatmap(
      sce, type_column = type_column, results,
      coords = coords, scale = scale, image_name = image_name,
      query_celltype = q, reference_celltype = reference,
      title = paste0(reference, "-", gsub("_", "-", q), " SES Score"),
      point_size = point_size)
  }
  
  # each query on its own, then the joint one
  heatmaps <- lapply(c(queries, paste(queries, collapse = "_")), heatmap_one)
  
  trendline <- plot_reference_query_interactions_multi_image_mean(
    results, sce, type_column = type_column, image = image_name)
  
  basic <- Plot_Basic(sce, type_column = type_column, highlight_celltypes = NULL,
                      coords = coords, title = NULL, point_size = point_size) +
    guides(colour = guide_legend(ncol = legend_ncol,
                                 override.aes = list(size = 2))) +
    theme(legend.title    = element_blank(),
          legend.text     = element_text(size = 6),
          legend.key.size = unit(0.3, "cm"))
  
  no_axes <- theme(axis.title = element_blank(), axis.text = element_blank(),
                   axis.ticks = element_blank(), panel.grid = element_blank())
  
  trendline | (wrap_plots(c(heatmaps, list(basic)), ncol = 2) & no_axes)
}