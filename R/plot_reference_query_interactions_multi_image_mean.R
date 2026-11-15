#' function to plot SES score trendline across permutation scales
#' 
#' @param multiscale_results results computed using `simulate_multiscale_mvlclq_fast_downsample`
#' @param sce_obj SpatialExperiment or SingleCellExperiment object
#' @param type_column Column name for cell type
#' @param image image to be plotted
#' @param cell_subset optional - reference cell types to keep
#' @param color_by optional - select which column used for the color grouping
#' @param significance_threshold significance threshold before correction. inside the function, this threshold is corrected for multiple testing
#' @param title title
#' @param facet facet by reference cell type with free y axes
#' @param show_trend currently unused
#' @param show_points overlay the individual reference-cell SES values
#' @param show_error draw the mean SES trend line with +/- standard error bars
#' @param custom_colors optional colors
#' @returns ggplot object
#' @export

plot_reference_query_interactions_multi_image_mean <- function(multiscale_results,
                                                               sce_obj,
                                                               type_column,
                                                               image = "image_1",
                                                               cell_subset = NULL,
                                                               color_by = NULL,
                                                               significance_threshold = 0.05,
                                                               title = "Multi-scale MVLCLQ Analysis",
                                                               facet = TRUE,
                                                               show_trend = TRUE,
                                                               show_points = FALSE,
                                                               show_error = TRUE,
                                                               custom_colors = NULL) {
  
  image_result <- multiscale_results[["results_by_image"]][[as.character(image)]]
  scales <- image_result[["scales"]]
  point_indices <- image_result[["point_indices"]]
  interaction_names <- image_result[["interaction_names"]]
  results_by_scale <- image_result[["results_by_scale"]]
  
  if (!is(sce_obj, "SpatialExperiment") && !is(sce_obj, "SingleCellExperiment")) {
    stop("Input must be a SpatialExperiment or SingleCellExperiment object")
  }
  cell_data <- as.data.frame(colData(sce_obj))
  
  # we make one data frame per scale x interaction, and then stacked them in the same order as before
  plot_data <- do.call(rbind, lapply(scales, function(scale) {
    do.call(rbind, lapply(interaction_names, function(interaction_name) {
      interaction_result <- results_by_scale[[as.character(scale)]]$interactions[[interaction_name]]
      #interaction_result <- results_by_scale$interactions[[interaction_name]]
      data.frame(scale = scale, point_index = point_indices,
        SES = interaction_result$ses_scores, cell_type = cell_data[point_indices, type_column],
        interaction = interaction_name)}))
  }))
  
  # revised label - reference-query label, e.g. A-B
  plot_data$color_group <- if (!is.null(color_by) && color_by %in% colnames(cell_data)) {
    cell_data[plot_data$point_index, color_by]
  } else {paste(plot_data$cell_type, plot_data$interaction, sep = "-")}
  
  if (!is.null(cell_subset)) plot_data <- plot_data[plot_data$cell_type %in% cell_subset, ]
  if (nrow(plot_data) == 0) stop("No data remaining after subsetting. Check your cell_subset parameter.")
  
  obs_summary <- plot_data %>%
    group_by(scale, cell_type, interaction, color_group) %>%
    summarise(mean_SES = mean(SES, na.rm = TRUE),
      std_SES = sd(SES, na.rm = TRUE), se_SES = std_SES / sqrt(length(SES)),
      n_points = n(), .groups = "drop")
  
  # corrected
  #n_test <- length(interaction_names)
  sig_threshold_z <- qnorm(1 - (significance_threshold / length(interaction_names)) / 2)
  
  p <- ggplot()
  
  if (show_error && nrow(obs_summary) > 0) {
    summary_data_clean <- obs_summary[is.finite(obs_summary$mean_SES) &
                                        is.finite(obs_summary$std_SES) &
                                        is.finite(obs_summary$se_SES) &
                                        obs_summary$std_SES >= 0, ]
    if (nrow(summary_data_clean) > 0) {
      p <- p +
        geom_errorbar(
          data = summary_data_clean,
          aes(x = scale, ymin = mean_SES - se_SES, ymax = mean_SES + se_SES, color = color_group),
          width = 0.1 * (max(scales) - min(scales)) / length(scales),
          position = position_dodge(width = 50), linewidth = 0.8, alpha = 0.7
        ) +
        geom_line(data = summary_data_clean, aes(x = scale, y = mean_SES, color = color_group), linewidth = 1.2) +
        geom_point(data = summary_data_clean, aes(x = scale, y = mean_SES, color = color_group), size = 2)
    }
  }
  
  if (show_points) {
    plot_data_clean <- plot_data[is.finite(plot_data$SES), ]
    if (nrow(plot_data_clean) > 0) {
      p <- p + geom_point(
        data = plot_data_clean, aes(x = scale, y = SES, color = color_group),
        alpha = 0.3, position = position_jitter(width = 0.1, height = 0)
      )
    }
  }
  
  p <- p +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.2) +
    geom_hline(yintercept = c(-sig_threshold_z, sig_threshold_z), color = "red", linewidth = 0.6, linetype = "dotted")
    #geom_hline(yintercept = c(-sig_threshold_z, sig_threshold_z), linewidth = 1)
  
  if (!is.null(custom_colors)) p <- p + scale_color_manual(values = custom_colors)
  if (facet && length(unique(plot_data$cell_type)) > 1) p <- p + facet_wrap(~cell_type, scales = "free_y")
  
  p +
    labs(title = title, x = "Scale", y = "SES Score", color = if (is.null(color_by)) "Reference-Query" else color_by) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 11, color = "black", angle = 90, vjust = 0.5, hjust = 0),
      axis.text.y = element_text(size = 11, color = "black"),
      axis.title = element_text(size = 13),
      plot.title = element_text(size = 14),
      panel.grid.major = element_line(linewidth = 0.1, colour = "gray90"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.1),
      legend.position = "bottom",
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11)
    ) +
    guides(color = guide_legend(nrow = 1))
}