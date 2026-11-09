
#' check cell count & cell proportion 
#'
#' @param sce_obj SpatialExperiment or SingleCellExperiment object
#' @param type_column Column name for cell type
#' @param image_column Column name for image index
#' @returns list of cell counts and proportions per image
#' @importFrom SummarizedExperiment colData
#' @export

calculate_occurrence_probability_multi_image <- function(sce_obj, type_column, image_column = "slide") {
  
  cell_data <- colData(sce_obj)
  imgs <- unique(cell_data[[image_column]])
  
  all_image_results <- list()
  for (i in 1:length(imgs)) {
    current_img <- imgs[i]
    
    cell_data_subset <- cell_data[cell_data[[image_column]] == current_img, ]
    
    type_cnts <- table(cell_data_subset[[type_column]])
    #total_count <- sum(type_cnts)
    probs <- type_cnts/sum(type_cnts)
    names(probs) <- names(type_cnts)
    
    all_image_results[[as.character(current_img)]] <- list(image_id = current_img,
      type_counts = type_cnts, total_count = sum(type_cnts), probabilities = probs)
  }
  
  all_image_results
}
