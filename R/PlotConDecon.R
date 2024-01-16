#' Plot ConDecon_obj
#' @description Plot ConDecon's predicted relative cell probability z-scores for each sample
#' @importFrom graphics plot
#' @import ggplot2
#' @param ConDecon_obj ConDecon object (output from RunConDecon)
#' @param umap 2-dimensional coordinates of reference single-cell data (default = first 2 dimensions of latent space)
#' @param samples Vector of query bulk sample(s) to plot (default = NULL, plot all samples)
#' @param pt.size size of the points plotted (default = 1)
#' @param cells Vector of cells from the reference single-cell data to plot (default = NULL, plot all cells)
#' @param title_names Title of each plot (default = samples, use NULL for no title)
#' @param relative Logical indicating whether to plot the relative abundance (Z-scores) or the raw abundance (default = TRUE)
#' @param upper_quant When relative = FALSE, plot the raw abundance removing the upper quartile or the long tails of the distribution (default = 0.95)
#' @param lower_quant When relative = FALSE, plot the raw abundance removing the lower quartile or the long tails of the distribution (default = 0.05)
#'
#' @return Plot the z-score cell probability values for each query bulk sample
#' @export
#'
#' @examples
#' data(counts_gps)
#' data(latent_gps)
#' data(bulk_gps)
#' data(variable_genes_gps)
#'
#' # For this example, we will reduce the training size to max.iter = 50 to reduce run time
#' ConDecon_obj = RunConDecon(counts = counts_gps, latent = latent_gps, bulk = bulk_gps,
#' variable.features = variable_genes_gps, max.iter = 50)
#'
#' # Plot ConDecon relative cell prob
#' PlotConDecon(ConDecon_obj)
PlotConDecon <- function(ConDecon_obj,
                         umap = ConDecon_obj$TrainingSet$latent[,1:2],
                         samples = NULL,
                         pt.size = 1,
                         cells = NULL,
                         title_names = samples,
                         relative = TRUE,
                         upper_quant = 0.95,
                         lower_quant = 0.05){

  umap <- data.frame(Dim_1 = umap[,1], Dim_2 = umap[,2], row.names = row.names(umap))

  #Identify the bulk samples to plot
  if(is.null(samples)){
    samples <- colnames(ConDecon_obj$Normalized_cell.prob)
    #Bulk samples name(s) are within column names of ConDecon Object
    #If samples represent indices:
  } else if(is.numeric(samples)){
    if(sum(samples %in% 1:ncol(ConDecon_obj$Normalized_cell.prob))>0){
      if(sum(samples %in% 1:ncol(ConDecon_obj$Normalized_cell.prob)) != length(samples)){
        warning("Some of the indices in 'samples' exceed the total number of samples in the ConDecon object\n")
      }
      samples <- samples[samples %in% 1:ncol(ConDecon_obj$Normalized_cell.prob)]
    } else{
      message("None of the 'samples' were found in ConDecon object")
      return(NULL)
    }
    #If samples represent sample names:
  } else if(is.character(samples)){
    if(sum(as.character(samples) %in% colnames(ConDecon_obj$Normalized_cell.prob))>0){
      if(sum(as.character(samples) %in% colnames(ConDecon_obj$Normalized_cell.prob)) != length(samples)){
        warning("Some of the names in 'samples' are not found in the ConDecon object\n")
      }
      samples <- samples[samples %in% colnames(ConDecon_obj$Normalized_cell.prob)]
    } else {
      message("None of the 'samples' were found in ConDecon object")
      return(NULL)
    }
  } else{
    message("'samples' must be a vector")
    return(NULL)
  }

  if(!is.vector(title_names) & !is.null(title_names)){
    message("'title_names' must be a character vector or NULL")
    return(NULL)
  } else if(length(title_names) == length(samples)){
    title_names <- as.character(title_names)
    names(title_names) <- as.character(samples)
  } else if(length(title_names) < length(samples) & !is.null(title_names)){
    warning("There are fewer items in vector 'title_names' than 'samples' being plotted\n")
    title_names <- NULL
  } else if(length(title_names) > length(samples)){
    warning("There are more items in vector 'title_names' than 'samples' being plotted\n")
    title_names <- NULL
  }

  #Identify the cells to plot
  if(is.null(cells)){
    cells <- row.names(ConDecon_obj$Normalized_cell.prob)
    #Cells name(s) are within row names of ConDecon Object
    #If cells represent indices:
  } else if(is.numeric(cells)){
    if(sum(cells %in% 1:nrow(ConDecon_obj$Normalized_cell.prob))>0){
      if(sum(cells %in% 1:nrow(ConDecon_obj$Normalized_cell.prob)) != length(cells)){
        warning("Some of the indices in 'cells' exceed the total number of cells in the ConDecon object\n")
      }
      cells <- cells[cells %in% 1:nrow(ConDecon_obj$Normalized_cell.prob)]
    } else{
      message("None of the 'cells' were found in ConDecon object")
      return(NULL)
    }
    #If cells represent cell names:
  } else if(is.character(cells)){
    if(sum(as.character(cells) %in% row.names(ConDecon_obj$Normalized_cell.prob))>0){
      if(sum(as.character(cells) %in% row.names(ConDecon_obj$Normalized_cell.prob)) != length(cells)){
        warning("Some of the names in 'cells' are not found in the ConDecon object\n")
      }
      cells <- cells[cells %in% row.names(ConDecon_obj$Normalized_cell.prob)]
    } else {
      message("None of the 'cells' were found in ConDecon object")
      return(NULL)
    }
  } else{
    message("'cells' must be a vector")
    return(NULL)
  }
  
  if(relative){
    for (i in samples){
      Dim_1 <- Dim_2 <- NULL # Setting the variables to NULL to appease R-CMD-check
      #https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
      ##scale_color_gradient2(low = "blue", high = "red", mid = "gray")
      graphics::plot(ggplot2::ggplot(umap[cells,], ggplot2::aes(x=Dim_1, y=Dim_2, color = ConDecon_obj$Relative_cell.prob[cells, i])) +
                       ggplot2::geom_point(size = pt.size) +
                       ggplot2::scale_color_gradient2(low = "#011F5B", high = "#990000", mid = "light gray") +
                       ggplot2::labs(color = "Cell Prob\nz-score", title = ifelse(is.null(title_names), "", title_names[as.character(i)])) +
                       ggplot2::theme_classic())
    }
  } else{
    for (i in samples){
      Dim_1 <- Dim_2 <- NULL # Setting the variables to NULL to appease R-CMD-check
      #https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
      ##scale_color_gradient2(low = "blue", high = "red", mid = "gray")
      norm_cell_prob <- ConDecon_obj$Normalized_cell.prob
      norm_cell_prob[norm_cell_prob < quantile(as.vector(ConDecon_obj$Normalized_cell.prob), probs = lower_quant)] <- quantile(as.vector(ConDecon_obj$Normalized_cell.prob), probs = lower_quant)
      norm_cell_prob[norm_cell_prob > quantile(as.vector(ConDecon_obj$Normalized_cell.prob), probs = upper_quant)] <- quantile(as.vector(ConDecon_obj$Normalized_cell.prob), probs = upper_quant)
      norm_cell_prob <- t(t(norm_cell_prob)/colSums(norm_cell_prob))
      graphics::plot(ggplot2::ggplot(umap[cells,], ggplot2::aes(x=Dim_1, y=Dim_2, color = norm_cell_prob[cells, i])) +
                       ggplot2::geom_point(size = pt.size) +
                       ggplot2::scale_color_gradient2(low = "#011F5B", high = "#990000", mid = "light gray", midpoint = stats::median(as.vector(norm_cell_prob))) +
                       ggplot2::labs(color = "Cell Prob", title = ifelse(is.null(title_names), "", title_names[as.character(i)])) +
                       ggplot2::theme_classic())
    }
  }

  invisible(ConDecon_obj)
}


#' Calculate Relative cell prob distribution
#' @description Calculate the cell probability z-score for each cell and bulk sample. Cell probability z-score associated with query sample_i:
#' z-score_i = (ConDecon_obj$Normalized_cell.prob[,i] - rowMeans(ConDecon_obj$Normalized_cell.prob)) / RowSds(ConDecon_obj$Normalized_cell.prob)
#' @importFrom grDevices boxplot.stats
#' @importFrom stats sd
#' @importFrom Matrix rowMeans
#' @importFrom matrixStats rowSds
#' @param ConDecon_obj ConDecon object (output from RunConDecon)
#' @param rm_outliers Boolean indicating whether to remove outlier cell probability values from the z-score calculation based on +- 1.5*IQR (default = FALSE)
#'
#' @return Relative_cell.prob matrix within ConDecon object containing the z-score cell probability values for each query bulk sample
#' @export
#'
#' @examples
#' data(counts_gps)
#' data(latent_gps)
#' data(bulk_gps)
#' data(variable_genes_gps)
#'
#' # For this example, we will reduce the training size to max.iter = 50 to reduce run time
#' ConDecon_obj = RunConDecon(counts = counts_gps, latent = latent_gps, bulk = bulk_gps,
#' variable.features = variable_genes_gps, max.iter = 50)
#'
#' # Recalc relative cell probability
#' ConDecon_obj = CalcRelativeCellProb(ConDecon_obj)
CalcRelativeCellProb <- function(ConDecon_obj,
                                 rm_outliers = FALSE){

  #Remove outliers defined as +-1.5*IQR
  if(rm_outliers){
    avg_prob <- apply(ConDecon_obj$Normalized_cell.prob, 1, function(cell_prob){
      mean(cell_prob[!(cell_prob %in% grDevices::boxplot.stats(cell_prob)$out)])
    })
    sd_prob <- apply(ConDecon_obj$Normalized_cell.prob, 1, function(cell_prob){
      stats::sd(cell_prob[!(cell_prob %in% grDevices::boxplot.stats(cell_prob)$out)])
    })
    #Ignore outliers
  } else{
    avg_prob <- Matrix::rowMeans(ConDecon_obj$Normalized_cell.prob)
    sd_prob <- matrixStats::rowSds(ConDecon_obj$Normalized_cell.prob)
  }

  #Calculate relative values
  ConDecon_obj$Relative_cell.prob <- (ConDecon_obj$Normalized_cell.prob - avg_prob)/sd_prob

  #Order list of lists in ConDecon object
  ConDecon_obj <- ConDecon_obj[c("Normalized_cell.prob", "Relative_cell.prob", "PredictCellProb", "Model", "TrainingSet")]
  class(ConDecon_obj) <- "ConDecon"
  return(ConDecon_obj)
}

