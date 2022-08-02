#' Plot ConDecon_obj
#' @description Calculate and plot the cell probability z-score for each cell and bulk sample. Cell probability z-score associated with query sample_i:
#' z-score_i = [ConDecon_obj$norm_cell.prob[,i] - rowMeans(ConDecon_obj$norm_cell.prob)] / RowSds(ConDecon_obj$norm_cell.prob)
#' @importFrom matrixStats rowSds
#' @importFrom grDevices boxplot.stats
#' @importFrom graphics plot
#' @importFrom stats sd
#' @import ggplot2
#' @import Matrix
#' @param ConDecon_obj ConDecon object (output from RunConDecon)
#' @param umap 2-dimensional coordinates of reference single-cell data (default = first 2 dimensions of latent space)
#' @param samples character vector of the name(s) of the query bulk sample(s) (default = NULL)
#' @param pt.size size of the points plotted (default = 1)
#' @param rm_outliers Boolean indicating whether to remove outlier cell probability values from the z-score calculation based on +- 1.5*IQR (default = FALSE)
#' @param cells vector of cell names/index indicating which cells from the reference single-cell data to plot (default = NULL, plot all cells)
#' @param title_names title of the plot (default = NULL, use bulk sample names from ConDecon object)
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
#' ConDecon_obj = RunConDecon(counts = counts_gps, latent = latent_gps, bulk = bulk_gps,
#' variable.features = variable_genes_gps)
#'
#' # add z-score cell abundances to ConDecon object: ConDecon$relative_cell.prob
#' ConDecon_obj = PlotConDecon(ConDecon_obj)
#'
#' # Plot the z-score cell abundances for the query bulk sample(s)
#' PlotConDecon(ConDecon_obj)
PlotConDecon <- function(ConDecon_obj,
                             umap = ConDecon_obj$TrainingSet$latent[,1:2],
                             samples = NULL,
                             pt.size = 1,
                             rm_outliers = FALSE,
                             cells = NULL,
                             title_names = NULL){

  umap <- data.frame(UMAP_1 = umap[,1], UMAP_2 = umap[,2])

  ##Which bulk samples to plot
  if(is.null(samples)){
    samples <- 1:ncol(ConDecon_obj$norm_cell.prob)
    if(is.null(title_names)){
      title_names <- colnames(ConDecon_obj$norm_cell.prob)
    }
  } else{
    samples <- as.character(samples)
    title_names <- samples
    names(title_names) <- samples
  }

  ##Which cells to plot
  if(is.null(cells)){
    norm_cell.prob <- ConDecon_obj$norm_cell.prob
  } else if (length(cells) > ncol(ConDecon_obj$norm_cell.prob)){
    message("The length of parameter cells is greater than the number of cells in
            ConDecon_obj$norm_cell.prob")
    return(NULL)
  }else{
    norm_cell.prob <- ConDecon_obj$norm_cell.prob[cells,]
    umap <- umap[cells,]
  }

  ##Remove outliers defined as +-1.5*IQR
  if(rm_outliers){
    avg_prob <- apply(norm_cell.prob, 1, function(cell_prob){
      mean(cell_prob[!(cell_prob %in% grDevices::boxplot.stats(cell_prob)$out)])
    })
    sd_prob <- apply(norm_cell.prob, 1, function(cell_prob){
      stats::sd(cell_prob[!(cell_prob %in% grDevices::boxplot.stats(cell_prob)$out)])
    })
  } else{
    avg_prob <- rowMeans(norm_cell.prob)
    sd_prob <- matrixStats::rowSds(norm_cell.prob)
  }

  ConDecon_obj$relative_prob <- NULL
  for (i in samples){
    subtract_avg <- (norm_cell.prob[,i]- avg_prob)/sd_prob
    ConDecon_obj$relative_prob <- cbind(ConDecon_obj$relative_prob, subtract_avg)

    graphics::plot(ggplot(umap, aes(x=UMAP_1, y=UMAP_2, color = subtract_avg)) +
           geom_point(size = pt.size) + scale_color_gradient2(low = "blue", high = "red",
                                                              mid = "gray") +
           labs(color = "Cell Prob\nz-score", title = title_names[i]) + theme_classic())
  }
  invisible(ConDecon_obj)
}
