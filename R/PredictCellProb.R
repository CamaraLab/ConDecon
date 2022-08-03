#' PredictCellProb
#' @description Using the fit model, predict the cell abundance distributions
#' @importFrom stats lm predict
#' @import Matrix
#' @param bulk matrix of the query bulk data (features x samples)
#' @param count single-cell count matrix (features x cells)
#' @param variable.features character vector of the most variable features
#' @param output ConDecon object with fit model
#'
#' @return ConDecon object with inferred cell probabilities for each query bulk sample
#' @export
#'
#' @examples
#' data(counts_gps)
#' data(variable_genes_gps)
#' data(bulk_gps)
#'
#' TrainingSet = BuildTrainingSet(count = counts_gps, latent = latent_gps)
#' ConDecon_obj = Map2Latent(output = TrainingSet, latent = latent_gps, count = counts_gps,
#' variable.features = variable_genes_gps)
#' ConDecon_obj = BuildModel(ConDecon_obj)
#'
#' ConDecon_obj = PredictCellProb(bulk = bulk_gps, count = counts_gps,
#' variable.features = variable_genes_gps, output = ConDecon_obj)
PredictCellProb <- function(bulk,
                            count,
                            variable.features,
                            output){

  k <- 5
  knn <- t(apply(output$TrainingSet$latent_distance, 1, function(i){
    order(i, decreasing = F)[1:k]
  }))

  #Infer bulk coefficients
  gene.index <- GeneIndex(knn, count, variable.features)
  bulk_nn <- CorBulk(bulk, variable.features, gene.index)
  bulk.coef <- stats::lm(bulk_nn ~ output$TrainingSet$latent + 0)$coefficients
  bulk.coef <- as.matrix(bulk.coef)

  #predict() needs at least 2 cases to predict
  if (ncol(bulk.coef) == 1){
    bulk.coef <- cbind(bulk.coef, 0)
  }

  row.names(bulk.coef) <- paste0("x",1:ncol(output$TrainingSet$latent))
  output$bulk_coefficients <- data.frame(t(bulk.coef))

  #Infer cell prob coefficients
  output$cell.prob_coefficients <- NULL
  for (j in 1:length(output$Model)){
    tmp <- stats::predict(output$Model[[j]], output$bulk_coefficients)
    output$cell.prob_coefficients <- rbind(output$cell.prob_coefficients, tmp)
  }

  #Infer cell prob
  output$infer_cell.prob <- NULL
  output$infer_cell.prob = output$TrainingSet$latent %*% output$cell.prob_coefficients

  #normalize cell prob st cell prob sum to 1
  colmin <- apply(output$infer_cell.prob, 2, min)
  cellprob_scale <- t(t(output$infer_cell.prob)+abs(colmin))
  output$norm_cell.prob = t(t(cellprob_scale)/colSums(cellprob_scale))

  return(output)
}
