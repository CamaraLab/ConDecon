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
#' # For this example, we will reduce the training size to max.iter = 50 to reduce run time
#' TrainingSet = BuildTrainingSet(count = counts_gps, latent = latent_gps, max.iter = 50)
#' ConDecon_obj = Map2Latent(TrainingSet = TrainingSet, latent = latent_gps, count = counts_gps,
#' bulk = bulk_gps, variable.features = variable_genes_gps)
#' ConDecon_obj = BuildModel(ConDecon_obj)
#'
#' ConDecon_obj = PredictCellProb(bulk = bulk_gps, count = counts_gps,
#' variable.features = variable_genes_gps, output = ConDecon_obj)
PredictCellProb <- function(bulk,
                            count,
                            variable.features,
                            output){

  #Variable features must be located in bulk AND single-cell data
  counts_genes <- unique(match(row.names(bulk), row.names(count)))
  count <- count[counts_genes[!is.na(counts_genes)],,drop=F]
  bulk_genes <- unique(match(row.names(count), row.names(bulk)))
  bulk <- bulk[bulk_genes[!is.na(bulk_genes)],,drop=F]
  bulk <- bulk[row.names(count),,drop=F]

  which_features <- unique(c(match(row.names(count), variable.features), match(row.names(bulk),
                                                                               variable.features)))
  which_features <- which_features[!is.na(which_features)]
  variable.features <- variable.features[which_features]

  #Find k=5 nearest neighbors of each cell
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
  output$PredictCellProb$bulk_coefficients <- data.frame(t(bulk.coef))

  #Infer cell prob coefficients
  output$PredictCellProb$cell.prob_coefficients <- NULL
  for (j in 1:length(output$Model)){
    tmp <- stats::predict(output$Model[[j]], output$PredictCellProb$bulk_coefficients)
    output$PredictCellProb$cell.prob_coefficients <- rbind(output$PredictCellProb$cell.prob_coefficients, tmp)
  }

  #Infer cell prob
  output$PredictCellProb$infer_cell.prob <- NULL
  output$PredictCellProb$infer_cell.prob = output$TrainingSet$latent %*% output$PredictCellProb$cell.prob_coefficients

  #normalize cell prob st cell prob sum to 1
  colmin <- apply(output$PredictCellProb$infer_cell.prob, 2, min)
  cellprob_scale <- t(t(output$PredictCellProb$infer_cell.prob)+abs(colmin))
  output$Normalized_cell.prob = t(t(cellprob_scale)/colSums(cellprob_scale))

  return(output)
}
