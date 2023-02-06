#' TransferFeatures
#' @description Using ConDecon's predicted cell probability distributions, transfer a numeric feature from the single-cell data to the bulk data
#' @importFrom Matrix colSums
#' @param ConDecon_obj ConDecon object (output from RunConDecon)
#' @param feature Numeric vector with support on the single-cell latent space (eg pseudotime, gene expression, etc)
#' @param feature_name String indicating the name where the transferred feature will be stored within the ConDecon object (default = the object name input into feature)
#' @param probs Numeric value indicating the bottom quartile of the distribution that should be removed for the purpose of smoothing long tails of the predicted cell probability distribuiton (default = 0.75)
#'
#' @return ConDecon object with a matrix called 'TransferFeatures' containing the transferred feature (rows) for each bulk sample (column)
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
#' # Transfer feature to ConDecon object: ConDecon_obj$TransferFeatures[feature_name,]
#' # For this example, randomly selected gene from the count matrix to transfer
#' random_gene = counts_gps[sample(x = 1:nrow(counts_gps), size = 1),]
#' ConDecon_obj = TransferFeatures(ConDecon_obj = ConDecon_obj, feature = random_gene)
TransferFeatures <- function(ConDecon_obj,
                             feature,
                             feature_name = deparse(substitute(feature)),
                             probs = 0.75){

  #Use invisible() to capture user's original name for feature
  invisible(paste0("Transferring ", feature_name, "... \n"))

  if(!inherits(ConDecon_obj, what = "ConDecon")){
    message("ConDecon_obj must be an object with class 'ConDecon' output from RunConDecon")
    return(ConDecon_obj)
  }

  if(is.null(ConDecon_obj$TransferFeatures)){
    ConDecon_obj$TransferFeatures <- NULL
  }

  if(!is.vector(feature)){
    message("feature must be a vector")
    return(ConDecon_obj)
  }
  if(length(feature) != nrow(ConDecon_obj$Normalized_cell.prob)){
    message("The length of 'feature' must be equal to the number of cells in the\n
            reference single-cell data")
    return(ConDecon_obj)
  }
  if(is.null(names(feature))){
    names(feature) <- rownames(ConDecon_obj$Normalized_cell.prob)
  } else if(sum(names(feature) %in% row.names(ConDecon_obj$Normalized_cell.prob)) != nrow(ConDecon_obj$Normalized_cell.prob)){
    message("Names of 'feature' must be the column/cell names in the reference single-cell data")
    return(ConDecon_obj)
  } else {
    feature = feature[row.names(ConDecon_obj$Normalized_cell.prob)]
  }

  ## If feature is calculated for a subset of the cells, remove cell with NA
  if(sum(is.na(feature)) == length(feature)){
    message("'feature' only contains NA values")
    return(ConDecon_obj)
  } else if(sum(is.na(feature))>0){
    Normalized_cell.prob <- ConDecon_obj$Normalized_cell.prob[!is.na(feature),]
    feature <- feature[!is.na(feature)]
  } else{
    Normalized_cell.prob <- ConDecon_obj$Normalized_cell.prob
  }

  cat(paste0("Transferring ", feature_name, "... \n"))
  ## Remove long tails (bottom quartile of distribution) from predicted cell probability distributions
  cell_prob_smooth = smooth_cell_prob(Normalized_cell.prob, probs = probs)
  ## Dot product between smoothed cell probability distribution and numeric feature
  transfer_feature = feature %*% cell_prob_smooth
  ## Normalize
  transfer_feature = transfer_feature/Matrix::colSums(cell_prob_smooth)
  row.names(transfer_feature) <- feature_name

  if(feature_name %in% row.names(ConDecon_obj$TransferFeatures)){
    ConDecon_obj$TransferFeatures[feature_name,] <- transfer_feature
  } else{
    ConDecon_obj$TransferFeatures <- rbind(ConDecon_obj$TransferFeatures, transfer_feature)
  }

  #Order list of lists in ConDecon object
  ConDecon_obj <- ConDecon_obj[c("Normalized_cell.prob", "Relative_cell.prob", "TransferFeatures", "PredictCellProb", "Model", "TrainingSet")]
  class(ConDecon_obj) = "ConDecon"
  return(ConDecon_obj)
}

#' Set values to zero
#' @param v Vector of predicted cell probability values for a single sample
#' @param threshold Numeric of the distribution that represents the determined quartile value
#'
#' @return Vector of cell probability values truncated based on threshold
set_zero <- function(v, threshold){
  v[v<threshold] <- 0
  return(v)
}

#' Smooth ConDecon's cell.prob matrix
#' @importFrom stats quantile
#' @param cell_prob Matrix of ConDecon's predicted cell probability distributions
#' @param probs Numeric value indicating the bottom quartile of the distribution that should be removed for the purpose of smoothing long tails of the predicted cell probability distribuiton (default = 0.75)
#'
#' @return Matrix with truncated cell probability values based on quartile value
smooth_cell_prob = function(cell_prob, probs = 0.75){

  threshold_samples = apply(cell_prob, 2, function(i){
    stats::quantile(i, probs = probs)
  })

  cell_prob_smooth <- NULL
  for(i in 1:ncol(cell_prob)){
    cell_prob_smooth <- cbind(cell_prob_smooth, set_zero(cell_prob[,i], threshold_samples[i]))
  }
  colnames(cell_prob_smooth) <- colnames(cell_prob)
  return(cell_prob_smooth)
}

