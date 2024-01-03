#' ComputeProximityScore
#' @description Compute the average distance between the query dataset and the nearest 10 training datasets normalized by the average distance between training datasets
#' @import pdist
#' @import Matrix
#' @param cond ConDecon object
#'
#' @return ConDecon object with computed proximity score
#' @export
ComputeProximityScore <- function(cond) {
  euclid_query_train <- as.matrix(pdist((cond$PredictCellProb$bulk_coefficients), t(cond$TrainingSet$bulk_coefficients)))
  k <- 10
  knn <- t(apply(euclid_query_train, 1, function(i){
    order(i, decreasing = F)[1:k]
  }))
  knn <- cbind(1:nrow(knn), knn)
  knn_dist <- t(apply(knn, 1, function(i){
    euclid_query_train[i[1],i[2:(k+1)]]
  }))
  
  euclid_train <- as.matrix(pdist(t(cond$TrainingSet$bulk_coefficients), t(cond$TrainingSet$bulk_coefficients)))
  k <- 11
  knn_train <- t(apply(euclid_train, 1, function(i){
    order(i, decreasing = F)[1:k]
  }))
  knn_dist_train <- t(apply(knn_train, 1, function(i){
    euclid_train[i[1],i[2:(k)]]
  }))
  
  cond$proximity_score <- rowMeans(knn_dist/mean(as.vector(knn_dist_train)))
  
  return(cond)
}
