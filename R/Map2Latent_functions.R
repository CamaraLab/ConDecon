#' Map2Latent
#' @description Use the latent space to expand the distributions
#' of cell abundances and correlations on a basis of functions
#' @importFrom stats lm cor
#' @param output Training set object
#' @param latent matrix of single-cell latent space (cells x dims)
#' @param count single-cell count matrix (features x cells)
#' @param variable.features character vector of variable features
#'
#' @return ConDecon object with low dimensional embedding of the space
#' of cell abundances and correlations
#' @export
#'
#' @examples
#' data(counts_gps)
#' data(latent_gps)
#' data(variable_genes_gps)
#'
#' TrainingSet = BuildTrainingSet(count = counts_gps, latent = latent_gps)
#'
#' ConDecon_obj = Map2Latent(output = TrainingSet, latent = latent_gps, count = counts_gps,
#' variable.features = variable_genes_gps)
Map2Latent <- function(output,
                       latent,
                       count,
                       variable.features){


  output$TrainingSet$latent <- as.matrix(latent[,1:output$TrainingSet$dims])

  #Coefficients for cell prob (no intercept)
  output$TrainingSet$cell.prob_coefficients <- stats::lm(output$TrainingSet$cell.prob ~ output$TrainingSet$latent + 0)$coefficients

  #Find index of KNN for dist matrix
  k <- 5
  knn <- t(apply(output$TrainingSet$latent_distance, 1, function(i){
    order(i, decreasing = F)[1:k]
  }))

  #Find the index of the variable features in each cell when ordered by expression
  #(averaged over its 5 nearest neighbors)
  gene.index <- GeneIndex(knn, count, variable.features)

  #Correlation btwn synthetic bulk and single-cell data
  output$TrainingSet$bulk_nn <- CorBulk(output$TrainingSet$synthetic_bulk, variable.features, gene.index)

  #Coefficients for synthetic bulk
  output$TrainingSet$bulk_coefficients <- stats::lm(output$TrainingSet$bulk_nn ~ output$TrainingSet$latent + 0)$coefficients

  return(output)
}

#' GeneIndex
#' @import Matrix
#' @importFrom tidyr gather
#' @param knn
#' @param count
#' @param variable.features
#'
#' @return
#'
#' @examples
GeneIndex <- function(knn,count,variable.features){

  knn.index <- tidyr::gather(as.data.frame(t(knn)))
  knn.index$new_key <- match(knn.index$key,row.names(knn))
  #colSums of knn.index will sum to k (the number of nearest neighbors per cell)
  #cells by cells matrix
  knn.index <- Matrix::sparseMatrix(i = knn.index$value, j = knn.index$new_key, dims = c(dim(knn)[1],
                                                                                         dim(knn)[1]))
  #normalize the count data before summing counts across cells
  norm_count <- t(t(count)/colSums(count))*10000
  #transform knn.index such that each column is the average expression of the the cell's nearest neighbors
  #genes by cell matrix
  count_nn <- (norm_count %*% knn.index)/5
  #genes_ID is the index of variable genes within the rownames of count_nn
  genes_id <- match(variable.features, row.names(count_nn))
  gene.index <- apply(count_nn, 2, GeneIndex_subfxn, geneID=genes_id)
  return(gene.index)
}

#' GeneIndex_subfxn
#'
#' @param i
#' @param geneID
#'
#' @return
#'
#' @examples
GeneIndex_subfxn <- function(i,geneID){
  #For each cell, order the genes by rank expression (from high to low)
  i.order <- rank(-i)
  return(i.order[geneID])
}

#' CorBulk
#'
#' @param bulk
#' @param genes
#' @param gene.index
#'
#' @return
#'
#' @examples
CorBulk <- function(bulk, genes, gene.index){
  genes_id <- match(genes,row.names(bulk))
  #Finds the rank of the differentially expressed genes when each bulk sample is ordered decreasingly by expression
  bulk.gene.index <- apply(bulk, 2, GeneIndex_subfxn, geneID=genes_id)
  cor.bulk <- stats::cor((gene.index), (bulk.gene.index), method = "pearson")
  return(cor.bulk)
}

