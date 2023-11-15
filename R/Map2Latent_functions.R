#' Map2Latent
#' @description Use the latent space to expand the distributions
#' of cell abundances and correlations on a basis of functions
#' @importFrom stats lm cor
#' @param TrainingSet Training set data (output from BuildTrainingSet)
#' @param latent matrix of single-cell latent space (cells x dims)
#' @param count single-cell count matrix (features x cells)
#' @param bulk matrix of query bulk data (features x samples)
#' @param variable.features character vector of variable features
#' @param k Number of nearest neighbor cells aggregated together when calculating rank correlation
#'
#' @return ConDecon object with low dimensional embedding of the space
#' of cell abundances and correlations
#' @export
#'
#' @examples
#' data(counts_gps)
#' data(latent_gps)
#' data(bulk_gps)
#' data(variable_genes_gps)
#'
#' # For this example, we will reduce the training size to max.iter = 50 to reduce run time
#' TrainingSet = BuildTrainingSet(count = counts_gps, latent = latent_gps, max.iter = 50)
#'
#' ConDecon_obj = Map2Latent(TrainingSet = TrainingSet, latent = latent_gps, count = counts_gps,
#' bulk = bulk_gps, variable.features = variable_genes_gps)
Map2Latent <- function(TrainingSet,
                       latent,
                       count,
                       bulk,
                       variable.features,
                       k = 1){

  #Variable features must be located in bulk AND single-cell data
  counts_genes <- unique(match(row.names(bulk), row.names(count)))
  count <- count[counts_genes[!is.na(counts_genes)],,drop=F]
  bulk_genes <- unique(match(row.names(count), row.names(bulk)))
  bulk <- bulk[bulk_genes[!is.na(bulk_genes)],,drop=F]
  bulk <- bulk[row.names(count),,drop=F]

  which_features <- unique(c(match(row.names(count), variable.features),
                             match(row.names(bulk), variable.features)))
  which_features <- which_features[!is.na(which_features)]
  variable.features <- variable.features[which_features]


  TrainingSet$TrainingSet$latent <- as.matrix(latent[,1:TrainingSet$TrainingSet$dims])

  #Coefficients for cell prob (no intercept)
  TrainingSet$TrainingSet$cell.prob_coefficients <-
    stats::lm(TrainingSet$TrainingSet$cell.prob ~ TrainingSet$TrainingSet$latent + 0)$coefficients

  #Find index of KNN for dist matrix
  knn <- as.matrix(apply(TrainingSet$TrainingSet$latent_distance, 1, function(i){
    order(i, decreasing = F)[1:k]
  }))
  
  if(k > 1){
    knn <- t(knn)
  }
  TrainingSet$TrainingSet$knn <- knn
  #Find the index of the variable features in each cell when ordered by expression
  #(averaged over its 5 nearest neighbors)
  gene.index <- GeneIndex(knn, count, variable.features, k = k)

  #Correlation btwn synthetic bulk and single-cell data
  TrainingSet$TrainingSet$bulk_nn <-
    CorBulk(TrainingSet$TrainingSet$synthetic_bulk, variable.features, gene.index)

  #Coefficients for synthetic bulk
  TrainingSet$TrainingSet$bulk_coefficients <-
    stats::lm(TrainingSet$TrainingSet$bulk_nn ~ TrainingSet$TrainingSet$latent + 0)$coefficients

  return(TrainingSet)
}

#' GeneIndex
#' @import Matrix
#' @importFrom tidyr gather
#' @param knn Matrix of cells with their nearest neighbor
#' @param count Matrix of single-cell count data
#' @param variable.features Character vector of variable features from single-cell data
#' @param k Number of nearest neighbor cells aggregated together when calculating rank correlation
#'
#' @return Matrix of rank indices for variable features ordered by expression
GeneIndex <- function(knn,count,variable.features, k){

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
  count_nn <- (norm_count %*% knn.index)/k
  #genes_ID is the index of variable genes within the rownames of count_nn
  genes_id <- match(variable.features, row.names(count_nn))
  gene.index <- apply(count_nn, 2, GeneIndex_subfxn, geneID=genes_id)
  return(gene.index)
}

#' GeneIndex_subfxn
#' @param i Bulk sample
#' @param geneID Vector with index of variable features
#'
#' @return Numeric vector of gene rank
GeneIndex_subfxn <- function(i,geneID){
  #For each cell, order the genes by rank expression (from high to low)
  i.order <- rank(-i)
  return(i.order[geneID])
}

#' CorBulk
#' @param bulk Matrix of query bulk data (features x samples)
#' @param genes Character vector of variable features from single-cell data
#' @param gene.index Numeric vector of gene rank
#'
#' @return Matrix of rank correlation for each bulk sample
CorBulk <- function(bulk, genes, gene.index){
  genes_id <- match(genes,row.names(bulk))
  #Finds the rank of the differentially expressed genes when each bulk sample is ordered decreasingly by expression
  bulk.gene.index <- apply(bulk, 2, GeneIndex_subfxn, geneID=genes_id)
  cor.bulk <- stats::cor((gene.index), (bulk.gene.index), method = "pearson")
  return(cor.bulk)
}

