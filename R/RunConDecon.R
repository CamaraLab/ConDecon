#' Run ConDecon
#' @description Run ConDecon (continuous deconvolution) to estimate the cell abundances associated with each query bulk sample
#' @importFrom methods canCoerce is
#' @import Matrix
#' @param counts single-cell counts matrix (features x cells)
#' @param latent matrix of single-cell latent space (cells x dims)
#' @param bulk matrix of query bulk data (features x samples)
#' @param variable.features character vector of the most variable features
#' @param max.iter size of the training dataset (default = 5,000)
#' @param max.cent max number of centers in the Gaussian (default = 5)
#' @param dims number of dimensions from latent (default = 10)
#' @param degree degree of the polynomial used
#' @param step manually parallelize building the training dataset
#' @param min.cent min number of centers in the Gaussian (default = 1)
#' @param n number of cells to be chosen to create the training dataset (default
#' is half the number of cells in the counts matrix)
#' @param trainingset pre-generated training dataset
#' @param sigma_min_cells min number of cells that should be captured by the standard
#' deviation of the Gaussian
#' @param sigma_max_cells max number of cells that should be captured by the standard
#' deviation of the Gaussian
#' @param verbose logical indicating whether to print progress (default = FALSE)
#'
#' @return ConDecon object with continuous deconvolution results
#' @export
#'
#' @examples
#' data(counts_gps)
#' data(latent_gps)
#' data(bulk_gps)
#' data(variable_genes_gps)
#'
#' # For this example, we will reduce the training size to max.iter = 50 to reduce run time
#' ConDecon_obj = RunConDecon(counts = counts_gps, latent = latent_gps,
#' bulk = bulk_gps, variable.features = variable_genes_gps, max.iter = 50)
RunConDecon <- function(counts,
                      latent,
                      bulk,
                      variable.features,
                      max.iter = 5000,
                      max.cent = 5,
                      dims = 10,
                      degree = 1,
                      step = ifelse(max.iter <= 10000, max.iter, 10000),
                      min.cent= 1,
                      n = round(ncol(counts)/2),
                      trainingset = NULL,
                      sigma_min_cells = NULL,
                      sigma_max_cells = NULL,
                      verbose = FALSE){


  ### CHECK INPUT ###

  if(is.null(step)){
    if(max.iter <= 10000){
      step <- max.iter
    } else{
      step <- 10000
    }
  }

  if(is.null(counts)){
    message("counts is required; see documentation")
    return(NULL)
  } else if(!methods::is(counts,"dgCMatrix")){
    if(methods::canCoerce(counts, "matrix")){
      counts <- Matrix(as.matrix(counts),sparse=TRUE)
    } else if(methods::canCoerce(counts, "data.frame")){
      counts <- Matrix(as.matrix(as.data.frame(counts)),sparse=TRUE)
    } else {
      message("counts must be a matrix")
      return(NULL)
    }
  }
  if (sum(Matrix::rowSums(counts) < 0) >0){
    message(paste("entries of counts must be non-negative"))
    return(NULL)
  } else if(sum(Matrix::rowSums(counts) == 0) > 0){
    counts <- counts[Matrix::rowSums(counts)>0,,drop=F]
    if(verbose == TRUE){
      cat("genes without expression will be removed from the counts data\n")
    }
  }
  counts_genes <- unique(match(row.names(bulk), row.names(counts)))
  counts <- counts[counts_genes[!is.na(counts_genes)],,drop=F]
  bulk_genes <- unique(match(row.names(counts), row.names(bulk)))
  bulk <- bulk[bulk_genes[!is.na(bulk_genes)],,drop=F]
  bulk <- bulk[row.names(counts),,drop=F]

  if (n < 1){
    message("n must be numeric value greater than one")
    return(NULL)
  } else if(n != round(n)){
    if(verbose == TRUE){
      cat("Rounding n to a whole number")
    }
    n <- round(n)
  }

  if(is.null(latent)){
    message("latent matrix is required; latent can be obtained by dimensionality reduction")
    return(NULL)
  } else if (methods::canCoerce(latent, "data.frame")){
    latent <- as.data.frame(latent)
  } else if (methods::canCoerce(latent, "matrix")){
    latent <- as.data.frame(as.matrix(latent))
  } else {
    message("latent must be a matrix or data.frame")
    return(NULL)
  }
  if (ncol(latent) < dims){
    if(verbose == TRUE){
      message(paste0("There are ", ncol(latent), " dimensions in the latent space, which is
                   less than the parameter dims\nOverriding the dims parameter such that
                   dims =  ", ncol(latent)))
    }
    dims = ncol(latent)
  }

  tot.cells <- ncol(counts)
  if(nrow(latent) != tot.cells | ncol(latent) <2){
    message("number of rows in latent must be equal to number of cells (columns) of counts;
            \nnumber of columns in latent must be greater than or equal to two")
    return(NULL)
  }

  if(degree > dims){
    message("degree must be less than or equal to dims")
    return(NULL)
  }

  if(is.vector(variable.features)){
    #Variable features must be located in bulk and single-cell data
    which_features <- unique(c(match(row.names(counts), variable.features), match(row.names(bulk),
                                                                                  variable.features)))
    which_features <- which_features[!is.na(which_features)]
    variable.features <- variable.features[which_features]
    if(verbose == TRUE){
      cat(paste0("There are ", length(variable.features), " variable.features in both the counts and bulk data\n"))
    }
  } else{
    message("variable.features must be a vector of features")
    return(NULL)
  }

  output <- NULL
  output <- vector(mode="list", length=0)

  #You can provide a trainingset object
  if (is.null(trainingset)){
    if(verbose == TRUE){
      message("Building training set")
    }
    output <- BuildTrainingSet(counts, latent, max.iter, max.cent, step, dims,
                                           min.cent, n, sigma_min_cells = sigma_min_cells,
                                           sigma_max_cells = sigma_max_cells, verbose = verbose)

  } else {
    if(!("TrainingSet" %in% names(trainingset))){
      message("training set does not exist")
      return(NULL)
    }
    if(trainingset$TrainingSet$dims != dims){
      warning("The number of latent dimensions used to create TrainingSet (TrainingSet$dims) is different from input dims\n")
    }
    output <- trainingset
  }

  if(verbose == TRUE){
    message("Map2Latent")
  }
  output <- Map2Latent(output, latent, counts, bulk, variable.features)

  if(verbose == TRUE){
    message("BuildModel")
  }
  output <- BuildModel(output, degree)

  if(verbose == TRUE){
    message("PredictCellProb")
  }
  output <- PredictCellProb(bulk, counts, variable.features, output)

  #Calc Relative cell probability
  output <- CalcRelativeCellProb(output)

  #Order list of lists in ConDecon object
  output <- output[c("Normalized_cell.prob", "Relative_cell.prob", "PredictCellProb", "Model", "TrainingSet")]
  class(output) = "ConDecon"
  return(output)
}
