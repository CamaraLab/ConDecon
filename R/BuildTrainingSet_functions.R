#' Build the training dataset
#' @description Use a Gaussian mixture model (with random parameters)
#' to generate a traning dataset from the reference single-cell data
#' @importFrom stats dist uniroot
#' @import Matrix
#' @param count single-cell count matrix (features x cells)
#' @param latent matrix of single-cell latent space (cells x dims)
#' @param max.iter size of the training dataset (default = 10,000)
#' @param max.cent max number of centers in the Gaussian (default = 5)
#' @param step manually parallelize building the training dataset
#' @param dims number of dimensions from latent (default = ncol(latent))
#' @param min.cent min number of centers in the Gaussian (default = 1)
#' @param n number of cells to be chosen to create the training dataset (default
#' is half the number of cells in the count matrix)
#' @param sigma_min_cells min number of cells that should be captured by the standard
#' deviation of the Gaussian
#' @param sigma_max_cells max number of cells that should be captured by the standard
#' deviation of the Gaussian
#' @param verbose logical indicating whether to print progress (default = TRUE)
#'
#' @return ConDecon object with training data
#' @export
#'
#' @examples
#' data(counts_gps)
#' data(latent_gps)
#'
#' TrainingSet = BuildTrainingSet(count = counts_gps, latent = latent_gps)
BuildTrainingSet <- function(count,
                             latent,
                             max.iter = 10000,
                             max.cent = 5,
                             step = 10000,
                             dims = ncol(latent),
                             min.cent = 1,
                             n = round(ncol(count)/2),
                             sigma_min_cells = NULL,
                             sigma_max_cells = NULL,
                             verbose = FALSE){

  output <- NULL
  output <- vector(mode="list", length=0)

  tot.cells <- ncol(count)

  output$latent_distance <- as.matrix(stats::dist(latent[,1:dims,drop=F]))
  output$dims <- dims

  N <- ncol(count)



  ### FIND SIGMA VALUES ###

  # Find sigma min
  if(is.null(sigma_min_cells)){
    N_min <- max(20, round(N/100))
  } else{
    N_min <- round(sigma_min_cells)
  }

  # Find sigma max
  if(is.null(sigma_max_cells)){
    N_max <- ifelse(round(N/5)<N_min, N, round(N/5))
  } else{
    N_max <- round(sigma_max_cells)
  }
  k <- N_max
  knn <- t(apply(output$latent_distance, 1, function(i){
    order(i, decreasing = F)[1:k]
  }))
  knn_dist <- t(apply(knn, 1, function(i){
    output$latent_distance[i[1],i]
  }))
  k_dist <- colMeans(knn_dist)
  error.sigma.max <- FALSE
  tryCatch( {sigma_max <- stats::uniroot(Find_Sigma, x_i=seq(0,max(k_dist), by = 0.01),
                                  width = 0.01, lower=min(knn_dist[knn_dist>0]),
                                  upper=max(knn_dist))$root }, error = function(e)
                                    {error.sigma.max <<- TRUE})
  # sigma_max <- uniroot(Find_Sigma, x_i=seq(0,max(k_dist), by = 0.01), width = 0.01, lower=min(knn_dist[knn_dist>0]), upper=max(knn_dist))$root

  # Find sigma min
  k <- N_min
  knn <- t(apply(output$latent_distance, 1, function(i){
    order(i, decreasing = F)[1:k]
  }))
  knn_dist <- t(apply(knn, 1, function(i){
    output$latent_distance[i[1],i]
  }))
  k_dist <- colMeans(knn_dist)
  # sigma_min <- uniroot(Find_Sigma, x_i=seq(0,max(k_dist), by = 0.01), width = 0.01, lower=min(knn_dist[knn_dist>0]), upper=max(knn_dist))$root
  error.sigma.min <- FALSE
  tryCatch( { sigma_min <- uniroot(Find_Sigma, x_i=seq(0,max(k_dist), by = 0.01),
                                   width = 0.01, lower=min(knn_dist[knn_dist>0]),
                                   upper=max(knn_dist))$root }, error = function(e)
                                     {error.sigma.min <<- TRUE})
  if(error.sigma.max & error.sigma.min){
    stop("Could not solve for sigma\nTry increasing sigma_max_cells and/or sigma_min_cells\nor decreasing the number of dimensions in the latent space to decrease sparsity")
    return(NULL)
  } else if(error.sigma.min & error.sigma.max == FALSE){
    sigma_min = sigma_max/5
  } else if(error.sigma.max & error.sigma.min == FALSE){
    sigma_max = sigma_min*5
  }
  # if(verbose == TRUE){
  #   message("Sigma found")
  # }



  ### SELECT TRAINING PARAMETERS ###

  output$parameters <- matrix(0, ncol = (1+max.cent*3), nrow = max.iter, dimnames =
                                list(NULL,c("num.centers",paste0("center.",1:max.cent),
                                            paste0("sigma.",1:max.cent),paste0("mix.",1:max.cent))))

  # Randomly choose the parameters of a Gaussian mixture model for each training set
  ## Choose number of centers
  output$parameters[,1] <- sample(min.cent:max.cent,max.iter, replace = TRUE)

  ## Choose centers, sigma, and constant
  ## Find the cell probabilities
  output$cell.prob <- matrix(0, ncol = max.iter, nrow = tot.cells)
  for(i in 1:max.iter){

    # Choose sigma
    output$parameters[i,(2+max.cent):(1+max.cent+output$parameters[i,1])] <-
      sample(10^seq(log10(sigma_min), log10(sigma_max), by = 0.0001), output$parameters[i,1],
             replace=TRUE)

    # Choose mixture between [1:100]
    mixture <- sample(1:100,output$parameters[i,1], replace=TRUE)
    # mixture <- rep(100/output$parameters[i,1], output$parameters[i,1])
    output$parameters[i,(2+max.cent*2):(1+max.cent*2+output$parameters[i,1])] <- mixture/sum(mixture)

    # Choose which cells will be the centers
    output$parameters[i,2:(1+output$parameters[i,1])] <-
      sample(1:tot.cells,output$parameters[i,1],replace = FALSE)

    # Cell Prob
    output$cell.prob[,i] <- CellProb(parameters = output$parameters[i,], latent_distance = output$latent_distance,
                                     max.cent = max.cent, tot.cells=tot.cells)
  }
  gc(verbose = FALSE)

  pick.cells <- PickCells(max.iter, step, output$cell.prob, max.cent, tot.cells, n, verbose)
  #rm(cell.prob)
  gc(verbose = FALSE)

  #Aggregate single-cell expression weighted by GMM to create synthetic bulk
  output$synthetic_bulk <- BulkCells(pick.cells, max.iter, step, max.cent, tot.cells, n, count, verbose)
  output$TrainingSet <- output
  return(output)
}

#' Find sigma
#'
#' @param sigma
#' @param x_i
#' @param width
#'
#' @return
#'
#' @examples
Find_Sigma <- function(sigma, x_i, width){
  sum((1/(sqrt(2*pi)*sigma))*exp(-(x_i)^2/(2*sigma^2)))*width-0.5
}

#' Calculate cell probability
#' @description  Calculate the probability of each cell based on a Gaussian distribution
#'
#' @param parameters
#' @param latent_distance
#' @param max.cent
#' @param tot.cells
#'
#' @return
#'
#' @examples
CellProb <- function(parameters,latent_distance,max.cent,tot.cells){
  # message("CellProb internal")
  num.cent <- parameters[1]
  g.mix.model <- t((parameters[(2+max.cent*2):(1+max.cent*2+num.cent)]/
                      (sqrt(2*pi)*parameters[(2+max.cent):(1+max.cent+num.cent)]))*
                     exp(-1*t(latent_distance[,parameters[2:(1+num.cent)]])^2/
                           (2*parameters[(2+max.cent):(1+max.cent+num.cent)]^2)))
  return(rowSums(g.mix.model)/sum(rowSums(g.mix.model)))
}

#' Pick cells based on the Gaussian kernal
#'
#' @param max.iter
#' @param step
#' @param cell.prob
#' @param max.cent
#' @param tot.cells
#' @param n
#' @param verbose
#'
#' @return
#'
#' @examples
PickCells <- function(max.iter,step,cell.prob,max.cent,tot.cells,n,verbose){
  # if(verbose == TRUE){
  #   message("PickCells")
  # }

  starting <- 1
  pick.cells <- NULL
  if ((max.iter/step) == round(max.iter/step)){
    for(i in seq(from = step, to = max.iter,by=step)){
      pick_cells <- lapply(starting:i,BulkCells_subfxn,cell.prob=cell.prob,
                           max.cent=max.cent, tot.cells=tot.cells,n=n)
      pick.cells <- rbind(pick.cells, do.call(rbind,pick_cells))
      starting <- i+1
    }
  } else if((max.iter/step) != round(max.iter/step)){
    for(i in c(seq(from = step, to = max.iter,by=step), max.iter)){
      pick_cells <- lapply(starting:i,BulkCells_subfxn,cell.prob=cell.prob,
                           max.cent=max.cent, tot.cells=tot.cells,n=n)
      pick.cells <- rbind(pick.cells, do.call(rbind, pick_cells))
      starting <- i+1
    }
  }
  return(pick.cells)
}

#' Sample cells
#'
#' @param i
#' @param cell.prob
#' @param max.cent
#' @param tot.cells
#' @param n
#'
#' @return
#'
#' @examples
BulkCells_subfxn <- function(i,cell.prob,max.cent,tot.cells,n){
  return(sample(1:tot.cells,n,replace = TRUE,prob = cell.prob[,i]))
}

#' Aggregate cells in bulk data based on the Gaussian kernal
#' @importFrom plyr count
#' @import Matrix
#' @importFrom tidyr gather
#' @param pick.cells
#' @param max.iter
#' @param step
#' @param max.cent
#' @param tot.cells
#' @param n
#' @param count
#'
#' @return
#'
#' @examples
BulkCells <- function(pick.cells, max.iter, step, max.cent, tot.cells, n, count, verbose){
  # if(verbose == TRUE){
  #   message("Bulk Cells")
  # }

  starting <- 1
  bulk.coef <- NULL
  bulk <- NULL
  if ((max.iter/step) == round(max.iter/step)){
    row.names(pick.cells) <- rep(1:step,(max.iter/step))
    for(i in seq(from = step, to = max.iter, by=step)){
      count.cells <- plyr::count(tidyr::gather(as.data.frame(t(pick.cells[starting:i,]))))
      gc(verbose = F)
      bulk_nn <- count %*% Matrix::sparseMatrix(i=count.cells[,2], j=as.numeric(count.cells[,1]),
                                        x=count.cells[,3], dims = c(tot.cells,step))
      bulk_nn1 <- t(t(bulk_nn)/colSums(bulk_nn))*1000000
      bulk <- c(bulk, bulk_nn1)
      rm(bulk_nn)
      starting <- i+1
    }
  } else if((max.iter/step) != round(max.iter/step)){
    row.names(pick.cells) <- c(rep(1:step,floor(max.iter/step)),1:(max.iter %% step))
    for(i in c(seq(from = step, to = max.iter,by=step), max.iter)){
      count.cells <- plyr::count(tidyr::gather(as.data.frame(t(pick.cells[starting:i,]))))
      gc(verbose = F)
      ##Change to the crossprod function
      bulk_nn <- count %*% Matrix::sparseMatrix(i=count.cells[,2], j=as.numeric(count.cells[,1]),
                                        x=count.cells[,3], dims = c(tot.cells,length(starting:i)))
      bulk_nn1 <- t(t(bulk_nn)/colSums(bulk_nn))*1000000
      bulk <- c(bulk, bulk_nn1)
      rm(bulk_nn)
      starting <- i+1
    }
  }
  bulk <- do.call(cbind, bulk)
  bulk <- as.matrix(bulk)
  return(bulk)
}
