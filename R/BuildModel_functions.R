#' Fit a regression model
#' @description Learn a model that maps between the space of cell abundances and correlations
#' @importFrom stats lm as.formula
#' @param output ConDecon object with low dimensional embedding of distributions of cell abundances and correlations
#' @param degree degree of the polynomial regression model (default = 1)
#'
#' @return ConDecon object with fit model
#' @export
#'
#' @examples
#' data(counts_gps)
#' data(latent_gps)
#' data(variable_genes_gps)
#' data(bulk_gps)
#'
#' # For this example, we will reduce the training size to max.iter = 50 to reduce run time
#' TrainingSet = BuildTrainingSet(count = counts_gps, latent = latent_gps, max.iter = 50)
#' ConDecon_obj = Map2Latent(output = TrainingSet, latent = latent_gps, count = counts_gps,
#' bulk = bulk_gps, variable.features = variable_genes_gps)
#'
#' ConDecon_obj = BuildModel(ConDecon_obj)
BuildModel <- function(output,
                       degree = 1){

  #Create a fit for each dimension
  model <- lapply(1:(output$TrainingSet$dims), MultiVarFit, cell.prob.coef = output$TrainingSet$cell.prob_coefficients,
                  bulk.coef = output$TrainingSet$bulk_coefficients, degree = degree, dims = output$TrainingSet$dims)

  if(sum(is.na(model[[1]])) > 0){
    cat("The training set is too small to support a polynomial of degree ",degree, "\n")
    cat("Consider increasing the size of the training set or decreasing the degree of the polynomial \n")
  }

  output$Model <- model
  names(output$Model) <- paste0("dim_", 1:output$TrainingSet$dims)

  return(output)

}

#' MultiVarFit
#' @importFrom stats lm as.formula
#' @param i The current dimension
#' @param cell.prob.coef The fit model coefficients for the space of cell abundances
#' @param bulk.coef The fit model coefficients for the space of rank correlations
#' @param degree The degree of the polynomial
#' @param dims The total number of dimensions of the latent space
#'
#' @return model from stats::lm
MultiVarFit <- function(i, cell.prob.coef, bulk.coef, degree, dims){

  row.names(bulk.coef) <- paste0("x",1:dims)
  tmp <- data.frame(t(bulk.coef), y = cell.prob.coef[i,])

  model <- stats::lm(stats::as.formula(paste0("y ~ polym(", paste(row.names(bulk.coef)[1:nrow(bulk.coef)],
                                                    collapse = ",")," ,degree =",degree,",
                                raw = T)")), data = tmp)

  return(model)
}
