# Tests for functions in RunConDecon.R
# ------------------------------------------------------------------------------

data("counts_gps")
data("latent_gps")
data("variable_genes_gps")
data("bulk_gps")

test_that("Check that ConDecon's predicted cell abundance matrix has the correct dims", {
  ConDecon_obj = RunConDecon(counts = counts_gps, latent = latent_gps, bulk = bulk_gps,
                  variable.features = variable_genes_gps, dims = 10, max.iter = 50,
                  verbose = TRUE)
  expect_equal(dim(ConDecon_obj$Normalized_cell.prob), c(ncol(counts_gps), ncol(bulk_gps)))
})

test_that("For degree = 2, check that ConDecon's predicted cell abundance matrix has the correct dims", {
  ConDecon_obj = RunConDecon(counts = counts_gps, latent = latent_gps, bulk = bulk_gps,
                             variable.features = variable_genes_gps, dims = 10, max.iter = 100,
                             verbose = FALSE, degree = 2)
  expect_equal(dim(ConDecon_obj$Normalized_cell.prob), c(ncol(counts_gps), ncol(bulk_gps)))
})

test_that("For n = 250.5, check that ConDecon's predicted cell abundance matrix has the correct dims", {
  ConDecon_obj = RunConDecon(counts = counts_gps, latent = latent_gps, bulk = bulk_gps,
                             variable.features = variable_genes_gps, dims = 10, max.iter = 50,
                             verbose = FALSE, n = 250.5)
  expect_equal(dim(ConDecon_obj$Normalized_cell.prob), c(ncol(counts_gps), ncol(bulk_gps)))
})


