# Tests for functions in RunConDecon.R
# ------------------------------------------------------------------------------

data("counts_gps")
data("latent_gps")
data("variable_genes_gps")
data("bulk_gps")

test_that("Check that dims of ConDecon's predicted cell abundances is correct", {
  x = RunConDecon(counts = counts_gps, latent = latent_gps, bulk = bulk_gps,
                  variable.features = variable_genes_gps, dims = 10, max.iter = 50,
                  verbose = FALSE)
  expect_equal(dim(x$norm_cell.prob), c(ncol(counts_gps), ncol(bulk_gps)))
})


