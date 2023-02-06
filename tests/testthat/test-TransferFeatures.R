# Tests for functions in TransferFeatures.R
# ------------------------------------------------------------------------------

data("counts_gps")
data("latent_gps")
data("variable_genes_gps")
data("bulk_gps")

ConDecon_obj = RunConDecon(counts = counts_gps, latent = latent_gps, bulk = bulk_gps,
                           variable.features = variable_genes_gps, dims = 10, max.iter = 50,
                           verbose = FALSE)
random_gene = counts_gps[sample(x = 1:nrow(counts_gps), size = 1),]

test_that("Check that ConDecon's feature transfer has the correct dims", {
  ConDecon_obj = TransferFeatures(ConDecon_obj = ConDecon_obj, feature = random_gene)
  expect_equal(ncol(ConDecon_obj$TransferFeatures), ncol(bulk_gps))
})

test_that("Check that feature = matrix returns NULL", {
  ConDecon_obj = TransferFeatures(ConDecon_obj = ConDecon_obj, feature = matrix(c(1,2,3,4), nrow = 2))
  expect_equal(ConDecon_obj$TransferFeatures, NULL)
})

test_that("Check that length(feature) != number of cells", {
  ConDecon_obj = TransferFeatures(ConDecon_obj = ConDecon_obj, feature = c(1,2,3,4))
  expect_equal(ConDecon_obj$TransferFeatures, NULL)
})

test_that("Check that when features does not have names, they are auto added", {
  ConDecon_obj = TransferFeatures(ConDecon_obj = ConDecon_obj, feature = as.vector(random_gene))
  expect_equal(ncol(ConDecon_obj$TransferFeatures), ncol(bulk_gps))
})

test_that("Check that if all features are NA, returns NULL", {
  ConDecon_obj = TransferFeatures(ConDecon_obj = ConDecon_obj, feature = rep(NA, length(random_gene)))
  expect_equal(ConDecon_obj$TransferFeatures, NULL)
})


test_that("Check that if some features are NA, returns the correct dims", {
  random_gene[c(1,33,55, 103, 250, 300)] <- NA
  ConDecon_obj = TransferFeatures(ConDecon_obj = ConDecon_obj, feature = random_gene)
  expect_equal(ncol(ConDecon_obj$TransferFeatures), ncol(bulk_gps))
})
