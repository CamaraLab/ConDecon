# Tests for functions in TransferFeatures.R
# ------------------------------------------------------------------------------

data("counts_gps")
data("latent_gps")
data("variable_genes_gps")
data("bulk_gps")

ConDecon_obj = RunConDecon(counts = counts_gps, latent = latent_gps, bulk = bulk_gps,
                           variable.features = variable_genes_gps, dims = 10, max.iter = 50,
                           verbose = FALSE)

test_that("Check that ConDecon's feature transfer has the correct dims", {
  random_gene = counts_gps[sample(x = 1:nrow(counts_gps), size = 1),]
  ConDecon_obj = TransferFeatures(ConDecon_obj = ConDecon_obj, feature = random_gene)
  expect_equal(ncol(ConDecon_obj$TransferFeatures), ncol(bulk_gps))
})
