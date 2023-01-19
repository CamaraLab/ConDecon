# Tests for functions in PlotConDecon.R
# ------------------------------------------------------------------------------

data("counts_gps")
data("latent_gps")
data("variable_genes_gps")
data("bulk_gps")

ConDecon_obj = RunConDecon(counts = counts_gps, latent = latent_gps, bulk = bulk_gps,
                variable.features = variable_genes_gps, dims = 10, max.iter = 50,
                verbose = FALSE)

data("umap_embedding_gps")
ConDecon_obj = PlotConDecon(ConDecon_obj = ConDecon_obj, umap = umap_embedding_gps)

test_that("Check that ConDecon's predicted relative cell abundance matrix has the correct dims", {
  expect_equal(dim(ConDecon_obj$Relative_cell.prob), c(nrow(ConDecon_obj$Normalized_cell.prob), ncol(ConDecon_obj$Normalized_cell.prob)))
})


