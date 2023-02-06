# Tests for functions in PlotConDecon.R
# ------------------------------------------------------------------------------

data("counts_gps")
data("latent_gps")
data("variable_genes_gps")
data("bulk_gps")

ConDecon_obj = RunConDecon(counts = counts_gps, latent = latent_gps, bulk = bulk_gps,
                variable.features = variable_genes_gps, dims = 10, max.iter = 50,
                verbose = FALSE)

data("meta_data_gps")

test_that("Check that ConDecon's predicted relative cell abundance matrix has the correct dims", {
  ConDecon_obj = PlotConDecon(ConDecon_obj = ConDecon_obj, umap = meta_data_gps[,c("UMAP_1", "UMAP_2")])
  expect_equal(dim(ConDecon_obj$Relative_cell.prob), c(nrow(ConDecon_obj$Normalized_cell.prob), ncol(ConDecon_obj$Normalized_cell.prob)))
  expect_equal(class(expect_invisible(PlotConDecon(ConDecon_obj = ConDecon_obj, umap = meta_data_gps[,c("UMAP_1", "UMAP_2")]))), "ConDecon")
})

test_that("Check character vector input into PlotConDecon samples parameter", {
  ConDecon_obj = PlotConDecon(ConDecon_obj = ConDecon_obj, umap = meta_data_gps[,c("UMAP_1", "UMAP_2")], samples = c("bulk1"))
  ConDecon_obj = CalcRelativeCellProb(ConDecon_obj, rm_outliers = TRUE)
  expect_equal(class(expect_invisible(PlotConDecon(ConDecon_obj = ConDecon_obj, umap = meta_data_gps[,c("UMAP_1", "UMAP_2")]))), "ConDecon")
})

test_that("Check numeric vector input into PlotConDecon samples parameter", {
  ConDecon_obj = PlotConDecon(ConDecon_obj = ConDecon_obj, umap = meta_data_gps[,c("UMAP_1", "UMAP_2")], samples = c(1))
  expect_equal(class(expect_invisible(PlotConDecon(ConDecon_obj = ConDecon_obj, umap = meta_data_gps[,c("UMAP_1", "UMAP_2")]))), "ConDecon")
})

test_that("Check character vector input into PlotConDecon title_names parameter", {
  ConDecon_obj = PlotConDecon(ConDecon_obj = ConDecon_obj, umap = meta_data_gps[,c("UMAP_1", "UMAP_2")], title_names = "bulk1", samples = "bulk1")
  expect_equal(class(expect_invisible(PlotConDecon(ConDecon_obj = ConDecon_obj, umap = meta_data_gps[,c("UMAP_1", "UMAP_2")]))), "ConDecon")
})

test_that("Check character vector input into PlotConDecon cells parameter", {
  ConDecon_obj = PlotConDecon(ConDecon_obj = ConDecon_obj, umap = meta_data_gps[,c("UMAP_1", "UMAP_2")], cells = c("Cell1", "Cell2", "Cell3", "Cell4", "Cell5"))
  expect_equal(class(expect_invisible(PlotConDecon(ConDecon_obj = ConDecon_obj, umap = meta_data_gps[,c("UMAP_1", "UMAP_2")]))), "ConDecon")
})

test_that("Check numeric vector input into PlotConDecon cells parameter", {
  ConDecon_obj = PlotConDecon(ConDecon_obj = ConDecon_obj, umap = meta_data_gps[,c("UMAP_1", "UMAP_2")], cells = c(1:5))
  expect_equal(class(expect_invisible(PlotConDecon(ConDecon_obj = ConDecon_obj, umap = meta_data_gps[,c("UMAP_1", "UMAP_2")]))), "ConDecon")
})
