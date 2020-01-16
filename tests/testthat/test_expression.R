library(spongeWeb)
context("Expression Values Parsed Correctly")

test_that("Gene Expression Values from API parsed correctly", {
  expect_error(get_geneExprValues(disease_name = "bla"), "API response is empty.*")
  expect_error(get_geneExprValues(disease_name = "kidney clear cell carcinoma", ensg_number = c("ENSG001bla")), "API response is empty.*")
  expect_error(get_geneExprValues(disease_name = "kidney clear cell carcinoma", gene_symbol = c("FooBar")), "API response is empty.*")
  expect_error(get_geneExprValues(disease_name = "kidney clear cell carcinoma", ensg_number = c("ENSG001bla"), gene_symbol = c("FooBar")), "More than one identification paramter is given.*")
  expect_error(get_geneExprValues(disease_name = NULL), "Required parameter disease_name is not given*")
  expect_error(get_geneExprValues(), "argument \"disease_name\" is missing, with no default*")
})

test_that("miRNA Expression Values from API parsed correctly", {
  expect_error(get_mirnaExprValues(disease_name = "bla"), "API response is empty.*")
  expect_error(get_mirnaExprValues(disease_name = "kidney clear cell carcinoma", mimat_number = c("MIR-bla")), "API response is empty.*")
  expect_error(get_mirnaExprValues(disease_name = "kidney clear cell carcinoma", hs_number = c("FooBar")), "API response is empty.*")
  expect_error(get_mirnaExprValues(disease_name = "kidney clear cell carcinoma", mimat_number = c("MIR-bla"), hs_number = c("FooBar")), "More than one identification paramter is given.*")
  expect_error(get_mirnaExprValues(disease_name = NULL), "Required parameter disease_name is not given*")
  expect_error(get_mirnaExprValues(), "argument \"disease_name\" is missing, with no default*")
})
