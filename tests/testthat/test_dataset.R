library(spongeWeb)
context("Check if dataset and run information are parsed correct")

test_that("Check if dataset information from API are parsed correct", {
  expect_equal(get_datasetInformation("kidney clear cell carcinoma"),
               data.frame(data_origin=c("TCGA"), dataset_ID=c(3), disease_name="kidney clear cell carcinoma",
                          disease_type=c("cancer"), download_url=c("https://portal.gdc.cancer.gov/projects/TCGA-CCSK"), stringsAsFactors=FALSE))
  expect_error(get_datasetInformation("foo"), "API response is empty.*")
  })

test_that("Check if run information from API are parsed correct", {
  expect_equal(get_runInformation("kidney clear cell carcinoma"),
               data.frame(coefficient_direction=c("<"), coefficient_threshold=c( -0.05), dataset.data_origin=c("TCGA"),dataset.dataset_ID=c(3),
                          dataset.disease_name="kidney clear cell carcinoma",f_test=c(0), f_test_p_adj_threshold=c(0.05), ks=c("seq(0.2, 0.9, 0.1)"),
                          log_level=c("ERROR"), m_max=c(8), min_corr=c(0.1), number_of_datasets=c(100000), number_of_samples=c(495), run_ID=c(8),
                          variance_cutoff=c(NA), stringsAsFactors=FALSE))
  expect_error(get_runInformation("foo"), "API response is empty.*")
})
