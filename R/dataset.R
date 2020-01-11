#' Dataset Information Function
#'
#' @param disease_name The name of the dataset of interest as string. If default is set, all available datasets with corresponding informations are shown. Fuzzy search is available (e.g. "kidney clear cell carcinoma" or just "kidney").
#'
#' @return Information about all or dataset disease_name disease_nameas data_frame.
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom dplyr mutate_at %>%
#'
#' @examples
#' datasetInformation("kidney clear cell carcinoma")
#' \dontrun{
#' datasetInformation(c("kidney clear cell carcinoma", "kidney papillary cell carcinoma"))
#' }
datasetInformation <- function(disease_name=NULL) {
  # Base URL path
  base_url = "http://10.162.163.20:5000/dataset"
  if (!is.null(disease_name)){
    full_url = paste0(base_url, "?disease_name=", disease_name)
  } else {
    full_url = base_url
  }

  # encode the URL with characters for each space.
  full_url <- URLencode(full_url)

  # Convert JSON to data frame
  dataset_df <- fromJSON(full_url)

  # turn columns to numeric and remove NA values
  dataset_df <- dataset_df %>%
    mutate_at(c( "dataset_ID"), as.numeric)

  return(dataset_df)
}

#' Run Information Function
#'
#' @description Retrieve all used parameters of the SPONGE method to create published results for the cancer type/dataset of interest.
#'
#' @param disease_name Name of the specific cancer type/dataset as string. Fuzzy search is available (e.g. "kidney clear cell carcinoma" or just "kidney").
#'
#' @return Run information about dataset disease_name as data_frame.
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom dplyr mutate_at %>%
#'
#' @examples
#' runInformation("kidney clear cell carcinoma")
#' \dontrun{
#' runInformation(c("kidney clear cell carcinoma", "kidney papillary cell carcinoma"))
#' }
runInformation <- function(disease_name){
  # check if required paramter is given
  try(if(is.null(disease_name)) stop("Required parameter disease_name is not given!"))

  # Base URL path
  base_url = "http://10.162.163.20:5000/dataset/runInformation"
  full_url = paste0(base_url, "?disease_name=", disease_name)

  # encode the URL with characters for each space.
  full_url <- URLencode(full_url)

  # Convert JSON to data frame
  run_df <- fromJSON(full_url)

  # turn columns to numeric and remove NA values
  #run_df <- run_df %>%
  #  mutate_at(c( "dataset_ID"), as.numeric)

  return(run_df)
}
