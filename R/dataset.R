#' Dataset Information Function
#'
#' @description Get information about all available datasets to start browsing or search for a specific cancer type/dataset.
#'
#' @param disease_name The name of the dataset of interest as string. If default is set, all available datasets with corresponding informations are shown.
#'                     Fuzzy search is available (e.g. "kidney clear cell carcinoma" or just "kidney").
#'
#' @return Information about all or dataset disease_name disease_nameas data_frame.
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom dplyr mutate_at %>%
#' @importFrom httr GET content headers
#'
#' @examples
#' get_datasetInformation("kidney clear cell carcinoma")
#'
#' \dontrun{
#' get_datasetInformation(c("kidney clear cell carcinoma", "kidney papillary cell carcinoma"))
#' }
get_datasetInformation <- function(disease_name=NULL) {
  # Base URL path
  base_url = paste(pkg.env$API.url, "/dataset", sep="")

  # Create full url
  if (!is.null(disease_name)){
    full_url = paste0(base_url, "?disease_name=", disease_name)
  } else {
    full_url = base_url
  }

  # Encode the URL with characters for each space.
  full_url <- URLencode(full_url)

  # Convert to text object using httr
  url_obj <- GET(full_url)

  # Parse url object
  raise <- content(url_obj, as="text", encoding = "UTF-8")

  #parse JSON
  new <- fromJSON(raise)

  # Determine if a url object returns '404 Not Found'
  if(headers(url_obj)$`content-type` == "application/problem+json")
    stop(paste("API response is empty. Reason: ", new$detail))
  else {
    # Turn columns to numeric and remove NA values
    new <- new %>%
      mutate_at(c( "dataset_ID"), as.numeric)
    return(new)
  }
}

#' Run Information Function
#'
#' @description Retrieve all used parameters of the SPONGE method to create published results for the cancer type/dataset of interest.
#'
#' @param disease_name Name of the specific cancer type/dataset as string.
#' Fuzzy search is available (e.g. "kidney clear cell carcinoma" or just "kidney").
#'
#' @return Run information about dataset disease_name as data_frame.
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom dplyr mutate_at %>%
#' @importFrom httr GET content headers
#'
#' @examples
#' get_runInformation("kidney clear cell carcinoma")
#'
#' \dontrun{
#' get_runInformation(c("kidney clear cell carcinoma", "kidney papillary cell carcinoma"))
#' }
get_runInformation <- function(disease_name){
  # check if required paramter is given
  if(is.null(disease_name))
    stop("Required parameter disease_name is not given!")

  # Base URL path
  base_url = paste(pkg.env$API.url, "/dataset/runInformation", sep="")
  full_url = paste0(base_url, "?disease_name=", disease_name)

  # Encode the URL with characters for each space.
  full_url <- URLencode(full_url)

  # Convert to text object using httr
  url_obj <- GET(full_url)

  # Parse url object
  raise <- content(url_obj, as="text", encoding = "UTF-8")

  #parse JSON
  new <- fromJSON(raise)

  # Determine if a url object returns '404 Not Found'
  if(headers(url_obj)$`content-type` == "application/problem+json")
    stop(paste("API response is empty. Reason: ", new$detail))
  else {
    # Flatten out nested elements
    new <- do.call("data.frame", new)

    # Turn columns to numeric and remove NA values
    new <- new %>%
      mutate_at(c("coefficient_threshold", "dataset.dataset_ID", "f_test", "f_test_p_adj_threshold", "m_max", "min_corr", "number_of_datasets", "number_of_samples", "run_ID"), as.numeric)
    new <- new %>%
      mutate_at(c("coefficient_direction", "ks", "log_level"), as.character)
    return(new)
  }
}
