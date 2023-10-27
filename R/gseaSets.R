#' Get gene sets
#'
#' @description Get all gene sets with results for the given diseases and conditions.
#'
#' @param disease_name_1 Name of the disease type/dataset for the first part of the comparison. Fuzzy search available.
#' @param disease_name_2 Name of the disease type/dataset for the second part of the comparison. Fuzzy search available.
#' @param condition_1 Condition for the first part of the comparison (e.g. disease or normal).
#' @param condition_2 Condition for the second part of the comparison (e.g. disease or normal).
#' @param disease_subtype_1 Name of the disease subtype for the first part of the comparison. Overtype is selected if it is not given.
#' @param disease_subtype_2 Name of the disease subtype for the second part of the comparison. Overtype is selected if it is not given.
#'
#' @return A data_frame containing all gene sets with results for the given disease and conditions.
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom httr GET content headers
#'
#' @examples
#' # Retrieve all gene sets with available results for disease names, subtypes
#' # and conditions
#' get_gseaSets(disease_name_1 = "liver", disease_name_2 = "thymoma",
#'              condition_1 = "disease", condition_2 = "disease")


get_gseaSets = function(disease_name_1, disease_name_2, condition_1, condition_2, disease_subtype_1 = NULL, disease_subtype_2 = NULL){
  # Check if required paramter is given
  if(is.null(disease_name_1))
    stop("Required parameter disease_name_1 is not given!")

  if(is.null(disease_name_2))
    stop("Required parameter disease_name_2 is not given!")

  if(is.null(condition_1))
    stop("Required parameter condition_1 is not given!")

  if(is.null(condition_2))
    stop("Required parameter condition_2 is not given!")


  # Base URL path
  base_url = paste(pkg.env$API.url, "/gseaSets", sep="")

  #Create full url
  full_url = paste(base_url, "?disease_name_1=", disease_name_1, "&disease_name_2=", disease_name_2, "&condition_1=", condition_1, "&condition_2=", condition_2, sep="")

  if (!is.null(disease_subtype_1)){
    full_url = paste(full_url, "&disease_subtype_1=", disease_subtype_1, sep="")
  }
  if (!is.null(disease_subtype_2)){
    full_url = paste(full_url, "&disease_subtype_2=", disease_subtype_2, sep="")
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
    # Flatten out nested elements
    new <- do.call("data.frame", new)

    return(new)
  }
}
