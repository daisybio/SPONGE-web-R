#' Overall Count for Interactions and miRNAs
#'
#' @description Function return current statistic about database - amount of shared miRNA, significant and insignificant interactions per dataset
#'
#' @param sponge_db_version The version of the SPONGE database to be used. Default is set $pkg.env$LATEST.
#' 
#' @return Overview of interaction counts about all or specific dataset as data_frame.
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom dplyr mutate_at %>%
#' @importFrom httr GET content headers
#'
#' @examples
#' get_overallCounts()
get_overallCounts <- function(sponge_db_version = pkg.env$LATEST) {
  # Base URL path
  full_url = paste(pkg.env$API.url, "/getOverallCounts", sep="")

  full_url <- paste(full_url, "?sponge_db_version=", sponge_db_version, sep="")

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
      mutate_at(c( "run_ID", "count_interactions","count_interactions_sign","count_shared_miRNAs"), as.numeric)
    return(new)
  }
}
