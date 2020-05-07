#' Find Associated Cancer hallmarks.
#'
#' @description Associated cancer hallmark for gene(s) of interest. Cancer Hallmark Genes (http://bio-bigdata.hrbmu.edu.cn/CHG/) is used as external source.
#'
#' @param gene_symbol A vector of gene symbol(s). Required parameter.
#'
#' @return A data_frame containing all associated cancer hallmarks.
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom dplyr mutate_at %>%
#' @importFrom httr GET content headers
#'
#' @examples
#' get_hallmark(gene_symbol=c("PTEN"))
get_hallmark <- function(gene_symbol){

  # Base URL path
  base_url = paste(pkg.env$API.url, "/getHallmark?", sep="")
  full_url = base_url

  # Create full url
  if(!is.null(gene_symbol)) {
    full_url = paste(full_url, "gene_symbol=", paste(gene_symbol, collapse=",", sep=""), "&", sep="")
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
    new <- do.call("data.frame", do.call("data.frame", new))
    return(new)
  }
}
