#' Find Associated GO Terms
#'
#' @description Associated GO terms for gene(s) of interest. QuickGO - a fast web-based browser of the Gene Ontology and Gene Ontology annotation data - is used as external source.
#'
#' @param gene_symbol A vector of gene symbol(s). Required parameter.
#'
#' @return A data_frame containing all associated GO terms.
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom dplyr mutate_at %>%
#' @importFrom httr GET content headers
#'
#' @examples
#' get_geneOntology(gene_symbol=c("PTEN","TIGAR"))
get_geneOntology <- function(gene_symbol){

  # Base URL path
  base_url = paste(pkg.env$API.url, "/getGeneOntology?", sep="")
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
