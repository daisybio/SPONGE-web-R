#' Number of Times Gene Involved in Complete Network and Significant Interactions.
#'
#' @param disease_name The name of the dataset of interest as string. If default is set, all available datasets with corresponding informations are shown.
#'                     Fuzzy search is available.
#' @param ensg_number A vector of ensg number(s). If ensg_number is set, gene_symbol must be NULL.
#' @param gene_symbol A vector of gene symbol(s). If gene_symbol is set, ensg_number must be NULL.
#' @param minCountAll Defines the minimal number of times a gene has to be involved in the complete network (e.g. the degree of the corresponding node must be greater than minCountAll).
#' @param minCountSign Defines the minimal number of times a gene has to be involved in significant (p.adj < 0.05) interactions in the network.
#'
#' @return A data_frame cotaining the amount of times a gene is involved in the complete network (equals to degree), column count_all, and in significant (FDR adjusted pValue < 0.05) interactions of the network, column count_sign.
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom dplyr mutate_at %>%
#' @importFrom httr GET content headers
#'
#' @examples
#' # Get all genes from specific cancer with a minimum number of 150 at significant ceRNA interactions
#' get_geneCount(disease_name = "kidney clear cell carcinoma", minCountSign = 150)
get_geneCount <- function(disease_name = NULL,
                          ensg_number = NULL,
                          gene_symbol = NULL,
                          minCountAll = NULL,
                          minCountSign = NULL){

  # all checks will be done from the API and its unit tests!

  # Base URL path
  base_url = "http://10.162.163.20:5000/getGeneCount?"
  full_url = base_url

  # Create full url
  if (!is.null(disease_name)){
    full_url = paste(full_url, "disease_name=", disease_name,"&", sep="")
  }
  if (!is.null(ensg_number)){
    full_url = paste(full_url, "ensg_number=", paste(ensg_number, collapse=",", sep=""), "&", sep="")
  }
  if(!is.null(gene_symbol)) {
    full_url = paste(full_url, "gene_symbol=", paste(gene_symbol, collapse=",", sep=""), "&", sep="")
  }
  if (!is.null(minCountAll)){
    full_url <- paste(full_url, "minCountAll=", minCountAll, "&", sep="")
  }
  if(!is.null(minCountSign)){
    full_url <- paste(full_url, "minCountSign=", minCountSign, "&", sep="")
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

    # Turn columns to numeric and remove NA values
     new <- new %>%
      mutate_at(c("count_all", "count_sign", "run.dataset.dataset_ID", "run.run_ID"), as.numeric) %>%
      mutate_at(c("gene.gene_symbol", "gene.ensg_number", "run.dataset.data_origin", "run.dataset.disease_name"), as.character)
    return(new)
  }
}
