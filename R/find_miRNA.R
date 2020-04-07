#' Find Sponged miRNAs
#'
#' @description Get all miRNAs that contribute to all interactions between the given identifiers (ensg_number or gene_symbol).
#'
#' @param disease_name Name of the specific cancer type/dataset. If default is set, all available datasets with corresponding informations are shown.
#'                     Fuzzy search available.
#' @param ensg_number A vector of ensg number(s). If ensg number is set, gene symbol and gene type must be NULL. One of the three identifiers must be provided.
#' @param gene_symbol A vector of gene symbol(s). If gene symbol is set, ensg number and gene type must be NULL. One of the three identifiers must be provided.
#' @param between If false (default), all interactions where one of the interaction partners fits the given genes of interest
#'                will be considered. If true, just interactions between the genes of interest will be considered.
#'
#' @return A data_frame containing all found miRNAs.
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom dplyr mutate_at %>%
#' @importFrom httr GET content headers
#'
#' @examples
#' get_sponged_miRNA(disease_name="kidney", gene_symbol = c("TCF7L1", "SEMA4B"), between=TRUE)
get_sponged_miRNA <- function(disease_name,
                       ensg_number = NULL,
                       gene_symbol = NULL,
                       between = FALSE){

  # all checks will be done from the API and its unit tests!

  # Base URL path
  base_url = paste(pkg.env$API.url, "/miRNAInteraction/findceRNA?", sep="")
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
  if(!is.logical(between)){
    stop(paste("Between parameter is not logical!"))
  } else {
    full_url <- paste(full_url, "between=", between, "&", sep="")
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
      mutate_at(c("coefficient", "run.run_ID"), as.numeric) %>%
      mutate_at(c("mirna.hs_nr", "mirna.mir_ID", "gene.ensg_number","gene.gene_symbol", "run.dataset"), as.character)
    return(new)
  }
}
