#' Gene Expression Function
#'
#' @description Get all expression values for gene(s) of interest.
#'
#' @param disease_name The name of the dataset of interest as string. If default is set, all available datasets with corresponding informations are shown.
#'                     Fuzzy search is available (e.g. "kidney clear cell carcinoma" or just "kidney").
#' @param ensg_number A vector of ensg number(s). If ensg_number is set, gene_symbol must be NULL.
#' @param gene_symbol A vector of gene symbol(s). If gene_symbol is set, ensg_number must be NULL.
#'
#' @return A data_frame with gene expression values.
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom dplyr mutate_at %>%
#' @importFrom httr GET content headers
#'
#' @examples
#' # Retrieve gene expression values for specif genes by ensg_numbers
#' get_geneExprValues(disease_name = "kidney clear cell carcinoma",
#'                    ensg_number = c("ENSG00000259090","ENSG00000217289"))
#'
#' # Retrieve gene expression values for specif genes by gene_symbols
#' get_geneExprValues(disease_name = "kidney clear cell carcinoma",
#'                    gene_symbol = c("SEPT7P1","TIGAR"))
#'
#'  \dontrun{
#' # Ensg_numbers and gene_symbols together.
#' get_geneExprValues(disease_name = "kidney clear cell carcinoma",
#'                    ensg_number = c("ENSG00000259090","ENSG00000217289"),
#'                    gene_symbol = c("SEPT7P1","TIGAR"))
#' }
get_geneExprValues = function(disease_name, ensg_number = NULL, gene_symbol = NULL){
  # Check if required paramter is given
  if(is.null(disease_name))
    stop("Required parameter disease_name is not given!")

  # Check if only one identifier is given.
  if(!is.null(gene_symbol) && !is.null(ensg_number))
    stop("More than one identification paramter is given. Please choose one out of ensg number or gene symbol.")

  # Base URL path
  base_url = "http://10.162.163.20:5000/exprValue/getceRNA"

  # Create full url
  if (!is.null(ensg_number)){
    full_url = paste(base_url, "?disease_name=", disease_name,"&ensg_number=", paste(ensg_number, collapse=",", sep=""), sep="")
  } else if(!is.null(gene_symbol)) {
   full_url = paste(base_url, "?disease_name=", disease_name,"&gene_symbol=", paste(gene_symbol, collapse=",", sep=""), sep="")
  } else
   full_url = paste(base_url, "?disease_name=", disease_name, sep="")

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
      mutate_at(c( "expr_value"), as.numeric) %>%
      mutate_at(c("dataset", "sample_ID"), as.character)
    return(new)
  }
}

#' miRNA Expression Function
#'
#' @description Get all expression values for miRNA(s) of interest.
#'
#' @param disease_name The name of the dataset of interest as string. If default is set, all available datasets with corresponding informations are shown.
#' Fuzzy search is available (e.g. "kidney clear cell carcinoma" or just "kidney").
#' @param mimat_number A vector of mimat_number(s). If mimat_number is set, hs_number must be NULL.
#' @param hs_number A vector of hs_number(s). If hs_number is set, mimat_number must be NULL.
#'
#' @return A data_frame with mirna expression values.
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom dplyr mutate_at %>%
#' @importFrom httr GET content headers
#'
#' @examples
#' # Retrieve gene expression values for specif miRNAs by mimat_numbers
#' get_mirnaExprValues(disease_name = "kidney clear cell carcinoma",
#'                     mimat_number = c("MIMAT0000076", "MIMAT0000261"))
#'
#' # Retrieve gene expression values for specif miRNAs by hs_numbers
#' get_mirnaExprValues(disease_name = "kidney clear cell carcinoma",
#'                     hs_number = c("hsa-miR-21-5p", "hsa-miR-183-5p"))
#'
#'  \dontrun{
#' # Mimat_numbers and hs_numbers together.
#' get_mirnaExprValues(disease_name = "kidney clear cell carcinoma",
#'                     mimat_number = c("MIMAT0000076", "MIMAT0000261"),
#'                     hs_number = c("hsa-miR-21-5p", "hsa-miR-183-5p"))
#' }
get_mirnaExprValues = function(disease_name, mimat_number = NULL, hs_number = NULL){
  # Check if required paramter is given
  if(is.null(disease_name))
    stop("Required parameter disease_name is not given!")

  # Check if only one identifier is given.
  if(!is.null(hs_number ) && !is.null(mimat_number))
    stop("More than one identification paramter is given. Please choose one out of mimat_number or hs_number.")

  # Base URL path
  base_url = "http://10.162.163.20:5000/exprValue/getmirNA"

  # Create full url
  if (!is.null(mimat_number)){
    full_url = paste(base_url, "?disease_name=", disease_name,"&mimat_number=", paste(mimat_number, collapse=",", sep=""), sep="")
  } else if(!is.null(hs_number )) {
    full_url = paste(base_url, "?disease_name=", disease_name,"&hs_number=", paste(hs_number , collapse=",", sep=""), sep="")
  } else
    full_url = paste(base_url, "?disease_name=", disease_name, sep="")

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
      mutate_at(c( "expr_value"), as.numeric) %>%
      mutate_at(c("dataset", "sample_ID"), as.character)
    return(new)
  }
}
