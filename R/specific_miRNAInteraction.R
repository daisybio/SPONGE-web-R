#' miRNA Induced ceRNA Interaction
#'
#' @description Get all ceRNA interactions where miRNA(s) of interest (different identifiers available - e.g. hs number or mimat number) contribute to.
#'
#' @param disease_name Name of the specific cancer type/dataset. If default is set, all available datasets with corresponding informations are shown.
#'                     Fuzzy search available.
#' @param mimat_number Mimat_number of interest. If mimat_number is set, hs_number must be NULL.
#' @param hs_number hs_number of interest. If hs_number is set, mimat_number must be NULL.
#' @param pValue Threshold of the FDR adjusted p-value. Default is 0.05.
#' @param pValueDirection Direction of the FDR adjusted p-value threshold (<, >). Must be set if pValue is set. Possible values are: "<", ">".
#' @param mscor Threshold of the 'multiple sensitivity correlation' (mscor).
#' @param mscorDirection Direction of the mscor threshold (<, >). Must be set if pValue is set. Possible values are: "<", ">".
#' @param correlation Threshold of the correlation.
#' @param correlationDirection Direction of the correlation threshold (<, >). Must be set if pValue is set. Possible values are: "<", ">".
#' @param limit Number of results that should be shown. Default value is 100 and can be up to 1000.
#'              For more results please use batches, the provided offset parameter or download the whole dataset.
#' @param offset Starting point from where results should be shown.
#' @param sponge_db_version The version of the SPONGE database to be used. Default is set $pkg.env$LATEST.
#'
#' @return A data_frame containing all ceRNA interactions fitting the parameters.
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom dplyr mutate_at %>%
#' @importFrom httr GET content headers
#'
#' @examples
#' # Retrieve all possible ceRNA interactions where miRNA of interest contribute to
#' get_specific_miRNAInteraction(disease_name = "kidney clear cell carcinoma",
#'                               mimat_number = "MIMAT0000076",
#'                               limit = 15)
#' \dontrun{
#' # Do not use both possible identifiers at once
#' get_specific_miRNAInteraction(disease_name = "kidney clear cell carcinoma",
#'                               mimat_number ="MIMAT0000076",
#'                               hs_number = "hsa-miR-21-5p",
#'                               limit = 15)
#'}
get_specific_miRNAInteraction <- function(disease_name = NULL,
                                  mimat_number = NULL,
                                  hs_number = NULL,
                                  pValue = 0.05,
                                  pValueDirection = "<",
                                  mscor = NULL,
                                  mscorDirection = "<",
                                  correlation = NULL,
                                  correlationDirection = "<",
                                  limit = 100,
                                  offset = NULL,
                                  sponge_db_version = pkg.env$LATEST) {

  # all checks will be done from the API and its unit tests!

  # Base URL path
  base_url = paste(pkg.env$API.url, "/miRNAInteraction/findSpecific?", sep="")
  full_url = base_url

  # Create full url
  if (!is.null(disease_name)){
    full_url = paste(full_url, "disease_name=", disease_name,"&", sep="")
  }
  if (!is.null(mimat_number)){
    full_url = paste(full_url, "mimat_number=", mimat_number, "&", sep="")
  }
  if(!is.null(hs_number)) {
    full_url = paste(full_url, "hs_number=", hs_number, "&", sep="")
  }
  if(!is.null(pValue)){
    full_url = paste(full_url, "pValue=", pValue,"&", sep="")
  }
  if(!is.null(pValueDirection)){
    if(pValueDirection %in% c("<",">"))
      full_url <- paste(full_url, "pValueDirection=", pValueDirection, "&", sep="")
    else
      stop(paste("pValueDirection:", pValueDirection," is not an allowed value. Please check the help page for further information."))
  }
  if(!is.null(mscor)){
    full_url = paste(full_url, "mscor=", mscor,"&", sep="")
  }
  if(!is.null(mscorDirection)){
    if(pValueDirection %in% c("<",">"))
      full_url <- paste(full_url, "mscorDirection=", mscorDirection, "&", sep="")
    else
      stop(paste("mscorDirection:", mscorDirection," is not an allowed value. Please check the help page for further information."))
  }
  if(!is.null(correlation)){
    full_url = paste(full_url, "correlation=", correlation,"&", sep="")
  }
  if(!is.null(correlationDirection)){
    if(pValueDirection %in% c("<",">"))
      full_url <- paste(full_url, "correlationDirection=", correlationDirection, "&", sep="")
    else
      stop(paste("correlationDirection:", correlationDirection," is not an allowed value. Please check the help page for further information."))
  }
  if (!is.null(limit)){
    full_url <- paste(full_url, "limit=", limit, "&", sep="")
  }
  if(!is.null(offset)){
    full_url <- paste(full_url, "offset=", offset, "&", sep="")
  }

  full_url <- paste(full_url, "sponge_db_version=", sponge_db_version, sep="")

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
      mutate_at(c("correlation", "p_value", "mscor", "run.dataset.dataset_ID", "run.run_ID"), as.numeric) %>%
      mutate_at(c("gene1.ensg_number", "gene1.gene_symbol", "gene2.ensg_number", "gene2.gene_symbol"), as.character)
    return(new)
  }
}
