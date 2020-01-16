#' miRNA Induced ceRNA Interaction
#'
#' @description Get all ceRNA interactions where miRNA(s) of interest (different identifiers available - e.g. hs number or mimat number) contribute to.
#'
#' @param disease_name Name of the specific cancer type/dataset. If default is set, all available datasets with corresponding informations are shown.
#'                     Fuzzy search available.
#' @param mimat_number A vector of mimat_number(s). If mimat_number is set, hs_number must be NULL.
#' @param hs_number A vector of hs_number(s). If hs_number is set, mimat_number must be NULL.
#' @param limit Number of results that should be shown. Default value is 100 and can be up to 1000.
#'              For more results please use batches, the provided offset parameter or download the whole dataset.
#' @param offset Starting point from where results should be shown.
#' @param information All available information about miRNA displayed (TRUE) or just ensg number (FALSE, default).
#'
#' @return A data_frame containing all ceRNA interactions fitting the paramters.
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
#'                               mimat_number = c("MIMAT0000076", "MIMAT0000261"),
#'                               limit = 15, information = FALSE)
#' \dontrun{
#' # Do not use both possible identifiers at once
#' get_specific_miRNAInteraction(disease_name = "kidney clear cell carcinoma",
#'                               mimat_number =c("MIMAT0000076", "MIMAT0000261"),
#'                               hs_number = c("hsa-miR-21-5p", "hsa-miR-183-5p"),
#'                               limit = 15, information = FALSE)
#'}
get_specific_miRNAInteraction <- function(disease_name = NULL,
                                  mimat_number = NULL,
                                  hs_number = NULL,
                                  limit = 100,
                                  offset = NULL,
                                  information = FALSE){

  # all checks will be done from the API and its unit tests!

  # Base URL path
  base_url = "http://10.162.163.20:5000/miRNAInteraction/findSpecific?"
  full_url = base_url

  # Create full url
  if (!is.null(disease_name)){
    full_url = paste(full_url, "disease_name=", disease_name,"&", sep="")
  }
  if (!is.null(mimat_number)){
    full_url = paste(full_url, "mimat_number=", paste(mimat_number, collapse=",", sep=""), "&", sep="")
  }
  if(!is.null(hs_number)) {
    full_url = paste(full_url, "hs_number=", paste(hs_number, collapse=",", sep=""), "&", sep="")
  }
  if (!is.null(limit)){
    full_url <- paste(full_url, "limit=", limit, "&", sep="")
  }
  if(!is.null(offset)){
    full_url <- paste(full_url, "offset=", offset, "&", sep="")
  }
  if(!is.logical(information)){
    stop(paste("Information is not logical!"))
  } else {
    full_url <- paste(full_url, "information=", information, "&", sep="")
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
    if(information){
      new <- new %>%
        mutate_at(c("interactions_genegene.run.run_ID"), as.numeric) %>%
        mutate_at(c("mirna.hs_nr", "mirna.id_type", "mirna.mir_ID", "mirna.seq"), as.character)
    } else {
      new <- new %>%
        mutate_at(c("interactions_genegene.run.run_ID"), as.numeric) %>%
        mutate_at(c("mirna.hs_nr", "mirna.mir_ID"), as.character)
    }
    return(new)
  }
}
