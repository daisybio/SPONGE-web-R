#' Occurences of Sponges miRNAs
#'
#' @description Get all mirna involved in cancer type/dataset of interest occuring a certain amout of times.
#'
#' @param disease_name Name of the specific cancer type/dataset. If default is set, all available datasets with corresponding informations are shown.
#'                     Fuzzy search available.
#' @param mimat_number A vector of mimat_number(s). If mimat_number is set, hs_number must be NULL.
#' @param hs_number A vector of hs_number(s). If hs_number is set, mimat_number must be NULL.
#' @param occurences Threshold of amount of contributions a miRNA should be considered.
#' @param sorting Possibilities for sorting of the results. Possible values are "pValue", "mscor" or "correlation".
#' @param descending Descending (TRUE, default) or ascending (FALSE) ordering of the results.
#' @param limit Number of results that should be shown. Default value is 100 and can be up to 1000.
#'              For more results please use batches, the provided offset parameter or download the whole dataset.
#' @param offset Starting point from where results should be shown.
#' @param sponge_db_version The version of the SPONGE database to be used. Default is set $pkg.env$LATEST.
#' 
#' @return A data_frame with all miRNAs occuring at leat "occurences" times.
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom dplyr mutate_at %>%
#' @importFrom httr GET content headers
#'
#' @examples
#' get_miRNAOccurences(disease_name="kidney clear cell carcinoma", occurences = 1000, limit = 10)
get_miRNAOccurences <- function(disease_name,
                                         mimat_number = NULL,
                                         hs_number = NULL,
                                         occurences = NULL,
                                         sorting = NULL,
                                         descending = TRUE,
                                         limit = 100,
                                         offset = NULL,
                                         sponge_db_version = pkg.env$LATEST) {
  # all checks will be done from the API and its unit tests!

  # Base URL path
  base_url = paste(pkg.env$API.url, "/miRNAInteraction/getOccurence?", sep="")
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
  if(!is.null(occurences)) {
    full_url = paste(full_url, "occurences=", occurences,"&", sep="")
  }
  if(!is.null(sorting)){
    if(sorting %in% c("pValue", "mscor", "correlation"))
      full_url <- paste(full_url, "sorting=", sorting, "&", sep="")
    else
      stop(paste("sorting:", sorting," is not an allowed value. Please check the help page for further information."))
  }
  if(!is.logical(descending)){
    stop(paste("Descending is not logical!"))
  } else {
    full_url <- paste(full_url, "descending=", descending, "&", sep="")
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
        mutate_at(c("occurences", "run.dataset.dataset_ID", "run.run_ID"), as.numeric) %>%
        mutate_at(c("mirna.hs_nr", "mirna.mir_ID", "run.dataset.data_origin", "run.dataset.disease_name"), as.character)
    return(new)
  }

}
