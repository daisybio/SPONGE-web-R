#' Survival Analysis Log Rank pValues
#'
#' @description Retrieve pValues from log rank test based on raw survival analysis data
#'
#' @param disease_name The name of the dataset of interest as string. If default is set, all available datasets with corresponding informations are shown.
#'                     Fuzzy search is available.
#' @param ensg_number A vector of ensg number(s). If ensg_number is set, gene_symbol must be NULL.
#' @param gene_symbol A vector of gene symbol(s). If gene_symbol is set, ensg_number must be NULL.
#'
#' @return A data_frame with gene information and corresponding log rank test pValue. For raw data use function \code{\link[spongeAPI]{survAna_rates}}
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom dplyr mutate_at %>%
#' @importFrom httr GET content headers
#'
#' @examples
#'  # Retrieve gene expression values for specif genes by ensg_numbers
#' get_survAna_pValues(disease_name = "kidney clear cell carcinoma",
#'                     ensg_number = c("ENSG00000259090","ENSG00000217289"))
#'
#' # Retrieve gene expression values for specif genes by gene_symbols
#' get_survAna_pValues(disease_name = "kidney clear cell carcinoma",
#'                     gene_symbol = c("SEPT7P1","TIGAR"))
#'
#'  \dontrun{
#' # Ensg_numbers and gene_symbols together.
#' get_survAna_pValues(disease_name = "kidney clear cell carcinoma",
#'                     ensg_number = c("ENSG00000259090","ENSG00000217289"),
#'                     gene_symbol = c("SEPT7P1","TIGAR"))
#' }
get_survAna_pValues <- function(disease_name, ensg_number = NULL, gene_symbol = NULL){
  # Check if required paramter is given
  if(is.null(disease_name))
    stop("Required parameter disease_name is not given!")

  # Check if only one identifier is given.
  if(!is.null(gene_symbol) && !is.null(ensg_number))
    stop("More than one identification paramter is given. Please choose one out of ensg number or gene symbol.")

  # Base URL path
  base_url = "http://10.162.163.20:5000/survivalAnalysis/getPValues"

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
      mutate_at(c("pValue"), as.numeric) %>%
      mutate_at(c("dataset"), as.character)
    return(new)
  }
}

#' Survival Analysis Raw Data
#'
#' @description Get all raw survival analysis data for kaplan meier plots and
#'
#' @param disease_name The name of the dataset of interest as string. If default is set, all available datasets with corresponding informations are shown.
#'                     Fuzzy search is available.
#' @param ensg_number A vector of ensg number(s). If ensg_number is set, gene_symbol must be NULL. One of the two identifiers must be provided.
#' @param gene_symbol A vector of gene symbol(s). If gene_symbol is set, ensg_number must be NULL.
#' @param sample_ID A vector of sample_ID of the patient/sample of interest.
#'
#' @return A data_frame with gene and patient/sample information and the "group information" encoded by column "overexpressed".
#' Information about expression value of the gene (FALSE = underexpression, gene expression <= mean gene expression over all samples,
#'                                                 TRUE = overexpression, gene expression >= mean gene expression over all samples)
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom dplyr mutate_at %>%
#' @importFrom httr GET content headers
#'
#' @examples
#'  # Retrieve survival rates for specif genes by ensg_numbers and specific patient/sample
#' get_survAna_rates(disease_name="kidney clear cell carcinoma",
#'                   ensg_number=c("ENSG00000259090", "ENSG00000217289"),
#'                   sample_ID = c("TCGA-BP-4968","TCGA-B8-A54F"))
#'  \dontrun{
#' # Ensg_numbers and gene_symbols together.
#' get_survAna_pValues(disease_name = "kidney clear cell carcinoma",
#'                     ensg_number = c("ENSG00000259090","ENSG00000217289"),
#'                     gene_symbol = c("SEPT7P1","TIGAR"))
#' }
get_survAna_rates <- function(disease_name, ensg_number = NULL, gene_symbol = NULL, sample_ID = NULL){
  # Check if required paramter is given
  if(is.null(disease_name))
    stop("Required parameter disease_name is not given!")

  # Check if only one identifier is given.
  if(!is.null(gene_symbol) && !is.null(ensg_number))
    stop("More than one identification paramter is given. Please choose one out of ensg number or gene symbol.")

  # Base URL path
  base_url = "http://10.162.163.20:5000/survivalAnalysis/getRates"

  # Create full url
  if (!is.null(ensg_number)){
    full_url = paste(base_url, "?disease_name=", disease_name,"&ensg_number=", paste(ensg_number, collapse=",", sep=""), sep="")
  } else if(!is.null(gene_symbol)) {
    full_url = paste(base_url, "?disease_name=", disease_name,"&gene_symbol=", paste(gene_symbol, collapse=",", sep=""), sep="")
  } else
    full_url = paste(base_url, "?disease_name=", disease_name, sep="")
  if(!is.null(sample_ID)){
    full_url = paste(full_url, "&sample_ID=", paste(sample_ID, collapse = ",", sep=""), sep="")
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

    # Turn columns to numeric and remove NA values
    new <- new %>%
      mutate_at(c("overexpression", "patient_information.disease_status"), as.logical) %>%
      mutate_at(c("patient_information.survival_time"), as.numeric) %>%
      mutate_at(c("dataset"), as.character)
    return(new)
  }
}

#' Clinical Data For Samples/Patients
#'
#' @param disease_name The name of the dataset of interest as string. If default is set, all available datasets with corresponding informations are shown.
#'                     Fuzzy search is available.
#' @param sample_ID A vector of sample_ID of the patient/sample of interest.
#'
#' @return A data_frame with clinal date of all available or specific sample/patient.
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom dplyr mutate_at %>%
#' @importFrom httr GET content headers
#'
#' @examples
#' get_survAna_sampleInformation(disease_name = "kidney clear cell carcinoma",
#'                              sample_ID = c("TCGA-BP-4968","TCGA-B8-A54F"))
#' \dontrun{
#' # Do not run function without disease.
#' get_survAna_sampleInformation(sample_ID = c("TCGA-BP-4968","TCGA-B8-A54F"))
#' }
get_survAna_sampleInformation <- function(disease_name, sample_ID = NULL){
  # Check if required paramter is given
  if(is.null(disease_name))
    stop("Required parameter disease_name is not given!")

  # Base URL path
  base_url = "http://10.162.163.20:5000/survivalAnalysis/sampleInformation"

  # Create full url
  full_url = paste(base_url, "?disease_name=", disease_name, sep="")
  if(!is.null(sample_ID)){
    full_url = paste(full_url, "&sample_ID=", paste(sample_ID, collapse = ",", sep=""), sep="")
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

    # Turn columns to numeric and remove NA values
    new <- new %>%
      mutate_at(c("disease_status"), as.logical) %>%
      mutate_at(c("survival_time"), as.numeric) %>%
      mutate_at(c("dataset", "sample_ID"), as.character)
    return(new)
  }


}
