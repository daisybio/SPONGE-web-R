#' Get differential expression information for transcripts
#'
#' @description Get differential expression information for the given comparison and transcript
#'
#' @param disease_name_1 Name of the disease type/dataset for the first part of the comparison. Fuzzy search available.
#' @param disease_name_2 Name of the disease type/dataset for the second part of the comparison. Fuzzy search available.
#' @param condition_1 Condition for the first part of the comparison (e.g. disease or normal).
#' @param condition_2 Condition for the second part of the comparison (e.g. disease or normal).
#' @param enst_number A vector of enst number(s). If enst number is not set, information for all transcripts is provided.
#' @param disease_subtype_1 Name of the disease subtype for the first part of the comparison. Overtype is selected if it is not given.
#' @param disease_subtype_2 Name of the disease subtype for the second part of the comparison. Overtype is selected if it is not given.
#'
#' @return A data_frame differential expression information for the given comparison and transcript.
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom dplyr mutate_at %>%
#' @importFrom httr GET content headers
#'
#' @examples
#' # Retrieve all differential expression information for the given
#' # comparison and transcript.
#' get_differentialExpressionTranscript(disease_name_1 = "liver",
#'                                disease_name_2 = "thymoma",
#'                                condition_1 = "disease",
#'                                condition_2 = "disease")


get_differentialExpressionTranscript = function(disease_name_1, disease_name_2, condition_1, condition_2, enst_number = NULL, disease_subtype_1 = NULL, disease_subtype_2 = NULL){
  # Check if required paramter is given
  if(is.null(disease_name_1))
    stop("Required parameter disease_name_1 is not given!")

  if(is.null(disease_name_2))
    stop("Required parameter disease_name_2 is not given!")

  if(is.null(condition_1))
    stop("Required parameter condition_1 is not given!")

  if(is.null(condition_2))
    stop("Required parameter condition_2 is not given!")

  # Base URL path
  base_url = paste(pkg.env$API.url, "/differentialExpressionTranscript", sep="")

  # Create full url
  if (!is.null(enst_number)){
    full_url = paste(base_url, "?disease_name_1=", disease_name_1, "&disease_name_2=", disease_name_2, "&condition_1=", condition_1, "&condition_2=", condition_2, "&enst_number=", paste(enst_number, collapse=",", sep=""), sep="")
  } else{
    full_url = paste(base_url, "?disease_name_1=", disease_name_1, "&disease_name_2=", disease_name_2, "&condition_1=", condition_1, "&condition_2=", condition_2, sep="")
  }

  if (!is.null(disease_subtype_1)){
    full_url = paste(full_url, "?disease_subtype_1=", disease_subtype_1, sep="")
  }
  if (!is.null(disease_subtype_1)){
    full_url = paste(full_url, "?disease_subtype_2=", disease_subtype_2, sep="")
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
      mutate_at(c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj", "stat"), as.numeric)

    return(new)
  }
}
