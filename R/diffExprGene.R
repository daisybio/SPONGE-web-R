#' Get differential expression information for genes
#'
#' @description Get differential expression information for the given comparison and gene.
#'
#' @param disease_name_1 Name of the disease type/dataset for the first part of the comparison. Fuzzy search available.
#' @param disease_name_2 Name of the disease type/dataset for the second part of the comparison. Fuzzy search available.
#' @param condition_1 Condition for the first part of the comparison (e.g. disease or normal).
#' @param condition_2 Condition for the second part of the comparison (e.g. disease or normal).
#' @param ensg_number A vector of ensg number(s). If ensg number is set, gene symbol must be NULL.
#' @param gene_symbol A vector of gene symbol(s). If gene symbol is set, ensg number must be NULL.
#' @param disease_subtype_1 Name of the disease subtype for the first part of the comparison. Overtype is selected if it is not given.
#' @param disease_subtype_2 Name of the disease subtype for the second part of the comparison. Overtype is selected if it is not given.
#'
#' @return A data_frame differential expression information for the given comparison and gene.
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom dplyr mutate_at %>%
#' @importFrom httr GET content headers
#'
#' @examples
#' # Retrieve all differential expression information for the given
#' # comparison and gene.
#' get_differentialExpressionGene(disease_name_1 = "liver",
#'                                disease_name_2 = "thymoma",
#'                                condition_1 = "disease",
#'                                condition_2 = "disease",
#'                                gene_symbol = "CYP2E1")

get_differentialExpressionGene = function(disease_name_1, disease_name_2, condition_1, condition_2, ensg_number = NULL, gene_symbol = NULL, disease_subtype_1 = NULL, disease_subtype_2 = NULL){
  # Check if required paramter is given
  if(is.null(disease_name_1))
    stop("Required parameter disease_name_1 is not given!")

  if(is.null(disease_name_2))
    stop("Required parameter disease_name_2 is not given!")

  if(is.null(condition_1))
    stop("Required parameter condition_1 is not given!")

  if(is.null(condition_2))
    stop("Required parameter condition_2 is not given!")


  # Check if only one identifier is given.
  if(!is.null(gene_symbol) && !is.null(ensg_number))
    stop("More than one identification paramter is given. Please choose one out of ensg number or gene symbol.")

  # Base URL path
  base_url = paste(pkg.env$API.url, "/differentialExpression", sep="")

  # Create full url
  if (!is.null(ensg_number)){
    full_url = paste(base_url, "?disease_name_1=", disease_name_1, "&disease_name_2=", disease_name_2, "&condition_1=", condition_1, "&condition_2=", condition_2, "&ensg_number=", paste(ensg_number, collapse=",", sep=""), sep="")
  } else if(!is.null(gene_symbol)) {
    full_url = paste(base_url, "?disease_name_1=", disease_name_1, "&disease_name_2=", disease_name_2, "&condition_1=", condition_1, "&condition_2=", condition_2, "&gene_symbol=", paste(gene_symbol, collapse=",", sep=""), sep="")
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
