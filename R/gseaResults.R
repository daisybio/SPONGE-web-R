#' Get gsea results
#'
#' @description Get gsea results for the given diseases, conditions, gene set and terms.
#'
#' @param disease_name_1 Name of the disease type/dataset for the first part of the comparison. Fuzzy search available.
#' @param disease_name_2 Name of the disease type/dataset for the second part of the comparison. Fuzzy search available.
#' @param condition_1 Condition for the first part of the comparison (e.g. disease or normal).
#' @param condition_2 Condition for the second part of the comparison (e.g. disease or normal).
#' @param gene_set Name of gene set that contains the terms which are retrieved.
#' @param disease_subtype_1 Name of the disease subtype for the first part of the comparison. Overtype is selected if it is not given.
#' @param disease_subtype_2 Name of the disease subtype for the second part of the comparison. Overtype is selected if it is not given.
#' @param term Names of the terms for which to get the results. Gets results for all terms if it is not given.
#'
#' @return A list of data_frames containing all results for the given diseases, conditions, gene set and terms.
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom dplyr %>%
#' @importFrom httr GET content headers
#'
#' @examples
#' # Retrieve gsea results for disease names, subtypes, conditions, gene set
#' # and terms
#' get_gseaResults(disease_name_1 = "liver", disease_name_2 = "thymoma",
#'              condition_1 = "disease", condition_2 = "disease",
#'              gene_set = "GO_Biological_Process_2023", term="GO:0001676")


get_gseaResults = function(disease_name_1, disease_name_2, condition_1, condition_2, gene_set, disease_subtype_1 = NULL, disease_subtype_2 = NULL, term = NULL){
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
  base_url = paste(pkg.env$API.url, "/gseaResults", sep="")

  #Create full url
  full_url = paste(base_url, "?disease_name_1=", disease_name_1, "&disease_name_2=", disease_name_2, "&condition_1=", condition_1, "&condition_2=", condition_2, "&gene_set=", gene_set, sep="")

  if (!is.null(disease_subtype_1)){
    full_url = paste(full_url, "&disease_subtype_1=", disease_subtype_1, sep="")
  }
  if (!is.null(disease_subtype_2)){
    full_url = paste(full_url, "&disease_subtype_2=", disease_subtype_2, sep="")
  }
  if (!is.null(term)){
    full_url = paste(full_url, "&term=", paste(term, collapse=",", sep=""), sep="")
  }


  # Encode the URL with characters for each space.
  full_url <- URLencode(full_url)

  # Convert to text object using httr
  url_obj <- GET(full_url)

  # Parse url object
  raise <- content(url_obj, as="text", encoding = "UTF-8")

  #parse JSON
  new <- fromJSON(raise)


  res <- new$res
  names(res) <- new$term
  res <- do.call("rbind", res)
  res["term"] <- sapply(strsplit(rownames(res), "\\."), "[[", 1)
  rownames(res) <- NULL
  new$res <- NULL

  lead_genes <- new$lead_genes
  names(lead_genes) <- new$term
  new$lead_genes <- NULL

  matched_genes <- new$matched_genes
  names(matched_genes) <- new$term
  new$matched_genes <- NULL

  # Determine if a url object returns '404 Not Found'
  if(headers(url_obj)$`content-type` == "application/problem+json")
    stop(paste("API response is empty. Reason: ", new$detail))
  else {
    # Flatten out nested elements
    new <- do.call("data.frame", new)

    new <- new %>%
      mutate_at(c("es", "nes", "pvalue", "fdr", "fwerp", "gene_percent"), as.numeric)

    return(list(result=new, res=res, lead_genes=lead_genes, matched_genes=matched_genes))
  }
}
