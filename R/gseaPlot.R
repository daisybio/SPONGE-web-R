#' Get gsea enrichment plot
#'
#' @description Get gsea enrichment plot for the given diseases, conditions, gene set and term.
#'
#' @param disease_name_1 Name of the disease type/dataset for the first part of the comparison. Fuzzy search available.
#' @param disease_name_2 Name of the disease type/dataset for the second part of the comparison. Fuzzy search available.
#' @param condition_1 Condition for the first part of the comparison (e.g. disease or normal).
#' @param condition_2 Condition for the second part of the comparison (e.g. disease or normal).
#' @param gene_set Name of gene set that contains the terms for which to get the enrichment plot.
#' @param term Name of the term for which to get the enrichment plot
#' @param disease_subtype_1 Name of the disease subtype for the first part of the comparison. Overtype is selected if it is not given.
#' @param disease_subtype_2 Name of the disease subtype for the second part of the comparison. Overtype is selected if it is not given.
#'
#' @return A ggplot object containing the gsea enrichment plot.
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom dplyr %>%
#' @importFrom httr GET content headers
#' @importFrom ggplot2 ggplot annotation_raster
#' @importFrom base64enc base64decode
#' @importFrom png readPNG
#'
#' @examples
#' # Retrieve gsea enrichment plot for disease names, subtypes, conditions,
#' # gene set and term
#' get_gseaPlot(disease_name_1 = "liver", disease_name_2 = "thymoma",
#'              condition_1 = "disease", condition_2 = "disease",
#'              gene_set = "GO_Biological_Process_2023", term="GO:0001676")


get_gseaPlot = function(disease_name_1, disease_name_2, condition_1, condition_2, gene_set, term, disease_subtype_1 = NULL, disease_subtype_2 = NULL){
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
  #base_url = paste(pkg.env$API.url, "/getGseaTerms", sep="")
  base_url = paste("localhost:5000", "/gseaPlot", sep="")

  #Create full url
  full_url = paste(base_url, "?disease_name_1=", disease_name_1, "&disease_name_2=", disease_name_2, "&condition_1=", condition_1, "&condition_2=", condition_2, "&gene_set=", gene_set, "&term=", term, sep="")

  if (!is.null(disease_subtype_1)){
    full_url = paste(full_url, "&disease_subtype_1=", disease_subtype_1, sep="")
  }
  if (!is.null(disease_subtype_2)){
    full_url = paste(full_url, "&disease_subtype_2=", disease_subtype_2, sep="")
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

    img <- base64decode(what=new)
    img <- readPNG(img)
    img <- ggplot() +
      annotation_raster(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)

    return(img)
  }
}
