#function for endpoint /findceRNA

#' Find Specific ceRNA (Gene)
#'
#' @description Get all ceRNAs in a disease of interest (search not for a specific ceRNA, but search for all ceRNAs satisfying filter functions).
#'
#' @param disease_name Name of the specific cancer type/dataset. If default is set, all available datasets with corresponding informations are shown.
#'                     Fuzzy search available.
#' @param gene_type Defines the type of gene of interest. One out of [3prime_overlapping_ncRNA, antisense, antisense_RNA, bidirectional_promoter_lncRNA, IG_C_gene, IG_C_pseudogene, IG_V_gene, IG_V_pseudogene, lincRNA, macro_lncRNA, miRNA, misc_RNA, Mt_rRNA, polymorphic_pseudogene, processed_pseudogene, processed_transcript, protein_coding, pseudogene, ribozyme, rRNA, rRNA_pseudogene, scaRNA, scRNA, sense_intronic, sense_overlapping, snoRNA, snRNA, TEC, TR_C_gene, TR_V_gene, TR_V_pseudogene, transcribed_processed_pseudogene, transcribed_unitary_pseudogene, transcribed_unprocessed_pseudogene, translated_processed_pseudogene, unitary_pseudogene, unprocessed_pseudogene, vaultRNA].
#' @param minBetweenness Threshold of the betweenness.
#' @param minNodeDegree Threshold of the degree.
#' @param minEigenvector Threshold of the eigenvektor.
#' @param sorting Possibilities for sorting of the results. Possible values are "degree", "betweenness" or "eigenvector".
#' @param descending Descending (TRUE, default) or ascending (FALSE) ordering of the results.
#' @param limit Number of results that should be shown. Default value is 100 and can be up to 1000.
#'              For more results please use batches, the provided offset parameter or download the whole dataset.
#' @param offset Starting point from where results should be shown.
#'
#' @return A data_frame containing all ceRNAs (genes) satisfying the given filters.
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom dplyr mutate_at %>%
#' @importFrom httr GET content headers
#'
#' @examples
#' get_ceRNA(disease_name = "kidney clear cell carcinoma",
#'           gene_type = "lincRNA", minBetweenness = 0.8)
get_ceRNA <- function(disease_name,
                      gene_type = NULL,
                      minBetweenness = NULL,
                      minNodeDegree = NULL,
                      minEigenvector = NULL,
                      sorting = NULL,
                      descending = TRUE,
                      limit = 100,
                      offset = NULL){
  # all checks will be done from the API and its unit tests!

  # Base URL path
  base_url = paste(pkg.env$API.url, "/findceRNA?", sep="")
  full_url = base_url

  # Create full url
  if (!is.null(disease_name)){
    full_url = paste(full_url, "disease_name=", disease_name,"&", sep="")
  }
  if(!is.null(gene_type)){
    types <- c("3prime_overlapping_ncRNA", "antisense", "antisense_RNA", "bidirectional_promoter_lncRNA", "IG_C_gene", "IG_C_pseudogene",
               "IG_V_gene", "IG_V_pseudogene", "lincRNA", "macro_lncRNA", "miRNA", "misc_RNA", "Mt_rRNA", "polymorphic_pseudogene",
               "processed_pseudogene", "processed_transcript", "protein_coding", "pseudogene", "ribozyme", "rRNA", "rRNA_pseudogene",
               "scaRNA", "scRNA", "sense_intronic", "sense_overlapping", "snoRNA", "snRNA", "TEC", "TR_C_gene", "TR_V_gene", "TR_V_pseudogene",
               "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene",
               "translated_processed_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene", "vaultRNA")
    if(gene_type %in% types)
      full_url = paste(full_url, "gene_type=", gene_type, "&", sep="")
    else
      stop(paste("Gene_type:", gene_type," is not an allowed value. Please check the help page for further information."))
  }
  if(!is.null(minBetweenness)){
    full_url <- paste(full_url, "minBetweenness=", minBetweenness, "&", sep="")
  }
  if(!is.null(minNodeDegree)){
    full_url <- paste(full_url, "minNodeDegree=", minNodeDegree, "&", sep="")
  }
  if(!is.null(minEigenvector)){
    full_url <- paste(full_url, "minEigenvector=", minEigenvector, "&", sep="")
  }
  if(!is.null(sorting)){
    if(sorting %in% c("degree", "betweenness", "eigenvector"))
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
      mutate_at(c("betweeness", "eigenvector", "node_degree", "run.dataset.dataset_ID", "run.run_ID"), as.numeric) %>%
      mutate_at(c("gene.ensg_number", "gene.gene_symbol"), as.character)
    return(new)
  }
}
