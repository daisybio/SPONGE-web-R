#' Find All Possible ceRNA Interactions
#'
#' @description Get all ceRNA interactions by given identifications (ensg_number, gene_symbol or gene_type), specific cancer type/dataset or different filter possibilities according different statistical values (e.g. FDR adjusted p-value).
#'
#' @param disease_name Name of the specific cancer type/dataset. If default is set, all available datasets with corresponding informations are shown.
#'                     Fuzzy search available.
#' @param ensg_number A vector of ensg number(s). If ensg number is set, gene symbol and gene type must be NULL. One of the three identifiers must be provided.
#' @param gene_symbol A vector of gene symbol(s). If gene symbol is set, ensg number and gene type must be NULL. One of the three identifiers must be provided.
#' @param gene_type Defines the type of gene of interest. One out of [3prime_overlapping_ncRNA, antisense, antisense_RNA, bidirectional_promoter_lncRNA, IG_C_gene, IG_C_pseudogene, IG_V_gene, IG_V_pseudogene, lincRNA, macro_lncRNA, miRNA, misc_RNA, Mt_rRNA, polymorphic_pseudogene, processed_pseudogene, processed_transcript, protein_coding, pseudogene, ribozyme, rRNA, rRNA_pseudogene, scaRNA, scRNA, sense_intronic, sense_overlapping, snoRNA, snRNA, TEC, TR_C_gene, TR_V_gene, TR_V_pseudogene, transcribed_processed_pseudogene, transcribed_unitary_pseudogene, transcribed_unprocessed_pseudogene, translated_processed_pseudogene, unitary_pseudogene, unprocessed_pseudogene, vaultRNA].
#' @param pValue Threshold of the FDR adjusted p-value. Default is 0.05.
#' @param pValueDirection Direction of the FDR adjusted p-value threshold (<, >). Must be set if pValue is set. Possible values are: "<", ">".
#' @param mscor Threshold of the 'multiple sensitivity correlation' (mscor).
#' @param mscorDirection Direction of the mscor threshold (<, >). Must be set if pValue is set. Possible values are: "<", ">".
#' @param correlation Threshold of the correlation.
#' @param correlationDirection Direction of the correlation threshold (<, >). Must be set if pValue is set. Possible values are: "<", ">".
#' @param sorting Possibilities for sorting of the results. Possible values are "pValue", "mscor" or "correlation".
#' @param descending Descending (TRUE, default) or ascending (FALSE) ordering of the results.
#' @param limit Number of results that should be shown. Default value is 100 and can be up to 1000.
#'              For more results please use batches, the provided offset parameter or download the whole dataset.
#' @param offset Starting point from where results should be shown.
#' @param information All available information about genes displayed (TRUE) or just ensg number (FALSE, default).
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
#' # Retrieve all possible ceRNAs for gene, identified by ensg_number,
#' # and threshold for pValue and mscor.
#' all_ceRNAInteractions(ensg_number=c("ENSG00000259090","ENSG00000217289"),
#'                       pValue=0.5, pValueDirection="<",
#'                       mscor=0.006, mscorDirection="<",
#'                       limit=15, information=FALSE)
all_ceRNAInteractions <- function(disease_name = NULL,
                                  ensg_number = NULL,
                                  gene_symbol = NULL,
                                  gene_type = NULL,
                                  pValue = 0.05,
                                  pValueDirection = "<",
                                  mscor = NULL,
                                  mscorDirection = "<",
                                  correlation = NULL,
                                  correlationDirection = "<",
                                  sorting = NULL,
                                  descending = TRUE,
                                  limit = 100,
                                  offset = NULL,
                                  information = FALSE){

  # all checks will be done from the API and its unit tests!

  # Base URL path
  base_url = "http://10.162.163.20:5000/ceRNAInteraction/findAll?"
  full_url = base_url

  # Create full url
  if (!is.null(disease_name)){
    full_url = paste(full_url, "disease_name=", disease_name,"&", sep="")
  }
  if (!is.null(ensg_number)){
    full_url = paste(full_url, "ensg_number=", paste(ensg_number, collapse=",", sep=""), "&", sep="")
  }
  if(!is.null(gene_symbol)) {
    full_url = paste(full_url, "gene_symbol=", paste(gene_symbol, collapse=",", sep=""), "&", sep="")
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
      mutate_at(c("correlation", "p_value", "mscor", "run.dataset.dataset_ID", "run.run_ID", "gene1.start_pos", "gene2.start_pos", "gene1.end_pos","gene2.end_pos", "gene1.chromosome_name", "gene2.chromosome_name"), as.numeric) %>%
      mutate_at(c("gene1.ensg_number", "gene1.gene_symbol", "gene2.ensg_number", "gene2.gene_symbol", "gene1.description", "gene2.description", "gene1.gene_type", "gene2.gene_type"), as.character)
    } else {
      new <- new %>%
        mutate_at(c("correlation", "p_value", "mscor", "run.dataset.dataset_ID", "run.run_ID"), as.numeric) %>%
        mutate_at(c("gene1.ensg_number", "gene1.gene_symbol", "gene2.ensg_number", "gene2.gene_symbol"), as.character)
    }
    return(new)
  }
}
