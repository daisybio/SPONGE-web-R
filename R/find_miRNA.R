#' Find Sponged miRNAs
#'
#' @description Get all miRNAs that contribute to all interactions between the given identifiers (ensg_number or gene_symbol).
#'
#' @param disease_name Name of the specific cancer type/dataset. If default is set, all available datasets with corresponding informations are shown.
#'                     Fuzzy search available.
#' @param ensg_number A vector of ensg number(s). If ensg number is set, gene symbol and gene type must be NULL. One of the three identifiers must be provided.
#' @param gene_symbol A vector of gene symbol(s). If gene symbol is set, ensg number and gene type must be NULL. One of the three identifiers must be provided.
#' @param gene_type Defines the type of gene of interest. One out of [3prime_overlapping_ncRNA, antisense, antisense_RNA, bidirectional_promoter_lncRNA, IG_C_gene, IG_C_pseudogene, IG_V_gene, IG_V_pseudogene, lincRNA, macro_lncRNA, miRNA, misc_RNA, Mt_rRNA, polymorphic_pseudogene, processed_pseudogene, processed_transcript, protein_coding, pseudogene, ribozyme, rRNA, rRNA_pseudogene, scaRNA, scRNA, sense_intronic, sense_overlapping, snoRNA, snRNA, TEC, TR_C_gene, TR_V_gene, TR_V_pseudogene, transcribed_processed_pseudogene, transcribed_unitary_pseudogene, transcribed_unprocessed_pseudogene, translated_processed_pseudogene, unitary_pseudogene, unprocessed_pseudogene, vaultRNA].
#'
#' @return A data_frame containing all found miRNAs.
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom dplyr mutate_at %>%
#' @importFrom httr GET content headers
#'
#' @examples
#' sponge_miRNA(disease_name="kidney", gene_symbol = c("TCF7L1", "SEMA4B"))
sponge_miRNA <- function(disease_name,
                       ensg_number = NULL,
                       gene_symbol = NULL,
                       gene_type = NULL){

  # all checks will be done from the API and its unit tests!

  # Base URL path
  base_url = "http://10.162.163.20:5000/miRNAInteraction/findceRNA?"
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
      mutate_at(c("interactions_genegene.run.run_ID"), as.numeric) %>%
      mutate_at(c("mirna.hs_nr", "mirna.mir_ID"), as.character)
    return(new)
  }
}
