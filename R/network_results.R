library("httr")
library("jsonlite")
library("utils")
library("ggplot2")
library(data.table)

#' Scores and coordinates within a 2D MDS map between Cancer Type network d other Cancer Type or Subtype networks.
#'
#' @param disease_name The name of the Cancer Type of interest as string. If default is set, only cancer types.
#' @param level Either "gene" or "transcript", which was introduced in SPONGEdb version 2. The default is "gene".
#'
#' @return list containing pairwise scores and 2D coordinates for all Cancer Types and, if available, the Types Subtypes.
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
#' @importFrom httr GET content headers
#'
#' @examples
#' # Get gene level scores and 2D coordinates for all Cancer Types and all Subtypes of "Breast invasive carcinoma".
#' get_network_results(disease_name = "Breast invasive carcinoma", level = "gene")
get_network_results <- function(disease_name = "Pan-cancer", level = "gene"){
  # all checks will be done from the API and its unit tests!
  # Base URL path
  base_url = paste(pkg.env$API.url, "/networkResults?", sep="")
  full_url = base_url
  # Create full url
  full_url = paste(full_url, "disease_name=", disease_name, "&", sep="")
  full_url = paste(full_url, "level=", level, "&", sep="")
  # Encode the URL with characters for each space.
  full_url <- URLencode(full_url)
  message(full_url)
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
    type_scores <- data.frame(new$type$scores$values)
    colnames(type_scores) <- new$type$scores$labels
    rownames(type_scores) <- new$type$scores$labels
    type_distances <- data.frame(x = new$type$euclidean_distances$x,
                                 y = new$type$euclidean_distances$y,
                                 labels = new$type$euclidean_distances$labels)
    type = list(scores=type_scores, distances=type_distances)
    if (!is.null(new$subtype$scores)){
      subtype_scores <- data.frame(new$subtype$scores$values)
      colnames(subtype_scores) <- new$subtype$scores$labels
      rownames(subtype_scores) <- new$subtype$scores$labels
      subtype_distances <- data.frame(x = new$subtype$euclidean_distances$x,
                                      y = new$subtype$euclidean_distances$y,
                                      labels = new$subtype$euclidean_distances$labels)
      subtype = list(scores=subtype_scores, distances=subtype_distances)
    } else {
      subtype <- list()
    }
    return(list(type=type, subtype=subtype))
  }
}

#' Shows a Heat map of hubness scores
#'
#' @param df A data.frame with hubness values returned from the get_network_results function.
#' @param title Allows the user to specify a plot title.
#' @param tri Toggles triangle shape for enhanced visibility of differences. Default is TRUE.
#'
#' @return
#' @export
#'
#' @import ggplot2
#' @importFrom data.table melt
#'
#' @examples
#' # Visualize gene level scores for "Breast invasive carcinoma" and  all its Subtypes.
#' res <- get_network_results(disease_name = "Breast invasive carcinoma", level = "gene")
#' plot_heatmap(res$subtype$scores)
plot_heatmap <- function(df, title = "Heatmap", tri = TRUE){
  if (tri) {
    df[upper.tri(df, diag = TRUE)] <- NA
  }
  m <- melt(df)
  df <- cbind(m, c(colnames(df)))
  colnames(df) <- c("data1", "value", "data2")
  df$data2 <- factor(df$data2, levels = unique(df$data2))
  ggplot(na.omit(df), aes(x = data1, y = data2, fill = value)) +
    geom_tile() +
    labs(title=title, x="", y="") +
    scale_x_discrete(position = "top") +
    scale_fill_gradient(low = munsell::mnsl("5P 7/12"), high = munsell::mnsl("5P 2/12"), space = 'Lab') +
    theme(axis.ticks=element_blank(),
          axis.line=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_line(color='#eeeeee'))
}

#' Shows a MDS plot of euclidean distances
#'
#' @param df A data.frame with cmdscaled values returned from the get_network_results function.
#' @param title Allows the user to specify a plot title.
#'
#' @return
#' @export
#'
#' @import ggplot2
#'
#' @examples
#' # Visualize gene level distances for "Breast invasive carcinoma" and  all its Subtypes.
#' res <- get_network_results(disease_name = "Breast invasive carcinoma", level = "gene")
#' plot_MDS(res$subtype$distances)
plot_MDS <- function(df, title = "MDS plot"){
  ggplot(df, aes(x = x, y = y)) +
    geom_point() +
    labs(title=title, x="", y="") +
    geom_text(aes(label=labels), size=3, nudge_y=1) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}
