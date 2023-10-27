#TODO: Explanation comment with params etc.


get_gseaSets = function(disease_name_1, disease_name_2, condition_1, condition_2, disease_subtype_1 = NULL, disease_subtype_2 = NULL){
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
  base_url = paste(pkg.env$API.url, "/getGseaTerms", sep="")

  #Create full url
  full_url = paste(base_url, "?disease_name_1=", disease_name_1, "&disease_name_2=", disease_name_2, "&condition_1=", condition_1, "&condition_2=", condition_2, sep="")

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
    new <- do.call("data.frame", new)

    return(new)
  }
}

get_gseaSets(disease_name_1 = "liver", disease_name_2 = "thymoma", condition_1 = "disease", condition_2 = "disease")


get_gseaTerms = function(disease_name_1, disease_name_2, condition_1, condition_2, gene_set, disease_subtype_1 = NULL, disease_subtype_2 = NULL){
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
  base_url = paste(pkg.env$API.url, "/getGseaTerms", sep="")

  #Create full url
  full_url = paste(base_url, "?disease_name_1=", disease_name_1, "&disease_name_2=", disease_name_2, "&condition_1=", condition_1, "&condition_2=", condition_2, "&gene_set=", gene_set, sep="")

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
    new <- do.call("data.frame", new)

    return(new)
  }
}

get_gseaTerms(disease_name_1 = "liver", disease_name_2 = "thymoma", condition_1 = "disease", condition_2 = "disease", gene_set = "GO_Biological_Process_2023")


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
  #base_url = paste(pkg.env$API.url, "/getGseaTerms", sep="")
  base_url = paste("localhost:5000", "/gseaResults", sep="")

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

get_gseaResults(disease_name_1 = "liver", disease_name_2 = "thymoma", condition_1 = "disease", condition_2 = "disease", gene_set = "GO_Biological_Process_2023")#, term="GO:0001676")

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

img <- get_gseaPlot(disease_name_1 = "liver", disease_name_2 = "thymoma", condition_1 = "disease", condition_2 = "disease", gene_set = "GO_Biological_Process_2023", term="GO:0001676")

