#' Join function for MG-RAST
#' 
#' Joins both the "Organism" dataset and "Ontology-KO" 
#' @param dfKEGG Dataframe containing the Ontology-KO data
#' @param dfORG Dataframe containing the Organism data
#' @export 

# library(stringr)
# library(dplyr)
# library(readr)

m5nr_join <- function(dfKEGG,dfORG){

  # rename the annotation columns
  colnames(dfKEGG)[13] <- "KEGGfunction"
  colnames(dfORG)[13] <- "organism"
  
  # keep only columns or interest
  dfKEGG <- dfKEGG[,c(1,2,13)]
  dfORG <- dfORG[,c(1,2,13)]
  
  # extract sequence id as its own column
  dfKEGG$id <- sapply(stringr::str_split(dfKEGG$`query sequence id`,"\\|"),function(x)x[2])
  dfORG$id <- sapply(stringr::str_split(dfORG$`query sequence id`,"\\|"),function(x)x[2])
  
  # join datasets by m5nr id and sequence read id
  df <- dplyr::left_join(dfKEGG,dfORG,by = c("hit m5nr id (md5sum)","id"))
  
  # only keep complete annotations
  df <- df[stats::complete.cases(df),]
  
  # extract KEGG accessions 
  df$KO <- apply(stringr::str_match(df$predicted.function,pattern = "accession\\=\\[(.*?)\\]"),1,function(x)x[2])
  
  df$organism <- stringr::str_trim(apply(stringr::str_match(df$organism,pattern = "\\[(.*?)[\\(|\\]]"),1,function(x)x[2]))
  
  
  return(df[,c(2,3,4,6,7)])
  
}