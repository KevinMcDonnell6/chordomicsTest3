#' Function for use with eggNOG emapper
#' 
#' Extracts KO accession and gene name from emapper.annotations file. Appends phylogeny to dataframe
#' @param df Dataframe downloaded from eggNOG emapper
#' @export



#library(httr)
# library(stringr)
# Name of eggnog file goes in these quotes eg. "test.annotations"
#               |
#               |
#               |
#              \ /
#               V
#df <- read_tsv("new.faa.emapper.annotations",col_names = F)


eggnog_tax <- function(df){
  
  # find column matching gene name
  index <- (1:ncol(df))[as.vector(!is.na(stringr::str_match(df[1,],"\\d{6}\\.")))]
  
  KOindex <- 7
  df <- df[,c(index,KOindex)]
  
  df$taxonomy_accession <- sapply(stringr::str_split(df[,1],"\\."),function(x)x[2])
  # df$taxonomy_accession <- sapply(str_split(df$X2,"\\."),function(x)x[2])
  
 
  len <- nrow(df)
  
  
  df_tax <- data.frame(SuperKingdom="",Kingdom="",Phylum="",
                       Class="",Order="",
                       Family="",Genus="",
                       Species="",stringsAsFactors = F)
  
  pb<- txtProgressBar(min = 0, max = len, style = 3)
  
  
  for(i in 1:len){
    setTxtProgressBar(pb, i)
    # print(i)
    LCA <- df$taxonomy_accession[i]
    # LCA <- str_replace_all(df$organism[i],pattern = " ",replacement = "%20")
    url  <- "https://www.uniprot.org"
    columns <- paste("id","entry%20name","lineage(SUPERKINGDOM)","lineage(KINGDOM)","lineage(PHYLUM)",
                     "lineage(CLASS)","lineage(FAMILY)","lineage(GENUS)","lineage(ORDER)","lineage(SPECIES)",sep=",")
    
    path <- paste("uniprot/?query=",LCA,"&sort=score&columns=",columns,"&limit=1&format=tab",sep = "")
    
    
    # send request
    uni.raw.result <- httr::GET(url = url, path = path)
    uni.raw.result
    
    if (uni.raw.result$status_code == 400 | uni.raw.result$status_code == 404 | rawToChar(uni.raw.result$content)==""){
      df_tax[i,"SuperKingdom"] <- NA
      df_tax[i,"Kingdom"] <- NA
      df_tax[i,"Phylum"] <- NA
      df_tax[i,"Class"] <- NA
      df_tax[i,"Order"] <- NA
      df_tax[i,"Family"] <- NA
      df_tax[i,"Genus"] <- NA
      df_tax[i,"Species"] <- NA
    }
    
    else{
      # convert returned data to charachters
      uni.this.raw.content <- rawToChar(uni.raw.result$content)
      uni.this.raw.content
      # convert data to dataframe
      res <- read.delim(textConnection(uni.this.raw.content),sep="\t",stringsAsFactors = F)
      
      df_tax[i,"SuperKingdom"] <- res$Taxonomic.lineage..SUPERKINGDOM.
      # Uniprot holds bacteria as Super kingdom not kingdom
      df_tax[i,"Kingdom"] <- ifelse(res$Taxonomic.lineage..SUPERKINGDOM.=="Bacteria" & is.na(res$Taxonomic.lineage..KINGDOM.),"Bacteria",res$Taxonomic.lineage..KINGDOM.)
      df_tax[i,"Phylum"] <- res$Taxonomic.lineage..PHYLUM.
      df_tax[i,"Class"] <- res$Taxonomic.lineage..CLASS.
      df_tax[i,"Order"] <- res$Taxonomic.lineage..ORDER.
      df_tax[i,"Family"] <- res$Taxonomic.lineage..FAMILY.
      df_tax[i,"Genus"] <- res$Taxonomic.lineage..GENUS.
      df_tax[i,"Species"] <- res$Taxonomic.lineage..SPECIES.
    }
    
    
  }
  
  df <- cbind(df,df_tax)
  
  return(df)
}

#d2 <- eggnog_tax(df)
#View(d2)

