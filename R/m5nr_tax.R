#' Function to access Taxonomy
#'
#' Get phlyogeny from a given species through the UniProt API.
#' @param df Dataframe containing a column named "organism"

# library(stringr)
# library(dplyr)
# library(readr)
# library(httr)
# library(jsonlite)


m5nr_tax <- function(df){
  
        len <- nrow(df)
        
        # df_tax <- data.frame(SuperKingdom=character(),Kingdom=character(),Phylum=character(),
        #                      Class=character(),Order=character(),
        #                      Family=character(),Genus=character(),
        #                      Species=character(),stringsAsFactors = F)
        df_tax <- data.frame(SuperKingdom="",Kingdom="",Phylum="",
                             Class="",Order="",
                             Family="",Genus="",
                             Species="",stringsAsFactors = F)
        
        pb<- txtProgressBar(min = 0, max = len, style = 3)
        
        
        for(i in 1:len){
          setTxtProgressBar(pb, i)
          # LCA <- df$Lowest.Common.Ancestor[i]
          LCA <- stringr::str_replace_all(df$organism[i],pattern = " ",replacement = "%20")
          url  <- "https://www.uniprot.org"
          columns <- paste("id","entry%20name","lineage(SUPERKINGDOM)","lineage(KINGDOM)","lineage(PHYLUM)",
                           "lineage(CLASS)","lineage(FAMILY)","lineage(GENUS)","lineage(ORDER)","lineage(SPECIES)",sep=",")
          
          path <- paste("uniprot/?query=",LCA,"&sort=score&columns=",columns,"&limit=1&format=tab",sep = "")
          
          
          # send request
          uni.raw.result <- httr::GET(url = url, path = path)
          uni.raw.result
          
          if (uni.raw.result$status_code == 400){
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
          
          # start_j <- str_which(as.character(df_tax[i,]),as.character(df$organism[i]))[1]#Lowest.Common.Ancestor[i]))[1]
          # if(!is.na(start_j) & start_j!=7){
          #   for(j in (start_j+1):7){df_tax[i,j]<- "Higher Taxa"}
          # }
        }

        return(df_tax)
        }

