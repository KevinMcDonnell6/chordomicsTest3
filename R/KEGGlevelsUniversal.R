#' Access KEGG Orthology
#' 
#' Gets functional information for MPA output csv using the KEGG Orthology API.
#' @param df Output csv from MPA
#' @export



#########################################################################
# Function gets functional information for MPA output from KEGG Orthology
#########################################################################
#library(httr)
#library(stringr)

KOLevels <- function(df,KOcolumnName,ECcolumnName=NULL){
  
  df <- data.frame(lapply(df, as.character), stringsAsFactors=FALSE)
  df$KO.Level1 <- NA
  df$predicted.function <- NA
  i<-0
  pb<- txtProgressBar(min = 0, max = length(df[,1]), style = 3)  
  #withProgress(message = 'Acquiring Data', value = 0, {
  n<-nrow(df)
  # for each entry send GET request
  for(k in paste(df[,KOcolumnName])){
    # View(test)
    i<-i+1
    setTxtProgressBar(pb, i)
    #incProgress(1/n,detail = paste("Part",i,"of",n))
    
    # if entry isnt empty
    if(k != "|" & k!="" & k!="NA"){
      
      # join KO ids for request
      k<-stringr::str_replace_all(paste(k),"\\|","+")
      
      # set url
      urlk<- "http://rest.kegg.jp"
      
      # set KOs for request
      pathk <- paste("get/ko:",k,sep = "")
      
      # send request
      kegg.raw.result <- httr::GET(url = urlk, path = pathk)
      
      # covert to characters
      kegg.this.raw.content <- rawToChar(kegg.raw.result$content)
      
      if(kegg.this.raw.content!=""){
        
        
        # extract BRITE section from returned data
        sent<- kegg.this.raw.content
        sent<- strsplit(sent,"ENTRY")
        sent_split<- stringr::str_extract_all(sent[[1]][2],"KEGG Orthology \\(KO\\) \\[BR\\:ko00001\\]\n([\\S\\s]*)")
        
        # extract levels from BRITE hierarchy
        names_ <-  stringr::str_split(sent_split,pattern = "\n")[[1]][1:4]
        
        # Sometimes results include extra level and so we need to go one level deeper (levels 3 and 4)
        if(stringr::str_match(names_[2],pattern = "\\d{5} *(.*)")[2] == "Brite Hierarchies" | stringr::str_match(names_[2],pattern = "\\d{5} *(.*)")[2] == "Not Included in Pathway or Brite"){
          df$KO.Level1[i] <- stringr::str_match(names_[3],pattern = "\\d{5} *(.*)")[2]
          df$predicted.function[i] <- stringr::str_match(names_[4],pattern = "\\d{5} *(.*)")[2]
          
        }
        
        # extract levels 2 and 3
        else{
          df$KO.Level1[i] <- stringr::str_match(names_[2],pattern = "\\d{5} *(.*)")[2]
          df$predicted.function[i] <- stringr::str_match(names_[3],pattern = "\\d{5} *(.*)")[2] 
          # print(df$KO.Level1[i])
        }
      }
    }
    
    # If KO is missing search for KO using EC number
    if(!is.null(ECcolumnName)){
    
      if(k == "|" & df[,ECcolumnName][i] != "|" & df[,ECcolumnName][i] != "UNKNOWN" & df[,ECcolumnName][i] != "|"){
        
          ec <- df[,ECcolumnName][i]
          ec <- stringr::str_split(ec,"\\|")[[1]][1]
          ec <- stringr::str_replace_all(ec,"\\.\\-","")
          
          # set url
          urlec<- "http://rest.kegg.jp"
          
          # set KOs for request
          pathec <- paste("find/ko/ec:",ec,sep = "")
          
          # send request
          ec.raw.result <- httr::GET(url = urlec, path = pathec)
          
          # covert to characters
          ec.this.raw.content <- rawToChar(ec.raw.result$content)
          
          k <- stringr::str_extract(ec.this.raw.content,"K\\d{5}")
          
          # repeat as above
          urlk<- "http://rest.kegg.jp"
          
          # set KOs for request
          pathk <- paste("get/ko:",k,sep = "")
          
          # send request
          kegg.raw.result <- httr::GET(url = urlk, path = pathk)
          
          # covert to characters
          kegg.this.raw.content <- rawToChar(kegg.raw.result$content)
          
          # extract BRITE section from returned data
          sent<- kegg.this.raw.content
          sent<-strsplit(sent,"ENTRY")
          sent_split<- stringr::str_extract_all(sent[[1]][2],"KEGG Orthology \\(KO\\) \\[BR\\:ko00001\\]\n([\\S\\s]*)")
          
          # extract levels from BRITE hierarchy
          names_ <-  stringr::str_split(sent_split,pattern = "\n")[[1]][1:4]
          
          # Sometimes results include extra level and so we need to go one level deeper (levels 3 and 4)
          if(!is.na(stringr::str_match(names_[2],pattern = "\\d{5} *(.*)")[2])){
            if(stringr::str_match(names_[2],pattern = "\\d{5} *(.*)")[2]=="Brite Hierarchies" | stringr::str_match(names_[2],pattern = "\\d{5} *(.*)")[2] == "Not Included in Pathway or Brite"){
              df$KO.Level1[i] <- stringr::str_match(names_[3],pattern = "\\d{5} *(.*)")[2]
              df$predicted.function[i] <- stringr::str_match(names_[4],pattern = "\\d{5} *(.*)")[2]
              
            }
            
            else{
              df$KO.Level1[i] <- stringr::str_match(names_[2],pattern = "\\d{5} *(.*)")[2]
              df$predicted.function[i] <- stringr::str_match(names_[3],pattern = "\\d{5} *(.*)")[2]  
            }
          }
        }
  
    }
  
  }
  
  df$KO.Level1 <- stringr::str_replace_all(tolower(df$KO.Level1),"protein families: ",replace="")
  df$KO.Level1 <- stringr::str_replace_all(tolower(df$KO.Level1),"unclassified: ",replace="")
  return(df)
}

# test <- read.csv("FAfoodWaste/Day1.csv")[1:10,]
# test <- read_tsv("new.faa.emapper.annotations",col_names = F)[1:10,]
# test <- read_tsv("new.faa.emapper.annotations",col_names = F)[1:10,]
# test <- read.csv("M5NRoutput(mgm4827488).csv")[1:10,]
# View(test)
# test <- KOLevels(test,KOcolumnName = "Meta.Protein.KO",ECcolumnName = "Meta.Protein.EC")
# test <- KOLevels(test,KOcolumnName = "X7")
# test <- KOLevels(test,KOcolumnName = "KO")
