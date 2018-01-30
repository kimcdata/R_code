genpept <- list()

j=1

  first <- gi[,j]
  first[which(first>0)] -> first
  
  curl = getCurlHandle()
  pages = list()
  for(u in 1:length(first)) {
     first[u] -> tar
     pages[[u]] = getURL(paste("http://www.ncbi.nlm.nih.gov/protein/",tar,sep=""), curl = curl)
     print(u)
  }

  genpept[[j]] <- c(1)

  print("done")

    for(i in 1:length(first)){
      acc <- c()
      print(i)
      tar <- first[i]
      tar <- sub("gi","",tar)

      url <- getURL(paste("http://www.ncbi.nlm.nih.gov/protein/",tar,sep=""))
      url <- unlist(strsplit(url, "\n"))

        ## IF SEQUENCE IS OUTDATED
        
        if(length(grep("Record removed", url))>0){
        genpept[[j]][i] <- ""
        } else {

        if(length(url[grep("This sequence has been replaced by", url)])>0){
           url[grep("This sequence has been replaced by", url)] -> tmp
           tar <- unlist(strsplit(unlist(strsplit(tmp, ">")),"<"))[5]
           url <- getURL(paste("http://www.ncbi.nlm.nih.gov/protein/",tar,sep=""))
           url <- unlist(strsplit(url, "\n"))
        }

        ## all sequences
        
        res <- unlist(strsplit(unlist(strsplit(url[grep("<p class=\"itemid\">", url)], ": ")), "<"))[3]
        #unlist(strsplit(unlist(strsplit(res," "))[5], "<"))[1]  -> res
        sub("\\.[0-9]","",res) -> acc
        genpept[[j]][i] <- acc

        }
        }


#        ### GenBank sequence
#
#        if(length(url[grep("<p class=\"itemid\">GenBank: ", url)])>0){
#          url[grep("<p class=\"itemid\">GenBank: ", url)] -> res
#        }
#
#        ## NCBI Ref Sequence
#
#        else if(length(url[grep("<p class=\"itemid\">NCBI Reference Sequence: ", url)])>0){
#          url[grep("NCBI Reference Sequence", url)] -> res
#        }
#
#        ##
#        else if(length(url[grep("<p class=\"itemid\">PDB: ", url)])>0){
#          url[grep("<p class=\"itemid\">PDB:", url)] -> res
#        }
#
#        ## REMOVE VERSION NUMBER
#
#        if(grep("\\.", res)>0){
#              acc <- unlist(strsplit(unlist(strsplit(res, ": ")), "\\."))[2]
#              genpept[[j]][i] <- acc
#        }


#    }

#}