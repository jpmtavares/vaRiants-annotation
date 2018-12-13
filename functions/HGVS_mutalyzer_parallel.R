##########################################################
#                                                        #
#                      Mutalyzer                         #
#                                                        #
##########################################################

################################################################
#             get HGVS annotation - main function              # 
################################################################
hgvs<-function(variants){
  ##__________________________________
  ## get variantPosition input (g.)
  ##__________________________________
    var<-ifelse(
      ##DELECTION
      nchar(as.character(variants$Ref))>1,
      paste(variants$Chr,":","g.",variants$Position,"del",sub('.','',variants$Ref),sep=""),
      ##INSERTION
      ifelse(nchar(as.character(variants$Alt))>1,
             paste(variants$Chr,":","g.",as.numeric(as.character(variants$Position))-1,"_",variants$Position,"ins",sub('.','',variants$Alt),sep=""),
             ##SNP
             paste(variants$Chr,":","g.",
                   variants$Position,variants$Ref,
                   ">",variants$Alt,sep=""))
    )
    ##__________________________________
    ## parallelize
    ##__________________________________
    registerDoParallel(cores=6)
    parallel<-foreach(x=var, y=variants$rs_ID, z=variants$Rank.Exons.Introns) %dopar% mutalyzer(z,x,y)
    
    ##__________________________________
    ## get refSeq_mRNA
    ##__________________________________
    NM<-lapply(parallel, function(x){
      filter(x, str_detect(V1, paste(unique(variants$refSeq_mRNA),collapse="|")) == TRUE) %>%
        filter(row_number()==n())
    })
    ##__________________________________
    ## get refSeq_protein
    ##__________________________________
    NP<-lapply(parallel, function(x){
      filter(x, str_detect(V1, paste(unique(variants$refSeq_protein),collapse="|")) == TRUE) %>%
        filter(row_number()==n())
    })
    ##__________________________________
    ## write output
    ##__________________________________
    hgvs<-data.frame(variants[!(names(variants) %in% c("refSeq_mRNA","refSeq_protein"))],
                     RefSeq_mRNA=do.call(rbind.data.frame,is.empty(NM))[,1],
                     HGVS_c=do.call(rbind.data.frame,is.empty(NM))[,2],
                     RefSeq_protein=do.call(rbind.data.frame,is.empty(NP))[,1],
                     HGVS_p=do.call(rbind.data.frame,is.empty(NP))[,2])
  return(hgvs)
}

################################################################
# Make a GET request to get snp HGVS annotation from mutalyzer #
################################################################
mutalyzer<-function(exonic,variantPosition,rs_ID){
  requests<-import("requests")
  #For exonic variants, look for refSeq_protein (p.) with rs_ID
  if(grepl("E",as.character(exonic))){
    parameters<-dict(rs_id=as.character(rs_ID))
    response<-requests$get("https://mutalyzer.nl/json/getdbSNPDescriptions?", params = parameters)
    results<-bind_cols(as.data.frame(str_split_fixed(response$json(),"\\:", 2)))
    
    ## If "Non existing rs_ID in the DB or no root element"
    if(any(grepl("Non existing|Incorrect",results[,1]))){
      parameters<-dict(build="hg19",variant=as.character(variantPosition))
      response<-requests$get("https://mutalyzer.nl/json/numberConversion?", params = parameters)
      results<-bind_cols(as.data.frame(str_split_fixed(response$json(),"\\:", 2)))
    }
  }else{ #look for refSeq_mRNA (c.) with variant Position (it's faster)
    parameters<-dict(build="hg19",variant=as.character(variantPosition))
    response<-requests$get("https://mutalyzer.nl/json/numberConversion?", params = parameters)
    results<-bind_cols(as.data.frame(str_split_fixed(response$json(),"\\:", 2)))
  }
  return(results)
}

################################################################
#              is some element of a list, empty?               #
################################################################
is.empty<-function(list){
  lapply(list,function(x){
    if(nrow(x)==0){
      x<-data.frame(V1=NA,V2=NA)
    }
    return(x)
  })
}
