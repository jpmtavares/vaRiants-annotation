##########################################################
#                                                        #
#             HGVS notation (with mutalyzer)             #
#                                                        #
##########################################################

##########################################################
#        access https://mutalyzer.nl to get HGVS         #
##########################################################
mutalyzer<-function(Chr,Position,rs_ID,Ref,Alt,refSeq_mRNA,refSeq_protein){
  ##########################################################
  #                     known variant                      #
  ##########################################################
  if(rs_ID!="."){
    doc.text <- try(read_html(paste("https://mutalyzer.nl/snp-converter?rs_id=",rs_ID,sep="")) %>%
                      html_nodes("p code") %>% 
                      html_text()
    )
    hgvs_mRna<-do.call('rbind',strsplit(doc.text[grepl(refSeq_mRNA,doc.text)],':',fixed=TRUE))
    hgvs_protein<-do.call('rbind',strsplit(doc.text[grepl(refSeq_protein,doc.text)],':',fixed=TRUE))
  }else{
    ##########################################################
    #                    new variant                        #
    ##########################################################
    ## INSERTION
    if(nchar(Alt)>1){
      #Read and parse html file
      doc.text <- try(read_html(paste("https://mutalyzer.nl/position-converter?assembly_name_or_alias=GRCh37&description=",Chr,"%3Ag.",as.numeric(as.character(Position)),"_",as.numeric(as.character(Position)),"ins",sub('.', '', Alt),sep="")) %>%
                        html_nodes("pre") %>% 
                        html_text()
      )
      #remove weird characters
      doc.text<-do.call('rbind',strsplit(doc.text,"\n\t\t",fixed=T))
      doc.text<-gsub(".*:\t","",doc.text)
      hgvs_mRna<-do.call('rbind',strsplit(doc.text[grepl(refSeq_mRNA,doc.text)],':',fixed=TRUE))
      hgvs_protein<-do.call('rbind',strsplit(doc.text[grepl(refSeq_protein,doc.text)],':',fixed=TRUE))
    }else{
      ## DELETION
      if(nchar(Ref)>1){
        #Read and parse html file
        doc.text <- try(read_html(paste("https://mutalyzer.nl/position-converter?assembly_name_or_alias=GRCh37&description=",Chr,"%3Ag.",as.numeric(as.numeric(as.character(Position))-(nchar(Ref)-2)),"_",as.numeric(as.character(Position)),"del",sub('.', '', Ref),sep="")) %>%
                          html_nodes("pre") %>% 
                          html_text()
        )
        #remove weird characters
        doc.text<-do.call('rbind',strsplit(doc.text,"\n\t\t",fixed=T))
        doc.text<-gsub(".*:\t","",doc.text)
        hgvs_mRna<-do.call('rbind',strsplit(doc.text[grepl(refSeq_mRNA,doc.text)],':',fixed=TRUE))
        hgvs_protein<-do.call('rbind',strsplit(doc.text[grepl(refSeq_protein,doc.text)],':',fixed=TRUE))
      }else{
        ## SNP
        #Read and parse html file
        doc.text <- try(read_html(paste("https://mutalyzer.nl/position-converter?assembly_name_or_alias=GRCh37&description=",Chr,"%3Ag.",as.numeric(as.character(Position)),Ref,"%3E",Alt,sep="")) %>%
                          html_nodes("pre") %>% 
                          html_text()
        )
        #remove weird characters
        doc.text<-do.call('rbind',strsplit(doc.text,"\n\t\t",fixed=T))
        doc.text<-gsub(".*:\t","",doc.text)
        hgvs_mRna<-do.call('rbind',strsplit(doc.text[grepl(refSeq_mRNA,doc.text)],':',fixed=TRUE))
        hgvs_protein<-do.call('rbind',strsplit(doc.text[grepl(refSeq_protein,doc.text)],':',fixed=TRUE))
      }
    }
  }
  ##########################################################
  #                   get correct output                   #
  ##########################################################
  ## INSERTION
  if(nchar(Alt)>1){
    var<-c(paste("ins",sub('.', '', Alt),sep=""),
           paste("ins",paste(rev(strsplit(sub('.', '', Alt),"")[[1]]),collapse = ""),sep=""),
           paste("ins",reverseDNA(sub('.', '', Alt)),sep=""), 
           paste("ins",reverseDNA(sub('.', '', Alt),complement = T),sep=""))
  }else{
    ## DELETION
    if(nchar(Ref)>1){
      var<-c(paste("del",sub('.', '', Ref),sep=""),
             paste("del",paste(rev(strsplit(sub('.', '', Ref),"")[[1]]),collapse = ""),sep=""),
             paste("del", reverseDNA(sub('.', '', Ref)),sep=""), 
             paste("del", reverseDNA(sub('.', '', Ref),complement = T),sep=""))
    }else{
      ## SNP
      var<-c(paste(Ref,">",Alt,sep=""),
             paste(Alt,">",Ref,sep=""),
             reverseDNA(paste(Ref,">",Alt,sep="")),
             reverseDNA(paste(Ref,">",Alt,sep=""),complement = T))
    }
  }
  ##########################################################
  #                    return output                       #
  ##########################################################
  if(is.null(hgvs_protein)){
    hgvs_protein<-matrix(".",nrow = 10, ncol = 2)
  }
  if(is.null(hgvs_mRna)){
    hgvs_mRna<-matrix(".",nrow = 10, ncol = 2)
    return(c(".",".",".","."))
  }else{
    return(c(#Chr,Position,rs_ID,Ref,Alt,
      rbind(hgvs_mRna[max(grep(paste(var,collapse="|"),hgvs_mRna[,2])),]),
      rbind(hgvs_protein[max(grep(paste(refSeq_protein,collapse="|"),hgvs_protein[,1])),])))
  }
}

##########################################################
#        apply mutalyzer function to data.frame          #
##########################################################
hgvs<-function(variants){
  mapply(mutalyzer,variants[,"Chr"],variants[,"Position"],variants[,"rs_ID"],
         variants[,"Ref"],variants[,"Alt"],variants[,"refSeq_mRNA"],
         variants[,"refSeq_protein"])%>% #applying mutalyzer function to a data.frame
    t()%>% #transpose data.frame
    cbind(variants,RefSeq_mRNA=.[,1],HGVS_c=.[,2],
          RefSeq_protein=.[,3],HGVS_p=.[,4]) %>% #add new columns
    select(-refSeq_mRNA,-refSeq_protein,-`1`,-`2`,-`3`,-`4`) #remove duplicated columns
}

##########################################################
#                    variant type                        #
##########################################################
variant_type<-function(HGVS_c, HGVS_p){
  ##__________________________________
  ## if variant is exonic
  ##__________________________________
  ifelse(str_detect(HGVS_c,"[\\*\\+\\-]")==FALSE, # if variant is exonic
         ifelse(grepl("=",HGVS_p)==TRUE | str_extract(HGVS_p,"[^p.](?:(?!\\d).)*") == str_extract(HGVS_p,".{3}$"), "silent", 
                ifelse(grepl("Ter",HGVS_p)==TRUE, "nonsense",
                       ifelse(grepl("fs",HGVS_p)==TRUE, "frame-shift", "missense"))),
         ##__________________________________
         ## if variant is intronic
         ##__________________________________
         ifelse(grepl("\\*",HGVS_c)==TRUE, "3-UTR",
                ifelse(grepl("\\.\\-",HGVS_c)==TRUE, "5-UTR",
                       ifelse(grepl("[0-9][-+][0-9][diATGC]|[0-9][-+][1][0-9][diATGC]",HGVS_c)==TRUE, "splicing", "intronic"))))
}

###########################################################
# correct HGVS protein annotation for synonymous variants #
###########################################################
synonymous<-function(HGVS_p){
  ifelse(grepl("=",HGVS_p)==TRUE,
         gsub("=",str_extract(HGVS_p,"[^p.](?:(?!\\d).)*"),HGVS_p),
         as.character(HGVS_p))
}