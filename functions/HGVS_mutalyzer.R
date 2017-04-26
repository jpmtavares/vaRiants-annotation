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
      doc.text <- try(read_html(paste("https://mutalyzer.nl/Position-converter?assembly_name_or_alias=GRCh37&description=",Chr,"%3Ag.",as.numeric(as.character(Position)),"_",as.numeric(as.character(Position)),"ins",sub('.', '', Alt),sep="")) %>%
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
        doc.text <- try(read_html(paste("https://mutalyzer.nl/Position-converter?assembly_name_or_alias=GRCh37&description=",Chr,"%3Ag.",as.numeric(as.numeric(as.character(Position))-(nchar(Ref)-2)),"_",as.numeric(as.character(Position)),"del",sub('.', '', Ref),sep="")) %>%
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
        doc.text <- try(read_html(paste("https://mutalyzer.nl/Position-converter?assembly_name_or_alias=GRCh37&description=",Chr,"%3Ag.",as.numeric(as.character(Position)),Ref,"%3E",Alt,sep="")) %>%
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
    var<-c(paste("ins",sub('.', '', Alt),sep=""), paste("ins",reverseDNA(sub('.', '', Alt),complement = T),sep=""))
  }else{
    ## DELETION
    if(nchar(Ref)>1){
      var<-c(paste("del",sub('.', '', Ref),sep=""), paste("del", reverseDNA(sub('.', '', Ref),complement = T),sep=""))
    }else{
      ## SNP
      var<-c(paste(Ref,">",Alt,sep=""),
             reverseDNA(paste(Ref,">",Alt,sep="")),
             reverseDNA(paste(Ref,">",Alt,sep=""),complement = T))
    }
  }
  ##########################################################
  #                    return output                       #
  ##########################################################
  if(is.null(hgvs_protein)){
    hgvs_protein<-matrix(".",nrow = 6, ncol = 2)
  }
  if(is.null(hgvs_mRna)){
    hgvs_mRna<-matrix(".",nrow = 6, ncol = 2)
    return(c(".",".",".","."))
  }else{
    return(c(#Chr,Position,rs_ID,Ref,Alt,
      rbind(hgvs_mRna[max(grep(paste(var,collapse="|"),hgvs_mRna[,2])),]),
      rbind(hgvs_protein[min(grep(paste(var,collapse="|"),hgvs_mRna[,2])),])))
  }
}

##########################################################
#        apply mutalyzer function to data.frame          #
##########################################################
hgvs<-function(variants){
  mapply(mutalyzer,variants[,"Chr"],variants[,"Position"],variants[,"rs_ID"],
         variants[,"Ref"],variants[,"Alt"],variants[,"refSeq_mRNA"],
         variants[,"refSeq_protein"])%>%
    t()%>%
    cbind(variants,RefSeq_mRNA=.[,1],HGVS_c=.[,2],
          RefSeq_protein=.[,3],HGVS_p=.[,4]) %>%
    select(-refSeq_mRNA,-refSeq_protein,-`1`,-`2`,-`3`,-`4`)
}
