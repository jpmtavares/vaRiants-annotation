##########################################################
#                                                        #
#                     MutationTaster                     #
#                                                        #
##########################################################

MT<-function(variants){
  ##########################################################
  #                        INDELS                          #
  ##########################################################
  if(nchar(variants["Ref"]) > 1 | nchar(variants["Alt"]) > 1){
    ## get DNA snippet
    snippet<-paste(getSeq(Hsapiens,variants["Chr"],
                          as.numeric(as.character(variants["Position"]))-20,as.numeric(as.character(variants["Position"]))-1),
                   paste("[",as.character(variants["Ref"]),"/",as.character(variants["Alt"]),"]",sep=""),
                   getSeq(Hsapiens,as.character(variants["Chr"]),as.numeric(as.character(variants["Position"]))+1,as.numeric(as.character(variants["Position"]))+20),
                   sep='')
    
    ## get reverse complement DNA snnipet if gene is located in negative strand
    if((variants["Strand"]=="-")&&(!is.na(variants["Strand"]))){
      snippet<-reverseDNA(snippet, complement = T)
    }
    
    # Read and parse HTML file
    session<-html_session("http://www.mutationtaster.org/")
    form<-set_values(html_form(session)[[1]],
                     transcript_stable_id_text=variants["ENSTranscript"],
                     sequence_type="gDNA",
                     sequence_snippet=snippet)
    doc.text <- submit_form(session, form) %>%
      read_html() %>% 
      html_nodes("h3") %>% 
      html_text()
    if(length(doc.text)==0){
      doc.text <- submit_form(session, form) %>%
        read_html() %>% 
        html_nodes("h2") %>% 
        html_text()
    }
    #read HTML tables
    tables <- submit_form(session, form) %>%
      read_html() %>% html_table(fill=T)
    #get precition result
    toMatch <- c("polymorphism","disease causing","wrong input format","data problem","annotation problem")
    MutationTaster<-unique(grep(paste(toMatch,collapse="|"), 
                                doc.text, value=TRUE))
    #get homozygous counts in 1000G and ExAC
    homozygous<-data.frame(homozygous_1000G=NA,homozygous_ExAC=NA)
    if((!is.null(MutationTaster))&&(!is.na(MutationTaster))){
      if(MutationTaster!=("data problem")&&MutationTaster!=("wrong input format")&&MutationTaster!=("annotation problem")){
        homozygous<-data.frame(homozygous_1000G=as.vector(tables[[3]][2,2]),homozygous_ExAC=as.vector(tables[[3]][3,2]))    }
    }
    
  }else{
    ##########################################################
    #                          SNP                           #
    ##########################################################
    ##___________________________________
    ## EXONIC SNP
    ##___________________________________
    if(str_detect(variants["HGVS_c"],"[\\*+-]")==FALSE){
      doc.text<-read_html(paste("http://www.mutationtaster.org/cgi-bin/MutationTaster/MutationTaster69.cgi?position_be=",str_extract(variants["HGVS_c"],"\\d.[\\d]"),"&new_base=",str_extract(variants["HGVS_c"],"(?<=\\>).*"),"&transcript_stable_id_text=",variants["ENSTranscript"],"&sequence_type=CDS",sep=""))%>%
        html_nodes("h3") %>% 
        html_text()
      if(length(doc.text)==0){
        doc.text <- read_html(paste("http://www.mutationtaster.org/cgi-bin/MutationTaster/MutationTaster69.cgi?position_be=",str_extract(variants["HGVS_c"],"\\d.[\\d]"),"&new_base=",str_extract(variants["HGVS_c"],"(?<=\\>).*"),"&transcript_stable_id_text=",variants["ENSTranscript"],"&sequence_type=CDS",sep="")) %>% 
          html_nodes("h2") %>% 
          html_text()
      }
      #read HTML tables
      tables <- read_html(paste("http://www.mutationtaster.org/cgi-bin/MutationTaster/MutationTaster69.cgi?position_be=",str_extract(variants["HGVS_c"],"\\d.[\\d]"),"&new_base=",str_extract(variants["HGVS_c"],"(?<=\\>).*"),"&transcript_stable_id_text=",variants["ENSTranscript"],"&sequence_type=CDS",sep="")) %>%
        html_table(fill=T)
      #get precition result
      toMatch <- c("polymorphism","disease causing","wrong input format","data problem","annotation problem")
      MutationTaster<-unique(grep(paste(toMatch,collapse="|"), 
                                  doc.text, value=TRUE))
      #get homozygous counts in 1000G and ExAC
      homozygous<-data.frame(homozygous_1000G=NA,homozygous_ExAC=NA)
      if((!is.null(MutationTaster))&&(!is.na(MutationTaster))){
        if(MutationTaster!=("data problem")&&MutationTaster!=("wrong input format")&&MutationTaster!=("annotation problem")){
          homozygous<-data.frame(homozygous_1000G=as.vector(tables[[3]][2,2]),homozygous_ExAC=as.vector(tables[[3]][3,2]))    }
      }
    }else{
      ##___________________________________
      ## INTRONIC SNP
      ##___________________________________
      ## get DNA snippet
      snippet<-paste(getSeq(Hsapiens,variants["Chr"],
                            as.numeric(as.character(variants["Position"]))-20,as.numeric(as.character(variants["Position"]))-1),
                     paste("[",as.character(variants["Ref"]),"/",as.character(variants["Alt"]),"]",sep=""),
                     getSeq(Hsapiens,as.character(variants["Chr"]),as.numeric(as.character(variants["Position"]))+1,as.numeric(as.character(variants["Position"]))+20),
                     sep='')
      
      ## get reverse complement DNA snnipet if gene is located in negative strand
      if((variants["Strand"]=="-")&&(!is.na(variants["Strand"]))){
        snippet<-reverseDNA(snippet, complement = T)
      }
      
      # Read and parse HTML file
      session<-html_session("http://www.mutationtaster.org/")
      form<-set_values(html_form(session)[[1]],
                       transcript_stable_id_text=variants["ENSTranscript"],
                       sequence_type="gDNA",
                       sequence_snippet=snippet)
      doc.text <- submit_form(session, form) %>%
        read_html() %>% 
        html_nodes("h3") %>% 
        html_text()
      if(length(doc.text)==0){
        doc.text <- submit_form(session, form) %>%
          read_html() %>% 
          html_nodes("h2") %>% 
          html_text()
      }
      #read HTML tables
      tables <- submit_form(session, form) %>%
        read_html() %>% html_table(fill=T)
      #get precition result
      toMatch <- c("polymorphism","disease causing","wrong input format","data problem","annotation problem")
      MutationTaster<-unique(grep(paste(toMatch,collapse="|"), 
                                  doc.text, value=TRUE))
      #get homozygous counts in 1000G and ExAC
      homozygous<-data.frame(homozygous_1000G=NA,homozygous_ExAC=NA)
      if((!is.null(MutationTaster))&&(!is.na(MutationTaster))){
        if(MutationTaster!=("data problem")&&MutationTaster!=("wrong input format")&&MutationTaster!=("annotation problem")){
          homozygous<-data.frame(homozygous_1000G=as.vector(tables[[3]][2,2]),homozygous_ExAC=as.vector(tables[[3]][3,2]))    }
      }
    }
  }
  return(data.frame(MutationTaster,
                    homozygous_1000G=as.numeric(as.character(homozygous[,1])),
                    homozygous_ExAC=as.numeric(as.character(homozygous[,2]))))
}
##########################################################
#            trying to make a better function            #
##########################################################
#mutationTaster<-function(Chr,Position,Ref,Alt,Strand,ENSTranscript){
#  url<-try(read_html(paste("http://www.mutationtaster.org/cgi-bin/MutationTaster/MT_ChrPos.cgi?chromosome=",sub(".{3}","",Chr),"&position=",Position,"&ref=",Ref,"&alt=",Alt,sep="")))%>%
#    html_nodes("td a") %>%
#    html_attr("href") %>% 
#    grep(ENSTranscript,.,value=T)
  
#  MutationTaster<-read_html(url[1]) %>%
#    html_nodes("h3") %>%
#    html_text()
#  if(length(MutationTaster)==0){
#    MutationTaster<-read_html(url[1]) %>%
#      html_nodes("h2") %>%
#      html_text()
#  }
  
#  tables<-read_html(url[1]) %>%
#    html_table(fill=T)
  
#  homozygous<-data.frame(homozygous_1000G=NA,homozygous_ExAC=NA)
#  if((!is.null(MutationTaster))&&(!is.na(MutationTaster))){
#    if(MutationTaster!=("data problem")&&MutationTaster!=("wrong input format")&&MutationTaster!=("annotation problem")){
#      homozygous<-data.frame(homozygous_1000G=as.vector(tables[[3]][2,2]),homozygous_ExAC=as.vector(tables[[3]][3,2]))    }
#  }
#  return(data.frame(MutationTaster[2],homozygous))
#}