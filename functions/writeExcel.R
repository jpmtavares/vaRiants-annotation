##########################################################
#                                                        #
#                  Output Excel file                     #
#                                                        #
##########################################################

writeExcel<-function(variants, notFound){
  
  novenove<-variants %>%
    filter(str_detect(HGVS_c,"[\\*\\+\\-][0-9]{1,2}[diATGC]")==TRUE | #all variants up to 99bp from SS
             str_detect(HGVS_c,"[\\*\\+\\-]")==FALSE) #add exonic
  
  non_novenove<-variants %>%
    filter(str_detect(HGVS_c,"[\\*\\+\\-][0-9]{1,2}[diATGC]")==FALSE) %>% #all variants beyond 99bp from SS
    filter(str_detect(HGVS_c,"[\\*\\+\\-]")==TRUE) #remove exonic
  
  cinquenta<-novenove %>%
    filter(str_detect(HGVS_c,"[\\*\\+\\-][5-9][0-9][diATGC]")==FALSE)
  
  off<-novenove %>%
    filter(str_detect(HGVS_c,"[\\*\\+\\-][5-9][0-9][diATGC]")==TRUE) %>%
    rbind(non_novenove)
  
  cov15<-cinquenta %>%
    filter(as.numeric(as.character(coverage)) > 14)
  
  covless15<-cinquenta %>%
    filter(as.numeric(as.character(coverage))<=14)
  
  #______________________________________________
  # WARNING: genes not found in vcf file in GRCh38 version
  #______________________________________________
  
  assign("last.warning", NULL, envir = baseenv())
  if(length(notFound)!=0){
    Warning<-data.frame(Warnings=paste("No variant was found in gene(s)",
                  paste(notFound, collapse = ", "), sep=" "))
  }else{
    Warning<-data.frame(Warnings="Finished with no warnings.")
  }
  
  return(list(cov15, covless15, off, Warning))
}