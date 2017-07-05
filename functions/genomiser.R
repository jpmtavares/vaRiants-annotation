##########################################################
#                                                        #
#               Genomiser (Exomiser tool)                #
#                                                        #
##########################################################
genomiser<-function(variants){
  #read genomiser files
  myfiles<-lapply(paste("../",list.files("../","*variants.tsv",recursive=TRUE),sep=""), function(x){
    file<-read.delim(x)
    return(file) #remove header from vcf file
  })
  
  #list2data.frame
  file<-do.call(rbind.data.frame, myfiles)
  
  #file<-read.delim(paste("../",list.files("../","*variants.tsv"),sep=""), header =T)
  ##correct coordinates for indels
  file$POS<-ifelse(nchar(as.character(file$REF))>1,
                   as.numeric(as.character(file$POS))+1,
                   ifelse(nchar(as.character(file$ALT))>1,
                          as.numeric(as.character(file$POS))+1,
                          as.numeric(as.character(file$POS))))
  
  geno<-file %>%
    setnames("X.CHROM","Chr") %>% #set column names
    setnames("POS","Position") %>%
    setnames("REF","Ref") %>%
    setnames("ALT","Alt")%>%
    mutate(Chr=paste("chr",.$Chr,sep="")) %>% #correct chr names
    left_join(variants,.) %>%
    setnames("EXOMISER_VARIANT_SCORE","Genomiser")%>%
    .[,c(names(variants),"Genomiser")] %>%
    unique()
  
  return(geno)
}