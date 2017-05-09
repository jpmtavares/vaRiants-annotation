##########################################################
#                                                        #
#               Genomiser (Exomiser tool)                #
#                                                        #
##########################################################
genomiser<-function(variants){
  file<-read.delim(paste("../",list.files("../","*variants.tsv"),sep=""), header =T)
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
    setnames("EXOMISER_GENE_COMBINED_SCORE","Genomiser")%>%
    select(Genomiser)
  
  return(geno)
}