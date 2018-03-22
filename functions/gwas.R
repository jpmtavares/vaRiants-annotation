##########################################################
#                                                        #
#       The NHGRI-EBI Catalog of published GWAS          #
#                                                        #
##########################################################
GWAS<-function(variants){
  Gstudies<-read.delim("./sources/gwas_catalog_v1.0.1-NHGRI-EBI.tsv",
                       header = T)
  gwas<-Gstudies %>%
    setnames("CHR_ID","Chr") %>% #set column names
    #setnames("CHR_POS","Position") %>%
    setnames("SNPS","rs_ID") %>%
    mutate(Chr=paste("chr",.$Chr,sep="")) %>% #correct chr names
    #mutate(CHR_POS=as.numeric(as.character(CHR_POS))) %>% #correct Position class
    mutate(rs_ID=as.character(rs_ID)) %>% #correct Position class
    left_join(variants,.) %>%
    setnames("LINK","GWAS")%>%
    .[,c(names(variants),"GWAS")]
  
  return(gwas)
}