##########################################################
#                                                        #
#             variants in Genes of Interest              #
#                                                        #
##########################################################
variants_genesofinterest<-function(variants,g,genes){
  
  variants_genes<-variants %>%
    separate_rows(ANN, sep=",") %>% # separate rows by ANN collumn separator ","
    filter(grepl(g, ANN) | grepl(g, CSQ)) %>% # find genes of interest in ANN and CSQ collumns
    mutate(ANNO=ifelse(!grepl(g, ANN), CSQ, ANN)) %>% # merge cells in collumns ANN and CSQ if they have genes of interest, in new collumn ANNO
    mutate(Gene.refGene=unlist(lapply(strsplit(ANNO,"|",fixed=T), "[[", 4))) %>% #change refGene names
    mutate(variant_type=unlist(lapply(strsplit(ANNO,"|",fixed=T), "[[", 2))) %>% #add variant_type
    select(-c(CSQ,ANN,ANNO)) %>%
    set_names(c("Chr","Position","rs_ID","Ref","Alt","coverage","coverage_ref",
                "coverage_alt","genotype","HGNC_symbol","variant_type")) %>% #set collumn names
    mutate(Chr=paste("chr",.$Chr,sep="")) %>% #correct chr names
    filter(HGNC_symbol%in%genes) %>% #select exact gene match
    #filter(!grepl("upstream|downstream", variant_type)) %>% # filter out intergenic variants
    unique %>%
    indelsCoordinates(.) # correct coordinates for indels
  
  return(variants_genes)
}