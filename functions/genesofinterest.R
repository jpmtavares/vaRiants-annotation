###################################################################
#                                                                 #
#                     Genes of Interest                           #
#                                                                 #
###################################################################

genesofinterest<-function(genesList){
  
  #read genes that changed names between GRCh37 and GRCh38
  grch37vs38<-read.delim("./sources/grch37vsgrch38.txt")
  
  #processing genes list
  genes<-unlist(strsplit(genesList,split="[, ]+|\t")) %>%
    gsub("[, ()=]+|\n","",.)
  
  #get gene names that changed between genome versions
  twonames<-join(data.frame(GRCh38=genes),grch37vs38)
  
  #add old version of gene names to the list of genes of interest
  GenesOfInterest<-unique(c(genes,as.character(twonames[!is.na(twonames[,2]),3])))
  
  grch37<-unlist(twonames[!is.na(twonames[,2]),3])
  grch38<-unlist(twonames[!is.na(twonames[,2]),1])
  
  #get gene list separated by "|"
  g<-paste(GenesOfInterest,collapse = "|")
  
  return(list(g,GenesOfInterest,grch37,grch38))
  
}