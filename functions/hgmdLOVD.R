##########################################################
#                                                        #
#                      HGMD e LOVD                       #
#                                                        #
##########################################################

hgmdLOVD<-function(variants){
  hgmd<-getURL(paste("http://genome.ucsc.edu/cgi-bin/hgBeacon/query?dataset=hgmd&chromosome=",
                     str_extract(as.character(variants[,"Chr"]),"[X,Y,\\d].*"),
                     "&position=", as.character(variants[,"Position"]),
                     "&alternateBases=", as.character(variants[,"Alt"]), 
                     "&format=text",sep=""))
  hgmd<-str_extract(hgmd,"[^\\n].*")
  
  lovd<-getURL(paste("http://genome.ucsc.edu/cgi-bin/hgBeacon/query?dataset=lovd&chromosome=",
                     str_extract(as.character(variants[,"Chr"]),"[X,Y,\\d].*"),
                     "&position=", as.character(variants[,"Position"]),
                     "&alternateBases=", as.character(variants[,"Alt"]), 
                     "&format=text",sep=""))
  lovd<-str_extract(lovd,"[^\\n].*")
  
  return(data.frame(variants,
                    HGMD=hgmd,
                    LOVD=lovd))
}