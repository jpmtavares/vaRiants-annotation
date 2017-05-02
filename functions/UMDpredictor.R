##########################################################
#                                                        #
#                      UMD-predictor                     #
#                                                        #
##########################################################
UMDpredictor<-function(genes){
  transcripts<-read.delim("./sources/GRCh37ENStranscripts.txt",
                          header = T)
  u<-semi_join(transcripts, as.data.frame(genes), by = c("HGNC_symbol" = "genes"))
  u$url<-paste("http://umd-predictor.eu/transcript_query2.php?name=",
               str_extract(u$ENSTranscript,"(?<=0{5}).*$"),sep="")
  
  umd<-lapply(u$url,function(x) {
    tryCatch({read.delim(x,header = F)},
             error=function(e){NULL})
  })
  u<-u[!sapply(umd, is.null),]
  umd<-umd[!sapply(umd, is.null)]
  
  umd<-mapply(cbind, "HGNC_symbol"=u$HGNC_symbol, "ENSTranscript"=u$ENSTranscript, umd,  SIMPLIFY=F)
  
  umdpredictor<-do.call("rbind", umd)
  names(umdpredictor)[c(3:8)]<-c("TranscriptPosition","Position","HGVS_c","HGVS_p",
                                 "score","UMD-predictor")
  return(umdpredictor)
}
