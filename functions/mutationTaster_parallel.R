##########################################################
#                                                        #
#                     MutationTaster                     #
#                                                        #
##########################################################
################################################################
#              get mutationTaster - main function              # 
################################################################
MT<-function(variants){
  ##__________________________________
  ## parallelize
  ##__________________________________
  registerDoParallel(cores=6)
  parallel<-foreach(v=variants$Chr, w=variants$Position, x=variants$Ref,
                    y=variants$Alt, z=variants$ENSTranscript) %dopar% mutationTaster(v,w,x,y,z)
  
  mt<-data.frame(variants,
                 MutationTaster=do.call(rbind.data.frame,parallel))
  return(mt)
}

################################################################
#           Make a GET request to get mutationTaster           #
################################################################
mutationTaster<-function(Chr,Position,Ref,Alt,ENSTranscript){
  requests<-import("requests")
  parameters<-dict(chromosome=gsub("chr","",Chr),
                   position=as.character(Position),
                   ref=as.character(Ref),
                   alt=as.character(Alt),
                   transcript_stable_id_text=as.character(ENSTranscript)) #is ensembl transcript working?
  response<-requests$get("http://www.mutationtaster.org/cgi-bin/MutationTaster/MT_ChrPos.cgi?", params = parameters)
  results<-data.frame(MutationTaster=str_extract(response$content,"polymorphism|disease causing|wrong input format|data problem|no suitable transcript"))
  
  return(results)
}