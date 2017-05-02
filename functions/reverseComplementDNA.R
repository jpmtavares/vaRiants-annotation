##########################################################
#                                                        #
#             reverse (and complement) DNA               #
#                                                        #
##########################################################
reverseDNA<-function(seq, complement=F){
  if(complement==T){
    split<-strsplit(as.character(chartr("TCGA[]", "AGCT][", seq)),"")[[1]]
    paste(rev(split),collapse = "")
  }else{
    chartr("TCGA[]", "AGCT][", seq)
  }
}