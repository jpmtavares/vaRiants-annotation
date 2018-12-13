##########################################################
#                                                        #
#                         Utils                          #
#                                                        #
##########################################################

################################################################
#                      variant type                            #
################################################################
variant_type<-function(HGVS_c, HGVS_p){
  ##__________________________________
  ## if variant is exonic
  ##__________________________________
  ifelse(str_detect(HGVS_c,"[\\*\\+\\-]")==FALSE, # if variant is exonic
         ifelse(grepl("=",HGVS_p)==TRUE | str_extract(HGVS_p,"[^p.](?:(?!\\d).)*") == str_extract(HGVS_p,".{3}$"), "silent", 
                ifelse(grepl("Ter",HGVS_p)==TRUE, "nonsense",
                       ifelse(grepl("fs",HGVS_p)==TRUE, "frame-shift", "missense"))),
         ##__________________________________
         ## if variant is intronic
         ##__________________________________
         ifelse(grepl("\\*",HGVS_c)==TRUE, "3-UTR",
                ifelse(grepl("\\.\\-",HGVS_c)==TRUE, "5-UTR",
                       ifelse(grepl("[0-9][-+][0-9][diATGC]|[0-9][-+][1][0-9][diATGC]",HGVS_c)==TRUE, "splicing", "intronic"))))
}

################################################################
#    correct HGVS protein annotation for synonymous variants   #
################################################################
synonymous<-function(HGVS_p){
  ifelse(grepl("=",HGVS_p)==TRUE,
         gsub("=",str_extract(HGVS_p,"[^p.](?:(?!\\d).)*"),HGVS_p),
         as.character(HGVS_p))
}