##########################################################
#                                                        #
#               correct indels coordinates               #
#                                                        #
##########################################################

indelsCoordinates<-function(variants){
  
  # if deletion or insertion, add 1 to Position
  # else, don't change anything
  variants$Position<-ifelse(nchar(as.character(variants$Ref))>1, # deletion
                            as.numeric(as.character(variants$Position))+1,
                            ifelse(nchar(as.character(variants$Alt))>1, # insertion
                                   as.numeric(as.character(variants$Position))+1,
                                   as.numeric(as.character(variants$Position))))
  
  return(variants)
  
}