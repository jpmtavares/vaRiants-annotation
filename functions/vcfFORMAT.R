##########################################################
#                                                        #
#              processing vcf FORMAT field               #
#                                                        #
##########################################################
# The final vcf is an ensemble of multiple variant callers
# which results in different FORMAT fields

vcfFORMAT<-function(format, somatic = FALSE){
  cov<-str_split_fixed(format, ":", 9)
  cov<-apply(cov, 2, function(x) gsub("^$|^ $", NA, x))
  cov<-as.data.frame(cov)
  cov$na_count<- apply(cov, 1, function(x) sum(is.na(x)))
  
  #genotype field - human-readable
  cov_original<-c("0/0","0/1","1/0","1/1","1/.","./1",
                  "0|0","0|1","1|0","1|1","1|.",".|1")
  cov_replacement<-c("homozygous reference","heterozygous","heterozygous","homozygous",
                     "heterozygous","heterozygous",
                     "homozygous reference","heterozygous","heterozygous","homozygous",
                     "heterozygous","heterozygous")
  
  cov$genotype<-mgsub(cov_original, cov_replacement, cov[,1])
  
  if(somatic==FALSE){
    cov$coverage<-cov[,3]
    cov$coverage_ref<-ifelse(cov$na_count==0, as.numeric(str_extract(cov[,4],"[^,]*")),
                             ifelse(cov$na_count==4, as.numeric(str_extract(cov[,2],"[^,]*")), NA))
    cov$coverage_alt<-ifelse(cov$na_count==0, as.numeric(str_extract(cov[,4],"(?<=,).*")),
                             ifelse(cov$na_count==4, as.numeric(str_extract(cov[,2],"(?<=,).*")), NA))
  }else{
    cov$coverage<-ifelse(cov$na_count==2,as.numeric(as.character(cov[,2])),
                         as.numeric(as.character(cov[,4])))
    cov$coverage_ref<-ifelse(cov$na_count==2, as.numeric(str_extract(cov[,4],"[^,]*")),
                             as.numeric(str_extract(cov[,2],"[^,]*")))
    cov$coverage_alt<-ifelse(cov$na_count==2, as.numeric(str_extract(cov[,4],"(?<=,).*")),
                             as.numeric(str_extract(cov[,2],"(?<=,).*")))
  }
  
  return(data.frame(cov[,c("coverage","coverage_ref","coverage_alt","genotype")]))
}