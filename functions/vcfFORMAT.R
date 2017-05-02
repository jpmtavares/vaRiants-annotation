##########################################################
#                                                        #
#              processing vcf FORMAT field               #
#                                                        #
##########################################################
# The final vcf is an ensemble of multiple variant callers
# which results in different FORMAT fields

vcfFORMAT<-function(format){
  cov<-str_split_fixed(format, ":", 9)
  cov<-apply(cov, 2, function(x) gsub("^$|^ $", NA, x))
  cov<-as.data.frame(cov)
  cov$na_count<- apply(cov, 1, function(x) sum(is.na(x)))
  
  #genotype field - human-readable
  cov_original<-c("0/0","0/1","1/0","1/1","1/.","./1",
                  "0|0","0|1","1|0","1|1")
  cov_replacement<-c("homozygous reference","heterozygous","heterozygous","homozygous",
                     "heterozygous","heterozygous",
                     "homozygous reference","heterozygous","heterozygous","homozygous")
  
  cov$genotype<-mgsub(cov_original, cov_replacement, cov[,1])
  cov$coverage<-cov[,3]
  cov$coverage_ref<-ifelse(cov$na_count==0, as.numeric(str_extract(cov[,4],"[^,]*")),
                           ifelse(cov$na_count==4, as.numeric(str_extract(cov[,2],"[^,]*")), NA))
  cov$coverage_alt<-ifelse(cov$na_count==0, as.numeric(str_extract(cov[,4],"(?<=,).*")),
                           ifelse(cov$na_count==4, as.numeric(str_extract(cov[,2],"(?<=,).*")), NA))
  
  return(data.frame(cov[,c("coverage","coverage_ref","coverage_alt","genotype")]))
}