##########################################################
#                                                        #
#              processing vcf FORMAT field               #
#                                                        #
##########################################################
# The final vcf is an ensemble of multiple variant callers
# which results in different FORMAT fields

vcfFORMAT<-function(format, somatic = FALSE){
  #split FORMAT information using ":"; outputs table with 9 collumns
  FORMAT<-str_split_fixed(format[[sample]], ":", 9)
  
  FORMAT<-apply(FORMAT, 2, function(x) gsub("^$|^ $", NA, x)) %>% #replace all empty cells (or cells with white-space) with NA
    as.data.frame() %>% #change matrix to data.frame
    mutate(rowNumber = as.numeric(rownames(.))) #keep row numbers
  #count number of NA row.wise
  FORMAT$na_count<- apply(FORMAT, 1, function(x) sum(is.na(x))) #germline na_count: 0, 4, 6
                                                                #somatic na_count: 2, ?
  #___________________________________________
  # get coverage
  #___________________________________________
  if(somatic==FALSE){ #for germline samples
    cov<-data.frame(coverage=FORMAT[,3])
    cov$coverage_ref<-ifelse(FORMAT$na_count==0, as.numeric(str_extract(FORMAT[,4],"[^,]*")),
                             ifelse(FORMAT$na_count==4, as.numeric(str_extract(FORMAT[,2],"[^,]*")), NA))
    cov$coverage_alt<-ifelse(FORMAT$na_count==0, as.numeric(str_extract(FORMAT[,4],"(?<=,).*")),
                             ifelse(FORMAT$na_count==4, as.numeric(str_extract(FORMAT[,2],"(?<=,).*")), NA))
  }else{ #for somatic samples
    cov$coverage<-ifelse(FORMAT$na_count==2,as.numeric(as.character(FORMAT[,2])),
                         as.numeric(as.character(FORMAT[,4])))
    cov$coverage_ref<-ifelse(FORMAT$na_count==2, as.numeric(str_extract(FORMAT[,4],"[^,]*")),
                             as.numeric(str_extract(FORMAT[,2],"[^,]*")))
    cov$coverage_alt<-ifelse(FORMAT$na_count==2, as.numeric(str_extract(FORMAT[,4],"(?<=,).*")),
                             as.numeric(str_extract(FORMAT[,2],"(?<=,).*")))
  }
  
  #___________________________________________
  # genotype field - human-readable
  #___________________________________________
  original<-c("0/0","0/1","1/0","1/1","1/.","./1",
              "0|0","0|1","1|0","1|1","1|.",".|1")
  replacement<-c("homozygous reference","heterozygous","heterozygous","homozygous",
                 "heterozygous","heterozygous",
                 "homozygous reference","heterozygous","heterozygous","homozygous",
                 "heterozygous","heterozygous")
  #replace original annotation by homozygous/heterozygous 
  cov$genotype<-mgsub(original, replacement, FORMAT[,1])
  
  #___________________________________________
  # get variants and combine coverage and genotype information
  #___________________________________________
  if(any(names(format)=="ANN")==FALSE){
    format$ANN=NA
  }else{
    if(any(names(format)=="CSQ")==FALSE){
      format$CSQ=NA
    }
  }
  variants<-data.frame(format[,c("CHROM","POS","avsnp147","REF","ALT")],
                       cov, format[,c("Gene.refGene","ANN","CSQ")])
  
  return(variants)
}