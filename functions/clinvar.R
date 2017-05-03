###################################################################
#                                                                 #
#              ClinVar: Download and processing file              #
#                                                                 #
###################################################################
##!!!!WARNING!!!!##
# it's imperative to replace all # \x2c and ' characters by _ (for example) except header
##!!!!WARNING!!!!##

#______________________________________________
# list files in NCBI Clinvar, and check dates
#______________________________________________
clinvarTab<-function(x){
  url<-"ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_[0-9]*.vcf.gz"
  filename<-grep("clinvar_[0-9]*.vcf.gz$",unlist(strsplit(getURL(url, dirlistonly = TRUE), split = "\n")),
                 value=T)
  
  date_updated<-ymd(str_extract(filename,"[0-9]{8}")) #get date of updated file in NCBI ClinVar
  date_last<-ymd(max(str_extract(list.files("./sources/", pattern = ".vcf.gz"),
                                 "[0-9]{8}"))) #get date of last downloaded file
  
  if(date_updated > date_last){ #if there is a newer file in NCBI Clinvar
    #then download it and process it
    download.file(url, paste("./sources/",filename,sep=""), method = "wget")
    
    #_________________________________________
    # preprocessing
    #_________________________________________
    system(paste("./utils/clinvar_preprocessing.sh '",str_extract(filename,".*[\\d]"),"'",sep=""))
    
    #_________________________________________
    # read vcf.gz files
    #_________________________________________
    clinvarVcf<-read.vcf(paste("./sources/",filename,sep=""), split.info = TRUE)
    
    #_________________________________________
    # processing vcf.gz files
    #_________________________________________
    
    ###list2df
    clinvarTab<-clinvarVcf$vcf[c("CHROM","POS","ID","REF","ALT","CLNSIG", "CLNDBN","CLNREVSTAT")]
    
    ###human-readable
    library(qdap)
    #CLNSIG
    clnsig_original<-c(0,1,2,3,4,5,6,7,255)
    clnsig_replacement<-c("Uncertain significance","not provided","Benign",
                          "Likely-benign","Likely-pathogenic","Pathogenic","drug-response",
                          "histocompatibility","other")
    clinvarTab$CLNSIG<-mgsub(clnsig_original, clnsig_replacement, clinvarTab$CLNSIG)
    
    #CLNREVSTAT
    clnrevstat_original<-c("no_assertion","no_criteria","single","mult","conf","exp","guideline")
    clnrevstat_replacement<-c("no stars","no stars",
                              "*","**",
                              "*","***",
                              "****")
    #clnrevstat_replacement<-c("No assertion provided","No assertion criteria provided",
    #                          "Criteria provided single submitter","Criteria provided multiple submitters no conflicts",
    #                          "Criteria provided conflicting interpretations","Reviewed by expert panel",
    #                          "Practice guideline")
    clinvarTab$CLNREVSTAT<-mgsub(clnrevstat_original, clnrevstat_replacement, clinvarTab$CLNREVSTAT)
    
    ###unifying file format 
    clinvar<-clinvarTab %>%
      mutate(CHROM=paste("chr",.$CHROM,sep=""),
             CLNDBN=gsub("#","_",as.character(.$CLNDBN)),
             CLNDBN=gsub("'","_",as.character(.$CLNDBN)),
             CLNDBN=gsub("\\\\x2c",",",as.character(.$CLNDBN))) %>%
      set_names(c("Chr","Start","rs_ID","Ref","Alt","clinvar_sig","clinvar_disease_name","clinvar_stars")) %>%
      mutate(Alt = strsplit(as.character(Alt), ",")) %>%
      unnest(Alt) %>%
      select(Chr, Start, rs_ID, Ref, Alt, clinvar_sig, clinvar_disease_name, clinvar_stars)
    
    #_________________________________________
    # write.file
    #_________________________________________
    write.table(clinvar,paste("./sources/",str_extract(filename,".*[\\d]"),".txt",sep=""),
                sep="\t",col.names=T,row.names=F,quote=F)
  }else{
    clinvar<-read.delim(paste("./sources/",grep(max(str_extract(list.files("./sources/", pattern = ".vcf.gz"), "[0-9]{8}")),
                                                list.files("./sources/", pattern = ".txt"),value = T),sep=""),
                        header = T)
  }
  return(clinvar)
}