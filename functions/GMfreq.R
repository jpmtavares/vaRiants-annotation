##########################################################
#                                                        #
#                    GenoMed Frequency                   #
#                                                        #
##########################################################
GMfrequency<-function(GM_freq){
  #___________________________________________________________
  # import ensemble.vcf.gz files
  #___________________________________________________________
  #get paths of ensemble.vcf.gz files
  temp<-system("ls -dR -1 /genomedarchive/Archive/Analysis/*/*/*/final/*/* | grep 'ensemble.vcf.gz$'",
               intern = T)
  temp<-temp[-grep("HCMcardio",temp)]
  
  if(length(temp) > GM_freq$N_samples[1]+8){ #if the number of samples in "Analysis" 
    #is 8x higher than the number of samples used in the last GMfreq file
    
    #read ensemble.vcf.gz files
    myfiles<-lapply(temp, function(x){
      file<-read.vcf(x)
      return(file$vcf) #remove header from vcf file
    }) 
    
    #___________________________________________________________
    # processing
    #___________________________________________________________
    #remove last column from vcf files
    myfiles<-lapply(myfiles, function(x) { x[,10] <- NULL; x })
    
    #get patients ID's
    samples<-temp %>% dirname() %>% dirname() %>% dirname() %>%
      basename()
    platform<-temp %>% dirname() %>% dirname() %>% dirname() %>% dirname() %>%
      basename()
    year<-temp %>% dirname() %>% dirname() %>% dirname() %>% dirname() %>% dirname() %>%
      basename()
    
    #add columns with patient ID, platform and year
    mysamples<-mapply(cbind, myfiles, "PatientID"=samples,
                      "Platform"=platform, "Year"=year, SIMPLIFY=F)
    
    #name each data.frame in mysamples list
    names(mysamples)<-samples
    #combine them in a single data.frame
    combine<-do.call("rbind", mysamples)
    
    #___________________________________________________________
    # get frequencies
    #___________________________________________________________
    count<-combine %>%
      select(CHROM, POS, REF, ALT) %>%
      group_by(CHROM, POS, REF, ALT) %>%
      summarise(inHouse_count=n()) %>%
      set_names(c("Chr","Position","Ref","Alt","inHouse_count"))
    
    freq<-data.frame(Chr=paste("chr",count$Chr, sep=""), count[,c(2:5)],
                     N_samples=length(samples),
                     inHouse_freq=count$inHouse_count/length(samples))
    
    #___________________________________________________________
    # write.table
    #___________________________________________________________
    output<-paste("GM_freq",format(Sys.time(), "%Y%m%d"),".txt",sep="")
    #write.table(freq, paste("./sources/",output,sep=""), row.names=F, col.names=T,
    #            quote=F,sep="\t")
    save(refSeqGenes,freq,file="./sources/bcbio_pipeline.Rdata")
  }else{
    freq<-GM_freq
  }
  return(freq)
}