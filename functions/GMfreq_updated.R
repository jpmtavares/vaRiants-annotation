##########################################################
#                                                        #
#                    GenoMed Frequency                   #
#                                                        #
##########################################################
GMfrequency<-function(GM_freq,path){
  
  #get date of the last todas_ file
  updated_todas<-max(str_extract(list.files(paste(path,"/Archive/VCF/",sep=""), pattern = "todas_"),
                                 "[0-9]{8}"))
  #get date of the last GMfreq file
  updated_GMfreq<-max(str_extract(list.files("./sources/", pattern = "GM_freq"),
                                  "[0-9]{8}"))
  
  if(updated_todas>updated_GMfreq){
    #___________________________________________________________
    # import todas.vcf.gz
    #___________________________________________________________
    todas<-read.vcf(paste(path,"/Archive/VCF/todas_",updated_todas,".vcf.gz",sep=""))
    todas<-todas$vcf
    #___________________________________________________________
    # no. of samples and alleles
    #___________________________________________________________
    nsamples<-todas %>%
      ncol()-9
    nalleles<-nsamples*2
    #___________________________________________________________
    # get number of mutated alleles per sample
    #___________________________________________________________
    alleles<-todas %>%
      mutate_at(vars(-CHROM, -POS, -ID, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT),
                funs(replace(., grepl("0/", .), 1))) %>%
      mutate_at(vars(-CHROM, -POS, -ID, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT),
                funs(replace(., grepl("1/", .), 2))) %>%
      select(-CHROM, -POS, -ID, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT) %>%
      replace(., is.na(.), 0) %>%
      mutate_if(sapply(., is.character), as.numeric)
    #___________________________________________________________
    # get number of mutated samples, homozygous and MAF
    #___________________________________________________________
    inHouse<-alleles %>%
      transmute(inHouse_samples=apply(alleles, 1, function(x) length(which(x>0)))) %>%
      mutate(inHouse_homozygous=apply(alleles, 1, function(x) length(which(x==2)))) %>%
      mutate(inHouse_maf=rowSums(alleles)/nalleles)
    
    #___________________________________________________________
    # get output table
    #___________________________________________________________
    freq<-data.frame(todas[,c("CHROM", "POS", "REF", "ALT")],
                     inHouse) %>%
      mutate(CHROM=paste("chr",CHROM, sep="")) %>%
      set_names(c("Chr","Position","Ref","Alt","inHouse_samples", "inHouse_homozygous", "inHouse_MAF"))
    #___________________________________________________________
    # write.table
    #___________________________________________________________
    output<-paste("GM_freq",format(Sys.time(), "%Y%m%d"),".txt",sep="")
    write.table(freq, paste("./sources/",output,sep=""), row.names=F, col.names=T,
                quote=F,sep="\t")
    #___________________________________________________________
    # save to bcbio_pipeline.Rdata
    #___________________________________________________________
    save(refSeqGenes,freq,file="./sources/bcbio_pipeline.Rdata")
  }else{
    freq<-GM_freq
  }
  return(freq)
}