###################################################################
#                                                                 #
#              ClinVar: Download and processing file              #
#                                                                 #
###################################################################
##!!!!WARNING!!!!##
# it's imperative to replace all # \x2c and ' characters by _ (for example) except header
##!!!!WARNING!!!!##

##########################################################
#                list and download file                  #
##########################################################
url<-"ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_[0-9]*.vcf.gz"
filename<-grep("clinvar_[0-9]*.vcf.gz$",unlist(strsplit(getURL(url, dirlistonly = TRUE), split = "\n")),
               value=T)

date<-ymd(str_extract(filename,"[0-9]{8}")) #get date of updated file in NCBI ClinVar
if(between(as.numeric(today()-date),0,7)){ #if there is a 7-days time lapse, download file
  temp<-tempfile()
  temp<-download.file(url, paste(temp,".vcf.gz",sep=""), method = "wget")
  
  #_________________________________________
  # read vcf.gz files
  #_________________________________________
  clinvarVcf<-read.vcf(paste(temp,".vcf.gz",sep=""), split.info = TRUE)
} 