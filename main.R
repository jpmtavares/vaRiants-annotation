###########################################################################
#                                                                         #
#                    Clinical Exome Analysis (exo)                        #
#                                                                         #
###########################################################################

#!/bin/Rscript
#______________________________________________
# set work directory, sample and genes
#______________________________________________
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("./")

#______________________________________________
# libraries
#______________________________________________
library(bedr)
library(stringr)
library(plyr)
library(dplyr)
library(tidyr)
library(qdap)
library(magrittr)
library(rvest)
library(readr)
library(XML)
library(RCurl)
library(lubridate)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
#______________________________________________
# R source files
#______________________________________________
source("./functions/HGVS_mutalyzer.R")
source("./functions/reverseComplementDNA.R")
source("./functions/UMDpredictor.R")
source("./functions/vcfFORMAT.R")
source("./functions/mutationTaster.R")
source("./functions/clinvar.R")
source("./functions/GMfreq.R")
source("./functions/genomiser.R")
source("./functions/hgmdLOVD.R")
source("./functions/gwas.R")
source("./functions/genesofinterest.R")
source("./functions/indelsCoordinates.R")
source("./functions/writeExcel.R")
source("./functions/variants_genesofinterest.R")
source("./functions/variants_transcripts.R")

#Samplename
sample<-str_extract(list.files("../",pattern="*multianno.vcf"),"[^\\.]*")
#interesting genes
gene.list<-"SCN5A"
#g<-paste(as.character(commandArgs(TRUE)[1]),sep="")
g<-unlist(genesofinterest(gene.list)[1])
#genes<-unlist(strsplit(g,"\\|"))
genes<-unlist(genesofinterest(gene.list)[2])

#______________________________________________
# preparing input files
#______________________________________________
load("./sources/bcbio_pipeline.Rdata") ## load previously saved .Rdata with reference transcripts, clinvar, and inHouse variant frequency
## save(refSeqGenes,freq,file="./sources/bcbio_pipeline.Rdata")

#UMD-predictor
umd<-UMDpredictor(genes)

#GM_freq
GM_freq<-GMfrequency(freq,str_extract(getwd(),"^.*(?=(/github))"))
#______________________________________________
# import tables
#______________________________________________
#variant calls (by at least two callers -> ensemble) & ANNOVAR annotation
anno<-read.vcf(paste("../",list.files("../",pattern="*multianno.vcf"),sep=""),
               split.info = T)
names(anno$vcf)[ncol(anno$vcf)]<-sample ##change last column name to SAMPLE

#______________________________________________
# processing basic transcript annotation
#______________________________________________
##get variants and combine coverage and genotype information
variants<-vcfFORMAT(anno$vcf)

##filter by genes of interest
variants_genes<-variants_genesofinterest(variants,g,genes)

## get transcript annotation from refSeqGenes file
trans_anno<-variants_transcripts(variants_genes,refSeqGenes)

#______________________________________________
# check if all genes of interest are present in refSeqGenes
#______________________________________________
if(any(is.na(trans_anno$Strand))){
  stop(paste("Add gene(s)",   # if not, stop analysis right here
             paste(as.character(unique(trans_anno$HGNC_symbol[is.na(trans_anno$Strand)])), collapse=", "),
             "to refSeqGenes.",sep=" "))
}else{ #______________________________________________ else, proceed
  
  #replace old gene nomenclature by GRCh38
  grch38<-as.character(unlist(genesofinterest(gene.list)[4]))
  grch37<-as.character(unlist(genesofinterest(gene.list)[3]))
  trans_anno<-data.frame(lapply(trans_anno, function(x) {
    mgsub(grch37, grch38, x) })) %>%
    unique()
  
  #______________________________________________
  # WARNING: genes not found in vcf file in GRCh38 version
  #______________________________________________
  notFound<-genes[!(genes%in%trans_anno$HGNC_symbol)]%>%
    .[!(.%in%grch37)]
  
  #______________________________________________
  # HGVS_mutalyzer
  #______________________________________________
  hgvs_anno<-hgvs(trans_anno) %>%
    cbind(.,type=mapply(variant_type, .[,"HGVS_c"], .[,"HGVS_p"]))
  
  #_________________________________________
  # NCBI ClinVar
  #_________________________________________
  clinvar<-clinvarTab() %>%
    join(hgvs_anno,.) %>%
    select(-Start)
  
  #_________________________________________
  # HGMD & LOVD (presence/absence)
  #_________________________________________
  #hgmd_lovd<-hgmdLOVD(clinvar)
  
  #_________________________________________
  # GWAS catalog
  #_________________________________________
  gwas<-GWAS(clinvar)
  
  #_________________________________________
  # HSF
  #_________________________________________
  HSF<-data.frame(gwas, HSF="error")
  
  #_________________________________________
  # MutationTaster
  #_________________________________________
  MT_anno<-HSF %>%
    apply(1,MT) %>% do.call(rbind.data.frame,.) %>%
    cbind(HSF, .)
  
  #correct HGVS
  MT_anno$HGVS_p<-mapply(synonymous,MT_anno$HGVS_p)
  
  #correct homozygous counting
  MT_anno$homozygous_ExAC[grep("\\.",MT_anno$homozygous_ExAC)]<-"-"
  MT_anno$homozygous_1000G[grep("\\.",MT_anno$homozygous_1000G)]<-"-"
  
  #_________________________________________
  # UMD-predictor
  #_________________________________________
  UMD<-MT_anno %>%
    join(umd) %>%
    select(-c(TranscriptPosition,score))
  
  UMD$HGVS_p<-ifelse(!is.na(UMD$HGVSp) & as.character(UMD$HGVSp) != as.character(UMD$HGVS_p),
                     as.character(UMD$HGVSp), as.character(UMD$HGVS_p))
  
  UMD<-select(UMD,-HGVSp) %>%
    mutate(type=variant_type(UMD$HGVS_c, UMD$HGVS_p)) %>%
    unique()
  
  #_________________________________________
  # Genomiser
  #_________________________________________
  if(length(list.files("../","*variants.tsv")) >0){
    Genomiser<-genomiser(UMD)
    UMD<-Genomiser
  }
  
  #_________________________________________
  # ANNOVAR predictions
  #_________________________________________
  predictions<-anno$vcf %>%
    select(CHROM,POS,avsnp147,REF,ALT,SIFT_score,SIFT_pred,PROVEAN_score,
           PROVEAN_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,
           MutationAssessor_score,MutationAssessor_pred,CADD13_PHRED,DANN_score,
           FATHMM_coding,FATHMM_noncoding,GWAVA_region_score,`GERP.._RS`,
           phyloP20way_mammalian,phastCons20way_mammalian,SiPhy_29way_logOdds,
           dpsi_max_tissue,PopFreqMax,X1000G_EUR,X1000G_ALL,gnomAD_genome_NFE,gnomAD_genome_ALL,
           ESP6500siv2_EA,ESP6500siv2_ALL) %>%
    set_names(c("Chr","Position","rs_ID","Ref","Alt","SIFT_score","SIFT_pred",
                "PROVEAN_score","PROVEAN_pred","Polyphen2_HVAR_score",
                "Polyphen2_HVAR_pred","MutationAssessor_score","MutationAssessor_pred",
                "CADD13_PHRED","DANN_score","FATHMM_coding","FATHMM_noncoding",
                "GWAVA_region_score","GERP++RS","phyloP20_mammalian",
                "phastCons20_mammalian","SiPhy","SPIDEX","PopFreqMax","1000G_EUR",
                "1000G_ALL","gnomAD_EUR","gnomAD_ALL","ESP_EUR","ESP_ALL")) %>%
    mutate(Chr=paste("chr",.$Chr,sep="")) %>%
    indelsCoordinates(.) %>% # correct indels coordinates
    join(UMD,.)
  
  #______________________________________________
  # human-readable predictions
  #______________________________________________
  #SIFT_pred
  predictions$SIFT_pred<-ifelse(as.character(predictions$SIFT_pred)=="T","tolerated",
                                ifelse(as.character(predictions$SIFT_pred)=="D","deleterious",as.character(predictions$SIFT_pred)))
  #Polyphen2_HVAR_pred
  predictions$Polyphen2_HVAR_pred<-ifelse(as.character(predictions$Polyphen2_HVAR_pred)=="D","probably damaging",
                                          ifelse(as.character(predictions$Polyphen2_HVAR_pred)=="P","possibly damaging",
                                                 ifelse(as.character(predictions$Polyphen2_HVAR_pred)=="B","benign",as.character(predictions$Polyphen2_HVAR_pred))))
  #MutationAssessor
  predictions$MutationAssessor_pred<-ifelse(as.character(predictions$MutationAssessor_pred)=="L","low",
                                            ifelse(as.character(predictions$MutationAssessor_pred)=="N","neutral",
                                                   ifelse(as.character(predictions$MutationAssessor_pred)=="M","medium",
                                                          ifelse(as.character(predictions$MutationAssessor_pred)=="H","high",as.character(predictions$MutationAssessor_pred)))))
  #PROVEAN
  predictions$PROVEAN_pred<-ifelse(as.character(predictions$PROVEAN_pred)=="N","neutral",
                                   ifelse(as.character(predictions$PROVEAN_pred)=="D","damaging",as.character(predictions$PROVEAN_pred)))
  #SPIDEX
  predictions$SPIDEX<-abs(as.numeric(predictions$SPIDEX))
  
  #______________________________________________
  # GM_freq
  #______________________________________________
  GMfreq<-GM_freq %>%
    indelsCoordinates(.) %>% # correct indels coordinates
    join(predictions,.) %>%
    select(-N_samples) %>%
    unique()
  
  GMfreq[is.na(GMfreq$HGVS_c),"HGVS_c"]<-"."
  
  #______________________________________________
  # Preparing output in excel
  #______________________________________________
  cov15<-as.data.frame(writeExcel(GMfreq,notFound)[1])
  covless15<-as.data.frame(writeExcel(GMfreq,notFound)[2])
  off<-as.data.frame(writeExcel(GMfreq,notFound)[3])
  Warning<-as.data.frame(writeExcel(GMfreq,notFound)[4])
  
  #______________________________________________
  # Writing table
  #______________________________________________
  # directory creation
  dir.create("./my_analysis", recursive=F, mode = "0777")
  library("WriteXLS")
  #write.table(GMfreq,paste("./my_analysis/",sample,"_all_variants.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
  
  WriteXLS(c("cov15","covless15","off","Warning"), ExcelFileName = paste("./my_analysis/",sample,"_genes_variants.xlsx",sep=""), SheetNames = c("High confidence","Low confidence","Off target","Warnings"), perl = "perl",
           row.names = FALSE, col.names = TRUE,
           envir = parent.frame())
}