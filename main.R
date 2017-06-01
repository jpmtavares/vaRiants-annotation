###########################################################################
#                                                                         #
#                    Clinical Exome Analysis (exo)                        #
#                                                                         #
###########################################################################

#!/bin/Rscript

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
#______________________________________________
# set work directory, sample and genes
#______________________________________________
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("./")
sample<-basename(getwd())
#interesting genes
g<-paste(as.character(commandArgs(TRUE)[1]),sep="")
genes<-unlist(strsplit(g,"\\|"))
#______________________________________________
# preparing input files
#______________________________________________
load("./sources/bcbio_pipeline.Rdata") ## load previously saved .Rdata with reference transcripts, clinvar, and inHouse variant frequency
## save(refSeqGenes,GM_freq,file="./sources/bcbio_pipeline.Rdata")

#UMD-predictor
umd<-UMDpredictor(genes)

#GM_freq
GM_freq<-GMfrequency(freq)
#______________________________________________
# import tables
#______________________________________________
#variant calls (by at least two callers -> ensemble) & ANNOVAR annotation
anno<-read.vcf(paste("../",list.files("../",pattern="*multianno.vcf",recursive=TRUE),sep=""),
               split.info = T)
names(anno$vcf)[ncol(anno$vcf)]<-sample ##change last column name to SAMPLE

#______________________________________________
# processing basic transcript annotation
#______________________________________________
##coverage and genotype
cov<-vcfFORMAT(anno$vcf[[sample]])

##get variants and combine coverage and genotype information
variants_tmp<-data.frame(anno$vcf[,c("CHROM","POS","avsnp147","REF","ALT")], 
                         cov, anno$vcf[,c("Gene.refGene","CSQ","ANN")]) %>%
  filter(grepl(g,.[,11]) | grepl(g,.[,12])) %>% #variants in genes of interest
  mutate(Gene.refGene=unlist(lapply(strsplit(.[,11],"|",fixed=T), "[[", 4))) %>%
  filter(.[,10] %in% genes) %>%
  select(-c(CSQ,ANN)) %>%
  set_names(c("Chr","Position","rs_ID","Ref","Alt","coverage","coverage_ref",
              "coverage_alt","genotype","HGNC_symbol")) %>% #set colnames
  mutate(Chr=paste("chr",.$Chr,sep="")) #correct chr names

variants<-data.frame(anno$vcf[,c("CHROM","POS","avsnp147","REF","ALT")], 
                     cov, anno$vcf[,c("Gene.refGene")]) %>%
  filter(.[,10] %in% genes) %>%
  set_names(c("Chr","Position","rs_ID","Ref","Alt","coverage","coverage_ref",
              "coverage_alt","genotype","HGNC_symbol")) %>% #set colnames
  mutate(Chr=paste("chr",.$Chr,sep="")) %>% #correct chr names
  rbind(variants_tmp,.) %>%
  unique()

rm(variants_tmp)

##correct coordinates for indels
variants$Position<-ifelse(nchar(as.character(variants$Ref))>1,
                          as.numeric(as.character(variants$Position))+1,
                          ifelse(nchar(as.character(variants$Alt))>1,
                                 as.numeric(as.character(variants$Position))+1,
                                 as.numeric(as.character(variants$Position))))

#______________________________________________
# get transcript annotation from refSeqGenes file
#______________________________________________
##variants within transcripts
transcripts<-left_join(variants,refSeqGenes) %>%
  filter(as.numeric(Start) <= as.numeric(Position) &
           as.numeric(Position) <= as.numeric(End)) %>%
  select(Chr, Position, rs_ID, Ref, Alt, coverage, coverage_ref, coverage_alt, genotype,
         HGNC_symbol, Rank.Exons.Introns, Strand, ENSGene, ENSTranscript, refSeq_mRNA,
         refSeq_protein)

##variants out of transcripts
inter_transcripts<-left_join(anti_join(variants,transcripts),refSeqGenes)%>%
  select(Chr, Position, rs_ID, Ref, Alt, coverage, coverage_ref, coverage_alt, genotype,
         HGNC_symbol, Rank.Exons.Introns, Strand, ENSGene, ENSTranscript, refSeq_mRNA,
         refSeq_protein) %>%
  mutate(Rank.Exons.Introns="intergenic") %>%
  unique()

##add both and sort by chromosome position
trans_anno<-rbind(transcripts,inter_transcripts) %>%
  arrange(Chr, Position, HGNC_symbol)

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
# HSF
#_________________________________________
HSF<-data.frame(clinvar, HSF="error")

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
  Genomiser<-data.frame(UMD,Genomiser=genomiser(UMD))
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
  mutate(Chr=paste("chr",.$Chr,sep=""),
         Position=ifelse(nchar(Ref)>1,
                         as.numeric(as.character(Position))+1,
                         ifelse(nchar(Alt)>1,
                                as.numeric(as.character(Position))+1,
                                as.numeric(as.character(Position))))) %>%
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
  mutate(Position=ifelse(nchar(as.character(Ref))>1,
                         as.numeric(Position)+1,
                         ifelse(nchar(as.character(Alt))>1,
                                as.numeric(Position)+1,
                                as.numeric(Position)))) %>%
  join(predictions,.) %>%
  select(-N_samples) %>%
  unique()

GMfreq[is.na(GMfreq$HGVS_c),"HGVS_c"]<-"."
#______________________________________________
# Preparing output in excel
#______________________________________________
novenove<-GMfreq %>%
  filter(str_detect(HGVS_c,"[\\*\\+\\-][0-9]{1,2}[diATGC]")==TRUE | #all variants up to 99bp from SS
           str_detect(HGVS_c,"[\\*\\+\\-]")==FALSE) #add exonic

non_novenove<-GMfreq %>%
  filter(str_detect(HGVS_c,"[\\*\\+\\-][0-9]{1,2}[diATGC]")==FALSE) %>% #all variants beyond 99bp from SS
  filter(str_detect(HGVS_c,"[\\*\\+\\-]")==TRUE) #remove exonic

cinquenta<-novenove %>%
  filter(str_detect(HGVS_c,"[\\*\\+\\-][5-9][0-9][diATGC]")==FALSE)

off<-novenove %>%
  filter(str_detect(HGVS_c,"[\\*\\+\\-][5-9][0-9][diATGC]")==TRUE) %>%
  rbind(non_novenove)

cov15<-cinquenta %>%
  filter(as.numeric(as.character(coverage)) > 14)

covless15<-cinquenta %>%
  filter(as.numeric(as.character(coverage))<=14)

#______________________________________________
# Writing table
#______________________________________________
# directory creation
dir.create("./my_analysis", recursive=F, mode = "0777")
library("WriteXLS")
write.table(GMfreq,paste("./my_analysis/",sample,"_all_variants.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)

WriteXLS(c("cov15","covless15","off"), ExcelFileName = paste("./my_analysis/",sample,"_genes_variants.xlsx",sep=""), SheetNames = c("High confidence","Low confidence","Off target"), perl = "perl",
         row.names = FALSE, col.names = TRUE,
         envir = parent.frame())