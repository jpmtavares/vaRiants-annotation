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
library(reticulate)
use_python("/usr/bin/python3.5")
library(doParallel)
#______________________________________________
# R source files
#______________________________________________
source("./functions/HGVS_mutalyzer_parallel.R")
source("./functions/reverseComplementDNA.R")
source("./functions/UMDpredictor.R")
source("./functions/vcfFORMAT.R")
source("./functions/utils.R")
source("./functions/mutationTaster_parallel.R")
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
gene.list<-"A2ML1, AARS2, ABCC9, ACAD9, ACADVL, ACTA1, ACTC1, ACTN2, AGK, AGL, AGPAT2, AKAP9, ALMS1, ANK2, ANKRD1, ANO5, ATP5E, ATPAF2, BAG3, BRAF, BSCL2, CACNA1C, CACNA1D, CACNA2D1, CACNB2, CALM1, CALM2, CALM3, CALR3, CAPN3, CASQ2, CAV3, CHRM2, COA5, COA6, COL7A1, COQ2, COX15, COX6B1, CRYAB, CSRP3, CTF1, CTNNA3, DES, DLD, DMD, DNAJC19, DOLK, DSC2, DSG2, DSG3, DSP, DTNA, ELAC2, EMD, EYA4, FAH, FHL1, FHL2, FHOD3, FKRP, FKTN, FLNC, FOXD4, FOXRED1, FXN, GAA, GATA4, GATA6, GATAD1, GFM1, GJA1, GJA5, GLA, GLB1, GNPTAB, GPD1L, GUSB, HCN4, HFE, HRAS, ILK, JPH2, JUP, KCND2, KCND3, KCNE1, KCNE2, KCNE3, KCNE5, KCNH2, KCNJ2, KCNJ5, KCNJ8, KCNK17, KCNQ1, KLF10, KRAS, LAMA2, LAMA4, LAMP2, LDB3, LDLR, LIAS, LMNA, LZTR1, MAP2K1, MAP2K2, MIB1, MLYCD, MRPL3, MRPL44, MRPS22, MTO1, MURC, MYBPC3, MYH6, MYH7, MYL2, MYL3, MYLK2, MYOM1, MYOT, MYOZ1, MYOZ2, MYPN, NEBL, NEXN, NF1, NKX2-5, NKX2-6, NNT, NRAS, OBSCN, PDHA1, PDLIM3, PHKA1, PITX2, PKP2, PLN, PMM2, PRDM16, PRKAG2, PSEN1, PSEN2, PTPN11, RAF1, RANGRF, RASA2, RBM20, RIT1, RRAS, RYR2, SCN10A, SCN1B, SCN2B, SCN3B, SCN4B, SCN5A, SCO2, SDHA, SGCA, SGCD, SHOC2, SLC22A5, SLC25A3, SLC25A4, SLMAP, SNTA1, SOS1, SOS2, SPEG, SPRED1, SURF1, SYNE1, SYNE2, TAZ, TBX20, TBX5, TCAP, TGFB3, TMEM43, TMEM70, TMPO, TNNC1, TNNI3, TNNI3K, TNNT2, TOR1AIP1, TPM1, TRDN, TRIM63, TRPM4, TSFM, TTN, TTR, TXNRD2, VCL, XK
"
#g<-paste(as.character(commandArgs(TRUE)[1]),sep="")
g<-unlist(genesofinterest(gene.list)[1])
#genes<-unlist(strsplit(g,"\\|"))
genes<-unlist(genesofinterest(gene.list)[2])

#______________________________________________
# preparing input files
#______________________________________________
load("./sources/bcbio_pipeline.Rdata") ## load previously saved .Rdata with reference transcripts, inHouse variant frequency
## save(refSeqGenes,freq,file="./sources/bcbio_pipeline.Rdata")

#UMD-predictor
umd<-UMDpredictor(genes)

#GM_freq
#GM_freq<-GMfrequency(freq,str_extract(getwd(),"^.*(?=(/github))"))
GM_freq<-freq
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
    mutate(type=mapply(variant_type, HGVS_c, HGVS_p))
  
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
  MT_anno<-MT(HSF)
  #MT_anno<-HSF
  #correct HGVS
  MT_anno$HGVS_p<-mapply(synonymous,MT_anno$HGVS_p)
  
  #correct homozygous counting
  #MT_anno$homozygous_ExAC[grep("\\.",MT_anno$homozygous_ExAC)]<-"-"
  #MT_anno$homozygous_1000G[grep("\\.",MT_anno$homozygous_1000G)]<-"-"
  
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
  # ANNOVAR predictions (dbnsfp35c - onde está Polyphen2_HVAR_score,Polyphen2_HVAR_pred,DANN_score?)
  #_________________________________________
  predictions<-anno$vcf %>%
    select(CHROM,POS,avsnp150,REF,ALT,SIFT_score,SIFT_pred,PROVEAN_score,
           PROVEAN_pred,
           MutationAssessor_score,MutationAssessor_pred,CADD13_PHRED,
           FATHMM_coding,FATHMM_noncoding,GWAVA_region_score,`GERP.._RS`,
           phyloP20way_mammalian,phastCons20way_mammalian,SiPhy_29way_logOdds,
           dpsi_max_tissue,PopFreqMax,X1000G_EUR,X1000G_ALL,gnomAD_genome_NFE,gnomAD_genome_ALL,
           ESP6500siv2_EA,ESP6500siv2_ALL,InterVar_automated) %>%
    set_names(c("Chr","Position","rs_ID","Ref","Alt","SIFT_score","SIFT_pred",
                "PROVEAN_score","PROVEAN_pred",
                "MutationAssessor_score","MutationAssessor_pred",
                "CADD13_PHRED","FATHMM_coding","FATHMM_noncoding",
                "GWAVA_region_score","GERP++RS","phyloP20_mammalian",
                "phastCons20_mammalian","SiPhy","SPIDEX","PopFreqMax","1000G_EUR",
                "1000G_ALL","gnomAD_EUR","gnomAD_ALL","ESP_EUR","ESP_ALL","InterVar")) %>%
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
  #predictions$Polyphen2_HVAR_pred<-ifelse(as.character(predictions$Polyphen2_HVAR_pred)=="D","probably damaging",
  #                                        ifelse(as.character(predictions$Polyphen2_HVAR_pred)=="P","possibly damaging",
  #                                               ifelse(as.character(predictions$Polyphen2_HVAR_pred)=="B","benign",as.character(predictions$Polyphen2_HVAR_pred))))
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
  # InterVar
  #______________________________________________
  intervar<-anno$vcf %>%
    select(CHROM, POS, avsnp150, REF, ALT,
           InterVar_automated, PVS1, PS1, PS2, PS3,
           PS4, PM1, PM2, PM3, PM4, PM5, PM6, PP1,
           PP2, PP3, PP4, PP5, BA1, BS1, BS2, BS3,
           BS4, BP1, BP2, BP3, BP4, BP5, BP6, BP7)
  # get criteria for intervar
  intervar_criteria<-apply(intervar, 1, function(x) paste(names(x[x == 1])[-1], collapse = ", "))
  # join output
  intervar<-data.frame(intervar, intervar_criteria) %>%
    select(CHROM, POS, avsnp150, REF, ALT, InterVar_automated, intervar_criteria) %>%
    mutate(CHROM = paste("chr", CHROM, sep="")) %>%
    set_names(c("Chr","Position","rs_ID","Ref","Alt", "InterVar", "InterVar_criteria")) %>%
    indelsCoordinates(.) %>% # correct indels coordinates
    join(predictions,.) %>%
    unique()
  
  #______________________________________________
  # GM_freq
  #______________________________________________
  GMfreq<-GM_freq %>%
    indelsCoordinates(.) %>% # correct indels coordinates
    join(intervar,.) %>%
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
