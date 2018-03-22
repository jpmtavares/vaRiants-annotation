##########################################################
#                                                        #
#                variants in Transcripts                 #
#                                                        #
##########################################################
variants_transcripts<-function(variants_genes,refSeqGenes){
  ##variants within transcripts
  transcripts<-left_join(variants_genes,refSeqGenes) %>%
    filter(as.numeric(as.character(Start)) <= as.numeric(as.character(Position)) &
             as.numeric(as.character(Position)) <= as.numeric(as.character(End))) %>%
    select(Chr, Position, rs_ID, Ref, Alt, coverage, coverage_ref, coverage_alt, genotype,
           HGNC_symbol, Rank.Exons.Introns, Strand, ENSGene, ENSTranscript, refSeq_mRNA,
           refSeq_protein)
  
  ##variants out of transcripts
  inter_transcripts<-left_join(anti_join(variants_genes,transcripts),refSeqGenes)%>%
    select(Chr, Position, rs_ID, Ref, Alt, coverage, coverage_ref, coverage_alt, genotype,
           HGNC_symbol, Rank.Exons.Introns, Strand, ENSGene, ENSTranscript, refSeq_mRNA,
           refSeq_protein) %>%
    mutate(Rank.Exons.Introns="intergenic") %>%
    unique()
  
  ##add both and sort by chromosome position
  trans_anno<-rbind(transcripts,inter_transcripts) %>%
    arrange(Chr, Position, HGNC_symbol) %>%
    unique
  
  return(trans_anno)
}