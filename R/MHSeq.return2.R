MHSeq.return2 <- function(MHQuant_out, MHQuant_seqs){
  # this function returns the bases of the MH sequence found
  # if it was found
  # it uses the "other sequences" .text file (not CRISPresso)
  # takes MHQuant_out from MHQuant function and 
  # MHQuant_seqs from MHQuant.prep3 as arguments
  
  max_MH <- MHQuant_out$max_MH
  MH_A <- MHQuant_out$MH_A
  MH_B <- MHQuant_out$MH_B
  
  if (max_MH==0) {
    MHSeq <- "No_MH"
  }
  
  if (max_MH>0&max_MH==MH_A){
    MHSeq <-  as.character(subseq(MHQuant_seqs$seq1A, start=MHQuant_seqs$MHSeqSize-max_MH+1, end=MHQuant_seqs$MHSeqSize))
  }
  
  if (max_MH>0&max_MH==MH_B){
    MHSeq <- as.character(subseq(MHQuant_seqs$seq1B, start=1, end=max_MH))
  }
  return(MHSeq)
}