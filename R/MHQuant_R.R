MHQuant_R <- function(MHQuant_seqs){
  # this function looks along the length of the sequences given by MHQuant.prep2
  # from L to R and if there is a match then give a score of + 1, if no match then score is 0 and loop breaks
  # outputs MH score for right hand MH (MH_B)
  MH_B <- 0
  for(i in 1:MHQuant_seqs$MHSeqSize){
    if (subseq(MHQuant_seqs$seq1B, start=i, end=i)==subseq(MHQuant_seqs$seq2B, start=i, end=i)){
      
      MH_B <- MH_B + 1
    } else {
      MH_B <- MH_B
      break
    }
  }
  return(MH_B)
}