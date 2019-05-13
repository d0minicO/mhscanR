MHQuant_L <- function(MHQuant_seqs){
  # this function looks along the length of the sequences given by MHQuant.prep2
  # from L to R and if there is a match then give a score of + 1, if no match then score is 0 and tries again
  # outputs MH score for left hand MH (MH_A)
  MH_A <- 0
  for(i in 1:MHQuant_seqs$MHSeqSize){
    if (subseq(MHQuant_seqs$seq1A, start=i, end=i)==subseq(MHQuant_seqs$seq2A, start=i, end=i)){
      MH_A <- MH_A + 1
    } else {
      MH_A <- 0
    }
  }
  return(MH_A)
}