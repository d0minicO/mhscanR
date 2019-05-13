
MHQuant.prep2 <- function(index_out, MHSeqSize){
  # this function extracts the sequences of the length defined in MHQuant.prep1
  # arguments are the output of index (index_out) and MHQuant.prep1 (MHQuant_prep)
  
  Seq <- index_out$delSeq                                 # deletion sequence
  Ref <- index_out$refSeq                                 # reference sequence
  Del_index_L <- min(index_out$delInd)                    # left hand deletion index
  Del_index_R <- max(index_out$delInd)                    # right hand deletion index
  L <- subseq(Seq, start=1, end=Del_index_L-1)            # Left hand deleted sequence
  R <- subseq(Seq, start=Del_index_R+1, end=nchar(Seq))   # Right hand deleted sequence
  
  # extract the sequences that remain intact from the aligned sequence, and for the missing sequence use the reference allele
  seq1A <- DNAString(subseq(Seq, start=Del_index_L-MHSeqSize, end=Del_index_L-1))
  seq1B <- DNAString(subseq(Ref, start=Del_index_L, end=Del_index_L+MHSeqSize-1))
  seq2A <- DNAString(subseq(Ref, start=Del_index_R-MHSeqSize+1, end=Del_index_R))
  seq2B <- DNAString(subseq(Seq, start=Del_index_R+1, end=Del_index_R+MHSeqSize))
  
  MHQuant_seqs <- list(seq1A=seq1A,
                       seq1B=seq1B,
                       seq2A=seq2A,
                       seq2B=seq2B,
                       MHSeqSize=MHSeqSize)
  return(MHQuant_seqs)
}