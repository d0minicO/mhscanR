MHQuant.prep1 <- function(index_out){
  # this function outputs the max allowable sequence to be used for MH quantification (MHSeqSize)
  # without going off the end of the read or exceeding the size of the deletion
  # it takes one argument which is the variable outputted from the index function
  # outputs MHSeqSize
  
  Seq <- index_out$delSeq                                 # deletion sequence
  Ref <- index_out$refSeq                                 # reference sequence
  Del_index_L <- min(index_out$delInd)                    # left hand deletion index
  Del_index_R <- max(index_out$delInd)                    # right hand deletion index
  sizeOfDel <- length(index_out$delInd)                   # size of the deletion
  L <- subseq(Seq, start=1, end=Del_index_L-1)            # Left hand deleted sequence
  R <- subseq(Seq, start=Del_index_R+1, end=nchar(Seq))   # Right hand deleted sequence
  
  # boolean for whether the deletion is too long to pull the entire sequence
  # from the remaining left or right sequences
  longDel <- (sizeOfDel >= nchar(L) | sizeOfDel >= nchar(R))
  shortSizeSeq <- (min(nchar(L), nchar(R))-1) # size of the short sequence to pull
  fullSizeSeq <- sizeOfDel # size of the long sequence to pull
  
  if (longDel){    # check whether to take the shorter size or longer size sequence for MH
    MHSeqSize <- shortSizeSeq
  } else {
    MHSeqSize <- fullSizeSeq
  }
  return(MHSeqSize)
}
