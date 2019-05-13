altMH <- function(MHSeq_out, index_out){
  # this function takes as arguments the MH sequence to look for
  # from MHSeq.return and index_out from index
  
  # extract the variables to compare
  testSeq <- index_out$delSeq
  refSeq <- index_out$refSeq
  MH <- MHSeq_out
  Del_index <- index_out$delInd
  
  # L and R read position end points
  Del_index_L <- min(index_out$delInd)
  Del_index_R <- max(index_out$delInd)
  
  # convert and get the DNA sequences using biostrings
  # subset the sequence to get deleted sequence from Reference
  refDel <- subseq(DNAString(refSeq), start=Del_index_L, end=Del_index_R)
  
  # count the MH pattern in the del sequence # subtract 1 to account for its presence once as the MH quantified initially
  altMH_out <- countPattern(MH,refDel)-1
  return(altMH_out)
}