is.mut.prep <- function(index_out){
  # this function prepares the deletion and reference sequences
  # for later use by function is.mutated
  # argument is index_out from the index function
  # it outputs 2 variables into seq_prep_out
  
  # 1 L hand del sequence +/- 10bp
  # 2 R hand del sequence +/- 10bp
  
  delInd <- index_out$delInd
  Del_index_L <- min(delInd)
  Del_index_R <- max(delInd)
  
  Seq <- DNAString(index_out$delSeq)
  Ref <- DNAString(index_out$refSeq)
  
  Del_test <- c(subseq(Seq, start=Del_index_L-10, end=Del_index_L-1),subseq(Seq, start=Del_index_R+1, end=Del_index_R+10))
  Ref_test <- c(subseq(DNAString(Ref), start=Del_index_L-10, end=Del_index_L-1),subseq(DNAString(Ref), start=Del_index_R+1, end=Del_index_R+10))
  
  mut_prep_out <- list(Del_test=Del_test,
                       Ref_test=Ref_test)
  return(mut_prep_out)
}