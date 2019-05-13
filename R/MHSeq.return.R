MHSeq.return <- function(MHQuant_out, index_out){
  # this function returns the bases of the MH sequence found in CRISPresso alleles
  # if it was found
  # takes MHQuant_out from MHQuant function and index_out as arguments
  
  delInd <- index_out$delInd
  Del_index_L <- min(delInd)
  Del_index_R <- max(delInd)
  Seq <- DNAString(index_out$delSeq)
  max_MH <- MHQuant_out$max_MH
  MH_A <- MHQuant_out$MH_A
  MH_B <- MHQuant_out$MH_B
  
  if (max_MH==0) {
    MHSeq <- "No_MH"
  }
  
  if (max_MH>0&max_MH==MH_A){
    MHSeq <-  as.character(subseq(Seq, start=Del_index_L-max_MH, end=Del_index_L-1))
  }
  
  if (max_MH>0&max_MH==MH_B){
    MHSeq <- as.character(subseq(Seq, start=Del_index_R+1, end=(Del_index_R+max_MH)))
  }
  return(MHSeq)
}