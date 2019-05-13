MHQuant.prep3 <- function(data2){
  # this function reads a .txt file containing aligned DNA sequences
  # and prepares them for MH analysis using the MHQuant function


  l1 <- DNAString(as.character(data2$seq1A))
  l2 <- DNAString(as.character(data2$seq1B))
  r1 <- DNAString(as.character(data2$seq2A))
  r2 <- DNAString(as.character(data2$seq2B))

  MHQuant_seqs <- list(seq1A=l1,
                       seq1B=l2,
                       seq2A=r1,
                       seq2B=r2,
                       MHSeqSize=as.numeric(nchar(as.character(data2$seq1A))))

  return(MHQuant_seqs)

}
