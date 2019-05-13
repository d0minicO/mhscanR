gather2 <- function(file){
  # Gathers the other .text file (not CRISPResso)
  # returns a cleaned up data frame for MH quantification
  dfx <- data.frame(read.delim(file, header=F, sep = "\t"))
  dfx <- cbind(dfx, MHSeqSize="")
  colnames(dfx) <- c("seqName", "seq1A", "seq1B", "seq2A", "seq2B", "MHSeqSize")
  
  return(dfx)
}