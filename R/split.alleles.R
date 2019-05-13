split.alleles <- function(input){
  # this function splits a CRISPResso read table based on reads and counts each read as a new allele
  # all the bases of microhomology are then extracted and returned as a DNAString object for GC quantification
  
  full_seq_MH <- list()
  # this for loop loops over the MHQuant output and produces a vector containing the exact number of MH bases found
  for (i in 1:nrow(input)) {
    # pull the i th row
    seq_MH <- input[i,]
    
    # replicate the MH sequence the muber of times it was read
    times_seq_MH <- rep(as.character(seq_MH$MH_seq), as.numeric(as.character(seq_MH$Reads)))
    
    # collapse the character vector into a single string
    times_seq2_MH <- paste(times_seq_MH, sep="", collapse="")
    # add this string to the row of the full_seq list
    full_seq_MH[[i]] <- times_seq2_MH
  }
  
  # collapse all the MH reads into a single MH string and make into DNAString
  input_MHs <- paste(full_seq_MH, sep="", collapse="")

  # make into DNA string object
  input_MHs <- DNAString(input_MHs)
  
  return(input_MHs)
  
}
