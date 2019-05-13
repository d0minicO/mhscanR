
index <- function(aligned, reference){
  # this function works on the aligned and reference sequences from each row of the crispresso text file
  # it uses the "-" characters used by crispresso to map the deletion
  # and outputs useful variables for later use
  # including index of the deletion, and deletion and reference sequence
  
  # using the "-" characters used be CRISPResso to define deleted regions and only keep alleles with a simple deletion
  Del_index1 <- (data.frame(stri_locate_all(pattern = '-', aligned, fixed = TRUE))[,1])
  Del_index2 <- split(Del_index1, cumsum(seq_along(Del_index1) %in% (which(diff(Del_index1)>1)+1)))
  
  # Trim and keep only the useful deletion DNA sequence
  I <- max(data.frame(Del_index2[1]))
  II <- (max(data.frame(Del_index2[3]))-min(data.frame(Del_index2[3])))
  DNA_1 <- stri_sub(aligned, I+1)
  DNA_2 <- substr(DNA_1, 1, nchar(DNA_1)-(II+1))
  
  # also trim and get the corresponding reference sequence
  ref1 <- stri_sub(reference, I+1)
  Ref <- substr(ref1, 1, nchar(ref1)-(II+1))
  
  # map the deletions
  Del_index3 <- data.frame(stri_locate_all(pattern = '-', DNA_2, fixed = TRUE))
  Del_index4 <- Del_index3[,1]
  
  index_out <- list(delSeq = DNA_2,
                    refSeq = Ref,
                    delInd = Del_index4)
  return(index_out)
}