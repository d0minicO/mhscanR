
is.mutated <- function(a, b, exclude = "-"){
  # function to compare mutations between two sequences a and b
  # adapted from https://www.r-bloggers.com/extract-different-characters-between-two-strings-of-equal-length/
  # if there is a mutation the function returns TRUE
  # if there is no mutation the function returns FALSE
  
  # first convert DNAString objects to characters for simple characters for matching
  a <- as.character(a)
  b <- as.character(b)
  
  if(nchar(a)!=nchar(b)) stop("Lengths of input strings differ. Please check your input.")
  split_seqs <- strsplit(c(a, b), split = "")
  only.diff <- (split_seqs[[1]] != split_seqs[[2]])
  only.diff[
    (split_seqs[[1]] %in% exclude) |
      (split_seqs[[2]] %in% exclude)
    ] <- NA
  diff.info<-data.frame(which(is.na(only.diff)|only.diff),
                        split_seqs[[1]][only.diff],split_seqs[[2]][only.diff])
  names(diff.info)<-c("position","poly.seq.a","poly.seq.b")
  is.mut.seq <- paste(diff.info$poly.seq.a, sep="", collapse="")
  is.mut <- is.mut.seq!=""
  return(is.mut)
}