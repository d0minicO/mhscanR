
gather <- function(file){
  # Gathers the files the specified directory
  # adds a buffer of "-" characters at beginning and end to allow indexing of the deletion
  # keeps only the deleted reads
  # returns a list of all the data tables in a list called data_all
  
  dfx <- data.frame(read.delim(file, header=T, sep = "\t"))
  
  dfx <- dfx[dfx$UNMODIFIED=="False"&dfx$n_deleted>0&dfx$n_inserted<1&dfx$n_mutated<1&!grepl("-", dfx$Reference_Sequence, fixed=T),]
  buffer <- as.character("-----")
  dfx$Aligned_Sequence <- paste0(buffer, dfx$Aligned_Sequence, buffer)
  dfx$Reference_Sequence <- paste0(buffer, dfx$Reference_Sequence, buffer)
  return(dfx)
}