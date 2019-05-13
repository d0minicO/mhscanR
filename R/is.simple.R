is.simple <- function(index_out){
  # this function takes the index of the deletion determined by the function index
  # and determines it is a simple deletion or not (i.e. if the deletion spans one contiguos region or not)
  # output is boolean
  runLen <- rle(diff(index_out$delInd))
  delInd <- index_out$delInd
  simple_out <- any(runLen$lengths>=(length(delInd)-1) & runLen$values==1) | length(delInd)==1
  return(simple_out)
}
