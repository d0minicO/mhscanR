
MHQuant_sub <- function(MHQuant_seqs){
  # this is a wrapper function for the two seperate MHQuant functions 
  # that quantify MHs in the two possible orientations for
  # MHs to occur either to the left or right of the deletion breakpoint
  
  # do the individual MH scoring and find the best MH out of both options
  MH_A <- MHQuant_L(MHQuant_seqs)
  MH_B <- MHQuant_R(MHQuant_seqs)
  max_MH <- max(MH_A, MH_B)  
  MHQuant_out <- list(max_MH=max_MH,
                      MH_A=MH_A,
                      MH_B=MH_B)
  return(MHQuant_out)
  
}