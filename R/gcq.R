gcq <- function(input, MH, equalTo=F, expected, CRISPResso=T){
  #' gcq
  #' 
  #' This function is for calculating GC content in CRISPR/Cas9 deletions that have already been analysed using the mhq function. 
  #' It filters the deletions for those with a given amount of microhomology and then performs 
  #' a chi sqaure test to compare the observed vs the expected GC content of the microhomologies.
  #' 
  #' \cr 
  #' Output is a dataframe with columns:
  #' 
  #' \itemize{
  #'   \item baseType = GC bases or AT bases
  #'   \item baseNum = number of bases of each type (in alleles with the given amount of microhomology)
  #'   \item baseProb = observed number of bases of each type in the microhomologies analysed (0 to 1 a.k.a 0 to 100%)
  #'   \item expectedProb = expected probability (0 to 1) of bases of each type (known background for the region of the deletions - determined by the user)
  #'   \item pval = chi sqaure test p value (chance of finding the observed vs expected probability)
  #' }
  #' 
  #' 
  #' 
  #' @examples
  #' gcq(mhqOutCRISPResso, MH=1, equalTo=F, expected=0.46, CRISPResso=T)
  #' gcq(mhqOutSanger, MH=2, equalTo=T, expected=0.51, CRISPResso=F)
  #' 
  #' @param input dataframe output after running mhq(yourData)
  #' @param MH microhomology amount to filter for
  #' @param equalTo if set to TRUE search ONLY for microhomologies equal to MH. If set to FALSE search for microhomologies greater than or equal to MH
  #' @param expected Eexpected background GC content over the region containing deletions (if 50% background, expected=0.5.) Determined by the user.
  #' @param CRISPResso are you analysing a dataframe containing analysed CRISPResso data?
  #' 
  #' @export


  if (CRISPResso==T){
    Results <- gcq_CRISPResso(input, MH, equalTo, expected)
  } else if (CRISPResso==F){
    Results <- gcq_Sanger(input, MH, equalTo, expected)
  } else {
    cat("Define whether you are analysing CRISPResso data or Sanger / Other data by setting CRISPResso=TRUE or FALSE")
  }
  return(Results)

}