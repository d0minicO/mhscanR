gcq_Sanger <- function(input, MH, equalTo=F, expected){
  # this cool function takes the output of mhq SANGER and tests using Chi Square if there is an enrichment for GC base pairing
  # input is dataframe from mhq CRISPResso = F ## will wrap in another function that splits based on CRISPResso like mhq does
  # MH is the amount of MH we wish to use as the threshold to look for GC quantification in
  # expected is the percentage (represented as value between 0-1, e.g. for an expected 50% GC content "expected=0.5")
  # this should be calculated elsewhere and will depend on each situation (e.g. take the background GC content over the region of where the deletions were found)
  # if equalTo=T then amount of MH EQUAL to number provided will be given (defaults to FALSE)
  # if equalTo=F (default) then amount of MH greater than or equal to the number provided will be given

  # e.g. we want to consider MHs of 1 bp only (MH=1, equalTo=T)
  # e.g. we want to consider all MHs of 2 bp or more (MH=2, equalTo=F)
  # e.g. we want to consider MHs of 2 bp only (MH=2, equalTo=T)

  # output is a dataframe with cols:
  # 1 baseType = GC bases or AT bases
  # 2 baseNum = number of bases of each type (in alleles with the given amount of microhomology)
  # 3 baseProb = observed number of bases of each type in the microhomologies analysed (0 to 1 a.k.a 0 to 100%)
  # 4 expectedProb = expected probability (0 to 1) of bases of each type (known background for the region of the deletions - determined by the user)
  # 5 pval = chi sqaure test p value (chance of finding the observed vs expected probability)

  #' @importFrom stats chisq.test
  #' @importFrom Biostrings letterFrequency
  NULL

  if(equalTo==T){
    # filter the alleles based on the given microhomology amount
    input$MH_filt <- input$MH_amount==MH
    input_MH_filt <- input[input$MH_filt==T,]
  } else if (equalTo==F){
    # filter the alleles based on the given microhomology amount
    input$MH_filt <- input$MH_amount>=MH
    input_MH_filt <- input[input$MH_filt==T,]
  } else {
    cat("You need to chose whether to search for GC content in alleles with microhomology *greater than* or *equal to* the given MH value")
    break
    }

  # calculate percentage of total alleles with that amount of microhomology
  Perc_input_MH_filt <-(nrow(input_MH_filt)/nrow(input))*100

  # collapse the bases of MH vector to allow GC quantification
  MH_bases <- paste(input_MH_filt$MH_sequence, sep="", collapse="")
  MH_bases_num <- nchar(MH_bases)

  # make into DNA string object for GC quantification
  MH_DNA <- DNAString(MH_bases)

  # biostrings to calculate probability and number of C or G and A or T bases
  MH_DNA_GCnum <- letterFrequency(MH_DNA, letters="CG", OR="|", as.prob=F)
  MH_DNA_ATnum <- letterFrequency(MH_DNA, letters="AT", OR="|", as.prob=F)

  # format the data for export and for chi square test
  bases <- c("GC","AT")
  givenSeqs_MH <- c(MH_DNA_GCnum,MH_DNA_ATnum)
  Probs <- c(expected,(1-expected))

  # now the chi square test using expected
  chi <- chisq.test(givenSeqs_MH,p=Probs)
  pval <- format.pval(chi$p.value)

  # put together and output a dataframe
  output <- data.frame(baseType=bases,
                       baseNum=givenSeqs_MH,
                       baseProb=givenSeqs_MH/sum(givenSeqs_MH),
                       expectedProb=Probs,
                       pval=pval)

  return(output)

}
