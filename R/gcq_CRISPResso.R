gcq_CRISPResso <- function(input, MH, equalTo=F,  expected){
  # this cool function uses the output folder of mhq CRISPResso files. Again, make sure there are no other .txt files inside this directory!
  #
  # and tests using Chi Square if there is an enrichment for GC base pairing
  # input is dataframe from mhq CRISPResso = T ## will wrap in another function that splits based on CRISPResso like mhq does
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

  #' @import stats
  #' @import Biostrings


  # gather the list of analysed files in the list of dataframes names
  files <- names(input)
  cat("Analysing the microhomology data from CRISPResso + mhq output...", "\n")

  Results <- list()
  for (file in files){ # start of for each file loop

    # gather the input
    data <- input[[file]]
    data <- data[-1,4:6]
    colnames(data) <- c("Reads", "MH", "MH_seq"
                         )
    # make sure MH amount is numeric
    data$MH <- as.numeric(as.character(data$MH))

    if(equalTo==T){
      # filter the alleles based on the given microhomology amount
      data$MH_filt <- data$MH==MH
      data_MH_filt <- data[data$MH_filt==T,]
    } else if (equalTo==F){
      # filter the alleles based on the given microhomology amount
      data$MH_filt <- data$MH>=MH
      data_MH_filt <- data[data$MH_filt==T,]
    } else {
      cat("You need to chose whether to search for GC content in alleles with microhomology *greater than* or *equal to* the given MH value")
    }

    # use split.allele function to split a CRISPResso read table based on reads and counts each read as a new allele
    # all the bases of microhomology are then extracted and returned as a DNAString object for GC quantification
    data_MHs <- split.alleles(data_MH_filt)

    # biostrings to calculate probability and numbers of "C|G" or "A|T" bases
    MH_DNA_GCnum <- letterFrequency(data_MHs, letters="CG", OR="|", as.prob=F)
    MH_DNA_ATnum <- letterFrequency(data_MHs, letters="AT", OR="|", as.prob=F)

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

    cat("Microhomology analysis complete for", file, "moving on to next...", "\n")          # tell us which file is now being analysed
    Results[[file]] <- output

  } # end of for each file loop

  names(Results) <- files # rename the list entries to all the files
  cat("Done!", "\n")
  return(Results)

}
