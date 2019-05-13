
mhq <- function(input, CRISPResso=TRUE){
  #' mhq
  #'
  #' This function is for quantifying microhomologies in pre-mapped DNA sequencing of CRISPR/Cas9 deletions.
  #' It accepts targeted amplicon next-generation sequencing data already analysed using CRISPResso (1.0.x) (Pinello et al., 2016).
  #' Alternatively, it accepts any sequencing data (eg Sanger sequencing) that has already been processed using a local alignment tool (e.g. MUSCLE).
  #'
  #' When analysing CRISPResso data (CRISPResso=TRUE): Input directory should contain ONLY one or more Allele_frequency_table.txt files that were outputted by CRISPResso (1.0.x) (Pinello et al., 2016). No other .txt files should be in this directory.
  #' Output is a data frame (or list of data frames if multiple files are being analysed). \cr
  #' \cr
  #' Also writes new .txt files to a new directory (/path/to/directory/MHQuant_out/MHQuant_out_Allele_Frequency_table.txt). \cr
  #' \cr
  #' Output contains these tab seperated columns: \cr
  #' \itemize{
  #'   \item MutantSequence - already determined by CRISPResso
  #'   \item ReferenceSequence - already determined by CRISPResso
  #'   \item SizeOfDeletion - already determined by CRISPResso
  #'   \item NumberOfReads - already determined by CRISPResso
  #'   \item MH_amount - microhomology amount found, if any, or 0
  #'   \item MH_sequence - microhomology sequence found, if any, or "No_MH"
  #'   \item altMH_count - alternative microhomologies found within the deleted sequence, if any, or "No_MH"
  #' }
  #' \cr
  #' \cr
  #' \cr
  #' When analysing Sanger sequencing data (CRISPResso=FALSE): Input should be a single .txt file containing deletion alleles that were already aligned using a local alignment tool such as MUSCLE (Edgar, 2004). Input text file should contain five tab-seperated columns containing the DNA sequences of the breakpoints to be analysed. The length of each of the DNA sequences can be varied between rows, but must be the same between columns. See this example: \cr
  #' \cr
  #' Example1   CGTGGCGAGG GCTGAGCTAT TGTTAGCACA GCTTCTCCA \cr
  #' Example2   CGTGGCGAGGCGTGG GCTGAGCTATGCTAT TGTTAGCACAGCACA GCTTCTCCACTCCA \cr
  #'
  #' \itemize{
  #'   \item column1 - Sequence name
  #'   \item column2 - 5' sequence not included in the deletion
  #'   \item column3 - 5' sequence included in the deletion
  #'   \item column4 - 3' sequence included in the deletion
  #'   \item column5 - 3' sequence not included in the deletion
  #' }
  #'
  #'
  #' @param input path to directory containing Allele_Frequency_tables.txt files or path to sequenceDataFile.txt
  #' @param CRISPResso Are you analysing CRISPResso Allele_frequency_table.txt files? (=TRUE/FALSE)
  #'
  #' @keywords CRISPR, Microhomology
  #' @export
  #'
  #' @examples
  #' mhq(input="~/exampleData/CRISPResso/", CRSISPresso=TRUE)
  #' mhq(input="~/exampleData/Sanger/sequenceDataFile.txt", CRSISPresso=FALSE)
  #'
  #' @importFrom Biostrings getSeq countPattern matchPattern reverseComplement
  #' @import stringi
  #' @import stringr
  NULL

  if (CRISPResso==T){
    directory <- input
    Results <- MHQuant_CRISPResso(directory)
  } else if (CRISPResso==F){
    file <- input
    Results <- MHQuant_Sanger(file)
  }
  return(Results)
}
