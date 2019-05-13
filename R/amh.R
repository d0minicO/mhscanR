
amh <- function(input, genome){
  #' amh
  #'
  #' This function is for counting the number of "alternative / more proximal" microhomologies that were bypassed during repair the repair of a CRISPR deletion.
  #' It should be used on Sanger sequenced alleles already analysed by running mhq(input, CRISPResso=F). In the other scenario, When analysing CRISPResso output data, mhq automatically counts alternative microhomologies so this script is not necessary.
  #'
  #' \cr
  #' Input is a tab seperated text file (without column names) containing deletions in the same format as output by the mhq function but with some additional columns added. Bed file coordinates are best obtained using a tool like UCSC BLAT.
  #' Also requires the desired genome to be preloaded as a BSGenome object
  #'
  #' \itemize{
  #' \item MutantSequence - same as mhq output
  #' \item ReferenceSequence - same as mhq output
  #' \item SizeOfDeletion - same as mhq output
  #' \item NumberOfReads - same as mhq output
  #' \item MH_amount - same as mhq output
  #' \item LD_chr - Additional: chromosome of deletion start
  #' \item LD_start - Additional: bed coordinate start of deletion span
  #' \item LD_stop - Additional: bed coordinate end of deletion span
  #' \item strand - Additional: strand that deletion was mapped to (or strand that microhomology is reported for)
  #' \item sg_start - Additional: bed coordinate start of sgRNAs span (if pairs of sgRNAs used, then start coordinate of 5' most sgRNA)
  #' \item sg_stop - Additional: bed coordinate end of sgRNAs span (if pairs of sgRNAs used, then end coordinate of 3' most sgRNA)
  #' }
  #'
  #' \cr
  #' Output is the original dataframe but with two additional columns:
  #'
  #' \itemize{
  #' \item altMH_count - Number of alternative microhomologies found
  #' \item add_delSize - Additional size of deletion beyond sgRNAs (same as deletion size when using one sgRNA)
  #' }
  #'
  #' @examples
  #' genome <- BSgenome.Mmusculus.UCSC.mm9
  #' amh(input="~/exampleData/Sanger/sequenceDataFile_altMH.txt", genome)
  #'
  #' @param input text file derived from output of mhq(yourData) with additional columns as specified below
  #' @param genome BSgenome object to pull sequences from. IMPORTANT: Has to be loaded before running the function. (i.e. genome <- BSgenome.Mmusculus.UCSC.mm9). Allows any BSgenome supported genome to be used.
  #' @export
  #'
  #'
  #' @import BSgenome
  #' @importFrom Biostrings getSeq countPattern matchPattern reverseComplement DNAString subseq
  #' @importFrom readr read_tsv
  #' @import dplyr
  #' @import utils

  NULL

  df <- read_tsv(input, col_names = F) %>% data.frame()
  colnames(df) <- c("ID", "L1", "L2", "R1", "R2", "MH_score", "MH", "LD_chr", "LD_start", "LD_stop", "strand", "sg_start", "sg_stop")

  # perform deletion classification to determine which sequences to search for MHs in
  df <- classify(df)

  # start the looping over rows of df, creating new list to populate with alternative MH counts + original data
  df_new <- list()

  for (i in 1:nrow(df)){ # start of for each row loop
    df_row <- df[i,]

    if (df_row$MH_score==0){
      df_row$altMH_count <- NULL
      df_row$add_delSize <- df_row$`y-b`+df_row$`a-x` # just calculate additional del size without calculating alternative MH counts for LDs with no MH
      df_new[[i]] <- df_row
    } else if (df_row$class1==T){
      df_new[[i]] <- class1(df_row)
    } else if (df_row$class2==T){
      df_new[[i]] <- class2(df_row)
    } else if (df_row$class3==T){
      df_new[[i]] <- class3(df_row)
    } else if (df_row$class4==T){
      df_new[[i]] <- class4(df_row)
    } else if (df_row$class5==T){
      df_new[[i]] <- class4(df_row)
    }

  } # end of for each row loop
  df_new <- bind_rows(df_new)
  return(df_new)
} # end of amh function

# sub functions
classify <- function(df){
  # this sub function takes a data frame in format used by amh function
  # it determines whether LDs extended in one or both directions past the sgRNA cut sites
  # outputs the same data frame with five new columns that are boolean for
  # for which deletion class we are dealing with

  # this is important in order to know which additional / alternative MHs to count
  # as we are only interested in MHs that lie between the sgRNA cut sites and the real deletion ends
  # not interested in counting MHs that lie between sgRNAs

  # determine if LD is L or R sided by performing crunching where sg coordinates are
  # relative to LD coordinates and what "class" of deletion this is

  df$"y-b" <- (df$LD_stop-df$sg_stop)
  df$"a-x" <- (df$sg_start-df$LD_start)

  # make a T F column for each class of LD to identify them
  df$class1 <- (df$`y-b`>50 & df$`a-x`>50) ## class 1 dels extend past both sgRNAs # these require searching for MH pattern in two sequences, one L and one R of sgRNAs
  df$class2 <- (df$`y-b`>50 & df$`a-x`<(-50)) ## class 2 dels extend only past downstream (right hand side / 3') sgRNA NOT removing sequence between sgRNAs
  df$class3 <- (df$`y-b`>50 & df$`a-x`< 50 & df$`a-x` > (-50)) ## class 3 dels extend only to downstream (right hand side / 3') sgRNA including sequence between sgRNAs within deletion
  df$class4 <- (df$`y-b`<(-50)) ## class 4 dels extend only past upstream (left hand side / 5') sgRNA NOT removing sequence between sgRNAs
  df$class5 <- (df$`y-b`< 50 &df$`y-b` > (-50)) ## class 5 dels extend only to upstream (left hand side / 5') sgRNA including sequence between sgRNAs within deletion

  return(df)

}

# Each class requires a different sequence to search for alternative MHs within
# and so gets its own sub function

class1 <- function(df_row){
  # this function searches for alternate MHs in class 1 LDs that extend past both sgRNAs

  if(df_row$strand == 1){
    seqa <- as.character(Biostrings::getSeq(genome, df_row$LD_chr, df_row$LD_start, df_row$sg_start)) # -1 from start coordinate otherwise 0-based indexing means original identified MH might not be found!
    seqb <- as.character(Biostrings::getSeq(genome, df_row$LD_chr, df_row$sg_stop, df_row$LD_stop))

    df_row$altMH_count <- (Biostrings::countPattern(df_row$MH,seqa)+Biostrings::countPattern(df_row$MH,seqb))-1 # minus one that is the original MH identified
    df_row$add_delSize <- nchar(as.character(paste0(seqa, seqb), collapse=""))

  } else if(df_row$strand == 0){
    seqa <- as.character(Biostrings::getSeq(genome, df_row$LD_chr, df_row$LD_start, df_row$sg_start))
    seqb <- as.character(Biostrings::getSeq(genome, df_row$LD_chr, df_row$sg_stop, df_row$LD_stop))

    df_row$altMH_count <- (Biostrings::countPattern(reverseComplement(DNAString(df_row$MH)),seqa)+Biostrings::countPattern(reverseComplement(DNAString(df_row$MH)),seqb))-1 # minus one that is the original MH identified
    df_row$add_delSize <- nchar(as.character(paste0(seqa, seqb), collapse=""))
  }
  return(df_row)
}

class2 <- function(df_row){
  # this function searches for MHs in class 2 LDs that extend past 3' sgRNA and does not include both sgRNAs in deletion

  if(df_row$strand== 1){
    seqa <- as.character(Biostrings::getSeq(genome, df_row$LD_chr, df_row$LD_start, df_row$LD_stop))
    df_row$altMH_count <- Biostrings::countPattern(df_row$MH,seqa)-1 # minus one that is the original MH identified
    df_row$add_delSize <- nchar(as.character(paste0(seqa), collapse=""))
  }else if(df_row$strand == 0){
    seqa <- as.character(Biostrings::getSeq(genome, df_row$LD_chr, df_row$LD_start, df_row$LD_stop))
    df_row$altMH_count <- Biostrings::countPattern(reverseComplement(DNAString(df_row$MH)),seqa)-1 # minus one that is the original MH identified
    df_row$add_delSize <- nchar(as.character(paste0(seqa), collapse=""))
  }
  return(df_row)
}

class3 <- function(df_row){
  # this function searches for MHs in class 3 LDs that extend past 3' sgRNA and includes both sgRNAs in deletion

  if(df_row$strand == 1){
    seqa <- as.character(Biostrings::getSeq(genome, df_row$LD_chr, df_row$sg_stop, df_row$LD_stop))
    df_row$altMH_count <- Biostrings::countPattern(df_row$MH,seqa) # DO NOT minus one because we did not scan over the part of the deletion that contained the original MH because it was contained within the sgRNA
    df_row$add_delSize <- nchar(as.character(paste0(seqa), collapse=""))
  }else if(df_row$strand == 0){
    seqa <- as.character(Biostrings::getSeq(genome, df_row$LD_chr, df_row$sg_stop, df_row$LD_stop))
    df_row$altMH_count <- Biostrings::countPattern(reverseComplement(DNAString(df_row$MH)),seqa) # DO NOT minus one because we did not scan over the part of the deletion that contained the original MH because it was contained within the sgRNA
    df_row$add_delSize <- nchar(as.character(paste0(seqa), collapse=""))
  }

  return(df_row)
}

class4 <- function(df_row){
  # this function searches for MHs in class 4 LDs that extend past 5' sgRNA and does not include both sgRNAs in deletion

  if(df_row$strand == 1){
    seqa <- as.character(Biostrings::getSeq(genome, df_row$LD_chr, df_row$LD_start, df_row$LD_stop))
    df_row$altMH_count <- Biostrings::countPattern(df_row$MH,seqa)-1 # minus one that is the original MH identified
    df_row$add_delSize <- nchar(as.character(paste0(seqa), collapse=""))
  }else if(df_row$strand == 0){
    seqa <- as.character(Biostrings::getSeq(genome, df_row$LD_chr, df_row$LD_start, df_row$LD_stop))
    df_row$altMH_count <- Biostrings::countPattern(reverseComplement(DNAString(df_row$MH)),seqa)-1 # minus one that is the original MH identified
    df_row$add_delSize <- nchar(as.character(paste0(seqa), collapse=""))
  }
  return(df_row)
}

class5 <- function(df_row){
  # this function searches for MHs in class 5 LDs that extend past 5' sgRNA and includes both sgRNAs in deletion

  if(df_row$strand == 1){
    seqa <- as.character(Biostrings::getSeq(genome, df_row$LD_chr, df_row$LD_start, df_row$sg_start))
    df_row$altMH_count <- Biostrings::countPattern(df_row$MH,seqa) # DO NOT minus one because we did not scan over the part of the deletion that contained the original MH because it was contained within the sgRNA
    df_row$add_delSize <- nchar(as.character(paste0(seqa), collapse=""))
  }else if(df_row$strand == 0){
    seqa <- as.character(Biostrings::getSeq(genome, df_row$LD_chr, df_row$LD_start, df_row$sg_start))
    df_row$altMH_count <- Biostrings::countPattern(reverseComplement(DNAString(df_row$MH)),seqa) # DO NOT minus one because we did not scan over the part of the deletion that contained the original MH because it was contained within the sgRNA
    df_row$add_delSize <- nchar(as.character(paste0(seqa), collapse=""))
  }
  return(df_row)
}
