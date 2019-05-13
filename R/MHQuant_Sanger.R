MHQuant_Sanger <- function(file){
  # this function quantifies microhomologies in any sequence data
  # that were aligned using software other than CRISPresso
  # e.g. larger deletions mapped with MUSCLE
  # argument is "file", a .txt file
  # containing the DNA sequences directly abutting the deletion
  # per the example below where
  #
  #  +  indicates sequences present in the deleted allele
  #  -  indicates sequences lost in the deleted allele, and
  # ... indicates variable size of the deletion
  #
  #       ++++++++++ ---------- ... ---------- +++++++++
  #       CGTGGCGAGG GCTGAGCTAT ... TGTTAGCACA GCTTCTCCA
  #


  data <- gather2(file)
  cat("Looking in the specified .txt file for premapped sequence data...", "\n")
  results_int <- NULL
  cat("Analysing...", "\n")
  for (i in 1:nrow(data)){ # start of for each row loop
    data2 <- data[i,]

    MHQuant_seqs <- MHQuant.prep3(data2)      # prepare the sequences for microhomology analysis
    MHQuant_out <- MHQuant_sub(MHQuant_seqs)  # look for MH in the corresponding left and hand sequences

    MHSeq_out <- MHSeq.return2(MHQuant_out,    # return MH the sequence, if any, or "No_MH"
                               MHQuant_seqs)

    # add the calculations to a new row of the output
    results_int <- rbind(results_int,
                         data.frame(seqName=data2$seqName,
                                    L1=as.character(MHQuant_seqs$seq1A),
                                    L2=as.character(MHQuant_seqs$seq1B),
                                    R1=as.character(MHQuant_seqs$seq2A),
                                    R2=as.character(MHQuant_seqs$seq2B),
                                    MH_amount=MHQuant_out$max_MH,
                                    MH_sequence=MHSeq_out))

  } # end of for each row loop

  Results <- results_int
  cat("Done!")
  return(Results) # output a data frame

}
