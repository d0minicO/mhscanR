# mhscanR
#### R package for quantifying microhomologies in pre-mapped DNA sequencing of CRISPR/Cas9 deletions

mhscanR contains three main functions for analysing microhomology sequences in CRISPR/Cas9 deletions.

> **mhq**

Quantifies micorhomologies in CRISPR/Cas9 deletions. Accepts targeted amplicon next-generation sequencing data analysed using CRISPResso (1.0.x, [Pinello et al., 2016](doi.org/10.1038/nbt.3583)), or accepts any other sequencing data (eg Sanger sequencing) processed using a local alignment tool (e.g. [MUSCLE](https://doi.org/10.1186/1471-2105-5-113)).

> **gcq** 

Quantifies GC content of microhomologies of different lengths. Performs statistical test of whether GC bases are enriched compared to an expected background GC content.

> **amh**

Counts alternative microhomologies in Sanger sequencing data of CRISPR/Cas9 deletions analysed using mhq. Note: when analysing CRISPResso data, alternative microhomology count is automatically calculated by the mhq function and this function does not need to be run separately.

* * *

### Example 1 -- Analysing microhomologies in deep sequenced CRISPR/Cas9 deletions analysed using CRISPResso
```
Results <- mhq(input=/path/to/directory/, CRSISPresso=TRUE)
```
- input=/path/to/directory/ containing **ONLY** one or more Allele_frequency_table.txt files that were outputted by CRISPResso (1.0.x) ([Pinello et al., 2016](doi.org/10.1038/nbt.3583)). No other .txt files should be in this directory.
- CRISPResso=TRUE

Output is a dataframe (or list of dataframes) containing tab seperated columns:
- MutantSequence (already determined by CRISPresso)
- ReferenceSequence (already determined by CRISPresso)
- SizeOfDeletion (already determined by CRISPresso)
- NumberOfReads (already determined by CRISPresso)
- MH_amount (microhomology amount found, if any, or "No_MH")
- MH_sequence (microhomology sequence found, if any, or "No_MH")
- altMH_count (alternative microhomologies found within the deleted sequence)

Writes analysed data files in /path/to/directory/MHQuant_out/MHQuant_out_Allele_Frequency_table.txt 
* * *
### Example 2 -- Analysing microhomologies in Sanger sequenced CRISPR/Cas9 deletions analysed using MUSCLE
```
Results <- mhq(input=/path/to/sequenceDataFile.txt, CRSISPresso=FALSE)
```
- Input= sequenceDataFile.txt
- CRISPResso=FALSE

 sequenceDataFile.txt should contain deletion alleles that were already aligned using a local alignment tool such as [MUSCLE](https://doi.org/10.1186/1471-2105-5-113). It needs five tab-seperated columns containing the DNA sequences of the breakpoints to be analysed. The length of each of the DNA sequences can be varied between rows, but must be the same between columns. For example:

```
Example1   CGTGGCGAGG GCTGAGCTAT TGTTAGCACA GCTTCTCCA
Example2   CGTGGCGAGGCGTGG GCTGAGCTATGCTAT TGTTAGCACAGCACA GCTTCTCCACTCCA
```
where:
- column1 = Sequence name
- column2 = 5' sequence **not included** in the deletion
- column3 = 5' sequence **included** in the deletion
- column4 = 3' sequence **included** in the deletion
- column5 = 3' sequence **not included** in the deletion

Output is a dataframe containing the original five columns as well as these additional columns:
- MH_amount (microhomology amount found, if any, or "No_MH")
- MH_sequence (microhomology sequence found, if any, or "No_MH")
* * *
### Example 3 -- Analysing GC content of microhomologies in deep sequenced CRISPR/Cas9 deletions or Sanger sequenced deletions analysed using mhq
```
GC_2orMoreBp_MHs <- gcq(mhqOutCRISPResso, MH=2, equalTo=F, expected=0.46, CRISPResso=T)
GC_3bp_MHs <- gcq(mhqOutSanger, MH=3, equalTo=T, expected=0.56, CRISPResso=F)
```

- input=list of dataframes output after running mhq(yourData, CRISPResso=T), or dataframe output after running mhq(yourData, CRISPResso=F)
- MH=Integer microhomology amount to filter for
- equalTo=TRUE or FALSE. If TRUE, filter for MHs of length=MH. If FALSE, filter for MHs of length>=MH
- expected=Expected background GC content over the region containing deletions. Determined by the user.
- CRISPResso=TRUE or FALSE. Are you analysing a dataframe containing analysed CRISPResso data or other / Sanger sequenced data?

Output is a dataframe with columns:
- baseType = GC bases or AT bases
- baseNum = number of bases of each type (in alleles with the given amount of microhomology)
- baseProb = observed number of bases of each type in the microhomologies analysed (0 to 1 a.k.a 0 to 100
- expectedProb = expected probability (0 to 1) of bases of each type (known background for the region of the deletions - determined by the user)
- pval = chi square test p value (chance of finding the observed vs expected probability)
* * *
### Example 4 -- Analysing alternative microhomologies in Sanger sequenced CRISPR/Cas9 deletions

```
library(BSgenome.Mmusculus.UCSC.mm9)
genome <- BSgenome.Mmusculus.UCSC.mm9
alt_MHs <- amh(input, genome)
```

- input= tab seperated text file containing data derived from output of mhq(yourData) **with additional columns** added (see below)
- genome= BSGenome object loaded before running the function.

Input text file columns required:
- MutantSequence - from mhq output
- ReferenceSequence - from mhq output
- SizeOfDeletion - from mhq output
- NumberOfReads - from mhq output
- MH_amount - from mhq output
- LD_chr - **Additional**: chromosome of deletion start
- LD_start - **Additional**: bed coordinate start of deletion span
- LD_stop - **Additional**: bed coordinate end of deletion span
- strand - **Additional**: strand that deletion was mapped to (or strand that microhomology is reported for)
- sg_start - **Additional**: bed coordinate start of sgRNAs span (if pairs of sgRNAs used, then start coordinate of 5' most sgRNA)
- sg_stop - **Additional**: bed coordinate end of sgRNAs span (if pairs of sgRNAs used, then end coordinate of 3' most sgRNA)

To gather the bed format information for the additional columns, a tool like UCSC BLAT can be used.

Output is a dataframe containing the original data and two additional columns:
- altMH_count - Number of alternative microhomologies found
- add_delSize - Additional size of deletion beyond sgRNAs (same as deletion size when using one sgRNA)
* * *
## To install mhscanR via RStudio
```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("d0minicO/mhscanR")
library(mhscanR)
```

Requires Biostrings, stringi, stringr, BSgenome, GenomicRanges and tidyverse. If the mhscanR installation fails,see http://bioconductor.org/ and https://www.tidyverse.org/ for help installing these packages first.

Queries, bugs, or discussions welcome: dominic.owens@balliol.ox.ac.uk

License
----

MIT
