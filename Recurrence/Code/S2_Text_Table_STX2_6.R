library(Biostrings)

#' Get the counts for homopolymers of any length in the genome for A, C, G or T
#' @param nucleotide: A, C, G or T
#' @param genome_dnastringset: genome sequence in dnastringset format(GRCh37/h19)
#' @return  counts for homopolymers of any length in the genome for A, C, G or T
countPolymersInGenome <- function(nucleotide, genome_dnastringset)
{
    base <- DNAString(paste(rep(x=nucleotide, times=1), collapse=""))
    base_matches <-vmatchPattern(base, genome_dnastringset)
    
    base_matches_reduced <- c()
    
    for(i in 1:length(base_matches))
    {
      base_matches_reduced <- c(base_matches_reduced, reduce(base_matches[[i]]))
    }
    
    base_matches_width <- c()
    
    for(i in 1:length(base_matches_reduced))
    {
      print(i)
      base_matches_width <- c(base_matches_width, width(base_matches_reduced[[i]]))
    }
    
    counts_width_base <- table(base_matches_width)
    
    return(counts_width_base)
} 

#' Get the combined counts for homopolymers of any length in the genome for A/T or C/G
#' @param pyr_bp: C or T
#' @param pur_bp: A or G
#' @param file_fastaGenome: location of the file with the genome sequence (GRCh37/h19)
#' @return combined counts for homopolymers of any length in the genome for A/T or C/G
getCountsPolymers_pyr_pur_pair <- function(pyr_bp, pur_bp, file_fastaGenome)
{
  genome_dnastringset <- readDNAStringSet(file_fastaGenome)
  genome_dnastringset <- genome_dnastringset[names(genome_dnastringset) %in% c(1:22, "X", "Y")]
  
  counts_polymers_pyr <- countPolymersInGenome(pyr_bp, genome_dnastringset)
  counts_polymers_pur <- countPolymersInGenome(pur_bp, genome_dnastringset)
  
  same_length <- names(counts_polymers_pur)[which(names(counts_polymers_pur) %in% names(counts_polymers_pyr))]
  counts_polymers_pyr_pur_pair <- counts_polymers_pur[same_length] + counts_polymers_pyr[same_length]
  counts_polymers_pyr_pur_pair <- c(counts_polymers_pyr_pur_pair, counts_polymers_pur[which(!(names(counts_polymers_pur) %in% names(counts_polymers_pyr)))], counts_polymers_pyr[which(!(names(counts_polymers_pyr) %in% names(counts_polymers_pur)))])
  
  return(counts_polymers_pyr_pur_pair)
}