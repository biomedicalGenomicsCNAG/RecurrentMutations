library(Biostrings)

#' Pair the reverse complements to each other
#' @param all_motifs: list of motifs
#' @return data.frame with all motifs and the reverse complements
getMotifPairs <- function(all_motifs)
{
  motif_pairs <- as.data.frame(matrix(nrow=length(all_motifs), ncol=3))
  colnames(motif_pairs) <- c("org_motif", "revCompl_motif", "paired_motif")
  
  motif_pairs[,"org_motif"] <- all_motifs
  motif_pairs[,"revCompl_motif"] <- as.character(reverseComplement(DNAStringSet(all_motifs)))
  motif_pairs[,"paired_motif"] <- apply(motif_pairs, 1, function(x) {
    paste(sort(c(x["org_motif"],x["revCompl_motif"])), collapse="_")
  })
  
  motif_pairs <- motif_pairs[!duplicated(motif_pairs$pair),]
  
  return(motif_pairs)
}

#'Estimate the frequence of a k-mer in the genome
#' @param length_kmer: length of k-mer
#' @param file_fastaGenome: location of the file with the genome sequence (GRCh37/h19)
#' @return data.frame with the estimation of the number of kmers in the genome
getEstimationFreqKmer <- function(length_kmer, fastafileGenome)
{
  genome_dnastringset <- readDNAStringSet(fastafileGenome)
  genome_dnastringset <- genome_dnastringset[names(genome_dnastringset) %in% c(1:22, "X", "Y")]
  
   #get counts of all i bp motifs per chromosome
  all_motifs_counts_perChr <- oligonucleotideFrequency(genome_dnastringset, width=k, step=k)

  #get counts of all i bp motifs in total
  all_motifs_counts_total <- colSums(all_motifs_counts_perChr)

  #pair the reverse complements to each other
  all_motifs <- names(all_motifs_counts_total)
  motif_pairs <- getMotifPairs(all_motifs)

  #sum the reverse complements
  estimFreqKmer <- as.data.frame(matrix(nrow=1, ncol=length(all_motifs_counts_total)))
  colnames(estimFreqKmer) <- names(all_motifs_counts_total)
  
  for(t in 1:length(all_motifs_counts_total))
  {
    print(t)
    curMotif <- names(all_motifs_counts_total)[t]
    curRevCompl <- as.character(reverseComplement(DNAStringSet(curMotif)))
    estimFreqKmer[1,curMotif] <- all_motifs_counts_total[curMotif] + all_motifs_counts_total[curRevCompl] 
  }
  
  return(estimFreqKmer)
}

#' Get the motifs with a specific base in position '?ndex_mut'.
#' @param index_mut: index of the mutation in the motif
#' @param base: base that is mutated
#' @param all_Motifs: all possible motifs 
#' @return list of specific motifs with 'base' in position 'index_mut'.
getSpecificMotifs <- function(index_mut, base, all_motifs)
{
  specific_motifs <- c()
  
  for(i in 1:length(all_motifs))
  {
    cur_motif_split <- unlist(strsplit(x = all_motifs[i], split = ""))
    
    if(cur_motif_split[index_mut] == base)
      specific_motifs <- c(specific_motifs, all_motifs[i])
  }
  
  return(specific_motifs)
}

#' Get the statistics for the enriched sequence motifs found in context of the genome and the cluster
#'
#' @param enriched_motifs: one or more enriched motifs found of equal length and ungapped
#' @param ref_base_mutation: the reference base that is mutated
#' @param loc_mutation_in_motif: the location of the base in the motif 
#' @param num_bp_context_available: number of bp retrieved at both sides of the mutation
#' @param mutations_recurrent: data.frame with the recurrent mutations in the cluster including sequence context 
#' @param mutations_nonRecurrent: data.frame with the non-recurrent mutations in the cluster including sequence context  
#' @param file_fastaGenome: location of the file with the genome sequence (GRCh37/h19)
#'
#' @example getEnrichmentStatistics("AACTT", "T", 4, mut_rec_clustL_TG, mut_nonRec_clustL_TG, file_fastaGenome)
#'
#' @return data.frame with the statistics
getEnrichmentStatistics <- function(enriched_motifs, ref_base_mutation, loc_mutation_in_motif, num_bp_context_available, mutations_recurrent, mutations_nonRecurrent, file_fastaGenome)
{
  statsEnrichedMotif <- as.data.frame(matrix(data=0, nrow=1, ncol=10))
  colnames(statsEnrichedMotif) <- c("perc_motif_genome_allKmer", "perc_motif_genome_KmerMatchLocSSM", "num_AllSSMs_withEnrichedMotif", "num_recSSMs_withEnrichedMotif", "num_nonRecSSMs_withEnrichedMotif", "num_AllSSMs_withoutEnrichedMotif", "num_recSSMs_withoutEnrichedMotif", "num_nonRecSSMs_withoutEnrichedMotif", "perc_AllSSMs_withEnrichedMotif", "perc_recSSMs_withEnrichedMotif", "perc_nonRecSSMs_withEnrichedMotif")
  
  length_motif <- unique(length(enriched_motifs))
  
  freq_allKmerInGenome <- getEstimationFreqKmer(length_motif, file_fastaGenome)
  
  all_kmers_MatchLocSSM <- getSpecificMotifs(loc_mutation_in_motif, ref_base, names(freq_allKmerInGenome))
  
  total_allKmer_genome <- sum(freq_allKmerInGenome)
  total_counts_KmerMatchLocSSM <- sum(total_allKmer_genome[which(names(total_allKmer_genome) %in% all_kmers_MatchLocSSM)])
  counts_enrichedMotif_genome <- total_allKmer_genome[which(names(total_allKmer_genome) %in% enriched_motifs)]
  
  # statistics in genome
  statsEnrichedMotif[1, "perc_motif_genome_KmerMatchLocSSM"] <- (counts_enrichedMotif_genome*100)/total_counts_KmerMatchLocSSM
  statsEnrichedMotif[1, "perc_motif_genome_allKmer"] <- (counts_enrichedMotif_genome*100)/total_counts_allKmer
  
  end_motif <- num_bp_context_available + 1 + (length_motif-loc_mutation_in_motif)
  start_motif <- end_motif - length_motif - 1
  
  mutations_recurrent$seq_context <- substr(x=mutations_recurrent[, paste("seq_context_", num_bp_context_available, "_bp_pyr",sep="")], start=start_motif, stop=end_motif)
  mutations_non_recurrent$seq_context <- substr(x=mutations_non_recurrent[, paste("seq_context_", num_bp_context_available, "_bp_pyr",sep="")], start=start_motif, stop=end_motif)
  
  counts_Kmers_allSSMs <- table(c(mutations_recurrent[, "seq_context"],mutations_non_recurrent[, "seq_context"]))
  counts_Kmers_recSSMs <- table(mutations_recurrent[, "seq_context"])
  counts_Kmers_nonRecSSMs <- table(mutations_non_recurrent[, "seq_context"])
  
  num_allSSMs_withEnrichedMotif <- sum(counts_Kmers_allSSMs[which(names(counts_Kmers_allSSMs) %in% enriched_motifs)])
  num_recSSMs_withEnrichedMotif <- sum(counts_Kmers_recSSMs[which(names(counts_Kmers_recSSMs) %in% enriched_motifs)])
  num_nonRecSSMs_withEnrichedMotif <- sum(counts_Kmers_nonRecSSMs[which(names(counts_Kmers_nonRecSSMs) %in% enriched_motifs)])
  
  statsEnrichedMotif[1, "num_allSSMs_withEnrichedMotif"] <- num_allSSMs_withEnrichedMotif
  statsEnrichedMotif[1, "num_recSSMs_withEnrichedMotif"] <- num_recSSMs_withEnrichedMotif
  statsEnrichedMotif[1, "num_nonRecSSMs_withEnrichedMotif"] <- num_nonRecSSMs_withEnrichedMotif
  
  statsEnrichedMotif[1, "num_allSSMs_withoutEnrichedMotif"] <- sum(counts_Kmers_allSSMs[which(!(names(counts_Kmers_allSSMs) %in% enriched_motifs))])
  statsEnrichedMotif[1, "num_recSSMs_withoutEnrichedMotif"] <- sum(counts_Kmers_recSSMs[which(!(names(counts_Kmers_recSSMs) %in% enriched_motifs))])
  statsEnrichedMotif[1, "num_nonRecSSMs_withoutEnrichedMotif"] <- sum(counts_Kmers_nonRecSSMs[which(names(!(counts_Kmers_nonRecSSMs) %in% enriched_motifs))])
  
  statsEnrichedMotif[1, "perc_allSSMs_withEnrichedMotif"] <- round((num_AllSSMs_withEnrichedMotif*100)/sum(counts_Kmers_allSSMs),1)
  statsEnrichedMotif[1, "perc_recSSMs_withEnrichedMotif"] <- round((num_recSSMs_withEnrichedMotif*100)/sum(counts_Kmers_recSSMs),1)
  statsEnrichedMotif[1, "perc_nonRecSSMs_withEnrichedMotif"] <- round((num_nonRecSSMs_withEnrichedMotif*100)/sum(counts_Kmers_nonRecSSMs),1)
  
  return(statsEnrichedMotif)
}