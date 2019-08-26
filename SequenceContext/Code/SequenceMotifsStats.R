library(Biostrings)

#' pair the reverse complements to each other
#' @param all_motifs: list of motifs
#' @return motif_pairs: data.frame with all motifs and the reverse complements
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
#' @return estimFreqKmer: estimation of the number of kmers in the genome
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

#' Get the motifs with a specific base in position 'índex_mut'.
#' @param index_mut: index of the mutation in the motif
#' @param base: base that is mutated
#' @param all_Motifs: all possible motifs 
#' @return specific_motifs: list of specific motifs with 'base' in position 'index_mut'.
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

#######################
# Cluster L: T>G SSMs #
#######################

freq_all5mer <- getEstimationFreqKmer(5, fastafileGenome)
  
# Cluster L: T>G
all_5mers_MatchLocSSM <- getSpecificMotifs(4, 'T', names(freq_all5mer))

clustL_total_counts_all5mer <- sum(freq_all5mer)
clustL_total_counts_5merMatchLocSSM <- sum(freq_all5mer[which(names(freq_all5mer) %in% all_5mers_MatchLocSSM)])
clustL_counts_matchingMotif <- freq_all5mer[which(names(freq_all5mer) == "AACTT")]

perc_clustL_motif_genome_matchingKmer <- (clustL_counts_matchingMotif*100)/clustL_total_counts_5merMatchLocSSM
perc_clustL_motif_genome_allKmer <- (clustL_counts_matchingMotif*100)/clustL_total_counts_all5mer

load("mutations_wSeqInfo_G_clustL_T_TG_S_rec_noDups.RData")
load("mutations_wSeqInfo_G_clustL_T_TG_S_nonRec.RData")

mutations_recurrent_noDups$motif_5bp <- substr(x=mutations_recurrent_noDups$seq_context_10_bp_pyr, start=248, stop=252)
mutations_unique$motif_5bp <- substr(x=mutations_unique$seq_context_10_bp_pyr, start=248, stop=252)

counts_in_clustL_allSSMs <- table(c(mutations_recurrent_noDups[, "motif_5bp"],mutations_unique[, "motif_5bp"]))
counts_in_clustL_recSSMs <- table(mutations_recurrent_noDups[, "motif_5bp"])
counts_in_clustL_nonRecSSMs <- table(mutations_unique[, "motif_5bp"])

perc_clustL_motif_AllSSMs <- round((counts_in_clustL_allSSMs["AACTT"]*100)/sum(counts_in_clustL_allSSMs),1) #13.1
perc_clustL_motif_recSSMs <- round((counts_in_clustL_recSSMs["AACTT"]*100)/sum(counts_in_clustL_recSSMs),1) # 38.9
perc_clustL_motif_nonRecSSMs <- round((counts_in_clustL_nonRecSSMs["AACTT"]*100)/sum(counts_in_clustL_nonRecSSMs),1) #12.1

counts_chiSquared <- matrix(data=0, nrow=2,ncol=2)
colnames(counts_chiSquared) <- c("counts_motif", "counts_notMotif")
rownames(counts_chiSquared) <- c("rec", "not_rec")

counts_chiSquared["rec","counts_motif"] <- sum(counts_in_clustL_recSSMs[which(names(counts_in_clustL_recSSMs) %in% c("AACTT"))])
counts_chiSquared["rec","counts_notMotif"] <- sum(counts_in_clustL_recSSMs[which(!(names(counts_in_clustL_recSSMs) %in% c("AACTT")))])
counts_chiSquared["not_rec","counts_motif"] <-  sum(counts_in_clustL_nonRecSSMs[which(names(counts_in_clustL_nonRecSSMs) %in% c("AACTT"))])
counts_chiSquared["not_rec","counts_notMotif"] <- sum(counts_in_clustL_nonRecSSMs[which(!(names(counts_in_clustL_nonRecSSMs) %in% c("AACTT")))])

chisq.test(counts_chiSquared)


# Cluster E: C>G
load(paste(resultsDir, "/both_strands_counts_motif_4bp.RData",sep=""))

both_strands_counts_motif_4bp <- both_strands_counts_motif

clustE_motifs_genome <- getSpecificMotifs(3, 'C', names(both_strands_counts_motif_4bp))

clustE_total_counts_all4mer <- sum(both_strands_counts_motif_4bp)
clustE_total_counts_4merMatchLocSSM <- sum(both_strands_counts_motif_4bp[which(names(both_strands_counts_motif_4bp) %in% clustE_motifs_genome)])
clustE_counts_matchingMotif <- sum(both_strands_counts_motif_4bp[which(names(both_strands_counts_motif_4bp) %in% c("CTCT", "CTCA"))])

perc_clustE_motif_genome_matchingKmer <- (clustE_counts_matchingMotif*100)/clustE_total_counts_4merMatchLocSSM
perc_clustE_motif_genome_allKmer <- (clustE_counts_matchingMotif*100)/clustE_total_counts_all4mer

load("/project/devel/PCAWG/somatic_consensus_aug2016/data/PCAWG_mutations_per_cluster/mutations_wSeqInfo_G_cluster_E_T_CG_S_rec_noDups.RData")
load("/project/devel/PCAWG/somatic_consensus_aug2016/data/PCAWG_mutations_per_cluster/mutations_wSeqInfo_G_cluster_E_T_CG_S_unique.RData")

mutations_recurrent_noDups$motif_4bp <- substr(x=mutations_recurrent_noDups$seq_context_10_bp_pyr, start=249, stop=252)
mutations_unique$motif_4bp <- substr(x=mutations_unique$seq_context_10_bp_pyr, start=249, stop=252)

counts_in_clustE_allSSMs <- table(c(mutations_recurrent_noDups$motif_4bp,mutations_unique$motif_4bp))
counts_in_clustE_recSSMs <- table(mutations_recurrent_noDups$motif_4bp)
counts_in_clustE_nonRecSSMs <- table(mutations_unique$motif_4bp)

perc_clustE_motif_AllSSMs <- round((sum(counts_in_clustE_allSSMs[which(names(counts_in_clustE_allSSMs) %in% c("CTCA", "CTCT"))])*100)/sum(counts_in_clustE_allSSMs),1) # 32.8
perc_clustE_motif_recSSMs <- round((sum(counts_in_clustE_recSSMs[which(names(counts_in_clustE_recSSMs) %in% c("CTCA", "CTCT"))])*100)/sum(counts_in_clustE_recSSMs),1)#55.0
perc_clustE_motif_nonRecSSMs <- round((sum(counts_in_clustE_nonRecSSMs[which(names(counts_in_clustE_nonRecSSMs) %in% c("CTCA", "CTCT"))])*100)/sum(counts_in_clustE_nonRecSSMs),1) #32.8

counts_chiSquared <- matrix(data=0, nrow=2,ncol=2)
colnames(counts_chiSquared) <- c("counts_motif", "counts_notMotif")
rownames(counts_chiSquared) <- c("rec", "not_rec")

counts_chiSquared["rec","counts_motif"] <- sum(counts_in_clustE_recSSMs[which(names(counts_in_clustE_recSSMs) %in% c("CTCA", "CTCT"))])
counts_chiSquared["rec","counts_notMotif"] <- sum(counts_in_clustE_recSSMs[which(!(names(counts_in_clustE_recSSMs) %in% c("CTCA", "CTCT")))])
counts_chiSquared["not_rec","counts_motif"] <-  sum(counts_in_clustE_nonRecSSMs[which(names(counts_in_clustE_nonRecSSMs) %in% c("CTCA", "CTCT"))])
counts_chiSquared["not_rec","counts_notMotif"] <- sum(counts_in_clustE_nonRecSSMs[which(!(names(counts_in_clustE_nonRecSSMs) %in% c("CTCA", "CTCT")))])

chisq.test(counts_chiSquared)

# Cluster H: C>A
load(paste(resultsDir, "logo/both_strands_counts_motif_6bp.RData",sep=""))

both_strands_counts_motif_6bp <- both_strands_counts_motif

clustH_motifs_genome <- getSpecificMotifs(3, 'C', names(both_strands_counts_motif_6bp))

clustH_total_counts_all6mer <- sum(both_strands_counts_motif_6bp)
clustH_total_counts_6merMatchLocSSM <- sum(both_strands_counts_motif_6bp[which(names(both_strands_counts_motif_6bp) %in% clustH_motifs_genome)])
clustH_counts_matchingMotif <- sum(both_strands_counts_motif_6bp[which(names(both_strands_counts_motif_6bp) %in% c("TTCTTT"))])

perc_clustH_motif_genome_matchingKmer <- (clustH_counts_matchingMotif*100)/clustH_total_counts_6merMatchLocSSM
perc_clustH_motif_genome_allKmer <- (clustH_counts_matchingMotif*100)/clustH_total_counts_all6mer


load("/project/devel/PCAWG/somatic_consensus_aug2016/data/PCAWG_mutations_per_cluster/mutations_wSeqInfo_G_clusterH_T_CA_S_rec_noDups.RData")
load("/project/devel/PCAWG/somatic_consensus_aug2016/data/PCAWG_mutations_per_cluster/mutations_wSeqInfo_G_clusterH_T_CA_S_unique.RData")

load("/project/devel/PCAWG/somatic_consensus_aug2016/data/PCAWG_mutations_clustH/mutations_wSeqInfo_G_clustH_T_CA_S_rec_noDups.RData")
load("/project/devel/PCAWG/somatic_consensus_aug2016/data/PCAWG_mutations_clustH/mutations_wSeqInfo_G_clustH_T_CA_S_nonRec.RData")

mutations_recurrent_noDups$motif_6bp <- substr(x=mutations_recurrent_noDups$seq_context_10_bp_pyr, start=249, stop=254)
mutations_unique$motif_6bp <- substr(x=mutations_unique$seq_context_10_bp_pyr, start=249, stop=254)

counts_in_clustH_allSSMs <- table(c(mutations_recurrent_noDups$motif_6bp,mutations_unique$motif_6bp))
counts_in_clustH_recSSMs <- table(mutations_recurrent_noDups$motif_6bp)
counts_in_clustH_nonRecSSMs <- table(mutations_unique$motif_6bp)

perc_clustH_motif_AllSSMs <- round((sum(counts_in_clustH_allSSMs[which(names(counts_in_clustH_allSSMs) %in% c("TTCTTT"))])*100)/sum(counts_in_clustH_allSSMs),1) #14.3
perc_clustH_motif_recSSMs <- round((sum(counts_in_clustH_recSSMs[which(names(counts_in_clustH_recSSMs) %in% c("TTCTTT"))])*100)/sum(counts_in_clustH_recSSMs),1) #32.2
perc_clustH_motif_nonRecSSMs <- round((sum(counts_in_clustH_nonRecSSMs[which(names(counts_in_clustH_nonRecSSMs) %in% c("TTCTTT"))])*100)/sum(counts_in_clustH_nonRecSSMs),1) #13.7


counts_chiSquared <- matrix(data=0, nrow=2,ncol=2)
colnames(counts_chiSquared) <- c("counts_motif", "counts_notMotif")
rownames(counts_chiSquared) <- c("rec", "not_rec")

counts_chiSquared["rec","counts_motif"] <- sum(counts_in_clustH_recSSMs[which(names(counts_in_clustH_recSSMs) %in% c("TTCTTT"))])
counts_chiSquared["rec","counts_notMotif"] <- sum(counts_in_clustH_recSSMs[which(!(names(counts_in_clustH_recSSMs) %in% c("TTCTTT")))])
counts_chiSquared["not_rec","counts_motif"] <-  sum(counts_in_clustH_nonRecSSMs[which(names(counts_in_clustH_nonRecSSMs) %in% c("TTCTTT"))])
counts_chiSquared["not_rec","counts_notMotif"] <- sum(counts_in_clustH_nonRecSSMs[which(!(names(counts_in_clustH_nonRecSSMs) %in% c("TTCTTT")))])

chisq.test(counts_chiSquared)

# Cluster H: T>G
load(paste(resultsDir, "both_strands_counts_motif_8bp.RData",sep=""))

both_strands_counts_motif_8bp <- both_strands_counts_motif

clustH_motifs_genome <- getSpecificMotifs(5, 'T', names(both_strands_counts_motif_8bp))

clustH_total_counts_all8mer <- sum(both_strands_counts_motif_8bp)
clustH_total_counts_8merMatchLocSSM <- sum(both_strands_counts_motif_8bp[which(names(both_strands_counts_motif_8bp) %in% clustH_motifs_genome)])
clustH_counts_matchingMotif <- sum(both_strands_counts_motif_8bp[which(names(both_strands_counts_motif_8bp) %in% c("AAATTTAT"))])

perc_clustH_motif_genome_matchingKmer <- (clustH_counts_matchingMotif*100)/clustH_total_counts_8merMatchLocSSM
perc_clustH_motif_genome_allKmer <- (clustH_counts_matchingMotif*100)/clustH_total_counts_all8mer



load("/project/devel/PCAWG/somatic_consensus_aug2016/data/PCAWG_mutations_clustH/mutations_wSeqInfo_G_clustH_T_TG_S_rec_noDups.RData")
load("/project/devel/PCAWG/somatic_consensus_aug2016/data/PCAWG_mutations_clustH/mutations_wSeqInfo_G_clustH_T_TG_S_nonRec.RData")


mutations_recurrent_noDups$motif_8bp <- substr(x=mutations_recurrent_noDups$seq_context_10_bp_pyr, start=247, stop=254)
mutations_unique$motif_8bp <- substr(x=mutations_unique$seq_context_10_bp_pyr, start=247, stop=254)

counts_in_clustH_allSSMs <- table(c(mutations_recurrent_noDups$motif_8bp,mutations_unique$motif_8bp))
counts_in_clustH_recSSMs <- table(mutations_recurrent_noDups$motif_8bp)
counts_in_clustH_nonRecSSMs <- table(mutations_unique$motif_8bp)

perc_clustH_motif_AllSSMs <- round((sum(counts_in_clustH_allSSMs[which(names(counts_in_clustH_allSSMs) %in% c("AAATTTAT"))])*100)/sum(counts_in_clustH_allSSMs),1) #1.6
perc_clustH_motif_recSSMs <- round((sum(counts_in_clustH_recSSMs[which(names(counts_in_clustH_recSSMs) %in% c("AAATTTAT"))])*100)/sum(counts_in_clustH_recSSMs),1) #12.5
perc_clustH_motif_nonRecSSMs <- round((sum(counts_in_clustH_nonRecSSMs[which(names(counts_in_clustH_nonRecSSMs) %in% c("AAATTTAT"))])*100)/sum(counts_in_clustH_nonRecSSMs),1) #1.6


counts_chiSquared <- matrix(data=0, nrow=2,ncol=2)
colnames(counts_chiSquared) <- c("counts_motif", "counts_notMotif")
rownames(counts_chiSquared) <- c("rec", "not_rec")

counts_chiSquared["rec","counts_motif"] <- sum(counts_in_clustH_recSSMs[which(names(counts_in_clustH_recSSMs) %in% c("AAATTTAT"))])
counts_chiSquared["rec","counts_notMotif"] <- sum(counts_in_clustH_recSSMs[which(!(names(counts_in_clustH_recSSMs) %in% c("AAATTTAT")))])
counts_chiSquared["not_rec","counts_motif"] <-  sum(counts_in_clustH_nonRecSSMs[which(names(counts_in_clustH_nonRecSSMs) %in% c("AAATTTAT"))])
counts_chiSquared["not_rec","counts_notMotif"] <- sum(counts_in_clustH_nonRecSSMs[which(!(names(counts_in_clustH_nonRecSSMs) %in% c("AAATTTAT")))])

chisq.test(counts_chiSquared)

# Cluster G: C>T
load(paste(resultsDir, "both_strands_counts_motif_6bp.RData",sep=""))
both_strands_counts_motif_6bp <- both_strands_counts_motif
clustG_motifs_genome <- getSpecificMotifs(4, 'C', names(both_strands_counts_motif_6bp))

clustG_total_counts_all6mer <- sum(both_strands_counts_motif_6bp)
clustG_total_counts_6merMatchLocSSM <- sum(both_strands_counts_motif_6bp[which(names(both_strands_counts_motif_6bp) %in% clustG_motifs_genome)])
clustG_counts_matchingMotif <- both_strands_counts_motif_6bp[which(names(both_strands_counts_motif_6bp) == "TTTCCT")]

perc_clustG_motif_genome_matchingKmer <- round((clustG_counts_matchingMotif*100)/clustG_total_counts_6merMatchLocSSM,1) #0.4
perc_clustG_motif_genome_allKmer <- round((clustG_counts_matchingMotif*100)/clustG_total_counts_all6mer,1) #0.1

load("/project/devel/PCAWG/somatic_consensus_aug2016/data/PCAWG_mutations_per_cluster/mutations_wSeqInfo_G_clusterG_T_CT_S_rec_noDups.RData")
load("/project/devel/PCAWG/somatic_consensus_aug2016/data/PCAWG_mutations_per_cluster/mutations_wSeqInfo_G_clusterG_T_CT_S_unique.RData")

mutations_recurrent_noDups$motif_6bp <- substr(x=mutations_recurrent_noDups$seq_context_10_bp_pyr, start=248, stop=253)
mutations_unique$motif_6bp <- substr(x=mutations_unique$seq_context_10_bp_pyr, start=248, stop=253)

counts_in_clustG_allSSMs <- table(c(mutations_recurrent_noDups$motif_6bp,mutations_unique$motif_6bp))
counts_in_clustG_recSSMs <- table(mutations_recurrent_noDups$motif_6bp)
counts_in_clustG_nonRecSSMs <- table(mutations_unique$motif_6bp)

perc_clustG_motif_AllSSMs <- round((counts_in_clustG_allSSMs["TTTCCT"]*100)/sum(counts_in_clustG_allSSMs),1)#5.2
perc_clustG_motif_recSSMs <- round((counts_in_clustG_recSSMs["TTTCCT"]*100)/sum(counts_in_clustG_recSSMs),1)#19.5 
perc_clustG_motif_nonRecSSMs <- round((counts_in_clustG_nonRecSSMs["TTTCCT"]*100)/sum(counts_in_clustG_nonRecSSMs),1) # 4.5 


counts_chiSquared <- matrix(data=0, nrow=2,ncol=2)
colnames(counts_chiSquared) <- c("counts_motif", "counts_notMotif")
rownames(counts_chiSquared) <- c("rec", "not_rec")

counts_chiSquared["rec","counts_motif"] <- sum(counts_in_clustG_recSSMs[which(names(counts_in_clustG_recSSMs) %in% c("TTTCCT"))])
counts_chiSquared["rec","counts_notMotif"] <- sum(counts_in_clustG_recSSMs[which(!(names(counts_in_clustG_recSSMs) %in% c("TTTCCT")))])
counts_chiSquared["not_rec","counts_motif"] <-  sum(counts_in_clustG_nonRecSSMs[which(names(counts_in_clustG_nonRecSSMs) %in% c("TTTCCT"))])
counts_chiSquared["not_rec","counts_notMotif"] <- sum(counts_in_clustG_nonRecSSMs[which(!(names(counts_in_clustG_nonRecSSMs) %in% c("TTTCCT")))])

chisq.test(counts_chiSquared)


# Cluster E: C>G
load(paste(resultsDir, "both_strands_counts_motif_4bp.RData",sep=""))

both_strands_counts_motif_4bp <- both_strands_counts_motif

clustM_motifs_genome <- getSpecificMotifs(3, 'C', names(both_strands_counts_motif_4bp))

clustM_total_counts_all4mer <- sum(both_strands_counts_motif_4bp)
clustM_total_counts_4merMatchLocSSM <- sum(both_strands_counts_motif_4bp[which(names(both_strands_counts_motif_4bp) %in% clustM_motifs_genome)])
clustM_counts_matchingMotif <- sum(both_strands_counts_motif_4bp[which(names(both_strands_counts_motif_4bp) %in% c("AGCT"))])

perc_clustM_motif_genome_matchingKmer <- (clustM_counts_matchingMotif*100)/clustM_total_counts_4merMatchLocSSM
perc_clustM_motif_genome_allKmer <- (clustM_counts_matchingMotif*100)/clustM_total_counts_all4mer

load("mutations_wSeqInfo_G_clustM_T_CG_S_rec_noDups.RData")
load("mutations_wSeqInfo_G_clustM_T_CG_S_nonRec.RData")

mutations_recurrent_noDups$motif_4bp <- substr(x=mutations_recurrent_noDups$seq_context_10_bp_pyr, start=249, stop=252)
mutations_unique$motif_4bp <- substr(x=mutations_unique$seq_context_10_bp_pyr, start=249, stop=252)

counts_in_clustM_allSSMs <- table(c(mutations_recurrent_noDups$motif_4bp,mutations_unique$motif_4bp))
counts_in_clustM_recSSMs <- table(mutations_recurrent_noDups$motif_4bp)
counts_in_clustM_nonRecSSMs <- table(mutations_unique$motif_4bp)

perc_clustM_motif_AllSSMs <- round((sum(counts_in_clustM_allSSMs[which(names(counts_in_clustM_allSSMs) %in% c("AGCT"))])*100)/sum(counts_in_clustM_allSSMs),1) # 3.2%
perc_clustM_motif_recSSMs <- round((sum(counts_in_clustM_recSSMs[which(names(counts_in_clustM_recSSMs) %in% c("AGCT"))])*100)/sum(counts_in_clustM_recSSMs),1)# 25.9%
perc_clustM_motif_nonRecSSMs <- round((sum(counts_in_clustM_nonRecSSMs[which(names(counts_in_clustM_nonRecSSMs) %in% c("AGCT"))])*100)/sum(counts_in_clustM_nonRecSSMs),1) #3%

counts_chiSquared <- matrix(data=0, nrow=2,ncol=2)
colnames(counts_chiSquared) <- c("counts_motif", "counts_notMotif")
rownames(counts_chiSquared) <- c("rec", "not_rec")

counts_chiSquared["rec","counts_motif"] <- sum(counts_in_clustM_recSSMs[which(names(counts_in_clustM_recSSMs) %in% c("AGCT"))])
counts_chiSquared["rec","counts_notMotif"] <- sum(counts_in_clustM_recSSMs[which(!(names(counts_in_clustM_recSSMs) %in% c("AGCT")))])
counts_chiSquared["not_rec","counts_motif"] <-  sum(counts_in_clustM_nonRecSSMs[which(names(counts_in_clustM_nonRecSSMs) %in% c("AGCT"))])
counts_chiSquared["not_rec","counts_notMotif"] <- sum(counts_in_clustM_nonRecSSMs[which(!(names(counts_in_clustM_nonRecSSMs) %in% c("AGCT")))])

chisq.test(counts_chiSquared)