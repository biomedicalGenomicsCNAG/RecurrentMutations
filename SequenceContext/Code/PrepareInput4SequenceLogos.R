library(rtracklayer)
library(parallel)
library(reshape2)
library(ggplot2)
library(Biostrings)
library(plyr)
library(Rsamtools)

source("main_utils.R")

#' Get sequence context for each mutation
#' @param mutations_df: data.frame with mutations
#' @param num_bp_context: number of bp to retrieve at both sides of the mutation
#' @param file_fastaGenome: location of the file with the genome sequence (GRCh37/h19)
#' @return mutations_df: data.frame with mutations with sequence context of 'num_bp_context'
getSeqContextSSMs <- function(mutations_df,num_bp_context, file_fastaGenome)
{
  genome_fasta <- open(FaFile(file_fastaGenome)) 
  
  chromosomes <- as.character(mutations_df$CHROM)
  chromosomes <- gsub("MT", "M", chromosomes)
  
  coordinates_bases <- GRanges(seqnames = chromosomes, ranges=IRanges(start=as.numeric(as.character(mutations_df$POS))-num_bp_context, end=as.numeric(as.character(mutations_df$POS))+num_bp_context))
  mutations_df$seq_context <- as.character(scanFa(genome_fasta, param=coordinates_bases))
  
  genome_fasta <- close(genome_fasta) 
  
  return(mutations_df)
}

#' Add to the mutations in the given list VCF files the sequence context of 'num_bp_context' in length at both sides
#' Store as tab-delimited file. 
#' @param sample_info_file: file with the mapping between the sample ID, the original VCF file with SSMs, the original VCF file with SIMs, the VCF file with recurrent SSMs and the VCF file with the recurrent SIMs 
#' @param isFiltered: boolean to indicate whether the VCF file still needs to be filtered or not based on the FILTER column
#' @param num_bp_context: number of bp to retrieve at both sides of the mutation
#' @param file_fastaGenome: location of the file with the genome sequence (GRCh37/h19)
#' @param dataDir: directory where the data is stored
#' @param samplesFolder: folder where the VCF files are stored
#' @param resultsDir: diretory to store mutations with the sequence context
addSequenceContexts2SSMs <- function(sample_info_file, filename_column, isFiltered, num_bp_context, file_fastaGenome,dataDir, samplesFolder, num_cores)
{
  tmp <- mclapply(1:rnow(sample_info_file), function(x)
  {

    mutations_df <- vcf2df(paste(dataDir, "/", samplesFolder, "/", sample_info_file[x, filename_column],sep=""), isFiltered)
    
    
    mutations_withContext <- getSeqContextSSMs(mutations_df,num_bp_context, file_fastaGenome)
    
    mutations_withContext$ssm_ref_strand <- paste(mutations_withContext$REF, ">", mutations_withContext$ALT, sep="")
    
    mutations_withContext[which(mutations_withContext$REF == "A"), "REF_nonRefStrand"] <- "T"
    mutations_withContext[which(mutations_withContext$REF == "C"), "REF_nonRefStrand"] <- "G"
    mutations_withContext[which(mutations_withContext$REF == "G"), "REF_nonRefStrand"] <- "C"
    mutations_withContext[which(mutations_withContext$REF == "T"), "REF_nonRefStrand"] <- "A"
    
    mutations_withContext[which(mutations_withContext$ALT == "A"), "ALT_nonRefStrand"] <- "T"
    mutations_withContext[which(mutations_withContext$ALT == "C"), "ALT_nonRefStrand"] <- "G"
    mutations_withContext[which(mutations_withContext$ALT == "G"), "ALT_nonRefStrand"] <- "C"
    mutations_withContext[which(mutations_withContext$ALT == "T"), "ALT_nonRefStrand"] <- "A"
    
    mutations_withContext$ssm_rev_compl <- paste(mutations_withContext$REF_nonRefStrand, ">", mutations_withContext$ALT_nonRefStrand, sep="")
    mutations_withContext$ssm_pyr_standard <- mutations_withContext$ssm_ref_strand
    mutations_withContext[which(mutations_withContext$ssm_ref_strand %in% c("A>C", "A>G", "A>T", "G>A", "G>C", "G>T")),"ssm_pyr_standard"] <- paste(mutations_withContext$REF_nonRefStrand, ">", mutations_withContext$ALT_nonRefStrand, sep="")
    
    mutations_withContext[,paste("seq_context_", num_bp_context, "_bp_refStrand",sep="")] <- mutations_withContext$seq_context
    mutations_withContext[,paste("seq_context_", num_bp_context, "_bp_rev_compl",sep="")] <- as.character(reverseComplement(DNAStringSet(mutations_withContext[,paste("seq_context_", num_bp_context, "_bp_refStrand",sep="")])))
    
    mutations_withContext[which(mutations_withContext$ssm_ref_strand %in% c("A>C", "A>G", "A>T", "G>A", "G>C", "G>T")),paste("seq_context_", num_bp_context, "_bp_pyr_standard",sep="")] <- as.character(reverseComplement(DNAStringSet(mutations_withContext[which(mutations_withContext$ssm_ref_strand %in% c("A>C", "A>G", "A>T", "G>A", "G>C", "G>T")),paste("seq_context_", num_bp_context, "_bp_refStrand",sep="")])))
    mutations_withContext[which(!(mutations_withContext$ssm_ref_strand %in% c("A>C", "A>G", "A>T", "G>A", "G>C", "G>T"))),paste("seq_context_", num_bp_context, "_bp_pyr_standard",sep="")] <- mutations_withContext[which(!(mutations_withContext$ssm_ref_strand %in% c("A>C", "A>G", "A>T", "G>A", "G>C", "G>T"))),paste("seq_context_", num_bp_context, "_bp_refStrand",sep="")]
    
    write.table(mutations_withContext, file=paste(resultsDir, "/ssms_",sample_info_file[x, "sample_id"],"_withSeqContext.txt",sep=""), quote = FALSE,sep="\t", row.names = FALSE, col.names = TRUE)
    
  }, mc.cores=num_cores)
}




#'Combine the indicated type of mutations ('ssm_subtype') of the samples in the list belonging to indicated 'group' in RData format. 
#' @param samplesDir: directory where the samples are stored with the sequence context included for each mutation
#' @param sample_ids: identifiers of the samples of which to combine their mutations
#' @param resultsDir: directory to store the combined set of mutations
#' @param resultsFolder: folder to store the combined set of mutations
#' @param ssm_subtype: CA, CG, CT, TA, TC or TG
#' @param group: name of the group of samples of which the mutations are combined
#' @param num_bp_context: number of bp extracted for sequence context 
#' @example sample_ids_clustH <- sample2clust[which(sample2clust$clust_char == "H"), "sample_id"]
#' @example combineMutations("TG","clustH", 10,samplesDir, sample_ids_clustH,resultsDir,"PCAWG_mutations_clustH")
combineMutationsOfGroup <- function(ssm_subtype, group, num_bp_context, samplesDir, sample_ids, resultsDir, resultsFolder)
{
  
  available_files <- list.files(samplesDir, full.names = TRUE)
  isFirst <- TRUE
  mutations_group <- NULL
  
  for(i in 1:length(sample_ids)) 
  {
    print(i)
    print(sample_ids[i])
    
    if(paste(samplesDir, "/ssms_", sample_ids[i], "_withSeqContext_", ssm_subtype, ".txt", sep="") %in% available_files){
      
      cur_mutations <- read.table(paste(samplesDir, "/ssms_", sample_ids[i], "_withSeqContext_", ssm_subtype, ".txt", sep=""), quote = "",sep="\t", header=TRUE, stringsAsFactors = FALSE, as.is=TRUE, colClasses = "character")
      
      print(nrow(cur_mutations))
      
      if(nrow(cur_mutations) > 0)
      {
        if(isFirst)
        {
          mutations_group <- cur_mutations[,c("CHROM", "POS", "REF","ALT", "mut_barcode", "subtype", paste("seq_context_", num_bp_context, "_bp_pyr",sep=""), paste("seq_context_", num_bp_context, "_bp_pur",sep=""))]
          isFirst <- FALSE
        } else {
          mutations_group <- rbind(mutations_group, cur_mutations[,c("CHROM", "POS", "REF","ALT", "mut_barcode", "subtype", paste("seq_context_", num_bp_context, "_bp_pyr", sep=""), paste("seq_context_", num_bp_context, "_bp_pur", sep=""))])
        }
      }
      
    } else {
      print("file not found")
      print(sample_ids[i])
    }
  }
  
  mutations_recurrent_noDups <- unique(mutations_group[duplicated(mutations_group$mut_barcode),])
  save(mutations_recurrent_noDups, file=paste(resultsDir, resultsFolder, "/mutations_wSeqInfo_G_", group, "_T_", ssm_subtype, "_S_rec_noDups.RData", sep=""))
  
  mutations_nonrecurrent <- mutations_group[which(!(mutations_group$mut_barcode %in% mutations_recurrent_noDups$mut_barcode)),]
  save(mutations_nonrecurrent, file=paste(resultsDir, resultsFolder, "/mutations_wSeqInfo_G_", group, "_T_", ssm_subtype, "_S_nonRec.RData", sep=""))
  
}