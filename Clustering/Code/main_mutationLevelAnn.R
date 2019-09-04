library(parallel)
library(vcfR)
library(Rsamtools)
library(GenomicFeatures)
library(rtracklayer)
library(plyr)

options(scipen=999)

#' Combine multiple bigWig files with the replication time data into a single file.
#' Data is downloaded from http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeUwRepliSeq - wavelet smoothed signal
#' @param filenames_replication_time: list of filenames with replication time data (one file per cell line) 
#' @param metadataDir: directory where the replication time data is stored
#' @return replication time scores of multiple bigWig files combined in an GRange object
readInBigWigFile <- function(filenames_replication_time, metadataDir)
{
  
  for(i in 1:length(filenames_replication_time))
  {
    print(i)
    cur_bigwig_file <- import(paste(metadataDir,"/", filenames_replication_time[i],sep=""))   
    colnames(elementMetadata(cur_bigwig_file))[1] <- paste("score_", sub(".bigWig", "",(sub("WaveSignalRep", "_", sub("wgEncodeUwRepliSeq", "",filenames_bigwig[i])))),sep="")
    
    if(i==1)
    {
      replicationTimeScores <- cur_bigwig_file
    }
    else
    {
      replicationTimeScores <- merge(replicationTimeScores, cur_bigwig_file)
    }
  }
  
  return(replicationTimeScores)
}

#' Collapse the GENCODE annotation into a single entry
#' @param mutations_annotated: GRanges object with the GENCODE annotation, multiple rows per mutation possible
#' @param num_cores: number of cores for mclapply
#' @return data.frame with mutations annotated with GENCODE, single row per mutation.
collapseGenCodeAnnotation <- function(mutations_annotated, num_cores)
{
  mut_barcodes <- unique(mutations_annotated$mut_barcode)
  
  mutations2gencode <- do.call(c, mclapply(mut_barcodes, function(x){
    
    cur_mutation <- mutations_annotated[which(mutations_annotated$mut_barcode == x), ]
    cur_mutationGenRegionsAnn <- cur_mutation[1,]
    
    annotation_cols <- as.data.frame(matrix(data=FALSE, nrow=1, ncol=9))
    colnames(annotation_cols) <- c("gene", "transcript", "exon", "UTR", "CDS", "stop_codon", "start_codon", "intergenic", "Selenocysteine")
    mcols(cur_mutationGenRegionsAnn) <- cbind(mcols(cur_mutation[1,]), annotation_cols)
    
    if(all(is.na(cur_mutation$type)))
      genRegions <- "intergenic"
    else    
      genRegions <- unique(as.character(cur_mutation$type))
    
    mcols(cur_mutationGenRegionsAnn)[,genRegions] <- TRUE
    
    mcols(cur_mutationGenRegionsAnn) <- mcols(cur_mutationGenRegionsAnn)[,-grep("type", colnames(mcols(cur_mutationGenRegionsAnn)))]
    
    return(cur_mutationGenRegionsAnn)
    
  }, mc.cores=num_cores))
  
  rm(mut_barcodes)
  gc()
  
  return(mutations2gencode)
}

#' Annotate mutations with GENCODE
#' @param mutations_granges: GRange object with the mutations
#' @param annotation_v19: GENCODE information (https://www.gencodegenes.org/human/releases.html)
#' @param num_cores: number of cores to use in mclapply
#' @return data.frame with mutations annotated with GENCODE
addGencodeAnnotation <- function(mutations_granges, annotation_v19, num_cores)
{

  print("findOverlaps")
  matches_annotation <- findOverlaps(query = mutations_granges, subject = annotation_v19)
  
  mutations_withMatch <- mutations_granges[queryHits(matches_annotation),]
  annotation_matches <- annotation_v19[subjectHits(matches_annotation),]
  
  rm(annotation_v19)
  gc()
  
  mutations_withMatch$type <- as.character(annotation_matches$type)

  mutation_notMatched <- mutations_granges[-queryHits(matches_annotation),]
  mutation_notMatched$type <- "intergenic"
  
  
  mutations_annotated <- c(mutations_withMatch, mutation_notMatched)

  rm(annotation_matches)
  rm(mutation_notMatched)
  rm(mutations_granges)
  rm(matches_annotation)
  gc()
  
  if(length(mutations_withMatch) > 1)
  {
      rm(mutations_withMatch)
    
      chromosomes <- sort(unique(as.character(mutations_annotated$CHROM)))
      
      print("collapseAnnotation")
      # annotate per chromsosome
      for(i in 1:length(chromosomes))
      {
          cur_collapsed <- collapseGenCodeAnnotation(mutations_annotated[which(mutations_annotated$CHROM == chromosomes[i]),], num_cores)
          
          if(i == 1)
            mutations2gencode <-cur_collapsed
          else 
            mutations2gencode <- c(mutations2gencode, cur_collapsed)
      }
  } else {
    mutations2gencode <- mutations_annotated
  }
  
  rm(mutations_annotated)
  gc()
  
  return(mutations2gencode)
}

#' Replication time is provided for windows of a 1000 bp and some SIMs overlap more than one window.
#' If so, take the average.
#' @param mutations2replTimeScores: data.frame with SIMs annotated with replication time scores
#' @param num_cores: number of cores to use in mclapply
#' @return data.frame with mutations annotated with replication time scores
collapseReplTimeScores <- function(mutations_replTimeScores, num_cores)
{
  mut_barcodes <- unique(mutations_replTimeScores$mut_barcode)
  
  mut2replTimeScores_collapsed <-  do.call(c, mclapply(mut_barcodes, function(x){
    
    cur_mutation <- mutations_replTimeScores[which(mutations_replTimeScores$mut_barcode == x), ]
    
    if(length(cur_mutation) > 1)
    {
        mcols(cur_mutation)[1,grep("score_", colnames(mcols(cur_mutation)))] <- as.data.frame(t(colMeans(as.data.frame(mcols(cur_mutation)[,grep("score_", colnames(mcols(cur_mutation)))]))))
    }
    
    return(cur_mutation[1,])
    
  }, mc.cores=num_cores))
  
  return(mut2replTimeScores_collapsed)
}

#' Annotate mutations with replication time data
#' @param mutations_granges: GRange object with the mutations
#' @param mutation_type: SSM or SIM
#' @param replicationTimeScores: replication time data
#' @param num_cores: number of cores to use in mclapply
#' @return data.frame with mutations annotated with replication time scores
addReplTimeScoresAnnotation <- function(mutations_granges, mutation_type, replicationTimeScores, num_cores)
{
    all_matches <- findOverlaps(query = mutations_granges, subject = replicationTimeScores)
    
    mutations_withMatch <- mutations_granges[queryHits(all_matches),]
    replTimeScores <- replicationTimeScores[subjectHits(all_matches),]
    
    mutations_replTimeScores <- mutations_withMatch
    mcols(mutations_replTimeScores) <- cbind(mcols(mutations_withMatch), mcols(replTimeScores)[,grep("score",colnames(mcols(replTimeScores)))])
    
    if(mutation_type == "sim")
    {
      mutations_replTimeScores <- collapseReplTimeScores(mutations_replTimeScores, num_cores)
    }
    
    mutation_notMatched <- mutations_granges[-queryHits(all_matches),]
    
    if(length(mutation_notMatched) > 0)
    {
        print("missing")
    }
    
    mutations_replTimeScores <- c(mutations_replTimeScores, mutation_notMatched)
    
    return(mutations_replTimeScores)
}

#' Annotate the mutations with the functional category retrieved from GENCODE and the replication time scores retrieved from ENCODE/University of Washington. #
#' @param sample_info: data.frame with the mapping between the sample ID, the original VCF file with SSMs, the original VCF file with SIMs, the VCF file with recurrent SSMs and the VCF file with the recurrent SIMs 
#' @param vcfIsFiltered: boolean to indicate whether or not the VCF file has been filtered based on the FILTER column
#' @param mutation_type: ssm or sim
#' @param annotation_v19: GENCODE information (https://www.gencodegenes.org/human/releases.html)
#' @param replicationTimeScores: replication time data
#' @param annSamplesDir: directory to store the tab-del sample files with the annotated mutations 
#' @param annSamplesFolder: folder to store the tab-del sample files with the annotated mutations 
#' @param num_cores: number of cores to use in mclapply
annotateAtMutationLevel <- function(sample_info, vcfIsFiltered, mutation_type, annotation_v19, replicationTimeScores, annSamplesDir, annSamplesFolder, num_cores)
{
  isFinished <- mclapply(1:nrow(sample_info), function(x)
  {
    print(x)
    
    # read in VCF file with mutations
    if(mutation_type == "ssm")
      mutations_df <- vcf2df(sample_info[x,"loc_all_ssm_vcf"],vcfIsFiltered)
    else 
      mutations_df <- vcf2df(sample_info[x,"loc_all_sim_vcf"],vcfIsFiltered)
    
    mutations_df$mut_barcode <-  paste(mutations_df$CHROM, mutations_df$POS, mutations_df$REF, mutations_df$ALT, sep="_")
    
    #convert to GRange format
    mutations_granges <- convert2GRanges(mutations_df, mutation_type, TRUE)
    
    sample_id <- sample_info[x,"sample_id"]
    
    # Add functional category from GenCode
    mutations_granges_genCode <- addGencodeAnnotation(mutations_granges, annotation_v19, num_cores)
    
    # Add replication time information
    mutations_granges_genCode_replTime <- addReplTimeScoresAnnotation(mutations_granges_genCode, mutation_type, replicationTimeScores, num_cores)
    
    rm(mutations_granges_genCode)
    gc()
    
    write.table(mutations_granges_genCode_replTime, file=paste(annSamplesDir, "/", annSamplesFolder, "/", sample_id, "_annotatedWithGenCode_ReplTime.txt",sep=""), quote = FALSE, sep="\t",row.names = FALSE, col.names = TRUE)
    
    
    rm(mutations_granges_genCode_replTime)
    gc()
    
    
    return(1)
  }, mc.cores=num_cores)
}