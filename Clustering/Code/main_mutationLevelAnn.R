library(parallel)
library(vcfR)
library(Rsamtools)
library(GenomicFeatures)
library(rtracklayer)
library(plyr)

#' Combine multiple bigWig files with the replication time data into a single file.
#' Data is downloaded from http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeUwRepliSeq - wavelet smoothed signal
#' @param filenames_replication_time: list of filenames with replication time data (one file per cell line) 
#' @param metadataDir: directory where the replication time data is stored
#' @return replication time scores from multiple bigWig files combined
readInBigWigFile <- function(filenames_replication_time, metadataDir)
{
  
  for(i in 1:length(file_names_replication_time))
  {
    print(i)
    cur_bigwig_file <- import(paste(metadataDir,"/", filenames_replication_time[i],sep=""))   
    colnames(elementMetadata(cur_bigwig_file))[1] <- paste("score_org_", sub(".bigWig", "",(sub("WaveSignalRep", "_", sub("wgEncodeUwRepliSeq", "",filenames_bigwig[i])))),sep="")
    
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
#' @param mutation_df: data.frame with mutations
#' @param num_cores: number of cores for mclapply
#' @return data.frame with mutations annotated with GENCODE
collapseGenCodeAnnotation <- function(mutations_df, num_cores)
{
  mut_barcodes <- unique(mutations_df$mut_barcode)
  
  mutations2gencode <-  mclapply(mut_barcodes, function(x){
    
    cur_mutation <- mutations_df[which(mutations_df$mut_barcode == x), ]
    
    annotation_cols <- as.data.frame(matrix(data=FALSE, nrow=1, ncol=9))
    colnames(annotation_cols) <- c("gene", "transcript", "exon", "UTR", "CDS", "stop_codon", "start_codon", "intergenic", "Selenocysteine")
    cur_mutationGenRegionsAnn <- cbind(cur_mutation[1,], annotation_cols)
    
    if(all(is.na(cur_mutation$type)))
      genRegions <- "intergenic"
    else    
      genRegions <- as.character(cur_mutation$type)
    
    cur_mutationGenRegionsAnn[1,genRegions] <- TRUE
    
    return(cur_mutationGenRegionsAnn[1,-grep("type", colnames(cur_mutationGenRegionsAnn))])
    
  }, mc.cores=num_cores)
  
  rm(mut_barcodes)
  rm(mutations_df)
  gc()
  
  mutations2gencode_formatted <- do.call(rbind.fill,mutations2gencode)
  
  rm(mutations2gencode)
  gc()
  
  return(mutations2gencode_formatted)
}

#' Annotate mutations with GENCODE
#' @param mutations_df: data.frame with mutations
#' @param mutations_granges: GRange object with the mutations
#' @param annotation_v19: GENCODE information (https://www.gencodegenes.org/human/releases.html)
#' @param num_cores: number of cores to use in mclapply
#' @return data.frame with mutations annotated with GENCODE
addGencodeAnnotation <- function(mutations_df, mutations_granges, annotation_v19, num_cores)
{

  print("findOverlaps")
  matches_annotation <- findOverlaps(query = mutations_granges, subject = annotation_v19)
  
  mutations_withMatch <- as.data.frame(mutations_df[queryHits(matches_annotation),])
  annotation_matches <- as.data.frame(annotation[subjectHits(matches_annotation),])
  
  rm(annotation)
  gc()
  
  mutations_annotated <- unique(cbind(mutations_withMatch, annotation_matches[,"type"]))
  colnames(mutations_annotated)[ncol(mutations_annotated)] <- "type"
  
  mutation_notMatched <- as.data.frame(mutations_df[-queryHits(matches_annotation),])
  
  if(nrow(mutation_notMatched) > 0)
      mutations_annotated[(nrow(mutations_annotated) + 1): (nrow(mutations_annotated) + nrow(mutation_notMatched)), colnames(mutation_notMatched)] <- mutation_notMatched
  
  rm(mutations_withMatch)
  rm(annotation_matches)
  rm(mutation_notMatched)
  rm(mutations_df)
  rm(mutations_granges)
  rm(matches_annotation)
  gc()
  
  mutations_annotated$mut_barcode <- paste(mutations_annotated$CHROM, mutations_annotated$POS, mutations_annotated$REF, mutations_annotated$ALT,sep="_")
  
  if(nrow(mutations_annotated) > 1)
  {
      chromosomes <- sort(unique(as.character(mutations_annotated$CHROM)))
      
      print("collapseAnnotation")
      # annotate per chromsosome
      for(i in 1:length(chromosomes))
      {
          print(chromosomes[i])
          cur_collapsed <- collapseGenCodeAnnotation(mutations_annotated[which(mutations_annotated$CHROM == chromosomes[i]),], num_cores)
          
          if(i == 1)
            mutations2gencode <-cur_collapsed
          else 
            mutations2gencode <- rbind(mutations2gencode, cur_collapsed)
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
  
  mut2replTimeScores_collapsed <-  mclapply(mut_barcodes, function(x){
    
    cur_mutation <- mutations_replTimeScores[which(mutations_replTimeScores$mut_barcode == x), ]
    
    if(nrow(cur_mutation) > 1)
    {
      cur_mutation[1,grep("score_", colnames(cur_mutation))] <- colMeans(t(apply(cur_mutation[,grep("score_", colnames(cur_mutation))], 1, as.numeric)))
    }
    
    return(cur_mutation[1,])
    
  }, mc.cores=num_cores)
  
  mutations2replTimeScores_formatted <- do.call(rbind.fill,mut2replTimeScores_collapsed)
  
  return(mutations2replTimeScores_formatted)
}

#' Annotate mutations with replication time data
#' @param mutations_df: data.frame with mutations
#' @param mutations_granges: GRange object with the mutations
#' @param mutation_type: SSM or SIM
#' @param replicationTimeScores: replication time data
#' @param num_cores: number of cores to use in mclapply
#' @return data.frame with mutations annotated with replication time scores
addReplTimeScoresAnnotation <- function(mutations_df, mutations_granges, mutation_type, replicationTimeScores, num_cores)
{
    all_matches <- findOverlaps(query = mutations_granges, subject = replicationTimeScores)
    
    mutations_withMatch <- as.data.frame(mutations_df[queryHits(all_matches),])
    replTimeScores <- as.data.frame(replicationTimeScores[subjectHits(all_matches),])
    
    mutation2replTimeScore <- cbind(mutations_withMatch, replTimeScores[,grep("score",colnames(replTimeScores))])
    
    if(mutation_type == "sim")
    {
      mutation2replTimeScore <- collapseReplTimeScores(mutation2replTimeScore, num_cores)
    }
    
    mutation_notMatched <- as.data.frame(mutations_df[-queryHits(all_matches),])
    
    
    if(nrow(mutation_notMatched) > 0)
    {
      mutation_notMatched$mut_barcode <- paste(mutation_notMatched$CHROM, mutation_notMatched$POS, mutation_notMatched$REF, mutation_notMatched$ALT,sep="_")
      mutation2replTimeScore[(nrow(mutation2replTimeScore) + 1): (nrow(mutation2replTimeScore) + nrow(mutation_notMatched)), colnames(mutation_notMatched)] <- mutation_notMatched
      
    }
    
    return(mutation2replTimeScore)
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
    mutations_granges <- convert2GRanges(mutations_df, mutation_type)
    
    sample_id <- sample_info[x,"sample_id"]
    
    # Add functional category from GenCode
    mutations_df_genCode <- addGencodeAnnotation(mutations_df, mutations_granges, annotation_v19, num_cores)
    
    # Add replication time information
    mutations_df_genCode_replTime <- addReplTimeScoresAnnotation(mutations_df_genCode, mutations_granges, mutation_type, replicationTimeScores, num_cores)
    
    rm(mutations_df_genCode)
    gc()
    
    write.table(mutations_df_genCode_replTime, file=paste(annSamplesDir, "/", annSamplesFolder, "/", sample_id, "_annotatedWithGenCode_ReplTime.txt",sep=""), quote = FALSE, sep="\t",row.names = FALSE, col.names = TRUE)
    
    
    rm(mutations_df_genCode_replTime)
    gc()
    
    
    return(1)
  }, mc.cores=num_cores)
}