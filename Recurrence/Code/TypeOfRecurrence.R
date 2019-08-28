library(parallel)

source("main_utils.R")

#' Add to each recurrent mutation in how many samples of each tumour type it is present 
#'
#' @param rec_mutations: data.frame with the recurrent mutations
#' @param mutation_type: ssm or sim
#' @param sample2ttype: data.frame linking the sample id to the tumor type
#' @param sample_info_file: file with the mapping between the sample ID, the original VCF file with SSMs, the original VCF file with SIMs, the VCF file with recurrent SSMs and the VCF file with the recurrent SIMs 
#'
#' @return data.frame with the recurrent mutations linked to the number of samples they are found in for each tumour type
addTtype2RecMut <- function(rec_mutations, mutation_type, sample2ttype, sample_info_file)
{
  sample_info <- read.table(file = sample_info_file, header=TRUE, sep = "\t", stringsAsFactors = FALSE)
  tumor_types <- sort(unique(as.character(sample2ttype$tumor_type)))
  
  ttype_info <- as.data.frame(matrix(data=0, nrow=nrow(rec_mutations), ncol=length(tumor_types)))
  colnames(ttype_info) <- tumor_types
  rec_mutations_2_ttype <- cbind(rec_mutations, ttype_info)
  
  for(i in 1:nrow(sample2ttype))
  {
    print(i)
    recurrent_mut_cur_sample <- vcf2df(sample_info[which(sample_info$sample_id == sample2ttype[i,"sample_id"]),paste("loc_recurrent_", mutation_type, "_vcf",sep="")],TRUE)
    recurrent_mut_cur_sample$mut_barcode <-  paste(recurrent_mut_cur_sample$CHROM, recurrent_mut_cur_sample$POS, recurrent_mut_cur_sample$REF, recurrent_mut_cur_sample$ALT, sep="_")
    rec_mutations_2_ttype[which(rec_mutations_2_ttype$mut_barcode %in% recurrent_mut_cur_sample$mut_barcode), sample2ttype[i,"tumor_type"]] <- rec_mutations_2_ttype[which(rec_mutations_2_ttype$mut_barcode %in% recurrent_mut_cur_sample$mut_barcode), sample2ttype[i,"tumor_type"]] + 1
  }
  
  return(rec_mutations_2_ttype)
}


#' Determine for each recurrent mutation whether it is recurrent only when considering the entire cohort ('pan-cancer'), recurrent within one tumour type or
#' recurrent in multiple tumour types.
#'
#' @param rec_mutations_2_ttype: data.frame with the recurrent mutations linked to the number of samples they are found in for each tumour type
#' @param ttype_cols: columns with counts for each tumour type in the data.frame 'rec_mutations_2_ttype' 
#' @param num_cores: number of cores to use in the mclapply
#'
#' @return data.frame with the recurrent mutations linked to the number of samples they are found in for each tumour type and the type of recurrence
getTypeOfRecurrence <- function(rec_mutations_2_ttype, ttype_cols, num_cores)
{
  
  rec_mut_2_typeOfRecurrence <- do.call(rbind, mclapply(1:nrow(rec_mutations_2_ttype), function(i)
  {
    print(i)
    num_rec_mut_per_ttype <- rec_mutations_2_ttype[i,ttype_cols]
    
    numSamples2NumTtypes <- as.data.frame(table(as.matrix(num_rec_mut_per_ttype)))
    colnames(numSamples2NumTtypes) <- c("num_samples", "num_ttypes")
    
    numSamples2NumTtypes <- numSamples2NumTtypes[which(numSamples2NumTtypes$num_samples != 0),]
    numSamples2NumTtypes$num_samples <- as.numeric(as.character(numSamples2NumTtypes$num_samples))
    
    if(max(numSamples2NumTtypes$num_samples) == 1)
    {
      rec_mutations_2_ttype[i,"rec_type"] <- "pancancer_only"
      
    } else if(!(1 %in% numSamples2NumTtypes$num_samples)){
      
      if(nrow(numSamples2NumTtypes) > 1){
        rec_mutations_2_ttype[i,"rec_type"] <- "multiple_cancers_unique"
        
      } else if(numSamples2NumTtypes$num_ttypes > 1){
        rec_mutations_2_ttype[i,"rec_type"] <- "multiple_cancers_unique"
        
      } else {
        rec_mutations_2_ttype[i,"rec_type"] <- "single_cancer_unique"
      }

    } else{
      
      if(nrow(numSamples2NumTtypes[which(numSamples2NumTtypes$num_samples > 1),]) > 1){
        rec_mutations_2_ttype[i,"rec_type"] <- "multiple_cancers_notUnique"
        
      } else if(numSamples2NumTtypes[which(numSamples2NumTtypes$num_samples > 1),"num_ttypes"] > 1){
        rec_mutations_2_ttype[i,"rec_type"] <- "multiple_cancers_notUnique"
        
      } else{
        rec_mutations_2_ttype[i,"rec_type"] <- "single_cancer_notUnique"
      }
    }
    
    return(rec_mutations_2_ttype[i,])
  }, mc.cores=num_cores))
  
  return(rec_mut_2_typeOfRecurrence)
}