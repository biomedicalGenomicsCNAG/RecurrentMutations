library(parallel)
library(vcfR)
library(Rsamtools)
library(GenomicFeatures)
library(rtracklayer)
library(plyr)


#' convert VCF to data.frame, filter if needed
#' @param vcf_filename: name of VCF file to convert to data.frame
#' @param isFiltered: boolean to indicate whether the VCF file still needs to be filtered or not based on the FILTER column
#' @return vcf_df: VCF file converted to data.frame
vcf2df <- function(vcf_filename, isFiltered)
{
  # read in VCF file
  vcf_file <- read.vcfR(vcf_filename, verbose=FALSE)
  
  # convert to data frame
  if(nrow(vcf_file@fix) != 1)
    vcf_df <- cbind(as.data.frame(getFIX(vcf_file), stringsAsFactors =FALSE), INFO2df(vcf_file))
  else if(nrow(vcf_file@fix) == 1)
    vcf_df <- cbind(t(as.data.frame(getFIX(vcf_file),stringsAsFactors =FALSE)), INFO2df(vcf_file))
  
  if(nrow(vcf_file@gt) > 0)
    vcf_df <- cbind(vcf_df,as.data.frame(vcf_file@gt,stringsAsFactors =FALSE))
  
  # only keep calls that passed all filters  
  if(!isFiltered) {
    vcf_df <- vcf_df[which(is.na(vcf_df$FILTER) | vcf_df$FILTER == "PASS"),]
  }
  
  
  vcf_df$REF <- as.character(vcf_df$REF)
  vcf_df$ALT <- as.character(vcf_df$ALT)
  
  return(vcf_df)
}

#' Convert the mutations in the data.frame to the GRanges format.
#' @param mutations_df: mutations in data.frame format
#' @param mutation_type: SSM or SIM
#' @return mutations_granges: mutations in GRange format
convert2GRanges <- function(mutations_df, mutation_type)
{   
  mut_locs <- unlist(IRangesList(apply(mutations_df, 1, function(x)
  {
    
    if(mutation_type == "sim")                        
    {
      if(substr(start=1, stop = 1, x= x["REF"]) !=  substr(start=1, stop = 1, x= x["ALT"]))
      {
        start_pos <- as.numeric(as.character(x["POS"]))
        end_pos <- start_pos + nchar(as.character(x["REF"]))-1
      }
      else if(nchar(x["REF"]) >  1)
      {
        start_pos <- as.numeric(as.character(x["POS"])) + 1
        end_pos <- start_pos + nchar(as.character(x["REF"]))-2
      }
      else if(nchar(x["ALT"]) >  1)
      {
        start_pos <- as.numeric(as.character(x["POS"]))
        end_pos <- start_pos
      }
      
      mut_loc <- IRanges(start=start_pos, end=end_pos)
    }
    else  
    {
      mut_loc <- IRanges(start=as.numeric(as.character(x["POS"])), end=as.numeric(as.character(x["POS"])))
    }
    
    return(mut_loc)
  })))
  
  mutations_df$CHROM <- gsub("MT", "M", mutations_df$CHROM)
  
  mutations_granges <- GRanges(seqnames=paste("chr", mutations_df$CHROM, sep=""), ranges=mut_locs, strand="*")
  names(mutations_granges) <- NULL
  
  return(mutations_granges)
}