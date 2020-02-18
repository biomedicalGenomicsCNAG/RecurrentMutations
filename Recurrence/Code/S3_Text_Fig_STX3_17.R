library(parallel)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(GenomicFeatures)
library(Rsamtools)
library(rtracklayer)

source("main_utils.R")
source("pcawg.colour.palette.R")

#' Tile the genome with overlapping windows, set the window size such that if all mutations would be equally spread there are 'avg_num_mut_in_tile' mutations per 'tile/window'
#' @param num_rec_mutations: total number of recurrent mutations used to determine window size
#' @param avg_num_mut_in_tile: average number of mutation in tile if mutations would be equally spread across the genome 
#' @param file_fastaGenome: location of the file with the genome sequence (GRCh37/h19)
#' @param num_cores: number of cores to use in mclapply
#' @return tiled_genome: genome divided into 'tiles/windows'
getTiledGenome <- function(rec_ssms_file, rec_sims_file, avg_num_mut_in_tile, file_fastaGenome, num_cores)
{
  rec_ssms <- read.table(file = rec_ssms_file, header=TRUE, comment.char="")
  colnames(rec_ssms)[1] <- "CHROM"
  rec_ssms$num_samples <- as.numeric(sub("rec_count=", "", rec_ssms$Count))
  
  rec_sims <- read.table(file = rec_sims_file, header=TRUE, comment.char="")
  colnames(rec_sims)[1] <- "CHROM"
  rec_sims$num_samples <- as.numeric(sub("rec_count=", "", rec_sims$Count))
  
  num_rec_mutations <- sum(rec_ssms$num_samples) + sum(rec_sims$num_samples)
  
  genome_fasta <- open(FaFile(file_fastaGenome)) 
  
  # get length of chromosomes
  length_chr <- seqlengths(genome_fasta)
  length_chr <- length_chr[c(as.character(1:22),"X","Y")]
  
  # get length of genome
  size_genome <- sum(as.numeric(length_chr))
  
  # determine tile size
  tile_size <- ceiling((avg_num_mut_in_tile*size_genome) / num_rec_mutations)
  
  # ensure the tile size is an even number
  if((tile_size %% 2) != 0)
    tile_size <- tile_size + 1
  
  # amount of overlap between tiles
  size_overlap_tiles <- tile_size/2
  
  #tile genome
  tiles <- tileGenome(seqlengths = length_chr, tilewidth = tile_size, cut.last.tile.in.chrom = TRUE)
  
  # construct the overlapping tiles
  overlapping_tiles <- tiles
  overlapping_tiles <- trim(shift(overlapping_tiles, size_overlap_tiles))
  
  # correct tiles that due to the 'shift' exceed the 5' and 3' ends of the chromosome
  overlapping_tiles_corrected <- do.call(c, mclapply(names(length_chr), function(x)
  {
    cur_tiles <- overlapping_tiles[which(seqnames(overlapping_tiles) %in% x),]
    
    if(any(start(cur_tiles) >= length_chr[x] | end(cur_tiles) > length_chr[x]))
      cur_tiles <- cur_tiles[-which(start(cur_tiles) >= length_chr[x] | end(cur_tiles) > length_chr[x]),]
    
    first_partial_tile <- cur_tiles[1,]
    start(first_partial_tile) <- 1
    end(first_partial_tile) <- size_overlap_tiles
    
    cur_tiles <- c(cur_tiles, first_partial_tile)
    
    return(cur_tiles)
  }, mc.cores=num_cores))  
  
  # sort as the two tile sets combined are not ordered according to chromosome and position
  tiled_genome <- sort(c(tiles, overlapping_tiles_corrected))
  
  tiled_genome$tile_index <- 1:length(tiled_genome)
  
  seqlevelsStyle(tiled_genome)  <- "UCSC"
  
  return(tiled_genome)
}

#' Collect the recurrent mutations of samples
#'
#' @param sample_ids: list of samples of interest
#' @param mutation_type: ssm or sim
#' @param recVcfIsFiltered: is the VCF file with the recurrent mutations already filtered or not? 
#' @param sample_info_file: path+filename of the file with the location of the VCF files 
#'
#' @return recurrent mutations in GRange format
getRecMutationsOfSamples <- function(sample_ids, mutation_type, recVcfIsFiltered, sample_info_file)
{
  sample_info <- read.table(file = sample_info_file, header=TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  rec_mut_granges <- do.call(c,unlist(lapply(sample_ids, function(i){
     
     recurrent_mut_cur_sample <- vcf2df(sample_info[which(sample_info$sample_id == i),paste("loc_recurrent_", mutation_type, "_vcf",sep="")],recVcfIsFiltered)
     
     if(nrow(recurrent_mut_cur_sample) > 0){
       
       cur_rec_mut_granges <- convert2GRanges(recurrent_mut_cur_sample, mutation_type, FALSE)
       return(cur_rec_mut_granges)
     }
   })))                           
  
  return(rec_mut_granges)
}


#' Find the tiles in the genome that contain a recurrent mutation
#'
#' @param rec_mutations_granges: recurrent mutations in GRange format  
#' @param tiled_genome: tiled genome in GRange format 
#'
#' @return data.frame with the tiles that contain recurrent mutations and how many
matchMutations2TiledGenome <- function(rec_mutations_granges, tiled_genome)
{
  hits_mutationsVStiles <- findOverlaps(rec_mutations_granges, subject = tiled_genome, type="any")
  
  mut_hits_df <- as.data.frame(rec_mutations_granges[queryHits(hits_mutationsVStiles), ])
  mut_hits_df$mut_index <- 1:nrow(mut_hits_df)
  
  tiles_withHits_df <- as.data.frame(tiled_genome[subjectHits(hits_mutationsVStiles), ])
 
  tile2Mut <- cbind(tiles_withHits_df, mut_index=mut_hits_df[,"mut_index"])
  
  tile_counts <- table(tile2Mut$tile_index)
  
  tile2numMut <- unique(tile2Mut[, c("tile_index", "seqnames")])
  rownames(tile2numMut) <- as.character(tile2numMut$tile_index)
  
  tile2numMut[names(tile_counts), "num_mutations"] <- tile_counts
  
  return(tile2numMut)
}  


#' If mutiple genes match the same tile collapse them into a single string
#'
#' @param tile2genAnn: tiles linked to genes (a tile may be linked to multiple genes)
#' @param num_cores: number of cores to be used in mclapply 
#'
#' @return tiles linked to genes (a tile is linked to a single string of genes)
collapseTileAnnotation <- function(tile2genAnn, start, end, num_cores)
{
   tile_indices <- unique(tile2genAnn$tile_index)
    
   tile2genAnn_collapsed <- do.call(rbind, mclapply(tile_indices[start:end], function(x){
     
     cur_tiles <- tile2genAnn[which(tile2genAnn$tile_index == x),]
     cur_genes <- paste(unique(cur_tiles$gene_name), collapse= " | ")
     
     cur_tile2genAnn_collapsed <- cur_tiles[1,]
     cur_tile2genAnn_collapsed$genes <- cur_genes
     
     return(cur_tile2genAnn_collapsed)
     
   },mc.cores=num_cores))
   
   return(tile2genAnn_collapsed)
}

#' Match the tiles that contain mutations to GENCODE annotation
#'
#' @param tiled_genome: tiled genome in GRange format  
#' @param annotationDir: directory of the GENCODE annotation  
#' @param num_cores: number of cores to be used in mclapply 
#' 
#' @return tiles linked to the overlapping genes 
matchGenAnn2TiledGenome <- function(tiled_genome, annotationDir, num_cores)
{
  annotation_v19 <- import(paste(annotationDir, "gencode.v19.annotation.gtf.gz", sep=""))
  
  hits_tilesVsgeneAnnotation <- findOverlaps(tiled_genome, subject = annotation_v19, type="any")
  
  tile_hits <- as.data.frame(tiled_genome[queryHits(hits_tilesVsgeneAnnotation), ])
  
  genAnn_hits <- as.data.frame(annotation_v19[subjectHits(hits_tilesVsgeneAnnotation), ])
  
  tile2genAnn <- unique(cbind(tile_hits, gene_name=genAnn_hits[,"gene_name"]))
  
  tile2genAnn_collapsed <- collapseTileAnnotation(tile2genAnn, num_cores)
  
  rm(annotation_v19)
  gc()
  
  return(tile2genAnn_collapsed)
}  

#' Find the tiles in the genome with recurrent mutations of the listed samples and annotate them with GENCODE.
#'
#' @param sample_ids: list of samples of interest 
#' @param recVcfIsFiltered: is the VCF file with the recurrent mutations already filtered or not?  
#' @param sample_info_file: path+filename of the file with the location of the VCF files  
#' @param tiled_genome: tiled genome in GRange format  
#' @param annotationDir: directory of the GENCODE annotation   
#' @param num_cores: number of cores to be used in mclapply  
#'
#' @return tiles in the genome with recurrent mutations of the listed samples and annotated with GENCODE.
getAnnotatedTiles2numMutations <- function(sample_ids,recVcfIsFiltered, sample_info_file, tiled_genome, annotationDir, num_cores)
{
  rec_ssms_granges <- getRecMutationsOfSamples(sample_ids, "ssm", recVcfIsFiltered, sample_info_file)
  rec_sims_granges <- getRecMutationsOfSamples(sample_ids, "sim", recVcfIsFiltered, sample_info_file)
  rec_mutations_granges <- c(rec_ssms_granges, rec_sims_granges)
  
  tile2numMut <- matchMutations2TiledGenome(rec_mutations_granges, tiled_genome)
  tilesWithHit2gene <- matchGenAnn2TiledGenome(tiled_genome[which(tiled_genome$tile_index %in% tile2numMut$tile_index),], annotationDir, num_cores)
  
  annotatedTiles2numMutations <- merge(tile2numMut[, c("tile_index", "seqnames", "num_mutations")], tilesWithHit2gene[,c("tile_index", "genes")], by="tile_index", all.x=TRUE)
  
  return(annotatedTiles2numMutations)
}

#' Plot the tiles against the number of recurrent mutations, colored according to chromosome and annotated with genes if applicable.
#'
#' @param annotatedTiles2numMutations:tiles in the genome with recurrent mutations of the listed samples and annotated with GENCODE. 
#' @param th_numMut_geneAnn: minimum number of mutations in a tile required to show the gene annotation 
#' @param filename_plot: filename for the plot 
#' @param resultsDir: directory to store the plot in 
plotTile2NumMut <- function(annotatedTiles2numMutations, th_numMut_geneAnn, filename_plot, resultsDir)
{
  colors_chr <- pcawg.colour.palette(sub("chr", "", unique(as.character(annotatedTiles2numMutations$seqnames))), "chromosomes")
  names(colors_chr) <-  unique(as.character(annotatedTiles2numMutations$seqnames))
 
  annotatedTiles2numMutations_filt <- annotatedTiles2numMutations[which(annotatedTiles2numMutations$num_mutations >= th_numMut_geneAnn), ]
  annotatedTiles2numMutations_filt <- annotatedTiles2numMutations_filt[which(!(is.na(annotatedTiles2numMutations_filt$genes))),]
  
  scatterplotTile2NumMut <- ggplot(aes(x=tile_index,y=num_mutations,  fill=seqnames),  data=annotatedTiles2numMutations) + scale_y_log10() +scale_fill_manual(values=colors_chr, name="chromosome") +  ylab("number of mutations") +   
    geom_point(size=2, shape=21) +  theme( axis.title.y=element_text(size=12), axis.title.x=element_blank(),axis.text.x =element_blank())  + annotation_logticks(sides = "l") + 
    geom_text(data=annotatedTiles2numMutations_filt, aes(x=tile_index, y=num_mutations, label=genes), size=2, vjust=-1)
  ggsave(file=paste(resultsDir, "/", filename_plot, ".pdf", sep=""))
}