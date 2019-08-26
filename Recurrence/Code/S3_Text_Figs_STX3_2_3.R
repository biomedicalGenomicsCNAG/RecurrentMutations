library(ggplot2)
library(Rsamtools)
library(GenomicRanges)


source("pcawg.colour.palette.R") #download from:

#' for the density plot of the replication time scores of the genome, remove the scores for 1000bp with > 'threshold_numN' N's.
#' @param replicationTimeScores: data.frame with replication time scores
#' @param threshold_numN: threshold for the maximum number of N's in the sequence of the 1000 bp window the replication time score is defined for
#' @param file_fastaGenome: location of the file with the genome sequence (GRCh37/h19)
#' @param num_cores: number of cores to use in mclapply
#' @return replicationTimeScores_filt: data.frame with replication time scores with removed windows with too many N's.
filterReplTimScores <- function(replicationTimeScores, threshold_numN, file_fastaGenome, num_cores)
{
    genome_fasta <- open(FaFile(file_fastaGenome)) 
    
    replicationTimeScores_formatted <- GRanges(seqnames = sub("chr", "", seqnames(replicationTimeScores)), ranges=IRanges(start=start(replicationTimeScores), end=end(replicationTimeScores)))
    
    replRegions_bases <- as.character(scanFa(genome_fasta, param=replicationTimeScores_formatted))
    
    replRegions_bases_numN <- as.data.frame(matrix(nrow=length(replRegions_bases), ncol=2))
    replRegions_bases_numN[,"bases"] <- replRegions_bases
    replRegions_bases_numN[,"num_N"] <- unlist(mclapply(1:length(replRegions_bases), function(x)
    {
      print(x)
      num_N <- str_count(replRegions_bases[x], pattern="N")
      return(num_N)
    }, mc.cores=num_cores))
    
    
    replicationTimeScores_filt <- replicationTimeScores[which(replRegions_bases_numN$num_N < threshold_numN),]
    
    return(replicationTimeScores_filt)
}

# each line will represent a sample of a certain tumor type
#' @param sample2ttype: data.frame with samples to plot linked to tumour type 
#' @param replicationTimeScores_df_filt: replication time scores filtered for scores for windows with >= threshold_numN of 'N'.
#' @param filename_plot: name of the file to store the plot in
#' @param dataDir: directory of the data
#' @param samplesFolder: folder of the annotated sample file with replication time scores
#' @param resultsDir: directory to store the plot
plotDensityReplTime <- function(sample2ttype, replicationTimeScores_df_filt, filename_plot,  dataDir, samplesFolder, resultsDir)
{
  firstPlot <- TRUE
  
  # colors of tumor types and overall replication time score for entire genome
  tumor_type_colors <- pcawg.colour.palette(as.character(sample2ttype$tumor_type), "tumour.subtype")
  names(tumor_type_colors) <- sample2ttype$tumor_type
  tumor_type_colors <- unique(c(tumor_type_colors, overall="black"))
  names(tumor_type_colors) <- c(unique(sample2ttype$tumor_type), "overall")
  
  # type of line, solid for samples and dashed for the genome
  tumor_type_line_type <- tumor_type_colors
  tumor_type_line_type[which(names(tumor_type_line_type) != "overall")] <- "solid"
  tumor_type_line_type["overall"] <- "dashed"
  
  # summarize scores of cell lines
  replicationTimeScores_df_filt$score_org_cancer_clines_median <- rowMedians(as.matrix(replicationTimeScores_df_filt[, c("score_org_Helas3_1", "score_org_Mcf7_1", "score_org_Sknsh_1", "score_org_Hepg2_1", "score_org_K562_1")]))
  
  # vertical line indicating median
  median_replTime_median_cancer_clines <- median(replicationTimeScores_df_filt$score_org_cancer_clines_median)

  for(i in 1:nrow(sample2ttype))
  {
    print(i)
    
    tryCatch({
      mutation_info <- read.table(paste(dataDir, "/", samplesFolder, "/", sample2ttype[i,"sample_id"], "_annotatedWithGenCode_ReplTime.txt",sep=""), quote = "",sep="\t", header=TRUE, stringsAsFactors = FALSE, as.is=TRUE, colClasses = "character")
      
     
      num_mutations <- nrow(mutation_info)
      
      mutation_info$tumor_type <- sample2ttype[i,"tumor_type"]
      
      
      mutation_info[, grep("score_org", colnames(mutation_info))] <- apply(mutation_info[, grep("score_org", colnames(mutation_info))], 2, as.numeric)  
      mutation_info$score_org_cancer_clines_median <- rowMedians(as.matrix(mutation_info[, c("score_org_Helas3_1", "score_org_Mcf7_1", "score_org_Sknsh_1", "score_org_Hepg2_1", "score_org_K562_1")]))
      mutation_info <- mutation_info[which(!(is.na(mutation_info$score_org_cancer_clines_median))),]
      
      if(firstPlot)
      {
        
        # only plot sample if it has at least 10 mutations
        if(nrow(mutation_info) >= 10)
        {
          total_plot_median <-  ggplot() + geom_density(data=mutation_info, aes(score_org_cancer_clines_median,  colour = tumor_type), alpha=0.1, kernel="gaussian") + xlab("median replication time score across cell lines") + guides(color=FALSE, linetype=FALSE) + scale_color_manual(values = tumor_type_colors) + scale_linetype_manual(values = tumor_type_line_type) + scale_y_continuous(breaks = seq(0,0.045, by=0.005), limits = c(0,0.045)) + scale_x_continuous(breaks = seq(0,100, by=10)) + theme(axis.title=element_text(size=20), axis.text=element_text(size=18, face="bold"), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
          
          firstPlot <- FALSE
          
        }
        
      } else{
        
        # only plot sample if it has at least 10 mutations
        if(nrow(mutation_info) >= 10)
        {
          total_plot_median <- total_plot_median + geom_density(data=mutation_info, aes(score_org_cancer_clines_median,  colour = tumor_type), alpha=0.1)
        } 
      }
    }, condition=function(ex) {
      print(ex);
    })
    
  }
  
  # plot density for all replication time scores from the entire genome
  replicationTimeScores_df_filt$tumor_type <- "overall"
  
  total_plot_median <-  total_plot_median + geom_density(data=replicationTimeScores_df_filt, aes(score_org_cancer_clines_median, linetype=tumor_type), size=1) + geom_vline(xintercept = median_replTime_median_cancer_clines, color = "black", size=0.5, linetype="dotted")
  ggsave(plot =total_plot_median, filename =  paste(resultsDir, "/density_plot_", filename_plot, ".pdf",sep=""))
}