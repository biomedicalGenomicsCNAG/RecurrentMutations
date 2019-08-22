
#################################
###    SEQUENCE LOGO PLOTS    ###
#################################

# R version (3.4.2)

# LIBRARIES
library(Biostrings)
library(ggplot2)
library(ggseqlogo)
prior_genome <- c(0.2952696, 0.2044550, 0.2045758, 0.2956996)
names(prior_genome) <- c("A","C","G","T")

# FUNCTION: "plotLogo21"

# Description: Plots sequence logo by using ggseqlogo (which in turn makes use of ggolot2)
 
# Arguments (this function is call in the next function and all the next parameters will be defined along the function)
##   data:character vector of sequences of the same length
##   sequence_length: length of the sequence to be plotted
##   yminrectsubtype: height of the bottom line for the box where the mutation subtype will be indicated.
##   ymaxrectsubtype: height of the top line for the box where the mutation subtype will be indicated.
##   ytext: height coordinate to plot the mutation subtype.
##   group: specific cluster wanted to be analysed (e.g. in case of Cluster A, provide: "cluster_A").
##   subtype: mutation type (e.g. for evaluating C>A mutations, provide: "C_A" ).
##   n_commas: number of sequences used to build the sequence logo (Formatted using comma for the thousand position)


plotLogo21 <- function(data, sequence_length, yminrectsubtype, ymaxrectsubtype, ytext, group, subtype, n_commas){
  
  basic_LogoPlot <- ggseqlogo(data, method="custom", seq_type='dna') +
    geom_hline(yintercept=0, size=0.3) +
    geom_vline(xintercept = seq(0.5,sequence_length+0.5,by=1), linetype="dotted",
               color = "gray", size=0.5) +
    scale_x_continuous(breaks = seq(1,sequence_length,by=1), labels = as.character(seq(-floor(sequence_length/2),floor(sequence_length/2),by=1))) +
    scale_y_continuous(limits = c(-0.5, 2.5), breaks=seq(-0.5,2.5,by=0.5)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 0, vjust=0.4, size = 20, hjust=0.5),
          axis.text.y = element_text(angle = 0, size = 20),
          axis.title.x=element_text(size=22),
          axis.title.y=element_text(size=22),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    annotate('rect', xmin = floor(sequence_length/2)+0.4, xmax = ceiling(sequence_length/2)+0.6, ymin = yminrectsubtype, ymax = ymaxrectsubtype, alpha = 1, col="white", fill="white") +
    annotate('rect', xmin = floor(sequence_length/2)+0.4, xmax = ceiling(sequence_length/2)+0.6, ymin = yminrectsubtype, ymax = ymaxrectsubtype, alpha = .1, col="black", fill="green") +
    annotate('text', x=ceiling(sequence_length/2), y=ytext, label=gsub("_",">",subtype), fontface="bold", size=5) +
    annotate('rect', xmin = sequence_length - 2.5, xmax = sequence_length+0.5, ymin = 2.40, ymax = 2.50, alpha = 1, col="white", fill="white") +
    annotate('text', x=sequence_length-1, y=2.45, label=paste0("n = ",n_commas), fontface="bold", size=5) +
    xlab("position with respect to the mutation") +
    ylab("relative entropy")
  
  return(basic_LogoPlot)
  
}



# FUNCTION: "SequenceLogo21"

# Description: Computes relative entropy per base in a position of the sequence and total entropy in that position from a matrix of counts of nucleotide per position,
#             builds the correspondig sequence Logo and highlights any position enriched.

# Arguments
##    file: file name of a matrix with the counts of each nucleotide pero position of the aligned sequences to stud. Needed a file in .RData format.
##    length_aroundMut: how many bases to extend next to the mutation (going in one direction, the same length will be used for the opposite direction).
##    threshold: minimum information content (measured in total entropy) to consider enrichment in a position of the sequence.
##    group: specific cluster wanted to be analysed (e.g. in case of Cluster A, provide: "cluster_A").
##    subtype: mutation type (e.g. for evaluating C>A mutations, provide: "C_A" ).
##    selection_case: "rec" if wanted to analyse recurrent mutations or "unique" if wanted to analyse non-recurrent mutations.
##    prior genome: corresponding proportion of each nucleotide in the genome.
##    Dir_toSave_logoplot: directory to save the sequence logo plot.



SequenceLogo21 <- function(file, length_aroundMut, threshold, group, subtype, selection_case, prior_genome, Dir_toSave_logoplot){

  # 1. LOAD DATA and assign some variables
    # data counts
  loaded <- load(file)
  data <- get(loaded)
  
  # total length of the sequence to plot
  sequence_length <- (length_aroundMut*2)+1
  
    # number of samples
  n <- data[which(data[,ceiling(sequence_length/2)] != 0), ceiling(sequence_length/2)]
  n <- max(n)  #n[[1]] #this is the number of sequences for this cluster and subtype (not anymore)
  n_commas <- formatC(n, big.mark=",")
  
    # mutation
  mutation <- gsub("_", ">", subtype)
    
    # some parameters for the actual plot
  if(subtype %in% c("T_A", "T_C", "T_G")){
    ymaxrect <- 1.76
    yminrectsubtype <- 1.78
    ymaxrectsubtype <- 1.88
    ytext <- 1.83
  }else{
    ymaxrect <- 2.3
    yminrectsubtype <- 2.32
    ymaxrectsubtype <- 2.42
    ytext <- 2.37
  }
  
  
  if(nrow(data) != 0){
    
    # 2. COMPUTE RELATIVE ENTROPIES AND TOTAL ENTROPY PER POSITION (Compute Relative Entropy per nucleotide per position and Total Entropy per position)
    
    matrix_total <- matrix(rep(n,sequence_length*4),nrow=4, ncol=sequence_length)
    prob_bases_seq_rev <- data/matrix_total
    prob_bases_seq_rev_prior <- prob_bases_seq_rev/prior_genome
    R_frequency_prior <- prob_bases_seq_rev * log2(prob_bases_seq_rev_prior)
    
    information_content <- colSums(R_frequency_prior, na.rm=TRUE)

    
    # 3. SELECT ENRICHED POSITIONS OVER DECIDED THRESHOLD 
    
    selected_positions <- which(information_content > threshold)
    

    # 4. PLOT LOGO
    
    if(length(selected_positions) > 1){
      
      # Check if the selection is gapped or ungapped, select parameters and plot logo
      continuous_pos <- seq(min(selected_positions), max(selected_positions),by=1)
      
      if(length(continuous_pos) == length(selected_positions)){    # (Case 1) Selected positions continuous -> UNGAPPED MOTIF

        # LOGO 21 positions. Function [5]
        basic_LogoPlot <- plotLogo21(R_frequency_prior, sequence_length, yminrectsubtype, ymaxrectsubtype, ytext, group, subtype, n_commas) #input data
        
        LogoPlot <- basic_LogoPlot + annotate('rect', xmin = selected_positions[1]-0.5, xmax = selected_positions[length(selected_positions)]+0.5, ymin = 0.0, ymax = ymaxrect, alpha = .1, col='black', fill='yellow')
        ggsave(paste0(Dir_toSave_logoplot,"LogoPlot_",group,"_",selection_case,"_",subtype,"_seq21.pdf"),  width = 30, height = 20, units = "cm", limitsize = FALSE)
        
        
      }else{                                                       # (Case 2) Selected positions non-continuous -> GAPPED MOTIF
        
        # Plot SequenceLogo of 21 positions
        basic_LogoPlot <- plotLogo21(R_frequency_prior, sequence_length, yminrectsubtype, ymaxrectsubtype, ytext, group, subtype, n_commas)
  
        # Selected positions
        Breaks <- c(0, which(diff(selected_positions) != 1), length(selected_positions)) 
        ind_selections <- sapply(seq(length(Breaks) - 1), function(i) selected_positions[(Breaks[i] + 1):Breaks[i+1]]) 
        
        # Adding rectangules for selected positions
        for (selection in ind_selections){
          basic_LogoPlot <- basic_LogoPlot + 
            annotate('rect', xmin = selection[1]-0.5, xmax = selection[length(selection)]+0.5, ymin = 0.0, ymax = ymaxrect, alpha = .1, col='black', fill='yellow')
        }
        ggsave(paste0(Dir_toSave_logoplot,"LogoPlot_",group,"_",selection_case,"_",subtype,"_seq21.pdf"),  width = 30, height = 20, units = "cm", limitsize = FALSE)
        
      }
      
      
    }else{                                                         # (Case 3) Normal plot with nothing selected
           
      # Plot SequenceLogo of 21 positions
      basic_LogoPlot <- plotLogo21(R_frequency_prior, sequence_length, yminrectsubtype, ymaxrectsubtype, ytext, group, subtype, n_commas)
        # save plot
      ggsave(paste0(Dir_toSave_logoplot,"LogoPlot_",group,"_",selection_case,"_",subtype,"_seq21.pdf"),  width = 30, height = 20, units = "cm", limitsize = FALSE)
      
    }
    
    
  }

  
}




## CALL THE FUNCTION: Example for recurrent C>T mutations in Cluster C.

length_aroundMut <- 10
threshold <- 0.25
group <- "cluster_G"
subtype <- "C_T"
selection_case <- "rec"
prior_genome <- c(0.2952696, 0.2044550, 0.2045758, 0.2956996)
names(prior_genome) <- c("A","C","G","T")
file <- paste0(path_to_file_directory,"/SeqContext_counts_",group,"_",selection_case,"_",subtype,"_21positions.RData")
Dir_toSave_logoplot <- getwd()

  # calling the function:
SequenceLogo21(file, length_aroundMut, threshold, group, subtype, selection_case, prior_genome, Dir_toSave_logoplot)

