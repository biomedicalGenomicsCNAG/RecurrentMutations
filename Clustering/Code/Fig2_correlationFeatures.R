library(corrplot)
library(Hmisc)

#' Get all correlations and correct for multiple-testing
#' @param sample2featuresAbsAndPerc: sample linked to the features in absolute number and in percentages
#' @param correlation_type: method to use to compute correlation
#' @param multiple_testing_correction: method to use to do multiple testing correction
getAllCorrelations <- function(sample2featuresAbsAndPerc, correlation_type, multiple_testing_correction){
 
  all_correlations <- as.data.frame(matrix(nrow=1010,ncol=5))
  colnames(all_correlations) <- c("set","feature_1", "feature_2", "correlation", "pval")
  
  #####################################################
  # Correlation between every pair of the 42 features #
  #####################################################
  
  poly_columns <- grep("Poly", colnames(sample2featuresAbsAndPerc), value=TRUE)
  poly_columns <- grep("perc_", poly_columns, value=TRUE)
  feature_names <- c("perc_rec_ssms", "perc_rec_sims", "perc_of_rec_mut_sim", "perc_rec_CA", "perc_rec_CG", "perc_rec_CT", "perc_rec_TA", "perc_rec_TC", "perc_rec_TG", "perc_rec_1bp_del_A_T", "perc_rec_1bp_del_C_G", "perc_rec_1bp_ins_A_T", "perc_rec_1bp_ins_C_G","num_ssms","num_sims", "perc_of_mut_sims", "perc_CA", "perc_CG", "perc_CT", "perc_TA", "perc_TC", "perc_TG", "perc_1bp_del_A_T", "perc_1bp_del_C_G", "perc_1bp_ins_A_T", "perc_1bp_ins_C_G", poly_columns)
  
  curRow <- 1
  
  for(i in 1:(length(feature_names)-1)){
    
    curStart <- i + 1
    
    for(j in curStart:length(feature_names)){
      
      corr_res <- rcorr(sample2featuresAbsAndPerc[,feature_names[i] ],sample2featuresAbsAndPerc[,feature_names[j] ], type=correlation_type)
      
      all_correlations[curRow,"feature_1"] <-feature_names[i] 
      all_correlations[curRow,"feature_2"] <-feature_names[j]
      all_correlations[curRow,"set"] <-"cohort"
      all_correlations[curRow,"pval"] <-corr_res$P[1,2]
      all_correlations[curRow,"correlation"] <- corr_res$r[1,2]
      curRow <- curRow + 1
    }
  }
  
  ############################################################################
  # Correlation between total and (number/percentage of) recurrent mutations #
  ############################################################################
  
  # Correlation between total number and recurrent SSMs
  correlation_num_rec_ssms <- rcorr(sample2featuresAbsAndPerc[,"num_ssms"],sample2featuresAbsAndPerc[,"num_rec_ssms"], type=correlation_type)
  pval_correlation_num_rec_ssms <- correlation_num_rec_ssms$P[1,2]
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"feature_1"] <- "num_ssms"
  all_correlations[curRow,"feature_2"] <- "num_rec_ssms"
  all_correlations[curRow,"pval"] <- pval_correlation_num_rec_ssms
  all_correlations[curRow,"correlation"] <- correlation_num_rec_ssms$r[1,2]
  curRow <- curRow + 1
  
  # Correlation between total number and recurrent SIMs
  correlation_num_rec_sims <- rcorr(sample2featuresAbsAndPerc[,"num_sims"],sample2featuresAbsAndPerc[,"num_rec_sims"], type=correlation_type)
  pval_correlation_num_rec_sims <- correlation_num_rec_sims$P[1,2]
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"feature_1"] <- "num_sims"
  all_correlations[curRow,"feature_2"] <- "num_rec_sims"
  all_correlations[curRow,"pval"] <- pval_correlation_num_rec_sims
  all_correlations[curRow,"correlation"] <- correlation_num_rec_sims$r[1,2]
  curRow <- curRow + 1
  
  # Correlation between total number and percentage of recurrent 1 bp SIMs
  sample2featuresAbsAndPerc$num_1bp_sims <- sample2featuresAbsAndPerc$num_1bp_del_A_T + sample2featuresAbsAndPerc$num_1bp_del_C_G + sample2featuresAbsAndPerc$num_1bp_ins_A_T + sample2featuresAbsAndPerc$num_1bp_ins_C_G
  sample2featuresAbsAndPerc$num_rec_1bp_sims <- sample2featuresAbsAndPerc$num_rec_1bp_del_A_T + sample2featuresAbsAndPerc$num_rec_1bp_del_C_G + sample2featuresAbsAndPerc$num_rec_1bp_ins_A_T + sample2featuresAbsAndPerc$num_rec_1bp_ins_C_G
  sample2featuresAbsAndPerc$perc_rec_1bp_sims <- (sample2featuresAbsAndPerc$num_rec_1bp_sims*100)/sample2featuresAbsAndPerc$num_1bp_sims 
  sample2featuresAbsAndPerc$perc_rec_1bp_sims[which(is.na(sample2featuresAbsAndPerc$perc_rec_1bp_sims))] <- 0 
  
  correlation_perc_rec_1bp_sims <- rcorr(sample2featuresAbsAndPerc[,"num_1bp_sims"],sample2featuresAbsAndPerc[,"perc_rec_1bp_sims"], type=correlation_type)
  pval_correlation_perc_rec_1bp_sims <- correlation_perc_rec_1bp_sims$P[1,2]
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"feature_1"] <- "num_1bp_sims"
  all_correlations[curRow,"feature_2"] <- "perc_rec_1bp_sims"
  all_correlations[curRow,"pval"] <- pval_correlation_perc_rec_1bp_sims
  all_correlations[curRow,"correlation"] <- correlation_perc_rec_1bp_sims$r[1,2]
  curRow <- curRow + 1
  
  #########################################################
  # number of SSM subtypes vs perc recurrent SSM subtypes #
  #########################################################
  
  # C>A SSMs
  correlation_perc_rec_CA <- rcorr(sample2featuresAbsAndPerc[,"num_CA"],sample2featuresAbsAndPerc[,"perc_rec_CA"], type=correlation_type)
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"feature_1"] <- "num_CA"
  all_correlations[curRow,"feature_2"] <- "perc_rec_CA"
  all_correlations[curRow,"pval"] <- correlation_perc_rec_CA$P[1,2]
  all_correlations[curRow,"correlation"] <- correlation_perc_rec_CA$r[1,2]
  curRow <- curRow + 1
  
  # C>G SSMs
  correlation_perc_rec_CG <- rcorr(sample2featuresAbsAndPerc[,"num_CG"],sample2featuresAbsAndPerc[,"perc_rec_CG"], type=correlation_type)
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"feature_1"] <- "num_CG"
  all_correlations[curRow,"feature_2"] <- "perc_rec_CG"
  all_correlations[curRow,"pval"] <- correlation_perc_rec_CG$P[1,2]
  all_correlations[curRow,"correlation"] <- correlation_perc_rec_CG$r[1,2]
  curRow <- curRow + 1
  
  # C>T SSMs
  correlation_perc_rec_CT <- rcorr(sample2featuresAbsAndPerc[,"num_CT"],sample2featuresAbsAndPerc[,"perc_rec_CT"], type=correlation_type)
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"feature_1"] <- "num_CT"
  all_correlations[curRow,"feature_2"] <- "perc_rec_CT"
  all_correlations[curRow,"pval"] <- correlation_perc_rec_CT$P[1,2]
  all_correlations[curRow,"correlation"] <- correlation_perc_rec_CT$r[1,2]
  curRow <- curRow + 1
  
  # T>A SSMs
  correlation_num_rec_TA <- rcorr(sample2featuresAbsAndPerc[,"num_TA"],sample2featuresAbsAndPerc[,"perc_rec_TA"], type=correlation_type)
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"feature_1"] <- "num_TA"
  all_correlations[curRow,"feature_2"] <- "perc_rec_TA"
  all_correlations[curRow,"pval"] <- correlation_num_rec_TA$P[1,2]
  all_correlations[curRow,"correlation"] <- correlation_num_rec_TA$r[1,2]
  curRow <- curRow + 1
  
  # T>C SSMs
  correlation_perc_rec_TC <- rcorr(sample2featuresAbsAndPerc[,"num_TC"],sample2featuresAbsAndPerc[,"perc_rec_TC"], type=correlation_type)
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"feature_1"] <- "num_TC"
  all_correlations[curRow,"feature_2"] <- "perc_rec_TC"
  all_correlations[curRow,"pval"] <- correlation_perc_rec_TC$P[1,2]
  all_correlations[curRow,"correlation"] <- correlation_perc_rec_TC$r[1,2]
  curRow <- curRow + 1
  
  # T>G SSMs
  correlation_perc_rec_TG <- rcorr(sample2featuresAbsAndPerc[,"num_TG"],sample2featuresAbsAndPerc[,"perc_rec_TG"], type=correlation_type)
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"feature_1"] <- "num_TG"
  all_correlations[curRow,"feature_2"] <- "perc_rec_TG"
  all_correlations[curRow,"pval"] <- correlation_perc_rec_TG$P[1,2]
  all_correlations[curRow,"correlation"] <- correlation_perc_rec_TG$r[1,2]
  curRow <- curRow + 1
  
  #########################################################
  # number of SIM subtypes vs perc recurrent SIM subtypes #
  #########################################################
  
  # 1 bp A/T deletions
  correlation_perc_rec_delAT <- rcorr(sample2featuresAbsAndPerc[,"num_1bp_del_A_T"],sample2featuresAbsAndPerc[,"perc_rec_1bp_del_A_T"], type=correlation_type)
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"feature_1"] <- "num_1bp_del_A_T"
  all_correlations[curRow,"feature_2"] <- "perc_rec_1bp_del_A_T"
  all_correlations[curRow,"pval"] <- correlation_perc_rec_delAT$P[1,2]
  all_correlations[curRow,"correlation"] <- correlation_perc_rec_delAT$r[1,2]
  curRow <- curRow + 1
  
  # 1 bp C/G deletions
  correlation_perc_rec_delCG <- rcorr(sample2featuresAbsAndPerc[,"num_1bp_del_C_G"],sample2featuresAbsAndPerc[,"perc_rec_1bp_del_C_G"], type=correlation_type)
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"feature_1"] <- "num_1bp_del_C_G"
  all_correlations[curRow,"feature_2"] <- "perc_rec_1bp_del_C_G"
  all_correlations[curRow,"pval"] <- correlation_perc_rec_delCG$P[1,2]
  all_correlations[curRow,"correlation"] <- correlation_perc_rec_delCG$r[1,2]
  curRow <- curRow + 1
  
  # 1 bp A/T insertions
  correlation_perc_rec_insAT <- rcorr(sample2featuresAbsAndPerc[,"num_1bp_ins_A_T"],sample2featuresAbsAndPerc[,"perc_rec_1bp_ins_A_T"], type=correlation_type)
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"feature_1"] <- "num_1bp_ins_A_T"
  all_correlations[curRow,"feature_2"] <- "perc_rec_1bp_ins_A_T"
  all_correlations[curRow,"pval"] <- correlation_perc_rec_insAT$P[1,2]
  all_correlations[curRow,"correlation"] <- correlation_perc_rec_insAT$r[1,2]
  curRow <- curRow + 1
  
  # 1 bp C/G insertions
  correlation_perc_rec_insCG <- rcorr(sample2featuresAbsAndPerc[,"num_1bp_ins_C_G"],sample2featuresAbsAndPerc[,"perc_rec_1bp_ins_C_G"], type=correlation_type)
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"feature_1"] <- "num_1bp_ins_C_G"
  all_correlations[curRow,"feature_2"] <- "perc_rec_1bp_ins_C_G"
  all_correlations[curRow,"pval"] <- correlation_perc_rec_insCG$P[1,2]
  all_correlations[curRow,"correlation"] <- correlation_perc_rec_insCG$r[1,2]
  curRow <- curRow + 1
  
  ###################################################################################################
  # Correlation between total number and number & percentage recurrent SSMs and SIMs per tumor type #
  ###################################################################################################
  
  ttype2corr_num_perc_rec <- as.data.frame(matrix(nrow=37, ncol=9))
  colnames(ttype2corr_num_perc_rec) <- c("tumor_type","corr_num_ssms", "pval_num_ssms","corr_num_sims", "pval_num_sims","corr_perc_ssms", "pval_perc_ssms","corr_perc_sims", "pval_perc_sims")
  
  tumor_types <- sort(unique(sample2featuresAbsAndPerc$tumor_type))
  
  
  for(i in 1:length(tumor_types)){
    
    print(i)
    cur_ttype <- sample2featuresAbsAndPerc[which(sample2featuresAbsAndPerc$tumor_type == tumor_types[i]),]
    
    ttype2corr_num_perc_rec[i, "tumor_type"] <- tumor_types[i]
    
    # if the tumor type has more than 4 samples (minimum for rcorr function)
    if(nrow(cur_ttype) > 4)
    {
      # correlation between number of SSMs and number of recurrent SSMs
      correlation_num_rec_ssms <- rcorr(cur_ttype[,"num_ssms"],cur_ttype[,"num_rec_ssms"], type=correlation_type)
      all_correlations[curRow,"set"] <-tumor_types[i]
      all_correlations[curRow,"feature_1"] <- "num_ssms"
      all_correlations[curRow,"feature_2"] <- "num_rec_ssms"
      all_correlations[curRow,"pval"] <- correlation_num_rec_ssms$P[1,2]
      all_correlations[curRow,"correlation"] <- correlation_num_rec_ssms$r[1,2]
      curRow <- curRow + 1
      
      # correlation between number of SSMs and percentage of recurrent SSMs
      correlation_perc_rec_ssms <- rcorr(cur_ttype[,"num_ssms"],cur_ttype[,"perc_rec_ssms"], type=correlation_type)
      all_correlations[curRow,"set"] <-tumor_types[i]
      all_correlations[curRow,"feature_1"] <- "num_ssms"
      all_correlations[curRow,"feature_2"] <- "perc_rec_ssms"
      all_correlations[curRow,"pval"] <- correlation_perc_rec_ssms$P[1,2]
      all_correlations[curRow,"correlation"] <- correlation_perc_rec_ssms$r[1,2]
      curRow <- curRow + 1
      
      # correlation between number of SIMs and number of recurrent SIMs
      correlation_num_sims <- rcorr(cur_ttype[,"num_sims"],cur_ttype[,"num_rec_sims"], type=correlation_type)
      all_correlations[curRow,"set"] <-tumor_types[i]
      all_correlations[curRow,"feature_1"] <- "num_sims"
      all_correlations[curRow,"feature_2"] <- "num_rec_sims"
      all_correlations[curRow,"pval"] <- correlation_num_sims$P[1,2]
      all_correlations[curRow,"correlation"] <- correlation_num_sims$r[1,2]
      curRow <- curRow + 1
      
      # correlation between number of SIMs and percentage of recurrent SIMs
      correlation_perc_sims <- rcorr(cur_ttype[,"num_sims"],cur_ttype[,"perc_rec_sims"], type=correlation_type)
      all_correlations[curRow,"set"] <-tumor_types[i]
      all_correlations[curRow,"feature_1"] <- "num_sims"
      all_correlations[curRow,"feature_2"] <- "perc_rec_sims"
      all_correlations[curRow,"pval"] <- correlation_perc_sims$P[1,2]
      all_correlations[curRow,"correlation"] <- correlation_perc_sims$r[1,2]
      curRow <- curRow + 1
      
      
    }
  }
  
  
  #########################################
  # Correct p-values for multiple testing #
  #########################################
  
  # replace all zero p-values with the smallest possible value
  all_correlations[which(all_correlations$pval == 0), "pval"] <- noquote(unlist(format(.Machine)))["double.eps"]
  
  # correct for multiple testing
  all_correlations[,"pval_adj"] <- p.adjust(all_correlations$pval, method =multiple_testing_correction)
  all_correlations$correlation_rounded <- round(all_correlations$correlation,2)
  
  return(all_correlations)
} 

#' Plot the correlation between the 42 features
#' @param all_correlations: all correlations computed for the manuscript with p-values adjusted for multiple testing
#' @param features_names: name of the 42 features
#' @param resultsDir: directory to store the plot
plotCorrelationBetweenFeatures <- function(all_correlations, features_names, resultsDir)
{
  
  # Get the correlations and adjusted p-values for the 42 features
  adjusted_pvals <- matrix(nrow=length(features_names), ncol=length(features_names))
  colnames(adjusted_pvals) <- features_names
  rownames(adjusted_pvals) <- features_names
  
  correlations <- matrix(nrow=length(features_names), ncol=length(features_names))
  colnames(correlations) <- features_names
  rownames(correlations) <- features_names
  
  for(i in 1:(length(features_names)-1)){
    
    curStart <- i + 1
    
    for(j in curStart:length(features_names)){
      
      adjusted_pvals[features_names[i], features_names[j]] <- all_correlations[which(all_correlations$feature_1 == features_names[i] & all_correlations$feature_2 == features_names[j]), "pval_adj"]
      correlations[features_names[i], features_names[j]] <- all_correlations[which(all_correlations$feature_1 == features_names[i] & all_correlations$feature_2 == features_names[j]), "correlation"]
    }
  }
  
  adjusted_pvals[lower.tri(adjusted_pvals)] <- t(adjusted_pvals)[lower.tri(adjusted_pvals)]
  
  # plot the correlations
  pdf(paste(resultsDir, "correlation_plot_features.pdf",sep=""))
  
  corrplot(correlations, order="original", p.mat = adjusted_pvals, sig.level = 0.05,pch.cex=0.5, tl.cex=0.5, tl.col="black")
  
  dev.off()
}