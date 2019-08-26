

getAllCorrelations <- function(sample2stats){
 
  all_correlations <- as.data.frame(matrix(nrow=1010,ncol=5))
  colnames(all_correlations) <- c("set","stat_1", "stat_2", "correlation", "pval")
  
  # Correlation between every pair of features
  poly_columns <- grep("Poly", colnames(sample2stats), value=TRUE)
  poly_columns <- grep("perc_", poly_columns, value=TRUE)
  feature_names <- c("perc_rec_ssms", "perc_rec_sims", "perc_of_rec_mut_sim", "perc_rec_CA", "perc_rec_CG", "perc_rec_CT", "perc_rec_TA", "perc_rec_TC", "perc_rec_TG", "perc_rec_1bp_del_A_T", "perc_rec_1bp_del_C_G", "perc_rec_1bp_ins_A_T", "perc_rec_1bp_ins_C_G","num_ssms","num_sims", "perc_of_mut_sims", "perc_CA", "perc_CG", "perc_CT", "perc_TA", "perc_TC", "perc_TG", "perc_1bp_del_A_T", "perc_1bp_del_C_G", "perc_1bp_ins_A_T", "perc_1bp_ins_C_G", poly_columns)
  
  curRow <- 1
  
  for(i in 1:(length(feature_names)-1)){
    
    curStart <- i + 1
    
    for(j in curStart:length(feature_names)){
      
      corr_res <- rcorr(sample2stats[,feature_names[i] ],sample2stats[,feature_names[j] ], type="spearman")
      
      all_correlations[curRow,"stat_1"] <-feature_names[i] 
      all_correlations[curRow,"stat_2"] <-feature_names[j]
      all_correlations[curRow,"set"] <-"cohort"
      all_correlations[curRow,"pval"] <-corr_res$P[1,2]
      all_correlations[curRow,"correlation"] <- corr_res$r[1,2]
      curRow <- curRow + 1
    }
  }
  
  # Correlation between total number and recurrent SSMs and SIMs
  correlation_num_rec_ssms <- rcorr(sample2stats[,"num_ssms"],sample2stats[,"num_rec_ssms"], type="spearman")
  pval_correlation_num_rec_ssms <- correlation_num_rec_ssms$P[1,2]
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"stat_1"] <- "num_ssms"
  all_correlations[curRow,"stat_2"] <- "num_rec_ssms"
  all_correlations[curRow,"pval"] <- pval_correlation_num_rec_ssms
  all_correlations[curRow,"correlation"] <- correlation_num_rec_ssms$r[1,2]
  curRow <- curRow + 1
  
  correlation_num_rec_sims <- rcorr(sample2stats[,"num_sims"],sample2stats[,"num_rec_sims"], type="spearman")
  pval_correlation_num_rec_sims <- correlation_num_rec_sims$P[1,2]
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"stat_1"] <- "num_sims"
  all_correlations[curRow,"stat_2"] <- "num_rec_sims"
  all_correlations[curRow,"pval"] <- pval_correlation_num_rec_sims
  all_correlations[curRow,"correlation"] <- correlation_num_rec_sims$r[1,2]
  curRow <- curRow + 1
  
  # Correlation between total number and percentage of recurrent 1 bp SIMs
  sample2stats$num_1bp_sims <- sample2stats$num_1bp_del_A_T + sample2stats$num_1bp_del_C_G + sample2stats$num_1bp_ins_A_T + sample2stats$num_1bp_ins_C_G
  sample2stats$num_rec_1bp_sims <- sample2stats$num_rec_1bp_del_A_T + sample2stats$num_rec_1bp_del_C_G + sample2stats$num_rec_1bp_ins_A_T + sample2stats$num_rec_1bp_ins_C_G
  sample2stats$perc_rec_1bp_sims <- (sample2stats$num_rec_1bp_sims*100)/sample2stats$num_1bp_sims 
  sample2stats$perc_rec_1bp_sims[which(is.na(sample2stats$perc_rec_1bp_sims))] <- 0 
  
  correlation_perc_rec_1bp_sims <- rcorr(sample2stats[,"num_1bp_sims"],sample2stats[,"perc_rec_1bp_sims"], type="spearman")
  pval_correlation_perc_rec_1bp_sims <- correlation_perc_rec_1bp_sims$P[1,2]
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"stat_1"] <- "num_1bp_sims"
  all_correlations[curRow,"stat_2"] <- "perc_rec_1bp_sims"
  all_correlations[curRow,"pval"] <- pval_correlation_perc_rec_1bp_sims
  all_correlations[curRow,"correlation"] <- correlation_perc_rec_1bp_sims$r[1,2]
  curRow <- curRow + 1
  
  # Correlation between total number and number/percentage recurrent SSMs and SIMs per tumor type
  ttype2corr_num_perc_rec <- as.data.frame(matrix(nrow=37, ncol=9))
  colnames(ttype2corr_num_perc_rec) <- c("tumor_type","corr_num_ssms", "pval_num_ssms","corr_num_sims", "pval_num_sims","corr_perc_ssms", "pval_perc_ssms","corr_perc_sims", "pval_perc_sims")
  
  tumor_types <- sort(unique(sample2stats$tumor_type))
  
  
  for(i in 1:length(tumor_types)){
    
    print(i)
    cur_ttype <- sample2stats[which(sample2stats$tumor_type == tumor_types[i]),]
    
    ttype2corr_num_perc_rec[i, "tumor_type"] <- tumor_types[i]
    
    # if the tumor type has more than 4 samples (minimum for rcorr function)
    if(nrow(cur_ttype) > 4)
    {
      correlation_num_rec_ssms <- rcorr(cur_ttype[,"num_ssms"],cur_ttype[,"num_rec_ssms"], type="spearman")
      ttype2corr_num_perc_rec[i, "corr_num_ssms"] <- correlation_num_rec_ssms$r[1,2]
      ttype2corr_num_perc_rec[i, "pval_num_ssms"] <- correlation_num_rec_ssms$P[1,2]
      
      all_correlations[curRow,"set"] <-tumor_types[i]
      all_correlations[curRow,"stat_1"] <- "num_ssms"
      all_correlations[curRow,"stat_2"] <- "num_rec_ssms"
      all_correlations[curRow,"pval"] <- correlation_num_rec_ssms$P[1,2]
      all_correlations[curRow,"correlation"] <- correlation_num_rec_ssms$r[1,2]
      curRow <- curRow + 1
      
      correlation_perc_rec_ssms <- rcorr(cur_ttype[,"num_ssms"],cur_ttype[,"perc_rec_ssms"], type="spearman")
      ttype2corr_num_perc_rec[i, "corr_perc_ssms"] <- correlation_perc_rec_ssms$r[1,2]
      ttype2corr_num_perc_rec[i, "pval_perc_ssms"] <- correlation_perc_rec_ssms$P[1,2]
      
      all_correlations[curRow,"set"] <-tumor_types[i]
      all_correlations[curRow,"stat_1"] <- "num_ssms"
      all_correlations[curRow,"stat_2"] <- "perc_rec_ssms"
      all_correlations[curRow,"pval"] <- correlation_perc_rec_ssms$P[1,2]
      all_correlations[curRow,"correlation"] <- correlation_perc_rec_ssms$r[1,2]
      curRow <- curRow + 1
      
      correlation_num_sims <- rcorr(cur_ttype[,"num_sims"],cur_ttype[,"num_rec_sims"], type="spearman")
      ttype2corr_num_perc_rec[i, "corr_num_sims"] <- correlation_num_sims$r[1,2]
      ttype2corr_num_perc_rec[i, "pval_num_sims"] <- correlation_num_sims$P[1,2]
      
      all_correlations[curRow,"set"] <-tumor_types[i]
      all_correlations[curRow,"stat_1"] <- "num_sims"
      all_correlations[curRow,"stat_2"] <- "num_rec_sims"
      all_correlations[curRow,"pval"] <- correlation_num_sims$P[1,2]
      all_correlations[curRow,"correlation"] <- correlation_num_sims$r[1,2]
      curRow <- curRow + 1
      
      correlation_perc_sims <- rcorr(cur_ttype[,"num_sims"],cur_ttype[,"perc_rec_sims"], type="spearman")
      ttype2corr_num_perc_rec[i, "corr_perc_sims"] <- correlation_perc_sims$r[1,2]
      ttype2corr_num_perc_rec[i, "pval_perc_sims"] <- correlation_perc_sims$P[1,2]
      
      all_correlations[curRow,"set"] <-tumor_types[i]
      all_correlations[curRow,"stat_1"] <- "num_sims"
      all_correlations[curRow,"stat_2"] <- "perc_rec_sims"
      all_correlations[curRow,"pval"] <- correlation_perc_sims$P[1,2]
      all_correlations[curRow,"correlation"] <- correlation_perc_sims$r[1,2]
      curRow <- curRow + 1
      
      
    }
  }
  
  # number of SSM subtypes vs perc recurrent SSM subtypes
  correlation_perc_rec_CA <- rcorr(sample2stats[,"num_CA"],sample2stats[,"perc_rec_CA"], type="spearman")
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"stat_1"] <- "num_CA"
  all_correlations[curRow,"stat_2"] <- "perc_rec_CA"
  all_correlations[curRow,"pval"] <- correlation_perc_rec_CA$P[1,2]
  all_correlations[curRow,"correlation"] <- correlation_perc_rec_CA$r[1,2]
  curRow <- curRow + 1
  
  correlation_perc_rec_CG <- rcorr(sample2stats[,"num_CG"],sample2stats[,"perc_rec_CG"], type="spearman")
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"stat_1"] <- "num_CG"
  all_correlations[curRow,"stat_2"] <- "perc_rec_CG"
  all_correlations[curRow,"pval"] <- correlation_perc_rec_CG$P[1,2]
  all_correlations[curRow,"correlation"] <- correlation_perc_rec_CG$r[1,2]
  curRow <- curRow + 1
  
  correlation_perc_rec_CT <- rcorr(sample2stats[,"num_CT"],sample2stats[,"perc_rec_CT"], type="spearman")
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"stat_1"] <- "num_CT"
  all_correlations[curRow,"stat_2"] <- "perc_rec_CT"
  all_correlations[curRow,"pval"] <- correlation_perc_rec_CT$P[1,2]
  all_correlations[curRow,"correlation"] <- correlation_perc_rec_CT$r[1,2]
  curRow <- curRow + 1
  
  correlation_num_rec_TA <- rcorr(sample2stats[,"num_TA"],sample2stats[,"perc_rec_TA"], type="spearman")
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"stat_1"] <- "num_TA"
  all_correlations[curRow,"stat_2"] <- "perc_rec_TA"
  all_correlations[curRow,"pval"] <- correlation_num_rec_TA$P[1,2]
  all_correlations[curRow,"correlation"] <- correlation_num_rec_TA$r[1,2]
  curRow <- curRow + 1
  
  correlation_perc_rec_TC <- rcorr(sample2stats[,"num_TC"],sample2stats[,"perc_rec_TC"], type="spearman")
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"stat_1"] <- "num_TC"
  all_correlations[curRow,"stat_2"] <- "perc_rec_TC"
  all_correlations[curRow,"pval"] <- correlation_perc_rec_TC$P[1,2]
  all_correlations[curRow,"correlation"] <- correlation_perc_rec_TC$r[1,2]
  curRow <- curRow + 1
  
  correlation_perc_rec_TG <- rcorr(sample2stats[,"num_TG"],sample2stats[,"perc_rec_TG"], type="spearman")
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"stat_1"] <- "num_TG"
  all_correlations[curRow,"stat_2"] <- "perc_rec_TG"
  all_correlations[curRow,"pval"] <- correlation_perc_rec_TG$P[1,2]
  all_correlations[curRow,"correlation"] <- correlation_perc_rec_TG$r[1,2]
  curRow <- curRow + 1
  
  # number of SIM subtypes vs perc recurrent SIM subtypes
  correlation_perc_rec_delAT <- rcorr(sample2stats[,"num_1bp_del_A_T"],sample2stats[,"perc_rec_1bp_del_A_T"], type="spearman")
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"stat_1"] <- "num_1bp_del_A_T"
  all_correlations[curRow,"stat_2"] <- "perc_rec_1bp_del_A_T"
  all_correlations[curRow,"pval"] <- correlation_perc_rec_delAT$P[1,2]
  all_correlations[curRow,"correlation"] <- correlation_perc_rec_delAT$r[1,2]
  curRow <- curRow + 1
  
  correlation_perc_rec_delCG <- rcorr(sample2stats[,"num_1bp_del_C_G"],sample2stats[,"perc_rec_1bp_del_C_G"], type="spearman")
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"stat_1"] <- "num_1bp_del_C_G"
  all_correlations[curRow,"stat_2"] <- "perc_rec_1bp_del_C_G"
  all_correlations[curRow,"pval"] <- correlation_perc_rec_delCG$P[1,2]
  all_correlations[curRow,"correlation"] <- correlation_perc_rec_delCG$r[1,2]
  curRow <- curRow + 1
  
  correlation_perc_rec_insAT <- rcorr(sample2stats[,"num_1bp_ins_A_T"],sample2stats[,"perc_rec_1bp_ins_A_T"], type="spearman")
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"stat_1"] <- "num_1bp_ins_A_T"
  all_correlations[curRow,"stat_2"] <- "perc_rec_1bp_ins_A_T"
  all_correlations[curRow,"pval"] <- correlation_perc_rec_insAT$P[1,2]
  all_correlations[curRow,"correlation"] <- correlation_perc_rec_insAT$r[1,2]
  curRow <- curRow + 1
  
  correlation_perc_rec_insCG <- rcorr(sample2stats[,"num_1bp_ins_C_G"],sample2stats[,"perc_rec_1bp_ins_C_G"], type="spearman")
  all_correlations[curRow,"set"] <-"cohort"
  all_correlations[curRow,"stat_1"] <- "num_1bp_ins_C_G"
  all_correlations[curRow,"stat_2"] <- "perc_rec_1bp_ins_C_G"
  all_correlations[curRow,"pval"] <- correlation_perc_rec_insCG$P[1,2]
  all_correlations[curRow,"correlation"] <- correlation_perc_rec_insCG$r[1,2]
  curRow <- curRow + 1
  
  ######
  #All p values
  
  all_correlations[which(all_correlations$pval == 0), "pval"] <- noquote(unlist(format(.Machine)))["double.eps"]
  all_correlations[,"pval_adj"] <- p.adjust(all_correlations$pval, method ="BY")
  all_correlations$correlation_rounded <- round(all_correlations$correlation,2)
  
  return(all_correlations)
} 

plotCorrelationBetweenFeatures <- function(all_correlations, features_names, resultsDir){
  
  
  # Adjusted with all correlations
  corrected_pval_all <- matrix(nrow=length(features_names), ncol=length(features_names))
  colnames(corrected_pval_all) <- features_names
  rownames(corrected_pval_all) <- features_names
  
  for(i in 1:(length(features_names)-1)){
    
    curStart <- i + 1
    
    for(j in curStart:length(features_names)){
      
      corrected_pval_all[features_names[i], features_names[j]] <- all_correlations[which(all_correlations$stat_1 == features_names[i] & all_correlations$stat_2 == features_names[j]), "pval_adj"]
    }
  }
  
  corrected_pval_all[lower.tri(corrected_pval_all)] <- t(corrected_pval_all)[lower.tri(corrected_pval_all)]
  
  pdf(paste(resultsDir, "correlation_plot_features.pdf",sep=""))
  
  corrplot(corr_features_spearman$r, order="original", p.mat = corrected_pval_all, sig.level = 0.05,pch.cex=0.5, tl.cex=0.5, tl.col="black")
  dev.off()
}