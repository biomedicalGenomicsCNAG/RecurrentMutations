library(rtracklayer)
library(parallel)

source("main_utils.R")

#' Get the genes with a SSM/SIM ('mut_type') affecting the coding sequence
#' @param mut_with_effect: mutations that change coding sequence
#' @param mut_type: SSM or SIM
#' @param annotation_v19: GENCODE annotation for GRCh37/h19
#' @return genes_affected: list of genes with mutation from the 'mut_with_effect' list.
getAffectedGenes <- function(mut_with_effect, mut_type, annotation_v19)
{
  # get matches of the list of mutation with the GENODE annotation
  mutations_granges <- convert2GRanges(mut_with_effect, mut_type)
  
  matches_annotation <- findOverlaps(query = mutations_granges, subject = annotation_v19)
  
  mutations_withMatch <- as.data.frame(mut_with_effect[queryHits(matches_annotation),])
  annotation_matches <- as.data.frame(annotation[subjectHits(matches_annotation),])
  
  mutations_annotated <- unique(cbind(mutations_withMatch, annotation_matches[,c("gene_name", "type")]))
  
  mutation_notMatched <- as.data.frame(mut_with_effect[-queryHits(matches_annotation),])
  
  # if a mutation hits a CDS,then  get the corresponding gene
  genesCDS2mut_barcode <- unique(mutations_annotated[which(mutations_annotated$CDS == TRUE & mutations_annotated$type == "CDS"), c("mut_barcode", "gene_name")])
  
  # get the mutations that do not hit a CDS according to GENCODE
  mut_not_matched_CDS <- mutations_annotated[which(!(mutations_annotated$mut_barcode %in% genesCDS2mut_barcode$mut_barcode)),]
  
  # if a mutation that does not hit a CDS, but does hit an exon, then get the corresponding gene
  genesExon2mut_barcode <- unique(mut_not_matched_CDS[which(mut_not_matched_CDS$exon == TRUE & mut_not_matched_CDS$type == "exon"), c("mut_barcode", "gene_name")])
  
  # get the mutations that do not hit a CDS nor an exon according to GENCODE
  mut_not_matched_CDS_exon <- mut_not_matched_CDS[which(!(mut_not_matched_CDS$mut_barcode %in% genesExon2mut_barcode$mut_barcode)),]
  
  # if a mutation that does not hit a CDS nor exon, but does hit a transcript, then get the corresponding gene
  genesTranscript2mut_barcode <- unique(mut_not_matched_CDS_exon[which(mut_not_matched_CDS_exon$transcript == TRUE & mut_not_matched_CDS_exon$type == "transcript"), c("mut_barcode", "gene_name")])
  
  # get the mutations that do not hit a CDS nor an exon nor a transcript according to GENCODE
  mut_not_matched_CDS_exon_tx <- mut_not_matched_CDS_exon[which(!(mut_not_matched_CDS_exon$mut_barcode %in% genesTranscript2mut_barcode$mut_barcode)),]
  
  # if a mutation that does not hit a CDS nor exon nor transcript, but does hit a gene, then get the corresponding gene
  genesGene2mut_barcode <- unique(mut_not_matched_CDS_exon_tx[which(mut_not_matched_CDS_exon_tx$gene == TRUE & mut_not_matched_CDS_exon_tx$type == "gene"), c("mut_barcode", "gene_name")])
  
  if(nrow(mutation_notMatched) > 0)
  {
    print("mutations not matching in GENCODE --> intergenic")
    print(mutation_notMatched)
  }
  
  # combine list of genes hit and return unique list
  genes_affected <- unique(c(genesCDS2mut_barcode$gene_name, genesExon2mut_barcode$gene_name,genesTranscript2mut_barcode$gene_name,genesGene2mut_barcode$gene_name))
  
  return(genes_affected)
}

#' Summarize the annotation to sample level: impact classification, functional category and replication time scores
#' @param sample2annotation
#' @param mut_type: SSM or SIM
#' @param annotation_v19: GENCODE 
#' @param replicationTimeScores_df: replication time scores
#' @param annSamplesDir: directory of the annotated samples at mutation level
#' @param annSamplesFolder: folder of the annotated samples at mutation level
#' @param num_cores: number of cores to use for mclapply
summarizeMutAnnotation2SampleLevel <- function(sample2annotation, mut_type, annotation_v19, replicationTimeScores_df, annSamplesDir, annSamplesFolder, num_cores)
{
  samples_annotatedMutations <- list.files(paste(annSamplesDir, "/", annSamplesFolder, "/",sep=""), pattern = "_annotatedWithGenCode_ReplTime.txt")

  # cut offs to be used to separate early from late replication regions.
  replicationTimeScores_df$score_org_cancer_clines_median <- rowMedians(as.matrix(replicationTimeScores_df[, c("score_org_Helas3_1", "score_org_Mcf7_1", "score_org_Sknsh_1", "score_org_Hepg2_1", "score_org_K562_1")]))
  median_replTime_median_cancer_clines <- median(replicationTimeScores_df$score_org_cancer_clines_median)
  
  # variant classification which affect coding sequence
  sim_var_class_with_effect <- c("De_novo_Start_InFrame","De_novo_Start_OutOfFrame", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site", "Start_Codon_Del", "Start_Codon_Ins", "Stop_Codon_Del", "Stop_Codon_Ins")
  ssm_var_class_with_effect <- c("De_novo_Start_InFrame","De_novo_Start_OutOfFrame", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "Start_Codon_SNP")
  

  # columns for variant classification information                                        
  if(mut_type == "ssm")
  {
    all_var_class <- c("3'UTR", "5'Flank","5'UTR","De_novo_Start_InFrame","De_novo_Start_OutOfFrame", "IGR", "Intron","lincRNA","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation", "RNA", "Silent","Splice_Site", "Start_Codon_SNP")

  } else {
    
    all_var_class <- c("3'UTR", "5'Flank","5'UTR","De_novo_Start_InFrame","De_novo_Start_OutOfFrame",  "Frame_Shift_Del","Frame_Shift_Ins", "IGR","Intron","In_Frame_Del", "In_Frame_Ins", "lincRNA","RNA","Splice_Site","Start_Codon_Del","Start_Codon_Ins","Stop_Codon_Del","Stop_Codon_Ins")
  }

  annotation_sampleLevel <-  as.data.frame(matrix(data=0, nrow=nrow(sample2annotation), ncol=(6+length(all_var_class)))) 
  colnames(annotation_sampleLevel) <- c( paste("genes_affected_by_", mut_type, sep=""),   paste("num_", mut_type,"_", sub("'", "",all_var_class),sep=""), paste("num_late_median_cancer_clines_", mut_type, sep=""), paste("perc_late_median_cancer_clines_", mut_type, sep=""), paste("num_early_median_cancer_clines_", mut_type, sep=""), paste("perc_early_median_cancer_clines_", mut_type, sep=""),paste("num_withReplTime_median_cancer_clines_", mut_type, sep=""))
  annotation_sampleLevel[,paste("genes_affected_by_", mut_type, sep="")] <- NA
  
  sample2annotation <- cbind(sample2annotation, annotation_sampleLevel)
  
  
  sample2annotation <- do.call(rbind, mclapply(1:nrow(sample2annotation), function(x)
  {
    
    if(paste(sample2annotation[x, "sample_id"], "_annotatedWithGenCode_ReplTime.txt", sep="") %in% samples_annotatedMutations)
    {
      cur_sample_mutLevelAnn <- read.table(file=paste(annSamplesDir, "/", annSamplesFolder, "/",sample2annotation[x, "sample_id"], "_annotatedWithGenCode_ReplTime.txt",sep=""), quote = "",sep="\t",header = TRUE, stringsAsFactors = FALSE)
      num_mutations <- nrow(cur_sample_mutLevelAnn)
      
      if(num_mutations > 0) 
      {
          # Impact classification
          cur_varClass_count <- table(cur_sample_mutLevelAnn$Variant_Classification)
          names(cur_varClass_count) <- paste("num_", mut_type,"_", sub("'", "",names(cur_varClass_count)),sep="")
          
          sample2annotation[x,names(cur_varClass_count)] <- cur_varClass_count
          
          # Functional category
          sample2annotation[x,paste("num_", mut_type,"_CDS",sep="")] <- nrow(cur_sample_mutLevelAnn[which(cur_sample_mutLevelAnn$CDS),])
          sample2annotation[x,paste("num_", mut_type,"_gene",sep="")] <- nrow(cur_sample_mutLevelAnn[which(cur_sample_mutLevelAnn$gene),])
          sample2annotation[x,paste("num_", mut_type,"_transcript",sep="")] <- nrow(cur_sample_mutLevelAnn[which(cur_sample_mutLevelAnn$transcript),])
          sample2annotation[x,paste("num_", mut_type,"_exon",sep="")] <- nrow(cur_sample_mutLevelAnn[which(cur_sample_mutLevelAnn$exon),])
          sample2annotation[x,paste("num_", mut_type,"_UTR",sep="")] <- nrow(cur_sample_mutLevelAnn[which(cur_sample_mutLevelAnn$UTR),])
          sample2annotation[x,paste("num_", mut_type,"_stop_codon",sep="")] <- nrow(cur_sample_mutLevelAnn[which(cur_sample_mutLevelAnn$stop_codon),])
          sample2annotation[x,paste("num_", mut_type,"_start_codon",sep="")] <- nrow(cur_sample_mutLevelAnn[which(cur_sample_mutLevelAnn$start_codon),])
          sample2annotation[x,paste("num_", mut_type,"_intergenic",sep="")] <- nrow(cur_sample_mutLevelAnn[which(cur_sample_mutLevelAnn$intergenic),])
          sample2annotation[x,paste("num_", mut_type,"_Selenocyste",sep="")] <- nrow(cur_sample_mutLevelAnn[which(cur_sample_mutLevelAnn$Selenocyste),])
      
        if(mut_type == "ssm"){
          
          mut_with_effect <- cur_sample_mutLevelAnn[which(cur_sample_mutLevelAnn$Variant_Classification %in% ssm_var_class_with_effect),]
          
        } else {
          
          mut_with_effect <- cur_sample_mutLevelAnn[which(cur_sample_mutLevelAnn$Variant_Classification %in% sim_var_class_with_effect),]  
        }
        
        if(nrow(mut_with_effect) > 0)
        {
          genes_affected <- getAffectedGenes(mut_with_effect, mut_type, annotation_v19)
          sample2annotation[[x,paste("genes_affected_by_", mut_type, sep="")]] <- list(genes_affected)
          
        }
        
        sample2annotation[x, paste("perc_", mut_type, "_CDS", sep="")] <- (sample2annotation[x, paste("num_", mut_type, "_CDS", sep="")] * 100)/num_mutations
        
        
        # replication time score annotation
        cur_sample_mutLevelAnn$score_org_cancer_clines_median <- rowMedians(as.matrix(cur_sample_mutLevelAnn[, c("score_org_Helas3_1", "score_org_Mcf7_1", "score_org_Sknsh_1", "score_org_Hepg2_1", "score_org_K562_1")]))
        
        sample2annotation[x,paste("num_late_median_cancer_clines_", mut_type,sep="")] <- nrow(cur_sample_mutLevelAnn[which(cur_sample_mutLevelAnn$score_org_cancer_clines_median < median_replTime_median_cancer_clines),])
        
        sample2annotation[x,paste("num_early_median_cancer_clines_", mut_type,sep="")] <- nrow(cur_sample_mutLevelAnn[which(cur_sample_mutLevelAnn$score_org_cancer_clines_median >= median_replTime_median_cancer_clines),])
        
        sample2annotation[x,paste("num_withReplTime_median_cancer_clines_", mut_type,sep="")] <- lateVsEarly_mutation_rate[x,"num_early_median_cancer_clines"] + lateVsEarly_mutation_rate[x,"num_late_median_cancer_clines"]
        
        sample2annotation[x,paste("perc_late_median_cancer_clines_", mut_type,sep="")] <- (lateVsEarly_mutation_rate[x,"num_late_median_cancer_clines"]*100)/lateVsEarly_mutation_rate[x,"num_withReplTime_median_cancer_clines"]
        
        sample2annotation[x,paste("perc_early_median_cancer_clines_", mut_type,sep="")] <- (lateVsEarly_mutation_rate[x,"num_early_median_cancer_clines"]*100)/lateVsEarly_mutation_rate[x,"num_withReplTime_median_cancer_clines"]
      }
      
    } else{
      sample2annotation[x,paste("num_", mut_type,"_CDS",sep="")] <- 0
      sample2annotation[x,paste("num_", mut_type,"_gene",sep="")] <- 0
      sample2annotation[x,paste("num_", mut_type,"_transcript",sep="")] <- 0
      sample2annotation[x,paste("num_", mut_type,"_exon",sep="")] <- 0
      sample2annotation[x,paste("num_", mut_type,"_UTR",sep="")] <- 0
      sample2annotation[x,paste("num_", mut_type,"_stop_codon",sep="")] <- 0
      sample2annotation[x,paste("num_", mut_type,"_start_codon",sep="")] <- 0
      sample2annotation[x,paste("num_", mut_type,"_intergenic",sep="")] <- 0
      sample2annotation[x,paste("num_", mut_type,"_Selenocyste",sep="")] <- 0
    }
    
    
    return(sample2annotation[x,])
    
  }, mc.cores=num_cores))
  
 return(sample2annotation)        
}

#' Combine the smoking status information of PCAWG with the one from TCGA (https://portal.gdc.cancer.gov/legacy-archive/search/f)
#' @param filename_pcawg_tobacco_status: location of the file with the PCAWG information on smoking status
#' @param metadataDir: directory of the PCAWG and TCGA annotation. To link to PCAWG: tcga_donor_uuid <--> bcr_patient_uuid  
#' @return data.frame with the smoking status of PCAWG and TCGA combined
getSmokingStatus <- function(filename_pcawg_tobacco_status, metadataDir) 
{
  ##############################################
  # Read in TCGA information on smoking status #
  ##############################################
  blca_tcga_clinical <- read.table(file = paste(metadataDir, "/nationwidechildrens.org_clinical_patient_blca.txt",sep=""), quote = "", sep = "\t",as.is = TRUE, na.strings = "", check.names = FALSE)
  colnames(blca_tcga_clinical) <- blca_tcga_clinical[2,]
  blca_tcga_clinical <- blca_tcga_clinical[-c(1:3),c("bcr_patient_uuid","bcr_patient_barcode","year_of_initial_pathologic_diagnosis","tobacco_smoking_history","stopped_smoking_year", "number_pack_years_smoked","age_began_smoking_in_years")]
  blca_tcga_clinical$year_of_tobacco_smoking_onset <- NA
  blca_tcga_clinical <- blca_tcga_clinical[,c("bcr_patient_uuid","bcr_patient_barcode","year_of_initial_pathologic_diagnosis","microsatellite_instability", "tobacco_smoking_history", "year_of_tobacco_smoking_onset","stopped_smoking_year", "number_pack_years_smoked","age_began_smoking_in_years")]
  
  cesc_tcga_clinical <- read.table(file = paste(metadataDir, "/nationwidechildrens.org_clinical_patient_cesc.txt",sep=""), quote = "", sep = "\t",as.is = TRUE, na.strings = "", check.names = FALSE)
  colnames(cesc_tcga_clinical) <- cesc_tcga_clinical[2,]
  cesc_tcga_clinical <- cesc_tcga_clinical[-c(1:3),c("bcr_patient_uuid","bcr_patient_barcode","year_of_initial_pathologic_diagnosis","tobacco_smoking_history","stopped_smoking_year", "number_pack_years_smoked","age_began_smoking_in_years")]
  cesc_tcga_clinical$year_of_tobacco_smoking_onset <- NA
  cesc_tcga_clinical <- cesc_tcga_clinical[,c("bcr_patient_uuid","bcr_patient_barcode","year_of_initial_pathologic_diagnosis","microsatellite_instability", "tobacco_smoking_history", "year_of_tobacco_smoking_onset","stopped_smoking_year", "number_pack_years_smoked","age_began_smoking_in_years")]
  
  dlbc_tcga_clinical <- read.table(file = paste(metadataDir, "/nationwidechildrens.org_clinical_patient_dlbc.txt",sep=""), quote = "", sep = "\t",as.is = TRUE, na.strings = "", check.names = FALSE)
  colnames(dlbc_tcga_clinical) <- dlbc_tcga_clinical[2,]
  dlbc_tcga_clinical <- dlbc_tcga_clinical[-c(1:3),c("bcr_patient_uuid","bcr_patient_barcode","year_of_initial_pathologic_diagnosis","tobacco_smoking_history", "stopped_smoking_year", "number_pack_years_smoked","age_began_smoking_in_years")]
  dlbc_tcga_clinical$year_of_tobacco_smoking_onset <- NA
  dlbc_tcga_clinical <- dlbc_tcga_clinical[,c("bcr_patient_uuid","bcr_patient_barcode","year_of_initial_pathologic_diagnosis","microsatellite_instability", "tobacco_smoking_history", "year_of_tobacco_smoking_onset","stopped_smoking_year", "number_pack_years_smoked","age_began_smoking_in_years")]
  
  hnsc_tcga_clinical <- read.table(file = paste(metadataDir, "/nationwidechildrens.org_clinical_patient_hnsc.txt",sep=""), quote = "", sep = "\t",as.is = TRUE, na.strings = "", check.names = FALSE)
  colnames(hnsc_tcga_clinical) <- hnsc_tcga_clinical[2,]
  hnsc_tcga_clinical <- hnsc_tcga_clinical[-c(1:3),c("bcr_patient_uuid","bcr_patient_barcode","year_of_initial_pathologic_diagnosis","tobacco_smoking_history", "year_of_tobacco_smoking_onset","stopped_smoking_year", "number_pack_years_smoked")]
  hnsc_tcga_clinical$age_began_smoking_in_years <- NA
  hnsc_tcga_clinical$age_began_smoking_in_years <- NA
  hnsc_tcga_clinical <- hnsc_tcga_clinical[,c("bcr_patient_uuid","bcr_patient_barcode","year_of_initial_pathologic_diagnosis","microsatellite_instability", "tobacco_smoking_history", "year_of_tobacco_smoking_onset","stopped_smoking_year", "number_pack_years_smoked","age_began_smoking_in_years")]
  
  kich_tcga_clinical <- read.table(file = paste(metadataDir, "/nationwidechildrens.org_clinical_patient_kich.txt",sep=""), quote = "", sep = "\t",as.is = TRUE, na.strings = "", check.names = FALSE)
  colnames(kich_tcga_clinical) <- kich_tcga_clinical[2,]
  kich_tcga_clinical <- kich_tcga_clinical[-c(1:3),c("bcr_patient_uuid","bcr_patient_barcode","year_of_initial_pathologic_diagnosis","tobacco_smoking_history", "year_of_tobacco_smoking_onset","stopped_smoking_year", "number_pack_years_smoked")]
  kich_tcga_clinical$age_began_smoking_in_years <- NA
  kich_tcga_clinical <- kich_tcga_clinical[,c("bcr_patient_uuid","bcr_patient_barcode","year_of_initial_pathologic_diagnosis","microsatellite_instability", "tobacco_smoking_history", "year_of_tobacco_smoking_onset","stopped_smoking_year", "number_pack_years_smoked","age_began_smoking_in_years")]
  
  kirc_tcga_clinical <- read.table(file = paste(metadataDir, "/nationwidechildrens.org_clinical_patient_kirc.txt",sep=""), quote = "", sep = "\t",as.is = TRUE, na.strings = "", check.names = FALSE)
  colnames(kirc_tcga_clinical) <- kirc_tcga_clinical[2,]
  kirc_tcga_clinical <- kirc_tcga_clinical[-c(1:3),c("bcr_patient_uuid","bcr_patient_barcode","year_of_initial_pathologic_diagnosis","tobacco_smoking_history", "year_of_tobacco_smoking_onset","stopped_smoking_year", "number_pack_years_smoked")]
  kirc_tcga_clinical$age_began_smoking_in_years <- NA
  kirc_tcga_clinical <- kirc_tcga_clinical[,c("bcr_patient_uuid","bcr_patient_barcode","year_of_initial_pathologic_diagnosis","microsatellite_instability", "tobacco_smoking_history", "year_of_tobacco_smoking_onset","stopped_smoking_year", "number_pack_years_smoked","age_began_smoking_in_years")]
  
  kirp_tcga_clinical <- read.table(file = paste(metadataDir, "/nationwidechildrens.org_clinical_patient_kirp.txt",sep=""), quote = "", sep = "\t",as.is = TRUE, na.strings = "", check.names = FALSE)
  colnames(kirp_tcga_clinical) <- kirp_tcga_clinical[2,]
  kirp_tcga_clinical <- kirp_tcga_clinical[-c(1:3),c("bcr_patient_uuid","bcr_patient_barcode","year_of_initial_pathologic_diagnosis","tobacco_smoking_history", "year_of_tobacco_smoking_onset","stopped_smoking_year", "number_pack_years_smoked")]
  kirp_tcga_clinical$age_began_smoking_in_years <- NA
  kirp_tcga_clinical <- kirp_tcga_clinical[,c("bcr_patient_uuid","bcr_patient_barcode","year_of_initial_pathologic_diagnosis","microsatellite_instability", "tobacco_smoking_history", "year_of_tobacco_smoking_onset","stopped_smoking_year", "number_pack_years_smoked","age_began_smoking_in_years")]
  
  luad_tcga_clinical <- read.table(file = paste(metadataDir, "/nationwidechildrens.org_clinical_patient_luad.txt",sep=""), quote = "", sep = "\t",as.is = TRUE, na.strings = "", check.names = FALSE)
  colnames(luad_tcga_clinical) <- luad_tcga_clinical[2,]
  luad_tcga_clinical <- luad_tcga_clinical[-c(1:3),c("bcr_patient_uuid","bcr_patient_barcode","year_of_initial_pathologic_diagnosis","tobacco_smoking_history", "year_of_tobacco_smoking_onset","stopped_smoking_year", "number_pack_years_smoked","age_began_smoking_in_years")]
  luad_tcga_clinical <- luad_tcga_clinical[,c("bcr_patient_uuid","bcr_patient_barcode","year_of_initial_pathologic_diagnosis","microsatellite_instability", "tobacco_smoking_history", "year_of_tobacco_smoking_onset","stopped_smoking_year", "number_pack_years_smoked","age_began_smoking_in_years")]
  
  lusc_tcga_clinical <- read.table(file = paste(metadataDir, "/nationwidechildrens.org_clinical_patient_lusc.txt",sep=""), quote = "", sep = "\t",as.is = TRUE, na.strings = "", check.names = FALSE)
  colnames(lusc_tcga_clinical) <- lusc_tcga_clinical[2,]
  lusc_tcga_clinical <- lusc_tcga_clinical[-c(1:3),c("bcr_patient_uuid","bcr_patient_barcode","year_of_initial_pathologic_diagnosis","tobacco_smoking_history", "year_of_tobacco_smoking_onset","stopped_smoking_year", "number_pack_years_smoked","age_began_smoking_in_years")]
  lusc_tcga_clinical <- lusc_tcga_clinical[,c("bcr_patient_uuid","bcr_patient_barcode","year_of_initial_pathologic_diagnosis","microsatellite_instability", "tobacco_smoking_history", "year_of_tobacco_smoking_onset","stopped_smoking_year", "number_pack_years_smoked","age_began_smoking_in_years")]
  
  tcga_tobacco_status <- rbind(blca_tcga_clinical,cesc_tcga_clinical,  dlbc_tcga_clinical, kich_tcga_clinical, kirc_tcga_clinical, kirp_tcga_clinical, luad_tcga_clinical, lusc_tcga_clinical,hnsc_tcga_clinical)
  
  #####################################
  # read in smoking status from PCAWG #
  #####################################
  pcawg_tobacco_status <- read.table(loc_file_pcawg_tobacco_status, sep="\t", quote = "", na.strings = "",as.is = TRUE, check.names = FALSE, header=TRUE, fill=TRUE)
  
  ##################################
  # merge info from TCGA and PCAWG #
  ##################################
  donor2smoking_status <- merge(pcawg_tobacco_status, tcga_tobacco_status, by.x="tcga_donor_uuid", by.y="bcr_patient_uuid", all.x=TRUE)
  
  ###############################
  # convert to same terminology #
  ###############################
  for(i in 1:nrow(donor2smoking_status))
  {
    if(!(is.na(donor2smoking_status[i,"tobacco_smoking_history"])) & (donor2smoking_status[i,"tobacco_smoking_history"] != "[Not Available]") & (donor2smoking_status[i,"tobacco_smoking_history"] != "[Unknown]"))
    {
      if(donor2smoking_status[i,"tobacco_smoking_history"] == 1)
        donor2smoking_status[i,"smoking_status"] <- "lifelong_non_smoker" 
      else if(donor2smoking_status[i,"tobacco_smoking_history"] == 2)
        donor2smoking_status[i,"smoking_status"] <- "cur_smoker" 
      else if(donor2smoking_status[i,"tobacco_smoking_history"] == 3)
        donor2smoking_status[i,"smoking_status"] <- "cur_ref_smoker_over15y" 
      else if(donor2smoking_status[i,"tobacco_smoking_history"] == 4)
        donor2smoking_status[i,"smoking_status"] <- "cur_ref_smoker_less15y"
      else if(donor2smoking_status[i,"tobacco_smoking_history"] == 5)
        donor2smoking_status[i,"smoking_status"] <- "cur_ref_smoker_unknown"
      else if(donor2smoking_status[i,"tobacco_smoking_history"] == 6)
        donor2smoking_status[i,"smoking_status"] <- "smoker_at_diagnosis"
      else if(donor2smoking_status[i,"tobacco_smoking_history"] == 7)
        donor2smoking_status[i,"smoking_status"] <- "smoking_hist_not_docu" 
    }
    else
    {
      if(!(is.na(donor2smoking_status[i,"tobacco_smoking_history_indicator"])))
      {
        if(donor2smoking_status[i,"tobacco_smoking_history_indicator"] == "Lifelong non-smoker (<100 cigarettes smoked in lifetime)")
          donor2smoking_status[i,"smoking_status"] <- "lifelong_non_smoker" 
        else if(donor2smoking_status[i,"tobacco_smoking_history_indicator"] == "Current smoker (includes daily smokers non-daily/occasional smokers)")
          donor2smoking_status[i,"smoking_status"] <- "cur_smoker" 
        else if(donor2smoking_status[i,"tobacco_smoking_history_indicator"] == "Current reformed smoker, duration not specified")
          donor2smoking_status[i,"smoking_status"] <- "cur_ref_smoker_unknown"
        else if(donor2smoking_status[i,"tobacco_smoking_history_indicator"] == "Current reformed smoker for <= 15 years")
          donor2smoking_status[i,"smoking_status"] <- "cur_ref_smoker_less15y"
        else if(donor2smoking_status[i,"tobacco_smoking_history_indicator"] == "Current reformed smoker for > 15 years")
          donor2smoking_status[i,"smoking_status"] <- "cur_ref_smoker_over15y"
        else if(donor2smoking_status[i,"tobacco_smoking_history_indicator"] == "Smoking history not documented")
          donor2smoking_status[i,"smoking_status"] <- "unknown"
      }
      else
      {
        donor2smoking_status[i,"smoking_status"] <- "unknown"
      }
      
    }
  }
  
  return(donor2smoking_status)
}

# Metadata: ID in PCAWG <--> ID in metadata file, location data
#'  @param sample2ttype_cluster: sample linked to tumor type, cluster and identifiers needed to link to annotation
#'  @param annotation_v19: GENCODE information (https://www.gencodegenes.org/human/releases.html)
#'  @param replicationTimeScores_df: data.frame with replication time scores
#'  @param metadataDir: directory in which the metadata is stored
#'  @param filename_drivers: filename of the file with predicted driver mutations. To link to PCAWG: tumor_wgs_aliquot_id <--> sample_id, https://www.biorxiv.org/content/10.1101/190330v2
#'  @param filename_IGHV_status: filename of the file with IGHV data. To link to PCAWG: submitted_donor_id <--> Case,  https://www.nature.com/articles/nature14666
#'  @param filename_MSI_classification_1: filename of the file with the MSI status according to MSI-Method 1. To link to PCAWG: tumor_wgs_submitter_sample_id <--> ID, https://docs.icgc.org/pcawg/data/
#'  @param filename_MSI_classification_2: filename of the file with the MSI status according to MSI-Method 2. To link to PCAWG:tumor_wgs_aliquot_id <--> tumor_wgs_aliquot_id, https://www.biorxiv.org/content/10.1101/208330v1
#'  @param filename_SBS_signatures: filename of the file with the SBS signatures contribution. To link to PCAWG: tumor_wgs_icgc_specimen_id <--> icgc_specimen_id,  https://www.biorxiv.org/content/10.1101/322859v2
#'  @param filename_DBS_signatures: filename of the file with the DBS signatures contribution. To link to PCAWG: tumor_wgs_icgc_specimen_id <--> icgc_specimen_id,  https://www.biorxiv.org/content/10.1101/322859v2
#'  @param filename_ID_signatures: filename of the file with the ID signatures contribution. To link to PCAWG: tumor_wgs_icgc_specimen_id <--> icgc_specimen_id,  https://www.biorxiv.org/content/10.1101/322859v2
#'  @param filename_pcawg_smoking_status: filename of the file with smoking status from PCAWG. Annotation PCAWG uses the icgc_donor_id, https://docs.icgc.org/pcawg/data/
#'  @param annSamplesDir: directory in which the annotated sample files on mutation levels are stored
#'  @param annSamplesFolder_ssms: folder of the annotated sample files on mutation level for SSMs
#'  @param annSamplesFolder_sims: folder of the annotated sample files on mutation level for SIMs  
#'  @param num_cores: number of cores to be used in mclapply              
#' @return data.frame with samples annotated with the tumour type, cluster assigned to, MSI status, IGHV status (Lymph-CLL only), tobacco smoking status (when available), contribution of each signature (SBS, DBS and ID) and the summary of the impact classification and functional category of its mutations. 
annotateAtSampleLevel <- function(sample2ttype_cluster, annotation_v19, replicationTimeScores_df, metadataDir, filename_drivers, filename_IGHV_status, filename_MSI_classification_1, filename_MSI_classification_2, filename_SBS_signatures, filename_DBS_signatures, filename_ID_signatures, filename_pcawg_smoking_status, annSamplesDir, annSamplesFolder_ssms, annSamplesFolder_sims, num_cores)
{
  
    ###################################################################
    # Impact classification, functional category and replication time #
    ###################################################################
  
    # summarize Impact classification, functional category and replication time annotation to sample level
    sample2annotation <- summarizeMutAnnotation2SampleLevel(sample2ttype_cluster, "ssm", annotation_v19, replicationTimeScores_df, annSamplesDir, annSamplesFolder_ssms, num_cores)
    sample2annotation <- summarizeMutAnnotation2SampleLevel(sample2annotation, "sim", annotation_v19, replicationTimeScores_df, annSamplesDir, annSamplesFolder_sims, num_cores)

    # combine the list of genes with a coding change due to a SIM and the list for SSMs
    sample2annotation$genes_affected_by_mut <- NA
    sample2annotation$genes_affected_by_mut <- lapply(1:nrow(sample2annotation), function(x) 
    {
      genes_affected_by_ssms <- unlist(sample2annotation[x,"genes_affected_by_ssm"])
      genes_affected_by_sims <- unlist(sample2annotation[x,"genes_affected_by_sim"])
      
      if(!any(is.na(genes_affected_by_sims)) & !any(is.na(genes_affected_by_ssms))) {
        affected_genes_all <- list(unique(c(genes_affected_by_ssms, genes_affected_by_sims)))
      } else if(!any(is.na(genes_affected_by_ssms))) {
        affected_genes_all <- list(genes_affected_by_ssms)
      } else if(!any(is.na(genes_affected_by_sims))) {
        affected_genes_all <- list(genes_affected_by_sims)
      } else{
        affected_genes_all <- NA
      }
      
      return(affected_genes_all)
    })
    
    # count the number of genes with a coding change due to either a SIM or SSM per sample
    sample2annotation$num_genes_affected_by_mut <- unlist(lapply(1:nrow(sample2annotation), function(x) 
    {
      cur_genes_affected_by_mut <- unlist(sample2annotation[x,"genes_affected_by_mut"])
      
      if(!(all(is.na(cur_genes_affected_by_mut))))
        num_genes_affected_by_mut <- length(cur_genes_affected_by_mut)
      else 
        num_genes_affected_by_mut <- 0
      
      return(num_genes_affected_by_mut)
    }
    ))
    
    # count the number of genes with a coding change due to a SSM per sample
    sample2annotation$num_genes_affected_by_ssm <- unlist(lapply(1:nrow(sample2annotation), function(x) 
    {
      cur_genes_affected_by_ssms <- unlist(sample2annotation[x,"genes_affected_by_ssm"])
      
      if(!(all(is.na(cur_genes_affected_by_ssms))))
        num_genes_affected_by_ssm <- length(cur_genes_affected_by_ssms)
      else 
        num_genes_affected_by_ssm <- 0
      
      return(num_genes_affected_by_ssm)
    }
    ))
    
    # count the number of genes with a coding change due to either a SIM per sample
    sample2annotation$num_genes_affected_by_sim <- unlist(lapply(1:nrow(sample2annotation), function(x) 
    {
      cur_genes_affected_by_sims <- unlist(sample2annotation[x,"genes_affected_by_sim"])
      
      if(!(all(is.na(cur_genes_affected_by_sims))))
        num_genes_affected_by_sim <- length(cur_genes_affected_by_sims)
      else 
        num_genes_affected_by_sim <- 0
      
      return(num_genes_affected_by_sim)
    }
    ))
    
    # number of SSMs in a gene
    sample2annotation$num_in_coding_ssms <- sample2annotation$num_ssm_De_novo_Start_InFrame + sample2annotation$num_ssm_De_novo_Start_OutOfFrame + sample2annotation$num_ssm_Missense_Mutation + sample2annotation$num_ssm_Nonsense_Mutation + sample2annotation$num_ssm_Nonstop_Mutation + sample2annotation$num_ssm_Splice_Site + sample2annotation$num_ssm_Start_Codon_SNP + sample2annotation$num_ssm_Silent
    
    # percentage of silent SSMs of the total number of SSMs in a gene
    sample2annotation$perc_silent_ssm_vs_in_coding <- (sample2annotation$num_ssm_Silent*100)/sample2annotation$num_in_coding_ssms
    sample2annotation[which(is.na(sample2annotation$perc_silent_ssm_vs_in_coding)), "perc_silent_ssm_vs_in_coding"] <- 0
    
    
    #####################
    # Predicted drivers #
    #####################
    
    drivers_per_donor <- read.table(paste(metadataDir, "/", filename_drivers, sep=""), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    
    sample2annotation$drivers <- NA
    
    for(i in 1:nrow(sample2annotation)){
       cur_drivers <-  list(unique(drivers_per_donor[which(drivers_per_donor$sample_id == sample2annotation[i,"tumor_wgs_aliquot_id"]),"gene"]))
       sample2annotation[[i,"drivers"]] <- cur_drivers
       sample2annotation[i,"num_drivers"] <- length(cur_drivers)
    }
    
    ##############
    # Signatures #
    ##############
    
    # convert absolute contribution of each SBS signature to perentage and add both to annotation
    sample2signContr_SBS <- read.table(paste(metadataDir, "/", filename_file_SBS_signatures,sep=""), sep=",", quote = "", header = TRUE)
    colnames(sample2signContr_SBS)[1:3] <- c("tumor_type", "icgc_specimen_id", "Accuracy_SBS")
    sign_SBS_cols <- grep("^SBS", colnames(sample2signContr_SBS),value=TRUE)
    sample2signContr_SBS[,paste("perc_", sign_SBS_cols,sep="")] <- (sample2signContr_SBS[,sign_SBS_cols]*100)/rowSums(sample2signContr_SBS[,sign_SBS_cols])
    sample2annotation <- merge(sample2signContr_SBS,sample2annotation, by.x="icgc_specimen_id", by.y="tumor_wgs_icgc_specimen_id")
    
    # convert absolute contribution of each DBS signature to perentage and add both to annotation
    sample2signContr_DBS <- read.table(paste(metadataDir, "/", filename_DBS_signatures,sep=""), sep=",", quote = "", header = TRUE)
    colnames(sample2signContr_DBS)[1:3] <- c("tumor_type", "icgc_specimen_id", "Accuracy_DBS")
    sign_DBS_cols <- grep("^DBS", colnames(sample2signContr_DBS), value=TRUE)
    sample2signContr_DBS[,paste("perc_", sign_DBS_cols,sep="")] <- (sample2signContr_DBS[,sign_DBS_cols]*100)/rowSums(sample2signContr_DBS[,sign_DBS_cols])
    sample2annotation <- merge(sample2signContr_DBS,sample2annotation, by.x="icgc_specimen_id", by.y="tumor_wgs_icgc_specimen_id")
    
    # convert absolute contribution of each ID signature to perentage and add both to annotation
    sample2signContr_ID <- read.table(paste(metadataDir, "/", filename_ID_signatures,sep=""), sep=",", quote = "", header = TRUE)
    colnames(sample2signContr_ID)[1:3] <- c("tumor_type", "icgc_specimen_id", "Accuracy_ID")
    sign_ID_cols <- grep("^ID", colnames(sample2signContr_ID), value=TRUE)
    sample2signContr_ID[,paste("perc_", sign_ID_cols,sep="")] <- (sample2signContr_ID[,sign_ID_cols]*100)/rowSums(sample2signContr_ID[,sign_ID_cols])
    sample2annotation <- merge(sample2signContr_ID,sample2annotation, by.x="icgc_specimen_id", by.y="tumor_wgs_icgc_specimen_id")
    
    # count per sample the number of signatures (SBS, DBS, ID) with a non-zero contribution
    for(i in 1:nrow(sample2annotation))
    {
      cur_sign_SBS_perc_contr <- sample2annotation[i, grep("^perc_SBS", colnames(sample2annotation), value=TRUE)]
      sample2annotation[i,"num_SBS_sign_nonZeroContr"] <- length(cur_sign_SBS_perc_contr[which(cur_sign_SBS_perc_contr >0)])
      
      cur_sign_DBS_perc_contr <- sample2annotation[i, grep("^perc_DBS", colnames(sample2annotation), value=TRUE)]
      sample2annotation[i,"num_DBS_sign_nonZeroContr"] <- length(cur_sign_DBS_perc_contr[which(cur_sign_DBS_perc_contr >0)])
      
      cur_sign_ID_perc_contr <- sample2annotation[i, grep("^perc_ID", colnames(sample2annotation), value=TRUE)]
      sample2annotation[i,"num_ID_sign_nonZeroContr"] <- length(cur_sign_ID_perc_contr[which(cur_sign_ID_perc_contr >0)])
    }

    # total number of signatures with a non-zero contribution
    sample2annotation[,"num_all_sign_nonZeroContr"] <- sample2annotation[,"num_SBS_sign_nonZeroContr"] + sample2annotation[,"num_DBS_sign_nonZeroContr"] + sample2annotation[,"num_ID_sign_nonZeroContr"]
    
    
    ###############
    # IGHV status #
    ###############
    
    #IGHV status is only available for Lymph-CLL samples
    ighv_info <- read.table(paste(metadataDir, "/", filename_IGHV_status,sep=""), sep="\t", quote = "", na.strings = "",as.is = TRUE, check.names = FALSE, header=TRUE, fill=TRUE)
    sample2annotation <- merge(sample2annotation, ighv_info, by.x="submitted_donor_id", by.y="Case", all.x=TRUE)
    
    
    ##################
    # Smoking status #
    ##################
    
    # Add to the smoking status information provided by PCAWG the information provided by TCGA
    donor2smokingStatus <- getSmokingStatus(filename_pcawg_tobacco_status, metadataDir)  
    sample2annotation <- merge(sample2annotation, donor2smokingStatus, by="icgc_donor_id", all.x=TRUE)
    
    
    ##############
    # MSI status #
    ##############
    
    
    # MSI-Method 1
    MSI_classification_1 <- read.table(paste(metadataDir, "/", filename_MSI_classification_1, sep=""), sep="\t", quote = "", na.strings = "",as.is = TRUE, check.names = FALSE, header=TRUE, fill=TRUE)
    sample2annotation <- merge(sample2annotation, MSI_classification_1, by.x="tumor_wgs_submitter_sample_id", by.y="ID", all.x=TRUE)
    
    # MSI-Method 2
    MSI_classification_2 <- read.table(paste(metadataDir, "/", filename_MSI_classification_2, sep=""), sep="\t", quote = "", na.strings = "",as.is = TRUE, check.names = FALSE, header=TRUE, fill=TRUE)
    sample2annotation <- merge(sample2annotation, MSI_classification_2, by="tumor_wgs_aliquot_id", all.x=TRUE)
 
       
    ###############################
    # set missing entries to zero #
    ###############################
    
    colsWithMissingData <- apply(sample2annotation, 2, function(x){ any(is.na(x))})
    colsWithMissingData <- names(colsWithMissingData[which(colsWithMissingData ==TRUE)])
    
    # for genes_affected_by_ssm, genes_affected_by_sim and genes_affected_by_mut if there are no values it will be NA instead 0. 
    colsWithMissingData <- colsWithMissingData[which(!colsWithMissingData %in% c("genes_affected_by_ssm","genes_affected_by_sim","genes_affected_by_mut"))]
    
    for(i in 1:length(colsWithMissingData))
    {
      sample2annotation[which(is.na(sample2annotation[,colsWithMissingData[i]])), colsWithMissingData[i]] <- 0
    }
    
    return(sample2annotation)
}