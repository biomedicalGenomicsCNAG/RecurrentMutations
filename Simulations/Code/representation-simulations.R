# Libraries 
library(ggplot2)
library(reshape2)
library(scales)
source("pcawg.colour.palette.R")
list_tumortypes <- c("Biliary-AdenoCA","Bladder-TCC","Bone-Benign","Bone-Epith","Bone-Osteosarc","Breast-AdenoCA","Breast-DCIS","Breast-LobularCA",
                     "Cervix-AdenoCA","Cervix-SCC","CNS-GBM","CNS-Medullo","CNS-Oligo","CNS-PiloAstro","ColoRect-AdenoCA","Eso-AdenoCA","Head-SCC","Kidney-ChRCC","Kidney-RCC",
                     "Liver-HCC","Lung-AdenoCA","Lung-SCC","Lymph-BNHL","Lymph-CLL","Myeloid-AML","Myeloid-MDS","Myeloid-MPN","Ovary-AdenoCA","Panc-AdenoCA","Panc-Endocrine",
                     "Prost-AdenoCA","Skin-Melanoma","SoftTissue-Leiomyo","SoftTissue-Liposarc","Stomach-AdenoCA","Thy-AdenoCA","Uterus-AdenoCA")


#########################################################################################################
#### 1. RECURRENCE WITHIN TUMOR TYPE AND ACROSS THE WHOLE PCAWG COHORT: SUMMARY OF 5000 SIMULATIONS  ####
#########################################################################################################

# Simulations data
load("recurrence_tumortype_pancancer_5000simul.RData")
all_recurrence_tt_pan <- recurrence_tumortype_pancancer
load("recurrence_tumortype_within_5000simul.RData")
all_recurrence_tt_wt <- recurrence_tumortype_within

# Observed data in the cohort
load("obs_ssms_numbers_per_tumor.RData")

# Melting simulated data
  # PAN
all_recurrence_tt_pan_melted <- melt(all_recurrence_tt_pan)
  # WT
for (column in 2:38){
  all_recurrence_tt_wt[,column] <- as.numeric(as.character(all_recurrence_tt_wt[,column]))
}
all_recurrence_tt_wt_melted <- melt(all_recurrence_tt_wt)

# Computing the ratio (how many times the observed value is higher that the simulated)
all_recurrence_tt_pan_melted_times <- all_recurrence_tt_pan_melted
all_recurrence_tt_pan_melted_times$real_value <- rep(0,37000)
all_recurrence_tt_wt_melted_times <- all_recurrence_tt_wt_melted
all_recurrence_tt_wt_melted_times$real_value <- rep(0,37000)

for(tumor in list_tumortypes){
  real_value_recurrence_pan <- ssms_numbers_per_tumor$num_rec_pan[which(ssms_numbers_per_tumor$tumor_type == tumor)]
  all_recurrence_tt_pan_melted_times$real_value[which(all_recurrence_tt_pan_melted_times$variable==tumor)] <- real_value_recurrence_pan
  
  real_value_recurrence_wt <- ssms_numbers_per_tumor$num_rec_wt[which(ssms_numbers_per_tumor$tumor_type == tumor)]
  all_recurrence_tt_wt_melted_times$real_value[which(all_recurrence_tt_wt_melted_times$variable==tumor)] <- real_value_recurrence_wt
}

# PAN
all_recurrence_tt_pan_melted_times$times_Real_over_Simul <- all_recurrence_tt_pan_melted_times$real_value/all_recurrence_tt_pan_melted_times$value

# WT
all_recurrence_tt_wt_melted_times$real_value <- as.numeric(all_recurrence_tt_wt_melted_times$real_value)
all_recurrence_tt_wt_melted_times$times_Real_over_Simul <- all_recurrence_tt_wt_melted_times$real_value/all_recurrence_tt_wt_melted_times$value

# [STXT1-1D] BOXPLOT Recurrence pancancer
  # colors
pcawg_colors <- pcawg.colour.palette(as.character(unique(all_recurrence_tt_pan_melted_times$variable)),"tumour.subtype")
names(pcawg_colors) <- as.character(unique(all_recurrence_tt_pan_melted_times$variable))

  # plot
boxplot_recurrence_pan <- ggplot(aes(x=variable,y=times_Real_over_Simul), data=all_recurrence_tt_pan_melted_times) + 
  geom_boxplot(color="black", fill=pcawg_colors) + 
  ylab("ratio of observed vs. simulated number of recurrent SSMs") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.4, size = 8, hjust=1),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12)) +
  scale_y_continuous(
    breaks = seq(floor(min(all_recurrence_tt_pan_melted_times$times_Real_over_Simul)), max(all_recurrence_tt_pan_melted_times$times_Real_over_Simul), by=1),
    limits = c(floor(min(all_recurrence_tt_pan_melted_times$times_Real_over_Simul)), max(all_recurrence_tt_pan_melted_times$times_Real_over_Simul)),
    labels = seq(floor(min(all_recurrence_tt_pan_melted_times$times_Real_over_Simul)), max(all_recurrence_tt_pan_melted_times$times_Real_over_Simul), by=1))

  # plot zooming in
boxplot_recurrence_pan <- ggplot(aes(x=variable,y=times_Real_over_Simul), data=all_recurrence_tt_pan_melted_times) + 
  geom_boxplot(color="black", fill=pcawg_colors) + 
  ylab("ratio of observed vs. simulated number of recurrent SSMs") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.4, size = 8, hjust=1),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12)) +
  scale_y_continuous(
    breaks = seq(floor(min(all_recurrence_tt_pan_melted_times$times_Real_over_Simul)), max(all_recurrence_tt_pan_melted_times$times_Real_over_Simul), by=1),
    limits = c(floor(min(all_recurrence_tt_pan_melted_times$times_Real_over_Simul)), max(all_recurrence_tt_pan_melted_times$times_Real_over_Simul)),
    labels = seq(floor(min(all_recurrence_tt_pan_melted_times$times_Real_over_Simul)), max(all_recurrence_tt_pan_melted_times$times_Real_over_Simul), by=1)) +
  coord_cartesian(xlim = NULL, ylim = c(1,8), expand = TRUE)

# [STXT1-1C] BOXPLOT Recurrence within tumor type
  # Filter "red": avoid this tumor types because for more than 40% of the simulations the recurrence is 0
tumortypes_to_avoid <- c("Bone-Benign", "Bone-Epith", "Breast-DCIS", "Breast-LobularCA",
                         "Breast-LobularCA", "Cervix-AdenoCA", "CNS-Oligo", "CNS-PiloAstro", "Kidney-ChRCC",
                         "Myeloid-AML", "Myeloid-MDS", "Myeloid-MPN", "SoftTissue-Leiomyo", "SoftTissue-Liposarc", 
                         "Thy-AdenoCA")

all_recurrence_tt_wt_melted_times_filter_red <- all_recurrence_tt_wt_melted_times[which(! all_recurrence_tt_wt_melted_times$variable %in% tumortypes_to_avoid),]
dim(all_recurrence_tt_wt_melted_times_filter_red) #115000      5

  # Filter "green": leave out "inf" cases (recurrent simulation = 0), not all tumortype because these cases were less than 40%
all_recurrence_tt_wt_melted_times_filter_red_green <- all_recurrence_tt_wt_melted_times_filter_red[which(! all_recurrence_tt_wt_melted_times_filter_red$times_Real_over_Simul == "Inf"),]
dim(all_recurrence_tt_wt_melted_times_filter_red_green) # 112820      5

  # colors
pcawg_colors <- pcawg.colour.palette(as.character(unique(all_recurrence_tt_wt_melted_times_filter_red_green$variable)),"tumour.subtype")
names(pcawg_colors) <- as.character(unique(all_recurrence_tt_wt_melted_times_filter_red_green$variable))

  # plot
boxplot_recurrence_wt <- ggplot(aes(x=variable,y=times_Real_over_Simul), data=all_recurrence_tt_wt_melted_times_filter_red_green) + 
  geom_boxplot(color="black", fill=pcawg_colors) + 
  ylab("ratio of observed vs. simulated number of recurrent SSMs") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.4, size = 8, hjust=1),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12)) +
  scale_y_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x, n=10),
                     limits = c(floor(min(all_recurrence_tt_wt_melted_times_filter_red_green$times_Real_over_Simul)), ceiling(max(all_recurrence_tt_wt_melted_times_filter_red_green$times_Real_over_Simul))),
                     labels = trans_format("log2", math_format(2^.x))) +
  annotation_logticks(base=2, sides="l")


############################################################
####  2. OVERALL RECURRENCE: SUMMARY 5000 SIMULATIONS   ####
############################################################

# Simulations data
load("recurrence_overall_5000simul.RData")
overall_dataframe_recurrence <- recurrence_overall

# Computing ratio
overall_dataframe_recurrence$obs_recurrence <- rep(1057935,5000)
overall_dataframe_recurrence$times_Real_over_Simul <- overall_dataframe_recurrence$obs_recurrence/overall_dataframe_recurrence$num_recurrent

# [STXT1-1A] Plot
overall_dataframe_recurrence$type <- as.factor(rep("Overall",5000))

boxplot_overall <- ggplot(aes(x=type,y=times_Real_over_Simul), data=overall_dataframe_recurrence) + 
  geom_boxplot(color="black", fill="#F8766D") + 
  ylab("ratio of observed vs. simulated number of recurrent SSMs") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size=12)) +
  scale_x_discrete(labels ="Overall") +
  scale_y_continuous( 
    breaks=seq(4.98,5.08,by=0.02),
    labels=format(round(seq(4.98,5.08,by=0.02),digits=2),nsmal=2),
    limits=c(4.98,5.08))


#####################################################################
####  3. RECURRENCE PER MUTATION TYPE: SUMMARY 5000 SIMULATIONS  ####
#####################################################################

# Simulations data
load("recurrence_subtype_5000simul.RData")

# Observed values of recurrence
median_rec_observed <- c(141747, 8856, 758411,25167, 38336, 85418)

# Adding observed values to the dataframe
c_a_obs <- rep(median_rec_observed[1], 5000)
c_g_obs <- rep(median_rec_observed[2], 5000)
c_t_obs <- rep(median_rec_observed[3], 5000)
t_a_obs <- rep(median_rec_observed[4], 5000)
t_c_obs <- rep(median_rec_observed[5], 5000)
t_g_obs <- rep(median_rec_observed[6], 5000)

all_recurrence_subtype_simul_obs <- as.data.frame(cbind(all_recurrence_subtype, c_a_obs, c_g_obs, c_t_obs, t_a_obs, t_c_obs, t_g_obs))
all_recurrence_subtype_simul_obs_ratio <- all_recurrence_subtype_simul_obs
all_recurrence_subtype_simul_obs_ratio$ratio_c_a <- all_recurrence_subtype_simul_obs_ratio$c_a_obs/all_recurrence_subtype_simul_obs_ratio$c_a
all_recurrence_subtype_simul_obs_ratio$ratio_c_g <- all_recurrence_subtype_simul_obs_ratio$c_g_obs/all_recurrence_subtype_simul_obs_ratio$c_g
all_recurrence_subtype_simul_obs_ratio$ratio_c_t <- all_recurrence_subtype_simul_obs_ratio$c_t_obs/all_recurrence_subtype_simul_obs_ratio$c_t
all_recurrence_subtype_simul_obs_ratio$ratio_t_a <- all_recurrence_subtype_simul_obs_ratio$t_a_obs/all_recurrence_subtype_simul_obs_ratio$t_a
all_recurrence_subtype_simul_obs_ratio$ratio_t_c <- all_recurrence_subtype_simul_obs_ratio$t_c_obs/all_recurrence_subtype_simul_obs_ratio$t_c
all_recurrence_subtype_simul_obs_ratio$ratio_t_g <- all_recurrence_subtype_simul_obs_ratio$t_g_obs/all_recurrence_subtype_simul_obs_ratio$t_g

all_recurrence_subtype_ratio_subset <- all_recurrence_subtype_simul_obs_ratio[,c(1,14:length(all_recurrence_subtype_simul_obs_ratio))]
all_recurrence_subtype_ratio_melted <- melt(all_recurrence_subtype_ratio_subset)

# [STXT1-1B] Plot

colors_sele <- c("blue", "black", "red", "pink", "yellow", "green")
subtypes <- c("C>A", "C>G", "C>T", "T>A"," T>C", "T>G")

boxplot_recurrence_substitutions <- ggplot(aes(x=variable,y=value), data=all_recurrence_subtype_ratio_melted) +  geom_boxplot(color="black", fill=colors_sele) + 
  ylab("ratio of observed vs. simulated number of recurrent SSMs") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust=0.4, size = 8, hjust=0.5),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12)) +
  scale_x_discrete(labels= subtypes) +
  scale_y_continuous(
    breaks = seq(1, 13, by=1),
    limits = c(min(1,floor(min(all_recurrence_subtype_ratio_melted$value))), ceiling(max(all_recurrence_subtype_ratio_melted$value))))

