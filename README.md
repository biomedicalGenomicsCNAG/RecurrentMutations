# Recurrent somatic mutations reveal new insights into consequences of mutagenic processes in cancer

<h4>Abstract</h4>
<p align="justify">The sheer size of the human genome makes it improbable that identical somatic mutations at the exact same position are observed in multiple tumours solely by chance. The scarcity of cancer driver mutations also precludes positive selection as the sole explanation. Therefore, recurrent mutations may be highly informative of characteristics of mutational processes. To explore the potential, we use recurrence as a starting point to cluster >2,500 whole genomes of a pan-cancer cohort. We describe each genome with 13 recurrence-based and 29 general mutational features. Using principal component analysis we reduce the dimensionality and create independent features. We apply hierarchical clustering to the first 18 principal components followed by k-means clustering. We show that the resulting 16 clusters capture clinically relevant cancer phenotypes. High levels of recurrent substitutions separate the clusters that we link to UV-light exposure and deregulated activity of POLE from the one representing defective mismatch repair, which shows high levels of recurrent insertions/deletions. Recurrence of both mutation types characterizes cancer genomes with somatic hypermutation of immunoglobulin genes and the cluster of genomes exposed to gastric acid. Low levels of recurrence are observed for the cluster where tobacco-smoke exposure induces mutagenesis and the one linked to increased activity of cytidine deaminases. Notably, the majority of substitutions are recurrent in a single tumour type, while recurrent insertions/deletions point to shared processes between tumour types. Recurrence also reveals susceptible sequence motifs, including TT[C>A]TTT and AAC[T>G]T for the POLE and ‘gastric acid-exposure’ clusters, respectively. Moreover, we refine knowledge of mutagenesis, including increased C/G deletion levels in general for lung tumours and specifically in midsize homopolymer sequence contexts for microsatellite instable tumours. Our findings are an important step towards the development of a generic cancer diagnostic test for clinical practice based on whole-genome sequencing that could replace multiple diagnostics currently in use.</p>

<br>

<h2>Main Figures and Table</h2>

<h4>Fig 1. Recurrence within each tumour type in absolute numbers and percentages</h4>

<p align="center">
  <img src="Figures/Fig1.png" width="350" title="">
</p>

<strong>DATA:</strong> [TXT table](Recurrence/recurrence_tumortype_ssms_sims.txt) <br>

<br>
<h4>Fig 2. Spearman’s rank correlation between the 42 mutational features</h4>

<p align="center">
  <img src="Figures/Fig2.png" width="350" title="">
</p>

<strong>DATA:</strong> The statistics used to compute the correlations in the manuscript, including those in S2 Text. ([RData](Clustering/Data/sample2featuresAbsAndPerc.RData) | [TXT file](Clustering/Data/sample2featuresAbsAndPerc.txt))<br>

<strong>CODE:</strong> Correct p-values for multiple testing and plot correlations for the 42 features ([RScript](Clustering/Code/Fig2_correlationFeatures.R)) <BR> 


<br>
<h4>Fig 3. Workflow of the recurrence-based approach to group cancer genomes</h4>

![alt text](Figures/Fig3.png)

<p><em>Note</em>: An interactive <strong>3D version of the PCA plot</strong> shown in Step 3 of the workflow is available <a href="https://plot.ly/~biomedicalGenomicsCNAG/1.embed">here</a>.</p>

<strong>DATA:</strong><br>

* Output of Step 1: 42 features ([RData](Clustering/Data/step1_sample2features.RData) | [TXT file](Clustering/Data/step1_sample2features.txt))<br>
* Output of Step 3 and 4: PCA object ([RData](Clustering/Data/res_pca.RData))<br>
* Output of Step 5 and 6: HCPC object ([RData](Clustering/Data/res_hcpc.RData)) and Samples linked to cluster and tumour type ([RData](Clustering/Data/sample2ttype2cluster.RData) | [TXT file](Clustering/Data/sample2ttype2cluster.txt))<br>
* Output of Step 7: Samples annotated with metadata ([RData](Clustering/Data/step7_sample2annotation.RData) | [TXT file](Clustering/Data/step7_sample2annotation.txt))<br>

<strong>CODE:</strong><br>
* Main workflow, step 1 through 7 ([RScript](Clustering/Code/main_workflow.R)) <BR>
   * Annotation of mutations ([RScript](Clustering/Code/main_mutationLevelAnn.R)) <BR>
   * Annotation of samples ([RScript](Clustering/Code/main_sampleLevelAnn.R)) <BR>
* Generic functions ([RScript](Clustering/Code/main_utils.R)) <BR> 

<br>
<h4>Fig 4: Key characteristics of the 16 clusters</h4>

<p align="center">
  <img src="Figures/Table1_Fig4.png" width="350" title="">
</p>

<strong>DATA:</strong><br>

* Piediagram ([RData](Clustering/Data/sample2ttype2cluster.RData) | [TXT file](Clustering/Data/sample2ttype2cluster.txt)) <br>
* Number of SSMs and SIMs ([RData](Clustering/Data/sample2num_mut.RData) | [TXT file](Clustering/Data/sample2num_mut.txt)) <br>
* Key characteristics and overall association to recurrence: HCPC object ([RData](Clustering/Data/res_hcpc.RData))

<br>
<h4>Fig 5. Overview of the 42 features and their association with each cluster</h4>
<p align="center">
  <img src="Figures/Fig4.png" width="350" title="">
</p>
<strong>DATA:</strong><br>

PCA object ([RData](Clustering/Data/res_pca.RData))
<br>
All significant associations between features and clusters ([RData](Clustering/Data/descr_var_hcpc.RData))
<br>
<br>
<strong>CODE:</strong><br>
Correct the p-values for the results of the v tests (features versus clusters) for multiple testing and plot associations ([RScript](Clustering/Code/Fig4_feature_association_2_clusters.R)) <BR> <BR> 
  
<h4>Fig 6. Enriched sequence motifs</h4> 
<p align="center">
  <img src="Figures/Fig6.png" width="350" title="">
</p>  

<strong>METHOD:</strong> [Description](SequenceContext/Readme.md) <br>

<strong>DATA:</strong>

<h5>Cluster A (C>A mutations):</h5>

* Nucleotide counts per position in the sequence:<br> 
    * Recurrent mutations ([RData](SequenceContext/Data/Counts/SeqContext_counts_cluster_A_rec_C_A_21positions.RData) | [TXT file](SequenceContext/Data/Counts/SeqContext_counts_cluster_A_rec_C_A_21positions.txt))<br>
    * Unique mutations ([RData](SequenceContext/Data/Counts/SeqContext_counts_cluster_A_unique_C_A_21positions.RData) | [TXT file](SequenceContext/Data/Counts/SeqContext_counts_cluster_A_unique_C_A_21positions.txt))<br>
* Relative Entropies per nucleotide and Total Entropy per position in the sequence:<br>
    * Recurrent mutations ([TXT file](SequenceContext/Data/Entropies/SeqContext_RelativeEntropies_and_total_cluster_A_rec_C_A_21positions.txt))<br>
    * Unique mutations ([TXT file](SequenceContext/Data/Entropies/SeqContext_RelativeEntropies_and_total_cluster_A_unique_C_A_21positions.txt))<br>

<h5>Cluster E (C>G mutations):</h5>

* Nucleotide counts per position in the sequence:<br> 
    * Recurrent mutations ([RData](SequenceContext/Data/Counts/SeqContext_counts_cluster_E_rec_C_G_21positions.RData) | [TXT file](SequenceContext/Data/Counts/SeqContext_counts_cluster_E_rec_C_G_21positions.txt))<br>
    * Unique mutations ([RData](SequenceContext/Data/Counts/SeqContext_counts_cluster_E_unique_C_G_21positions.RData) | [TXT file](SequenceContext/Data/Counts/SeqContext_counts_cluster_E_unique_C_G_21positions.txt))<br>
* Relative Entropies per nucleotide and Total Entropy per position in the sequence:<br>
    * Recurrent mutations ([TXT file](SequenceContext/Data/Entropies/SeqContext_RelativeEntropies_and_total_cluster_E_rec_C_G_21positions.txt))<br>
    * Unique mutations ([TXT file](SequenceContext/Data/Entropies/SeqContext_RelativeEntropies_and_total_cluster_E_unique_C_G_21positions.txt))<br>
    
<h5>Cluster G (C>T mutations):</h5>

* Nucleotide counts per position in the sequence:<br> 
    * Recurrent mutations ([RData](SequenceContext/Data/Counts/SeqContext_counts_cluster_G_rec_C_T_21positions.RData) | [TXT file](SequenceContext/Data/Counts/SeqContext_counts_cluster_G_rec_C_T_21positions.txt))<br>
    * Unique mutations ([RData](SequenceContext/Data/Counts/SeqContext_counts_cluster_G_unique_C_T_21positions.RData) | [TXT file](SequenceContext/Data/Counts/SeqContext_counts_cluster_G_unique_C_T_21positions.txt))<br>
* Relative Entropies per nucleotide and Total Entropy per position in the sequence:<br>
    * Recurrent mutations ([TXT file](SequenceContext/Data/Entropies/SeqContext_RelativeEntropies_and_total_cluster_G_rec_C_T_21positions.txt))<br>
    * Unique mutations ([TXT file](SequenceContext/Data/Entropies/SeqContext_RelativeEntropies_and_total_cluster_G_unique_C_T_21positions.txt))<br>    

<h5>Cluster H (C>A mutations):</h5>

* Nucleotide counts per position in the sequence:<br> 
    * Recurrent mutations ([RData](SequenceContext/Data/Counts/SeqContext_counts_cluster_H_rec_C_A_21positions.RData) | [TXT file](SequenceContext/Data/Counts/SeqContext_counts_cluster_H_rec_C_A_21positions.txt))<br>
    * Unique mutations ([RData](SequenceContext/Data/Counts/SeqContext_counts_cluster_H_unique_C_A_21positions.RData) | [TXT file](SequenceContext/Data/Counts/SeqContext_counts_cluster_H_unique_C_A_21positions.txt))<br>
* Relative Entropies per nucleotide and Total Entropy per position in the sequence:<br>
    * Recurrent mutations ([TXT file](SequenceContext/Data/Entropies/SeqContext_RelativeEntropies_and_total_cluster_H_rec_C_A_21positions.txt))<br>
    * Unique mutations ([TXT file](SequenceContext/Data/Entropies/SeqContext_RelativeEntropies_and_total_cluster_H_unique_C_A_21positions.txt))<br>

<h5>Cluster L (T>G mutations):</h5>

* Nucleotide counts per position in the sequence:<br> 
    * Recurrent mutations ([RData](SequenceContext/Data/Counts/SeqContext_counts_cluster_L_rec_T_G_21positions.RData) | [TXT file](SequenceContext/Data/Counts/SeqContext_counts_cluster_L_rec_T_G_21positions.txt))<br>
    * Unique mutations ([RData](SequenceContext/Data/Counts/SeqContext_counts_cluster_L_unique_T_G_21positions.RData) | [TXT file](SequenceContext/Data/Counts/SeqContext_counts_cluster_L_unique_T_G_21positions.txt))<br>
* Relative Entropies per nucleotide and Total Entropy per position in the sequence:<br>
    * Recurrent mutations ([TXT file](SequenceContext/Data/Entropies/SeqContext_RelativeEntropies_and_total_cluster_L_rec_T_G_21positions.txt))<br>
    * Unique mutations ([TXT file](SequenceContext/Data/Entropies/SeqContext_RelativeEntropies_and_total_cluster_L_unique_T_G_21positions.txt))<br>

<strong>CODE:</strong> <BR>
* Add sequence context to SSMs and collect mutation of a specific subtype of a specific group of sample ([RScript](SequenceContext/Code/PrepareInput4SequenceLogos.R)) <br>
* Compute and plot sequence logos ([RScript](SequenceContext/Code/SequenceLogos.R)) <br>
* Compute statistics for sequence motifs ([RScript](SequenceContext/Code/SequenceMotifsStats.R)) <br>

<br>
<h4>Fig 7. Factors impacting on recurrence in the context of the clusters</h4> 
<p align="center">
  <img src="Figures/Fig7.png" width="350" title="">
</p>

<br>
<h2>Supporting Information</h2>

<h4>S1 Figure. Clustering tree showing tumour type distribution for 2 to 20 clusters.</h4>
<p align="center">
  <img src="Figures/S1Fig.png" width="350" title="">
</p>

<strong>DATA:</strong> <br>
Table with the samples and the clusters they are assigned to at each cluster resolution (from 2 to 20 clusters) [(TXT file)](Clustering/Data/HierarchClust20_sample2clust.txt) 
HCPC objects containing the association of the features to the clusters at a specific cluster resolution. ([Folder](Clustering/Data/HCPC_clustertree))


<br>
<h4>S2 Figure. PCA and clustering with and without the recurrence-related features.</h4>
<p align="center">
  <img src="Figures/S2Fig.png" width="350" title="">
</p>

<strong>DATA:</strong> <br>
Table with the samples and the clusters they are assigned to after clustering without and with the recurrence-related features (C and D) [(TXT file)](Clustering/Data/Clustering_sample2clustRec2clustGen.txt) 



<br>
<h4>S3 Figure. Enriched sequence motifs for C>G SSMs in cluster M.</h4>
<p align="center">
  <img src="Figures/S3Fig.png" width="350" title="">
</p>

<strong>METHOD:</strong> [Description](SequenceContext/Readme.md) <br>

<strong>DATA:</strong>

<h5>Cluster M (C>G mutations):</h5>

* Nucleotide counts per position in the sequence:<br> 
    * Recurrent mutations ([RData](SequenceContext/Data/Counts/SeqContext_counts_cluster_M_rec_C_G_21positions.RData) | [TXT file](SequenceContext/Data/Counts/SeqContext_counts_cluster_M_rec_C_G_21positions.txt))<br>
    * Unique mutations ([RData](SequenceContext/Data/Counts/SeqContext_counts_cluster_M_unique_C_G_21positions.RData) | [TXT file](SequenceContext/Data/Counts/SeqContext_counts_cluster_M_unique_C_G_21positions.txt))<br>
* Relative Entropies per nucleotide and Total Entropy per position in the sequence:<br>
    * Recurrent mutations ([TXT file](SequenceContext/Data/Entropies/SeqContext_RelativeEntropies_and_total_cluster_M_rec_C_G_21positions.txt))<br>
    * Unique mutations ([TXT file](SequenceContext/Data/Entropies/SeqContext_RelativeEntropies_and_total_cluster_M_unique_C_G_21positions.txt))<br>

<strong>CODE:</strong> <BR>
* Add sequence context to SSMs and collect mutation of a specific subtype of a specific group of sample ([RScript](SequenceContext/Code/PrepareInput4SequenceLogos.R)) <br>
* Compute and plot sequence logos ([RScript](SequenceContext/Code/SequenceLogos.R)) <br>
* Compute statistics for sequence motifs ([RScript](SequenceContext/Code/SequenceMotifsStats.R)) <br>

<br>

<h4>S4 Figure. Enriched sequence motifs for T>G SSMs in cluster H.</h4>
<p align="center">
  <img src="Figures/S4Fig.png" width="350" title="">
</p>

<strong>METHOD:</strong> [Description](SequenceContext/Readme.md) <br>

<strong>DATA:</strong>

<h5>Cluster H (T>G mutations):</h5>

* Nucleotide counts per position in the sequence:<br> 
    * Recurrent mutations ([RData](SequenceContext/Data/Counts/SeqContext_counts_cluster_H_rec_T_G_21positions.RData) | [TXT file](SequenceContext/Data/Counts/SeqContext_counts_cluster_H_rec_T_G_21positions.txt))<br>
    * Unique mutations ([RData](SequenceContext/Data/Counts/SeqContext_counts_cluster_H_unique_T_G_21positions.RData) | [TXT file](SequenceContext/Data/Counts/SeqContext_counts_cluster_H_unique_T_G_21positions.txt))<br>
* Relative Entropies per nucleotide and Total Entropy per position in the sequence:<br>
    * Recurrent mutations ([TXT file](SequenceContext/Data/Entropies/SeqContext_RelativeEntropies_and_total_cluster_H_rec_T_G_21positions.txt))<br>
    * Unique mutations ([TXT file](SequenceContext/Data/Entropies/SeqContext_RelativeEntropies_and_total_cluster_H_unique_T_G_21positions.txt))<br>

<strong>CODE:</strong> <BR>
* Add sequence context to SSMs and collect mutation of a specific subtype of a specific group of sample ([RScript](SequenceContext/Code/PrepareInput4SequenceLogos.R)) <br>
* Compute and plot sequence logos ([RScript](SequenceContext/Code/SequenceLogos.R)) <br>
* Compute statistics for sequence motifs ([RScript](SequenceContext/Code/SequenceMotifsStats.R)) <br>

---

<h4>S1 Table. Tumour type abbreviation, full name and number of samples.</h4>

<h4>S2 Table. Recurrence in pan-cancer context and within tumour type(s).</h4>

<strong>CODE:</strong><BR> Determine the type of recurrence for each recurrent mutation ('pan-cancer only', 'single tumour type', 'multiple tumour types') ([RScript](Recurrence/Code/TypeOfRecurrence.R)) <br>
  
---

<h4>S1 Text. Estimation of the levels of recurrence when purely driven by chance.</h4>

<strong>METHOD:</strong> [Description](Simulations/README.md) <br>

<strong>DATA:</strong><br>

* Summary tables of 5000 simulations:<br> 
    * Recurrence pan-cancer per tumour type ([RData](Simulations/Data/recurrence_tumortype_pancancer_5000simul.RData) | [TXT file](Simulations/Data/recurrence_tumortype_pancancer_5000simul.txt))<br>
    * Recurrence within tumour type ([RData](Simulations/Data/recurrence_tumortype_within_5000simul.RData) | [TXT file](Simulations/Data/recurrence_tumortype_within_5000simul.txt))<br>
    * Recurrence per substitution type ([RData](Simulations/Data/recurrence_subtype_5000simul.RData) | [TXT file](Simulations/Data/recurrence_subtype_5000simul.txt))<br>
    * Overall recurrence ([RData](Simulations/Data/recurrence_overall_5000simul.RData) | [TXT file](Simulations/Data/recurrence_overall_5000simul.txt))<br>

<strong>CODE:</strong><br>

* [Performance of simulations](Simulations/Code/sim-all-cancers.R) <br>
* [Summarise per tumour type](Simulations/Code/summarise.R)<br>
* [Summarise per substitution type](Simulations/Code/summarise-mut-type.R)<br>
* [Boxplot representation](Simulations/Code/representation-simulations.R)<br>

<br>
<h4>S2 Text. Recurrence versus general mutational characteristics.</h4>

<strong>DATA: </strong>Features in absolute terms and percentages ([RData](Clustering/Data/sample2featuresAbsAndPerc.RData) | [TXT file](Clustering/Data/sample2featuresAbsAndPerc.txt))<br>

<strong>CODE:</strong> Estimate for frequency of homopolymers in the genome ([RScript](Recurrence/Code/S2_Text_Table_STX2_6.R)) <br>
<BR>   
<h4>S3 Text. Detailed cluster-specific descriptions.</h4>

<strong>DATA: </strong>Features in absolute terms and percentages plus metadata ([RData](Clustering/Data/sample2features2annotation.RData) | [TXT file](Clustering/Data/sample2features2annotation.txt))<br>

<strong>CODE:</strong> <BR>
* Density plots of the replication time scores ([RScript](Recurrence/Code/S3_Text_Figs_STX3_2_3.R))
* Plot the number of recurrent mutations along the genome for Lymph-BNHL and Lymph-CLL in cluster M. ([RScript](Recurrence/Code/S3_Text_Fig_STX3_17.R)) <br>

<BR>
<h4>S4 Text. Smoking history and related mutational subtypes.</h4>

---

<h4>S1 File. Characteristic plots summarising each of the 42 features.</h4>

<strong>DATA: </strong>Features in absolute terms and percentages ([RData](Clustering/Data/sample2featuresAbsAndPerc.RData) | [TXT file](Clustering/Data/sample2featuresAbsAndPerc.txt))<br>
<BR>
<h4>S2 File. Sample distribution per tumour type across the 16 clusters.</h4>

<strong>DATA: </strong>Samples linked to cluster and tumour type ([RData](Clustering/Data/sample2ttype2cluster.RData) | [TXT file](Clustering/Data/sample2ttype2cluster.txt))<br>
