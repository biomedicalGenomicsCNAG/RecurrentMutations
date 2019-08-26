
library(RMySQL)
library(parallel)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(GenomicFeatures)
library(Rsamtools)

load(paste("/project/devel/PCAWG/landscape_recurrence/results/stats_pcawg_annotated2clinical.RData",sep=""))

stats <- stats_pcawg_annotated2clinical
stats$clust_char <- factor(as.character(stats$clust_char), levels = sort(unique(as.character(stats$clust_char))))

t_sample_indices_clustM_CLL <- stats[which(stats$cluster == "M" & stats$tumor_type =="Lymph-CLL"),"sample_id"]
t_sample_indices_clustM_BNHL <- stats[which(stats$cluster == "M" & stats$tumor_type =="Lymph-BNHL"),"sample_id"]


getStatsTilePerTtype("tile_rec_mutations", c(t_sample_indices_clustM_CLL,t_sample_indices_clustM_BNHL),"lymph_bnhl_cll_clustM")

getTiledGenome <- function(num_mutations, mc_cores)
{
  genome_fasta <- open(FaFile(fastafileGenome)) 
  
  length_chr <- seqlengths(genome_fasta)
  length_chr <- length_chr[c(as.character(1:22),"X","Y")]
  size_genome <- sum(as.numeric(length_chr))
  tile_size <- ceiling((10*size_genome) / num_mutations)
  
  if((tile_size %% 2) != 0)
    tile_size <- tile_size + 1
  
  size_overlap_tiles <- tile_size/2
  
  tiles <- tileGenome(seqlengths = length_chr, tilewidth = tile_size, cut.last.tile.in.chrom = TRUE)
  tmp_tiles <- tiles
  
  tmp_tiles <- trim(shift(tmp_tiles, size_overlap_tiles))
  
  tmp_tiles_corrected <- do.call(c, mclapply(names(length_chr), function(x)
  {
    cur_tiles <- tmp_tiles[which(seqnames(tmp_tiles) %in% x),]
    print(cur_tiles)
    
    if(any(start(cur_tiles) >= length_chr[x] | end(cur_tiles) > length_chr[x]))
      cur_tiles <- cur_tiles[-which(start(cur_tiles) >= length_chr[x] | end(cur_tiles) > length_chr[x]),]
    
    first_partial_tile <- cur_tiles[1,]
    start(first_partial_tile) <- 1
    end(first_partial_tile) <- size_overlap_tiles
    
    cur_tiles <- c(cur_tiles, first_partial_tile)
    
    return(cur_tiles)
  }, mc.cores=mc_cores))  
  
  overlapping_tiles <- sort(c(tiles, tmp_tiles_corrected))
  
  overlapping_tiles$tile_index <- 1:length(overlapping_tiles)
  
  return(overlapping_tiles)
}


matchTile2Mut <- function(mutations_gr, tiled_genome, mut_type, table_name)
{
  # overlap between mutations and tiled genome
  hits_mutationsVStiles <- findOverlaps(mutations_gr, subject = tiled_genome)
  
  mut_hits <- as.data.frame(mutations_gr[queryHits(hits_mutationsVStiles), ])
  colnames(mut_hits) <- paste("mut_", colnames(mut_hits), sep="")
  
  tiled_genome_hits <- as.data.frame(tiled_genome[subjectHits(hits_mutationsVStiles), ])
  colnames(tiled_genome_hits) <- paste("tile_", colnames(tiled_genome_hits), sep="")
  mut2tile_info <- cbind(mut_hits, tiled_genome_hits)
  
  mut2tile <- mut2tile_info[, c("mut_mut_index", "tile_tile_index")]
  colnames(mut2tile) <- c(paste(mut_type, "_index",sep=""), "tile_index")
  
  numsamples2mut <- dbGetQuery(con, paste("SELECT ", mut_type, "_index, num_samples FROM ", mut_type, "_filtered2stats_v2;", sep=""))
  
  mut2tile <- merge(mut2tile, numsamples2mut, by=paste(mut_type, "_index",sep=""), all.x=TRUE)
  
  if(mut_type == "ssm")
    mut2class <- dbGetQuery(con, paste("SELECT ", mut_type, "_index, substitution_class FROM ", mut_type, "2signature_info_v2;", sep=""))
  else
    mut2class <- dbGetQuery(con, paste("SELECT ", mut_type, "_index, type FROM ", mut_type, "_v2;", sep=""))
  
  mut2tile <- merge(mut2tile, mut2class, by=paste(mut_type, "_index",sep=""), all.x=TRUE)
  dbWriteTable(conn = con, name = table_name, value = mut2tile, overwrite=FALSE, append=TRUE, row.names=0)
}


mainDir <- "/project/devel/PCAWG/somatic_consensus_aug2016/"
setwd(mainDir)

dataDir <- paste(mainDir, "data/",sep="")
resultsDir <-  paste(mainDir, "results/",sep="")

con = dbConnect(MySQL(), group='client', dbname="mstobbe_schema")
source("/project/devel/PCAWG/somatic_consensus_aug2016/scripts/pcawg.colour.palette.R")


getStatsTilePerTtype <- function(name_tiled_genome, t_sample_indices, group_name)
{
  load(paste(dataDir,"t_sample2tumortype.RData", sep=""))
  
  con = dbConnect(MySQL(), group='client', dbname="mstobbe_schema")
  
  # all tiles
  all_tile_indices <- dbGetQuery(con, paste("SELECT tile_index FROM ", name_tiled_genome, ";", sep=""))[,1]  
  
  # all tiles linked to mutations
  tile2ssm <- dbGetQuery(con, paste("SELECT tile_index, ssm_index FROM ssm2", name_tiled_genome, ";",sep=""))
  tile2sim <- dbGetQuery(con, paste("SELECT tile_index, sim_index FROM sim2", name_tiled_genome,";", sep="")) 
  
  # mutations linked to samples from specific tumor type
  ssm2tsample_ttype <- dbGetQuery(con, paste("SELECT ssm_index, t_sample_index FROM ssm_filtered2tumor_sample_v2 WHERE t_sample_index IN ('", paste(t_sample_indices, collapse="', '"), "');", sep=""))
  sim2tsample_ttype <- dbGetQuery(con, paste("SELECT sim_index, t_sample_index FROM sim_filtered2tumor_sample_v2 WHERE t_sample_index IN ('", paste(t_sample_indices, collapse="', '"), "');", sep=""))
  
  # mutations from specific tumor type
  ssm_indices_ttype <- unique(ssm2tsample_ttype$ssm_index)
  sim_indices_ttype <- unique(sim2tsample_ttype$sim_index)
  
  # number of samples per SSM in specific tumor type
  ssm2num_samples <- table(ssm2tsample_ttype$ssm_index)
  ssm2num_samples_df <- as.data.frame(ssm2num_samples)
  ssm2num_samples_df <- ssm2num_samples_df[which(as.character(ssm2num_samples_df$Var1) %in% as.character(unique(tile2ssm$ssm_index))),]
  rownames(ssm2num_samples_df) <- as.character(ssm2num_samples_df$Var1)
  
  sim2num_samples <- table(sim2tsample_ttype$sim_index)
  sim2num_samples_df <- as.data.frame(sim2num_samples)
  sim2num_samples_df <- sim2num_samples_df[which(as.character(sim2num_samples_df$Var1) %in% as.character(unique(tile2sim$sim_index))),]
  rownames(sim2num_samples_df) <- as.character(sim2num_samples_df$Var1)
  
  tile_indices_with_ssm <- tile2ssm[which(tile2ssm$ssm_index %in% ssm_indices_ttype), "tile_index"]
  tile_indices_with_sim <- tile2sim[which(tile2sim$sim_index %in% sim_indices_ttype), "tile_index"]
  
  tile_indices_with_hits <- unique(c(tile_indices_with_ssm, tile_indices_with_sim))
  tile_indices_without_hits <- all_tile_indices[which(!(all_tile_indices %in% tile_indices_with_hits))]
  
  create_table_query <- paste("CREATE TABLE IF NOT EXISTS `mstobbe_schema`.`", name_tiled_genome,"_", group_name, "2stats` (`tile_rec_mutations2stats_index` INT(4) NOT NULL AUTO_INCREMENT,`tile_index` INT(4) NULL,`num_mutations` SMALLINT(2) NULL,PRIMARY KEY (`tile_rec_mutations2stats_index`), UNIQUE INDEX `unique_tile_index` (`tile_index` ASC));", sep="")
  tmp <- dbSendQuery(con,create_table_query)
  
  tile_stats_without_stats <- as.data.frame(matrix(data=0, nrow=length(tile_indices_without_hits), ncol=2))
  colnames(tile_stats_without_stats) <- c("tile_index", "num_mutations")
  tile_stats_without_stats[, "tile_index"] <- tile_indices_without_hits
  
  dbWriteTable(conn = con, name = paste(name_tiled_genome,"_", group_name, "2stats",sep=""), value = tile_stats_without_stats, overwrite=FALSE, append=TRUE, row.names=0)
  
  dbDisconnect(con) 
  
  tmp <- mclapply(1:length(tile_indices_with_hits), function(i)
  {
    con = dbConnect(MySQL(), group='client', dbname="mstobbe_schema")
    
    cur_stats <- as.data.frame(matrix(data=0, nrow=1, ncol=2))
    colnames(cur_stats) <- c("tile_index", "num_mutations")
    
    print(i)
    ssm_hits <- dbGetQuery(con, paste("SELECT * FROM ssm2tile_rec_mutations WHERE tile_index=", tile_indices_with_hits[i] ,";", sep=""))  
    sim_hits <- dbGetQuery(con, paste("SELECT * FROM sim2tile_rec_mutations WHERE tile_index=",tile_indices_with_hits[i] ,";", sep=""))   
    
    ssm_hits <- ssm_hits[which(ssm_hits$ssm_index %in% ssm_indices_ttype), ]
    sim_hits <- sim_hits[which(sim_hits$sim_index %in% sim_indices_ttype), ]
    
    #ssm_hits$num_samples <- ssm2num_samples[as.character(ssm_hits$ssm_index)]
    #sim_hits$num_samples <- sim2num_samples[as.character(sim_hits$sim_index)]
    
    ssm_hits$num_samples <- ssm2num_samples_df[as.character(ssm_hits$ssm_index), "Freq"]
    sim_hits$num_samples <- sim2num_samples_df[as.character(sim_hits$sim_index), "Freq"]
    
    ssm_hits_rep <- ssm_hits[rep(seq_len(nrow(ssm_hits)), ssm_hits$num_samples),]
    
    sim_hits_rep <- sim_hits[rep(seq_len(nrow(sim_hits)), sim_hits$num_samples),]
    
    cur_stats[1,"tile_index"] <- tile_indices_with_hits[i]
    cur_stats[1,"num_ssms"] <- nrow(ssm_hits_rep)
    cur_stats[1,"num_sims"] <- nrow(sim_hits_rep)
    
    cur_stats[1,"num_mutations"] <-  cur_stats[1,"num_sims"] +  cur_stats[1,"num_ssms"]
    
    tmp <- dbSendQuery(con, paste("INSERT INTO ", name_tiled_genome,"_", group_name, "2stats (tile_index, num_mutations) VALUES ('", cur_stats$tile_index, "', '", cur_stats$num_mutations, "');",  sep=""))
    dbClearResult(tmp)
    
    dbDisconnect(con) 
    
    return(0)
  }, mc.cores=4)
  
}


load(paste("/project/devel/PCAWG/landscape_recurrence/results/stats_pcawg_annotated2clinical.RData",sep=""))

stats <- stats_pcawg_annotated2clinical
stats$clust_char <- factor(as.character(stats$clust_char), levels = sort(unique(as.character(stats$clust_char))))

t_sample_indices_clustM_CLL <- stats[which(stats$clust_char == "M" & stats$tumor_type =="Lymph-CLL"),"sample_id"]
t_sample_indices_clustM_BNHL <- stats[which(stats$clust_char == "M" & stats$tumor_type =="Lymph-BNHL"),"sample_id"]

getStatsTilePerTtype("tile_rec_mutations",t_sample_indices_clustM_CLL,"lymph_cll_clustM")
getStatsTilePerTtype("tile_rec_mutations", t_sample_indices_clustM_BNHL,"lymph_bnhl_clustM")
getStatsTilePerTtype("tile_rec_mutations", c(t_sample_indices_clustM_CLL,t_sample_indices_clustM_BNHL),"lymph_bnhl_cll_clustM")


getGeneLabels <- function(tile_stats, tile_stats_chr, threshold, table_name)
{
  tile_minTh <- tile_stats[which(tile_stats$num_mutations >= threshold),"tile_index"]
  
  
  tile2gene <- dbGetQuery(con, paste("SELECT tile_index, gene_name FROM ", table_name, ";",sep="")) 
  tile_minTh2gene <- tile2gene[which(tile2gene$tile_index %in% tile_minTh),]
  
  #temp filter
  tile_minTh2gene <- tile_minTh2gene[-grep("^AL",tile_minTh2gene$gene_name),]
  tile_minTh2gene <- tile_minTh2gene[-grep("^AC0",tile_minTh2gene$gene_name),]
  tile_minTh2gene <- tile_minTh2gene[-grep("^RP11",tile_minTh2gene$gene_name),]
  tile_minTh2gene <- tile_minTh2gene[-grep("^hsa-mir-",tile_minTh2gene$gene_name),]
  tile_minTh2gene <- tile_minTh2gene[-grep("^MIR",tile_minTh2gene$gene_name),]
  tile_minTh2gene <- tile_minTh2gene[-grep("^LINC",tile_minTh2gene$gene_name),]
  if(length(grep("^LL",tile_minTh2gene$gene_name))> 0)
    tile_minTh2gene <- tile_minTh2gene[-grep("^LL",tile_minTh2gene$gene_name),]
  
  
  tile_indices_minTh <- unique(tile_minTh2gene$tile_index)
  tile_minTh_summary <- as.data.frame(matrix(nrow=length(tile_indices_minTh), ncol=2))
  colnames(tile_minTh_summary) <- c("tile_index", "gene_name")
  
  for(i in 1:length(tile_indices_minTh)) 
  {
    tile_minTh_summary[i,"tile_index"] <-  tile_indices_minTh[i]
    tile_minTh_summary[i,"gene_name"] <- paste(tile_minTh2gene[which(tile_minTh2gene$tile_index == tile_indices_minTh[i]),"gene_name"], collapse=" / ")
  }
  
  
  tile_minTh_summary <- merge(tile_minTh_summary,tile_stats_chr, by="tile_index", all.x=TRUE)
  
  genes <- unique(tile_minTh_summary$gene_name)
  tile_minTh_summary_filt <- as.data.frame(matrix(nrow=length(genes), ncol=5))
  colnames(tile_minTh_summary_filt) <- c("tile_index", "mid_pos", "num_mutations", "chr", "gene_name")
  
  for(i in 1:length(genes))
  {
    tile_minTh_summary_filt[i,"gene_name"] <- genes[i]
    cur_tiles <- tile_minTh_summary[which(tile_minTh_summary$gene_name == genes[i]), ]
    max_num_mut <- max(cur_tiles$num_mutations)
    tile_minTh_summary_filt[i,"num_mutations"] <- max_num_mut
    tile_minTh_summary_filt[i,"tile_index"] <- cur_tiles[which(cur_tiles$num_mutations == max_num_mut), "tile_index"][1]
    tile_minTh_summary_filt[i,"mid_pos"] <- cur_tiles[which(cur_tiles$num_mutations == max_num_mut), "mid_pos"][1]
    tile_minTh_summary_filt[i,"chr"] <- cur_tiles[which(cur_tiles$num_mutations == max_num_mut), "chr"][1]
  }
  
  return(tile_minTh_summary_filt)
}

plotStatsSpecTtype <- function(group_name, threshold)
{
  con = dbConnect(MySQL(), group='client', dbname="mstobbe_schema")  
  tile_stats <- dbGetQuery(con, paste("SELECT * FROM tile_rec_mutations_", group_name, "2stats;", sep="")) 
  
  tile2chr <- dbGetQuery(con, "SELECT * FROM tile_rec_mutations;") 
  tile2chr$mid_pos <- floor(tile2chr$end_pos/2)
  
  chromosomes <- unique(as.character(tile2chr$chr),sep="")
  chromosomes <- sub("23", "X",chromosomes)
  chromosomes <- sub("24", "Y",chromosomes)
  colors_chr <- pcawg.colour.palette(chromosomes, "chromosomes")
  names(colors_chr) <- unique(as.character(tile2chr$chr),sep="")
  
  
  
  tile_stats_chr <- merge(tile_stats, tile2chr, by="tile_index", all.x=TRUE)
  
  tile_min10_summary_filt <- getGeneLabels(tile_stats, tile_stats_chr, threshold,"tile_rec_mutations2gene")
  #tile_min10_summary_filt <- tile_min10_summary_filt[which(!(tile_min10_summary_filt$tile_index %in% c(59290,92967, 190740,336052,404012,404094,404069, 404177, 404157, 404171, 467753, 499849,499854, 499817,499826, 499843,462945,499698,499851,404148,499848,499813, 404130))),]
  # tile_min10_summary_filt[which(tile_min10_summary_filt$tile_index=="499850"), "gene_name"] <- "IGLL5 / IGLC1 / IGLJ2 / IGLC2 / IGLJ3 / IGLC3 / IGLJ4 / IGLC4"
  # tile_min10_summary_filt[which(tile_min10_summary_filt$tile_index=="499847"), "gene_name"] <- "IGLV3-1 / IGLL5 / IGLJ1"
  # tile_min10_summary_filt[which(tile_min10_summary_filt$tile_index=="499812"), "gene_name"] <- "IGLVVI-25-1 / IGLV3-25 / IGLV3-24"
  # 
  # tile_min10_summary_filt[which(tile_min10_summary_filt$tile_index=="404129"), "gene_name"] <- "IGHVII-22-1 / IGHVIII-22-2 / IGHV3-23 / IGHV1-24"
  # tile_min10_summary_filt[which(tile_min10_summary_filt$tile_index=="404147"), "gene_name"] <- "IGHVII-33-1 / IGHV3-33-2 / IGHV4-34 / IGHV7-34-1"
  # 
  # tile_min10_summary_filt[which(tile_min10_summary_filt$tile_index=="404070"), "gene_name"] <- "IGHD3-3 / IGHD2-2 / KIAA0125 / IGHD1-1 / IGHD1-7 / IGHD6-6 / IGHD5-5 / IGHD4-4"
  # 
  # tile_min10_summary_filt[which(tile_min10_summary_filt$tile_index=="499846"), "gene_name"] <- "IGLV3-2 / IGLV3-1 / IGLL5 / IGLJ1 / IGLC1 / IGLJ2"
  # 
  tile_min10_summary_filt$chr <- factor(tile_min10_summary_filt$chr)
  
  
  tile_stats_chr <- tile_stats_chr[which(tile_stats_chr$num_mutations>0),]
  tmp_tile_stats_chr <- tile_stats_chr
  tmp_tile_stats_chr$chr <- as.numeric(as.character(tmp_tile_stats_chr$chr))
  tmp_tile_stats_chr[which(tmp_tile_stats_chr$chr == 23), "chr"] <- "X"
  tmp_tile_stats_chr[which(tmp_tile_stats_chr$chr == 24), "chr"] <- "Y"
  tmp_tile_stats_chr$chr <- factor(tmp_tile_stats_chr$chr, levels=c(1:22, "X","Y"))
  tmp_colors_chr <- colors_chr
  names(tmp_colors_chr)[23:24] <- c("X","Y")
  tmp_tile_min10_summary_filt <- tile_min10_summary_filt
  tmp_tile_min10_summary_filt$chr <- as.numeric(as.character(tmp_tile_min10_summary_filt$chr))
  tmp_tile_min10_summary_filt[which(tmp_tile_min10_summary_filt$chr == 23), "chr"] <- "X"
  tmp_tile_min10_summary_filt[which(tmp_tile_min10_summary_filt$chr == 24), "chr"] <- "Y"
  
  tile_stats_chr <- tile_stats_chr[which(tile_stats_chr$num_mutations>0),]
  scatterplot_absNumSIMvsabsSSM <- ggplot(aes(x=tile_index,y=num_mutations,  fill=chr),  data=tmp_tile_stats_chr) + scale_y_log10() +scale_fill_manual(values=tmp_colors_chr, name="chromosome") +  ylab("number of mutations") +   #scale_fill_manual(values=colors_cust) + ylab("absolute number of SSMs") + #+ #ylab(paste("num ", mutation_type, sep="")) +  
    geom_point(size=2, shape=21) +  theme( axis.title.y=element_text(size=12), axis.title.x=element_blank(),axis.text.x =element_blank())  + annotation_logticks(sides = "l") + geom_text(data=tmp_tile_min10_summary_filt, aes(x=tile_index, y=num_mutations, label=gene_name), size=2, vjust=-1)
  ggsave(file=paste(resultsDir, "tiling/scatterplot_num_mut_rec_mut_", group_name, "th", threshold,"_noXLabs_withlegend.pdf", sep=""))
  
  scatterplot_absNumSIMvsabsSSM <- ggplot(aes(x=tile_index,y=num_mutations,  fill=chr),  data=tmp_tile_stats_chr) + guides(fill=FALSE) + scale_y_log10() +scale_fill_manual(values=tmp_colors_chr, name="chromosome") +  ylab("number of mutations") +   #scale_fill_manual(values=colors_cust) + ylab("absolute number of SSMs") + #+ #ylab(paste("num ", mutation_type, sep="")) +  
    geom_point(size=2, shape=21) +  theme( axis.title.y=element_text(size=12), axis.title.x=element_blank(),axis.text.x =element_blank())  + annotation_logticks(sides = "l") + geom_text(data=tmp_tile_min10_summary_filt, aes(x=tile_index, y=num_mutations, label=gene_name), size=2, vjust=-1)
  ggsave(file=paste(resultsDir, "tiling/scatterplot_num_mut_rec_mut_", group_name, "th", threshold,"_noXLabs_nolegend.pdf", sep=""))
  
  
  
  dbDisconnect(con)
}