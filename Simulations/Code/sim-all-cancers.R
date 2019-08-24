library(plyr)
library(parallel)

cores <- 3

source("sim-snp-lib.R")

if (interactive()){
    outfn <- "olap-test.rds"
    outfn2 <- "allg.rds"
    rndseed <- 1
} else{
    args  <- commandArgs(trailingOnly = TRUE)
    outfn <- args[1]
    outfn2 <- args[2]
    rndseed <- as.integer(args[3])
    # for example :
    # Rscript sim-all-cancers.R olap-test.rds allg.rds
}
set.seed(rndseed)
load("stats_ssms2tumor_type.RData")
muts <- stats_ssms2tumor_type
g <- list()
st_time <- proc.time()[3]
istart <- 1 
istop <- nrow(muts) 
aux_make_genome <- function(i){
    # here i is the identifier of a sample 
    # out of ~ 2500
    # or better said identifies a row in muts
    # which is the original table passed on
    # by Miranda.
    cat(sprintf("%d-start\n",i))
    nmut <- c(muts[i, "num_CT"], muts[i,"num_CA"], muts[i,"num_CG"],
    muts[i, "num_TA"], muts[i, "num_TC"], muts[i, "num_TG"])
    legr <- c(1144530852, 2861327131 - 1144530852)
    totmut <- sum(nmut)
    res <- make_genome(i, legr, nmut)
    # res is a data frame with nmut rows
    res$sid <- rep(muts[i, "sample_id"], totmut)
    res$tumor_type <- rep(muts[i, "tumor_type"], totmut)
    cat(sprintf("%d-end\n", i))
    current_time <- proc.time()[3] - st_time
    cat(sprintf("%d\t%.1f\t%.1f\t%d\n", i, current_time, 
    current_time/(i - istart + 1), totmut))
    return(res) 
}
idx <- as.list(seq(istart, istop))
st_time <- proc.time()[3]
g <- parallel::mclapply(idx, aux_make_genome, mc.cores=cores)
allg <- plyr::rbind.fill(g)
current_time <- proc.time()[3] - st_time
cat(sprintf("computing intersection after %.1f secs\n", current_time))
olap <- find_overlapping_pos(allg)
olap2  <- olap[olap$freq >= 2, ]
cat("computing allg2...\n", file=stderr())
allg2 <- allg[allg$pos %in% olap2$pos, ]
cat(dim(olap2), "\n", file=stderr())
# in this merge lines will be repeated
# because there are mutations of the same type
# hitting the same position in different genomes.
tbl <- merge(olap2, allg2, by.x=c("type", "pos"), by.y=c("type","pos"))
cat("aggregating to find tumor of origin of the som mut...\n", file=stderr())
out_tbl <- aggregate(tumor_type ~ type + pos, data=tbl, 
    FUN=function(v){
        return(paste0(v, collapse=":"))
        })
nr <- nrow(out_tbl)
filler <- rep(0, nr)
ttype_counts <- data.frame(
  "Bladder-TCC" = filler ,        "Breast-AdenoCA"=filler,      "Ovary-AdenoCA"=filler,      
  "Panc-Endocrine"=filler,      "Prost-AdenoCA"=filler,       "Kidney-RCC"=filler ,        
  "Skin-Melanoma"=filler  ,     "Stomach-AdenoCA"=filler,     "Thy-AdenoCA"=filler,        
  "Liver-HCC"=filler ,          "CNS-Medullo"=filler  ,       "Cervix-SCC"=filler   ,      
  "Cervix-AdenoCA"=filler,      "Panc-AdenoCA"=filler  ,      "Myeloid-AML"=filler  ,      
  "Eso-AdenoCA"=filler ,        "Lymph-CLL"=filler ,          "Head-SCC"=filler  ,         
  "Bone-Epith"=filler   ,       "Bone-Benign"=filler ,        "Bone-Osteosarc"=filler,     
  "Lymph-BNHL"=filler,          "Myeloid-MPN"=filler,         "Myeloid-MDS"=filler,        
  "Biliary-AdenoCA"=filler,     "Breast-LobularCA"=filler,    "Breast-DCIS"=filler ,       
  "ColoRect-AdenoCA"=filler,    "SoftTissue-Liposarc"=filler, "SoftTissue-Leiomyo"=filler, 
  "Kidney-ChRCC"=filler,        "CNS-GBM"=filler,             "CNS-Oligo"=filler,          
  "Lung-AdenoCA"=filler,        "Lung-SCC"=filler ,           "CNS-PiloAstro"=filler,      
  "Uterus-AdenoCA"=filler, check.names=FALSE )
for (i in seq(1,nr)){
    tts <- strsplit(out_tbl$tumor_type[i], ":")[[1]]
    for (j in seq(1, length(tts))){
       ttype_counts[i, tts[j]] <- ttype_counts[i, tts[j]] + 1 
    }
}
out_tbl <- cbind(out_tbl, ttype_counts)
cat(sprintf("writing to disk...\n"))
saveRDS(object = out_tbl, file=outfn)
saveRDS(object = allg, file=outfn2)
