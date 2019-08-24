library(plyr)

if (interactive()){
    allfn <- "all1.rds"
} else {
    args <- commandArgs(trailingOnly=TRUE)
    allfn <- args[1]
}
cat("reading allg\n", file=stderr())
if (!exists("allg")){
    allg <- readRDS(allfn)
}
utypes <- c("c_t", "c_a", "c_g", "t_a", "t_c", "t_g")
allgl <- list()
for (ty in utypes){
    allgl[[ty]] <- allg[allg$type==ty,]
    cat(dim(allgl[[ty]]), "\n", file=stderr())
    cnts        <- count(allgl[[ty]], vars="pos")
    totuniq     <- nrow(cnts)
    rec         <- nrow(cnts[cnts$freq>=2,])
    cat(ty, totuniq, rec, "\n")
}
