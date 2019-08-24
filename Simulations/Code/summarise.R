if (interactive()){
    olapfn <- "olap1.rds"
    allfn <- "all1.rds"
} else {
    args <- commandArgs(trailingOnly=TRUE)
    olapfn <- args[1]
    allfn <- args[2]
}
cat("reading olap\n", file=stderr())
if (!exists("olap")){
    olap <- readRDS(olapfn)
    cat(dim(olap),"\n", file=stderr())
}
cat("reading allg\n", file=stderr())
if (!exists("allg")){
    allg <- readRDS(allfn)
    cat(dim(allg), "\n", file=stderr())
}

for (abbr in unique(allg$tumor_type)){
#for (abbr in c("Bladder-TCC")){
    all_sel <- allg[allg$tumor_type==abbr, ]
    # non unique rows of course
    # for example a recurrent mutation which
    # appears twice in the same tumor type (Breast-TCC)
    # corresponds to two lines with the same 
    # type and position.
    cat("hashing...\n", file=stderr())
    hash_sel <- paste0(all_sel$type, ":", all_sel$pos)
    cat("hashing:table...\n", file=stderr())
    tbl_sel <- table(hash_sel)
    ssm <- length(tbl_sel)
    cat("merging with olap...\n", file=stderr())
    tmp     <- merge(all_sel, olap, 
        by.x=c("type", "pos"), by.y=c("type", "pos"))
    # it would have been nice to write a comment here
    # in due time.
    tmp <- tmp[, -c(3, 4)]
    tmp <- tmp[!duplicated(tmp),]
    rec_ssm <- sum(tmp[, abbr] > 0)
    tumor_specific_rec_ssm <- nrow(tmp[tmp[, abbr]>=2,])
    cat(sprintf("%s\t%d\t%d\t%s\t%d\t%s\n", 
        abbr, ssm, rec_ssm, 
        format(rec_ssm * 100 / ssm, digits=4),
        tumor_specific_rec_ssm, 
        format(tumor_specific_rec_ssm *100 / ssm, digits=4)))
}
tot_ssm <- nrow(allg) - sum(olap$freq - 1)
tot_rec_ssm <- nrow(olap)
cat(sprintf("%s\t%d\t%d\t%s\t%d\t%s\n", 
        "Overall", tot_ssm, tot_rec_ssm, 
        format(tot_rec_ssm * 100 / tot_ssm, digits=4),
        NA, NA))


