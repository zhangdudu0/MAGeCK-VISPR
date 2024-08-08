# redirect all output to the log file
log <- file(snakemake@log[[1]])
sink(log, append = TRUE)
sink(log, append = TRUE, type = "message")

library(sva)

# read count
count <- read.table(snakemake@input[["counts"]], sep = "\t", header = TRUE, row.names=1)
samples <- colnames(count)[2:ncol(count)]
# read batch matrix
info <- read.table(snakemake@input[["batchmatrix"]], sep = "\t", header = TRUE, row.names=1)
info <- info[samples, , drop = FALSE]
# get batch info
batch <- factor(info[, "batch"])
# design matrix
if(ncol(info) >= 3) {
  # use covariates
  mod <- model.matrix(~info[, 3:ncol(info)])
} else {
  # no covariates
  mod <- model.matrix(~1, data=info)
}
# log-transform counts
log_count <- as.matrix(log(count[, 2:ncol(count)] + 1, 2))
sgrnas <- rownames(log_count)
nonzero <- apply(log_count, 1, var) != 0
# remove batch with combat
combat_data <- ComBat(dat = log_count[nonzero, ], batch = batch, mod = mod)
# regularize corrected counts to get rid of minor numerical instabilities
combat_data <- pmax(combat_data, 0)
# generate counts again
combat_count <- 2 ^ combat_data - 1
combat_count <- rbind(combat_count, log_count[!nonzero, ])[sgrnas, ]
# replace original counts with batch-corrected ones
combat_count <- data.frame(sgRNA = rownames(count), Gene = count[, "Gene"], combat_count)
# write results
write.table(combat_count, file = snakemake@output[[1]], sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
