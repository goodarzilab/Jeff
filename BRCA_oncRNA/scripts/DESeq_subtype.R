library(DESeq2)
library(EnhancedVolcano)
library(argparse)
library(data.table)
#Script used to for BRCA_smRNA project.

parser = ArgumentParser()
parser$add_argument("-m", "--cm", nargs=1, help="Path to count matrix", type="character")
parser$add_argument("-d", "--coldata", nargs=1, help="Path to colData", type="character")
parser$add_argument("-c1", "--condition1", nargs=1, help="Condition 1", type="character")
parser$add_argument("-c2", "--condition2", nargs=1, help="Condition 2", type="character")
parser$add_argument("-o", "--output", nargs=1, help="Path to output directory", type="character")
args <- parser$parse_args()

#Load Data
cm <- read.table(args$cm, sep=",", as.is=TRUE, row.names=1, header=TRUE)
coldata <- read.table(args$coldata, sep=",", as.is=TRUE, row.names=1, header=TRUE)

dds <- DESeqDataSetFromMatrix(countData = cm, colData = coldata, design = ~ condition)
dds <- DESeq(dds)

dds_res <- as.data.frame(results(dds, contrast=c("condition", args$condition1, args$condition2), tidy=TRUE))
dds_res <- dds_res[!is.na(dds_res$padj),]
norm_counts <- counts(dds, normalized=TRUE)

#Output tables
name <- paste(args$output,args$condition1, "v", args$condition2, sep="") 
res_output <- paste(name, "_dds_res.csv", sep="")
counts_output <- paste(name, "_norm_counts.csv", sep="")
write.table(dds_res, file=res_output, sep=",", col.names=TRUE, quote=F)
write.table(norm_counts, file=counts_output, sep=",", col.names=TRUE, quote=F)

