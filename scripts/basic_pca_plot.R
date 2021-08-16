#!/usr/bin/env Rscript

list.of.packages <- c("optparse", "data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library("optparse")
library(data.table)

option_list = list(
  make_option(c("-e", "--eigenvec"), type="character", default=NULL, 
              help="plink eigenvec file", metavar="character"),
  make_option(c("-v", "--eigenval"), type="character", default=NULL, 
              help="plink eigenval file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output PCA pdf", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

pca <- fread(opt$eigenvec, header=F)
eigenval <- scan(opt$eigenval)

pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

pdf(opt$output)
plot(pca$V3, pca$V4, xlab=paste0("PC1 (", signif(pve$pve[1], 3), "%)"), ylab=paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
dev.off()