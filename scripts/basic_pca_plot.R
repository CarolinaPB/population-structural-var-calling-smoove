#!/usr/bin/env Rscript

list.of.packages <- c("optparse", "data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library("optparse")
library(data.table)
library(ggplot2)

option_list = list(
  make_option(c("-e", "--eigenvec"), type="character", default=NULL, 
              help="plink eigenvec file", metavar="character"),
  make_option(c("-v", "--eigenval"), type="character", default=NULL, 
              help="plink eigenval file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output PCA pdf", metavar="character"),
  make_option(c("-s", "--sample_list"), type="character", default=NULL, 
              help="sample list", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

pca <- fread(opt$eigenvec, header=F)
eigenval <- scan(opt$eigenval)

pve <- data.frame(pve = eigenval/sum(eigenval)*100)

samples <- fread(opt$sample_list, col.names = c("sample", "bam", "population"), header = F)

pca_samples <- merge(pca, samples, by.x="V1", by.y="sample")
setnames(pca_samples, old = c("V3","V4"),new = c("PC1", "PC2") )

b <- ggplot(pca_samples, aes(PC1, PC2, shape=population, col=population)) + geom_point()
b <- b + scale_shape_manual(values = LETTERS[1:26])
b <- b + theme_light()
b <- b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

ggsave(filename = opt$output,
plot = b, 
device = "pdf", 
)
