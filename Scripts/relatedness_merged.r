#!/usr/bin/Rscript
 
 
args <- commandArgs(trailingOnly = TRUE)
print(args)
 
name=args[1]
 
 
related <- read.table(paste(name, ".genome", sep=''), header = TRUE)
dim(related)
write.table(unique(related[,1:2]), paste("merged_related_IDs", sep=''), append = FALSE, quote = FALSE, row.names = FALSE, col.names = FALSE)
