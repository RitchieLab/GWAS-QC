#!/usr/bin/Rscript

if (!require(tidyverse)) install.packages('tidyverse')
library(tidyverse)
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)
 
args <- commandArgs(trailingOnly = TRUE)
if (length(args)!=2) {
  stop("Two arguments must be supplied (input file and output prefix)", call.=FALSE)
}

print(args)

file=args[1]
 
 
het <- read.csv(paste(file, "_het.het", sep=""), sep = "")
miss <- read.csv(paste(file, "_miss.imiss", sep=""), sep = "")

################################################################################
### If you want to run this yourself to change parameter
#setwd("~/Desktop/GWAS/")

#het <- read.csv("PMBB-Release-2020-2.0_genetic_genotype_het.het", sep="")
#miss <- read.csv("PMBB-Release-2020-2.0_genetic_genotype_miss.imiss", sep="")
################################################################################

x <- miss %>% select(IID, F_MISS)
y <- het %>% mutate(HR = (`N.NM.`-`O.HOM.`)/`N.NM.`) %>% select(IID,HR)
to_plot <- inner_join(x,y, by = "IID")

# Calculate mean and sd of HR
HR_mean=mean(y$HR)
HR_SD=sd(y$HR)

# Calculate mean of Missingness
MISS_mean=mean(x$F_MISS)

rbPal_1 <- colorRampPalette(c('#4292c6','#6baed6'))
rbPal_2 <- colorRampPalette(c('#6baed6','#deebf7'))

to_plot$Col1 <- rbPal_1(100)[as.numeric(cut(to_plot$F_MISS,breaks = 100))]
to_plot$Col2 <- rbPal_2(100)[as.numeric(cut(to_plot$F_MISS,breaks = 100))]

pdf(paste(file, "_HM_plot.pdf", sep=""))
plot(to_plot$F_MISS,to_plot$HR, 
     #col = to_plot$Col,
     col = ifelse(to_plot$HR <= HR_mean-(3*HR_SD), "#deebf7", 
                  ifelse(to_plot$HR >= HR_mean+(3*HR_SD), "#deebf7", to_plot$Col2)),
     las = 1, 
     #xlim=c(0, 0.12), # change parameters if you want to adjust the x-axis limits
     #ylim=c(0.1, 0.25), # change parameters if you want to adjust the y-axis limits
     xlab = "Proportion of Missing Genotypes",
     ylab = "Heterozygosity Rate", 
     #main = "NAME OF YOUR PLOT", # Uncomment if you want to include a title for your plot
     pch = 19,
     cex = 0.5,
     log = "x")
    # These commands make the threshold lines for HR and Missingness 
    # according to the H3ABioNet tutorial
    abline(h=HR_mean+(3*HR_SD),col=2,lty=2, lwd = 1)
    abline(h=HR_mean-(3*HR_SD),col=2,lty=2, lwd = 1)
    abline(v=MISS_mean,col=2,lty=2, lwd = 1)
    abline(v=MISS_mean-(0.05*MISS_mean),col=1,lty=1, lwd = 1.25)
    abline(v=MISS_mean+(0.05*MISS_mean),col=4,lty=1, lwd = 1.25)
dev.off()

print("Finished making:", quote = FALSE)
print(paste(file, "_HM_plot.pdf", sep=""), quote = FALSE)
