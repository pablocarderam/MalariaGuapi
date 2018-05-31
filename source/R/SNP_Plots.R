# SNP Analysis for vcf data

# Troubleshoot with Pablo at pablocarderam@gmail.com

rm(list = ls()) # Daoist pu

# Imports
library(readr)
library(ggplot2)
library(plyr)
library(gridExtra)
library(stringr)

# Function gets legend, stolen from http://stackoverflow.com/a/12041779/914024
gglegend <- function(x){ 
  tmp <- ggplot_gtable(ggplot_build(x)) 
  leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box") 
  tmp$grobs[[leg]]
}

# Rerun SNP graphs from here:
snp_cov <- read_csv("~/Documents/GitHub/MalariaGuapi/dat/full_snp_cov_1000bp_windows_long_form.csv")

plotChrom <- function(snp_cov_tab,legend_pos="none") {
  chrom = snp_cov_tab$Chromosome[1]
  
  alpha_val=1
  if(legend_pos=="none") {
    alpha_val=0.25
  }
  
  p1 = qplot(x=snp_cov_tab$Position,y=snp_cov_tab$Count, color=snp_cov_tab$Sample, size=I(0.25), geom="point", alpha=I(alpha_val)) + 
    theme(legend.position=legend_pos) + guides(col = guide_legend(title="Sample",ncol = 4))
  p = p1 + labs(list(title=chrom, x= "Base Position", y="Density"))
  return(p)
}

chromPlts = dlply(snp_cov, .(snp_cov$Chromosome), function(x)plotChrom(x))
chromPlts = c(chromPlts[2:16],chromPlts[1]) # reorder plots
leg = gglegend(plotChrom(snp_cov,"bottom"))

snp_plt = do.call("arrangeGrob", c(chromPlts, ncol=2))

ggsave("plt/snp_density_all.png", plot=snp_plt, width =8, height=10, dpi = 600)
ggsave("plt/snp_density_all_legend.png", plot=leg, height=10, dpi = 600)

