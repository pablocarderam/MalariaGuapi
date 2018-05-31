# SNP Analysis for vcf data

# Troubleshoot with Pablo at pablocarderam@gmail.com
# Based off of Boqiang Hu's code at http://huboqiang.cn/2016/03/03/RscatterPlotPCA

rm(list = ls()) # Daoist pu

# Imports
library(readr)
library(ggplot2)
library(gridExtra)

# Get data
samples <- read_csv("~/Documents/GitHub/MalariaGuapi/dat/full_snp_cov_1000bp_windows_wide_form.csv")

row.names(samples) = samples$Sample # row names needed for pca
samples$Sample = NULL

samples_pca = prcomp(samples) # PCA

# Get output data
dat_out = as.data.frame(samples_pca$x)
# Add sample "type"
dat_out$Sample_type = sapply( strsplit(as.character(row.names(samples)), split = "(?<=[a-zA-Z])\\s*(?=[0-9])", perl = TRUE), "[[", 1 )

# Get percentage of each component
percentage = round(samples_pca$sdev / sum(samples_pca$sdev) * 100, 2)
percentage = paste( colnames(dat_out), "(", paste( as.character(percentage), "%", ")", sep="") )

# Plot
pca_plt = ggplot(dat_out,aes(x=PC1,y=PC2,color=Sample_type)) + 
  geom_point(alpha = 0.5, size = 0.75) + 
  xlab(percentage[1]) + ylab(percentage[2]) +
  ggtitle("PCA based on genomic SNP counts")

# Plot contributions
dat_out_r = as.data.frame(samples_pca$rotation)
dat_out_r$feature = row.names(dat_out_r)

cont_plt = ggplot(dat_out_r,aes(x=PC1,y=PC2,label=feature,color="red" )) + 
  geom_point(alpha = 0.5, size = 0.75) + geom_text(size=1) + theme(legend.position="none") +
  ggtitle("Feature contribution to each component")

pca_plts = arrangeGrob(pca_plt,cont_plt, ncol=2)
ggsave("plt/snp_pca.png", plot=pca_plts, width =8, height=3, dpi = 600)
 