# SNP Analysis for vcf data

# Troubleshoot with Pablo at pablocarderam@gmail.com

rm(list = ls()) # Daoist pu

# Imports
library(readr)
library(plyr)
library(stringr)
library(tidyr)

# Set up database
snp_cov = data.frame(Chromosome=character(), # initialize data frame
                     Start=integer(),
                     End=integer(),
                     Count=integer()
)
files_path = "dat/snp_coverage/"
files = dir(files_path, pattern =".txt") # get all snp count table files

for(i in 1:length(files)) { # for each table file (sample),
  snp_cov_i = read_delim(paste(files_path,files[i],sep=""), "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE) # read it
  colnames(snp_cov_i) = c("Chromosome","Start","End","Count") # name columns
  snp_cov_i$Sample = strsplit(files[i],"_SNP_coverage")[[1]][1] # get this sample's name
  snp_cov_i$Position = (snp_cov_i$Start+snp_cov_i$End)/2 # save window midpoint for graphing
  
  snp_cov = rbind.data.frame(snp_cov,snp_cov_i) # add to main data frame
}

snp_cov$Chromosome = data.frame(str_split_fixed(snp_cov$Chromosome, "_", 3))$X2 # get chromosome numbers
snp_cov$Chromosome = paste("Chromosome", snp_cov$Chromosome) # add word to numbers
snp_cov[snp_cov=="Chromosome M76611"] = "Mitochondrial" # special case mitochondrial chromosome
snp_cov[snp_cov=="Chromosome API"] = "Apicoplast" # special case apicoplast chromosome

write_csv(snp_cov,"dat/full_snp_cov_1000bp_windows_long_form.csv") # save main data frame

# Reshape main data frame for clustering analyses and so on
samples = cbind.data.frame(snp_cov$Sample,snp_cov$Count) # keep only sample and count data
colnames(samples) = c("Sample", "Count") # add column names
samples$Range = paste(snp_cov$Chromosome, paste(snp_cov$Start,snp_cov$End,sep = "-"), sep = ": ")
# save complete positional data as single column
samples=spread(samples, Range, Count) # reshape long to wide form using tidyr package
write_csv(samples,"dat/full_snp_cov_1000bp_windows_wide_form.csv") # save wide data frame

