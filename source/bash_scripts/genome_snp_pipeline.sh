#!/bin/bash
#SBATCH -c 8                               # Request 8 cores
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-05:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=24000                        # Memory total in MB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE                   # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=pablo_cardenasramirez@hms.harvard.edu   # Email to which notifications will be sent

# Code used in P. falciparum genome analysis for a single genome cram file
# Pablo Cardenas, pablocarderam@gmail.com

# Requires python, perl, gcc, samtools, awk, fastqc, seqtk, bwa (or bowtie2), \
# bcftools, bedtools

# Called from file multi_genome_snp_caller.sh

hostname

# ARGUMENTS
# First argument is path to bash script, so we skip 0
sample_file_name=$1 # 1: Path to cram file to be processed
dat=$2 # 2: Path to folder with input files
out=$3 # 3: Path to output folder
ref_gen=$4 # 4 Path to reference genome fasta file

sample="${sample_file_name%%.*}" # Sample name without extension

# save number of cores
cores=8

# ESTIMATE COVERAGE
# Get coverage of first 2 million bp of each file
samtools depth -a $dat/$sample_file_name | head -n 2000000 | \
awk '{sum+=$3} END {print sum/NR}' > \
$out/output/coverage/${sample}_coverage.txt

echo Coverage estimated.

# GET PAIR END FILES, CONVERT TO SORTED FASTQ
# This doesn't work, for some reason samtools's sorting wasn't doingits job here.
# Switched to converting to FastQ first and sorting after that, as shown below
# in the code block after this one
# 
# Convert cram to bam, cores core parallelization
# samtools view -@ $cores -b $dat/$sample_file_name > $out/temp/${sample}.bam
# Sort reads
# samtools sort -@ $cores -n $out/temp/${sample}.bam -o $out/temp/${sample}_sort.bam
# Convert bam to fastq split by pair ends
# bedtools bamtofastq -i $out/temp/${sample}_sort.bam \
# -fq $out/temp/${sample}_end1.fq -fq2 $out/temp/${sample}_end2.fq    
# 
# Interleave pair end reads into organized single zipped file
# seqtk mergepe $out/temp/${sample}_end1.fq $out/temp/${sample}_end2.fq | \
# gzip > $out/temp/${sample}.fq.gz

# This code does work:
# converts to bam | converts to fastq | sorts according to seq ID as per 
# https://edwards.sdsu.edu/research/sorting-fastq-files-by-their-sequence-identifiers/
samtools fastq -@ $cores $dat/$sample_file_name | \
cat | paste - - - - | sort -k1,1 -t " " | \
tr "\t" "\n" | gzip > $out/temp/${sample}.fq.gz

echo Got pair-end files.

# QUALITY CONTROL with seqtk
# Controversial (especially for RNA-seq) but still considered a good idea. 
# https://doi.org/10.1371/journal.pone.0085024

# FastQC Quality Analysis report
fastqc $out/temp/${sample}.fq.gz -o $out/output/qc/

# Fixed length (redundant with default phred below)
# seqtk trimfq -b 10 -e 10 $out/temp/${sample}_end1.fq > \
#$out/temp/${sample}_trim.fq
# echo Trimmed.

# Quality control and trimming using modified Mott algorithm (phred, 
# http://www.phrap.org/phredphrap/phred.html), deletes read if shorter than 
# 50 bp after QC, then
# | drop pair-end reads if one of pair lost
seqtk trimfq -l 50 $out/temp/${sample}.fq.gz | \
seqtk dropse > $out/temp/${sample}_qc.fq
# Separate pair-end reads again
seqtk seq -1 $out/temp/${sample}_qc.fq > $out/temp/${sample}_qc_end1.fq 
seqtk seq -2 $out/temp/${sample}_qc.fq > $out/temp/${sample}_qc_end2.fq

# FastQC Quality Analysis report after QC
fastqc $out/temp/${sample}_qc.fq -o $out/output/qc/

echo Quality controlled.

# ALIGNMENT 

# Using bowtie2 (abandoned in favor of BWA for P. falciparum, 
# https://doi.org/10.1016/j.ygeno.2017.03.001)

# Global alignment
# bowtie2 -x $dat/Pfal_3D7 -U $dat/${sample}_qc.fq -S $out/temp/${sample}.sam
# echo Global alignment done.

# Local alignment
# bowtie2 --local -x $dat/Pfal_3D7 -U $out/temp/${sample}.fq \
#-S $out/temp/${sample}_local.sam
# echo Local alignment done.
# samtools view -bS $out/temp/${sample}_local.sam > \
#$out/temp/${sample}_local.bam
# echo Local alignment bammed.
# samtools sort $out/temp/${sample}_local.bam \
#-o $out/temp/${sample}_local.sort.bam
# echo Local alignment sorted.

# Using bwa mem, 8 cores
bwa mem -t $cores $ref_gen $out/temp/${sample}_qc_end1.\
fq $out/temp/${sample}_qc_end2.fq > $out/temp/${sample}_global.sam

echo Global alignment done.

# sam to bam and sorting using samtools
samtools view -@ 8 -bS $out/temp/${sample}_global.sam > \
$out/temp/${sample}_global.bam
echo Global alignment bammed.
samtools sort -@ 8 $out/temp/${sample}_global.bam \
-o $out/temp/${sample}_global_sort.bam
echo Global alignment sorted.

# Create bcf file
samtools mpileup -uf $ref_gen $out/temp/${sample}_global_sort.bam | \
bcftools view -Ov - > $out/temp/${sample}_global.raw.bcf
echo Global alignment bcf created

# SNP Calling
bcftools call --skip-variants indels --multiallelic-caller --variants-only \
-O v $out/temp/${sample}_global.raw.bcf \
-o $out/output/snp_call/${sample}_global_snpcall.vcf

# generate txt file with SNP counts per window
bedtools coverage -a ${ref_gen}.windows.bed \
-b $out/output/snp_call/${sample}_global_snpcall.vcf \
-counts > $out/output/snp_coverage/${sample}_SNP_coverage.txt

echo End of SNP analysis.
