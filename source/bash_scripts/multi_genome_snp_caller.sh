#!/bin/bash
#SBATCH -c 8                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-00:10                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=16000                        # Memory total in MB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=pablo_cardenasramirez@hms.harvard.edu   # Email to which notifications will be sent

# Code used in P. falciparum genome analysis for multiple genome cram files
# Pablo Cardenas, pablocarderam@gmail.com

# Requires python, perl, gcc, samtools, awk, fastqc, seqtk, bwa (or bowtie2), \
# bcftools, bedtools  

# Usage: sbatch bash_scripts/multi_genome_snp_caller.sh /n/scratch2/pc182/Test/dat test \
#/n/scratch2/pc182/Test/Pfal3D7/PlasmoDB-37_Pfalciparum3D7_Genome.fasta
# File genome_snp_pipeline.sh must be in subdirectory bash_scripts accessible
# to this script under 

hostname
cd $HOME/projects/MalariaGuapi

# ARGUMENT HANDLING, DIRECTORIES, AND PATH VARIABLES
# First argument is path to bash script, so we skip 0
dat=$1 # 1: Path to folder with input files
out_dir_name=$2 # 2: Name of output folder
# 3: Path to genome fasta file
ref_gen=${3:-$HOME/projects/MalariaGuapi/dat/PlasmoDB-37_Pfalciparum3D7_Genome.fasta} 
# 4: Path to SNP pipeline script
pipeline=${4:-bash_scripts/genome_snp_pipeline.sh} 

# Make output directory, set variable
mkdir $dat/../$out_dir_name
out=$dat/../$out_dir_name 

# Make directory for output files
mkdir $out/outputt
# Make directory for snp call vcf output files
mkdir $out/output/snp_call
# Make directory for snp coverage output files
mkdir $out/output/snp_cov
# Make directory for read coverage files
mkdir $out/output/coverage
# Make directory for quality control files
mkdir $out/output/qc
# Directory for intermediate files
mkdir $out/temp
# Directory for output and error files produced by single jobs
mkdir $out/output/logs

# PREPARE REFERENCE FILES

# Make list of all sample files within directory
ls $dat | grep cram > $out/sample_list.txt

# Generate genome indexes for alignment with bwa using bwtsw algorithm for \
# large genomes
bwa index -a bwtsw $ref_gen
# Bowtie2 version (abandoned in favor of bwa)
# bowtie2-build $dat/PlasmoDB-37_Pfalciparum3D7_Genome.fasta $dat/Pfal_3D7

# Generate fai file to get chromosome sizes
samtools faidx $ref_gen
# generates file with chromosome sizes
cut -f1,2 ${ref_gen}.fai > ${ref_gen}.sizes
# make windows in which to bin SNPs
bedtools makewindows -g ${ref_gen}.sizes -w 1000 > ${ref_gen}.windows.bed

# 3: P > ${refgen}.sizesath to genome fasta file

# Embarrassingly parallel for all samples
for sample_file_name in $(cat $out/sample_list.txt) ; do

    # Call separate bash file with individual genome analysis pipeline
    sbatch -o $out/output/logs/${sample_file_name}.o \
-e $out/output/logs/${sample_file_name}.e \
$pipeline $sample_file_name $dat $out $ref_gen

    echo $sample_file_name

done

