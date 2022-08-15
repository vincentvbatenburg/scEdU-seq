#!/bin/bash 
#SBATCH --nodes 1
#SBATCH -t 8:00:00
#SBATCH -c 6

# to run on multiple directories, run:
# find . -name "JvB-*" -exec sbatch 2_trim.sh {} \;

DATADIR=$1

cd $DATADIR


cutadapt -o trimmed.R1.fastq.gz -p trimmed.R2.fastq.gz demultiplexedR1.fastq.gz demultiplexedR2.fastq.gz -m 3 -a "IlluminaSmallAdapterConcatBait=GGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTT" -a "IlluminaIndexAdapter=GGAATTCTCGGGTGCCAAGGAACTCCAGTCACN{6}ATCTCGTATGCCGTCTTCTGCTTG"  -A "IlluminaPairedEndPCRPrimer2.0=AGATCGGAAGAGCGN{6}CAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG;min_overlap=5" -A "universalPrimer=GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT;min_overlap=5" -a  "IlluminaGEX=TTTTTAATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATC;min_overlap=5" -a "IlluminaMultiplexingPCRPrimer=GGAACTCCAGTCACN{6}TCTCGTATGCCGTCTTCTGCTTG;min_overlap=5" -A "Aseq=TGGCACCCGAGAATTCCA" -a "Aseq=TGGCACCCGAGAATTCCA"  -a "illuminaSmallRNAAdapter=TCGTATGCCGTCTTCTGCTTGT"

