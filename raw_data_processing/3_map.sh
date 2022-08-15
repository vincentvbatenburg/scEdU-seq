#!/bin/bash 
#SBATCH --nodes 1
#SBATCH -t 8:00:00
#SBATCH -c 6

# to run on multiple directories, run:
# find . -name "JvB-*" -exec sbatch 3_map.sh {} \;

DATADIR=$1

cd $DATADIR

bwa mem -t 6 -M /hpc/hub_oudenaarden/group_references/ensembl/97/homo_sapiens/primary_assembly_NOMASK_ERCC92.fa trimmed.R1.fastq.gz trimmed.R2.fastq.gz | samtools view -b - > ./unsorted.bam; samtools sort -T ./temp_sort -@ 4 ./unsorted.bam > ./sorted.unfinished.bam; mv ./sorted.unfinished.bam ./sorted.bam; samtools index ./sorted.bam; rm ./unsorted.bam;


