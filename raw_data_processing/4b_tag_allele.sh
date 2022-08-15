#!/bin/bash 
#SBATCH --nodes 1
#SBATCH -t 8:00:00
#SBATCH -c 6

# to run on multiple directories, run:
# find . -name "JvB-*" -exec sbatch 4b_tag_allele.sh {} \;

DATADIR=$1

cd $DATADIR
bamtagmultiome.py sorted.bam -alleles ../RPE-1.phased_monosomy.vcf.gz -allele_samples A,B -blacklist ../blacklist.bed -o tagged_allelic_overhang.bam -method nla_no_overhang --multiprocess -tagthreads 6 -mapfile /hpc/hub_oudenaarden/group_references/ensembl/97/homo_sapiens/simulated_CATG_single_69.mappability.safe.bgzf

