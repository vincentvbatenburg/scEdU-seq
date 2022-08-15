#!/bin/bash 
#SBATCH --nodes 1
#SBATCH -t 8:00:00
#SBATCH -c 1

# to run on multiple directories, run:
# find . -name "JvB-*" -exec sbatch 5a_cuts_tsv.sh {} \;

DATADIR=$1

cd $DATADIR

bamTabulator.py --dedup --noqcfail bl-tagged.bam SM,RX,is_proper_pair,is_read1,reference_name,mapping_quality,XA,reference_start,reference_end,DS,flag,dt,RS,ms |awk '$3 == "True" && $4 == "True" && $6>30 && $7 == "None"'| awk '{print $1"\t"$5"\t"$8"\t"$13"\t"$14}' > "$(pwd | grep -oE '[^/]*$')".tsv; gzip *.tsv


