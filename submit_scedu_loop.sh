#!/bin/bash

files=(${@})

while getopts c: flag
do
    case "${flag}" in
        c) config=${OPTARG};;
    esac
done

for file in "${files[@]}"; do
       
       if echo $file | grep -q -- "-c"; then
         continue
       fi
       
        if echo $file | grep -q "config"; then
         continue
       fi
       
       echo $file 
      # echo $config
       
       sbatch -t 50:00:00 --mem 5G --nodes 1 -c 1 -o "$file"_%j.out /hpc/hub_oudenaarden/vincentvb/forkdotV2/scripts/scedu_hpc_analysis.R $file $config
       
       #--gres=tmpspace:100G
       
       done
       
       