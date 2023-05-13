#!/bin/bash

sbatch -t 1:00:00 --mem 20G --nodes 1 -c 2 -o plot_$1_%j.out /hpc/hub_oudenaarden/vincentvb/forkdotV2/scripts/scedu_plotter.R $@

#--gres=tmpspace:100G