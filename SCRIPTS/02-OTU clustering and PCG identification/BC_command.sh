#!/bin/bash

# Download BacterialCore at https://git.io/Je5V3
# This command runs BacterialCore with a cutoff of 0.5%
# This call can substitute the QIIME call in processing_dada2_cluster.R
python BacterialCore.py -f /home/silviatm/micro/tomate_subset.fa -o /home/silviatm/micro/bc/tomate_bc_0.005/ -p 3 -initial_level 0.99 -min_core_percentage 1 -cutoff 0.005 -taxo_p 0.9 -tree_level 99 -tree_type_analysis 3 -threads 16
