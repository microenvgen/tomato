This repository contains scripts used in the analyses of "Coupled phylogenetic and functional enrichment in the tomato rhizosphere microbiome" by Talavera-Marcos S, Gallego R, Chaboy R, Rastrojo A, and Aguirre de Cárcer D.


### Initial analyses

These include primer removal and dada2.

1. `remove_341F_805R.py`. Remove primer sequences 

2. `dada2multithreads.R`. Initial sample processing 
	
	
3. `Tomatos.Rmd`. Process sequences and initial analyses (which requires internally also executing blast.sh and cdpcoa.R)

### OTU clustering and PCG identification

1. `processing_dada2_cluster.R`. Process to final OTU table.

It first creates a FASTA of processed 16S sequences.

Then [QIIME](https://github.com/biocore/qiime/releases/tag/1.9.1) scripts are called to cluster the ASVs into OTUs.

The final OTU table is filtered to remove chromosome, non-bacterial and mitochondrial sequences. At this point we removed low quality samples from the input folder and re-run the script with 40 samples in total.

2. `BC_command.sh`. Includes the [BacterialCore.py](https://git.io/Je5V3) arguments. Internally, uses QIIME as well. Uses the FASTA from the previous script.

### PICRUSt2 and KO analysis

- `Pipeline_PCGs_Picrust2.sh`. Main script. Obtains the metagenome of each PCG separately. Uses PICRUSt2's `picrust2_pipèline.py` and `core.py`. The input files are the BIOM and FASTA files obtained with BacterialCore.

- `Pipeline_Picrust2.sh`. Does the same but for all the community.

- `create_PCG_combinations.R`. Generate N combinations of OTUs, including one per PCG. Ensures the selected OTUs do coexist in the actual samples.

- `create_random_combinations.R`. Generate N combinations of OTUs, without taking into the account the PCGs. Ensures the selected OTUs do coexist in the actual samples.

- `explore_combinations.R`. Create a report.

- `final_evaluation.R`. Student's T-test

- `analysisKO.R`. Obtain KO terms exclusive to each node and part of the minimal metagenome. Also creates plots for KO annotations from PICRUSt2. It also generates Supplementary Table 2.
