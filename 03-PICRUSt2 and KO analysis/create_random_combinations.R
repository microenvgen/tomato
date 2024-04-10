## This script generates random combinations of OTUs present at a given
## community type with multiple samples. Outputs tables for explore_combinations.R
## - All combinations of OTUs are known to coexists in at least one sample.
## - They belong to random PCGs - we do not check for this
## These combinations are used to save a subset of the input abundance table to
## a file that will be used by PICRUSt's predict_metagenome.py. The subsequent
## commands are saved to a file specified as out_script.
## 
## Output:
## - A set of temporal abundance tables, each one for each combination created
## - A metadata file containing the random combinations of indices generated for
##   each iteration. The file is saved in the working directory with the name "metadata_n.txt", where "n" is the number of iterations specified in the script.
## 
## Usage:
##   1. Set the number of iterations per sample (N) and any other variables
##      like input files
##   2. Run the script in the command line
## 
## Notes:
##   With PICRUSt1, we subset the OTUs from a precalculated table for the
##   translation from OTU to KO. This was however limited to a set of 99%
##   identity-clustered OTUs. With PICRUSt2, we use the sequence itself (or a close
##   relative; but not limited), and draw from a precalculated table, but the IDs
##   are different. Se we can't use the PICRUSt2 database directly. But that's
##   why there's a generated file which keeps the OTU -> KO results!
##   - PICRUSt1: OTU2KO is the precalculated file (or a subset)
##   - PICRUSt2: OTU2KO is a more complete results file
################################################################################

## Load necessary packages
library(tidyverse)

## Working dir
setwd("/home/silvia/AAA/2023-05-05_Metagenoma_minimo")

####VARIABLES###################################################################
## Set number of combinations
N = 100

## Set MIN_ABUN
MIN_ABUN <- 0 # only look at OTUs with an abundance higher than MIN_ABUN

## Set number of OTUs to be selected
N_OTUs = 12 # 12 == 12 PCGs coherentes para cutoff = 0.005

## Set output folder
output_folder <- paste0("/home/silvia/AAA/2023-05-05_Metagenoma_minimo/results_PCG_evaluation/minim_metag_combis_random", N_OTUs)
if (!file.exists(output_folder)) {system(paste("mkdir -p", output_folder))}

# Set results file
results_file <- paste0("/home/silvia/AAA/2023-05-05_Metagenoma_minimo/results_PCG_evaluation/results_random_combis_", N_OTUs, ".csv")

# Set whole core metagenome file (for explore_combinations.R)
metag_file <- "results_picrust2/KOs/core_metag_All.txt"

# File with KO predictions for each OTU
OTU2KOf <- "/home/silvia/AAA/2023-05-05_Metagenoma_minimo/metagenome_predictions_picrust2/KO_predicted.tsv.gz"

####FILES#######################################################################
## Load dataset
otu_table <- read.csv("tomate_bc_0.01/Tree/0.99/table.from_biom_0.99.txt", sep = "\t", row.names = 1, skip = 1)
otu_table <- otu_table[1:(ncol(otu_table) - 1)]

## Load OTU-to-KO table
OTU2KO <- read.csv(OTU2KOf, sep ="\t", row.names = 1)

## Initialize metadata vector
metadata <- list()
metadata_file <- paste0("/home/silvia/AAA/2023-05-05_Metagenoma_minimo/results_PCG_evaluation/metadata_random", N_OTUs, ".txt")

# # all 20, coherent of not
# N_OTUs = 20 # 20 == 20 PCGs TOTALES para cutoff = 0.005
# output_folder <- paste0("/home/silvia/AAA/2023-05-05_Metagenoma_minimo/results_PCG_evaluation_all20/minim_metag_combis_random", N_OTUs)
# results_file <- paste0("/home/silvia/AAA/2023-05-05_Metagenoma_minimo/results_PCG_evaluation_all20/results_random_combis_", N_OTUs, ".csv")
# metadata_file <- paste0("/home/silvia/AAA/2023-05-05_Metagenoma_minimo/results_PCG_evaluation_all20/metadata_random", N_OTUs, ".txt")
# if (!file.exists(output_folder)) {system(paste("mkdir -p", output_folder))}

################################################################################
## Get list of present non-zero abundance OTUs for each sample
present_otus <- apply(otu_table, 2, function(x) {
  p <- names(x[x > MIN_ABUN])
}
)

## For each sample, we will obtain N random combinations of (coexisting) OTUs
for (sa in names(present_otus)) {
  ## Loop N times (per sample)
  for (n in 1:N) {
    ## Get random combis
    random_combi <- sample(present_otus[[sa]], N_OTUs)
    
    # Save random_combi to metadata
    metadata[[paste0(sa, n)]] <- c(paste0(sa, n), unname(random_combi))
    
    ## Use random_combi as indeces to subset rows from OTU2KO
    temp_file <- paste0(output_folder, "/temp_", sa, "_", n, ".txt")
    write.table(OTU2KO[random_combi, ],
                file = paste0(output_folder, '/', sa, '_', n, '_METAGENOME_TABLE.txt'),
                row.names = T, col.names = T, sep = "\t")
  }     
}

# Save metadata
write.table(metadata, file = metadata_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Get final results
system(paste0("Rscript explore_combinations.R ", output_folder, " ", results_file, " ", metag_file, " ", getwd()))
