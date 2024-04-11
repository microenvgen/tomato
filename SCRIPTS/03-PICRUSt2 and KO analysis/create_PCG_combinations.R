## This script generates random combinations of OTUs present at a given
## community type with multiple samples. Outputs tables for explore_combinations.R
## - All combinations of OTUs are known to coexists in at least one sample.
## - Each OTUs belongs to a different PCG, thus being a random minimal
## representation of the minimal metagenome.
## These combinations are used to subset the input abundance table, obtained
## by PICRUSt's predict_metagenome.py and metagenome_contributions.py. The subsequent
## commands are saved to a file specified as out_script.
## 
## Warning: less than 5% of OTUs have a KO equivalent in PICRUSt's table.
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

## Set if we want to consider an extra non-PCG node ("others")
others <- FALSE
if (others) {suffix <- "test"} else { suffix <- ""}

## Set output folder
output_folder <- paste0("/home/silvia/AAA/2023-05-05_Metagenoma_minimo/results_PCG_evaluation/minim_metag_combis", suffix, "")
if (!file.exists(output_folder)) {system(paste("mkdir -p", output_folder))}

# Set results file
results_file <- paste0("/home/silvia/AAA/2023-05-05_Metagenoma_minimo/results_PCG_evaluation/results_PCG_combis_", suffix, ".csv")

# Set whole core metagenome file (for explore_combinations.R)
metag_file <- "results_picrust2/KOs/core_metag_All.txt"

# File with KO predictions for each OTU
OTU2KOf <- "/home/silvia/AAA/2023-05-05_Metagenoma_minimo/metagenome_predictions_picrust2/KO_predicted.tsv.gz"

# PCG table
pcgtablef <- "results_0005_coherent.txt"

# OTUtoKO
OTU2KOf <- "/home/silvia/AAA/2023-05-05_Metagenoma_minimo/metagenome_predictions_picrust2/KO_predicted.tsv.gz"


# # all 20, coherent of not
# pcgtablef <- "tomate_bc_0.005/Tree/results.txt"
# N_OTUs = 20 # 20 == 20 PCGs TOTALES para cutoff = 0.005
# output_folder <- paste0("/home/silvia/AAA/2023-05-05_Metagenoma_minimo/results_PCG_evaluation_all20/minim_metag_combis", suffix, "")
# results_file <- paste0("/home/silvia/AAA/2023-05-05_Metagenoma_minimo/results_PCG_evaluation_all20/results_PCG_combis_", suffix, ".csv")
# if (!file.exists(output_folder)) {system(paste("mkdir -p", output_folder))}

####FILES#######################################################################
## Load dataset
otu_table <- read.csv("tomate_bc_0.01/Tree/0.99/table.from_biom_0.99.txt", sep = "\t", row.names = 1, skip = 1)
otu_table <- otu_table[1:(ncol(otu_table) - 1)]

## Load OTU-to-KO table
OTU2KO <- read.csv(OTU2KOf, sep ="\t", row.names = 1)

## Initialize metadata vector
metadata <- list()
metadata_file <- paste0("/home/silvia/AAA/2023-05-05_Metagenoma_minimo/results_PCG_evaluation/metadata_PCG", suffix, ".txt")

## Load pcg_table
# pcg_table <- read.csv("tomate_bc_0.01/Tree/results.txt", sep ="\t")[c("Core", "Leaves")]
pcg_table <- read.csv(pcgtablef, sep ="\t")[c("Core", "Leaves")]
pcg_table <- pcg_table[1:(nrow(pcg_table) - 1),]
################################################################################
## Get list of present non-zero abundance OTUs for each sample
present_otus <- apply(otu_table, 2, function(x) {
  p <- names(x[x > MIN_ABUN])
  }
)

## Create a named vector of OTU to core assignments
otu_to_core <- pcg_table %>%
  separate_rows(Leaves, sep = ";")
# otu_to_core <- otu_to_core[otu_to_core$Leaves %in% rownames(OTU2KO), ]

## (!) if we want to consider an extra non-PCG node
if (others) {
  warning(paste("Adding", (length(unique(unlist(present_otus))) - length(otu_to_core$Leaves)), "non-PCG leaves as 'others'"))
  new_leaves <- tibble("Core" = "others",
                       "Leaves" = unique(unlist(present_otus))[!(unique(unlist(present_otus)) %in% otu_to_core$Leaves)])
  otu_to_core <- bind_rows(otu_to_core, new_leaves)
}

## For each sample, we will obtain N random combinations of (coexisting) OTUs
## from each of the different Cores
for (sa in names(present_otus)) {
  ## Select those OTU that coexist in sample sa
  otu_to_core_sa <- otu_to_core[otu_to_core$Leaves %in% present_otus[[sa]], ]
  if (length(unique(otu_to_core_sa$Core)) < length(pcg_table$Core)) {
    warning(paste("Can't obtain combinations for sample", sa, "since there are only", length(otu_to_core_sa), "coexisting OTUs and we need at least", length(pcg_table$Core)))
    next
  }
  ## Split $Core by unique values
  split_OTUs <- split(otu_to_core_sa$Leaves[seq_along(otu_to_core_sa$Core)],
                         otu_to_core_sa$Core)
  
  ## Loop N times (per sample)
  for (n in 1:N) {
    ## Get random combis
    random_combi <- sapply(split_OTUs, function(x) sample(x, 1))
    
    # Save random_combi to metadata
    metadata[[paste0(sa, n)]] <- unname(random_combi)
    
    ## Use random_combi as indeces to subset rows from OTU2KO
    temp_file <- paste0(output_folder, "/temp_", sa, "_", n, ".txt")
    
    write.table(OTU2KO[random_combi, ],
                file = paste0(output_folder, '/', sa, '_', n, '_METAGENOME_TABLE.txt'),
                row.names = T, col.names = T, sep = "\t")
  }     
}

# Save metadata
write.csv(as.data.frame(metadata_file) %>% t, file = "~/test", quote = FALSE, row.names = T)

# Get final results
system(paste0("Rscript explore_combinations.R ", output_folder, " ", results_file, " ", metag_file, " ", getwd()))