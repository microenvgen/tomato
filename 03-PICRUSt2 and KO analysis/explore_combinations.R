# Similarly to core.R, this script parses the KOs present in a KOs abundance table.
# Creates a report.
library(stringr)

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3 || length(args) > 4) {
  stop("Usage: Rscript explore_combinations.R <combi_folder> <results_file> <total_core_kos_file> [working_directory]")
}

# Input + output
combi_folder <- args[1] # "/home/silvia/AAA/2023-05-05_Metagenoma_minimo/results_PCG_evaluation/minim_metag_combis"
results_file <- args[2] # "/home/silvia/AAA/2023-05-05_Metagenoma_minimo/results_PCG_evaluation/results_PCG_combis.csv"
metag_file   <- args[3]

# Set the working directory if provided
if (length(args) == 4) {
  setwd(args[4])
} else {
  setwd("~/AAA/2023-05-05_Metagenoma_mínimo")
}

# Initialize table of results
results <- data.frame(matrix(ncol = 3))
colnames(results) <- c("Sample", "Reactions in minimal metagenome", "% Coverage minimal metagenome")

# Read files
filenames <- list.files(combi_folder,
                         pattern = "*METAGENOME_TABLE.txt")
results_file <- args[2] # "/home/silvia/AAA/2023-05-05_Metagenoma_minimo/results_PCG_evaluation/results_PCG_combis.csv"

# Core in all, add to results table
all_core_kos <- read.table(metag_file, header = T, sep = "\t")
total_kos <- nrow(all_core_kos)
results[1, ] <- c("All", total_kos, 100) # add to result table

# Core covered by combinations
for (file in filenames) {
  # Read present KOs
  combi_name <- paste(str_split(file, "_")[[1]][1:2], collapse = "_")
  combi_abunds <- read.csv(paste0(combi_folder, "/", file), sep= "\t", row.names = 1, header = 1)
  combi_kos <- colnames(combi_abunds[colSums(combi_abunds) > 0])
  
  # Core kos of samples
  kos_in_core   <- sum(combi_kos %in% unlist(all_core_kos))
  results <- rbind(results,
                   c(combi_name, kos_in_core, paste0(100*round((kos_in_core/total_kos), 4))))
}

# Save results and remove temporal files
write.csv(results, results_file, row.names = F)
