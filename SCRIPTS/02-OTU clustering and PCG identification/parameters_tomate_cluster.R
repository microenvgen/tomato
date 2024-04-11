setwd("~/micro/tomate/")

# PRIMERS (460, 422 amplicon length)
# 341F	CCTACGGGNGGCWGCAG
# 805R	GACTACHVGGGTATCTAATCC

## Data options
data <- "tomate_23sep"
main = "~/micro/tomate"
setwd(main)

logs_output_f <- paste0(main,"/logs")

fastaF_input_f  <- paste0("/lustre/proyectos/microbioma/resources/Datos/tomate_2022-09/fastQ/R1/") #
fastaR_input_f  <- paste0("/lustre/proyectos/microbioma/resources/Datos/tomate_2022-09/fastQ/R2/") #

fastaF_output_f <- paste0(main,"/",data,"_dada_filter_F_reads/")
fastaR_output_f <- paste0(main,"/",data,"_dada_filter_R_reads/")

dada2cleaned_output_f  <- paste0(main,"/",data,"_dada2cleaned/")
species_assignment_f <- '/home/silviatm/micro/ratones_ugr_2022-06-29/silva_species_assignment_v123.fa.gz' #
taxa_train_set <- '/home/silviatm/micro/ratones_ugr_2022-06-29/silva_nr_v123_train_set.fa' #
  
OTU_table_loc <- paste0(main,"/", data, "_otu_table.csv")
map_table_loc <- paste0(main,"/", data, "_map_table.csv")

## Plot and .RData options
quality_plots  = FALSE # TODO doesn't work for some reason.
plot_errors    = T
bar_plot_top   = T
         top   = 20
plotkrona      = FALSE
read_len_plot  = FALSE # TODO check if it works with 0-read samples

saveallRData   = TRUE
savefinalRData = TRUE

## Run options
multithread <- 4 # dada function loads the whole fastq data to each core;
                 # make sure you have (fastq size) * multithread space available
rarefy <- TRUE
assign_tax_to_ASVs <- FALSE
save_ASVs <- TRUE
  ASV_table_loc     <- paste0(main,"/", data, "_asv_table.csv")
  ASV_tax_table_loc <- paste0(main,"/", data, "_asv_tax_table.csv")

remove_chloroplasts <- TRUE
remove_mitochondria <- TRUE
remove_nonbacteria  <- TRUE

# trimming+merging troubleshooting https://github.com/benjjneb/dada2/issues/1440
# I took a look at the samples and it looks to me like @fanli-gcb suggested: The relatively low quality of the sequences combined with the relatively long amplicon is challenging. High frequency variants are still picked out fine, but low frequency stuff is often not getting identified correctly in both F and R reads, and thus failing to merge.  This is one dataset where I would consider allowing maxMismatch=1 in mergePairs, which got another 10-20% of the reads through. But there are limits to what the method can do with long low-quality amplicons.
# --> allow for not so good quality if merging goes poorly
cutadapt_args1 <- "module load cutadapt; " #any preceeding args to cutadapt
FWD_PRIMER_LEN <- 17
REV_PRIMER_LEN <- 21
ADAPTER_FWD <- 'CCTACGGGNGGCWGCAG' # or just primers
ADAPTER_REV <- 'GACTACHVGGGTATCTAATCC'
cutadapt_args2 <- "--discard-untrimmed --report=minimal "#-m 100" # any trailing args to cutadapt
fastaF_trim_f <- paste0(main,"/",data,"_dada_trim_F_reads/")
fastaR_trim_f <- paste0(main,"/",data,"_dada_trim_R_reads/")

# prepare cutadapt script
cat(paste0("#!/bin/bash\n",
           "cutoutputF=$1\n",
           "cutoutputR=$2\n",
           "fnFs=$3\n",
           "fnRs=$4\n",
           cutadapt_args1, 
           "cutadapt -g ", ADAPTER_FWD, " -G ", ADAPTER_REV," -o $cutoutputF -p $cutoutputR $fnFs $fnRs --action=trim ",
           cutadapt_args2, " ",
           "| tail -1 | cut -f2,7"), # only get the initial and final (= with primers) number of sequences
    file=paste0("cutadapt_",data,".sh"))

truncLen <- c(280-17, 220-21) 

maxmismatches <- 1 
removeBimeraDenovo_method <- "consensus"
allowoneoff <- TRUE

verbose <- TRUE
width   <- 20000 # how many chars before a linebreak when writing final .fasta
