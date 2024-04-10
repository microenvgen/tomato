# =====PARAMETERS===============================================================
# source("~/AAA/2022-09-13_tomate/parameters_tomate.R")
source("/home/silviatm/micro/tomate/parameters_SOIL_cluster.R")
# ==============================================================================

# =====LIBRARIES================================================================
r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)

library('dada2')
library('logr')
library('phyloseq')
library('ggplot2')
library("ShortRead")
library("stringr")
library("optparse")
library("data.table")
require("gsubfn")
require("tidyverse")
# ==============================================================================

# =====FUNCTIONS================================================================
my_transpose <- function(df){
  t_df=data.table::transpose(df)
  colnames(t_df) <- rownames(df)
  rownames(t_df) <- colnames(df)
  return(t_df)
}
get_abundance_and_tax_from_table <- function (exp_file, # required
                                              species_are_rows=TRUE,
                                              sep="\t", # extra internal options you can change
                                              row.names=1,
                                              skip=1,
                                              taxa_fields=NULL,
                                              tax_col_name="taxonomy",
                                              tax_sep=";",# 4 chars or longer "breaks" renamer (will include this as if it were a real name)
                                              NA_option="___",
                                              check.names=FALSE) { 
  # This function reads a .csv/.tsv speciesXsamples file, where the last "sample" 
  # column (or row) is the taxonomy. Then returns the abundances and taxa data in
  # separate data.frames
  require(gsubfn)
  require(tidyverse)
  if (is.null(taxa_fields)) {taxa_fields = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")}
  
  exp=read.csv(exp_file,sep = sep,skip = skip,row.names=row.names,check.names = check.names)
  if (!species_are_rows) {exp<-my_transpose(exp)}
  
  tax<-exp["taxonomy"]; exp<-exp[1:dim(exp)[2]-1]
  tax <- tax %>% separate("taxonomy",sep = tax_sep,taxa_fields)
  tax[is.na(tax)]<- NA_option # avoid na-related errors
  
  return(list(exp,tax))
}
# ==============================================================================

## Parse arguments
## ---------------
option_list <- list(
  make_option(c("-p", "--params"), type="character", default=NULL,
              help="Parameters file, see examples", metavar="character")
              )
parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

source(opt$params)

## Start log
## ---------
log_open(file_name = paste0(logs_output_f,"/",data),
         logdir = FALSE,
         autolog = TRUE)

log_print(paste0("Data: ",data,".\nCreating necessary directories..."))
for (f in c(logs_output_f, dada2cleaned_output_f, fastaF_trim_f, fastaR_trim_f)) {
  if (!file.exists(f)) {system(paste("mkdir -p", f))}
}

## Initial data exploration
## ------------------------
## Quality plots
if (quality_plots) {
  log_print("Creating quality plots...")
  # > F
  pdf(paste0(logs_output_f, "/", data, "_qualityProfile_F.pdf")) # fixed
  p<-plotQualityProfile(fastaF_input_f)
  plot(p)
  dev.off()
  # > R
  pdf(paste0(logs_output_f, "/", data, "_qualityProfile_R.pdf"))
  p<-plotQualityProfile(fastaR_input_f)
  plot(p)
  dev.off()
}


## Filter and trim
## ---------------
log_print("Removing primers...")

# prepare input and output files
fnFs <- sort(list.files(fastaF_input_f, full.names = TRUE))
fnRs <- sort(list.files(fastaR_input_f, full.names = TRUE))
X <- my_transpose(data.frame(
  cutoutputF <- str_replace(fnFs, fastaF_input_f, fastaF_trim_f),
  cutoutputR <- str_replace(fnRs, fastaR_input_f, fastaR_trim_f),
  fnFs <- fnFs,
  fnRs <- fnRs)
)

# call cutadapt
trimoutput <-
  mcmapply(X=X, 
    function(X) {
       print("Running the following")
       print(paste0("bash ", "cutadapt_",data,".sh ",
                X[1], " ", X[2], " ", X[3], " ", X[4]))
       system(command = paste0("bash ", "cutadapt_",data,".sh ",
                X[1], " ", X[2], " ", X[3], " ", X[4]), intern = T)
    }, mc.cores = multithread
  ) %>% data.frame %>% separate(., col=".", into= c("input.seq", "noprimers"), sep="\t")

log_print(paste0("Trimming to ",truncLen, "..."))
trimmed <- filterAndTrim(fastaF_trim_f, fastaF_output_f,
                         fastaR_trim_f, fastaR_output_f,
                         truncLen = truncLen,
                         verbose=verbose,
                         multithread=multithread)
log_print(head(trimmed))

if (saveallRData){
  save.image("trimmed.RData")
}

## Learn the error rates
## ---------------------
log_print("Learning the error rates...")
errF <- learnErrors(fastaF_output_f, multithread = multithread, verbose = verbose)
errR <- learnErrors(fastaR_output_f, multithread = multithread, verbose = verbose)

log_print("Convergence for F and R error learning:")
log_print(dada2:::checkConvergence(errF)) 
log_print(dada2:::checkConvergence(errR))
# (Even if it doesn't converge, it's ok if the value is way smaller than 
# at the beggining. Here's an example:)
# [1] 3.252136e+01 4.968894e-01 4.556446e-02 2.076082e-02 3.130193e-03
# [6] 3.386694e-03 3.383264e-03 4.794636e-03 2.790897e-03 1.037385e-05

if (saveallRData){
  save.image("after_learnErrors.RData")
}

## Plot Errors
## -----------
if (plot_errors) {
  log_print("Creating error ratio plots...")
  for (i in 1:2) {
    err <- list(errF, errR)[[i]]
    pdf(file = paste0(logs_output_f, "/", data, "_errorPlot_", c("F", "R")[i]))
    p1<-plotErrors(err, nominalQ = TRUE)
    plot(p1)
    dev.off()
  }
  log_print("Done")
}

## Implement the main algorithm with dada()
## ----------------------------------------
log_print("RUNNING DADA...")
dadaF <- dada(fastaF_output_f, err = errF, multithread = multithread)
dadaR <- dada(fastaR_output_f, err = errR, multithread = multithread)
log_print(dadaF[[1]])
log_print(dadaR[[1]])
log_print("Done")

if (saveallRData){
  save.image("after_dada.RData")
}

## Merge forward and reverse reads
## -------------------------------
log_print("Merging forward and reverse reads...")
mergeFR <- mergePairs(dadaF, fastaF_output_f, dadaR, fastaR_output_f, verbose = verbose,
                      maxMismatch = maxmismatches, minOverlap = 12)
log_print(head(mergeFR[[1]]))
log_print("Done")

## Merged read length distribution
if (read_len_plot) {
  pdf(file = paste0(logs_output_f, "/", data, "_readlenPlot"))
  lens <- c()
  for (name in names(mergeFR)) {
    lens <- c(lens, as.numeric(lapply(mergeFR[[name]][[1]], nchar)))
  }
  hist(lens,breaks = 100)
  dev.off()
}

if (saveallRData){
  save.image("after_merging.RData")
}

## Sequence table
log_print("Creating sequence table...")
seq_table <- makeSequenceTable(mergeFR)
log_print(dim(seq_table))
log_print("Done")

## Remove chimeras
log_print("Removing chimeras...")
seq_table_nochim <- removeBimeraDenovo(seq_table, 
                                       method = removeBimeraDenovo_method, 
                                       multithread = multithread, 
                                       verbose = verbose,
                                       allowOneOff = allowoneoff)
log_print(dim(seq_table_nochim))

if (saveallRData){
  save.image("after_chim_removal.RData")
}

getN <- function(x) sum(getUniques(x))
track <- cbind(trimoutput,
               trimmed[,2,drop=F], 
               sapply(dadaF, getN)[rownames(trimmed)],  # does not work if it's only one sample; do getN(dadaF/dadaR/mergeFR)
               sapply(dadaR, getN)[str_replace(rownames(trimmed), "_R1_", "_R2_")], 
               sapply(mergeFR, getN)[rownames(trimmed)], 
               rowSums(seq_table_nochim)[rownames(trimmed)])
colnames(track) <- c('input', "no.primers", 'filtered', 'dadaF', 'dadaR', 'merged', 'nochim')
log_print("Done")


if (assign_tax_to_ASVs) {
  ## Taxonomy assignation
  log_print("Assigning taxons to ASVs...")
  taxaASV <- assignTaxonomy(seq_table_nochim, taxa_train_set, multithread = multithread)
  taxaASV <- addSpecies(taxaASV, species_assignment_f)
  taxa_print <- taxaASV
  rownames(taxa_print) <- NULL
  
  log_print(head(taxa_print))
  log_print("Done")
}

if (save_ASVs) {
  write.csv(seq_table_nochim, file = ASV_table_loc, quote = FALSE)
  if (assign_tax_to_ASVs) {
    write.csv(taxaASV, file = ASV_tax_table_loc)
  }
  log_print("Saved abundance and taxonomy tables for unclustered ASV")
}

## dada2 to fasta (re-replication, necessary for OTU assignment, before rarification)
log_print("Creating .fa files...")
for (i in seq(nrow(seq_table_nochim))) {  
  ids <- paste0(row.names(seq_table_nochim)[i],
                '_',
                seq(rowSums(seq_table_nochim)[i]))  
  
  seqs.re_replicated <- rep(colnames(seq_table_nochim),times=seq_table_nochim[i,]) 
  
  writeFasta(object=ShortRead(sread=DNAStringSet(seqs.re_replicated),
                              id=BStringSet(ids)),
             file=paste0(dada2cleaned_output_f,"/",
                         str_remove(row.names(seq_table_nochim)[i],
                                    ".fastq.gz"),"_dada2cleaned",".fasta"),
             width=width)
}
ASV_fa <- paste0(dada2cleaned_output_f,'/',data,'_ASV.fa')
system(paste0('cat ',dada2cleaned_output_f,'/* > ', ASV_fa))
system(paste0("rm ", dada2cleaned_output_f, "/*_dada2cleaned*"))
log_print("Done")
save.image("before_qiime.RData") # DEBUG

# cat(... = paste0("#!/bin/bash
# # >>> conda initialize >>>
# # !! Contents within this block are managed by 'conda init' !!
# __conda_setup=\"$('/home/silvia/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)\"
# if [ $? -eq 0 ]; then
#     eval \"$__conda_setup\"
# else
#     if [ -f \"/home/silvia/anaconda3/etc/profile.d/conda.sh\" ]; then
#         . \"/home/silvia/anaconda3/etc/profile.d/conda.sh\"
#     else
#         export PATH=\"/home/silvia/anaconda3/bin:$PATH\"
#     fi
# fi
# unset __conda_setup
# # <<< conda initialize <<<
# conda activate qiime1

cat(... = paste0("#!/bin/bash
module load qiime/1.9.1
INPUT=",ASV_fa,"
NAME=$(basename -s '_ASV.fa' $INPUT)
mkdir -p ./qiime

#clustering
pick_otus.py -i $INPUT -o ./qiime -m usearch61_ref -C -z -s 0.99 -x -v --threads ", cores, "
pick_rep_set.py -i ./qiime/$NAME'_ASV_otus.txt' -f $INPUT -o ./qiime/rep_set.fasta
# tax table
assign_taxonomy.py -i ./qiime/rep_set.fasta -o ./qiime/rep_set.assignedTax -m rdp

#otu table
make_otu_table.py -i qiime/$NAME'_ASV_otus.txt' -t qiime/rep_set.assignedTax/rep_set_tax_assignments.txt -o qiime/otu_table_$NAME.biom
biom convert -i qiime/otu_table_$NAME.biom -o qiime/table.from_biom_$NAME.txt --to-tsv --header-key taxonomy
cp qiime/table.from_biom_$NAME.txt ", OTU_table_loc),
file="temp.sh")
system('sh temp.sh')
system("rm temp.sh")


## Phyloseq object creation (for rarefication)
log_print("Loading the new abundance+tax table...")
list[ot, tt]<- get_abundance_and_tax_from_table(exp_file = OTU_table_loc,
                                                species_are_rows = TRUE,
                                                sep = "\t",
                                                check.names=FALSE)
 
newn <- unlist(lapply(rownames(trimmed), FUN=function(x){str_split(x[1], "_")[[1]][1]}))
track <- cbind(track, "OTU.clustering" = rowSums(my_transpose(ot)[newn,]))

## Remove mitoc and chlor OTUs
print("Without removing chloro and mitoc OTUs:")
print(dim(ot))
if (remove_chloroplasts) {
  print("Removing chloroplast OTUs...")
  is.chloro <- tt[,"Class"] %in% " c__Chloroplast"
  ot <- ot[!is.chloro,]
  tt <- tt[!is.chloro,]
  print("After removing chloroplast OTUs:")
  print(dim(ot))
  track <- cbind(track, "no.chloro" = rowSums(my_transpose(ot)[newn,]))
}


if (remove_mitochondria) {
  print("Removing mitochondria OTUs...")
  is.mito <- tt[,"Family"] %in% " f__mitochondria"
  ot <- ot[!is.mito,]
  tt <- tt[!is.mito,]
  print("After removing mitochondria OTUs:")
  print(dim(ot))
  track <- cbind(track, "no.mito" = rowSums(my_transpose(ot)[newn,]))
}

if (remove_nonbacteria) {
  print("Removing non-bacterial OTUs...")
  is.bact <- tt[,"Kingdom"] %in% "k__Bacteria"
  ot <- ot[is.bact,]
  tt <- tt[is.bact,]
  print("After removing non-bacterial OTUs:")
  print(dim(ot))
  track <- cbind(track, "only.bact" = rowSums(my_transpose(ot)[newn,]))
}

log_print(track)

log_print("Creating phyloseq object...")
ps <- phyloseq(otu_table(my_transpose(ot),
                         taxa_are_rows = FALSE),
               tax_table(as(tt, "matrix"),
                         errorIfNULL = FALSE))
sn <- sample_names(ps)

## rarefy
if (rarefy) {
  log_print("Rarefying...")
  pruned <- prune_samples(sample_sums(ps) > 0, ps) #remove those samples with zero reads (NOTE: this changes rownames starting with a number n to Xn)
  ps_final <- rarefy_even_depth(pruned) # TODO more params to tweak here. to the smallest number of seqs by default
} else {
  pruned <- prune_samples(sample_sums(ps) > 0, ps) #remove those samples with zero reads
  ps_final <- pruned
}


## Bar plot top X
if (bar_plot_top) {
  if (is.null(top)) {
    top <- length(taxa_names(ps_final) )
    log_print("Creating bar plot for all taxa...")
  } else {
    log_print(paste0("Creating bar plot for top ", top, " taxa..."))
  }
  toptop <- names(sort(taxa_sums(ps_final), decreasing = TRUE))[1:top]
  ps_top <- transform_sample_counts(ps_final, function(OTU) OTU/sum(OTU))
  pdf(file = paste0(logs_output_f, "/", data, "_phyloseqBarPlot_"), width = 20)
  p1<-plot_bar(ps_top, fill = 'Order') + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p2<-plot_bar(ps_top, fill = 'Family') + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p3<-plot_bar(ps_top, fill = 'Genus') + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  plot(p1); plot(p2); plot(p3)
  dev.off()
  log_print("Done")
}

## Plot Krona
if (plotkrona) {
  log_print("Creating Krona plot...")
  plot_krona(ps_final,'GP-krona','Type')
  log_print("Done")
}

log_print(paste0("Proccessing finished correctly for ",data,"."))

if (savefinalRData) {
  save.image(paste0("FINAL",data,".RData"))
}

## Save RAREFIED, WITH NO CLOROPLASTS/MITOCONDRIA files
if (rarefy || remove_chloroplasts || remove_mitochondria) {
  OTU1 <- my_transpose(as.data.frame(otu_table(ps_final)))
  TAX1 <- as.data.frame(tax_table(ps_final))
  OTU1$taxonomy <- do.call(paste, c(TAX1[1:ncol(TAX1)], sep=";"))
  write.table(OTU1, file = OTU_table_loc, quote = FALSE, sep = "\t", col.names=NA)
}
