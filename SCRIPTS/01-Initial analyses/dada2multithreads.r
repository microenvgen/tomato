

truncLen <- c(268, 218)
maxmismatches <- 1 
removeBimeraDenovo_method <- "consensus"
allowoneoff <- TRUE

# =====LIBRARIES================================================================
r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)


library(dada2)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(doParallel)
library(readxl)
library(DECIPHER)
library(phangorn)


# File parsing
pathF <- "R1/" 
pathR <- "R2/" 
filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered")
fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
trimmed <- filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(268,218), rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)


# File parsing
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
    derepF <- derepFastq(filtFs[[sam]])
    ddF <- dada(derepF, err=errF, multithread=TRUE)
    derepR <- derepFastq(filtRs[[sam]])
    ddR <- dada(derepR, err=errR, multithread=TRUE)
    merger <- mergePairs(ddF, derepF, ddR, derepR,maxMismatch = maxmismatches, minOverlap = 12)
    mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)

save.image("dada2multithreads_intermediate.RData")

# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)

# Remove chimeras
seqtab_nochim <- removeBimeraDenovo(seqtab, method=removeBimeraDenovo_method, multithread=TRUE,allowOneOff = allowoneoff)

# Assign taxonomy
tax <- assignTaxonomy(seqtab_nochim, "silva_nr_v123_train_set.fa.gz", multithread=TRUE)

# Write to disk
#saveRDS(seqtab, "path/to/study/seqtab_final.rds") # CHANGE ME to where you want sequence table saved
#saveRDS(tax, "path/to/study/tax_final.rds") # CHANGE ME ...

getN <- function(x) sum(getUniques(x))
track <- cbind(trimmed, 
               sapply(mergers, getN), 
               rowSums(seqtab_nochim))
colnames(track) <- c('input', 'filtered', 'dada+merged', 'nochim')        
                 
                    
save.image("dada2multithreads.RData")





