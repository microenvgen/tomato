#> This script plots a list of KO terms for each node in a list of nodes. Before
#> doing so however, it subsets only those terms present in 90% of the OTUs in
#> the node or (2) EXCLUSIVE to each node and part of the minimal metagenome.

library(dplyr)
library(ggplot2)

# === PARAMETERS ===============================================================
nodes <- c("Node119365", "Node126898", "Node147305",
           "Node15327", "Node19398", "Node21166",
           "Node28120", "Node45985", "Node52227",
           "Node53994", "Node55420",
           "Node55904")

# nota: porque quiero ver el pangenoma exclusivo de cada PCG, pero exclusivo no entre ellos mismos sino en toda la comunidad.
nodes_all20 <- c("Node13271", "Node13287", "Node15327", "Node19398", "Node21166", "Node28120", "Node45985", "Node52227", "Node53994", "Node55420", "Node55904", "Node70899", "Node71253", "Node101358", "Node115605", "Node117321", "Node119365", "Node126898", "Node137058", "Node147305")
  
input_folder <- "/home/silvia/AAA/2023-05-05_Metagenoma_minimo/results_picrust2_PCG_all20/KOs" # ídem nota
save_dir <- "/home/silvia/AAA/2024-01-03_plot_pangenomes"
prefix_pan <- "pan90_"
prefix_mm <- "metagminim_"
input_suffix <- "_core_metag_All.txt"
colorplot1 <- "#048c7f"
colorplot2 <- "#045c9f"
titleplot1 <- "Pangenoma del PCG al 90%"
titleplot2 <- "Contribución del PCG al metagenoma mínimo"
imagename <- "analysisKOfinal.RData"

# ignored bc not informative at all
fields_to_ignore <- c("Metabolic pathways", "Microbial metabolism in diverse environments")

# params for pangenome 90%
threshold <- 0.9
contributions_f <- "/home/silvia/AAA/2023-05-05_Metagenoma_minimo/metagenome_predictions_picrust2_PCG_all20/"

# === READ DBs =================================================================
ko2pathway <- "/home/silvia/mambaforge/lib/python3.10/site-packages/picrust2/default_files/pathway_mapfiles/KEGG_pathways_to_KO.tsv"
lines <- readLines(ko2pathway)
ko2pathway <- c(); namesko2pathway <- c()
for (line in lines) {
  elements <- strsplit(line, "\t")[[1]]
  pw  <- elements[1]
  kos <- elements[-1]
  ko2pathway <- c(ko2pathway, rep(pw, length(kos)))
  namesko2pathway <- c(namesko2pathway, kos)
}
names(ko2pathway) <- namesko2pathway

pathway2description <- "/home/silvia/mambaforge/lib/python3.10/site-packages/picrust2/default_files/description_mapfiles/KEGG_pathways_info.tsv.gz"
lines <- readLines(pathway2description)
pathway2description <- c(); namespathway2description <-c()
for (line in lines) {
  elements <- strsplit(line, "\t")[[1]]
  pw  <- elements[1]
  descs <- elements[-1]
  namespathway2description <- c(namespathway2description, rep(pw, length(descs)))
  pathway2description <- c(pathway2description, descs)
}
names(pathway2description) <- namespathway2description

ko2description <- "/home/silvia/mambaforge/lib/python3.10/site-packages/picrust2/default_files/description_mapfiles/ko_info.tsv.gz"
lines <- readLines(ko2description)
ko2description <- c(); namesko2description <- c()
for (line in lines) {
  elements <- strsplit(line, "\t")[[1]]
  pw  <- elements[-1]
  kos <- elements[1]
  ko2description <- c(ko2description, rep(pw, length(kos)))
  namesko2description <- c(namesko2description, kos)
}
names(ko2description) <- namesko2description

# === FUNCTIONS ================================================================
# Function to read KO identifiers from a file
read_ko_file <- function(node, input_folder) {
  file_path <- paste0(input_folder, "/", node, input_suffix)
  ko_list <- readLines(file_path)
  return(ko_list)
}

# === Pangenomas al 90% y KO exclusivas ========================================
all_pans <- list()
otus <- list()
for (node in nodes_all20) { # ídem nota
  contribs <- read.table(paste0(contributions_f, node, "/KO_metagenome_out/pred_metagenome_contrib.tsv.gz"),
                         header = 1)[c("function.", "taxon")]
  # OTU/KO combinations
  unique_combinations <- contribs %>%
    distinct(function., taxon, .keep_all = TRUE) # reduces from 500k+ to 40k+
  
  # list of OTUs of this PCG
  otus[[node]] <- unique(contribs$taxon)
  
  # number of OTUs in which each KO appears
  function_counts <- contribs %>%
    group_by(function.) %>%
    summarize(unique_taxon_count = n_distinct(taxon))
  
  # this number should be 90% of total (TODO: 20/22, pero he visto alguno con 19/22)
  # 1463 out of 2710 for Node55904
  all_pans[[node]] <- filter(function_counts, unique_taxon_count >= (threshold * length(otus[[node]])))[["function."]]
}

pans90 <- list(); pans90_FULL <- list()
for (node in nodes) {
  kos90 <- all_pans[[node]][!all_pans[[node]] %in% unlist(all_pans[-which(names(all_pans) == node)])]
  pathways90 <- ko2pathway[kos90] # there will be NAs
  descriptions90 <- pathway2description[pathways90]
  not_na <- !is.na(pathways90)
  pans90[[node]] <- data.frame(
    KO = kos90[not_na],
    pathways = pathways90[not_na],
    Description = descriptions90[not_na],
    KODescription = ko2description[kos90[not_na]],
    row.names = NULL
  )
  pans90_FULL[[node]] <- data.frame(
    KO = kos90,#[not_na],
    pathways = pathways90,#[not_na],
    Description = descriptions90,#[not_na],
    KODescription = ko2description[kos90],#[not_na]],
    row.names = NULL
  )
}

plots1 <- list()
for (node in names(pans90)) { # no hago node in nodes por si falta alguno!
  summary_df <- pans90[[node]][order(pans90[[node]]$Description), ]
  summary_df <- summary_df[!(summary_df$Description %in% fields_to_ignore), ]
  if (nrow(summary_df)) {
    # Create a bar plot using ggplot2
    theplot <- ggplot(summary_df,
                      aes(x = `Description`)) + #Description)) +
      geom_bar(width = 0.8,
               fill = colorplot1, # either here or in aes() !
               position = position_dodge(width= 1)) +
      scale_x_discrete(position = "bottom") +
      coord_flip(clip = "off") +  # Horizontal bars
      geom_text(aes(label =  Description),
               color = "black",
               stat="count",
               size = 4,
               position = position_dodge(width= 1), hjust=-0.05) +
      labs(title = titleplot1,
           x = NULL,
           y = NULL) +
      theme(legend.position = "none",
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text = element_text(size = 15),
            plot.title = element_text(hjust = 0),
            panel.background = element_rect("white"),
            plot.margin = margin(0,
                                 6*max(nchar(summary_df$Description)),
                                 0,0)
      ) # Hide legend
    plots1[[node]] <- theplot
  }
}

# adjust size and save
for (node in names(pans90)) { # no hago node in nodes por si falta alguno!
  ggsave(file.path(save_dir, paste0(prefix_pan, node, ".png")), plots1[[node]], width = 7, height = 12)
}

# === Metagenoma mínimo + KOs exclusivos =======================================
# First we put all of them together
all_kos <- list()
for (node in nodes_all20) { # ídem nota
  # Fetch pathway information
  all_kos[[node]]  <- read_ko_file(node, input_folder)
}

pathways_dfs <- list(); pathways_dfs_FULL <- list()
for (node in nodes) { # FIX
  ko_list <- all_kos[[node]][!all_kos[[node]] %in% unlist(all_kos[-which(names(all_kos) == node)])]
  if (length(ko_list)== 0) {
    warning(paste("Node", node, "does not have any exclusive KO terms."))
  } else {
    pathways <- ko2pathway[ko_list] # there will be NAs
    descriptions <- pathway2description[pathways]
    not_na <- !is.na(pathways)
    pathways_dfs[[node]] <- data.frame(
      KO = ko_list[not_na],
      pathways = pathways[not_na],
      Description = descriptions[not_na],
      KODescription = ko2description[ko_list[not_na]],
      row.names = NULL
    )
    not_na <- !is.na(pathways)
    pathways_dfs_FULL[[node]] <- data.frame(
      KO = ko_list,#[not_na],
      pathways = pathways,#[not_na],
      Description = descriptions,#[not_na],
      KODescription = ko2description[ko_list],#[not_na]],
      row.names = NULL
    )
  }
}

plots2 <- list()
for (node in names(pathways_dfs)) { # no hago node in nodes por si falta alguno!
  summary_df <- pathways_dfs[[node]][order(pathways_dfs[[node]]$Description), ]
  summary_df <- summary_df[!(summary_df$Description %in% fields_to_ignore), ]
  if (nrow(summary_df)) {
    # Create a bar plot using ggplot2
    theplot <- ggplot(summary_df,
                      aes(x = `Description`)) + #Description)) +
      geom_bar(width = 0.8,
               fill = colorplot2, # either here or in aes() !
               position = position_dodge(width = 1)) +
      coord_flip(clip = "off") +  # Horizontal bars
      geom_text(aes(label =  Description),
                color = "black",
                stat="count",
                size = 4,
                position = position_dodge(width= 1), hjust=-0.05) +
      labs(title = titleplot2,
           x = NULL,
           y = NULL) +
      theme(legend.position = "none",
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text = element_text(size = 15),
            plot.title = element_text(hjust = 0),
            panel.background = element_rect("white"),
            plot.margin = margin(0,
                                 6*max(nchar(summary_df$Description)),
                                 0,0)
            ) # Hide legend
    plots2[[node]] <- theplot
  }
}

# adjust size and save
for (node in names(pathways_dfs)) { # no hago node in nodes por si falta alguno!
  ggsave(file.path(save_dir, paste0(prefix_mm, node, ".png")), plots2[[node]], width = 7, height = 12)
}

save.image(paste0(save_dir, "/", imagename))










# see tax for otus.
plots_g <- list()
tax_f <- "~/AAA/2023-05-05_Metagenoma_minimo/tomate_bc_0.005/Tree/0.99/table.from_biom_0.99.txt"
tax   <- read.table(tax_f, sep="\t", row.names = 1)[41]
for (node in nodes_all20) {
  taxa_list <- tax[as.character(otus[[node]]),] %>% stringr::str_split("; ")
  genera_counts <- table(sapply(taxa_list, function(x) {
    genus <- ifelse(length(x) >= 6, gsub("g__", "", x[6]), "Unknown") # Extract genus or use "Unknown"
    if (genus == "") genus <- "Unknown" # Handle empty genera
    genus
  }))
  genera_df <- data.frame(Genera = names(genera_counts), Count = as.numeric(genera_counts))
  plots_g[[node]] <- ggplot(genera_df, aes(x = Genera, y = Count, fill = Genera)) +
    geom_bar(stat = "identity", color = "black") +
    labs(title = paste0("Nodo ", stringr::str_sub(node, 5, nchar(node)), ": OTUs de cada género"), x = "Género", y = "Número de OTUs") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.text = element_text(size = 7))
}
pdf("~/borrarpdf.pdf", width = 6, height = 6)
plots_g
dev.off()


if (FALSE) {
  png(paste0("~/Documents/tesis_final/figure/", node, "_con.png"))
  # node = "Node19398"
  # node = "Node28120"
  node = "Node53994"
  # abc<-c("K10013", "K10108", "K10109", "K10110")
  # abc<-pathways_dfs$Node28120 [pathways_dfs$Node28120$Description== "ABC transporters",]$KO
  contribs <- read.table(paste0(contributions_f, node, "/KO_metagenome_out/pred_metagenome_contrib.tsv.gz"),
                         header = 1)[c("function.", "taxon")]
  abc <- pathways_dfs[[node]][pathways_dfs[[node]]$Description=="Cell cycle - Caulobacter",]$KO
  select <- contribs$function. %in% abc
  contribsabc <- contribs[select,]
  
  taxa_list <- tax[as.character(unique(contribsabc$taxon)),] %>% stringr::str_split("; ")
  genera_counts <- table(sapply(taxa_list, function(x) {
    genus <- ifelse(length(x) >= 6, gsub("g__", "", x[6]), "Unknown") # Extract genus or use "Unknown"
    if (genus == "") genus <- "Unknown" # Handle empty genera
    genus
  }))
  genera_df <- data.frame(Genera = names(genera_counts), Count = as.numeric(genera_counts))
  ggplot(genera_df, aes(x = Genera, y = Count, fill = Genera)) +
    geom_bar(stat = "identity", color = "black") +
    labs(title = paste0("Nodo ", stringr::str_sub(node, 5, nchar(node)),
                        ": OTUs de cada género"), x = "Género", y = "Número de OTUs") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  dev.off()
  
  png(paste0("~/Documents/tesis_final/figure/", node, "_sin.png"))
  contribsabc <- contribs[!select,]
  taxa_list <- tax[as.character(unique(contribsabc$taxon)),] %>% stringr::str_split("; ")
  genera_counts <- table(sapply(taxa_list, function(x) {
    genus <- ifelse(length(x) >= 6, gsub("g__", "", x[6]), "Unknown") # Extract genus or use "Unknown"
    if (genus == "") genus <- "Unknown" # Handle empty genera
    genus
  }))
  genera_df <- data.frame(Genera = names(genera_counts), Count = as.numeric(genera_counts))
  ggplot(genera_df, aes(x = Genera, y = Count, fill = Genera)) +
    geom_bar(stat = "identity", color = "black") +
    labs(title = paste0("Nodo ", stringr::str_sub(node, 5, nchar(node)), ": OTUs de cada género"), x = "Género", y = "Número de OTUs") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  dev.off()
}