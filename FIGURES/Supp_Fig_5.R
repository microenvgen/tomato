library(phyloseq)
library(tidyverse)

load("figures.RData")
theme_set(theme_bw())


KO_table<-read.table("../SCRIPTS/01-Initial analyses/pred_metagenome_unstrat_KO.tsv", sep="\t", header=TRUE)
rownames(KO_table)=KO_table$`function`
KO_table=KO_table[,-which(colnames(KO_table)=="function.")]
SampleN_KO <- colnames(KO_table) %>% stringr::str_remove("sample_") %>% as.numeric
colnames(KO_table)=sample_data(ps_si)$Muestra[match(SampleN_KO, sample_data(ps_si)$Juntar)]

KO<-phyloseq(otu_table(KO_table, taxa_are_rows=T),sample_data(ps_si))
KO.nmds.bray <- ordinate(KO, method="NMDS", distance="bray")

setEPS()
postscript("Suppl_Figure_5.eps", width = 16, height = 9) # Adjusted size for EPS (16x9 inches)
# tiff("Figure_3.tiff", width = 1280, height = 784)
plot_ordination(KO,
                KO.nmds.bray,
                color="Soil",
                shape="Type") +
  geom_point(size = 6) +  # Double the point size (default is 3)
  theme(
    plot.title = element_text(size = 30),   # Double the title size
    axis.title = element_text(size = 24),   # Double axis title size
    axis.text = element_text(size = 18),    # Double axis labels size
    legend.text = element_text(size = 18),  # Double the legend text size
    legend.title = element_text(size = 20), # Double the legend title size
    legend.spacing.y = unit(.8, "cm")        # Increase vertical spacing between legend items
  ) +
  guides(color = guide_legend(byrow = TRUE), # Ensure legend items are spaced by row
         shape = guide_legend(byrow = TRUE)) 
dev.off()