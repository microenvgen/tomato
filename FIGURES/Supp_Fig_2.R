# from 'Tomato 1.Rmd'

library(phyloseq)
library(tidyverse)

load("supp_fig_2.rdata")
theme_set(theme_bw())

# normalize
plot(sort(sample_sums(ps_si))) | abline(h=45000)

# figure
nivel<-"Family"
ps_si_nivel <- tax_glom(ps_si, nivel)

Top10nivel = names(sort(taxa_sums(ps_si_nivel), TRUE)[1:10])
ps_si_nivel = prune_taxa(Top10nivel, ps_si_nivel)
# p <- plot_bar(ps_si_nivel, fill=nivel)
# ggsave(plot = p, width = 8, height = 6,
#        filename = "Suppl_Figure_2.eps",
#        device = "eps")
top10_abundance <- sample_sums(ps_si_nivel)
Others <- 45000 - top10_abundance

# Add the "non-top 10" taxon to the phyloseq object
otu_table_extended <- cbind(otu_table(ps_si_nivel),
                            data.frame(Others))
tax_table_extended <- as.data.frame(tax_table(ps_si_nivel))
tax_table_extended["Others", ] <- c(rep("Others", ncol(tax_table_extended)))
ps_si_extended <- phyloseq(
  otu_table(as.matrix(otu_table_extended), taxa_are_rows = F),
  tax_table(as.matrix(tax_table_extended)),
  sample_data(ps_si)
)

# color + plot
color_levels <- c(
  "Bacillaceae",
  "Chitinophagaceae",
  "Comamonadaceae",
  "Cytophagaceae",
  "Flavobacteriaceae",
  "Oxalobacteraceae",
  "Pseudomonadaceae",
  "Rhizobiaceae",
  "Sphingobacteriaceae",
  "Streptomycetaceae",
  "Others"
)
custom_colors <- c(
  # "#D3D3B0",
  scales::hue_pal()(10),
  "#e8e8df"
)
names(custom_colors) <- color_levels

p <- plot_bar(ps_si_extended, fill = nivel) +
  scale_fill_manual(values = custom_colors)
# order families
p$data$Family <- factor(p$data$Family, levels = c("Others", color_levels[1:10]))
p$data$Sample <- factor(p$data$Sample, levels = names(top10_abundance[order(top10_abundance, decreasing = T)]))
ggsave(plot = p, width = 8, height = 6,
       filename = "Suppl_Figure_2.eps",
       device = "eps")
ggsave(plot = p, width = 8, height = 6,
       filename = "Suppl_Figure_2.tiff",
       device = "tiff")
