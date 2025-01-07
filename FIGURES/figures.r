load("figures.RData")

library(phyloseq)
library(tidyverse)
library(gridExtra)
library(cowplot)
# ggplot2 3.5 breaks some packages (https://github.com/jtlandis/ggside/issues/56. cowplot seems to be one of these 

theme_set(theme_bw())

tiff("Figure_1a.tiff", width = 1280, height = 784)
plot_ordination(ps_si, ord.nmds.bray,
                shape="Type",
                color="Soil",
                title="A") +
  # title="Bray NMDS (Stress=0.171)") +
  geom_point(size = 6) +  # Double the point size (default is 3)
  theme(
    plot.title = element_text(size = 30),  # Double the title size
    axis.title = element_text(size = 24),  # Double axis title size
    axis.text = element_text(size = 18),   # Double axis labels size
    legend.text = element_text(size = 18), # Double the legend text size
    legend.title = element_text(size = 20) # Double the legend title size
  )
dev.off()




tiff("Figure_1b.tiff", width = 1280, height = 784)
plot_ordination(ps_si3, ps.dpcoa, color="Soil", shape="Type") +
  ggtitle("B") +
  ggplot2::scale_colour_discrete() +
  geom_point(size = 6) +  # Double the point size (default is 3)
  theme(
    plot.title = element_text(size = 30),  # Double the title size
    axis.title = element_text(size = 24),  # Double axis title size
    axis.text = element_text(size = 18),   # Double axis labels size
    legend.text = element_text(size = 18), # Double the legend text size
    legend.title = element_text(size = 20) # Double the legend title size
  )
dev.off()




tiff("Figure_2.tiff", width = 1280, height = 784)
plot(Plot_by_10_classes) +
  geom_point(size = 3) +  # Double the point size (default is 3)
  theme(
    # plot.title = element_text(size = 30),  # Double the title size
    axis.title = element_text(size = 24),  # Double axis title size
    axis.text = element_text(size = 18),   # Double axis labels size
    legend.text = element_text(size = 18), # Double the legend text size
    legend.key.height = unit(1.2, "cm"), # Separation between legend lines
    legend.title = element_text(size = 20) # Double the legend title size
  )
dev.off()




EC_table<-read.table("../SCRIPTS/01-Initial analyses/pred_metagenome_unstrat_EC.tsv", sep="\t", header=TRUE)
rownames(EC_table)=EC_table$`function`
EC_table=EC_table[,-which(colnames(EC_table)=="function.")]
SampleN_EC <- colnames(EC_table) %>% stringr::str_remove("sample_") %>% as.numeric # old colnames in EC_table correspond to the sample numbers in our phyloseq objects
colnames(EC_table)=sample_data(ps_si)$Muestra[match(SampleN_EC, sample_data(ps_si)$Juntar)]

colnames(sample_data(ps_si))[4] <- "Type"
EC<-phyloseq(otu_table(EC_table, taxa_are_rows=T), sample_data(ps_si))
EC.nmds.bray <- ordinate(EC, method="NMDS", distance="bray")
setEPS()
postscript("Figure_3.eps", width = 16, height = 9) # Adjusted size for EPS (16x9 inches)
# tiff("Figure_3.tiff", width = 1280, height = 784)
plot_ordination(EC,
                EC.nmds.bray,
                # title="Bray NMDS (Stress=0.167)",
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



# --------------

scale = 1.2
{
  # Figure 1 (A)
  plot_A <- plot_ordination(ps_si, ord.nmds.bray,
                            shape = "Type", 
                            color = "Soil") +
    annotate("text", x = Inf, y = Inf, label = "A", 
             hjust = 22.5, vjust = 2, size = 12, fontface = "bold") +
    geom_point(size = 8) + # even bigger points
    theme(
      axis.title = element_text(size = 24), 
      axis.text = element_text(size = 18), 
      legend.position = "none",
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10))
    )
  
  # Figure 2 (B)
  plot_B <- plot_ordination(ps_si3, ps.dpcoa,
                            shape = "Type", 
                            color = "Soil") +
    annotate("text", x = Inf, y = Inf, label = "B", 
             hjust = 21.6, vjust = 2, size = 12, fontface = "bold") +
    geom_point(size = 8) +
    theme(
      axis.title = element_text(size = 24), 
      axis.text = element_text(size = 18), 
      legend.position = "none",
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10))
    )
  
  # common legend
  legend <- get_legend(plot_A + theme(legend.position = "bottom",
                                      legend.text = element_text(size = 18),
                                      legend.title = element_text(size = 20),    
                                      plot.background = element_rect(fill = "white", color = NA), # white bg for legend
                                      panel.background = element_rect(fill = "white", color = NA),
                                      # legend.spacing.x = unit(.5, 'cm'),
                                      # legend.key.width = unit(2, "cm",
                                      legend.text.align = 1
  ))
  
  # plot A + B
  combined_plot <- plot_grid(plot_A, NULL, plot_B,
                             nrow = 1,
                             rel_widths = c(1, 0.05, 1))
  final_plot <- plot_grid(combined_plot, legend, 
                          ncol = 1,
                          rel_heights = c(1, 0.2)) +
    theme(plot.background = element_rect(fill = "white", color = NA), # white bg for legend
          panel.background = element_rect(fill = "white", color = NA))
  ggsave("Figure_1_2_combined.tiff", plot = final_plot, width = 16*scale, height = 8*scale)

}

