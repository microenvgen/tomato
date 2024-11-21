# from 'Tomato 1.Rmd'

library(phyloseq)
library(tidyverse)
library(picante)
library(biomeUtils) # devtools::install_github("RIVM-IIV-Microbiome/biomeUtils")

theme_set(theme_bw())

load("supp_fig_3.rdata")
colnames(sample_data(ps_si3))<-c("Sample", "Soil", "Type", "Ojo", "Juntar")
p_faith <- calculatePD(ps_si3)
df <- data.frame(x=sample_data(p_faith)$PD, Soil=sample_data(p_faith)$Soil, Type=sample_data(p_faith)$Type)
df2 <- data.frame(x=sample_data(p_faith)$SR, Soil=sample_data(p_faith)$Soil, Type=sample_data(p_faith)$Type)

df %>%
  ggplot (aes (x = Soil, y = x)) + 
  geom_boxplot() +
  geom_point (aes(color = Type)) + 
  ggtitle("A") +
  theme (axis.title.y = element_blank()) -> p

df2  %>%
  ggplot (aes (x = Soil, y = x)) + 
  geom_boxplot() +
  geom_point (aes(color = Type)) + 
  ggtitle("B") +
  theme (axis.title.y = element_blank())-> q

final <- p + q + patchwork::plot_layout(guides = "collect")

ggsave(plot = final, width = 8, height = 6,
       filename = "Suppl_Figure_3.eps",
       device = "eps")
ggsave(plot = final, width = 8, height = 6,
       filename = "Suppl_Figure_3.jpeg",
       device = "jpeg")
# ggsave(plot = final, width = 8, height = 6,
#        filename = "Suppl_Figure_3.tiff",
#        device = "tiff")