library(mapSpain)
library(ggrepel)
library(tidyverse)

# projection used: EPSG:3857 (https://epsg.io/3857)

set.seed(123)

esp_can <- esp_get_country() 
can_prov <- esp_get_can_provinces()
can_box <- esp_get_can_box() 
munic <- esp_get_munic() 

muestra = c("A", "B", "C", "D", "E", "F", "G", "H")
origen = c("Guadalajara (Guadalajara)", "Fuencarral (Madrid)", "Valdelatas (Madrid)", "Vigo (Pontevedra)", "Ruiseñada (Cantabria)", "Moraira (Alicante)", "Lobón (Badajoz)", "Algete (Madrid)")
Location = paste(muestra, origen, sep = " - ")

mapa <- data.frame(
  muestra = muestra,
  before = c(7, 8, 7, 6, 8, 9, 5, 8),
  after = c(6, 8, 4, 0, 6, 8, 5, 3),
  origen = origen,
  Location = Location,
  LAU_CODE = c("19130", "28079", "28006", "36057", "39024", "03128",
               "06072", #"06498",
               "28009"))

mymunic <-    
  munic |>  
  select(geometry, LAU_CODE)

mapa <- mymunic |> 
  left_join(mapa, by = "LAU_CODE") |> 
  filter(!is.na(Location))

colorsc <- "Set2"
p <- ggplot(esp_can) +
  geom_sf() +
  geom_sf(data = can_prov) +
  geom_sf(data = can_box) + 
  geom_sf(data = mapa, aes(fill = Location), color = "black", size = 1) +
  geom_label_repel(data = mapa, aes(label = muestra, geometry = geometry, colour = Location),
                   stat = "sf_coordinates", size = 4,
                   show.legend = FALSE,
                   # nudge_x = .3,
                   nudge_y = .2) +
  ggtitle("Location for each topsoil") +
  scale_fill_brewer(palette = colorsc) +
  scale_color_brewer(palette = colorsc) +
  # theme_void() +
  theme(legend.position = c(0.18, 0.5),    
        legend.justification = "center",   
        legend.background = element_rect(fill = "white", colour = "black"),
        legend.margin = margin(6, 6, 6, 6),
        legend.text = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.background = element_rect("white", linetype = "blank")) + 
  guides(fill = guide_legend(ncol = 1)) + 
  coord_sf()

ggsave(plot = p, width = 8, height = 6,
       filename = "Suppl_Figure_1.eps",
       device = "eps")

ggsave(plot = p, width = 8, height = 6,
       filename = "Suppl_Figure_1.jpeg",
       device = "jpeg")
