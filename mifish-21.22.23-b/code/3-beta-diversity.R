
## Set seed ====
set.seed(123)

## Load libraries ====
library(microViz)
#library(pairwiseAdonis)
library(phyloseq)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(sf)

#taxa_level <- "Species"
#taxa_level <- "Genus"
taxa_level <- "Family"


pcoa.colors <- function(pcoa){
  min0 <- min(pcoa)
  rn0 <- diff(range(pcoa))
  scale.to.color <- function(x, min0, rn0){
    format.hexmode(round((x-min0)/rn0*255))
  }
  sc.out <- apply(pcoa, 2, scale.to.color, 
                  min0 = min0, rn0 = rn0)
  col0 <- apply(sc.out, 1, 
                function(x) paste("#", 
                                  paste(x, collapse = ""),
                                  sep = ""))
  
  return(col0)
}

#Geospatial layers:
moorea <- st_read('./project-data/map-data/moorea.shp')

shore_colors <- c('#1f78b4', '#fb9a99', '#33a02c')
habitat_colors <- brewer.pal(5, "Set1")

#set theme parameters
mytheme <- theme(legend.position="bottom",
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 legend.text=element_text(size=12),
                 legend.title=element_text(size=12))

physeq.rare <- 
  readRDS("./mifish-21.22.23-b/data/physeq-rare.rds")

physeq.rare <- tax_glom(physeq.rare, taxrank=taxa_level) 

data.rare <- as(sample_data(physeq.rare), "data.frame")
# 
# #data.r.ra %>% View()
# 
# bc.pcoa <- phyloseq::ordinate(physeq.r.ra, "PCoA", "bray")
# 
# #saveRDS(bc.pcoa, 'test.pcoa.rds')
# #saveRDS(physeq.r.ra, 'test.physeq.r.ra.rds')
# 
# 
# plot_ordination(physeq.r.ra, bc.pcoa, 
#                          #type = "Site", 
#                          color = "Habitat", 
#                          title = "Bray-Curtis PCoA") +
#   geom_point(size = 4, aes(shape = Island_shore)) +
#   scale_color_manual(values=habitat_colors) +
#   theme_bw()

#PCoA - Bray distance______________________________________________________

ord <- ordinate(physeq.rare, method="PCoA", distance = "bray")

plot_ordination(physeq.rare, ord, 
                #type = "Site", 
                color = "Habitat",
                shape = "Island_shore") +
  geom_point(size = 4, aes(shape = Island_shore)) +
  scale_color_manual(values=habitat_colors) +
  facet_wrap(~Year) +
  theme_bw()

ggsave(paste0('./mifish-21.22.23-b/plots/', 
              taxa_level,
              '_pcoa_bray_21.22.23.png'), 
       width = 10,
       height = 3)

df <- data.frame(data.rare,
                 PCoA1 = ord$vectors[,1],
                 PCoA2 = ord$vectors[,2],
                 PCoA3 = ord$vectors[,3]) %>% 
  dplyr::select(Site, PCoA1, PCoA2, PCoA3,
                Latitude, Longitude, Year) %>% 
  dplyr::filter(!is.na(Longitude)) %>% 
  mutate(erich.hex = 
           dplyr::select(., PCoA1, PCoA2, PCoA3) %>% 
           pcoa.colors()) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs=4326)

bbox <- st_bbox(df)

ggplot() +
  geom_sf(data=moorea, fill='white', color='black') +
  geom_sf(data=df, aes(color=erich.hex), size=3, show.legend = F, alpha=0.75) +
  scale_colour_identity() +
  theme_bw() + mytheme + 
  #guides(color = guide_colourbar(barwidth = 15))+
  facet_wrap(~Year, ncol = 3)

ggsave(paste0('./mifish-21.22.23-b/maps/fish_',
              taxa_level,
              '_PCoA_bray_RGB_21.22.23.png'), 
       width = 14, height = 4, dpi=600)

# presence/absence --------------------------------------------------------

#physeq.pa <- phyloseq_standardize_otu_abundance(physeq.rare, method = "pa")
physeq.pa <- tax_transform(physeq.rare, trans = "binary")
ord = ordinate(physeq.pa, method="PCoA", distance = "jaccard")

plot_ordination(physeq.pa, ord, 
                #type = "Site", 
                color = "Habitat",
                shape = "Island_shore") +
  geom_point(size = 4, aes(shape = Island_shore)) +
  scale_color_manual(values=habitat_colors) +
  facet_wrap(~Year) +
  theme_bw()

ggsave(paste0('./mifish-21.22.23-b/plots/pres_abs_',
              taxa_level,
              '_pcoa_jaccard_21.22.23.png'), 
       width = 10,
       height = 3)

df <- data.frame(data.rare,
                 PCoA1 = ord$vectors[,1],
                 PCoA2 = ord$vectors[,2],
                 PCoA3 = ord$vectors[,3]) %>% 
  dplyr::select(Site, PCoA1, PCoA2, PCoA3,
                Latitude, Longitude, Year) %>% 
  dplyr::filter(!is.na(Longitude)) %>% 
  mutate(erich.hex = 
           dplyr::select(., PCoA1, PCoA2, PCoA3) %>% 
           pcoa.colors()) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs=4326)

ggplot() +
  geom_sf(data=moorea, fill='white', color='black') +
  geom_sf(data=df, aes(color=erich.hex), size=3, show.legend = F, alpha=0.75) +
  scale_colour_identity() +
  theme_bw() + mytheme + 
  #guides(color = guide_colourbar(barwidth = 15))+
  facet_wrap(~Year, ncol = 3)

ggsave(paste0('./mifish-21.22.23-b/maps/pres_abs_',
              taxa_level,
              '_pcoa_jaccard_RGB_21.22.23.png'), 
       width = 9, height = 4, dpi=600)


