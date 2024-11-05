
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

#set theme parameters
mytheme <- theme(legend.position="bottom",
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 legend.text=element_text(size=12),
                 legend.title=element_text(size=12))
habitat_colors <- brewer.pal(5, "Set1")

## Load data ====
physeq.rare <- 
  readRDS("./mifish-21.22/data/physeq-rare.rds")

sample_data(physeq.rare)$Year <- as.factor(sample_data(physeq.rare)$Year)

## Validate phyloseq object and fix if necessary
physeq.rare <- phyloseq_validate(physeq.rare)

physeq.rare <- tax_fix(
  physeq.rare,
  min_length = 4,
  unknowns = NA,
  suffix_rank = "classified",
  sep = " ",
  anon_unique = TRUE,
  verbose = TRUE) 

physeq.rare <- tax_glom(physeq.rare, taxrank="Genus") 

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

ggsave('./mifish-21.22/plots/genus_pcoa_bray_21.22.png', width = 10)

df <- data.frame(data.rare,
                 PCoA1 = ord$vectors[,1],
                 PCoA2 = ord$vectors[,2],
                 PCoA3 = ord$vectors[,3]) %>% 
  dplyr::select(Site, PCoA1, PCoA2, PCoA3,
                Latitude, Longitude, Year) %>% 
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

ggsave(filename='./mifish-21.22/maps/fish_genus_PCoA_bray_RGB_21.22.png', 
       width = 9, height = 4, dpi=600)

#Weight_unifrac-----------------------------------------------------------------
# take into account differences in abundance of taxa between samples, but takes longer to calculate
physeq.rare.trans <- microbiome::transform(physeq.rare, "compositional")

ord <- ordinate(physeq.rare.trans , "PCoA", "unifrac", weighted=T)

plot_ordination(physeq.rare.trans, 
                ord, 
                #type = "Site", 
                color = "Habitat",
                shape = "Island_shore") +
  geom_point(size = 4, aes(shape = Island_shore)) +
  scale_color_manual(values=habitat_colors) +
  facet_wrap(~Year) +
  theme_bw()

ggsave('./mifish-21.22/plots/genus_pcoa_unifrac_21.22.png', width = 10)

df <- data.frame(data.rare,
                 PCoA1 = ord$vectors[,1],
                 PCoA2 = ord$vectors[,2],
                 PCoA3 = ord$vectors[,3]) %>% 
  dplyr::select(Site, PCoA1, PCoA2, PCoA3,
                Latitude, Longitude, Year) %>% 
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

ggsave(filename='./mifish-21.22/maps/fish_genus_PCoA_unifrac_RGB_21.22.png', 
       width = 9, height = 4, dpi=600)


#Unweight Frac Consumer_pressure_________________________________________________________________________-
# considers the presence/absence of taxa between sample pairs.

#physeq.rare.trans <- microbiome::transform(physeq.rare, "compositional")
physeq.pa <- tax_transform(physeq.rare, trans = "binary")

ord <- ordinate(physeq.pa , "PCoA", "unifrac", weighted=F)

plot_ordination(physeq.pa, 
                ord, 
                color = "Habitat",
                shape = "Island_shore") +
  geom_point(size = 4, aes(shape = Island_shore)) +
  scale_color_manual(values=habitat_colors) +
  theme_bw() +
  facet_wrap(~Year)

ggsave('./mifish-21.22/plots/pres_abs_genus_pcoa_unifrac_21.22.png', width = 10)

df <- data.frame(data.rare,
                 PCoA1 = ord$vectors[,1],
                 PCoA2 = ord$vectors[,2],
                 PCoA3 = ord$vectors[,3]) %>% 
  dplyr::select(Site, PCoA1, PCoA2, PCoA3,
                Latitude, Longitude, Year) %>% 
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

ggsave(filename='./mifish-21.22/maps/pres_abs_genus_pcoa_unifrac__RGB_21.22.png', 
       width = 9, height = 4, dpi=600)

# presence/absence --------------------------------------------------------

#physeq.pa <- phyloseq_standardize_otu_abundance(physeq.rare, method = "pa")

ord = ordinate(physeq.pa, method="PCoA", distance = "jaccard")

plot_ordination(physeq.pa, ord, 
                #type = "Site", 
                color = "Habitat",
                shape = "Island_shore") +
  geom_point(size = 4, aes(shape = Island_shore)) +
  scale_color_manual(values=habitat_colors) +
  facet_wrap(~Year) +
  theme_bw()

ggsave('./mifish-21.22/plots/pres_abs_genus_pcoa_jaccard_21.22.png', width = 10)

df <- data.frame(data.rare,
                 PCoA1 = ord$vectors[,1],
                 PCoA2 = ord$vectors[,2],
                 PCoA3 = ord$vectors[,3]) %>% 
  dplyr::select(Site, PCoA1, PCoA2, PCoA3,
                Latitude, Longitude, Year) %>% 
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

ggsave(filename='./mifish-21.22/maps/pres_abs_genus_pcoa_jaccard_RGB_21.22.png', 
       width = 9, height = 4, dpi=600)


# 
# df.colors <- 
#   df  %>% 
#   dplyr::select(c(PCoA1, PCoA2, PCoA3))
# 
# erich.hex <- pcoa.colors(df.colors)
# 
# df <- cbind(df, erich.hex) %>% 
#   st_as_sf(coords = c("Longitude", "Latitude"), crs=4326)
# 
# 
# 
#  
# 
# m15 <- ggplot() +
#   geom_sf(data=moorea, fill='white', color='black') +
#   geom_sf(data=df, aes(color=erich.hex), size=3, show.legend = F, alpha=0.75) +
#   scale_colour_identity() +
#   theme_bw() + mytheme + 
#   #guides(color = guide_colourbar(barwidth = 15))+
#   facet_wrap(~Year, ncol = 3)
# m15
# ggsave(m15, filename='./maps/fish_family_PCoA_RGB_21.22.23.png', 
#        width = 9, height = 4, dpi=600)
# 
# 
