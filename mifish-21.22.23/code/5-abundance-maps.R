#load libraries
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(sf)
library(RColorBrewer)
library(ggspatial)
library(colorspace)

#Geospatial layers:
moorea <- st_read('./project-data/map-data/moorea.shp')

topn <- 10
taxa_level <- 'Genus'

pal <- brewer.pal(12,"Paired") 
# Add more colors to this palette :
pal <- colorRampPalette(pal)(topn)

dark <- darken(pal, 0.85)
light <- lighten(pal, 0.99)

#write_rds(top, "./mifish-21.22.23/data/top20_genus.rds")
top.names <- 
  read_rds(paste0("./mifish-21.22.23/data/top",topn,"_",taxa_level,".rds"))

top.names <- top.names[order(top.names)]

## Load data ====
physeq.norare <- 
  readRDS("./mifish-21.22.23/data/physeq-norare-21.22.23.rds") 

# sample_data(physeq.rare)$Habitat <- 
#   as.factor(sample_data(physeq.rare)$Habitat)

# sample_data(physeq.rare)$hab_year <- 
#   paste0(sample_data(physeq.rare)$Habitat, sample_data(physeq.rare)$Year)

relabund <-
  physeq.norare %>% 
  tax_glom(taxrank = "Genus") %>% 
  transform_sample_counts(function(x) {x/sum(x)}) %>% 
  psmelt() %>% 
  dplyr::select(Site, Genus, Year, Habitat, 
                Abundance, Latitude, Longitude) %>% 
  mutate(
    #Family = stringr::str_remove(Family, "F_")
    Genus = stringr::str_remove(Genus, "G_")
  ) %>% 
  group_by(Site, Genus, Year) %>% 
  slice(which.max(Abundance))

tax_top <- relabund %>% 
  #filter(Family %in% top)
  filter(Genus %in% top.names) %>% 
  group_by(Site, Genus, Year) %>% 
  slice(which.max(Abundance)) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"))

st_crs(tax_top) <- st_crs(moorea)

bbox <- st_bbox(tax_top)

#set theme parameters
mytheme <- theme(legend.position="bottom",
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 legend.text=element_text(size=12),
                 legend.title=element_text(size=12)) 

for(i in 1:length(top.names)){
  midpoint <- tax_top %>% 
    st_drop_geometry() %>% 
    filter(Genus == top.names[i]) %>% 
    group_by(Genus) %>% 
    summarise(mn = min(Abundance, na.rm = T),
              mx = max(Abundance, na.rm = T)) %>% 
    mutate(midpoint = (mx - mn)/3) %>% 
    pull(midpoint)
  
  tmp <- 
    tax_top %>% 
    filter(Genus == top.names[i]) %>% 
    arrange(-Abundance)
  
  m <- ggplot() +
    geom_sf(data=moorea, fill='white', color='black') +
    geom_sf(data=tmp, 
            aes(color=Abundance), 
            size=2) +
    theme_bw() +
    mytheme + 
    scale_color_gradient2(name = top.names[i],
                         low = light[i],
                         mid = pal[i],
                         high = dark[i],
                         midpoint = midpoint) +
    #scale_color_distiller(palette='Blues', top[i], direction = 1) +
    guides(color = guide_colourbar(barwidth = 15)) +
    facet_wrap(~Year, ncol = 3)
  
  ggsave(paste0("./mifish-21.22.23/maps/relabund_",top.names[i],"_21.22.23.png"), 
         width = 10.75, height = 4, dpi=600)
}


chlorurus <- 
  relabund %>% 
  filter(Genus == "Chlorurus") %>% 
  select(Site, Genus, Abundance, Year) %>%
  mutate(site_year = paste0(Site, "_", Year)) %>% 
  left_join(
    read_rds("./mifish-21.22.23/data/alpha-diversity-metrics-12s-rare1000.rds") %>% 
      dplyr::select(-Site, -Year),
    join_by(site_year)
  )

write_rds(chlorurus, 
          './mifish-21.22.23/data/chlorurus-rel-abund-biodiversity-metrics-21.22.23.rds')
