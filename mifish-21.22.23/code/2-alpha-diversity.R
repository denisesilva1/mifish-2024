#Load all the libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyverse)
library(cowplot)
library(vegan)
library(microViz)
library(ggpubr)
library(rstatix) 

## Load data ====
erich <- read_rds(file = "./mifish-23/data/12s-erich_rare1000.RDS") %>% 
  mutate(site_year = paste0(Site, "_", Year)) %>% 
  select(Site, 
         site_year,
         Year,
         Latitude, 
         Longitude, 
         Habitat,
         Island_shore,
         Fish_Species_Richness,
         Fish_Shannon_Diversity,
         Fish_Phylogenetic_Diversity,
         Fish_Chao1_Diversity,
         Fish_Evenness)

erich.tmp <- read_rds(file = "./mifish-21.22/data/12s-erich_rare1000.RDS") %>% 
  mutate(site_year = paste0(Site, "_", Year)) %>% 
  select(Site, 
         site_year,
         Year,
         Latitude, 
         Longitude, 
         Habitat,
         Island_shore,
         Fish_Species_Richness,
         Fish_Shannon_Diversity,
         Fish_Phylogenetic_Diversity,
         Fish_Chao1_Diversity,
         Fish_Evenness) 

erich <- rbind(erich, erich.tmp)

write_rds(erich, "./mifish-21.22.23/data/alpha-diversity-metrics-12s-rare1000.rds")

ggplot() +
  geom_boxplot(data = erich,
               aes(x = Habitat,
                   y = Fish_Shannon_Diversity, 
                   fill = Habitat)) +
  #scale_y_continuous(trans='log10') +
  theme_bw() +
  theme(axis.line = element_line(colour = "gray"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  xlab("Habitat") +
  ylab("Shannon Diversity - Fish") +
  guides(fill=guide_legend()) +
  facet_wrap(~Year)
#legend.position = "none")

ggsave(file = './mifish-21.22.23/plots/Fish_Shannon_Diversity_habitat.png', 
       dpi = 600, width = 35, height = 8, units = "cm")


ggplot() +
  geom_boxplot(data = erich,
               aes(x = Habitat,
                   y = Fish_Phylogenetic_Diversity, 
                   fill = Habitat)) +
  #scale_y_continuous(trans='log10') +
  theme_bw() +
  theme(axis.line = element_line(colour = "gray"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  xlab("Habitat") +
  ylab("Phylogenetic Diversity - Fish") +
  guides(fill=guide_legend()) +
  facet_wrap(~Year)
#legend.position = "none")

ggsave(file = './mifish-21.22.23/plots/Fish_Phylogenetic_Diversity.png', 
       dpi = 600, width = 35, height = 8, units = "cm")


ggplot() +
  geom_boxplot(data = erich,
               aes(x = Habitat,
                   y = Fish_Species_Richness, 
                   fill = Habitat)) +
  #scale_y_continuous(trans='log10') +
  theme_bw() +
  theme(axis.line = element_line(colour = "gray"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  xlab("Habitat") +
  ylab("Observed Richness - Fish") +
  guides(fill=guide_legend()) +
  facet_wrap(~Year)
#legend.position = "none")

ggsave(file = './mifish-21.22.23/plots/Fish_Species_Richness.png', 
       dpi = 600, width = 35, height = 8, units = "cm")


ggplot() +
  geom_boxplot(data = erich,
               aes(x = Habitat,
                   y = Fish_Chao1_Diversity, 
                   fill = Habitat)) +
  #scale_y_continuous(trans='log10') +
  theme_bw() +
  theme(axis.line = element_line(colour = "gray"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  xlab("Habitat") +
  ylab("Chao Index - Fish") +
  guides(fill=guide_legend()) +
  facet_wrap(~Year)
#legend.position = "none")

ggsave(file = './mifish-21.22.23/plots/Fish_Chao1_Diversity.png', 
       dpi = 600, width = 35, height = 8, units = "cm")



ggplot() +
  geom_boxplot(data = erich,
               aes(x = Habitat,
                   y = Fish_Evenness, 
                   fill = Habitat)) +
  #scale_y_continuous(trans='log10') +
  theme_bw() +
  theme(axis.line = element_line(colour = "gray"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  xlab("Habitat") +
  ylab("Evenness - Fish") +
  guides(fill=guide_legend()) +
  facet_wrap(~Year)
#legend.position = "none")

ggsave(file = './mifish-21.22.23/plots/Fish_Evenness.png', 
       dpi = 600, width = 35, height = 8, units = "cm")

