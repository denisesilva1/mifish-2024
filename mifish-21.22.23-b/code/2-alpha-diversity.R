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
physeq_rrfy <- 
  readRDS("./mifish-21.22.23-b/data/physeq-rare.rds")
data.rrfy <- as(sample_data(physeq_rrfy), "data.frame") 

OTU <- 
  as(phyloseq::otu_table(physeq_rrfy), 
     "matrix") %>% 
  t()

# Calculate alpha diversity metrics ---------------------------------------
rich <- 
  vegan::estimateR(OTU) %>% 
  t() %>% 
  data.frame() %>% 
  dplyr::select(S.obs, S.chao1) 

faith.pd <- 
  picante::pd(samp = OTU,
              tree = phyloseq::phy_tree(physeq_rrfy), 
              include.root = FALSE) %>% 
  dplyr::select(PD)

shannon <- 
  vegan::diversity(OTU, 
                   index = "shannon") %>% 
  as_tibble() %>% 
  dplyr::rename(shannon = value)

# Combine alpha diversity metrics -----------------------------------------
lookup <- c('Fish_Species_Richness'='S.obs',
            'Fish_Shannon_Diversity'='shannon',
            'Fish_Phylogenetic_Diversity'='PD',
            'Fish_Chao1_Diversity' = 'S.chao1',
            'Fish_Evenness'='evenness')

alpha_div <- 
  cbind(rich, shannon, faith.pd) %>% 
  dplyr::mutate(evenness = shannon / log(S.obs)) %>% 
  cbind(data.rrfy) %>% 
  dplyr::rename(all_of(lookup)) %>% 
  dplyr::relocate(Site, Year) %>% 
  dplyr::filter(!is.na(Latitude)) %>% 
  as_tibble()

write_rds(alpha_div, 
          file = './mifish-21.22.23-b/data/ati-alpha-diversity-metrics-species.rds', 
          compress = "xz")

ggplot() +
  geom_boxplot(data = alpha_div,
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

ggsave(file = './mifish-21.22.23-b/plots/Fish_Shannon_Diversity_habitat.png', 
       dpi = 600, width = 25, height = 8, units = "cm")


ggplot() +
  geom_boxplot(data = alpha_div,
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

ggsave(file = './mifish-21.22.23-b/plots/Fish_Phylogenetic_Diversity.png', 
       dpi = 600, width = 25, height = 8, units = "cm")


ggplot() +
  geom_boxplot(data = alpha_div,
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

ggsave(file = './mifish-21.22.23-b/plots/Fish_Species_Richness.png', 
       dpi = 600, width = 25, height = 8, units = "cm")


ggplot() +
  geom_boxplot(data = alpha_div,
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

ggsave(file = './mifish-21.22.23-b/plots/Fish_Chao1_Diversity.png', 
       dpi = 600, width = 25, height = 8, units = "cm")



ggplot() +
  geom_boxplot(data = alpha_div,
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

ggsave(file = './mifish-21.22.23-b/plots/Fish_Evenness.png', 
       dpi = 600, width = 25, height = 8, units = "cm")
