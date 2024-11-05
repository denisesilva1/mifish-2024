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
physeq.rare <- 
  readRDS("./mifish-21.22/data/physeq-rare.rds")

sample_data(physeq.rare)$sample_id <- 
  sample_data(physeq.rare) %>% 
  row.names()

#sample_data(physeq.rare) %>% View()

#physeq.rare = subset_samples(physeq.rare, Site != "LTER 4B")
#physeq.rare = subset_samples(physeq.rare, Site != "56")
#saveRDS(physeq.rare, "./physeq-rare.rds")

#Fix missing data in 2022
data.rare <- as(sample_data(physeq.rare), "data.frame")

#View(data.rare)

## Validate phyloseq object ====

physeq.rare <- phyloseq_validate(physeq.rare) # NAs detected

## Make taxa names uniquely identifiable ====
physeq.rare <- 
  tax_fix(
    physeq.rare,
    min_length = 4,
    unknowns = NA,
    suffix_rank = "classified",
    sep = " ",
    anon_unique = TRUE,
    verbose = TRUE
  ) 

#Filter the Ncs
#physeq.rare <- 
  #subset_samples(physeq.rare, genus != "Negative Control")

OTU <- as(otu_table(physeq.rare), "matrix")
#if( taxa_are_rows(physeq_rare) ){ OTU <- t(OTU) }
rich <- data.frame(t(estimateR(t(OTU))))['S.obs']
names(rich) <- "Observed"

FaithPD = picante::pd(samp = t(OTU), 
                      tree = phy_tree(physeq.rare), 
                      include.root = T)["PD"]
names(FaithPD) <- "FaithPD"

shannon = diversity(t(OTU), index="shannon") %>% 
  data.frame()
names(shannon) <- "Shannon"

chao1 = estimateR(t(OTU)) %>% 
  t() %>% 
  data.frame() %>% 
  select(S.chao1)

names(chao1) <- "Chao1"


erich <- cbind(rich, shannon, FaithPD, chao1) 
  
head(erich)

##Add evenness as a function of shannon
erich$evenness <- erich$Shannon/log(erich$Observed) 
erich <- erich %>% 
  mutate(evenness = ifelse(is.na(evenness), 1, evenness))


#Fill out the rest of the data from the phyloseq object; missing Year = 2022
erich <- cbind(erich, data.rare) 
rownames(erich) <- rownames(data.rare)

#Remove and include new calculations
lookup <- c('Fish_Species_Richness'='Observed',
            'Fish_Shannon_Diversity'='Shannon',
            'Fish_Phylogenetic_Diversity'='FaithPD',
            'Fish_Chao1_Diversity' = 'Chao1',
            'Fish_Evenness'='evenness')

erich <- dplyr::rename(erich, all_of(lookup)) 

#View(erich)

saveRDS(erich, file = "./mifish-21.22/data/12s-erich_rare1000.RDS", compress = TRUE)


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

ggsave(file = './mifish-21.22/plots/Fish_Shannon_Diversity_habitat.png', 
       dpi = 600, width = 25, height = 8, units = "cm")


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

ggsave(file = './mifish-21.22/plots/Fish_Phylogenetic_Diversity.png', 
       dpi = 600, width = 25, height = 8, units = "cm")


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

ggsave(file = './mifish-21.22/plots/Fish_Species_Richness.png', 
       dpi = 600, width = 25, height = 8, units = "cm")


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

ggsave(file = './mifish-21.22/plots/Fish_Chao1_Diversity.png', 
       dpi = 600, width = 25, height = 8, units = "cm")



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

ggsave(file = './mifish-21.22/plots/Fish_Evenness.png', 
       dpi = 600, width = 25, height = 8, units = "cm")
