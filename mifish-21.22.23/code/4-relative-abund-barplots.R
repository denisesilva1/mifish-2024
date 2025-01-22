library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(tidyverse)
library(readr)

taxa_level <- "Species"
#taxa_level <- "Genus"
#taxa_level <- "Family"

#topn <- 10
topn <- 20

pal <- brewer.pal(12,"Paired") 
# Add more colors to this palette :
pal <- colorRampPalette(pal)(topn)

## Load data ====
physeq.norare23 <- 
  read_rds("./mifish-23/data/physeq_norare.rds")
physeq.norare23@phy_tree <- NULL

physeq.norare <- 
  read_rds("./mifish-21.22/data/physeq_norare.rds")
physeq.norare@phy_tree <- NULL

physeq.norare <- 
  phyloseq::merge_phyloseq(physeq.norare, physeq.norare23)

#write_rds(physeq.norare, "./mifish-21.22.23/data/physeq-norare-21.22.23.rds")

sample_data(physeq.norare)$Habitat <- 
  as.factor(sample_data(physeq.norare)$Habitat)

sample_data(physeq.norare)$hab_year <- 
  paste0(sample_data(physeq.norare)$Habitat, 
         sample_data(physeq.norare)$Year)

sample_data(physeq.norare)$common <- "common"

# top2 <- 
#   physeq.norare %>% 
#   merge_samples("common") %>% 
#   tax_glom(taxrank = taxa_level) %>% 
#   transform_sample_counts(function(x) {(x/sum(x))*100}) %>% 
#   psmelt() %>% 
#   as_tibble() %>% 
#   dplyr::select(taxa_level, 
#                 Abundance) %>% 
#   rename(Taxa = taxa_level) %>% 
#   mutate(Taxa = str_sub(Taxa, start = 3))
  
relabund <-
  physeq.norare %>% 
  merge_samples("hab_year") %>% 
  #tax_glom(taxrank = "Genus") %>% 
  tax_glom(taxrank = taxa_level) %>% 
  transform_sample_counts(function(x) {(x/sum(x))*100}) %>% 
  psmelt() %>% 
  as_tibble() %>% 
  #dplyr::select(Genus, Year, Habitat, Abundance) %>% 
  dplyr::select(taxa_level, Year, Habitat, 
                Abundance, Latitude, Longitude) %>% 
  mutate(Habitat = as.character(Habitat)) %>% 
  mutate(Habitat = case_when(
    Habitat == '1' ~ 'Bay',
    Habitat == '2' ~ 'Fringing Reef',
    Habitat == '3' ~ 'Mid-Lagoon',
    Habitat == '4' ~ 'Reef Crest',
    Habitat == '5' ~ 'Reef Pass')
  ) %>% 
  rename(Taxa = taxa_level) %>% 
  mutate(Taxa = str_sub(Taxa, start = 3))

relabund %>% 
  dplyr::group_by(Taxa) %>% 
  #dplyr::group_by(Genus) %>% 
  dplyr::summarise(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  arrange(desc(Abundance)) %>% 
  View()

top <- relabund %>% 
  dplyr::group_by(Taxa) %>% 
  #dplyr::group_by(Genus) %>% 
  dplyr::summarise(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  arrange(desc(Abundance)) #%>% View()
  slice_head(n = topn) %>% 
  pull(Taxa)
#pull(Genus)

# Write out names of top 10
write_rds(top, paste0("./mifish-21.22.23/data/top",topn,"_",taxa_level,".rds"))

tax_top <- relabund %>% 
  filter(Taxa %in% top)
#filter(Genus %in% top)

plot_RelAbund <- 
  ggplot(data = tax_top, 
         aes(x = Habitat, y = Abundance, 
             #fill = Family)) +
             fill = Taxa)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = pal, name = taxa_level) +
  # scale_fill_manual(values = c("#644194", "#CAB2D6",
  #                              "#4A77AF", "#B3CDE1",
  #                              "#5D9D40", "#BDDD94",
  #                              "#E31A1C", "#FB9A99",
  #                              "#FF7F00", "#FDBF6F")) + 
  # scale_x_discrete(labels = c("Reef Pass", "Reef Crest", "Mid-Lagoon",
  #                             "Fringing Reef", "Bay")) +
  ggh4x::facet_wrap2(~ Year, ncol = 1,
                     axes = "y", remove_labels = "y",
                     scale = "free_y") +
  coord_flip() +
  #facet_wrap(~ Year, ncol = 3, scale="free") +
  #labs(title = "Relative Abundance") + 
  ylab("Relative Abundance (%)") + 
  xlab("Habitat") + 
  theme_bw() + 
  theme(axis.text = element_text(face = "bold", size = 11.5), 
        axis.title = element_text(face = "bold", size = 12), 
        title = element_text(face = "bold")) 

plot_RelAbund

ggplot2::ggsave(#"./mifish-21.22.23/plots/fish_ati_relabund_family_top20_21.22.23.png", 
  paste0("./mifish-21.22.23/plots/fish_ati_relabund_",
         taxa_level,
         "_top",
         as.character(topn),
         "_21.22.23.png"),
  plot_RelAbund, 
  height = 300, 
  width = 700,
  units = "mm",
  scale = 0.5, 
  dpi = 1000)
