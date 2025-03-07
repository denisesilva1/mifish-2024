#load libraries
library(qiime2R)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(gridExtra)
library(vegan)
library(decontam)
library(data.table)


# Make Phyloseq Object ----------------------------------------------------

# metadata (Just to see how the table looks)
fread("./mifish-21.22/data/metadata_mifish_2024.txt", sep="\t") %>%
  View()

#Build a phyloseq object
physeq <- qza_to_phyloseq(features = './mifish-21.22/data/mifish_12s-nobact_table.qza',
                          tree = './mifish-21.22/data/mifish-rooted-tree.qza', 
                          taxonomy = './mifish-21.22/data/mifish_12s-tax.qza',
                          metadata = './mifish-21.22/data/metadata_mifish_2024.txt')


#Remove unassigned reads (unassigned at the kingdom level) & for good measure the euks
#physeq <- subset_taxa(physeq,  Kingdom != "Unassigned") 
physeq <- subset_taxa(physeq, Order != "Chloroplast")
#physeq <- subset_taxa(physeq, Kingdom != "Bacteria")

#make a sample data frame
sample.data <- as(sample_data(physeq), "data.frame") 
#View(sample.data)

# Decontamination ----------------------------------------------------

#Start by looking at the library size
sample.data$LibrarySize <- sample_sums(physeq)
sample.data <- sample.data[order(sample.data$LibrarySize),]
sample.data$Index <- seq(nrow(sample.data))

#visualize library size 
ggplot() +
  geom_point(data = sample.data, 
             aes(x=Index, y=LibrarySize)) +
  theme_bw()

#Remove singletons
physeq <- prune_taxa(taxa_sums(physeq) > 1, physeq)

#Check final numbers
print(microbiome::summarize_phyloseq(physeq)) 
#Total reads = /mifish-21.22
# Max. number of reads = 93311
#Avg reads/sample = 10298.7406417112
#Min No. of reads = 0

# Make two finalized phyloseq objects, one non-rarefied and one rarefied
# physeq.final is non rarefied
write_rds(physeq, file = "./mifish-21.22/data/physeq_norare.rds")

# Make a rarefied phyloseq object for alpha diversity analyses
otu <- otu_table(physeq)
class(otu) <- "matrix"
otu <- t(otu)
rarecurve(otu, step=50, label = F)
abline(v = 1000, col = 'red')


physeq_rare <- 
  rarefy_even_depth(physeq, 
                    sample.size = 1000, 
                    rngseed = 1020) #Set seed to be reproducible

write_rds(physeq_rare, file = "./mifish-21.22/data/physeq-rare.rds")

##Export the raw otu and taxonomy file for use in downstream analyses
otu <- otu_table(physeq)
tax <- tax_table(physeq)
otu.rare <- otu_table(physeq_rare)
tax.rare <- tax_table(physeq_rare)
write.csv(otu, "./mifish-21.22/data/12s-asv-table.csv", row.names = F)
write.csv(tax, "./mifish-21.22/data/12s-tax-table.csv", row.names = F)
write.csv(otu, "./mifish-21.22/data/12s-rarefied-asv-table.csv", row.names = F)
write.csv(tax, "./mifish-21.22/data/12s-rarefied-tax-table.csv", row.names = F)

