library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("patchwork")
library("rbiom")
library(readxl)

set.seed(0)

data <- read_excel("../data/metadata_LATINBIOTA_MEXICO.xlsx", sheet = "Data")
metadata <- data[, c("Lane", "Lifestyle", "Gender", "Age group")]
metadata <- as.data.frame(metadata)
rownames(metadata) <- data$Lane




merged_metagenomes <- import_biom("../data/bracken_final.biom")



ids <- merged_metagenomes@sam_data$Id
metadata2 <- as.data.frame(ids)

metadata_final <- metadata[match(ids, row.names(metadata)), , drop = FALSE]

merged_metagenomes@sam_data <- sample_data(metadata_final)

p <- plot_richness(physeq = merged_metagenomes, measures = c("Observed","Chao1","Shannon"), x="Lifestyle", color = 'Age.group', title = "Alpha diversity")
ggsave('../imgs/observed_diversity.png', dpi=500, units="px", p)

p <- plot_richness(physeq = merged_metagenomes, measures = c("Observed"), color = 'Lifestyle', sortby = "Observed" ,title = "Alpha diversity")+theme(axis.text.x = element_text(size = 3.5)) 
ggsave('../imgs/observed_diversity2.png', dpi=500, width = 5000, units="px", p)


merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4)
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
percentages <- transform_sample_counts(merged_metagenomes, function(x) x*100 / sum(x) )



meta_ord <- ordinate(physeq = percentages, method = "NMDS", distance = "bray")
p <- plot_ordination(physeq = percentages, ordination = meta_ord, color = "Lifestyle", label="Lane", shape="Gender", title="Bracken taxonomy")
ggsave('../imgs/bracken_woutliers.png', dpi=500, units="px", p)


#37082_3#16
#37082_2#4
#37035_2#22
#37035_7#18
#36703_3#20
#36703_3#4

percentages2 <- subset_samples(percentages, !(Lane %in% c("37082_3#16", "36703_3#4", "37035_2#22", "37035_7#18", "36703_3#20", "37082_2#4")))
meta_ord <- ordinate(physeq = percentages2, method = "NMDS", distance = "bray")
p2 <- plot_ordination(physeq = percentages2, ordination = meta_ord, color = "Lifestyle", label = "Lane", shape="Gender", title = 'Bracken taxonomy without outliers') 
ggsave('../imgs/bracken.png', dpi=500, units="px", p2)                            
                          


percentages3 <- subset_taxa(percentages2, Kingdom == "Bacteria")
meta_ord <- ordinate(physeq = percentages3, method = "NMDS", distance = "bray")
p3 <- plot_ordination(physeq = percentages3, ordination = meta_ord, color = "Lifestyle", label = "Lane", shape="Gender", title = 'Taxonomy (Bacteria)') 
ggsave('../imgs/bracken_onlyBacteria.png', dpi=500, units="px", p3)  


percentages3 <- subset_taxa(percentages2, Kingdom == "Archaea")
meta_ord <- ordinate(physeq = percentages3, method = "NMDS", distance = "bray")
p3 <- plot_ordination(physeq = percentages3, ordination = meta_ord, color = "Lifestyle", label = "Lane", shape="Gender", title = 'Taxonomy (Archaea)') 
ggsave('../imgs/bracken_onlyArchaea.png', dpi=500, units="px", p3) 


percentages3 <- subset_taxa(percentages2, Kingdom == "Eukaryota")
meta_ord <- ordinate(physeq = percentages3, method = "NMDS", distance = "bray")
p3 <- plot_ordination(physeq = percentages3, ordination = meta_ord, color = "Lifestyle", label = "Lane", shape="Gender", title = 'Taxonomy (Eukaryota)') 
ggsave('../imgs/bracken_onlyEukaryota.png', dpi=500, units="px", p3)  


percentages3 <- subset_taxa(percentages2, Kingdom == "Viruses")
meta_ord <- ordinate(physeq = percentages3, method = "NMDS", distance = "bray")
p3 <- plot_ordination(physeq = percentages3, ordination = meta_ord, color = "Lifestyle", label = "Lane", shape="Gender", title = 'Taxonomy (Virus)') 
ggsave('../imgs/bracken_onlyVirus.png', dpi=500, units="px", p3)  


TopNGenus <- names(sort(taxa_sums(merged_metagenomes), TRUE)[1:5])
Top5Genus <- prune_taxa(TopNGenus,merged_metagenomes)
plot_heatmap(Top5Genus)
