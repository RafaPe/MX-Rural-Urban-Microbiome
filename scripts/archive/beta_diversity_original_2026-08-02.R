library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("patchwork")
library("rbiom")
library(dplyr)
library(readxl)
library(ggrepel)
setwd('C:/Users/rafap/Documents/tesis_maestria/tesis_maestria/scripts')

set.seed(0)

data <- read_excel("../data/metadata_LATINBIOTA_MEXICO.xlsx", sheet = "Data")
metadata <- data[, c("Lane", "Lifestyle", "Gender", "Age group")]
metadata <- as.data.frame(metadata)
rownames(metadata) <- data$Lane




merged_metagenomes <- import_biom("../data/bracken_taxonomy_final.biom")



ids <- merged_metagenomes@sam_data$Id
metadata2 <- as.data.frame(ids)

metadata_final <- metadata[match(ids, row.names(metadata)), , drop = FALSE]

metadata_final <- metadata_final %>%
  rename(Age.group = `Age group`) %>% # Asegurar que el nombre de la columna sea 'Age.group'
  mutate(
    # Crear la nueva columna recodificada
    Age.group.Recoded = case_when(
      # Fusionar adolescentes con adultos
      Age.group %in% c("Teenager", "Adult") ~ "Adult",
      # Fusionar niÃ±os con infantes
      Age.group %in% c("Child", "Infant") ~ "Child",
      # Elderly se queda solo
      Age.group == "Elderly" ~ "Elderly",
      # Mantener cualquier otro valor
      TRUE ~ as.character(Age.group) 
    )
  )

merged_metagenomes@sam_data <- sample_data(metadata_final)

samples_to_exclude <- c("37082_3#16", "36703_3#4", "37035_2#22", "37035_7#18", "36703_3#20", "37082_2#4")


merged_metagenomes_filtered <- subset_samples(
  merged_metagenomes,
  !sample_names(merged_metagenomes) %in% samples_to_exclude
)

# 6. GENERAR EL PLOT DE DIVERSIDAD ALFA CON DATOS FILTRADOS Y RECODIFICADOS
p <- plot_richness(
  physeq = merged_metagenomes_filtered, # Usamos el objeto FILTRADO
  measures = c("Observed","Chao1","Shannon"), 
  x = "Lifestyle",
  color = 'Age.group.Recoded', # Usamos la nueva columna recodificada
  title = "Diversidad Alfa (Muestras Filtradas y Grupos de Edad Recodificados)" # Nuevo tÃ­tulo
) + 
  geom_boxplot() # AÃ±adir boxplots para mejor visualizaciÃ³n

print(p)

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

plot_ordination(physeq = percentages, ordination = meta_ord, color = "Lifestyle", label="Lane", shape="Gender", title="Bracken taxonomy")

p <- plot_ordination(
  physeq = percentages, 
  ordination = meta_ord, 
  color = "Lifestyle", 
  shape = "Gender", 
  title = "Bracken taxonomy"
)

p_final <- p +
  # CAPAS DE MEJORA
  geom_text_repel(aes(label = Lane), size = 3, show.legend = FALSE) +
  geom_point(size = 3.5, alpha = 0.8) +
  
  # TRADUCCIÃ“N DE LEYENDAS Y ETIQUETAS DE VALOR
  scale_color_manual(
    values = c("Rural" = "#72B2A2", "Urban" = "#E59572"), # AsegÃºrate de que los colores sean los que quieres
    labels = c("Rural" = "Rural", "Urban" = "Urbano") 
  ) +
  scale_shape_manual(
    values = c("Female" = 17, "Male" = 15, "NA"=16), # AsegÃºrate de que las formas sean las que quieres
    labels = c("Female" = "Mujer", "Male" = "Hombre") 
  ) +
  
  # TRADUCCIÃ“N DE TÃTULOS Y ETIQUETAS DE EJES
  labs(
    title = "Escalamiento Multidimensional No-MÃ©trico (NMDS)",
    subtitle = paste("Stress del NMDS:", round(meta_ord$stress, 3), 
                     "â€“ MÃ©trica de Distancia: Bray-Curtis"),
    color = "Estilo de Vida", # Nuevo tÃ­tulo de leyenda
    shape = "Sexo",           # Nuevo tÃ­tulo de leyenda
    x = "NMDS1",              # Etiqueta de eje
    y = "NMDS2"               # Etiqueta de eje
  ) +
  theme_bw()

print(p_final)

# 2. AÃ±ade la capa de texto *sin* que aparezca en la leyenda (show.legend = FALSE)
p_final <- p +
  # Usar geom_text_repel requiere la librerÃ­a ggrepel para evitar superposiciones
  geom_text_repel(
    aes(label = Lane), 
    size = 3, # TamaÃ±o de la etiqueta
    show.legend = FALSE # <--- Â¡ESTO ES CRUCIAL! Evita que aparezca en la leyenda.
  ) +
  # Opcional: Ajustar el tamaÃ±o de los puntos
  geom_point(size = 3)

print(p_final)


#37082_3#16
#37082_2#4
#37035_2#22*
#37035_7#18
#36703_3#20
#36703_3#4

percentages2 <- subset_samples(percentages, !(Lane %in% c("37082_3#16", "36703_3#4", "37035_2#22", "37035_7#18", "36703_3#20", "37082_2#4")))
#percentages2 <- subset_samples(percentages, !(Lane %in% c("37082_3#16", "37082_2#4")))
meta_ord <- ordinate(physeq = percentages2, method = "NMDS", distance = "bray")
p2 <- plot_ordination(physeq = percentages2, ordination = meta_ord, color = "Lifestyle", label = "Lane", shape="Gender", title = 'Bracken taxonomy without outliers') 
ggsave('../imgs/bracken.png', dpi=500, units="px", p2) 

p2 <- plot_ordination(
  physeq = percentages2, 
  ordination = meta_ord, 
  color = "Lifestyle", 
  shape = "Gender", 
  title = "Bracken taxonomy"
)

p_final <- p2 +
  # CAPAS DE MEJORA
  geom_text_repel(aes(label = Lane), size = 3, show.legend = FALSE) +
  geom_point(size = 3.5, alpha = 0.8) +
  
  # TRADUCCIÃ“N DE LEYENDAS Y ETIQUETAS DE VALOR
  scale_color_manual(
    values = c("Rural" = "#72B2A2", "Urban" = "#E59572"), # AsegÃºrate de que los colores sean los que quieres
    labels = c("Rural" = "Rural", "Urban" = "Urbano") 
  ) +
  scale_shape_manual(
    values = c("Female" = 17, "Male" = 15, "NA"=16), # AsegÃºrate de que las formas sean las que quieres
    labels = c("Female" = "Mujer", "Male" = "Hombre") 
  ) +
  
  # TRADUCCIÃ“N DE TÃTULOS Y ETIQUETAS DE EJES
  labs(
    title = "Escalamiento Multidimensional No-MÃ©trico (NMDS) sin outliers",
    subtitle = paste("Stress del NMDS:", round(meta_ord$stress, 3), 
                     "â€“ MÃ©trica de Distancia: Bray-Curtis"),
    color = "Estilo de Vida", # Nuevo tÃ­tulo de leyenda
    shape = "Sexo",           # Nuevo tÃ­tulo de leyenda
    x = "NMDS1",              # Etiqueta de eje
    y = "NMDS2"               # Etiqueta de eje
  ) +
  theme_bw()

print(p_final)
                          


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

