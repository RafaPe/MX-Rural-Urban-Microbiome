library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("patchwork")
library("rbiom")
library(dplyr)
library(readxl)
library(ggrepel)
# Buscar la raiz del proyecto desde el directorio de trabajo actual hacia arriba.
# Esto evita rutas NA y funciona desde la raiz, scripts/ u otra subcarpeta.
find_project_root <- function(start = getwd()) {
  current_dir <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    metadata_candidate <- file.path(
      current_dir,
      "data",
      "metadata_LATINBIOTA_MEXICO.xlsx"
    )

    if (file.exists(metadata_candidate)) {
      return(current_dir)
    }

    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      break
    }
    current_dir <- parent_dir
  }

  stop(
    "No se encontro la raiz del proyecto desde: ", getwd(),
    ". Abre el proyecto tesis_maestria - copia y ejecuta el script de nuevo."
  )
}

project_root <- find_project_root()
data_dir <- file.path(project_root, "data")
img_dir <- file.path(project_root, "img")
dir.create(img_dir, recursive = TRUE, showWarnings = FALSE)

metadata_path <- file.path(data_dir, "metadata_LATINBIOTA_MEXICO.xlsx")
biom_path <- file.path(data_dir, "bracken_taxonomy_final.biom")

if (!file.exists(biom_path)) {
  stop("No se encontro el archivo BIOM en: ", biom_path)
}

set.seed(0)

data <- read_excel(
  metadata_path,
  sheet = "Data"
)
metadata <- data[, c("Lane", "Lifestyle", "Gender", "Age group")]
metadata <- as.data.frame(metadata)
rownames(metadata) <- data$Lane




merged_metagenomes <- import_biom(
  biom_path
)



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

samples_to_exclude <- c("37082_3#16", "36703_3#4", "37035_7#18", "36703_3#20", "37082_2#4", "37035_2#22")


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
ggsave(file.path(img_dir, 'observed_diversity.png'), dpi=500, units="px", p)

p <- plot_richness(physeq = merged_metagenomes, measures = c("Observed"), color = 'Lifestyle', sortby = "Observed" ,title = "Alpha diversity")+theme(axis.text.x = element_text(size = 3.5)) 
ggsave(file.path(img_dir, 'observed_diversity2.png'), dpi=500, width = 5000, units="px", p)


merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4)
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
percentages <- transform_sample_counts(merged_metagenomes, function(x) x*100 / sum(x) )



# NMDS con todas las muestras. Los outliers se conservan y se destacan
# con un color y un simbolo diferentes a los grupos de estilo de vida.
meta_ord_with_outliers <- ordinate(
  physeq = percentages,
  method = "NMDS",
  distance = "bray"
)

p_with_outliers_base <- plot_ordination(
  physeq = percentages,
  ordination = meta_ord_with_outliers,
  color = "Lifestyle",
  shape = "Gender"
)

outlier_data <- p_with_outliers_base$data %>%
  filter(Lane %in% samples_to_exclude)

p_with_outliers <- p_with_outliers_base +
  geom_point(size = 3.5, alpha = 0.8) +
  geom_point(
    data = outlier_data,
    aes(color = "Outlier"),
    size = 3,
    stroke = 1.3,
    show.legend = TRUE
  ) +
  geom_text_repel(
    data = outlier_data,
    aes(label = Lane),
    color = "#C2185B",
    size = 3,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Rural" = "#72B2A2",
      "Urban" = "#E59572",
      "Outlier" = "#C2185B"
    ),
    labels = c(
      "Rural" = "Rural",
      "Urban" = "Urbano",
      "Outlier" = "Outlier"
    )
  ) +
  scale_shape_manual(
    values = c("Female" = 17, "Male" = 15),
    labels = c("Female" = "Mujer", "Male" = "Hombre"),
    na.value = 16
  ) +
  labs(
    title = "Escalamiento Multidimensional No Metrico (NMDS) con outliers",
    subtitle = paste(
      "Stress del NMDS:", round(meta_ord_with_outliers$stress, 3),
      "- Distancia Bray-Curtis"
    ),
    color = "Grupo",
    shape = "Sexo",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_bw()

print(p_with_outliers)
ggsave(
  file.path(img_dir, 'bracken_woutliers.png'),
  plot = p_with_outliers,
  dpi = 600,
  width = 8,
  height = 6,
  units = "in"
)

#37082_3#16
#37082_2#4
#37035_2#22*
#37035_7#18
#36703_3#20
#36703_3#4

percentages2 <- subset_samples(percentages, !(Lane %in% samples_to_exclude))
#percentages2 <- subset_samples(percentages, !(Lane %in% c("37082_3#16", "37082_2#4")))
meta_ord <- ordinate(physeq = percentages2, method = "NMDS", distance = "bray")
p2 <- plot_ordination(physeq = percentages2, ordination = meta_ord, color = "Lifestyle", label = "Lane", shape="Gender", title = 'Bracken taxonomy without outliers') 
ggsave(file.path(img_dir, 'bracken.png'), dpi=600, units="px", p2)

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
    title = "Escalamiento Multidimensional No-Métrico (NMDS) sin outliers",
    subtitle = paste("Stress del NMDS:", round(meta_ord$stress, 3),
                     "- Distancia Bray-Curtis"),
    color = "Estilo de Vida", # Nuevo tÃ­tulo de leyenda
    shape = "Sexo",           # Nuevo tÃ­tulo de leyenda
    x = "NMDS1",              # Etiqueta de eje
    y = "NMDS2"               # Etiqueta de eje
  ) +
  theme_bw()

print(p_final)
ggsave(file.path(img_dir, 'bracken_final.png'), dpi=600, units="px", p_final)


percentages3 <- subset_taxa(percentages2, Kingdom == "Bacteria")
meta_ord <- ordinate(physeq = percentages3, method = "NMDS", distance = "bray")
p3 <- plot_ordination(physeq = percentages3, ordination = meta_ord, color = "Lifestyle", label = "Lane", shape="Gender", title = 'Taxonomy (Bacteria)') 
ggsave(file.path(img_dir, 'bracken_onlyBacteria.png'), dpi=500, units="px", p3)


percentages3 <- subset_taxa(percentages2, Kingdom == "Archaea")
meta_ord <- ordinate(physeq = percentages3, method = "NMDS", distance = "bray")
p3 <- plot_ordination(physeq = percentages3, ordination = meta_ord, color = "Lifestyle", label = "Lane", shape="Gender", title = 'Taxonomy (Archaea)') 
ggsave(file.path(img_dir, 'bracken_onlyArchaea.png'), dpi=500, units="px", p3)


percentages3 <- subset_taxa(percentages2, Kingdom == "Eukaryota")
meta_ord <- ordinate(physeq = percentages3, method = "NMDS", distance = "bray")
p3 <- plot_ordination(physeq = percentages3, ordination = meta_ord, color = "Lifestyle", label = "Lane", shape="Gender", title = 'Taxonomy (Eukaryota)') 
ggsave(file.path(img_dir, 'bracken_onlyEukaryota.png'), dpi=500, units="px", p3)


percentages3 <- subset_taxa(percentages2, Kingdom == "Viruses")
meta_ord <- ordinate(physeq = percentages3, method = "NMDS", distance = "bray")
p3 <- plot_ordination(physeq = percentages3, ordination = meta_ord, color = "Lifestyle", label = "Lane", shape="Gender", title = 'Taxonomy (Virus)') 
ggsave(file.path(img_dir, 'bracken_onlyVirus.png'), dpi=500, units="px", p3)


TopNGenus <- names(sort(taxa_sums(merged_metagenomes), TRUE)[1:5])
Top5Genus <- prune_taxa(TopNGenus,merged_metagenomes)
plot_heatmap(Top5Genus)
