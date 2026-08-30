# Diversidad alfa por estilo de vida y grupo de edad
#
# Ejecutar desde la raiz del proyecto con, por ejemplo:
#   & 'C:/Program Files/R/R-4.5.1/bin/Rscript.exe' scripts/alpha_diversity.R
#
# Salidas:
#   results/figures/alpha_diversity/alpha_diversity_by_lifestyle_age.png
#   results/tables/alpha_diversity/alpha_diversity_values.csv
#   results/tables/alpha_diversity/alpha_diversity_kruskal_wallis.csv
#   results/tables/alpha_diversity/alpha_diversity_lifestyle_by_age_wilcoxon_holm.csv

required_packages <- c("phyloseq", "ggplot2", "readxl", "dplyr", "tidyr", "tibble")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) {
  stop(
    "Faltan los siguientes paquetes de R: ",
    paste(missing_packages, collapse = ", "),
    ". Instalalos y vuelve a ejecutar el script.",
    call. = FALSE
  )
}

library(phyloseq)
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)

find_project_root <- function(start = getwd()) {
  current_dir <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    if (file.exists(file.path(current_dir, "data", "metadata_LATINBIOTA_MEXICO.xlsx"))) {
      return(current_dir)
    }

    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) break
    current_dir <- parent_dir
  }

  stop(
    "No se encontro la raiz del proyecto desde: ", getwd(),
    ". Abre el proyecto y ejecuta el script otra vez.",
    call. = FALSE
  )
}

project_root <- find_project_root()
data_dir <- file.path(project_root, "data")
figures_dir <- file.path(project_root, "results", "figures", "alpha_diversity")
tables_dir <- file.path(project_root, "results", "tables", "alpha_diversity")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

metadata_path <- file.path(data_dir, "metadata_LATINBIOTA_MEXICO.xlsx")
biom_path <- file.path(data_dir, "bracken_taxonomy_final.biom")

if (!file.exists(biom_path)) {
  stop("No se encontro el archivo BIOM en: ", biom_path, call. = FALSE)
}

# Se conservan las decisiones del bloque original de beta_diversity.R.
samples_to_exclude <- c(
  "37082_3#16", "36703_3#4", "37035_7#18",
  "36703_3#20", "37082_2#4", "37035_2#22"
)

# estimate_richness() convierte, por ejemplo, 36703_3#10 en X36703_3.10.
# Esta clave permite unir esos nombres de R con los identificadores originales.
normalise_sample_id <- function(x) {
  x <- sub("^X", "", as.character(x))
  gsub("\\.", "#", x)
}

metadata <- read_excel(metadata_path, sheet = "Data") |>
  select(Lane, Lifestyle, `Age group`) |>
  as.data.frame()

physeq <- import_biom(biom_path)
sample_ids <- sample_names(physeq)

metadata <- metadata[match(sample_ids, metadata$Lane), , drop = FALSE]
if (anyNA(metadata$Lane)) {
  missing_ids <- sample_ids[is.na(metadata$Lane)]
  stop(
    "Hay muestras del BIOM que no aparecen en los metadatos: ",
    paste(missing_ids, collapse = ", "),
    call. = FALSE
  )
}
rownames(metadata) <- sample_ids

metadata <- metadata |>
  rename(Age_group_original = `Age group`) |>
  mutate(
    Age_group = case_when(
      # La figura de referencia agrupaba los registros sin edad como adultos.
      Age_group_original %in% c("Teenager", "Adult") | is.na(Age_group_original) ~ "Adult",
      Age_group_original %in% c("Child", "Infant") ~ "Child",
      Age_group_original == "Elderly" ~ "Elderly",
      TRUE ~ as.character(Age_group_original)
    ),
    Lifestyle = trimws(as.character(Lifestyle)),
    Lifestyle = recode(Lifestyle, "Urban" = "Urbano", "Rural" = "Rural"),
    Age_group = recode(
      Age_group,
      "Child" = "Niño",
      "Adult" = "Adulto",
      "Elderly" = "Adulto mayor"
    ),
    Lifestyle = factor(Lifestyle, levels = c("Rural", "Urbano")),
    Age_group = factor(Age_group, levels = c("Niño", "Adulto", "Adulto mayor"))
  )

sample_data(physeq) <- sample_data(metadata)
physeq <- prune_samples(
  !(normalise_sample_id(sample_names(physeq)) %in% samples_to_exclude),
  physeq
)

# La figura de referencia solo compara Niños y Adultos; se excluye Elderly.
physeq <- prune_samples(
  as.character(sample_data(physeq)$Age_group) != "Adulto mayor",
  physeq
)

# `sample_data` es un objeto S4. Se extraen sus columnas explícitamente para
# obtener un data.frame ordinario, compatible con dplyr en todas las versiones.
sample_metadata <- data.frame(
  SampleID_key = normalise_sample_id(sample_names(physeq)),
  Lifestyle = as.character(sample_data(physeq)$Lifestyle),
  Age_group = as.character(sample_data(physeq)$Age_group),
  stringsAsFactors = FALSE
)

# estimate_richness calcula los mismos tres índices usados antes en plot_richness.
alpha_values <- estimate_richness(physeq, measures = c("Observed", "Chao1", "Shannon")) |>
  tibble::rownames_to_column("SampleID") |>
  mutate(SampleID_key = normalise_sample_id(SampleID)) |>
  left_join(sample_metadata, by = "SampleID_key") |>
  pivot_longer(
    cols = c(Observed, Chao1, Shannon),
    names_to = "Index",
    values_to = "Value"
  ) |>
  mutate(
    Index = factor(Index, levels = c("Observed", "Chao1", "Shannon")),
    Group = interaction(Lifestyle, Age_group, sep = " - ", drop = TRUE)
  )

if (anyNA(alpha_values$Lifestyle) | anyNA(alpha_values$Age_group)) {
  incomplete_samples <- alpha_values |>
    filter(is.na(Lifestyle) | is.na(Age_group)) |>
    distinct(SampleID, Lifestyle, Age_group)
  stop(
    "Hay muestras sin estilo de vida o grupo de edad tras el filtrado: ",
    paste(incomplete_samples$SampleID, collapse = ", "),
    call. = FALSE
  )
}

# Prueba omnibus no parametrica: compara los grupos Lifestyle × Age_group.
kruskal_results <- alpha_values |>
  group_by(Index) |>
  summarise(
    n = n(),
    statistic = unname(kruskal.test(Value ~ Group)$statistic),
    df = unname(kruskal.test(Value ~ Group)$parameter),
    p_value = kruskal.test(Value ~ Group)$p.value,
    .groups = "drop"
  ) |>
  mutate(p_value_holm_3_indices = p.adjust(p_value, method = "holm"))

# Comparaciones planeadas: Rural vs. Urbano por separado en cada grupo de edad
# y para cada índice. Hay seis pruebas; Holm controla el error familiar entre ellas.
lifestyle_by_age_results <- alpha_values |>
  filter(Age_group %in% c("Niño", "Adulto")) |>
  group_by(Index, Age_group) |>
  group_modify(function(data, key) {
    group_counts <- table(data$Lifestyle)
    if (!all(c("Rural", "Urbano") %in% names(group_counts))) {
      stop(
        "Falta Rural o Urbano para la comparación de ",
        as.character(key$Index), " en ", as.character(key$Age_group),
        call. = FALSE
      )
    }

    test <- wilcox.test(Value ~ Lifestyle, data = data, exact = FALSE)
    medians <- tapply(data$Value, data$Lifestyle, median)

    tibble(
      comparison = "Rural vs. Urbano",
      n_rural = unname(group_counts["Rural"]),
      n_urbano = unname(group_counts["Urbano"]),
      median_rural = unname(medians["Rural"]),
      median_urbano = unname(medians["Urbano"]),
      wilcoxon_W = unname(test$statistic),
      p_value = test$p.value
    )
  }) |>
  ungroup() |>
  mutate(
    p_value_holm_6_comparisons = p.adjust(p_value, method = "holm"),
    significance = case_when(
      p_value_holm_6_comparisons < 0.001 ~ "***",
      p_value_holm_6_comparisons < 0.01 ~ "**",
      p_value_holm_6_comparisons < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

group_summary <- alpha_values |>
  group_by(Index, Lifestyle, Age_group) |>
  summarise(
    n = n(),
    median = median(Value),
    q1 = quantile(Value, 0.25),
    q3 = quantile(Value, 0.75),
    .groups = "drop"
  )

write.csv(
  alpha_values,
  file.path(tables_dir, "alpha_diversity_values.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  group_summary,
  file.path(tables_dir, "alpha_diversity_group_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  kruskal_results,
  file.path(tables_dir, "alpha_diversity_kruskal_wallis.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  lifestyle_by_age_results,
  file.path(tables_dir, "alpha_diversity_lifestyle_by_age_wilcoxon_holm.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

age_colours <- c("Niño" = "#F4A300", "Adulto" = "#0066CC")

alpha_plot <- ggplot(alpha_values, aes(x = Lifestyle, y = Value, colour = Age_group)) +
  geom_boxplot(
    aes(group = interaction(Lifestyle, Age_group)),
    position = position_dodge(width = 0.72),
    width = 0.62,
    linewidth = 0.45,
    outlier.shape = NA
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.72
    ),
    size = 1.45,
    alpha = 0.8
  ) +
  facet_wrap(~Index, scales = "free_y") +
  scale_colour_manual(values = age_colours, drop = TRUE) +
  labs(
    title = "Diversidad Alfa por Estilo de Vida y Grupo de Edad",
    x = "Estilo de vida",
    y = "Valor del índice de diversidad",
    colour = "Grupo de edad"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0, face = "plain"),
    strip.background = element_rect(fill = "grey90", colour = "black"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "right"
  )

ggsave(
  filename = file.path(figures_dir, "alpha_diversity_by_lifestyle_age.png"),
  plot = alpha_plot,
  width = 9,
  height = 5.7,
  units = "in",
  dpi = 600,
  bg = "white"
)

message("Figura guardada en: ", file.path(figures_dir, "alpha_diversity_by_lifestyle_age.png"))
message("Pruebas estadisticas guardadas en: ", tables_dir)
print(kruskal_results)
print(lifestyle_by_age_results)
