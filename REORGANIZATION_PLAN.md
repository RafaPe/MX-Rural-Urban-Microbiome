# Plan de reorganización y reproducibilidad

## Objetivo

Mantener una ruta analítica corta y reproducible para la tesis, separar los
resultados de los datos y conservar todos los experimentos históricos sin
mezclarlos con el trabajo activo.

## Estructura objetivo

```text
notebooks/
  01_otu_eda.ipynb
  02_clr_pca_exploration.ipynb
  03_bin_exploration.ipynb
  04_clustering_experiments.ipynb
  05_resistance_genes_by_group.ipynb
  06_gbm_analysis.ipynb
  07_taxonomic_indices_predictive_comparison.ipynb
notebooks_archive/2026-07-pre-refactor/
  duplicates/
  legacy/
  exploratory/
src/
  paths.py
  data_loading.py
  compositional.py
  metadata.py
config/
  analysis_profiles.yaml
results/
  figures/<analysis>/
  tables/<analysis>/
  models/<analysis>/
  logs/
  tmp/
```

## Decisiones aplicadas

1. Los siete notebooks priorizados son los canónicos y se ordenan con prefijos
   numéricos.
2. Las variantes, copias y análisis no priorizados se conservan bajo
   `notebooks_archive/2026-07-pre-refactor`; no se eliminan.
3. Las dos variantes de MetaPhlAn se consolidan en
   `02_clr_pca_exploration.ipynb`. El archivo base procede de
   `metaphlan_EDA_genus-genus.ipynb`; la variante `Copy2` queda archivada.
4. Las diferencias de fuente, total de abundancia y metadatos se expresan como
   perfiles en `config/analysis_profiles.yaml`, no como copias de notebook.
5. Las funciones reutilizables de carga, composición y metadatos viven en
   `src/`. Los notebooks deben usarlas y contener solamente narrativa,
   parámetros, análisis y visualizaciones.
6. Los datos fuente se leen exclusivamente de `datasets/`; las salidas se
   escriben en `results/`, separando figuras, tablas, modelos, registros y
   temporales.

## Convención de salidas

Cada notebook define `ANALYSIS` con su nombre de archivo y crea sus directorios
mediante `src.paths.analysis_output_dirs(ANALYSIS)`. Las figuras se guardan en
`results/figures/<analysis>/`; los CSV/TSV/XLSX/JSON en
`results/tables/<analysis>/`; los modelos en `results/models/<analysis>/`.

Las salidas regenerables y pesadas pertenecen en `results/tmp/` y no deben
versionarse. Las figuras y tablas finales se revisan antes de versionarse.

## Refactorización de MetaPhlAn

El notebook canónico debe seleccionar un perfil con
`load_profile("metaphlan_genus_percentage")` o
`load_profile("metaphlan_genus_relative")`. `load_table`,
`add_remainder_category`, `clr_transform` y `merge_metadata` sustituyen las
celdas repetidas. Cualquier nuevo experimento se añade como un perfil y se
documenta, en lugar de clonarse el notebook.

## Verificación pendiente

Antes de retirar los archivos archivados de una rama futura, ejecutar cada
notebook canónico desde cero, comparar sus tablas y figuras con la versión
histórica, y registrar las diferencias metodológicas en el propio notebook.
