# Experimentos de microbioma intestinal

Este repositorio reúne experimentos de tesis para caracterizar y comparar
poblaciones rurales, urbanas y conjuntos públicos internacionales mediante
abundancias taxonómicas, rutas funcionales, genes de resistencia y módulos del
microbioma intestinal.

No es un pipeline de producción: cada libreta representa un experimento y puede
abrirse desde la raíz del repositorio. El código reutilizable se concentra en
`src/`, los insumos permitidos en `data/` y las salidas en `results/`.

## Experimentos

| Libreta | Experimento |
| --- | --- |
| `01_otu_eda.ipynb` | Comparación PCA de abundancias OTU crudas, escaladas y CLR |
| `02_clr_pca_exploration.ipynb` | PCA composicional de taxonomía, funciones y CAMDA 2026 |
| `03_bin_exploration.ipynb` | Calidad y clasificación taxonómica de bins |
| `04_clustering_experiments.ipynb` | Clustering y optimización de representaciones taxonómicas y funcionales |
| `05_resistance_genes_by_group.ipynb` | Comparación de genes de resistencia entre estilos de vida |
| `06_gbm_analysis.ipynb` | Análisis de módulos del microbioma intestinal |
| `07_taxonomic_indices_predictive_comparison.ipynb` | Evaluación predictiva de índices PCA con y sin CLR |

La libreta 04 contiene búsquedas de Optuna extensas y puede tardar varios
minutos u horas según los parámetros y el equipo.

## Instalación

```bash
conda env create -f environment.yml
conda activate clustering-env
jupyter lab
```

Abre Jupyter desde la raíz del repositorio. Las libretas también detectan cuando
se ejecutan desde `notebooks/` y resuelven las rutas del proyecto.

## Estructura

```text
config/       parámetros declarativos de análisis
data/         insumos externos públicos y tablas procesadas desidentificadas
notebooks/    experimentos numerados
results/      figuras y tablas finales
scripts/      análisis auxiliares que no viven en una libreta
src/          funciones compartidas por los experimentos
```

Consulta `data/README.md` y `data/manifest.csv` antes de añadir un nuevo
dataset. Los archivos crudos, encuestas, metadatos completos y cualquier dato
potencialmente identificable no deben versionarse.

## Resultados y reproducibilidad

Cada experimento escribe en una subcarpeta con su nombre dentro de
`results/figures/` o `results/tables/`. Se conservan las figuras finales y las
tablas resumen necesarias para interpretar la tesis; temporales, modelos y logs
se ignoran.

Los experimentos aleatorios declaran una semilla en su bloque inicial de
parámetros. Cambiar una semilla o un parámetro científico debe quedar registrado
en la libreta y en el commit correspondiente.

## Privacidad

Los identificadores de las tablas públicas son identificadores técnicos de
secuenciación. Aun así, antes de publicar una versión del repositorio debe
revisarse la autorización de redistribución de cada entrada del manifiesto. La
eliminación de un archivo en la rama actual no lo elimina automáticamente del
historial anterior de Git.
