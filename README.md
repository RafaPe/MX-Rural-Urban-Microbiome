# Gut Microbiome Thesis Analyses

This repository contains the analysis code, notebooks, and selected derived
results for a master's thesis investigating variation in gut microbiome
composition across rural, urban, and publicly available international cohorts.

The project covers taxonomic composition, functional pathways, antimicrobial
resistance genes, genome bins, gut-brain modules (GBMs), alpha and beta
diversity, clustering, co-occurrence networks, and predictive taxonomic
indices. It is a research repository rather than a production pipeline: each
notebook documents a self-contained experiment, while reusable utilities live
in `src/`.

## Repository layout

```text
config/       declarative analysis profiles and shared parameters
data/         dataset documentation and local input locations
notebooks/    numbered analysis notebooks
results/      selected final figures and summary tables
scripts/      standalone R analyses
src/          shared Python utilities
```

## Analyses

| Component | Description |
| --- | --- |
| `01_otu_eda.ipynb` | Exploratory PCA comparison of raw, scaled, and CLR-transformed taxonomic abundances |
| `02_clr_pca_exploration.ipynb` | CLR/PCA exploration of taxonomic and functional profiles, including CAMDA 2026 data |
| `03_bin_exploration.ipynb` | Quality assessment and taxonomic characterization of genome bins |
| `04_clustering_experiments.ipynb` | Clustering and representation optimization for taxonomic and functional profiles |
| `05_resistance_genes_by_group.ipynb` | Comparison of antimicrobial-resistance genes across lifestyle groups |
| `06_gbm_analysis.ipynb` | Gut-brain module analysis |
| `07_taxonomic_indices_predictive_comparison*.ipynb` | Predictive comparison of taxonomic indices, including the BHC variant |
| `scripts/alpha_diversity.R` | Alpha-diversity analysis by lifestyle and age group |
| `scripts/beta_diversity.R` | Beta-diversity ordination and visualization workflow |
| `scripts/networks.R` | Taxonomic co-occurrence network analysis |

Notebook 04 includes extensive Optuna searches and can take minutes or hours,
depending on the selected parameters and available hardware.

## Getting started

Create the Python environment from the repository root:

```bash
conda env create -f environment.yml
conda activate clustering-env
jupyter lab
```

Launch Jupyter from the repository root. The notebooks resolve project paths
when they are started from either the root directory or `notebooks/`. The R
scripts require a local R installation and the packages listed at the top of
each script.

## Data access

The code and selected derived results are publicly available. Taxonomic and
functional datasets—including input profiles, pathway tables, GBM modules,
genome-bin classifications, and resistance annotations—are restricted and are
intentionally not distributed in this repository.

Restricted datasets are available upon request from the repository maintainer,
subject to the appropriate authorization and data-use conditions. The expected
local paths and the analyses that consume each dataset are documented in
[`data/manifest.csv`](data/manifest.csv). Once access is granted, place the
files in the documented locations and keep them out of version control.

Only dataset documentation, public CAMDA metadata, and a deidentified sample
metadata table are versioned under `data/`. Do not commit restricted inputs,
raw data, surveys, full metadata, identifiable information, or regenerable
intermediate files.

## Results and reproducibility

Final figures and summary tables are stored under `results/figures/` and
`results/tables/`. Models, logs, and temporary outputs are reproducible and
ignored by Git.

Where applicable, notebooks declare their random seed and analysis parameters.
Notebooks 04 and 07 use `RANDOM_STATE = 42`; changes to a seed or scientific
parameter should be documented in the notebook and its corresponding commit.
Full numerical reproduction requires the restricted input data and the
software environment defined in `environment.yml`.

## License

See [LICENSE](LICENSE).
