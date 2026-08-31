# Data access

This directory documents the input locations used by the analyses. It does not
distribute restricted research inputs.

## Versioned files

The repository tracks only the following data-related files:

- `manifest.csv`, which documents required inputs, their access level, and the
  analyses that use them;
- `external/camda_2026_metadata.tsv`, a public CAMDA 2026 metadata file; and
- `processed/latinbiota_sample_metadata.csv`, a deidentified table containing
  only the technical sequencing identifier and the categories required by the
  notebooks.

## Restricted inputs

Taxonomic profiles, functional pathway tables, GBM modules, genome-bin
classifications, resistance annotations, and related derived matrices are local
only and ignored by Git. They are available upon request from the repository
maintainer, subject to the appropriate authorization and data-use conditions.

After receiving authorization, place the files at the paths listed in
`manifest.csv`. Do not commit them, raw data, surveys, full metadata, or any
potentially identifiable information. The `raw/`, `interim/`, and `private/`
directories are likewise intended for local use only.
