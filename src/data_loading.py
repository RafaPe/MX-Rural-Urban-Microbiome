"""Lectura validada de tablas de abundancia y metadatos."""

from pathlib import Path

import pandas as pd


def load_table(path: str | Path, *, index_col: int | str | None = 0, **kwargs) -> pd.DataFrame:
    """Carga CSV, TSV o Excel y reporta extensiones no admitidas."""
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path, index_col=index_col, **kwargs)
    if suffix == ".tsv":
        return pd.read_csv(path, sep="\t", index_col=index_col, **kwargs)
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path, index_col=index_col, **kwargs)
    raise ValueError(f"Formato no admitido: {path}")


def validate_sample_ids(abundance: pd.DataFrame, metadata: pd.DataFrame) -> None:
    """Falla temprano si las muestras de abundancia no tienen metadatos."""
    missing = abundance.index.difference(metadata.index)
    if len(missing):
        raise ValueError(f"{len(missing)} muestras no tienen metadatos: {missing[:5].tolist()}")
