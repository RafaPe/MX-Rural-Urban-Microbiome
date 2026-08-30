"""Armonización explícita de metadatos de las muestras."""

import pandas as pd


def merge_metadata(abundance: pd.DataFrame, metadata: pd.DataFrame, *, sample_id: str = "Lane") -> pd.DataFrame:
    """Une metadatos al índice de muestras conservando su orden."""
    metadata = metadata.copy()
    if sample_id in metadata.columns:
        metadata = metadata.set_index(sample_id)
    return abundance.join(metadata, how="left")


def recode_age_group(values: pd.Series) -> pd.Series:
    """Agrupa Infant/Child y Teenager/Adult de forma documentada."""
    mapping = {"Infant": "Child", "Child": "Child", "Teenager": "Adult", "Adult": "Adult"}
    return values.map(mapping).fillna(values)
