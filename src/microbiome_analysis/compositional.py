"""Transformaciones para tablas de abundancia composicional."""

import numpy as np
import pandas as pd


def add_remainder_category(data: pd.DataFrame, *, total: float, name: str = "otros") -> pd.DataFrame:
    """Añade la abundancia restante a una categoría, sin modificar la entrada."""
    result = data.copy()
    remainder = total - result.sum(axis=1)
    result[name] = result.get(name, 0) + remainder.clip(lower=0)
    return result


def clr_transform(data: pd.DataFrame, *, pseudocount: float = 1e-6) -> pd.DataFrame:
    """Aplica CLR con pseudoconteo a una matriz de abundancias no negativas."""
    if (data < 0).any().any():
        raise ValueError("CLR requiere abundancias no negativas.")
    values = data.astype(float) + pseudocount
    log_values = np.log(values)
    return log_values.sub(log_values.mean(axis=1), axis=0)
