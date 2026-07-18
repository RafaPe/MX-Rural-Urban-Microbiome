"""Rutas reproducibles del proyecto y de sus productos analíticos."""

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"
DATASETS_DIR = DATA_DIR / "processed"
EXTERNAL_DATA_DIR = DATA_DIR / "external"
RESULTS_DIR = PROJECT_ROOT / "results"


def analysis_output_dirs(analysis: str) -> dict[str, Path]:
    """Crea y devuelve los directorios de salida para un análisis."""
    locations = {
        "figures": RESULTS_DIR / "figures" / analysis,
        "tables": RESULTS_DIR / "tables" / analysis,
        "models": RESULTS_DIR / "models" / analysis,
        "logs": RESULTS_DIR / "logs",
        "tmp": RESULTS_DIR / "tmp" / analysis,
    }
    for location in locations.values():
        location.mkdir(parents=True, exist_ok=True)
    return locations
