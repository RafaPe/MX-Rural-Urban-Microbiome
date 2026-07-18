"""Carga de perfiles declarativos para análisis que comparten un notebook."""

from pathlib import Path

import yaml

from .paths import PROJECT_ROOT


def load_profile(name: str, path: Path | None = None) -> dict:
    """Devuelve un perfil y convierte sus rutas a ubicaciones del proyecto."""
    path = path or PROJECT_ROOT / "config" / "analysis_profiles.yaml"
    with path.open(encoding="utf-8") as stream:
        profiles = yaml.safe_load(stream)
    if name not in profiles:
        options = ", ".join(profiles)
        raise KeyError(f"Perfil desconocido: {name}. Opciones: {options}")
    profile = profiles[name].copy()
    for key in ("abundance_path", "metadata_path"):
        profile[key] = PROJECT_ROOT / profile[key]
    return profile
