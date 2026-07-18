"""Añade una cabecera reproducible a cada notebook canónico.

Se ejecuta de forma idempotente y utiliza únicamente la librería estándar para
no depender de Jupyter en tareas de mantenimiento.
"""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOKS = {
    "01_otu_eda.ipynb": "01_otu_eda",
    "02_clr_pca_exploration.ipynb": "02_clr_pca_exploration",
    "03_bin_exploration.ipynb": "03_bin_exploration",
    "04_clustering_experiments.ipynb": "04_clustering_experiments",
    "05_resistance_genes_by_group.ipynb": "05_resistance_genes_by_group",
    "06_gbm_analysis.ipynb": "06_gbm_analysis",
}
MARKER = "# Canonical notebook setup"


def source(text: str) -> list[str]:
    return [line + "\n" for line in text.strip().splitlines()]


def prepare(path: Path, analysis: str) -> bool:
    notebook = json.loads(path.read_text(encoding="utf-8"))
    if any(MARKER in "".join(cell.get("source", [])) for cell in notebook["cells"]):
        return False
    markdown = {
        "cell_type": "markdown",
        "metadata": {},
        "source": source(
            f"# {analysis.replace('_', ' ').title()}\n\n"
            "Notebook canónico. Los datos se leen de `datasets/` y los productos "
            "derivados se guardan en `results/`."
        ),
    }
    code = {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": source(
            f"{MARKER}\n"
            "from pathlib import Path\n"
            "import sys\n\n"
            "PROJECT_ROOT = Path.cwd().resolve()\n"
            "if PROJECT_ROOT.name == 'notebooks':\n"
            "    PROJECT_ROOT = PROJECT_ROOT.parent\n"
            "if str(PROJECT_ROOT) not in sys.path:\n"
            "    sys.path.insert(0, str(PROJECT_ROOT))\n\n"
            "from src.paths import DATASETS_DIR, analysis_output_dirs\n\n"
            f"ANALYSIS = '{analysis}'\n"
            "OUTPUTS = analysis_output_dirs(ANALYSIS)\n"
            "# Guardar figuras nuevas en OUTPUTS['figures'] y tablas en OUTPUTS['tables']."
        ),
    }
    notebook["cells"] = [markdown, code, *notebook["cells"]]
    path.write_text(json.dumps(notebook, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    return True


if __name__ == "__main__":
    updated = [name for name, analysis in NOTEBOOKS.items() if prepare(ROOT / "notebooks" / name, analysis)]
    print(f"Notebooks updated: {', '.join(updated) or 'none'}")
