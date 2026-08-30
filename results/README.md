# Resultados

Las libretas escriben aquí sus productos derivados, nunca dentro de `notebooks/`
ni de `data/`.

- `figures/`: figuras finales en PNG, SVG o PDF.
- `tables/`: tablas resumen necesarias para interpretar o continuar un
  experimento.
- `models/`, `logs/` y `tmp/`: productos regenerables ignorados por Git.

Cada experimento usa una subcarpeta con el nombre de su libreta, por ejemplo
`02_clr_pca_exploration`.

No se deben versionar copias completas de los datos de entrada, tablas que
contengan identificadores sensibles ni resultados intermedios que la propia
libreta pueda regenerar. Una tabla consumida como entrada por un experimento
independiente pertenece en `data/processed/`, no en `results/`.
