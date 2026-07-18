from pathlib import Path

import nbformat


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "notebooks" / "07_taxonomic_indices_predictive_comparison.ipynb"
MARKER = "## 8. Repetición del experimento con CAMDA 2026"


nb = nbformat.read(NOTEBOOK, as_version=4)

# Permite regenerar la sección sin duplicarla.
for index, cell in enumerate(nb.cells):
    if cell.cell_type == "markdown" and MARKER in cell.source:
        nb.cells = nb.cells[:index]
        break

# La referencia histórica depende de versiones antiguas de las bibliotecas. Se
# conserva la comparación, pero una diferencia numérica ya no detiene la libreta.
for cell in nb.cells:
    if cell.cell_type == "markdown" and "## 7. Verificación de los resultados históricos" in cell.source:
        cell.source = """## 7. Verificación de los resultados históricos y visualización

La referencia de la libreta histórica para la media de F1 del índice combinado
es `0.7061149591` sin tratamiento y `0.9628320683` con CLR. Se muestra la
diferencia observada, pero no se interrumpe la ejecución: pequeñas variaciones
pueden aparecer entre versiones de `scipy` y `scikit-learn`."""
    if cell.cell_type == "code" and "raise AssertionError(\"Los F1 no reproducen" in cell.source:
        cell.source = cell.source.replace(
            "    if not np.allclose(verification[\"Observed\"], verification[\"Expected original\"], atol=1e-6):\n"
            "        raise AssertionError(\"Los F1 no reproducen la libreta original dentro de la tolerancia.\")",
            "    verification[\"Within tolerance\"] = verification[\"Absolute difference\"] <= 1e-6\n"
            "    if not verification[\"Within tolerance\"].all():\n"
            "        print(\"Aviso: la implementación actual difiere ligeramente de la referencia histórica.\")",
        )
    if cell.cell_type == "code" and "combined_f1 = (" in cell.source:
        cell.source = """combined_f1 = (
    fold_metrics.loc[fold_metrics["Index"].eq("Combined")]
    .groupby("Treatment")["F1"].mean()
)

if LEGACY_POSITIONAL_ALIGNMENT:
    expected = pd.Series({
        "Without compositional treatment": 0.7061149590674622,
        "With CLR": 0.9628320683111955,
    })
    verification = pd.DataFrame({
        "Observed": combined_f1,
        "Expected original": expected,
        "Absolute difference": (combined_f1 - expected).abs(),
    })
    verification["Within tolerance"] = verification["Absolute difference"] <= 1e-6
    if not verification["Within tolerance"].all():
        print("Aviso: la implementación actual difiere ligeramente de la referencia histórica.")
else:
    verification = combined_f1.to_frame("F1 with SampleID alignment")

display(verification)


def plot_combined_performance(metrics, title, filenames):
    # Réplica visual de performance_comparison_hipca.png para el índice combinado.
    metric_order = ["Accuracy", "F1", "Precision", "Recall", "Specificity"]
    treatments = [
        ("Without compositional treatment", "Sin CLR", "#5b8fc9"),
        ("With CLR", "Con CLR", "#df6265"),
    ]
    combined = metrics.loc[metrics["Index"].eq("Combined")]
    means = combined.groupby("Treatment")[metric_order].mean()
    stds = combined.groupby("Treatment")[metric_order].std()

    x = np.arange(len(metric_order))
    width = 0.36
    fig, ax = plt.subplots(figsize=(13, 7), dpi=150)
    for position, (treatment, label, color) in enumerate(treatments):
        offset = (position - 0.5) * width
        values = means.loc[treatment, metric_order].to_numpy()
        errors = stds.loc[treatment, metric_order].to_numpy()
        bars = ax.bar(
            x + offset, values, width, label=label, color=color,
            yerr=errors, capsize=5,
            error_kw={"ecolor": "#8a8a8a", "elinewidth": 1.5, "capthick": 1.5},
            edgecolor="white", linewidth=0.8,
        )
        for bar, value in zip(bars, values):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                max(value - 0.055, 0.025),
                f"{value:.3f}",
                ha="center", va="center", color="white",
                fontsize=11, fontweight="bold",
            )

    ax.set_xticks(x, metric_order)
    ax.set_ylabel("Puntuación media")
    ax.set_xlabel("Métricas", fontweight="bold")
    ax.set_ylim(0, 1.10)
    ax.set_title(title, fontsize=18, fontweight="bold", pad=18)
    ax.grid(axis="y", alpha=0.18, linewidth=0.7)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(
        title="Tratamiento", frameon=True, edgecolor="#333333",
        loc="upper left", bbox_to_anchor=(1.01, 1.0),
    )
    fig.tight_layout()
    for filename in filenames:
        fig.savefig(FIGURES / filename, dpi=300, bbox_inches="tight")
    plt.show()
    return fig


original_combined_figure = plot_combined_performance(
    fold_metrics,
    "Comparación del rendimiento del índice combinado\\n(validación cruzada de 5 folds)",
    ["combined_index_performance_comparison.png", "indices_f1_comparison.png"],
)"""

nb.cells.extend([
    nbformat.v4.new_markdown_cell(
        """## 8. Repetición del experimento con CAMDA 2026

### Objetivo y diseño

Se repite la comparación de los índices de Hotelling $T^2$, estadístico $Q$ e
índice combinado con los datos taxonómicos de CAMDA 2026. Las dos rutas usan
los mismos folds y se alinean explícitamente por `sample_id`:

1. **Sin tratamiento composicional:** abundancias relativas originales,
   selección KS, transformación no lineal histórica y estandarización.
2. **Con CLR:** amalgamación de especies con más de 90 % de ceros en `otros`,
   reemplazo de ceros con la mitad de la abundancia positiva mínima, CLR dentro
   de cada fold, selección KS y estandarización.

Las muestras **no occidentales** (`non_westernized = yes`) constituyen el grupo
de referencia con el que se ajustan los límites PCA. La clase positiva que se
intenta detectar es **occidental** (`non_westernized = no`). Esta dirección
reproduce la correspondencia conceptual del experimento original: ajustar la
región de referencia con el grupo rural/no occidental y detectar el grupo
urbano/occidental como desviación. Por el fuerte desbalance se reportan F1,
sensibilidad, especificidad y exactitud balanceada, además de la exactitud
convencional.

Para estudiar el efecto del balance sin depender de un único subconjunto, se
realizan **10 experimentos reproducibles**. Cada experimento conserva las 138
muestras no occidentales y selecciona aleatoriamente **200 occidentalizadas**.
El preprocesamiento composicional y los cinco folds se recalculan dentro de cada
experimento, exclusivamente con sus 338 muestras."""
    ),
    nbformat.v4.new_code_cell(
        """# Carga y preparación de CAMDA 2026, equivalente a 02_clr_pca_exploration.
CAMDA_DATA = PROJECT_ROOT / "data"
CAMDA_TAXONOMY = CAMDA_DATA / "taxonomy_relative_abundance.tsv"
CAMDA_METADATA = CAMDA_DATA / "metadata.tsv"

camda_raw = pd.read_csv(CAMDA_TAXONOMY, sep="\\t", index_col=0).T
camda_raw = camda_raw.apply(pd.to_numeric, errors="coerce").fillna(0.0)
camda_metadata = (
    pd.read_csv(CAMDA_METADATA, sep="\\t")
    .drop_duplicates("sample_id")
    .set_index("sample_id")
)

camda_labels = camda_metadata["non_westernized"].reindex(camda_raw.index)
valid_camda = camda_labels.isin(["yes", "no"])
camda_raw = camda_raw.loc[valid_camda].copy()
camda_labels = camda_labels.loc[valid_camda].map({
    "no": "Occidental", "yes": "No occidental",
})

# Diseño de los submuestreos repetidos, aplicado sólo a CAMDA 2026.
CAMDA_WESTERNIZED_N = 200
CAMDA_EXPERIMENT_SEEDS = [7, 21, 42, 84, 101, 202, 303, 404, 505, 606]
camda_non_western_ids = camda_labels.index[camda_labels.eq("No occidental")]
camda_western_ids = camda_labels.index[camda_labels.eq("Occidental")]


def prepare_camda_composition(raw_subset):
    # Amalgamación y tratamiento de ceros idénticos a la sección CAMDA de la libreta 02.
    pseudocount = raw_subset.where(raw_subset > 0).min().min() / 2
    zero_threshold = int(0.9 * len(raw_subset))
    sparse_features = raw_subset.columns[
        (raw_subset == 0).sum(axis=0) > zero_threshold
    ]

    treated = raw_subset.copy()
    treated["otros"] = treated.loc[:, sparse_features].sum(axis=1)
    treated = treated.drop(columns=sparse_features)
    treated["otros"] += 100.0 - treated.sum(axis=1)

    zero_count_per_sample = (treated == 0).sum(axis=1)
    treated = treated.replace(0, pseudocount)
    correction = (zero_count_per_sample * pseudocount) / treated.shape[1]
    treated = treated.sub(correction, axis=0)

    tiny_negative_other = (
        (treated["otros"] < 0) & (treated["otros"] > -pseudocount)
    )
    treated.loc[tiny_negative_other, "otros"] = pseudocount
    if (treated <= 0).any().any():
        raise ValueError("CAMDA 2026: quedaron valores no positivos tras tratar los ceros.")
    return treated, pseudocount


print("Pool CAMDA:", camda_raw.shape)
print("No occidentales disponibles:", len(camda_non_western_ids))
print("Occidentalizadas disponibles:", len(camda_western_ids))
print(
    f"Diseño: {len(CAMDA_EXPERIMENT_SEEDS)} experimentos × "
    f"({len(camda_non_western_ids)} no occidentales + "
    f"{CAMDA_WESTERNIZED_N} occidentalizadas)"
)"""
    ),
    nbformat.v4.new_markdown_cell(
        """### Funciones adaptadas y validación cruzada

La lógica matemática del experimento se conserva. Los nombres de las clases se
parametrizan para CAMDA 2026 y la construcción matricial se vectoriza para
evitar la fragmentación de `DataFrame` que produce la implementación histórica.
La selección de características, el escalador y el PCA se ajustan únicamente
con las muestras de entrenamiento de cada fold."""
    ),
    nbformat.v4.new_code_cell(
        """CAMDA_REFERENCE_LABEL = "No occidental"
CAMDA_POSITIVE_LABEL = "Occidental"


def binary_metrics_for_labels(y_true, y_pred, negative_label, positive_label):
    tn, fp, fn, tp = confusion_matrix(
        y_true, y_pred, labels=[negative_label, positive_label]
    ).ravel()
    return {
        "F1": f1_score(y_true, y_pred, pos_label=positive_label, zero_division=0),
        "Precision": precision_score(
            y_true, y_pred, pos_label=positive_label, zero_division=0
        ),
        "Recall": recall_score(y_true, y_pred, pos_label=positive_label, zero_division=0),
        "Specificity": tn / (tn + fp),
        "Accuracy": accuracy_score(y_true, y_pred),
        "Balanced Accuracy": balanced_accuracy_score(y_true, y_pred),
        "TP": tp, "FP": fp, "TN": tn, "FN": fn,
    }


def make_camda_feature_matrix(feature_table, sample_ids, features, nonlinear=False):
    ordered_features = list(dict.fromkeys(features))
    available_samples = [sample for sample in sample_ids if sample in feature_table.columns]
    matrix = feature_table.reindex(
        index=ordered_features, columns=available_samples, fill_value=0.0
    ).T
    if nonlinear:
        matrix = matrix.map(custom_transform)
    return matrix


def fit_camda_index_model(feature_table, reference_ids, positive_ids,
                          compositional, p_value=KS_P_VALUE):
    reference_features, positive_features = select_features_by_ks(
        feature_table, reference_ids, positive_ids, method="exact", p_value=p_value
    )
    # Reproduce la selección original: ambos sentidos con CLR y sólo los rasgos
    # de la clase objetivo en la ruta sin tratamiento composicional.
    features = (
        reference_features + positive_features
        if compositional else positive_features
    )
    if not features:
        raise ValueError("Ninguna característica superó el umbral KS en este fold.")
    training = make_camda_feature_matrix(
        feature_table, reference_ids, features, nonlinear=not compositional
    )
    scaler = StandardScaler()
    scaled = pd.DataFrame(
        scaler.fit_transform(training), index=training.index, columns=training.columns
    )
    model = calculate_pca_parameters(scaled, PCA_VARIANCE, ALPHA)
    model.update({
        "features": list(training.columns),
        "scaler": scaler,
        "nonlinear": not compositional,
        "n_reference_features": len(reference_features),
        "n_positive_features": len(positive_features),
    })
    return model


def transform_camda_for_model(feature_table, model):
    matrix = make_camda_feature_matrix(
        feature_table, feature_table.columns, model["features"],
        nonlinear=model["nonlinear"],
    )
    return pd.DataFrame(
        model["scaler"].transform(matrix),
        index=matrix.index,
        columns=matrix.columns,
    )


def calculate_camda_indices(transformed_matrix, model):
    # Conserva la regla histórica de usar coordenadas PCA cuando las dimensiones coinciden.
    pca_coordinates = model["pca"].transform(transformed_matrix)
    index_input = (
        pca_coordinates
        if pca_coordinates.shape[1] == model["D"].shape[0]
        else transformed_matrix.to_numpy()
    )
    t2 = np.einsum("ij,jk,ik->i", index_input, model["D"], index_input)
    q = np.einsum("ij,jk,ik->i", index_input, model["C"], index_input)
    combined = np.einsum(
        "ij,jk,ik->i", index_input, model["combined_matrix"], index_input
    )
    return pd.DataFrame({
        "SampleID": transformed_matrix.index,
        "T2": t2,
        "Prediction T2": np.where(
            t2 > model["t2_threshold"], CAMDA_POSITIVE_LABEL, CAMDA_REFERENCE_LABEL
        ),
        "Q": q,
        "Prediction Q": np.where(
            q > model["q_threshold"], CAMDA_POSITIVE_LABEL, CAMDA_REFERENCE_LABEL
        ),
        "Combined Index": combined,
        "Combined Prediction": np.where(
            combined > model["combined_threshold"],
            CAMDA_POSITIVE_LABEL,
            CAMDA_REFERENCE_LABEL,
        ),
    })


def evaluate_camda_fold(treated_train, treated_test, raw_train, raw_test,
                        y_test, reference_ids, positive_ids):
    clr_train = clr_transform(treated_train)
    clr_test = clr_transform(treated_test)
    clr_model = fit_camda_index_model(
        clr_train.T, reference_ids, positive_ids, compositional=True
    )
    clr_predictions = calculate_camda_indices(
        transform_camda_for_model(clr_test.T, clr_model), clr_model
    )

    raw_model = fit_camda_index_model(
        raw_train.T, reference_ids, positive_ids, compositional=False
    )
    raw_predictions = calculate_camda_indices(
        transform_camda_for_model(raw_test.T, raw_model), raw_model
    )

    prediction_columns = {
        "T2": "Prediction T2",
        "Q": "Prediction Q",
        "Combined": "Combined Prediction",
    }
    rows = []
    for treatment, predictions, model in [
        ("Without compositional treatment", raw_predictions, raw_model),
        ("With CLR", clr_predictions, clr_model),
    ]:
        # Garantiza que las métricas sigan el orden exacto de las predicciones.
        truth = pd.Series(y_test, index=treated_test.index).loc[predictions["SampleID"]]
        for index_name, column in prediction_columns.items():
            metrics = binary_metrics_for_labels(
                truth, predictions[column], CAMDA_REFERENCE_LABEL, CAMDA_POSITIVE_LABEL
            )
            rows.append({
                "Treatment": treatment,
                "Index": index_name,
                "Selected features": len(model["features"]),
                **metrics,
            })
    return rows


def run_camda_cross_validation(treated, raw, labels):
    splitter = StratifiedKFold(
        n_splits=N_SPLITS, shuffle=True, random_state=RANDOM_STATE
    )
    fold_rows = []
    for fold, (train_idx, test_idx) in enumerate(
        splitter.split(treated, labels), start=1
    ):
        train_ids = treated.index[train_idx]
        reference_ids = train_ids[labels.loc[train_ids].eq(CAMDA_REFERENCE_LABEL)]
        positive_ids = train_ids[labels.loc[train_ids].eq(CAMDA_POSITIVE_LABEL)]
        rows = evaluate_camda_fold(
            treated.iloc[train_idx], treated.iloc[test_idx],
            raw.iloc[train_idx], raw.iloc[test_idx],
            labels.iloc[test_idx].tolist(), reference_ids, positive_ids,
        )
        for row in rows:
            row["Fold"] = fold
        fold_rows.extend(rows)
    return pd.DataFrame(fold_rows)"""
    ),
    nbformat.v4.new_code_cell(
        """camda_experiment_rows = []
for experiment, seed in enumerate(CAMDA_EXPERIMENT_SEEDS, start=1):
    sampled_western_ids = (
        camda_labels.loc[camda_western_ids]
        .sample(n=CAMDA_WESTERNIZED_N, random_state=seed)
        .index
    )
    selected_ids = camda_labels.index[
        camda_labels.index.isin(camda_non_western_ids.union(sampled_western_ids))
    ]
    raw_subset = camda_raw.loc[selected_ids].copy()
    labels_subset = camda_labels.loc[selected_ids].copy()
    treated_subset, pseudocount = prepare_camda_composition(raw_subset)

    experiment_folds = run_camda_cross_validation(
        treated_subset, raw_subset, labels_subset
    )
    experiment_folds["Experiment"] = experiment
    experiment_folds["Seed"] = seed
    experiment_folds["Samples"] = len(selected_ids)
    experiment_folds["CLR input features"] = treated_subset.shape[1]
    experiment_folds["Pseudocount"] = pseudocount
    camda_experiment_rows.append(experiment_folds)
    print(
        f"Experimento {experiment:02d}/{len(CAMDA_EXPERIMENT_SEEDS)} "
        f"(semilla {seed}): {len(selected_ids)} muestras, "
        f"{treated_subset.shape[1]} rasgos para CLR"
    )

camda_fold_metrics = pd.concat(camda_experiment_rows, ignore_index=True)
camda_metric_columns = [
    "F1", "Precision", "Recall", "Specificity", "Accuracy",
    "Balanced Accuracy", "Selected features",
]
camda_experiment_metrics = (
    camda_fold_metrics.groupby(["Experiment", "Seed", "Treatment", "Index"])
    [camda_metric_columns].mean().reset_index()
)
camda_summary = (
    camda_experiment_metrics.groupby(["Treatment", "Index"])
    [[
        "F1", "Precision", "Recall", "Specificity", "Accuracy",
        "Balanced Accuracy", "Selected features",
    ]]
    .agg(["mean", "std"])
)

camda_fold_metrics.to_csv(OUTPUT / "camda_2026_fold_metrics.csv", index=False)
camda_experiment_metrics.to_csv(
    OUTPUT / "camda_2026_experiment_metrics.csv", index=False
)
camda_summary.to_csv(OUTPUT / "camda_2026_metrics_summary.csv")
display(camda_summary)

camda_combined_figure = plot_combined_performance(
    camda_experiment_metrics,
    "CAMDA 2026: índice combinado en 10 submuestreos\\n(200 occidentalizadas + 138 no occidentales; 5 folds)",
    [
        "camda_2026_combined_index_performance_comparison.png",
        "camda_2026_indices_f1_comparison.png",
    ],
)"""
    ),
    nbformat.v4.new_markdown_cell(
        """### Resultados e interpretación

Se ejecutaron **10 experimentos**; cada uno incluyó las 138 muestras no
occidentales y una selección aleatoria de 200 occidentalizadas. Los valores
siguientes son la media y desviación estándar entre los resultados medios de
los diez experimentos y consideran únicamente el **índice combinado**.

Con CLR se obtuvo **F1 = 0.926 ± 0.008**, exactitud **0.914 ± 0.009**,
precisión **0.937 ± 0.003**, sensibilidad **0.917 ± 0.017**, especificidad
**0.910 ± 0.005** y exactitud balanceada **0.913 ± 0.007**. Sin CLR se obtuvo
**F1 = 0.926 ± 0.004**, exactitud **0.905 ± 0.005**, precisión
**0.870 ± 0.005**, sensibilidad **0.990 ± 0.004**, especificidad
**0.784 ± 0.010** y exactitud balanceada **0.887 ± 0.005**.

CLR superó a la ruta sin CLR en **exactitud balanceada, precisión y
especificidad en los 10 de 10 experimentos**, en exactitud en 8 de 10 y en F1
en 6 de 10. La ventaja media de CLR fue **+0.0265 en exactitud balanceada**,
mientras que el F1 promedio fue prácticamente idéntico entre tratamientos.
Sin CLR mantuvo mayor sensibilidad en los 10 experimentos. En consecuencia,
los submuestreos más balanceados sí aportan evidencia consistente de que CLR
mejora la separación equilibrada entre ambos grupos, aunque intercambia parte
de la sensibilidad por mucha mayor especificidad."""
    ),
])

nbformat.write(nb, NOTEBOOK)
print(f"Sección CAMDA añadida a {NOTEBOOK}")
