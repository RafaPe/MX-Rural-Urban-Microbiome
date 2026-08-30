# Redes de coocurrencia del microbioma rural y urbano.
#
# Directorio de trabajo fijo. Si el proyecto se mueve, solo hay que cambiar esta línea.
setwd("C:/Users/rafap/Documents/tesis_maestria/tesis_maestria - copia/scripts")

required_packages <- c("NetCoMi", "SpiecEasi")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    paste0(
      "Faltan paquetes de R necesarios: ",
      paste(missing_packages, collapse = ", "),
      ".\nInstálalos antes de ejecutar este script. Consulta: ",
      "https://github.com/stefpeschel/NetCoMi"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages(library(NetCoMi))

data_dir <- "../datasets"
output_dir <- "../results/figures/08_networks"

abundance_path <- file.path(data_dir, "metaphlan_genus.tsv")
metadata_path <- file.path(data_dir, "metadata_final.csv")

for (path in c(abundance_path, metadata_path)) {
  if (!file.exists(path)) {
    stop("No se encontró el archivo requerido: ", path, call. = FALSE)
  }
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Leyendo abundancias y metadatos...")
abund <- read.delim(
  abundance_path,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)
metadata <- read.csv(
  metadata_path,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

if (anyDuplicated(rownames(abund))) {
  stop("La tabla de abundancias contiene identificadores de muestra duplicados.", call. = FALSE)
}
if (anyDuplicated(rownames(metadata))) {
  stop("Los metadatos contienen identificadores de muestra duplicados.", call. = FALSE)
}
if (!"Lifestyle" %in% colnames(metadata)) {
  stop("Los metadatos no contienen la columna 'Lifestyle'.", call. = FALSE)
}

missing_metadata <- setdiff(rownames(abund), rownames(metadata))
if (length(missing_metadata) > 0L) {
  stop(
    "Faltan metadatos para ", length(missing_metadata),
    " muestras de la tabla de abundancias.",
    call. = FALSE
  )
}

# Reordenar explícitamente evita depender del orden original de ambos archivos.
metadata <- metadata[rownames(abund), , drop = FALSE]
group <- trimws(as.character(metadata$Lifestyle))
keep_samples <- group %in% c("Rural", "Urban")

if (!all(keep_samples)) {
  message("Se excluyeron ", sum(!keep_samples), " muestras sin grupo Rural/Urban.")
  abund <- abund[keep_samples, , drop = FALSE]
  group <- group[keep_samples]
}

abund_rural <- as.matrix(abund[group == "Rural", , drop = FALSE])
abund_urban <- as.matrix(abund[group == "Urban", , drop = FALSE])
storage.mode(abund_rural) <- "double"
storage.mode(abund_urban) <- "double"

if (nrow(abund_rural) == 0L || nrow(abund_urban) == 0L) {
  stop("Se necesitan muestras de los grupos Rural y Urban.", call. = FALSE)
}
if (any(!is.finite(abund_rural)) || any(!is.finite(abund_urban))) {
  stop("La tabla de abundancias contiene valores no numéricos o no finitos.", call. = FALSE)
}

message(
  "Construyendo redes (",
  nrow(abund_rural), " muestras rurales y ",
  nrow(abund_urban), " urbanas)..."
)
set.seed(123456)
net_final <- netConstruct(
  data = abund_rural,
  data2 = abund_urban,
  measure = "spieceasi",
  measurePar = list(
    method = "glasso",
    nlambda = 50,
    lambda.min.ratio = 0.01,
    pulsar.params = list(rep.num = 20)
  ),
  filtTax = "numbSamp",
  filtTaxPar = list(numbSamp = 15),
  verbose = TRUE
)

message("Calculando propiedades de las redes...")
props_final <- netAnalyze(
  net_final,
  centrLCC = TRUE,
  clustMethod = "cluster_louvain"
)

palette <- c(
  "#4A90E2", "#E85D75", "#50C878", "#9B59B6", "#F39C12",
  "#1ABC9C", "#E67E22", "#3498DB", "#95A5A6", "#D35400"
)

common_plot_args <- list(
  x = props_final,
  layout = "spring",
  sameLayout = TRUE,
  layoutGroup = "union",
  repulsion = 0.7,
  shortenLabels = "none",
  labelScale = FALSE,
  labelFont = 3,       # 3 = cursiva en los dispositivos gráficos de R
  hubLabelFont = 3,    # conserva la cursiva también en los nodos hub
  nodeSize = "degree",
  nodeSizeSpread = 4,
  nodeColor = "cluster",
  colorVec = palette,
  sameClustCol = TRUE,
  nodeTransp = 35,
  borderWidth = 2,
  borderCol = "white",
  highlightHubs = TRUE,
  hubBorderWidth = 3,
  hubBorderCol = "#2C3E50",
  edgeWidth = 1,
  negDiffCol = TRUE,
  posCol = "#27AE60",
  negCol = "#E74C3C",
  edgeTranspLow = 20,
  edgeTranspHigh = 0,
  cexNodes = 1.2
)

save_png <- function(path, width, height, draw) {
  grDevices::png(
    filename = path,
    width = width,
    height = height,
    units = "px",
    res = 300,
    bg = "white",
    type = "cairo-png"
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  draw()
}

important_path <- file.path(output_dir, "redes_generos_importantes.png")
message("Generando ", important_path)
save_png(important_path, width = 3600, height = 1800, draw = function() {
  set.seed(123456)
  do.call(plot, c(
    common_plot_args,
    list(
      nodeFilter = "highestDegree",
      nodeFilterPar = 15,
      rmSingles = "inboth",
      cexLabels = 0.75,
      cexHubLabels = 0.75,
      cexTitle = 1.5,
      # Márgenes laterales mayores para separar visualmente ambas redes.
      mar = c(1, 3, 4, 3),
      title1 = "Red rural",
      title2 = "Red urbana"
    )
  ))
})

# NetCoMi fija internamente mfrow = c(1, 2) cuando grafica dos redes. Para la
# figura completa se preparan primero ambos objetos qgraph sin dibujarlos y luego
# se trazan en un dispositivo con dos filas y una columna.
set.seed(123456)
all_networks <- do.call(plot, c(
  common_plot_args,
  list(
    nodeFilter = "none",
    rmSingles = "none",
    cexLabels = 0.40,
    cexHubLabels = 0.40,
    showTitle = FALSE,
    doPlot = FALSE
  )
))

all_path <- file.path(output_dir, "redes_todos_los_generos.png")
message("Generando ", all_path)
save_png(all_path, width = 2400, height = 4000, draw = function() {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(2, 1), mar = c(1, 1, 4, 1), xpd = NA)

  plot(all_networks$q1)
  graphics::title(main = "Red rural: todos los géneros", cex.main = 1.5)

  plot(all_networks$q2)
  graphics::title(main = "Red urbana: todos los géneros", cex.main = 1.5)
})

message("Listo. Figuras guardadas en: ", output_dir)
