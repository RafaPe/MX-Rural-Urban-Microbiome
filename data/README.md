# Datos de los experimentos

El repositorio sólo versiona insumos externos públicos y productos procesados
sin identificadores directos.

- `external/`: fuentes públicas obtenidas de terceros, conservadas sin cambios
  sustantivos.
- `processed/`: matrices y metadatos derivados que consumen las libretas.
- `private/`: datos sensibles disponibles únicamente de forma local e ignorados
  por Git.
- `raw/` e `interim/`: entradas crudas y productos temporales locales, también
  ignorados.

`manifest.csv` documenta la procedencia y las libretas consumidoras. Antes de
hacer público el repositorio se debe revisar nuevamente la autorización de
redistribución de cada archivo. Los libros de encuestas, nutrición, EEG y el
metadato LATINBIOTA completo no forman parte de la versión pública.

La tabla `processed/latinbiota_sample_metadata.csv` se deriva del metadato
privado. Conserva únicamente el identificador técnico de secuenciación y las
categorías requeridas por los experimentos; elimina nombres de proveedor,
identificadores públicos, edad y BMI exactos.
