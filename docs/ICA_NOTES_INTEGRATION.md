# Integración de apuntes ICA

Fecha: 2026-05-23

## Fuentes revisadas

Los apuntes se revisaron como material interno de apoyo teórico, no como fuente bibliográfica formal del informe.

| Archivo | Páginas | Extracción de texto | Uso en el informe |
|---|---:|---|---|
| `20251210114845334.pdf` | 1 | Sin OCR útil | Pre-ICA, centrado, whitening/sphering, FastICA, convergencia y relación con EEGLAB/MNE. |
| `20251210114832765.pdf` | 1 | Sin OCR útil | Formulación matricial, covarianza, PCA, matriz sphering, datos esferizados y matriz de desmezcla. |
| `20251210114902727.pdf` | 2 | Sin OCR útil | Pasos de ICA en BrainVision Analyzer: sphering, unmixing, control de calidad, evaluación de ICs y back-projection. |
| `Apuntes ICA.pdf` | 3 | Sin OCR útil | BSS, separación de fuentes de ruido/actividad cerebral, topografías, PSD, Infomax/FastICA y reconstrucción limpia. |

## Ideas aprovechadas

- ICA se presenta como separación ciega de fuentes para limpieza y auditoría, no como localización cortical.
- El preprocesado pre-ICA incluye centrado por canal, matriz de covarianza, PCA/SVD y whitening/sphering.
- La matriz de datos se interpreta como canales por muestras, con número de componentes normalmente igual o cercano al número de canales válidos.
- FastICA se formula como estimación de componentes independientes mediante maximización de no gaussianidad y control de convergencia.
- La revisión de ICs combina topografía, PSD, trazas temporales y heurísticos estadísticos.
- La reconstrucción limpia se entiende como back-projection al espacio de canales tras retirar componentes artefactuales.

## Cambios aplicados

- Ampliada la subsección `Fase 2: ICA y limpieza de componentes` en `report/sections_es/methods.tex`.
- Añadida una caja `Decisión metodológica` para separar el papel de ICA como herramienta de limpieza de una interpretación neuroanatómica directa.
- Añadida una tabla operativa sobre centrado, whitening/sphering, FastICA y back-projection.
- Añadida una caja `Control QC` con supuestos y cautelas: mezcla lineal instantánea, independencia aproximada, no gaussianidad e indeterminación de escala/signo/orden.

## Cautelas

- Los PDF son escaneos sin OCR; la lectura se hizo visualmente a partir de páginas renderizadas en `/private/tmp/neuro_ica_notes/`.
- No se copiaron imágenes manuscritas al informe para mantener estética homogénea y evitar depender de material escaneado.
- No se modificaron algoritmos Julia ni parámetros de ejecución.
