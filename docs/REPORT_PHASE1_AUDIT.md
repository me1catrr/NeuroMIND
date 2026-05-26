# Auditoría Fase 1 — Integración informe, dashboard y resultados NeuroMIND

Fecha: 2026-05-23  
Ejecución auditada: `results/subjects/sub-M05/ses-T2/eyesclosed/`  
Informe destino: `report/main_es.tex`  
PDF principal: `report/build/pdf/main_es.pdf`

## Propósito

Esta auditoría prepara la ampliación del informe oficial de NeuroMIND. Su objetivo es mapear, de forma trazable, qué produce el pipeline Julia, qué muestra el dashboard, qué debe entrar en el informe LaTeX y qué artefactos deben permanecer como material auxiliar.

No introduce conclusiones clínicas. La ejecución revisada corresponde a un único sujeto/sesión/condición y debe tratarse como validación técnica del pipeline.

## Estado posterior de ejecución

Esta auditoría ya fue usada como guía de implementación. Las fases derivadas quedaron ejecutadas así:

- Fase 2: apéndice técnico del dashboard incorporado al informe.
- Fase 3: metodología ampliada con contratos operativos por fase.
- Fase 4: sección 5 rediseñada como resultados fase a fase de NeuroMIND.
- Fase 5: documentación externa consolidada en `README.md`, `report/README.md`, `CHANGELOG.md` y `CODEX.md`.

Los apartados siguientes se conservan como registro de la auditoría original y como referencia para futuras ampliaciones con más sujetos, CSD activo y análisis grupales.

## Hallazgos principales

- NeuroMIND ya genera suficientes outputs para ampliar la sección 5 del informe desde una lectura resumida hacia una lectura fase a fase.
- El dashboard cubre más granularidad que el informe metodológico: el informe agrupa el pipeline en Fase 0--8, mientras que el dashboard usa paneles Fase 0--14.
- La sección 5 actual ya integra QC, PSD, potencia por banda, wPLI alfa/theta y surrogates/FDR, pero todavía no explota ICA, segmentación, rechazo de artefactos, métricas de red, matrices p/q completas ni paneles del dashboard.
- La ampliación debe separar con claridad tres niveles: resultado observado, interpretación técnica y limitación metodológica.

## Equivalencia de fases

| Bloque conceptual | Informe actual | Dashboard | Módulos Julia principales | Estado para informe |
|---|---:|---:|---|---|
| Proyecto, dataset y trazabilidad | Resumen / Datos | 0 | `src/io/Config.jl`, `src/types.jl` | Añadir tabla de ejecución |
| BIDS y metadatos | Fase 0 | 1 | `src/io/BIDSLoader.jl` | Completar con validación |
| Señal raw | Fase 0--1 | 2 | `src/io/BIDSLoader.jl`, `src/qc/QualityControl.jl` | Incluir vista + métricas |
| QC inicial | Fase 0 | 3 | `src/qc/QualityControl.jl` | Ampliar sección 5.2 |
| Filtrado | Fase 1 | 4 | `src/preprocessing/Filtering.jl` | Añadir respuesta/impacto |
| ICA | Fase 2a--2b | 5 | `src/ica/ICACore.jl`, `src/ica/ICAClassification.jl`, `src/ica/ICAInspection.jl` | Nueva subsección |
| Segmentación | Fase 3 | 6 | `src/segmentation/Epochs.jl` | Nueva tabla de épocas |
| Rechazo de artefactos | Fase 5 | 7 | `src/segmentation/Epochs.jl` | Nueva subsección o unir con segmentación |
| Espectral | Fase 6 | 8 | `src/spectral/PowerSpectrum.jl` | Ya existe, ampliar con tablas |
| wPLI | Fase 7 | 9 | `src/connectivity/wPLI.jl`, `src/connectivity/CSD.jl`, `src/connectivity/GraphMetrics.jl` | Ya existe, añadir métricas de red |
| Surrogates/FDR | Fase 8 | 10 | `src/statistics/Surrogates.jl`, `src/statistics/FDR.jl` | Ya existe, ampliar trazabilidad |
| Resumen final | Seguimiento | 11 | `src/report/HTMLReport.jl`, `src/webapp/App.jl` | Usar como síntesis |
| Exportación/reporte | Apéndices / CODEX | 12 | `src/report/HTMLReport.jl`, `src/webapp/App.jl` | Documentar en apéndice dashboard |
| Transversal | Diseño previsto | 13 | `src/statistics/GroupStats.jl`, `src/visualization/*` | Dejar como preparado si no hay cohorte |
| Longitudinal | Diseño previsto | 14 | `src/longitudinal/LongitudinalAnalysis.jl` | Dejar como preparado si no hay cohorte |

## Endpoints dashboard relevantes

| Panel | Endpoint principal | Función | Uso propuesto en informe |
|---:|---|---|---|
| 0 | `/api/project_overview` | Resumen del proyecto, dataset y estado | Apéndice dashboard |
| 1 | `/api/bids_metadata` | Metadatos BIDS y validación | Apéndice + metodología |
| 2 | `/api/phase2_info` | Información de señal raw | Resultados fase raw/QC |
| 3 | `/api/bids/channel_stats` | Estadísticas por canal | Resultados QC |
| 4 | `/api/filter_config`, `/api/filter_response` | Parámetros y respuesta del filtrado | Metodología + resultados de filtrado |
| 5 | `/api/phase5_ica_info`, `/api/ica_features` | ICA, componentes y scores | Nueva subsección ICA |
| 6 | `/api/phase6_segmentation` | Segmentación y cobertura | Nueva subsección segmentación |
| 7 | `/api/phase7_ar` | Rechazo de artefactos por época | Nueva subsección rechazo |
| 8 | `/api/phase8_spectral` | PSD, bandas, regiones e índices | Ampliar espectral |
| 9 | `/api/phase9_wpli` | Matrices wPLI, red y aristas | Ampliar conectividad |
| 10 | `/api/phase10_surrogates` | Observado, nula, p/q, significativas | Ampliar inferencia |
| 11 | `/api/phase11_summary` | KPIs integrados | Síntesis técnica |
| 12 | `/api/phase12_files` | Inventario de outputs | Apéndice dashboard/exportación |
| 13 | `/api/phase13_transversal` | Comparación transversal | Pendiente hasta cohorte |
| 14 | `/api/phase14_longitudinal` | Seguimiento T1/T2 | Pendiente hasta cohorte |

## Inventario de outputs auditados

### Trazabilidad y configuración

| Archivo | Contenido | Decisión |
|---|---|---|
| `overview.csv` | Sujeto, sesión, tarea, frecuencia, canales, duración, canales malos | Incluir como tabla compacta |
| `config_snapshot.toml` | Configuración completa de ejecución | Resumir; conservar como trazabilidad |
| `pipeline_log.txt` | Registro textual de ejecución | Citar en apéndice/exportación |

### QC inicial y señal

| Archivo | Contenido | Decisión |
|---|---|---|
| `qc_summary.csv` | Media, RMS, desviación, rango, z-RMS, flags por canal | Incluir top canales y resumen |
| `channel_statistics.csv` | Duplicado funcional de QC por canal | Usar como fuente o archivar como redundante |
| `channel_coverage.csv` | Cobertura por canal tras segmentación/rechazo | Incluir gráfico/tabla breve |
| `channel_artifact_summary.csv` | Canales implicados en artefactos | Resumir si aporta detalle |
| `filtered_signal_preview.png` | Vista raw/filtrada | Incluir o sustituir por figura actual recortada |

### ICA

| Archivo | Contenido | Decisión |
|---|---|---|
| `ica_summary.json` | Componentes, aceptados/rechazados, duración, tipos | Incluir tabla resumen |
| `ica_component_features.csv` | Scores ocular, muscular, línea, jump por IC | Incluir top componentes por score |
| `ica_components.csv` | Matriz/componentes espaciales | No incluir completa; fuente auxiliar |
| `ica_activations.csv` | Activaciones temporales | No incluir completa; posible figura de ejemplo |
| `ica_mixing_matrix.csv`, `ica_unmixing_matrix.csv` | Matrices internas ICA | No incluir; conservar trazabilidad |
| `ica_signal_before.csv`, `ica_signal_after.csv` | Señal antes/después | No incluir completa; usar si se crea figura comparativa |

### Segmentación y rechazo de artefactos

| Archivo | Contenido | Decisión |
|---|---|---|
| `segmentation_summary.json` | Épocas totales/válidas/rechazadas, longitud, solape | Incluir tabla |
| `segments_table.csv` | Estado de 100 épocas | Resumir con conteos; no incluir completa |
| `artifact_rejection_summary.json` | Rechazo por amplitud/gradiente | Incluir tabla breve |
| `rejected_segments.csv` | Épocas rechazadas | Incluir si son pocas; útil para trazabilidad |

### Espectral

| Archivo | Contenido | Decisión |
|---|---|---|
| `spectral_summary.json` | Método, ventana, nfft, número de épocas | Incluir en tabla de parámetros |
| `psd_by_channel.csv` | PSD por canal y frecuencia | No incluir completa; usar para gráficos |
| `psd_all_channels.png` | PSD multicanal | Ya incluida; mantener |
| `band_power_summary.csv` | Potencia por banda/canal | Incluir resumen por banda/región |
| `band_power_summary.png` | Mapa/figura de potencia por banda | Ya incluida; mantener |
| `regional_psd.csv` | PSD por región | Buena candidata para tabla/figura nueva |
| `spectral_indices.csv` | Índices espectrales por canal | Candidato secundario |

### Conectividad wPLI y red

| Archivo | Contenido | Decisión |
|---|---|---|
| `connectivity_summary.json` | Medias, máximos, densidad por banda | Incluir tabla principal |
| `connectivity_edges.csv` | Aristas por banda | Incluir top aristas por banda |
| `network_metrics.csv` | Strength, degree, norm_strength por canal | Añadir gráfico/tabla de hubs exploratorios |
| `wpli_*.csv` | Matrices wPLI por banda | No incluir completas; usar para figuras |
| `wpli_*.png` | Heatmaps por banda | Mantener alfa/theta; delta/beta/gamma como material auxiliar |

### Surrogates/FDR

| Archivo | Contenido | Decisión |
|---|---|---|
| `surrogate_summary.json` | Método, n, alpha, FDR, totales significativos | Ya incluido; reforzar trazabilidad |
| `surrogate_quality.csv` | Resumen por banda de observado/nula/p/FDR | Ya incluido; mantener |
| `surrogate_null_stats_*.csv` | Nula por arista y banda | No incluir completa; usar para figura de distribución si procede |
| `wpli_observed_*.csv` | Matrices observadas para contraste | No incluir completa |
| `wpli_pvalues_*.csv` | p-valores por banda | Resumir; posible figura p-value distribution |
| `wpli_qvalues_*.csv` | q-valores por banda | Resumir |
| `wpli_significant_*.csv` | Máscaras significativas | Resumir; todas vacías tras FDR |
| `significant_connections.csv` | Aristas significativas finales | Incluir como resultado negativo: 0 aristas |

## Propuesta de nueva sección 5

Estructura recomendada para la Fase 4 de trabajo:

```text
5 Resultados fase a fase de NeuroMIND
  5.1 Ejecución analizada y trazabilidad
  5.2 Fase 0-1: BIDS, señal raw y QC inicial
  5.3 Fase 2: filtrado y respuesta espectral
  5.4 Fase 3: ICA y clasificación de componentes
  5.5 Fase 4-5: segmentación, baseline y rechazo de artefactos
  5.6 Fase 6: análisis espectral y potencia por banda
  5.7 Fase 7: conectividad wPLI y métricas de red
  5.8 Fase 8: surrogates, FDR e inferencia
  5.9 Lectura integrada de la ejecución
  5.10 Limitaciones y próximos relanzamientos
```

## Propuesta de apéndice dashboard

Estructura recomendada para la Fase 2 de trabajo:

```text
Apéndice técnico: dashboard NeuroMIND
  Propósito del dashboard
  Arquitectura Genie.jl + HTML/JS
  Relación results/ → API → paneles
  Tabla de paneles y endpoints
  Paneles implementados y paneles preparatorios
  Uso práctico para QC y depuración
  Limitaciones
```

## Criterios de selección para el informe

Incluir directamente:

- Tablas compactas con menos de 15 filas o agregados por banda/fase.
- Figuras que resumen decisiones QC o resultados interpretables.
- Resultados negativos relevantes, como ausencia de aristas significativas tras FDR.
- Captions explicativas con interpretación técnica prudente.

Resumir:

- CSV largos por canal, frecuencia, arista o época.
- Matrices completas repetidas por banda.
- Logs y snapshots de configuración.

Conservar fuera del flujo narrativo:

- Matrices ICA completas.
- Activaciones temporales completas.
- Matrices p/q por banda completas.
- Figuras redundantes de bandas no discutidas en texto.

## Riesgos y cautelas

- La ejecución auditada es individual; no permite conclusiones clínicas ni comparaciones grupales.
- Las matrices actuales se interpretan en espacio sensor con `use_csd=false`.
- `significant_connections.csv` está vacío tras FDR; esto debe presentarse como resultado inferencial negativo, no como fallo del pipeline.
- La numeración del dashboard no coincide exactamente con la numeración del informe; debe explicarse en el apéndice para evitar confusión.
- Las figuras del dashboard son herramientas de inspección interactiva; el informe debe usar versiones estáticas seleccionadas o reconstrucciones vectoriales.

## Siguiente fase

La Fase 2 debe incorporar al informe un apéndice técnico del dashboard usando esta auditoría como fuente. No requiere modificar algoritmos Julia. Los cambios esperados se concentrarán en:

- `report/sections_es/appendix.tex`
- `report/assets/figures/generated/neurosmart_generated_figures.tex` si se crea un diagrama TikZ nuevo
- `report/CHANGELOG.md`
- `CODEX.md`
