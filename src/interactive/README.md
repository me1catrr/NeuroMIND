# README — `src/interactive`

Guía operativa de lanzamiento de visores interactivos, en orden recomendado: primero
single-subject (`sub-M05/ses-T2/eyesclosed`), después cohorte transversal y finalmente
cohorte longitudinal.

## Contexto rápido

- Ruta base single-subject: `results/subjects/sub-M05/ses-T2/eyesclosed/`
- Soporte compartido (no se lanza): `viewer_support.jl` y `viewer_common.js`

## Paso a paso recomendado

### Paso n.º 1 — `plot_raw.jl` (`:8765`)

Muestra la señal EEG cruda por canal (selector de canales, ventana temporal y escala).

````
julia --project=. src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_raw.jl
````

- **Archivos que lee:** `tables/raw_signal.csv`

### Paso n.º 2 — `plot_raw_butterfly.jl` (`:8766`)

Muestra la señal cruda en modo butterfly/montaje apilado para comparar canales a la vez.

````
julia --project=. src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_raw_butterfly.jl
````

- **Archivos que lee:** `tables/raw_signal.csv`

### Paso n.º 3 — `plot_raw_histogram.jl` (`:8768`)

Muestra histogramas/densidades de amplitud (por canal o combinados).

````
julia --project=. src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_raw_histogram.jl
````

- **Archivos que lee:** `tables/raw_signal.csv`

### Paso n.º 4 — `plot_filtered_vs_raw.jl` (`:8769`)

Muestra la comparación de la señal raw frente a cada etapa de filtrado (notch, bandreject, highpass, lowpass).

````
julia --project=. src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_filtered_vs_raw.jl
````

- **Archivos que lee:**
  - `tables/raw_signal.csv`
  - `tables/filtered_signal_notch.csv`
  - `tables/filtered_signal_bandreject.csv`
  - `tables/filtered_signal_highpass.csv`
  - `tables/filtered_signal_lowpass.csv`

### Paso n.º 5 — `plot_ica_components.jl` (`:8770`)

Muestra un panel de inspección ICA (activaciones de IC, topomap, PSD de IC, features y decisión keep/reject/review).

````
julia --project=. src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_ica_components.jl
````

- **Archivos que lee:**
  - `tables/ica/ica_components.csv`
  - `tables/ica/ica_component_features.csv`
  - `tables/ica/ica_activations.csv`
  - `json/ica_summary.json`
  - `figures/ica/ica_topomap_*.png`

### Paso n.º 6 — `plot_ica_before_after.jl` (`:8771`)

Muestra la comparación temporal de señal antes y después de ICA por canal.

````
julia --project=. src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_ica_before_after.jl
````

- **Archivos que lee:**
  - `tables/ica/ica_signal_before.csv`
  - `tables/ica/ica_signal_after.csv`
  - `json/ica_summary.json`

### Paso n.º 7 — `plot_epochs.jl` (`:8772`)

Muestra calidad de épocas (histograma), superposición de épocas válidas y visualización de épocas rechazadas.

````
julia --project=. src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_epochs.jl
````

- **Archivos que lee:**
  - `tables/segments_table.csv`
  - `tables/rejected_segments.csv`
  - `cache/ica_result.jls`
  - `ica_labels_auto.csv`
  - `config_snapshot.toml`

### Paso n.º 8 — `plot_baseline.jl` (`:8773`)

Muestra el efecto de la corrección de baseline (before/after) por canal y época, con resumen de offsets.

````
julia --project=. src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_baseline.jl
````

- **Archivos que lee:**
  - `cache/ica_result.jls`
  - `ica_labels_auto.csv`
  - `tables/segments_table.csv`
  - `config_snapshot.toml`

### Paso n.º 9 — `plot_spectral.jl` (`:8767`)

Muestra el análisis espectral completo en cuatro modos (`bands_abs`, `psd`, `bands_rel`, `asd`) e índices espectrales.

````
julia --project=. src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_spectral.jl
````

- **Archivos que lee:**
  - `tables/psd_by_channel.csv`
  - `tables/band_power_summary.csv`
  - `tables/spectral_indices.csv`
  - `tables/regional_psd.csv`

### Paso n.º 10 — `plot_spectral_PSD.jl` (`:8767`)

Muestra una versión simplificada centrada solo en curvas PSD (canal seleccionado vs promedio).

````
julia --project=. src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_spectral_PSD.jl
````

- **Archivos que lee:**
  - `tables/psd_by_channel.csv`
  - `tables/band_power_summary.csv`
  - `tables/spectral_indices.csv`
  - `tables/regional_psd.csv`
- **Nota:** comparte puerto con `plot_spectral.jl`; úsalo como vista rápida, dejando `plot_spectral.jl` como visor principal.

### Paso n.º 11 — `plot_connectivity.jl` (`:8774`)

Muestra conectividad wPLI (heatmaps, grafo 10-20, hubs y métricas de red por banda).

````
julia --project=. src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_connectivity.jl
````

- **Archivos que lee:**
  - `tables/connectivity/wpli_{BAND}.csv`
  - `tables/connectivity/connectivity_edges.csv`
  - `tables/connectivity/network_metrics.csv`
  - `json/connectivity_summary.json`

### Paso n.º 12 — `plot_surrogate.jl` (`:8775`)

Muestra inferencia estadística por surrogates (mapas p/q/z/sig, volcano, comparación observado vs nulo, red FDR).

````
julia --project=. src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_surrogate.jl
````

- **Archivos que lee:**
  - `tables/surrogate/*`
  - `json/surrogate_summary.json`
  - `cache/surrogate/null_distribution_{band}.jls` (si existe)

### Paso n.º 13 — `plot_transversal.jl` (`:8781`)

Muestra la comparación de cohorte MS vs Control (overview, conectividad, volcano, espectral y métricas globales).

````
julia --project=. src/interactive/plot_transversal.jl
````

- **Archivos que lee:**
  - `results/transversal/eyesclosed/*`
  - `results/transversal/eyesopen/*`
  - `statistics_contract.toml`
  - `transversal_summary.json`

### Paso n.º 14 — `plot_longitudinal.jl` (`:8780`)

Muestra cambios pareados T1->T2 por condición (EC/EO), con conectividad, potencia y métricas globales.

````
julia --project=. src/interactive/plot_longitudinal.jl
````

- **Archivos que lee:**
  - `results/longitudinal/eyesclosed/*`
  - `results/longitudinal/eyesopen/*`
  - `statistics_contract.toml`
  - `longitudinal_summary.json`

## Notas prácticas

- Si un puerto está ocupado, cierra el visor anterior antes de lanzar otro.
- El flujo recomendado single-subject es: raw -> filtrado -> ICA -> epochs/baseline -> espectral -> conectividad -> surrogate.
