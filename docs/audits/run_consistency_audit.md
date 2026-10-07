# Auditoría de trazabilidad — corrida sub-M05/ses-T2/eyesclosed

**Fecha de auditoría:** 2026-07-14  
**Modo:** solo lectura (sin modificar código, documentos, resultados ni presentación)  
**Auditor:** agente técnico (Cursor)  
**Alcance:** corrida canónica usada en el manual LaTeX (`Report_Pre/`) y la presentación Slidev (`Report_Pre/slidev/`)

---

## Veredictos por bloque

| Bloque | Veredicto | Motivo principal |
|--------|-----------|------------------|
| Identidad de corrida | **ÁMBAR** | Corrida principal identificada (2026-07-09), pero mezcla fases de mayo (caché ICA, figuras `_EC.png`) |
| QC | **ÁMBAR** | Métricas CSV/JSON coherentes con julio; figuras del manual en `Report_Pre/figures/plots/` desactualizadas (2026-07-03) |
| Filtrado | **VERDE** | `config_snapshot.toml`, `pipeline_log.txt` y señales alineados |
| ICA | **ÁMBAR** | Clasificación, figuras y exports de julio; **matriz FastICA de mayo** (caché) |
| Segmentación | **VERDE** | `segmentation_summary.json`, `rejected_segments.csv` y figuras coherentes |
| Espectral | **VERDE** | `spectral_summary.json`, `band_power_summary.csv` regenerados 2026-07-09 |
| Conectividad | **VERDE** | `connectivity_summary.json`, `wpli_*.csv/png` (sin sufijo) de julio |
| Surrogates | **VERDE** | Única corrida cohorte con `surrogate_summary.json` reciente; verificación puntual documentada |
| Manual (`Report_Pre/`) | **ÁMBAR** | Cifras mayormente correctas; **figuras QC copiadas sin resincronizar** tras corrida julio |
| Presentación (Slidev) | **ÁMBAR** | Manifiesto y `subject-data.json` bien trazados; **`signal_preview_EC.png` es de mayo** |

**Veredicto global: ÁMBAR** — reproducible con advertencias explícitas sobre caché ICA, figuras `_EC`/`figures/plots` antiguas y surrogates como verificación puntual.

---

## 1. Resumen ejecutivo

La corrida de referencia para el informe técnico es **`sub-M05 / ses-T2 / eyesclosed`**, ejecutada el **2026-07-09** entre las **10:49:50** y **~10:58:32** (521.7 s), según `pipeline_log.txt`. El identificador de manifiesto es **`sub-M05_ses-T2_eyesclosed`**.

La mayoría de tablas y JSON de las fases 5–8 (segmentación, espectral, conectividad, surrogates) pertenecen de forma demostrable a esa ejecución. Las cifras numéricas del manual y de Slidev coinciden con los archivos de resultados en los casos verificados (Fp2 z=4.63, 99 épocas válidas, 447 conexiones significativas, etc.).

Los problemas de trazabilidad más graves son:

1. **ICA:** la descomposición FastICA proviene de una caché del **2026-05-24** (`results/M05/T2/cache/ica_result.jls`), 46 días antes del resto de la corrida.
2. **Figuras desactualizadas en el manual:** `Report_Pre/figures/plots/qc_raw_butterfly.png` y `qc_psd_raw_average.png` datan del **2026-07-03**, mientras las copias en `results/.../figures/` son del **2026-07-09**.
3. **Figura desactualizada en Slidev:** `figures/signal_preview_EC.png` (y su copia en `slidev/public/`) es del **2026-05-24**; el pipeline de julio escribió `filtered_signal_preview.png` en la raíz pero no actualizó la ruta del manifiesto.
4. **10 archivos `figures/*_EC.png`** en la carpeta BIDS siguen siendo de mayo y no fueron regenerados en julio.
5. **`channel_statistics_compare.csv`** existe (generado 2026-07-09 11:07) pero **no tiene generador en `src/` ni `scripts/`** del repositorio.

El manual y Slidev documentan varias de estas limitaciones; el PDF del manual compilado el **2026-07-14** puede seguir mostrando figuras QC anteriores a la corrida canónica.

---

## 2. Alcance

### Incluido

- Corrida canónica `results/subjects/sub-M05/ses-T2/eyesclosed/`
- Manual LaTeX: `Report_Pre/main.tex` → `Report_Pre/main.pdf`
- Presentación Slidev: `Report_Pre/slidev/slides.md` → `slides-export.pdf`
- Manifiesto YAML, inventario SHA-256, `subject-data.json`
- Caché ICA legacy: `results/M05/T2/cache/`
- Documentación operativa: `AGENTS.md`, `README.md` (anexo M05), auditorías previas
- Logs batch: `logs/batch_run_2026-05-26_11-33.csv`

### Excluido

- Contenido teórico del manual sin afirmaciones verificables sobre M05
- Cohorte completa (199 runs del 2026-05-26) salvo como contexto
- Validación BrainVision Analyzer (sin CSV en repo)
- Modificación de ningún artefacto

---

## 3. Fuentes revisadas

| Categoría | Ruta | Rol |
|-----------|------|-----|
| Instrucciones agentes | `AGENTS.md`, `CLAUDE.md` | Contexto operativo |
| README anexo M05 | `README.md` §Anexo | Trazabilidad documentada |
| Inventario salidas | `Report_Pre/INVENTARIO_SALIDAS_M05_2026-07-02.md` | Mapa figuras/tablas |
| Mapa resultados Slidev | `Report_Pre/slidev/docs/m05-results-map.md` | Narrativa verificada |
| Auditoría legacy | `docs/REPORT_PHASE1_AUDIT.md` | Referencia histórica (`report/` legacy) |
| Log ejecución | `results/subjects/sub-M05/ses-T2/eyesclosed/pipeline_log.txt` | Cronología paso a paso |
| Config efectiva | `results/subjects/sub-M05/ses-T2/eyesclosed/config_snapshot.toml` | Única fuente reproducible |
| Config original | `config/_scratch_surrogates_verification.toml` | **Ausente** (citada en log) |
| Manifiesto | `Report_Pre/slidev/manifests/runs/sub-M05_ses-T2_eyesclosed.yml` | Congelación corrida |
| Inventario hashes | `Report_Pre/slidev/generated/run-manifest/sub-M05_ses-T2_eyesclosed.json` | SHA-256 por archivo |
| Datos agregados slides | `Report_Pre/slidev/generated/subjects/M05/subject-data.json` | Cifras en diapositivas |
| Manual fuente | `Report_Pre/chapters/*.tex` | Afirmaciones textuales |
| Presentación | `Report_Pre/slidev/slides.md` | Afirmaciones y rutas de figuras |
| PDF manual | `Report_Pre/main.pdf` (mtime 2026-07-14 13:10) | Compilación |
| PDF slides | `Report_Pre/slidev/slides-export.pdf` (mtime 2026-07-14 12:59) | Exportación |
| Caché ICA | `results/M05/T2/cache/ica_result.jls` | FastICA serializado |
| Índice cohorte | `results/subjects_index.csv` | `processed_at` por grabación |
| Log batch | `logs/batch_run_2026-05-26_11-33.csv` | Contexto M05 SKIPPED |

---

## 4. Identidad de la corrida principal

| Campo | Valor | Evidencia |
|-------|-------|-----------|
| `manifest_id` | `sub-M05_ses-T2_eyesclosed` | Manifiesto YAML |
| Sujeto | M05 | `overview.csv`, log |
| Sesión | T2 | `overview.csv`, log |
| Condición | eyesclosed (EC) | `overview.csv`, log |
| Inicio | `2026-07-09T10:49:50.695` | `pipeline_log.txt` L1 |
| Fin | ~`2026-07-09T10:58:32` (521.7 s) | `pipeline_log.txt` L92 |
| Comando | `scripts/run_single_subject.jl` | Manifiesto, log |
| Config usada | `_scratch_surrogates_verification.toml` | Log L3 — **archivo no existe** |
| Config reproducible | `config_snapshot.toml` | Snapshot en export_dir |
| `git_commit` | `null` | Manifiesto (no registrado en corrida) |
| `config_hash` | `null` | Manifiesto (no registrado en corrida) |
| Canales | 31 (Fp2 incluido) | `overview.csv`, `[montage] exclude_fp2=false` |
| Muestras / fs / duración | 50180 / 500 Hz / 100.36 s | `overview.csv`, log |
| Carga datos | TSV legacy (no solo .vhdr) | `README.md` §Fase 0 |

**Nota:** No existe campo `run_id` en el pipeline Julia. La trazabilidad depende del manifiesto Slidev + `config_snapshot.toml` + `pipeline_log.txt`.

---

## 5. Cronología reconstruida

Zona horaria de los logs del pipeline: **hora local del sistema** (España, UTC+2 en julio). El manifiesto JSON usa **UTC** (`generated_at: 2026-07-14T10:39:19.452Z`).

| Fecha/hora | Tipo de evento | Evidencia | Ruta | Confianza | Observaciones |
|------------|----------------|-----------|------|-----------|---------------|
| 2026-05-24 17:34:57 | Ejecución ICA temprana (cohorte) | `timestamp` mínimo en `ica_summary.json` cohorte | Varios sujetos | Alta | Contexto batch inicial |
| 2026-05-24 18:47:00 | Serialización caché ICA M05 | mtime FS; manifiesto `cache_date` | `results/M05/T2/cache/ica_result.jls` | Alta | **Origen real de FastICA** |
| 2026-05-24 18:46 | Generación figuras `*_EC.png` | mtime FS | `results/.../eyesclosed/figures/*_EC.png` (10 archivos) | Alta | No sobrescritas en julio |
| 2026-05-26 11:33–12:19 | Batch masivo Phase C | `logs/batch_run_2026-05-26_11-33.csv` | 199 runs `processed_at` 2026-05-26 | Alta | M05 EC **SKIPPED** |
| 2026-05-26 11:46:18 | Run rápido M05 eyesopen | log / batch CSV | `eyesopen/pipeline_log.txt` (3.2 s) | Alta | Sin `ica_summary.json` |
| 2026-07-03 11:49–18:38 | Copia manual figuras QC | mtime FS | `Report_Pre/figures/plots/qc_*` | Alta | Copia manual, no pipeline |
| 2026-07-06 09:44 | Walkthrough raw regenerado | mtime FS | `Report_Pre/figures/plots/walkthrough_M05_raw_*.png` | Media | Anterior a corrida julio |
| 2026-07-08 | Decisión montaje Fp2 | comentario snapshot | `config_snapshot.toml [montage]` | Alta | `exclude_fp2=false` |
| **2026-07-09 10:49:50** | **Inicio ejecución principal** | log L1 | `pipeline_log.txt` | **Alta** | Config scratch (borrada) |
| 2026-07-09 10:49:52 | Carga BIDS completada | log L6–8 | log | Alta | 31 ch, 50180 muestras |
| 2026-07-09 10:49:53 | QC canales | log L11–12 | log, `channel_statistics.csv` | Alta | Fp2 z>3, σ̄=14.4 µV |
| 2026-07-09 10:50:11 | Figuras QC generadas | log L13 | `figures/qc_*.png` (jul-09) | Alta | 18 s de generación |
| 2026-07-09 10:50:11 | **Restauración caché ICA** | log L23 | cache + log | Alta | `duration_s: 0.9` en JSON |
| 2026-07-09 10:50:12 | Rechazo auto IC 5,9,13,14,21 | log L25 | log | Alta | Umbral 1.5 documentado en log |
| 2026-07-09 10:50:19 | Export ICA completado | `ica_summary.json` timestamp | `ica_summary.json` | Alta | Timestamp = export, no FastICA |
| 2026-07-09 10:50:21 | Segmentación + AR | `segmentation_summary.json` | JSON | Alta | 100→99 épocas |
| 2026-07-09 10:50:24–10:58:28 | Espectral + wPLI | summaries JSON | `spectral/connectivity_summary.json` | Alta | |
| 2026-07-09 10:50:27–10:58:26 | Surrogates (8 bandas) | log L65–74 | `surrogate_summary.json` | Alta | 447 sig. total |
| 2026-07-09 10:58:30 | Índice cohorte actualizado | `processed_at` | `subjects_index.csv` | Alta | |
| 2026-07-09 11:07 | Post-proceso walkthrough | mtime FS | `channel_statistics_compare.csv`, `results/walkthrough_*.png` | Media | Fuera del pipeline estándar |
| 2026-07-13 14:24 | Sync assets Slidev | mtime FS | `slidev/public/subjects/M05/figures/` | Alta | Copió `signal_preview_EC.png` de mayo |
| 2026-07-14 10:39:19Z | Inventario SHA-256 | `generated_at` | `run-manifest/*.json` | Alta | Metadatos inventario |
| 2026-07-14 12:59 | Export PDF Slidev | mtime FS | `slides-export.pdf` | Media | Compilación, no ejecución |
| 2026-07-14 13:10 | Compilación PDF manual | mtime FS | `main.pdf` | Media | Compilación, no ejecución |

### Distinción de tipos de fecha

| Tipo | Ejemplo | Interpretación |
|------|---------|----------------|
| Ejecución principal | 2026-07-09 10:49:50 | Inicio pipeline completo |
| Ejecución parcial | 2026-05-24 18:47 | Solo FastICA → caché |
| Restauración caché | 2026-07-09 10:50:11 | Carga `ica_result.jls`, no recalcula ICA |
| Regeneración figuras | 2026-07-09 10:50:08–17 | QC/ICA figuras en export_dir |
| Copia manual/post-hoc | 2026-07-03, 2026-07-13 | `figures/plots/`, sync Slidev |
| Compilación documento | 2026-07-14 | PDF manual/slides |
| Fecha declarada manualmente | 2026-07-09 en capítulos | Debe contrastarse con log |

---

## 6. Inventario de outputs por fase

Leyenda de **Estado:**

- **confirmado** — pertenece a la corrida 2026-07-09
- **caché válida** — reutilizado por hash de config ICA
- **probable** — inferido, no demostrado al 100 %
- **ejecución anterior** — mtime o contenido de mayo/antes
- **copiado posterior** — generado fuera del pipeline estándar
- **inconsistente** — conflicto entre fuentes
- **no localizable** — procedencia de código no encontrada

Hashes SHA-256 del inventario `run-manifest` (generado 2026-07-14).

### Fase 1 — Carga BIDS

| Fase | Output | Ruta | Fecha interna | Fecha FS | Config | Hash (prefijo) | Estado |
|------|--------|------|---------------|----------|--------|----------------|--------|
| 1 | overview | `overview.csv` | — | 2026-07-09 10:58 | snapshot | `af75ff80…` | confirmado |
| 1 | raw signal | `raw_signal.csv` | — | 2026-07-09 10:50 | snapshot | (en manifiesto) | confirmado |
| 1 | log | `pipeline_log.txt` | 2026-07-09T10:49:50 | 2026-07-09 10:58 | scratch† | `c30ce7ab…` | confirmado |

### Fase 2 — QC inicial

| Fase | Output | Ruta | Fecha interna | Fecha FS | Config | Hash (prefijo) | Estado |
|------|--------|------|---------------|----------|--------|----------------|--------|
| 2 | channel stats | `channel_statistics.csv` | — | 2026-07-09 10:58 | `[qc]` | `ea1968b3…` | confirmado |
| 2 | qc summary | `qc_summary.csv` | — | 2026-07-09 10:58 | `[qc]` | (manifiesto) | confirmado |
| 2 | butterfly QC | `figures/qc_raw_butterfly.png` | — | 2026-07-09 10:50 | `[qc.figures]` | `846d55a5…` | confirmado |
| 2 | PSD raw avg | `figures/qc_psd_raw_average.png` | — | 2026-07-09 10:50 | `[qc.figures]` | (manifiesto) | confirmado |
| 2 | variance topomap | `figures/qc_variance_topomap.png` | — | 2026-07-09 10:50 | `[qc.figures]` | (manifiesto) | confirmado |
| 2 | compare raw/filt | `channel_statistics_compare.csv` | — | 2026-07-09 11:07 | — | — | copiado posterior / no localizable |
| 2 | copia manual QC | `Report_Pre/figures/plots/qc_*.png` | — | **2026-07-03** | — | — | **ejecución anterior** (respecto a jul-09) |

### Fase 3 — Filtrado

| Fase | Output | Ruta | Fecha interna | Fecha FS | Config | Hash (prefijo) | Estado |
|------|--------|------|---------------|----------|--------|----------------|--------|
| 3 | señal pre-ICA | `ica_signal_before.csv` | — | 2026-07-09 10:50 | `[filtering] eeg_julia` | (manifiesto) | confirmado |
| 3 | preview filtrado | `filtered_signal_preview.png` | — | 2026-07-09 10:58 | `[filtering]` | — | confirmado (no en manifiesto) |
| 3 | preview EC stale | `figures/signal_preview_EC.png` | — | **2026-05-24** | — | `a6367127…` | **ejecución anterior** |
| 3 | walkthrough plots | `Report_Pre/figures/plots/walkthrough_M05_*` | — | 2026-07-06 / 07-09 | — | — | probable (mezcla fechas) |

### Fase 4 — ICA

| Fase | Output | Ruta | Fecha interna | Fecha FS | Config | Hash (prefijo) | Estado |
|------|--------|------|---------------|----------|--------|----------------|--------|
| 4 | caché FastICA | `results/M05/T2/cache/ica_result.jls` | — | **2026-05-24** | `[ica] eeg_julia` | — | **caché válida** |
| 4 | hash config ICA | `results/M05/T2/cache/ica_config.hash` | — | **2026-05-24** | hash `-26e3c03259100db8` | — | caché válida |
| 4 | summary | `ica_summary.json` | 2026-07-09T10:50:19 | 2026-07-09 10:50 | `[ica]` | `b7cb542b…` | confirmado (export julio, ICA mayo) |
| 4 | components | `ica_components.csv` | — | 2026-07-09 10:50 | umbral 1.5 | (manifiesto) | confirmado |
| 4 | señal post-ICA | `ica_signal_after.csv` | — | 2026-07-09 10:50 | — | (manifiesto) | confirmado |
| 4 | topomaps / butterfly | `figures/ica_*.png` | — | 2026-07-09 10:50 | `[ica.figures]` | (manifiesto) | confirmado |

**Componentes rechazados (verificado en `ica_components.csv`):** IC5, IC9, IC13, IC14, IC21 — tipos jump×2, line_noise×3. Varianza retenida: 82.6%.

### Fase 5 — Segmentación, baseline, rechazo

| Fase | Output | Ruta | Fecha interna | Fecha FS | Config | Hash (prefijo) | Estado |
|------|--------|------|---------------|----------|--------|----------------|--------|
| 5 | segmentation | `segmentation_summary.json` | 2026-07-09T10:50:21 | 2026-07-09 10:50 | epoch 1.0 s, 0 overlap | `e1d788a5…` | confirmado |
| 5 | segments table | `segments_table.csv` | — | 2026-07-09 10:50 | — | (manifiesto) | confirmado |
| 5 | rejected epoch | `rejected_segments.csv` | — | 2026-07-09 10:50 | época 12, C4 72.56 µV | (manifiesto) | confirmado |
| 5 | AR summary | `artifact_rejection_summary.json` | 2026-07-09T10:50:21 | 2026-07-09 10:50 | ±70 µV, 31 ch | (manifiesto) | confirmado |
| 5 | baseline fig | `figures/baseline_before_after_C4_ep089.png` | — | 2026-07-09 10:50 | 0–0.10 s, 2 passes | (manifiesto) | confirmado |

### Fase 6 — Espectral

| Fase | Output | Ruta | Fecha interna | Fecha FS | Config | Hash (prefijo) | Estado |
|------|--------|------|---------------|----------|--------|----------------|--------|
| 6 | spectral summary | `spectral_summary.json` | 2026-07-09 10:58:28 | 2026-07-09 10:58 | fft_hamming, nfft=1024 | `40580758…` | confirmado |
| 6 | band power | `band_power_summary.csv` | — | 2026-07-09 10:58 | 7 bandas | (manifiesto) | confirmado |
| 6 | PSD fig (raíz) | `psd_all_channels.png` | — | 2026-07-09 10:58 | — | (manifiesto) | confirmado |
| 6 | PSD fig EC stale | `figures/psd_all_channels_EC.png` | — | **2026-05-24** | — | — | ejecución anterior |
| 6 | regional PSD | `regional_psd.csv` | — | 2026-07-09 10:58 | 5 regiones × 7 bandas | (manifiesto) | confirmado |

### Fase 7 — Conectividad wPLI

| Fase | Output | Ruta | Fecha interna | Fecha FS | Config | Hash (prefijo) | Estado |
|------|--------|------|---------------|----------|--------|----------------|--------|
| 7 | connectivity | `connectivity_summary.json` | 2026-07-09T10:58:29 | 2026-07-09 10:58 | hilbert, 31 ch | `6b75e81f…` | confirmado |
| 7 | network metrics | `network_metrics.csv` | — | 2026-07-09 10:58 | strength top: C4 5.86 | (manifiesto) | confirmado |
| 7 | wPLI ALPHA | `wpli_ALPHA.csv` / `.png` | — | 2026-07-09 10:58 | 465 aristas | `1b1598f4…` (png) | confirmado |
| 7 | wPLI EC stale | `figures/wpli_*_EC.png` (7) | — | **2026-05-24** | — | — | ejecución anterior |

### Fase 8 — Surrogates e inferencia

| Fase | Output | Ruta | Fecha interna | Fecha FS | Config | Hash (prefijo) | Estado |
|------|--------|------|---------------|----------|--------|----------------|--------|
| 8 | surrogate summary | `surrogate_summary.json` | 2026-07-09T10:58:26 | 2026-07-09 10:58 | circular_shift, N=200, seed=42 | `d2b97df5…` | confirmado |
| 8 | significant edges | `significant_connections.csv` | — | 2026-07-09 10:58 | 447 total | (manifiesto) | confirmado |
| 8 | per-band p/q/sig | `wpli_pvalues_{band}.csv` etc. | — | 2026-07-09 10:58 | 7 bandas | (manifiesto) | confirmado |

† Config scratch original eliminada; `config_snapshot.toml` es la fuente reproducible.

---

## 7. Resultados reutilizados desde caché

| Artefacto | Fecha real ICA | Fecha corrida julio | Mecanismo | Validez |
|-----------|----------------|---------------------|-----------|---------|
| `results/M05/T2/cache/ica_result.jls` | 2026-05-24 18:47:49 | 2026-07-09 10:50:11 | `cache_valid` si `ica_config.hash` coincide | Válida por hash de config `[ica]` |
| Matrices W, A, activaciones S | Mayo | Exportadas julio | Deserialización + clasificación/rechazo julio | Coherente si config ICA no cambió |
| `ica_summary.json` `duration_s: 0.9` | — | Julio | Confirma carga, no cálculo FastICA | Alta confianza |

**Evidencia de carga desde caché:**

```
[10:50:11]   ICA cargado desde caché
```

(`pipeline_log.txt` L23)

**Implicación:** La narrativa «una sola ejecución de extremo a extremo el 2026-07-09» es correcta para filtrado→surrogates, pero **incorrecta si se omite que FastICA es de mayo**. El manual (`05_ica.tex` L166–168) y Slidev (diapositiva pipeline) documentan la caché.

---

## 8. Resultados pertenecientes a ejecuciones anteriores

### En `results/subjects/sub-M05/ses-T2/eyesclosed/figures/` (mtime 2026-05-24)

1. `signal_preview_EC.png`
2. `psd_all_channels_EC.png`
3. `band_power_summary_EC.png`
4. `wpli_ALPHA_EC.png`
5. `wpli_BETA_HIGH_EC.png`
6. `wpli_BETA_LOW_EC.png`
7. `wpli_BETA_MID_EC.png`
8. `wpli_DELTA_EC.png`
9. `wpli_GAMMA_EC.png`
10. `wpli_THETA_EC.png`

Documentado en `Report_Pre/INVENTARIO_SALIDAS_M05_2026-07-02.md` (nota 2026-07-09): las copias con sufijo `_EC` en `figures/` **no se regeneran** en cada corrida; las copias canónicas sin sufijo en la raíz sí.

### En `Report_Pre/figures/plots/` (copias para LaTeX)

| Archivo | mtime | vs results jul-09 |
|---------|-------|-------------------|
| `qc_raw_butterfly.png` | 2026-07-03 18:38 | results: 2026-07-09 10:50 |
| `qc_psd_raw_average.png` | 2026-07-03 18:38 | results: 2026-07-09 10:50 |
| `qc_variance_topomap.png` | 2026-07-03 18:17 | results: 2026-07-09 10:50 |
| `walkthrough_M05_raw_Fp2_Cz.png` | 2026-07-06 09:44 | `results/walkthrough_*`: 2026-07-09 11:07 |

### Legacy `results/M05/T2/`

Árbol paralelo para dashboard; marcado `legacy_paths` en manifiesto. Caché ICA vive aquí por diseño (`SingleSubjectPipeline.jl` L399–407).

### Cohorte batch 2026-05-26

199 grabaciones con `processed_at` 2026-05-26; config distinta (Fp2 excluido por defecto, sin surrogates en la mayoría). **No representan la config del informe julio.**

---

## 9. Figuras o tablas desactualizadas

| Artefacto | Usado en | Fecha en uso | Fecha canónica | Acción recomendada |
|-----------|---------|--------------|----------------|-------------------|
| `Report_Pre/figures/plots/qc_raw_butterfly.png` | Manual Cap. 4 (`03_raw.tex` L171) | 2026-07-03 | 2026-07-09 | Resincronizar desde results |
| `Report_Pre/figures/plots/qc_psd_raw_average.png` | Manual Cap. 4 (`03_raw.tex` L200) | 2026-07-03 | 2026-07-09 | Resincronizar |
| `figures/signal_preview_EC.png` | Slidev + manifiesto | 2026-05-24 | `filtered_signal_preview.png` 2026-07-09 | Actualizar manifiesto y sync |
| `slidev/public/.../signal_preview_EC.png` | Diapositiva «Señal raw» | 2026-07-13 (copia mayo) | jul-09 | Re-sync tras corregir fuente |
| `figures/*_EC.png` (10) | Dashboard / posible confusión | 2026-05-24 | raíz sin sufijo jul-09 | Archivar o regenerar |
| `channel_statistics_compare.csv` | Manifiesto `use: no_usar` | 2026-07-09 11:07 | N/A | No citar; localizar generador |

---

## 10. Inconsistencias del manual

| Afirmación del manual | Ubicación | Fuente real | Coincide | Clasificación |
|----------------------|-----------|-------------|----------|---------------|
| Fp2 único canal z>3 (z=4.63) | `03_raw.tex` Cuadro 5.1 | `channel_statistics.csv` | Sí | Exacta |
| 31 canales, Fp2 incluido | `config_snapshot` citado; `09_espectral.tex` L452–453 | `config_snapshot.toml [montage]` | Sí | Exacta |
| σ̄ raw 14.4 µV, sin amplitude_warning | `03_raw.tex` | `overview.csv` | Sí | Exacta |
| Figuras QC butterfly y PSD | `03_raw.tex` L171, L200 | `Report_Pre/figures/plots/qc_*` | Parcial | **Desactualizada** (copias jul-03) |
| Topografía varianza Fp2 aislado | `03_raw.tex` L560 | `figures/qc_variance_topomap.png` en results (jul-09) | Sí* | Correcta con contexto (*manual puede usar copia jul-03) |
| 5 IC rechazados, 82.6% varianza | `05_ica.tex` | `ica_summary.json` | Sí | Correcta con contexto (caché ICA mayo) |
| Reducción σ Fp2 −50.1% (100 s) | `05_ica.tex` Cuadro 8.6 L665 | `ica_headplot` / registro completo | Sí | Exacta |
| Caché ICA documentada | `05_ica.tex` L166–168 | `ica_result.jls` mayo | Sí | Exacta |
| 99 épocas válidas, retención 99% | `10_conclusiones.tex` L109 | `segmentation_summary.json` | Sí | Exacta |
| Época 12 rechazada, C4 peor canal | Cap. segmentación/AR | `rejected_segments.csv` | Sí | Exacta |
| Ventana Hamming (no Hanning) | `09_espectral.tex` L394 | `spectral_summary.json` `method: fft_hamming_taper` | Sí | Exacta |
| DELTA 30.99%, ALPHA 28.69% | `09_espectral.tex` | `spectral_summary.json` | Sí | Exacta |
| 447 conexiones significativas | `11_surrogates.tex` | `surrogate_summary.json` `n_sig_total: 447` | Sí | Correcta con contexto (verificación puntual) |
| Surrogates no en pipeline estándar | `10_conclusiones.tex` L171–174 | `config/single_subject.toml` `enabled=false` | Sí | Exacta |
| Validación BVA C4 DELTA ×18.48 | `09_espectral.tex` L927+ | Colaborador externo (J. Espuny) | — | **No verificable** |
| 78 sujetos con config anterior | `10_conclusiones.tex` L158–162 | batch 2026-05-26, `subjects_index.csv` | Sí | Exacta |

---

## 11. Inconsistencias de la presentación Slidev

| Diapositiva / bloque | Afirmación o cifra | Fuente declarada | Fuente real | Estado | Corrección propuesta |
|---------------------|-------------------|------------------|-------------|--------|---------------------|
| Portada / frontmatter | 2026-07-09 | manifiesto | `pipeline_log.txt` | Confirmado | — |
| Pipeline 8 fases | ICA «cargado desde caché» | log | log + cache mayo | Confirmado | — |
| Señal raw | preview 31 canales | `figures/signal_preview_EC.png` | mtime **2026-05-24** | **Inconsistente** | Usar `filtered_signal_preview.png` o regenerar |
| QC butterfly/PSD | rutas `results/.../figures/qc_*` | declarado | mtime 2026-07-09 | Confirmado | — |
| Fp2 z-score | 4.63 | `channel_statistics.csv` | verificado | Confirmado | — |
| Reducción ICA Fp2 | 42.0% (10 s) | `ica_signal_before/after.csv` | cálculo: 11.50→6.67 µV | Confirmado | — |
| Tabla σ canales | Fp2 −50.1% (100 s) | informe Cuadro 8.6 | cita manual, no recomputada | Correcta (cita) | Mantener caveat ventana |
| 447 conexiones sig. | ALPHA 258 + GAMMA 189 | `surrogate_summary.json` | verificado | Confirmado | — |
| Archivos con SHA-256 | N archivos | run-manifest | generado 2026-07-14 | Confirmado | Distinguir inventario vs ejecución |
| Limitación surrogates | verificación puntual | `subject-data.json` warnings | `config_snapshot` `enabled=true` | Confirmado | — |
| Advertencia 74.7% | no usar compare CSV | comentario slide ICA | `channel_statistics_compare` = raw→filt | Confirmado | — |

---

## 12. Diferencias manual – presentación – repositorio

| Tema | Manual | Presentación | Evidencia repositorio | Veredicto |
|------|--------|--------------|----------------------|-----------|
| Fecha última ejecución | Implícita cap. 11–12 | 2026-07-09 explícita | `pipeline_log.txt` 10:49:50 | Consistente |
| Fecha ICA real | Documentada en cap. ICA | Caché 2026-05-24 declarada | `ica_result.jls` 2026-05-24 | Consistente si se declara |
| Config activa | `config_snapshot.toml` | manifiesto YAML | snapshot julio | Consistente |
| Política Fp2 | exclude_fp2=false | igual | snapshot | Consistente |
| Figuras QC | `figures/plots/` jul-03 | `results/figures/` jul-09 | divergencia mtime | **Inconsistente (manual)** |
| signal_preview | no central | mayo 2026 | `filtered_signal_preview` jul-09 | **Inconsistente (slides)** |
| Reducción ICA Fp2 | 50.1% (100 s) | 42% (10 s) + tabla 50.1% | ambas reproducibles | Consistente con caveat |
| Duración épocas | 1.0 s | 1.0 s | `segmentation_summary.json` | Consistente |
| Número épocas | 99/100 | 99/100 | JSON + CSV | Consistente |
| Baseline | 0–0.10 s, 2 passes | igual | snapshot | Consistente |
| Ventana espectral | Hamming | Hamming | `spectral_summary.json` | Consistente |
| Método wPLI | Hilbert | Hilbert | snapshot + JSON | Consistente |
| Surrogates enabled | verificación puntual | advertencia | snapshot `true` vs default `false` | Consistente |
| Método surrogates | circular_shift | circular_shift | log + JSON; TOML `method` no leído | Consistente (código fijo) |
| Conexiones sig. | 447 | 447 | `surrogate_summary.json` | Consistente |
| Resultado verificación puntual | cap. surrogates | diapositiva fase 8 | único `surrogate_summary` reciente en cohorte | Consistente |

---

## 13. Hallazgos críticos

1. **Mezcla temporal ICA:** FastICA de **2026-05-24** + resto del pipeline de **2026-07-09**. Válido por hash, pero debe declararse siempre en informes.
2. **Manual PDF con figuras QC anteriores a la corrida canónica:** `\graphicspath{{figures/plots/}}` sirve copias del **2026-07-03**, no del **2026-07-09**.
3. **Slidev muestra `signal_preview_EC.png` de mayo** en la diapositiva «Señal raw», sincronizada el 2026-07-13 sin validar mtime de origen.
4. **10 figuras `*_EC.png` obsoletas** permanecen en la carpeta BIDS junto a outputs actuales — riesgo de uso accidental (dashboard).
5. **`channel_statistics_compare.csv` sin procedencia en código** — la columna `std_uv_filt` produce un 74.7% que NO es reducción ICA (es raw→filtrado).
6. **Config scratch eliminada** — solo `config_snapshot.toml` permite reproducir parámetros; `git_commit` y `config_hash` del manifiesto son `null`.
7. **Desfase código ↔ log ICA (posible):** el log del 2026-07-09 registra «Componentes rechazados (auto, umbral=1.5)»; el código vigente en `SingleSubjectPipeline.jl` usa `load_ica_labels()` (solo CSV manual) y el mensaje de log no incluye «auto». `README.md` L1174 documenta esta discrepancia. **Re-ejecutar hoy podría no reproducir el mismo rechazo** si no existe `ica_labels.csv`.

---

## 14. Hallazgos menores

1. `filtered_signal_preview.png` (jul-09) no está en el manifiesto YAML ni en el inventario SHA-256.
2. `walkthrough_M05_raw_*.png` en `Report_Pre/figures/plots/` (jul-06) vs `results/` (jul-09 11:07) — fechas distintas, contenido probablemente similar pero no verificado byte a byte.
3. `qc_amplitude_histograms.png` en results tiene mtime **2026-07-03** (único QC no regenerado jul-09); subpáginas `_p01/_p02` sí son jul-09/07.
4. `docs/REPORT_PHASE1_AUDIT.md` referencia rutas legacy `report/main_es.tex` — no aplica a `Report_Pre/`.
5. M05 EC fue **SKIPPED** en batch 2026-05-26; la entrada en `subjects_index.csv` solo refleja la corrida julio.
6. `band_topomap_grid.png` en `figures/` tiene mtime **2026-07-08** (un día antes de la corrida principal).

---

## 15. Elementos no verificables

| Elemento | Motivo |
|----------|--------|
| Validación BrainVision Analyzer (correlaciones DELTA, C4 ×18.48) | Sin CSV/export BVA en `results/subjects/sub-M05/` |
| Generador de `channel_statistics_compare.csv` | No hay referencia en `src/` ni `scripts/` |
| `git_commit` de la corrida | `null` en manifiesto; no registrado en pipeline |
| Identidad byte-a-byte de figuras walkthrough jul-06 vs jul-09 | No se calculó diff de imágenes |
| Contenido visual exacto del PDF manual compilado 2026-07-14 | No se realizó revisión visual página a página |
| Si la caché ICA de mayo usó exactamente la misma señal filtrada que julio | Misma config declarada; señal filtrada no almacenada en mayo para diff directo |

---

## 16. Correcciones mínimas recomendadas

### Documentales (prioridad alta)

1. Resincronizar `Report_Pre/figures/plots/qc_raw_butterfly.png` y `qc_psd_raw_average.png` desde `results/subjects/sub-M05/ses-T2/eyesclosed/figures/` (jul-09).
2. Actualizar manifiesto Slidev: sustituir `figures/signal_preview_EC.png` por `filtered_signal_preview.png` o regenerar preview con nombre estable.
3. Ejecutar `sync-subject-assets.mjs` tras corregir manifiesto; re-exportar `slides-export.pdf`.
4. Recompilar `main.pdf` tras sync de figuras QC.
5. Añadir nota al pie en Cap. 4: «Figuras sincronizadas desde corrida 2026-07-09».

### Organización

1. Introducir `run_id` (UUID o timestamp) escrito en cada `*_summary.json`.
2. Script `sync-report-figures.mjs` espejo de Slidev para el manual LaTeX.
3. Política post-corrida: eliminar o mover a `legacy/` los `figures/*_EC.png` no regenerados.
4. Registrar `git_commit` y hash de config al finalizar pipeline.

### Código (solo si se confirma desfase ICA)

1. `SingleSubjectPipeline.jl`: reintegrar rechazo automático por `artifact_threshold` si debe coincidir con log 2026-07-09, o corregir README si el log es histórico.
2. Al guardar preview filtrado, actualizar también `figures/signal_preview_EC.png` o deprecar ese nombre.

---

## 17. Comprobaciones que requerirían nueva ejecución

| Comprobación | Motivo | Autorización |
|--------------|--------|--------------|
| Invalidar caché ICA y re-ejecutar fase 4 | Alinear FastICA con fecha corrida | **Sí** |
| Regenerar `figures/*_EC.png` o eliminarlas | Eliminar ambigüedad dual-layout | **Sí** |
| Re-ejecutar pipeline M05 completo con config snapshot | Verificar reproducibilidad end-to-end | **Sí** |
| Diff visual figuras QC jul-03 vs jul-09 | Cuantificar impacto en manual PDF | No (solo lectura posible) |
| Batch Phase C cohorte completa | 78 sujetos con config vigente | **Sí** (estratégica) |

**No se ejecutó ninguna de estas comprobaciones durante esta auditoría.**

---

## 18. Evidencias y rutas exactas

### Corrida canónica

```
results/subjects/sub-M05/ses-T2/eyesclosed/
├── pipeline_log.txt
├── config_snapshot.toml
├── overview.csv
├── channel_statistics.csv
├── ica_summary.json
├── segmentation_summary.json
├── spectral_summary.json
├── connectivity_summary.json
├── surrogate_summary.json
├── ica_signal_before.csv
├── ica_signal_after.csv
├── rejected_segments.csv
└── figures/  (mezcla jul-09 + 10 archivos *_EC.png de mayo)
```

### Caché ICA

```
results/M05/T2/cache/
├── ica_result.jls      (2026-05-24 18:47:49, 12.5 MB)
└── ica_config.hash     (contenido: -26e3c03259100db8)
```

### Publicaciones

```
Report_Pre/main.tex
Report_Pre/main.pdf
Report_Pre/figures/plots/          ← copias manual (varias desactualizadas)
Report_Pre/slidev/slides.md
Report_Pre/slidev/slides-export.pdf
Report_Pre/slidev/manifests/runs/sub-M05_ses-T2_eyesclosed.yml
Report_Pre/slidev/generated/run-manifest/sub-M05_ses-T2_eyesclosed.json
Report_Pre/slidev/generated/subjects/M05/subject-data.json
Report_Pre/slidev/public/subjects/M05/figures/
```

### Métricas verificadas (cálculo directo sobre CSV)

| Métrica | Valor | Fuente |
|---------|-------|--------|
| Fp2 z-score RMS | 4.63 | `channel_statistics.csv` |
| Fp2 std raw | 54.58 µV | `channel_statistics.csv` |
| Reducción raw→filtrado Fp2 | 74.7% | `channel_statistics_compare.csv` (raw 54.58 → filt 13.83) |
| Reducción ICA aislada Fp2 (10 s) | 42.0% | `ica_signal_before/after.csv` (11.50 → 6.67 µV) |
| Reducción σ Fp2 registro completo | 50.1% | `05_ica.tex` Cuadro 8.6 (13.83 → 6.91 µV) |
| Épocas válidas | 99/100 | `segmentation_summary.json` |
| Conexiones significativas | 447/465 | `surrogate_summary.json` |
| ALPHA significativas | 258 | `surrogate_summary.json` |
| GAMMA significativas | 189 | `surrogate_summary.json` |
| Strength máximo | C4 = 5.86 | `network_metrics.csv` |

---

## Plan de corrección (Fase 7)

### A. Correcciones documentales

| Acción | Motivo | Riesgo | Coste | Prioridad | Requiere autorización |
|--------|--------|--------|-------|-----------|----------------------|
| Resincronizar QC plots al manual | PDF usa figuras jul-03 | Bajo | 15 min | Alta | No |
| Corregir fuente signal_preview en Slidev | Figura de mayo en diapositiva raw | Bajo | 30 min | Alta | No |
| Notas de pie fecha corrida + caché ICA | Trazabilidad | Nulo | 1 h | Media | No |
| Distinguir 42% vs 50.1% en captions | Evitar confusión | Nulo | 30 min | Media | No |

### B. Correcciones de organización

| Acción | Motivo | Riesgo | Coste | Prioridad | Requiere autorización |
|--------|--------|--------|-------|-----------|----------------------|
| `run_id` en summaries | Anti-mezcla | Bajo | 2–4 h dev | Alta | No |
| `sync-report-figures.mjs` | Automatizar copias manual | Bajo | 2 h | Media | No |
| Archivar `*_EC.png` obsoletos | README advierte pero no enforced | Bajo | 30 min | Alta | Sí (borrado) |
| `git_commit` en manifiesto al ejecutar | Procedencia | Bajo | 1 h | Media | No |

### C. Correcciones de código

| Archivo | Función / zona | Problema | Prioridad |
|---------|----------------|----------|-----------|
| `src/SingleSubjectPipeline.jl` | ICA rechazo ~L434 | Log jul-09 dice auto; código usa solo `load_ica_labels` | Alta |
| `src/SingleSubjectPipeline.jl` | `_save_ica_results` ~L1552 | `artifact_thresh = 1.5` hardcoded, no lee `cfg` | Media |
| `src/SingleSubjectPipeline.jl` | guardado figuras | No actualiza `figures/signal_preview_EC.png` | Media |
| `src/statistics/Surrogates.jl` | `surrogate_test` | `cfg.surrogates["method"]` no se lee | Baja (documentado) |

### D. Ejecuciones necesarias

| Acción | Motivo | Riesgo | Coste | Prioridad | Requiere autorización |
|--------|--------|--------|-------|-----------|----------------------|
| Borrar caché ICA + re-run fase 4 M05 | Alinear ICA con julio | Medio | ~5 min | Media | **Sí** |
| Pipeline completo M05 desde snapshot | Regenerar todos los outputs | Bajo | ~9 min | Media | **Sí** |
| Re-export PDF manual + slides | Tras sync figuras | Nulo | ~15 min | Alta | No |
| Batch Phase C | Cohorte con config vigente | Bajo científico | 6–34 h | Baja | **Sí** |

---

## Preguntas que requieren decisión humana

1. ¿Se acepta la **caché ICA de mayo** como válida para el informe, o se exige re-ejecución de FastICA?
2. ¿Figura canónica de preview: `filtered_signal_preview.png` (jul-09) o regenerar `signal_preview_EC.png`?
3. ¿Eliminar los 10 `figures/*_EC.png` de mayo o archivarlos bajo `legacy/`?
4. ¿Quién generó `channel_statistics_compare.csv` y debe incorporarse al repo o descartarse definitivamente?
5. ¿El manual debe referenciar figuras vía symlink a `results/` o mantener copias con script de sync obligatorio?

---

*Fin del informe. Generado en modo solo lectura. No se modificaron manual, Slidev, código, resultados ni timestamps.*
