# ROADMAP — Análisis Completo MINDEM-IMIBIC

## Dataset: 41 pacientes EM + 37 controles sanos  
**Última actualización:** 2026-05-24  
**Estado general:** Fase A (Auditoría) → Fase B (BIDS) en marcha

---

## 1. Cohort Summary

| Grupo | N | Sesiones | EC+EO | Análisis |
|-------|---|----------|-------|---------|
| MS (PAC) | 41 (M4–M44) | T1, T2 (≈29 pares) | ✓ | Transversal + Longitudinal |
| HC (CON) | 37 (MC1–MC40) | Única (sin T1/T2) | ✓ | Transversal |

**Sujetos excluidos:**
- M4: solo Oddball en T1, sin T2 → excluido
- M5: T1 solo Oddball → excluido de longitudinal; T2 resting state válido
- M6: solo Oddball en T1, sin T2 → excluido
- MC5, MC6, MC11: en carpeta EXCLUIDOS
- M3_T2: en carpeta EXCLUIDOS

**Sujetos solo T1 (sin longitudinal):** M10, M17, M26–M29, M32, M35, M38, M41  
**Duplicado T2:** M37 → usar T2_1; T2_2 como respaldo  

---

## 2. Fases del Pipeline

### Fase A — Auditoría del Dataset ✅ IMPLEMENTADO
**Script:** `scripts/audit_full_dataset.jl`  
**Salida:** `data/full_data/inventory.csv`

- Parseo de todos los `.vhdr` en `data/full_data/Pacientes MINDEM_IMIBIC_27 03 25/`
- Normalización de nombres de condición (>15 variantes ortográficas)
- Clasificación: EC / EO / ODDBALL / EXCLUIDO
- Detección de pares T1/T2 para análisis longitudinal
- Generación de `participants.tsv`, `groups.csv`, `longitudinal_pairs.csv`

---

### Fase B — Conversión BIDS (Ligera) ✅ IMPLEMENTADO
**Script:** `scripts/build_bids_full.jl`  
**Salida:** `data/bids/raw/`, `data/bids/electrodes/`

Crea estructura BIDS **sin copiar datos binarios** (referencia al archivo original):
- `sub-{id}/ses-{sess}/eeg/metadata.json` con `data_format: "brainvision"`, `vhdr_path`
- `sub-{id}/ses-{sess}/electrodes.tsv` (copia del template con 31 canales 10-20)
- `participants.tsv`, `groups.csv`, `longitudinal_pairs.csv`
- `dataset_description.json`

---

### Fase C — Pipeline Individual por Sujeto ✅ DISPONIBLE
**Script:** `scripts/run_batch_pipeline.jl`  
**Config:** `config/batch_pipeline.toml`

Para cada sujeto/sesión/condición válida en `inventory.csv`:
1. Lee metadata JSON de BIDS (con `vhdr_path`)
2. Carga señal BrainVision binario via `BrainVisionLoader`
3. Ejecuta `run_single_subject_pipeline` (pasos 1-8: QC, filtrado, ICA, segmentación, AR, espectral, wPLI, surrogates)
4. Guarda en `results/subjects/sub-{id}/ses-{sess}/{task}/`
5. Registra en `logs/batch_run_YYYY-MM-DD.csv`

**Estimación de tiempo:** ~5-15 min por sujeto × ~210 sujetos/sesiones = 18-50 horas  
**Paralelización recomendada:** procesamiento por lotes (ECvsEO, T1 vs T2)

---

### Fase D — Análisis Transversal (MS vs Controles) 
**Script existente:** `scripts/run_transversal_analysis.jl`  
**Entrada:** `data/bids/groups.csv` + resultados individuales

Para cada banda y condición:
- Welch t-test con corrección BH-FDR (q < 0.05)
- Cohen's d por par de electrodos
- Mapas de conectividad por grupo
- Salida: `results/group/transversal/{EC,EO}/`

Dashboard: Panel 13

---

### Fase E — Análisis Longitudinal (T1 → T2 en EM)
**Script existente:** `scripts/run_longitudinal_analysis.jl`  
**Entrada:** `data/bids/longitudinal_pairs.csv`

Para cada paciente con T1+T2:
- Paired t-test + Cohen's dz + BH-FDR
- Evolución temporal de conectividad por banda
- Salida: `results/group/longitudinal/{EC,EO}/`

Dashboard: Panel 14

---

### Fase F — Replicación MNE-Python
**Carpeta:** `EEG_MNE_Python/`  
**Propósito:** Replicar el pipeline NeuroMIND en Python/MNE para validación cruzada

Pasos a implementar:
1. `EEG_MNE_Python/pipeline_mne.py` — pipeline completo MNE
2. `EEG_MNE_Python/compare_results.py` — comparación Julia vs Python
3. Panel 15 en dashboard: comparativa de métricas wPLI Julia/Python

---

### Fase G — Enriquecimiento del Informe
- Demografías completas en informe PDF/HTML
- Estadísticas grupo (edad, sexo, educación)
- Correlaciones clínicas (si se dispone de escalas EDSS, cognitivas)
- Figuras de publicación para `wPLI_paper/`

---

## 3. Estructura de Archivos Generados

```
data/
├── full_data/
│   ├── inventory.csv                    ← Fase A: inventario completo
│   └── Pacientes MINDEM_IMIBIC_27 03 25/  ← datos originales (no subir al repo)
└── bids/
    ├── dataset_description.json
    ├── participants.tsv
    ├── groups.csv
    ├── longitudinal_pairs.csv
    ├── electrodes/
    │   └── sub-{id}_ses-{sess}_electrodes.tsv   ← posiciones 10-20
    └── raw/
        └── sub-{id}_ses-{sess}_task-{task}_run-01_eeg_metadata.json

results/
├── subjects/
│   └── sub-{id}/ses-{sess}/{task}/
│       ├── overview.csv
│       ├── wpli_{band}.csv
│       ├── surrogate_summary.json
│       └── ...  (todos los archivos del pipeline individual)
└── group/
    ├── transversal/{EC,EO}/
    │   ├── band_statistics.csv
    │   ├── group_statistics_{band}.csv
    │   └── transversal_summary.json
    └── longitudinal/{EC,EO}/
        ├── band_statistics_longitudinal.csv
        ├── longitudinal_statistics_{band}.csv
        └── longitudinal_summary.json

logs/
└── batch_run_YYYY-MM-DD.csv            ← registro de procesamiento en lote
```

---

## 4. Orden de Ejecución

```bash
# 1. Auditar datos
julia --project=. scripts/audit_full_dataset.jl

# 2. Crear estructura BIDS
julia --project=. scripts/build_bids_full.jl

# 3. Lanzar pipeline en lote (puede dividirse por sesiones)
julia --project=. scripts/run_batch_pipeline.jl --condition EC
julia --project=. scripts/run_batch_pipeline.jl --condition EO

# 4. Análisis transversal
julia --project=. scripts/run_transversal_analysis.jl

# 5. Análisis longitudinal
julia --project=. scripts/run_longitudinal_analysis.jl

# 6. Dashboard para visualizar resultados
julia --project=. scripts/launch_dashboard.jl
```

---

## 5. Control de Calidad del Dataset

| Sujeto | T1_EC | T1_EO | T2_EC | T2_EO | Estado |
|--------|-------|-------|-------|-------|--------|
| M4 | ✗ | ✗ | — | — | EXCLUIDO (solo Oddball) |
| M5 | ✗ | ✗ | ✓ | ✓ | Solo T2 |
| M6 | ✗ | ✗ | — | — | EXCLUIDO |
| M7-M9 | ✓ | ✓ | ✓ | ✓ | T1+T2 |
| M10 | ✓ | ✓ | — | — | Solo T1 |
| M11-M25 | ✓ | ✓ | ✓* | ✓* | T1+T2 (*ver notas) |
| M26-M35 | ✓ | ✓ | — | — | Solo T1 |
| M36-M44 | ✓ | ✓ | ✓ | ✓ | T1+T2 |
| MC* (excl MC5,6,11) | ✓ | ✓ | — | — | Solo sesión única |

*M37 tiene T2_1 y T2_2 duplicados → usar T2_1 como T2 primario

---

## 6. Notas de Implementación

### BrainVision Loader
- Archivo binario `.eeg`: IEEE_FLOAT_32, MULTIPLEXED (31ch × N_muestras × 4 bytes)
- No se escala (resolución 0.0488281 µV/bit sólo para INT_16; IEEE_FLOAT_32 ya en µV)
- Canal de referencia: Cz (implícito)
- 31 canales EEG + FCz (REF) + Fpz (GND) en el `.vhdr`, pero el `.eeg` tiene 31

### Normalización de Condiciones
Las siguientes variantes se normalizan a EC/EO:
- EC: ojoscerrados, Ojoscerrados, OJOSCERRADOS, Ojoscerradps, ojoscerreados, 
      ojoscerradpsb, OjoscCerrados, "Ojos cerrados", "ojos cerrados", ojosCerrados,
      OjosCerrados, ojocerrados
- EO: ojosabiertos, Ojosabiertos, OJOSABIERTOS, OJOSACIERTOS, ojosabioertos,
      ojosabiertes, "Ojos abiertos", "OJOS ABIERTOS", "OJOS ABIERTTOS",
      ojoscabiertos, OjosAbiertos, "ojos abietos"

### Longitudinal Pairs (N≈28 pares MS con T1+T2)
M7, M8, M9, M11, M12, M13, M14, M15, M16, M18, M19, M20, M21, M22,
M23, M24, M25, M30, M31, M33, M34, M36, M37, M39, M40, M42, M43, M44

### Seguridad
- Datos EEG reales NO van al repo (`data/`, `results/`, `logs/` en `.gitignore`)
- Los scripts producen artefactos en `data/bids/` y `results/` (ambos en `.gitignore`)
- Los scripts y configs SÍ van al repo
