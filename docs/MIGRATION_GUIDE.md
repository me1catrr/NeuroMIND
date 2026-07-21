# Guía de Migración: EEG_Julia → NeuroMIND

> **Nota (2026-07-21):** este documento describe la migración histórica
> EEG_Julia → NeuroMIND. El ejemplo de la Fase 1 usa `config/pipeline.toml`
> y `load_config`/`run_pipeline!`, que pertenecen a la ruta legacy archivada
> en `deprecated/code/config/` y `deprecated/code/src/Pipeline.jl`. La ruta activa actual usa
> `config/single_subject.toml` y `SingleSubjectPipeline.jl` — ver README.md.

## Estrategia

La migración es **incremental**: EEG_Julia sigue funcionando intacto mientras
NeuroMIND se desarrolla en paralelo. Se migra módulo a módulo, validando que
los resultados numéricos son idénticos antes de retirar el código antiguo.

---

## Fase 1: Validación de equivalencia numérica

Antes de cualquier migración, confirmar que NeuroMIND reproduce exactamente
los resultados de EEG_Julia para el dataset de referencia (M05/T2/EC).

```julia
# En EEG_Julia (resultados de referencia)
# Archivos en: EEG_Julia/data/Connectivity/wPLI/dict_wpli.bin

# En NeuroMIND
using NeuroMIND
cfg  = load_config("config/pipeline.toml")
# ... ejecutar pipeline ...

# Comparar matrices wPLI ALPHA
W_old = load_legacy_wpli("../data/Connectivity/wPLI/dict_wpli.bin", "ALPHA")
W_new = conn.matrices["ALPHA"]
@assert maximum(abs.(W_old .- W_new)) < 1e-6 "wPLI ALPHA no coincide"
```

---

## Fase 2: Migración módulo a módulo

### Prioridad 1 — I/O y configuración
- [ ] Migrar `config/default_config.jl` → `config/pipeline.toml`
- [ ] Migrar `src/modules/paths.jl` → `src/io/Config.jl`
- [ ] Verificar que `load_eeg_bids` lee el mismo formato TSV

### Prioridad 2 — Preprocessing
- [ ] Migrar `src/Preprocessing/filtering.jl` → `src/preprocessing/Filtering.jl`
- [ ] Comparar señales filtradas sample a sample

### Prioridad 3 — ICA
- [ ] Migrar `src/ICA/ICA.jl` + `ICA_cleaning.jl`
- [ ] Exportar labels ICA existentes como CSV para `load_ica_labels`

### Prioridad 4 — Spectral
- [ ] Migrar `src/Spectral/FFT.jl` → `src/spectral/PowerSpectrum.jl`
- [ ] Verificar potencia por banda contra `data/Processing/FFT/dict_FFT_power.bin`

### Prioridad 5 — Conectividad (ya portada)
- [ ] `src/Connectivity/wPLI.jl` → `src/connectivity/wPLI.jl` ✓
- [ ] `src/Connectivity/CSD.jl` → `src/connectivity/CSD.jl` (opcional) ✓

### Prioridad 6 — Estadística
- [ ] Migrar lógica surrogate de `Surrogate.jl`
- [ ] Añadir FDR que faltaba en el sistema antiguo

---

## Fase 3: Dashboards Pluto

Los notebooks Pluto de EEG_Julia se reemplazan por dashboards NeuroMIND:

| EEG_Julia Pluto | → | NeuroMIND Dashboard |
|----------------|---|---------------------|
| `Pluto/BIDS/BIDS.jl` | → | QC implícito en pipeline |
| `Pluto/Preprocessing/Preprocessing.jl` | → | Ejecutado automáticamente |
| `Pluto/ICA/ICA.jl` | → | `dashboards/ica.jl` (pendiente) |
| `Pluto/Spectral/Spectral.jl` | → | `dashboards/spectral.jl` ✓ |
| `Pluto/Connectivity/Connectivity.jl` | → | `dashboards/connectivity.jl` ✓ |
| `Pluto/Surrogate/Surrogate.jl` | → | `dashboards/statistics.jl` (pendiente) |

---

## Fase 4: Eliminar EEG_Julia

Solo cuando:
1. Todos los tests de equivalencia pasan
2. Los dashboards están validados con datos reales
3. Los resultados publicados siguen siendo reproducibles

---

## Qué NO migrar

- `tools/pages/` — herramientas de GitHub Pages: son independientes
- `Javier_results/` — datos de referencia: copiar como fixtures de test
- `data/` — cache de resultados: se regeneran con NeuroMIND

---

## Notas de diseño

### Por qué cambiar la serialización

EEG_Julia usa `Serialization.serialize` directamente, sin versionado.
Un cambio de versión de Julia rompe todos los `.bin`.
NeuroMIND envuelve cada resultado en un `Dict` con versión del framework,
detectando incompatibilidades antes de cargar datos corruptos.

### Por qué TOML en lugar de `default_config.jl`

Un archivo `.jl` de configuración es código ejecutable:
difícil de parametrizar, testear y versionar.
Un TOML es datos puros: parseado sin efectos secundarios,
editable sin conocer Julia, y difable claramente en git.

### Por qué estructuras en lugar de Dicts

`dict_csd["eeg_csd"]` no es autoexplicativo.
`epochs.data` sí lo es. Los tipos Julia permiten al compilador
generar código eficiente y al IDE ofrecer autocompletado.
