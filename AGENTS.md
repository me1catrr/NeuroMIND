# AGENTS.md — NeuroMIND (OpenAI Codex CLI)

> Leído automáticamente por Codex CLI al iniciar una sesión en este directorio.
> Contexto completo del proyecto para agentes OpenAI. Última actualización: 2026-05-21.

---

## Proyecto

**NeuroMIND** — Framework Julia para análisis de conectividad funcional EEG en Esclerosis Múltiple.

- **Repo**: https://github.com/me1catrr/NeuroMIND
- **Desarrollador**: Rafael Castro Triguero · me1catrr@uco.es · https://github.com/me1catrr
- **Lenguaje**: Julia 1.9+ (no Python, no R)
- **Dashboard**: web local Genie.jl · `http://localhost:8080`
- **Proyecto relacionado**: `../EEG_Julia/` — original del que NeuroMIND es refactorización

## Reglas estrictas — NUNCA hacer

```
❌ No subir data/, results/, reports/, exports/ al repo
❌ No subir .claude/, .cursor/, .vscode/, .DS_Store
❌ No subir señales EEG reales, derivados por sujeto, logs clínicos
❌ No commitear directamente en main (salvo cambios triviales de 1 línea)
❌ No usar MultivariateStats — ICA implementado en puro Julia (LinearAlgebra + Random)
❌ No cambiar lógica científica (ICA, wPLI, PSD) sin comparar contra EEG_Julia/
```

## Git identity

```bash
git config --global user.name "Rafael Castro Triguero"
git config --global user.email "me1catrr@uco.es"
```

## Pipeline — 8 pasos

```
[1/8] Carga BIDS       → EEGRecording
[2/8] QC canales       → flag por z-score
[3/8] Filtrado         → HP 0.5 Hz + LP 48 Hz + Notch 50 Hz (Butterworth filtfilt)
[4/8] ICA              → FastICA sobre señal CONTINUA filtrada (ANTES de segmentar)
[5/8] Segmentación     → sobre señal limpiada por ICA + baseline + rechazo AR
[6/8] Espectral PSD    → Hanning + potencia por banda
[7/8] wPLI             → Hilbert analítica across-segments (CSD opcional, off por defecto)
[8/8] Guardado         → results/subjects/sub-{id}/ses-{sess}/{task}/
```

**Regla crítica**: ICA (paso 4) SIEMPRE antes de segmentar (paso 5).

## Archivos clave

| Archivo | Rol |
|---------|-----|
| `src/NeuroMIND.jl` | Entry point del módulo |
| `src/types.jl` | Tipos: EEGRecording, EpochSet, ICAResult, SpectralResult, ConnectivityMatrix |
| `src/SingleSubjectPipeline.jl` | Pipeline completo 8 pasos |
| `src/ica/ICACore.jl` | FastICA puro Julia (PCA whitening + tanh) |
| `src/webapp/App.jl` | Servidor Genie + rutas API |
| `web/views/dashboard.html` | SPA dashboard (~6300 líneas, HTML+CSS+JS inline) |
| `config/single_subject.toml` | Parámetros del pipeline |
| `scripts/run_single_subject.jl` | CLI del pipeline |

## Comandos frecuentes

```bash
# Verificar sintaxis
julia --project=. -e 'include("src/NeuroMIND.jl"); println("OK")'

# Tests
julia --project=. tests/runtests.jl

# Pipeline
julia --project=. scripts/run_single_subject.jl

# Dashboard
julia --project=. scripts/launch_dashboard.jl

# Verificar antes de push (debe devolver vacío)
git ls-files | grep -E '(^data/|^results/|\.DS_Store$|^\.claude/)'
```

## Dashboard — 13 paneles (Phase 0–12)

Implementados: 0, 1, 2, 3, 4, 5. Paneles 6–12: placeholders vacíos.

## Dependencias Julia

```
Base64, CSV, CairoMakie, DSP, DataFrames, Dates, FFTW,
Genie, LinearAlgebra, Random, Serialization, Statistics, StatsBase, TOML
```
Sin: MultivariateStats, JSON3, Pluto.

## Estado actual (2026-05-21)

### Hecho
- Pipeline 8 pasos con ICA propio (puro Julia)
- Dashboard paneles 0–5
- API ICA: `/api/phase5_ica_info`, `/api/ica_activation`, `/api/ica_signal`
- README.md, CLAUDE.md, AGENTS.md, .cursorrules

### Pendiente
- Dashboard paneles 6–12
- Clasificación automática de componentes ICA
- Modo reproducibilidad EEG_Julia (CSD + mismos parámetros)
- GitHub Actions CI

## Workflow Git

```bash
git fetch origin
git checkout -b feat/<nombre>
# ... trabajar ...
git add <archivos específicos>
git commit -m "feat: descripción"
git push -u origin feat/<nombre>
gh pr create --title "..." --body "..."
```

> Para contexto completo, ver también `CLAUDE.md` en este mismo directorio.
