# Plan de continuación — Report_Pre + revisión de código/resultados

> Creado: 2026-10-07 (cierre de sesión en Claude Code). Pensado para seguir en
> **Cursor** (o cualquier agente): lee también `AGENTS.md` (§10 estado, §12 reglas).

---

## 0. Estado al cerrar

| Elemento | Estado |
|---|---|
| 🟡 PR [#5](https://github.com/me1catrr/NeuroMIND/pull/5) `feat/graphify-setup` → `main` | Abierto, *mergeable*. 21 commits (julio + graphify). Revisar y fusionar. |
| 🟡 PR [#6](https://github.com/me1catrr/NeuroMIND/pull/6) `feat/ci` → `feat/graphify-setup` | Abierto. Primera ejecución del CI en curso al cerrar; comprobar resultado en la pestaña *Checks*. Fusionar **después** de #5. |
| 🟢 Tests locales | 449/449 pass (`julia --project=. test/runtests.jl`). |
| 🟢 graphify | Grafo versionado en `graphify-out/`; hooks git activos en el iMac. |
| ⚪ Mac Studio | Pendiente: aplicar `graphify-out/patches/graphify_julia_returntype.patch` y `graphify hook install`. |
| 🟡 `Report_Pre/` (repo git propio) | Cambios sin commitear en `main.tex`, `appendices/A_referencia.tex`, caps. 00–02. |
| 🔴 `Report_Pre/figures/plots/` | Desincronizado de `results/` desde 2026-07-27 (sync manual, ver AGENTS.md §9). |

> [!IMPORTANT]
> Si el CI de #6 falla, el log de GitHub Actions dice en qué paso. Lo más probable
> es un test que dependa de archivos locales no versionados: hacerlo *skip* si falta
> el archivo (mismo patrón que las regresiones ALPHA–EC en `test/runtests.jl`).

---

## 1. Cómo trabajar en Cursor

- **Reglas ya activas:** `AGENTS.md` (contexto del proyecto) y `.cursor/rules/graphify.mdc`
  (el agente consulta el grafo antes de leer código).
- **Orientarse sin gastar contexto:** antes de abrir archivos grandes,
  ```bash
  graphify query "¿qué funciones calculan la potencia por banda?"
  graphify explain "segment_recording()"
  graphify path "run_single_subject_pipeline()" "compute_wpli()"
  ```
  Limitación: no ve llamadas a funciones exportadas vía `using NeuroMIND`
  (p. ej. desde `test/`): *sin aristas ≠ código muerto*, confirmar con búsqueda.
- **Grafo al día:** el hook post-commit re-extrae el código solo. Si cambian `.md`
  de NeuroMIND, la parte semántica requiere `/graphify . --update` (Claude Code) o
  pedírselo al agente de Cursor.
- **Ritmo acordado para el informe:** capítulo a capítulo; primero auditar y
  presentar hallazgos (tabla con 🔴🟡🟢), **sin tocar nada**; decidir; y solo
  entonces editar. Si hay que dejar constancia, `.md` de auditoría dentro de la
  carpeta revisada.

---

## 2. Bucle por capítulo (Report_Pre ↔ código ↔ resultados)

Para cada capítulo, en este orden:

1. **Leer el capítulo** (`Report_Pre/chapters/NN_*.tex`) y listar qué afirma:
   parámetros, figuras, tablas, listings.
2. **Contrastar con el código** del módulo correspondiente (tabla §3), usando
   `graphify explain` para orientarse y `config/pipeline.toml` como fuente de
   parámetros.
3. **Contrastar con los resultados** (`results/subjects/sub-M05/ses-T2/eyesclosed/`
   para el caso único; `results/{transversal,longitudinal}/eyesclosed/` para grupo).
4. **Clasificar discrepancias:** 🔴 error en código/resultado · 🟡 texto o figura
   desactualizada · 🟢 coherente.
5. **Corregir** en la rama adecuada:
   - código NeuroMIND → `feat/<tema>` + syntax check + tests;
   - informe → commit en el repo `Report_Pre/`.
6. **Figuras:** regenerar en NeuroMIND, luego copiar a `Report_Pre/figures/plots/`.

> [!WARNING]
> Nunca cambiar lógica científica (ICA, wPLI, PSD, estadística) sin contrastar con
> `../EEG_Julia/` (AGENTS.md §12).

---

## 3. Mapa capítulo → código → salidas

| Cap. | Tema | Código principal | Salidas a revisar |
|---|---|---|---|
| 04 | Diseño y adquisición | `scripts/audit_full_dataset.jl`, `src/io/BrainVisionLoader.jl` | `qc_decision_table.csv`, inventario |
| 05 | Señal raw | `src/io/`, `src/qc/QualityControl.jl` | `raw_signal.csv`, `qc_summary.csv` |
| 06 | Filtrado | `src/preprocessing/Filtering.jl` | figuras de filtrado |
| 07–08 | Artefactos / ICA | `src/ica/ICACore.jl`, `ICAClassification.jl` | `ica_summary.json`, `ica_signal_*.csv` |
| 09–11 | Segmentación / baseline / AR | `src/segmentation/Epochs.jl` | nº épocas válidas, rechazo ±70 µV |
| 12 | Espectral | `src/spectral/PowerSpectrum.jl` | `band_power_summary.csv` |
| 13 | Conectividad | `src/connectivity/wPLI.jl` | `wpli_{band}.csv`, `connectivity_edges.csv` |
| 14 | Surrogates | `src/statistics/Surrogates.jl` | `surrogate_summary.json`, `wpli_qvalues_*` |
| 15 | Transversal | `src/transversal/Transversal.jl` | `results/transversal/eyesclosed/` |
| 16 | Longitudinal | `src/longitudinal/Longitudinal.jl` | `results/longitudinal/eyesclosed/` |

---

## 4. Cola de trabajo (orden sugerido)

1. 🟡 Fusionar #5 y #6 (tras CI verde); configurar el Mac Studio (§0).
2. 🟡 Commitear o descartar los cambios pendientes de `Report_Pre/` (caps. 00–02).
3. 🔴 **Sincronizar figuras** `results/` → `Report_Pre/figures/plots/` y recompilar
   (`make`); revisar caps. 15–16, que citan resultados de grupo.
4. 🔴 **Decidir `GroupStats.jl` / `FDR.jl`:** la producción (`Transversal.jl`,
   `Longitudinal.jl`) no los usa; solo `test/runtests.jl`. Opciones: cablearlos
   (comparar numéricamente antes/después) o marcarlos como referencia/eliminar.
   Afecta a lo que el informe puede afirmar sobre la estadística validada.
5. 🟡 Recorrer capítulos 04 → 16 con el bucle §2.
6. ⚪ Tests de integración del pipeline completo; reportar el fallo Julia a graphify.

### Prompt de arranque para Cursor

```
Lee AGENTS.md y docs/PLAN_continuacion_2026-10.md. Vamos con el punto <N> de la
cola §4. Antes de leer código usa `graphify query/explain`. Presenta hallazgos en
tabla (🔴🟡🟢) sin modificar nada hasta que yo decida.
```
