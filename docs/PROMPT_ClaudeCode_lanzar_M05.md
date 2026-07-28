# Prompt para Claude Code — Validar M05 tras la unificación

## Contexto
Durante esta semana se ha reestructurado el proyecto: configuración unificada en
`config/pipeline.toml`, salida **única** en el árbol BIDS (se eliminó el doble
árbol `results/{ID}/{SES}/` con sufijos `_EC`), montaje de 31 canales
(`exclude_fp2 = false`, `n_channels_used = 31`), wPLI clásico (`use_dwpli = false`)
y carga diferida del dashboard (Genie). Nada de esto se ha ejecutado todavía.

## Objetivo
Lanzar el pipeline del sujeto de referencia **M05** y **verificar** que el
comportamiento es el esperado. **No refactorices ni "mejores" nada todavía** —
solo ejecutar, comprobar y reportar. Si algo falla, repórtalo pero no lo
arregles sin confirmación.

## Paso 0 — Preparación (antes de lanzar)
1. Confirma que el módulo carga sin errores y **sin compilar Genie**:
   ```bash
   julia --project=. -e 'using NeuroMIND; println("OK: módulo cargado")'
   ```
   Debe imprimir OK. Comprueba en la salida que NO aparece compilación de Genie
   (la carga diferida debe evitarlo).
2. Haz una copia de seguridad del resultado actual de M05 (para poder comparar
   la corrida nueva con la del 2026-07-09):
   ```bash
   cp -r results/subjects/sub-M05/ses-T2/eyesclosed /tmp/M05_prev_2026-07-09
   ```
3. Anota el estado de `results/` antes de lanzar:
   ```bash
   ls -1 results/                                   # no debe existir results/M05
   ls results/M05 2>/dev/null && echo "OJO: existe results/M05 antes de lanzar"
   ```

## Paso 1 — Lanzar M05
```bash
julia --project=. scripts/run_single_subject.jl --config config/pipeline.toml --force
```
(El sujeto por defecto de `config/pipeline.toml` es M05 / ses-T2 / eyesclosed.
`--force` fuerza el re-cómputo aunque exista caché.)

Durante la ejecución, observa la salida por terminal de los 8 pasos y anota
tiempos y cualquier `@warn`/`⚠`/`✗`.

## Paso 2 — Verificaciones (checklist)

**A. Salida ÚNICA en BIDS — no debe reaparecer el árbol heredado:**
```bash
test -d results/M05 && echo "❌ FALLO: reapareció results/M05/ (árbol heredado)" \
                    || echo "✅ OK: no hay results/M05/"
ls -1 results/subjects/sub-M05/ses-T2/eyesclosed/ | head
```

**B. 31 canales / 465 aristas en las matrices wPLI:**
```bash
# nº de columnas de la matriz (1 col 'channel' + 31 canales = 32)
head -1 results/subjects/sub-M05/ses-T2/eyesclosed/wpli_ALPHA.csv | awk -F, '{print NF" columnas (esperado 32)"}'
# nº de aristas en connectivity_edges.csv por banda (esperado 465)
grep -c ',ALPHA,' results/subjects/sub-M05/ses-T2/eyesclosed/connectivity_edges.csv
```

**C. El log del paso 8/8 NO debe usar nombres con sufijo `_EC`:**
```bash
grep -iE 'qc_channels_EC|_EC\.csv|_EC\.png' results/subjects/sub-M05/ses-T2/eyesclosed/pipeline_log.txt \
  && echo "❌ quedan _log con sufijo _EC (revisar)" \
  || echo "✅ OK: log sin sufijos _EC"
```

**D. Caché ICA dentro del árbol BIDS:**
```bash
ls results/subjects/sub-M05/ses-T2/eyesclosed/cache/ica_result.jls \
  && echo "✅ OK: caché en el árbol BIDS"
```

**E. `config_snapshot.toml` coincide con la config vigente (decisiones clave):**
```bash
grep -E 'exclude_fp2|use_dwpli|n_channels_used|method *=|wpli_method' \
  results/subjects/sub-M05/ses-T2/eyesclosed/config_snapshot.toml
# esperado: exclude_fp2=false, use_dwpli=false, n_channels_used=31,
#           method="first_window_mean" (baseline), wpli_method="hilbert"
```

**F. Índices globales actualizados:**
```bash
grep '^M05,' results/subjects_index.csv
grep '^M05,' results/qc/qc_decision_table.csv
```

**G. (Reproducibilidad) Comparar con la corrida previa:**
```bash
# Las tablas numéricas deberían coincidir (misma config). Diferencias esperadas
# solo en timestamps y en spectral_summary (antes reportaba epoch_length 2.0,
# ahora 1.0). Señala cualquier OTRA diferencia numérica.
diff <(sort /tmp/M05_prev_2026-07-09/overview.csv) \
     <(sort results/subjects/sub-M05/ses-T2/eyesclosed/overview.csv)
```

## Paso 3 — Reportar
Resume en un mensaje:
- ✅/❌ de cada check A–G.
- Tiempos por paso y total; si ICA se recomputó (debería, con `--force`).
- Cualquier `@warn`/error durante la corrida (esperado al menos el aviso de
  DELTA: 0.5 ciclos/época < 4.0).
- Nº de ficheros generados (CSV / PNG / JSON) y si falta alguno respecto a la
  corrida previa.
- **No modifiques código.** Si un check falla, descríbelo y espera instrucciones.

## Después
Con M05 validado, el siguiente paso (en otra sesión) será montar la capa de
presentación en terminal (Console/progreso/resumen, TDD fases 2-3), usando
`Crayons`/`ProgressMeter`/`PrettyTables` que ya están en el `Manifest.toml`.
