# Auditoría experta — Phase C (batch completo) + análisis transversal/longitudinal + M05 post-unificación

**Fecha de auditoría:** 2026-07-25
**Modo:** solo lectura sobre código y resultados existentes (los únicos cambios aplicados
son los 4 fixes de bajo riesgo listados en §9, autorizados explícitamente)
**Auditor:** agente técnico (Claude Code), con verificación cruzada de 3 sub-agentes de
exploración independientes + verificación directa de código para los hallazgos críticos
**Alcance:** corrida individual `sub-M05/ses-T2/eyesclosed` tras la unificación de salida
BIDS, batch completo Phase C (201/206 grabaciones), y los análisis de grupo transversal
(MS vs. control) y longitudinal (T1 vs. T2) en sus 4 combinaciones banda×condición
**Continúa de:** `docs/audits/run_consistency_audit.md` (2026-07-14), que auditó una única
corrida de M05 previa a la unificación BIDS. Los hallazgos de esa auditoría que siguen
vigentes o se resolvieron se referencian en §6.

---

## Veredictos por bloque

| Bloque | Veredicto | Motivo principal |
|--------|-----------|-------------------|
| M05 post-unificación (checks A–F) | **ÁMBAR** | Árbol BIDS único y config correctos; pero el montaje real cae a 30/31 canales por exclusión QC dinámica no documentada en AGENTS.md |
| Batch Phase C | **VERDE** | 201/5/0 confirmado exacto; los 5 SKIP tienen causa identificada (4 duplicados de ejecución, 1 duplicado de dato sin curar) |
| Integridad de datos de cohorte | **ÁMBAR** | Duplicado M16 T1 EC sin curar; conteos de AGENTS.md (78 sujetos) no coinciden con el roster real (77) ni con los datos procesados (75) |
| Estadística de grupo — implementación | **ROJO** | El código que genera todos los resultados de `results/transversal/` y `results/longitudinal/` no está cubierto por ningún test; el módulo que sí tiene tests (`GroupStats.jl`) no se usa en producción |
| Estadística de grupo — corrección matemática | **ÁMBAR** | Los tests de decisión (Mann-Whitney, Wilcoxon) están bien implementados; las columnas "t-test" de referencia usan una aproximación incorrecta (Z en vez de t) |
| Emparejamiento longitudinal | **ÁMBAR** | Bug de orden de filtro que puede excluir sujetos válidos del análisis EC-only o EO-only |
| Resultados científicos | **VERDE (con matices obligatorios)** | Cifras internamente consistentes y trazables a sus CSV/JSON fuente; requieren las salvedades de §4 antes de citarse |
| Reproducibilidad/trazabilidad | **ÁMBAR** | Sin regresión de julio 14, pero manifiesto de ejecución (`git_commit`/`config_hash`) sigue sin implementarse |
| Sincronización con `Report_Pre/` | **ÁMBAR** | El propio informe se auto-impuso una condición de citabilidad que Phase C probablemente satisface, pero el texto no se actualizó |

**Veredicto global: ÁMBAR**, con un hallazgo que merece ROJO aislado (estadística de
grupo sin cobertura de test) por su relevancia para la validez del resultado científico
central de este ciclo de trabajo.

---

## 0. Precondición de citabilidad — ¿puede citarse ya este batch?

`Report_Pre/chapters/10_conclusiones.tex:151–161` (el informe/tesis que consume estos
resultados, no versionado en git) contiene una condición explícita escrita por el propio
autor:

> "Los 78 sujetos ya procesados no reflejan la configuración actual. Se generaron antes
> del cierre del circuito de rechazo ICA (2026-07-06) y antes de que Fp2 pasara a
> incluirse por defecto (2026-07-08) [...]. **Ningún resultado de grupo (transversal o
> longitudinal) puede citarse todavía sin reprocesar el dataset completo con la
> configuración vigente.**"

Y en "Próximos pasos" (línea 196–199): *"Lanzar el batch completo [...] con la
configuración vigente (ICA con rechazo real, Fp2 incluido, 31 canales)"*.

**Evaluación:** Phase C se ejecutó el 2026-07-25, posterior a ambas fechas (07-06 y
07-08), con `exclude_fp2=false` confirmado en `config_snapshot.toml` de M05. Por fechas,
**es probable que esta condición ya esté satisfecha** — este batch parece ser exactamente
el reprocesamiento que el propio informe pedía. Pero:

1. No verifiqué el diff de código exacto del "cierre del circuito de rechazo ICA" del
   07-06 (la auditoría de julio 14 sí documentó un desfase entre el log y
   `load_ica_labels()` en ese momento; no repetí esa comprobación específica ahora).
2. El hallazgo §2.3 (exclusión dinámica de Fp2 por QC pese a `exclude_fp2=false`) matiza
   "Fp2 incluido" de forma que el propio capítulo de conclusiones no contempla.
3. `Report_Pre/` no se tocó en esta auditoría (fuera del alcance acordado) — el capítulo
   sigue afirmando que el resultado de grupo no es citable.

**Recomendación:** antes de citar cualquier cifra de §4 en un documento académico,
confirmar explícitamente que Phase C cumple ambas condiciones y actualizar
`10_conclusiones.tex` §Limitaciones y §Próximos pasos en consecuencia (no se hizo aquí).

---

## 1. Resumen ejecutivo

El ciclo de trabajo auditado (M05 individual → batch de 206 grabaciones → transversal +
longitudinal) es **coherente y en gran parte correcto**: los conteos del batch cuadran
exactamente con lo documentado (201 OK/5 SKIP/0 ERR), los resultados de grupo son
trazables a sus fuentes, y todos los problemas de trazabilidad de la auditoría anterior
(julio 14) están resueltos.

Hay, sin embargo, tres problemas que sí importan para la validez del resultado científico:

1. **El código que calcula las 24 tablas de significancia de `results/transversal/` y
   `results/longitudinal/` no tiene ningún test.** Existe un módulo (`GroupStats.jl`) con
   tests reales, pero los scripts de producción no lo usan — reimplementan la estadística
   por su cuenta. Los "241+ tests passing" que cita `AGENTS.md` no cubren este código.
2. **Un bug de orden de filtro en el emparejamiento longitudinal** puede excluir sujetos
   válidos incluso del análisis EC-only, contradiciendo el diseño documentado del propio
   script.
3. **La exclusión de canales varía por sujeto** (QC dinámico por z-score, siempre activo,
   independiente de `exclude_fp2`), lo que explica por qué el análisis transversal EC
   termina usando solo 24 de 31 canales tras intersectar toda la cohorte (~23% de aristas
   perdidas) — mecanismo real y documentado en el código, pero no reflejado en `AGENTS.md`.

El resultado central —21 conexiones significativas en transversal EC (dominadas por
ALPHA), prácticamente nada en EO, y **cero** en ambos longitudinales— es internamente
consistente y probablemente refleja un patrón real, pero el n longitudinal efectivo
(15–18 pares, no los 27 "completos" que cita `AGENTS.md`) exige matizar cualquier
conclusión de "no hay progresión detectable" como "no hay progresión detectable con este
poder estadístico".

---

## 2. Corrida individual — `sub-M05/ses-T2/eyesclosed` post-unificación

### 2.1 Checklist de validación (`docs/PROMPT_ClaudeCode_lanzar_M05.md`, checks A–F)

| Check | Resultado | Evidencia |
|-------|-----------|-----------|
| A — sin árbol legacy `results/M05/` | ✅ PASA | `find results -iname "*M05*"` solo devuelve `results/subjects/sub-M05` |
| B — 31 canales / 465 aristas | ❌ **FALLA** | `wpli_ALPHA.csv` tiene 31 columnas (`channel` + 30 canales, Fp2 ausente); `connectivity_edges.csv` tiene 435 filas/banda, no 465 |
| C — logs sin sufijo `_EC` | ✅ PASA | `grep "_EC" pipeline_log.txt` → 0 coincidencias |
| D — caché ICA dentro del árbol BIDS | ✅ PASA | `cache/ica_result.jls` (11.9 MB) vive en `results/subjects/sub-M05/ses-T2/eyesclosed/cache/` |
| E — `config_snapshot.toml` coherente | ✅ PASA (literal) | `exclude_fp2=false`, `use_dwpli=false`, `n_channels_used=31`, `method="first_window_mean"`, `wpli_method="hilbert"` — todos correctos como *declaración*; ver 2.3 para el efecto práctico |
| F — índices globales actualizados | ✅ PASA | `subjects_index.csv` y `qc/qc_decision_table.csv` tienen fila M05 con timestamp `2026-07-25T10:43:01–02`, coincidente con el fin del pipeline |

### 2.2 Cronología de la corrida

Ejecución única `10:20:58 → 10:43:02` (1324.0 s totales). Fases 1–7 (carga → wPLI): ~23 s.
Surrogates (`10:21:21 → 10:43:01`, ~21.6 min): ALPHA 135.6 s, BETA_HIGH 136.2 s,
BETA_LOW 131.3 s, BETA_MID 132.5 s, **DELTA 495.5 s (~4× el resto)**, GAMMA 131.6 s,
THETA 131.9 s. La banda DELTA domina el tiempo total sin explicación evidente en el log;
no confirmado como sistémico (no se comprobó contra los otros 200 sujetos — ver §7.2).

`grep -inE "warn|⚠|✗"` sobre `pipeline_log.txt` → 0 coincidencias persistidas. El aviso de
DELTA (0.5 ciclos/época < `min_cycles_for_wpli=4.0`) sí está implementado
(`src/connectivity/wPLI.jl:145-152`, `@warn`) y documentado (`AGENTS.md §7`), pero
`@warn` va a stderr mientras `pipeline_log.txt` se construye por impresión manual —
probablemente se emitió en consola sin quedar persistido en el log del sujeto.

### 2.3 Mecanismo real de exclusión de canales (hallazgo central de este bloque)

El check B "falla" no por un bug de código sino por un **mecanismo intencional y
documentado en el propio código, pero no reflejado en la documentación de más alto
nivel**. Mecanismo exacto (`src/SingleSubjectPipeline.jl:975-990`):

```julia
# A la lista efectiva se unen siempre bad_ch del QC (z-score).
exclude_fp2 = Bool(get(montage_cfg, "exclude_fp2", true))
...
all_excl = unique(vcat(montage_excl, bad_ch))   # bad_ch SIEMPRE se une, sin importar exclude_fp2
```

`config/pipeline.toml:717,725-726` incluso deja constancia de que el campo
`n_channels_analysis` es **"INERTE (documental) — informativo, no se lee"**: no gobierna
nada, es solo una nota.

En M05, Fp2 tiene z-score = 4.63 (por encima del umbral 3.0 de QC — ver
`AGENTS.md §7`), por lo que cae en `bad_ch` y se excluye del análisis **aunque
`exclude_fp2=false` diga que se debe conservar**. El propio log resultante es engañoso a
primera lectura:

```
[10:21:15]   Montaje análisis: 30 ch · Fp2=incluido · excluidos=Fp2
```

Ambas afirmaciones son ciertas bajo definiciones distintas ("incluido" = no excluido por
la política estática; "excluidos=Fp2" = excluido por QC dinámico), pero leídas de corrido
parecen contradecirse. **Consecuencia documental:** `AGENTS.md §7` afirma sin matizar
que el montaje vigente es "31×31 con 465 edges per band" — cierto solo si ningún canal
resulta sospechoso por QC, lo cual no se cumple en M05 y probablemente tampoco en una
fracción no despreciable de la cohorte (ver cascada en §3.3). Corregido en `AGENTS.md`
como parte de los fixes de esta sesión (§9.1).

### 2.4 `figures/aux/` — los 12 visores interactivos de verificación

`results/subjects/sub-M05/ses-T2/eyesclosed/figures/aux/` contiene 12 scripts Julia
standalone (`plot_raw.jl`, `plot_raw_butterfly.jl`, `plot_raw_histogram.jl`,
`plot_filtered_vs_raw.jl`, `plot_ica_components.jl`, `plot_ica_before_after.jl`,
`plot_epochs.jl`, `plot_baseline.jl`, `plot_spectral.jl`, `plot_spectral_PSD.jl`,
`plot_connectivity.jl`, `plot_surrogate.jl`), fechados 23–25 de julio, más ~22 PNG
estáticos generados por ellos como snapshots de verificación puntual.

Inspección de `plot_raw.jl` (representativo del patrón de los 12): cada uno monta su
propio servidor HTTP mínimo sobre `Sockets` puro (sin `HTTP.jl`/Genie —
`listen`/`accept`/`readline`/`write` manuales), sirve una UI HTML/JS/canvas autocontenida,
y calcula su ruta de datos de entrada como `normpath(joinpath(@__DIR__, "..", "..",
"tables", "raw_signal.csv"))` — es decir, relativa a su propia posición dos niveles bajo
la carpeta de resultados del sujeto. Puertos dedicados 8765–8775, uno por aspecto,
hardcodeados a `sub-M05/ses-T2/EC`.

El patrón es **arquitectónicamente idéntico** al de los visores de cohorte oficiales
(`src/transversal/plot_transversal.jl` :8781, `src/longitudinal/plot_longitudinal.jl`
:8780): mismo estilo de servidor HTTP a mano, misma factura de UI. Diferencia clave: los
visores de cohorte están versionados y documentados en `AGENTS.md`; estos 12 **viven
dentro de `results/`**, que `.gitignore:8` excluye por completo — no están en el
repositorio, no se referencian desde `src/`/`scripts/`/`NeuroMIND.jl`, y se perderían si
se borra `results/` o se cambia de máquina. Si en el futuro se implementa el
"clean-slate por unidad hoja" que propone `docs/TDD_ejecucion_rutinas.md` (borrar y
regenerar la carpeta del sujeto antes de cada corrida), estos 12 scripts se borrarían sin
que git los proteja. Versionados en esta sesión — ver §9.3.

---

## 3. Batch completo — Phase C (2026-07-25)

### 3.1 Conteo real vs. declarado

El log `results/logs/batch_run_2026-07-25_11-26.csv` (11:26:02.858 → 12:17:20.707,
**51 min 18 s**) es la corrida de Phase C real. Existen además 5 logs diminutos de
ensayos manuales previos de un solo sujeto (`batch_run_2026-07-25_{11-16,11-17,11-20,
11-23,11-24}.csv`; uno de ellos registra un `ERROR` de una prueba aislada —
`MethodError(open; (Base.DevNull(); "w"))`— sin relación con Phase C).

**Conteo confirmado: 201 OK / 5 SKIPPED / 0 ERROR** — coincide exactamente con lo
documentado en `AGENTS.md`.

### 3.2 Los 5 SKIP, con causa verificada

| Sujeto/sesión/cond. | Timestamp SKIP | Causa real |
|---|---|---|
| M10 T1 EC | 11:26:02.858 | Ya procesado en ensayo previo (OK 11:17:32) |
| M10 T1 EO | 11:26:02.869 | Ya procesado en ensayo previo (OK 11:20:40) |
| M11 T1 EC | 11:26:02.871 | Ya procesado en ensayo previo (OK 11:24:27) |
| **M16 T1 EC** | 11:30:14.912 | **Duplicado dentro de la misma corrida** — 0.08 s antes, la misma clave ya se había marcado OK |
| M05 T2 EC | 11:52:38.671 | Ya procesado antes del batch (10:43:02, la corrida individual de §2) |

4 de los 5 SKIP son benignos (reprocesamiento evitado de sujetos ya corridos ese mismo
día). El quinto (M16) tiene una causa raíz distinta y real:

`data/full_data/inventory.csv` tiene **dos archivos `.vhdr` para la misma clave**
(M16, T1, EC): `ojoscerrados.vhdr` y `ojoscerradpsb.vhdr` (línea 24 y 25), **ambos con
`excluded=false`** y sin nota de curación. `load_jobs()` en
`scripts/run_batch_pipeline.jl:140-146` no deduplica por clave `(subject, session,
condition)`, solo respeta el flag `excluded` — así que generó dos jobs para la misma
combinación, y el segundo se auto-descartó como "ya procesado" en vez de fallar o
alertar.

Contraste: `inventory.csv` tiene un caso análogo *sí* curado correctamente — M37 T2 con
dos archivos por condición, donde el segundo está marcado
`excluded=true, note="DUPLICADO T2: usar T2_1"`. M16 es el único duplicado sin ese flag.
**No se puede determinar desde el código ni los nombres de archivo cuál de los dos
`.vhdr` de M16 es el válido** — se deja como pregunta abierta para el usuario (§10).

### 3.3 Conteos de cohorte — discrepancias con `AGENTS.md`

| Métrica | `AGENTS.md` | Real en disco |
|---|---|---|
| Sujetos totales | 78 (41 MS + 37 control) | **77** (41 MS + **36** control — faltan MC5/6/11/12 del rango declarado) |
| Sujetos con datos en `results/subjects/` | 78 (implícito) | **75** (39 MS + 36 control; M4 y M6 son grabaciones `ODDBALL`, correctamente excluidas por no ser resting-state) |
| Grabaciones válidas (`excluded=false`) | 206 | 206 ✓ en `inventory.csv`, pero **205** filas únicas en `qc_decision_table.csv`/`subjects_index.csv` (la diferencia es el duplicado M16 de §3.2) |
| Desglose QC (`qc_decision_table.csv`, `final_decision`) | — | `include`=126, `include_with_warning`=53, `exclude`=19, `manual_review`=7 (total 205) |

No es una regresión de esta sesión — son conteos que probablemente ya estaban desfasados
antes del batch y que Phase C simplemente hizo visibles al ser la primera pasada completa
con la configuración vigente. Corregible con una actualización de texto en `AGENTS.md`
(fuera del alcance de los fixes de bajo riesgo de esta sesión — ver §9).

---

## 4. Análisis de grupo — resultados reales

### 4.1 Frescura de los resultados

Todos los CSV/JSON centrales de `results/transversal/` y `results/longitudinal/` están
fechados **después** de Phase C (fin 12:17:20): CSV de aristas 12:19–12:20,
`*_summary.json` 12:42:19, figuras EC 12:41. Reflejan el batch más reciente, no una
corrida anterior. Única excepción menor: 2 CSV auxiliares de "viewer"
(`viewer_edges_top_dz_ALPHA.csv` 12:10:52, `viewer_edges_top_d_ALPHA.csv` 12:18:55) rozan
el cierre de Phase C, pero son artefactos secundarios de exploración, no resúmenes
centrales.

### 4.2 Cifras — tabla única de referencia

| Análisis | n incluidos | n_total_sig | best_band | Desglose bandas con q<0.05 |
|---|---|---|---|---|
| Transversal EC | MS=32, control=36 | **21** | ALPHA | ALPHA=12, THETA=8, BETA_MID=1, resto=0 |
| Transversal EO | MS=27, control=36 | **3** | BETA_LOW | BETA_LOW=3, resto=0 |
| Longitudinal EC | pares=15/27 candidatos | **0** | "" (vacío) | 0 en las 7 bandas |
| Longitudinal EO | pares=18/27 candidatos | **0** | "" (vacío) | 0 en las 7 bandas |

`paired_subjects.csv` (EC y EO) confirma 27 candidatos MS con T1+T2 presentes —
coincide con lo declarado en `AGENTS.md`, pero solo como conteo de candidatos crudos.
Tras aplicar QC por sesión, quedan `included=true` **15/27 en EC** y **18/27 en EO**
(12 y 9 excluidos respectivamente, todos por `QC T1/T2=exclude` o `manual_review`).

### 4.3 Lectura científica — con las salvedades obligatorias

- **Transversal EC (21 sig., dominado por ALPHA)** es el resultado más fuerte y el único
  que se alinea con parte de la hipótesis oficial del proyecto (alteraciones en banda
  alfa). El hallazgo secundario en THETA (8 conexiones) **no** estaba anticipado por la
  hipótesis declarada (que menciona solo alfa/beta) — mencionarlo como hallazgo
  exploratorio, no confirmatorio.
- **Transversal EO (solo 3 sig., en BETA_LOW)** es sustancialmente más débil que EC —
  consistente con que la conectividad en reposo con ojos cerrados suele ser más estable y
  con mayor SNR que con ojos abiertos, pero la caída es grande y merece mención explícita
  como limitación, no omitirse.
- **Longitudinal EC y EO — cero conexiones significativas en ambos.** Esto **no debe
  presentarse como "no hay progresión de conectividad T1→T2"** sin matizar: el n efectivo
  tras QC (15–18 pares) es muy inferior a los "27 pares completos" que cita `AGENTS.md`,
  y con ese tamaño muestral la ausencia de significancia tras FDR-BH es compatible tanto
  con "no hay cambio real" como con "no hay poder estadístico suficiente para detectarlo".
  La afirmación correcta es: *no se detectó ningún cambio significativo con el tamaño
  muestral disponible tras QC*.
- El hard-intersect de canales (§2.3 → cascada de grupo: transversal EC usa 24/31
  canales, 276/465 aristas posibles = 59%; longitudinal EC usa 27/31) reduce la potencia
  estadística real por debajo de lo que "31 canales" sugiere, sin que ningún `@warn`
  activo lo señale durante la ejecución — solo queda registrado pasivamente en
  `band_statistics.csv`, la consola, y la UI de los visores (`n_channels`).

---

## 5. Código de los análisis de grupo — revisión experta

### 5.1 Metodología estadística (verificada matemáticamente)

- **Transversal**: Mann-Whitney U (`mannwhitney_p`,
  `scripts/run_transversal_analysis.jl:166-195`) es el test de decisión, con corrección
  de empates en la varianza. `welch_t` (líneas 151-163) es solo columna de referencia.
- **Longitudinal**: Wilcoxon signed-rank pareado (`wilcoxon_p`,
  `scripts/run_longitudinal_analysis.jl:160-176`) es el test de decisión — diseño
  correcto, pareo real T1/T2 por sujeto, no tratado como muestras independientes.
  `paired_t` es solo referencia.
- **Elección de test apropiada** para ambos diseños (no paramétrico, apto para
  distribuciones no gaussianas de wPLI y para el tamaño muestral disponible).

**Imprecisión confirmada (impacto bajo):** `welch_t` y `paired_t` convierten su
estadístico a p-valor vía aproximación normal estándar (`_norm_cdf` —
`run_transversal_analysis.jl:158-159`, `run_longitudinal_analysis.jl:186-188`), no con
distribución t de Student ni grados de libertad de Welch-Satterthwaite. Son, en la
práctica, aproximaciones Z etiquetadas como "t-test". No afecta la decisión real (que usa
Mann-Whitney/Wilcoxon), pero el nombre de columna puede inducir a error a quien las use
de forma aislada.

**Asimetría de rigor:** Mann-Whitney corrige empates en σ²
(`run_transversal_analysis.jl:176-188`); Wilcoxon (longitudinal) no tiene corrección
análoga en σW (`run_longitudinal_analysis.jl:165-174`). Ninguna de las dos aplica
corrección de continuidad (±0.5).

### 5.2 Corrección de comparaciones múltiples (FDR-BH)

Aplicado **por banda separada** sobre el conjunto de aristas de esa banda
(`bh_qvalues(p_mw)` dentro de `for band in BANDS`, `run_transversal_analysis.jl:605`;
análogo en `run_longitudinal_analysis.jl:684`), no pooled a través de bandas×aristas.
Existe una segunda familia FDR más pequeña (7 tests, uno por banda) para
`global_mean_wpli_statistics.csv` (línea 724).

El algoritmo BH en sí (`bh_qvalues`, líneas 116-124) está **correctamente implementado**:
orden ascendente de p-valores, paso monótono desde el p mayor hacia abajo, cap a 1.0 —
verificado línea por línea. La frase "full cohort" del changelog de `AGENTS.md` es
ambigua y puede leerse como "una sola familia global de corrección", cuando en realidad
son 7 familias independientes (una por banda); no es incorrecto, pero merece
clarificación.

### 5.3 Política de intersección de canales ("hard intersect")

`common_set = reduce(intersect, all_ch_sets)` (`run_transversal_analysis.jl:578-581` y
`:734-738`; `run_longitudinal_analysis.jl:650-653`, `:778-782`): no se descarta el
sujeto completo si le falta un canal, se descarta **ese canal para toda la cohorte** de
esa banda/condición. Cuantificado con datos reales (§4.3): 24/31 canales en transversal
EC, 27/31 en longitudinal EC — y la variabilidad es aún mayor de lo que sugiere solo EC:
al regenerar las figuras EO (§9.4) se confirmó **21/31 canales en transversal EO** y
**28/31 en longitudinal EO**, es decir la condición EO pierde casi un tercio del montaje
en el análisis transversal. No hay `@warn` activo si `n_channels` cae por debajo de 31;
la información solo queda disponible pasivamente (CSV, consola, UI del visor). Esto ya
está identificado como limitación conocida en `AGENTS.md §10` ("channel-intersection
policy... currently hard intersect") — esta auditoría aporta la cuantificación real que
faltaba.

### 5.4 Bug de emparejamiento longitudinal

`scripts/audit_full_dataset.jl:411` define:

```julia
include_longitudinal = t1_ec && t1_eo && t2_ec && t2_eo   # exige las 4 condiciones
```

`run_longitudinal_analysis.jl:426-429` filtra por ese flag **antes** de aplicar la
lógica fina por-condición (líneas 432-441), que fue diseñada explícitamente para aceptar
un sujeto con EC completo aunque le falte EO (o viceversa) — así lo indica el propio
header del script ("EC y EO en paralelo", línea 9). Efecto: un sujeto con par EC completo
pero EO incompleto queda excluido **incluso del análisis EC-only**, contradiciendo el
diseño documentado. El fix de *formato* del booleano (changelog 2026-05-26, `_as_bool`)
es correcto — el problema actual es de *orden/semántica* del filtro combinado, no de
parseo.

No se cuantificó cuántos de los 12 (EC) / 9 (EO) pares excluidos tras QC (§4.2) se deben
a este bug específico frente a exclusión QC genuina — requeriría instrumentar el script,
fuera del alcance de esta sesión de diagnóstico.

### 5.5 Dos implementaciones estadísticas paralelas — la testeada no es la de producción

Confirmado directamente (no solo por inferencia): `src/statistics/GroupStats.jl` +
`FDR.jl` tienen tests reales (`test/runtests.jl:708` `@testset "GroupStats"`, línea
727 testea `wilcoxon_signed_rank`) y están incluidos en `NeuroMIND.jl:114-115`.

`grep -n "^using\|^import" scripts/run_transversal_analysis.jl
scripts/run_longitudinal_analysis.jl` confirma que ninguno de los dos scripts importa ese
módulo:

```
scripts/run_transversal_analysis.jl:46:  using CSV, DataFrames, Statistics, LinearAlgebra, Dates, TOML, Printf, CairoMakie
scripts/run_transversal_analysis.jl:49:  using .GroupVizCommon
scripts/run_longitudinal_analysis.jl:50: using CSV, DataFrames, Statistics, LinearAlgebra, Dates, TOML, Printf, CairoMakie
scripts/run_longitudinal_analysis.jl:53: using .GroupVizCommon
```

Ambos scripts reimplementan localmente `mannwhitney_p`, `wilcoxon_p`, `bh_qvalues`,
`_erf_approx` — código sin ningún test, y es el que efectivamente genera **todas** las
tablas de `results/transversal/` y `results/longitudinal/`. Ningún archivo bajo `test/`
referencia `run_transversal_analysis.jl`, `run_longitudinal_analysis.jl`,
`mannwhitney_p` ni `wilcoxon_p` — la cobertura cero es total, no parcial. Los "241+ tests
passing" que `AGENTS.md` cita como logro **no cubren el código real de producción para el
análisis de grupo**. Las dos implementaciones ya han divergido: la corrección de empates
en Mann-Whitney existe en los scripts pero está ausente en `GroupStats.jl:13-40`.

Este es, junto con §5.4, el hallazgo con mayor impacto potencial sobre la validez
científica del resultado — no porque haya evidencia de que la implementación actual esté
mal (el algoritmo BH y el diseño de los tests de decisión se verificaron correctos línea
por línea, §5.1–5.2), sino porque **no hay red de seguridad automatizada** que lo
garantice hacia adelante.

### 5.6 Calidad de código — hallazgos menores

- `load_wpli` sin try/catch (`run_transversal_analysis.jl:235`,
  `run_longitudinal_analysis.jl:237`): una celda `missing`/mal formada en un solo CSV de
  un solo sujeto puede abortar todo el batch de análisis de grupo (EC+EO) — contrasta con
  el try/catch sí presente para la generación de figuras (`:664-678`, `:680-696`).
  Corregido en esta sesión (§9.2).
- ~300+ líneas duplicadas verbatim entre los dos scripts (`load_wpli`, `load_band_power`,
  `load_qc_table`/`qc_decision`/`qc_ok`, `realign_matrix`, `nodal_strength_degree`,
  `_erf_approx`/`bh_qvalues` — confirmado por checksum idéntico). `GroupVizCommon.jl`
  demuestra que el equipo sabe factorizar cuando lo considera prioritario; no se aplicó a
  la capa de carga de datos/estadística.
- `0.05` (umbral FDR) hardcodeado ~4 veces por script (`:625`,`:815` transversal;
  `:701`,`:854` longitudinal) en vez de leer `cfg.statistics.fdr_q` — que sí es un campo
  esperado por los tests (`test/runtests.jl:802`) y usado en `GroupStats.jl:132`, pero
  `config/pipeline.toml` **no tiene sección `[statistics]`** en absoluto.
- `top_edges_by_effect_$(band).csv` existe solo en transversal (`:636-639`); sin
  equivalente en longitudinal — deriva estructural entre los dos scripts.
- Redundancia menor: `paired_t` devuelve un `dz` que se descarta
  (`pp_vec[k], _ = paired_t(v1, v2)`, `run_longitudinal_analysis.jl:681`) y se recalcula
  acto seguido vía `cohen_dz` (línea 682).
- `N_PAIRED_DESIGN=30` (`run_longitudinal_analysis.jl:60`, comentario "Fig. 3.1 — pares EM
  T1+T2 comparables"): confirmado con grep que aparece en 7 sitios (líneas 60, 448, 473,
  518, 547, 594, 600, 899), todos dentro de `println`/interpolación de string o campos
  JSON de metadata — **nunca en aritmética real**. Es una tercera cifra (junto a 27
  candidatos y 15/18 incluidos) que solo aparece en texto/logs; cosmético, no afecta
  ningún cálculo, pero puede confundir a quien lea los JSON de salida.
- Cita bibliográfica incorrecta: `AGENTS.md`/changelog citan la aproximación de `erf`
  como "A&S 26.2.17"; la fórmula realmente implementada (coeficientes
  `0.254829592…`, `p=0.3275911`) es **A&S 7.1.26** — correctamente citada solo dentro de
  `GroupStats.jl:161` (el módulo no usado en producción). Corregido en esta sesión (§9.1).

---

## 6. Regresión vs. auditoría de julio 14 — todo resuelto

| Problema (auditoría 2026-07-14) | Estado ahora |
|---|---|
| Caché ICA desincronizada 46 días | **Mejorado**: ahora 3 días de diferencia (caché 22 jul, corrida 25 jul); reuso intencional y trazable (`"ICA cargado desde caché"` en el log) |
| 10 figuras `*_EC.png` obsoletas de mayo mezcladas | **Resuelto** — 0 coincidencias bajo `sub-M05` |
| Árbol legacy `results/M05/` | **Resuelto** — confirmado ausente |
| `channel_statistics_compare.csv` sin generador | **Resuelto** — el archivo ya no existe |
| `git_commit`/`config_hash` nulos en manifiesto | **No resuelto** — sigue sin implementarse; consistente con que `docs/TDD_ejecucion_rutinas.md` (22 julio) es una propuesta de diseño aún no implementada |

---

## 7. Buenas prácticas confirmadas

Para que esta auditoría no quede solo como lista de problemas:

- Regla `erf()` respetada en todo el repo: `grep -rn "erf("` en `scripts/` + `src/` = 0
  resultados; solo existe `_erf_approx(`, sin dependencia de `SpecialFunctions`.
- `honest_best_band` (`GroupVizCommon.jl:166-174`) verificado end-to-end contra datos
  reales: transversal EC (`n_total_sig=21`) → `"ALPHA"` (correcto, 12 sig. vs 8 THETA);
  longitudinal EC (`n_total_sig=0`) → `""`. Funciona exactamente como promete el
  changelog — no hay riesgo de reportar una "banda ganadora" cuando no hubo hallazgos.
- Separación limpia entre cálculo y visualización: los visores de cohorte
  (`plot_transversal.jl`, `plot_longitudinal.jl`) nunca recalculan estadística, solo leen
  columnas `p_mannwhitney`/`q_value` ya calculadas, y sí muestran `n_channels` —
  cumpliendo la mitigación de transparencia que promete `AGENTS.md §10` para el hard
  intersect.
- Fallbacks defensivos: `groups.csv` ausente → `@warn` + `exit(1)` explícito (líneas
  388-395); `longitudinal_pairs.csv` ausente → auto-detección por escaneo de directorio
  (líneas 451-474) en vez de fallo silencioso.
- Exclusiones auditables por sujeto con motivo explícito, confirmado en los CSV reales
  (`results/transversal/EC/subject_inclusion.csv`, `paired_subjects.csv`).
- El algoritmo BH y los tests de decisión (Mann-Whitney, Wilcoxon pareado) están
  matemáticamente bien implementados y son la elección apropiada para el diseño (§5.1–5.2).
- Todos los problemas de trazabilidad de la auditoría previa quedaron resueltos (§6).

---

## 8. Deuda de sincronización con `Report_Pre/`

No se editó ningún archivo de `Report_Pre/` en esta sesión (fuera de alcance acordado),
pero quedan señalados para una revisión posterior:

- `chapters/10_conclusiones.tex` §Limitaciones (líneas 151-161) y §Próximos pasos
  (líneas 196-199): la condición de citabilidad que impone probablemente ya está
  satisfecha por Phase C — ver §0. Requiere confirmación explícita y actualización de
  texto.
- Cualquier cifra de "78 sujetos" citada en `Report_Pre/` debería revisarse contra los
  conteos reales de §3.3 (77 en el roster, 75 con datos procesados).
- El caso de trabajo verificado de `Report_Pre/` sigue siendo únicamente
  `sub-M05/ses-T2/eyesclosed` — el hallazgo de exclusión dinámica de Fp2 (§2.3) afecta
  directamente a los capítulos que citan "31 canales/465 aristas" para ese sujeto
  (`03_raw`, `05_ica`, `09_espectral`, `10_conectividad`, `11_surrogates`, según el
  propio `AGENTS.md §13` referenciado desde `Report_Pre`).

---

## 9. Fixes de bajo riesgo aplicados en esta sesión

Autorizados explícitamente por el usuario; ninguno altera resultados científicos ya
generados.

1. **`AGENTS.md`**: aclarado §7 (exclusión de canales = `exclude_fp2` estático + QC
   dinámico por z-score, que siempre se aplica; "31×31/465 aristas" matizado como el caso
   sin canales sospechosos) y corregida la cita "A&S 26.2.17" → "A&S 7.1.26" en el
   changelog.
2. **`scripts/run_transversal_analysis.jl`** y **`run_longitudinal_analysis.jl`**:
   `load_wpli` blindado con try/catch (mensaje con sujeto/banda/causa en vez de abortar
   todo el batch de análisis de grupo).
3. **`figures/aux/` versionado**: los 12 visores copiados a
   `scripts/aux_viewers/sub-M05_ses-T2_eyesclosed/` con rutas de datos corregidas para
   apuntar a `results/subjects/sub-M05/ses-T2/eyesclosed/` desde la nueva ubicación.
   Originales en `results/` intactos.
4. **Figuras EO regeneradas**: `heatmap_triplet_*`/`sig_network_*`/`explore_network_topN_*`
   completadas para EO en transversal y longitudinal (antes solo existían para EC).

---

## 10. Preguntas que requieren decisión humana

1. **M16 T1 EC**: ¿cuál de los dos `.vhdr` (`ojoscerrados.vhdr` o `ojoscerradpsb.vhdr`)
   es la grabación válida? No se puede inferir del código ni del nombre de archivo.
   Necesario para marcar `excluded=true` en el otro dentro de `data/full_data/inventory.csv`.
2. ¿Se confirma que Phase C satisface la condición de citabilidad de
   `Report_Pre/chapters/10_conclusiones.tex` (§0)? Si sí, ¿se actualiza ese capítulo en
   una sesión posterior?
3. ¿Se autoriza reconectar `run_transversal_analysis.jl`/`run_longitudinal_analysis.jl` a
   `src/statistics/GroupStats.jl` (§5.5)? Cambiaría potencialmente p/q-valores exactos
   (aunque el algoritmo verificado es equivalente) — requiere diff explícito contra
   `results/` actual antes de sobrescribir nada ya citado.
4. ¿Se autoriza corregir el orden de filtro del emparejamiento longitudinal (§5.4)?
   Cambiaría el n de pares incluidos — mismo requisito de diff explícito.
5. ¿Se añade una sección `[statistics]` a `config/pipeline.toml` y se conecta a los
   scripts (§5.6), o se deja el `0.05` hardcodeado documentado como decisión deliberada?

---

## 11. Plan de remediación por fases (referencia para sesiones futuras)

- **Fase 0 (decisión humana, sin código)**: preguntas de §10.
- **Fase 1 (aplicada hoy, §9)**: aclaraciones de documentación, blindaje `load_wpli`,
  versionado de `figures/aux/`, figuras EO.
- **Fase 2 (cambia cifras — requiere aprobación explícita antes de ejecutar)**:
  reconectar la estadística de grupo a `GroupStats`/`FDR` (§5.5) y corregir el
  emparejamiento longitudinal (§5.4). Ambos deben re-ejecutarse juntos, diferenciarse
  contra los resultados actuales, y mostrarse el diff antes de sobrescribir nada.
- **Fase 3 (mayor alcance, diferible)**: manifiesto de ejecución (`run_manifest.json`,
  `git_commit`, `config_hash` — ya diseñado en `docs/TDD_ejecucion_rutinas.md`),
  deduplicación de las ~300 líneas compartidas entre scripts de grupo, política explícita
  sobre variabilidad de canales por sujeto (documentar vs. forzar montaje uniforme
  descartando sujetos con canales sospechosos en vez de canales).

---

*Fin del informe. Los únicos cambios de código realizados durante esta auditoría son los
4 fixes de bajo riesgo listados en §9, aplicados solo tras autorización explícita.*
