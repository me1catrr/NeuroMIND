# CHANGELOG — NeuroMIND

> Historial detallado de cambios, movido desde `AGENTS.md` §13 el 2026-10-13 para que
> no se cargue en cada sesión de los agentes. Las entradas nuevas van arriba.
> `AGENTS.md` §13 apunta aquí (y `Report_Pre/chapters/12_espectral.tex` cita esa sección).

### 2026-10-15 — `build_clinical_master.jl --check`

- Nueva opción `--check`: tras escribir el libro y el CSV cruza 22 comprobaciones (BIDS en disco ↔
  inventario ↔ `raw_selection.csv` ↔ `participants.tsv` ↔ libro/CSV releídos ↔ base original de
  Alejandro ↔ cifras del capítulo 4 y de sus figuras), imprime las filas de los Cuadros 4.1/4.2 en
  formato LaTeX y sale con código 1 si algo falla. Probada con un fichero sobrante en el BIDS.
  Las cifras del informe se siguen corrigiendo a mano; no se genera LaTeX.

### 2026-10-15 — M07 recupera su EO de T1 y se cierran las fichas aplicadas

- **M07 T1 EO:** `config/raw_selection.csv` declara `M7_T1_*_Cog-10-1` como EO (acción `use`,
  confirmación de Alejandro del 2026-09-28, ficha ADQ-03). `audit_full_dataset.jl` admite ahora
  `use` sobre un archivo sin sufijo de ojos y lo reetiqueta en lugar de dejarlo como ODDBALL.
  Inventario: 215 filas, 210 utilizables, 5 excluidas. `build_bids_full.jl --subjects M07`
  añade el registro (señal idéntica al original) y el dataset pasa a 210 grabaciones
  (105 EC + 105 EO); pares longitudinales: 28 EC y 28 EO. `participants.tsv` no cambia.
- **Validación con `bids-validator-deno` 3.0.2:** 0 errores; los 4 avisos son
  `EVENTS_TSV_MISSING`, `TOO_FEW_AUTHORS`, `JSON_KEY_RECOMMENDED` y `SIDECAR_KEY_RECOMMENDED`
  (todos documentados en `data/bids/README.md`). Además, 204 `.eeg` completos idénticos byte a
  byte a su origen; los 6 recortes comprobados antes contra su `SourceWindow`.
- **Cruce completo de consistencia (15/10):** BIDS ↔ `inventory.csv` (210 registros, 0 sobrantes,
  0 faltantes, 105 sesiones) ↔ `raw_selection.csv` (15 filas, 12 archivos) ↔ `participants.tsv`
  (`has_t1`/`has_t2` frente a carpetas `ses-*`) ↔ libro y CSV unificados (308 celdas `eeg_*`, hoja
  Participantes idéntica al CSV, 215 filas en Grabaciones EEG) ↔ base de Alejandro (462 valores
  EDSS/D-FIS/BDI-II sin diferencias) ↔ Cuadros 4.1 y 4.2 del informe.
- **Fichas trasladadas** (cuaderno e informe): MONT-07, ADQ-01, ADQ-02, ADQ-03, ADQ-04, ADQ-05,
  ADQ-06, ADQ-09, ADQ-10, PRE-03 e INF-01. PRE-03 recibe su respuesta (dos grabaciones BIDS
  independientes). El test de regresión de ADQ-09 sigue pendiente.
- **`.gitignore`** excluye `*.xlsx`, `*.vhdr`, `*.vmrk` y `*.eeg` en cualquier carpeta (ADQ-10).
  `AGENTS.md`: «Hanning» corregido a Hamming-taper (INF-01) y cohorte longitudinal a 28 pares.
- **Informe:** cap. 4 con Cuadro 4.1, secciones 4.2.2, 4.3 y 4.5, figuras 4.5, 4.6 y 4.8
  actualizadas; el Cuadro 4.4 pasa a «dataset NeuroMIND».

### 2026-10-15 — Documento clínico unificado para Javier y Alejandro

- **`scripts/build_clinical_master.jl`** genera `data/sourcedata/clinical/datos_clinicos_unificados.xlsx`
  (Léeme · Participantes · Grabaciones EEG · Control de calidad · Diccionario) y `.csv`
  (hoja Participantes). Parte de `participants.tsv` (ya fusiona la base MIND-EM de
  28-09-2026 con EDSS/D-FIS/BDI-II numéricos), `inventory.csv` (decisiones de Javier sobre
  raws) y las observaciones de Javier. Las fuentes originales pasan a `data/sourcedata/clinical/anteriores/` (rutas de `audit_full_dataset.jl` actualizadas; salida idéntica). Verificado contra `data/bids/`: 308 combinaciones sujeto/sesión/condición y los 6 tramos de M04/M05/M06 T1 coinciden. Sin fecha de nacimiento ni nombres de archivo originales.
- **Control de calidad frente a las categorías de Javier** (hoja del 10-09, anterior a la base
  del 28-09): EDSS y D-FIS transversal son coherentes con los valores numéricos; quedan por
  confirmar M12 (D-FIS longitudinal Baja con T1=18, mediana del subgrupo 18), M12 vs M24 (mismos
  D-FIS T1/T2, distinta categoría longitudinal) y M18 vs M25 (mismos D-FIS, respuesta Mayor vs
  Menor), además de las edades de M13, M29 y M35.
- El `.xlsx` se escribe como SpreadsheetML + `zip` del sistema (sin dependencias nuevas).

### 2026-10-15 — `data/` reorganizado y dataset BIDS-EEG real (validado)

- **Nueva estructura** (todo bajo `.gitignore`): `data/sourcedata/{brainvision,clinical,docs}`
  (solo lectura), `data/bids/` (dataset BIDS-EEG), `data/tables/inventory.csv` y
  `data/_archive/` (copias previas, borrables tras verificar). `data/full_data/`,
  `data/clinical/`, `diagnostico_dataset.xlsx` y el BIDS ligero (`raw/`, `electrodes/`)
  salen de su sitio; nada se borró.
- **`build_bids_full.jl` reescrito:** genera `data/bids/sub-*/ses-*/eeg/` con el tríptico
  renombrado (se reescriben `DataFile=`/`MarkerFile=`; desaparecen las iniciales),
  `electrodes.tsv` + `coordsystem.json` por sesión y sidecars heredados en la raíz
  (`task-*_eeg.json`, `task-*_channels.tsv`). Plantillas en `config/bids/` (versionadas).
  M04/M05/M06 T1 se cortan **antes de filtrar** en dos grabaciones con intervalos
  semiabiertos `[tmin, tmax)` (decisión de Javier/Alejandro, ADQ-06): cierra PRE-03,
  el cargador no necesita ventana temporal. Comprueba que los 209 `.vhdr` comparten canales
  y orden (último canal Fp2). Sidecars corregidos: `EEGReference=FCz`, tierra = electrodo
  `Gnd` del gorro ECI (no Fpz), equipo actiCHamp (no BrainAmp).
- **Validado** con `bids-validator-deno` 3.0.2 (`uv tool install bids-validator-deno`):
  0 errores, 4 avisos documentados en `data/bids/README.md`. Señal BIDS ≡ `sourcedata`
  muestra a muestra en las 209 grabaciones (comparación con MNE).
- **`audit_full_dataset.jl`:** nuevas rutas; `inventory.csv` con `filepath` relativo (portable
  entre Macs); fusiona Excel MIND-EM, CSV demográfico y `Hoja1` de Javier en
  `participants.tsv` (24 columnas, `participants.json`): fechas de nacimiento e inclusión,
  EDSS/D-FIS/BDI-II T1 y T2, diagnóstico, tratamiento, categorías de Javier. Ya no genera
  `groups.csv` ni `longitudinal_pairs.csv`; `groups.example.csv` y `subgroups.csv` eliminados
  (absorbidos). Imprime discrepancias entre fuentes (edad M13, M29, M35).
- **Código:** `load_single_subject` y `detect_first_subject` leen el árbol BIDS (una sola vía
  de carga; se elimina la ruta TSV heredada, `load_eeg_bids` y `_parse_bids_json`).
  `Transversal.load_t1_cohort` y `Longitudinal.load_pairs` calculan cohorte y pares desde
  `inventory.csv`. Helpers de ruta `bids_eeg_dir/bids_vhdr_path/bids_electrodes_path` en
  `Config.jl`. `Project.toml`: dependencia nueva `XLSX` (lectura de la base clínica).
- **mne_brain:** lee directamente `../data/bids`; se elimina `scripts/build_bids_full.py`
  y la copia `mne_brain/data/BIDS` (archivada). Tests reescritos.
- **Privacidad:** `config/raw_selection.csv` (versionado en un repositorio público) ya no
  lleva iniciales ni fechas de grabación: `M4_T1_<iniciales>_<fecha>_COG-10-1` se escribe
  `M4_T1_*_COG-10-1` y la auditoría compara con `_anon_stem`. Inventario y `participants.tsv`
  resultan idénticos (mismo hash) antes y después.
- **Pendiente:** reprocesar M04/M05/M06 T1 y relanzar los análisis de grupo; las rutas del
  dashboard que leían señal desde `*_eeg_data.tsv` (visor de señal cruda, PSD cruda) ya no
  funcionaban con el BIDS ligero y siguen pendientes de portar al lector BrainVision;
  coordenadas de electrodos sin sistema de referencia documentado (declarado `Other`).

### 2026-10-14 — Preparación de datos según la hoja maestra de Javier (sin reprocesar)

- Nuevo `config/raw_selection.csv` (fuente única de decisiones de selección de raws;
  acciones `segment|use|discard|exclude`). Recoge: M4/M5/M6 T1 (raw único con abiertos
  + cerrados, ventanas en segundos de Javier; coinciden con los marcadores `New Segment`
  del `.vmrk`), M16 T1 EC (usar `_ojoscerrados`), M31 T2 (no analizar), M37 T2 (usar T2_2;
  el inventario previo estaba invertido).
- `audit_full_dataset.jl` lee ese CSV: se retiran las reglas que excluían M4/M5/M6 como
  Oddball; `inventory.csv` gana `tmin_s,tmax_s` (cols. 18–19, al final para no desplazar
  columnas que leen otros scripts). Efecto: 209 grabaciones válidas (antes 206); M5 pasa a
  par longitudinal completo, M31 sale del longitudinal, M4/M6 entran en transversal T1.
- `build_bids_full.jl`: nueva opción `--rebuild`; los JSON incluyen `tmin_s/tmax_s` cuando
  hay corte; `EEGReference` corregido de `Cz` a `FCz` (referencia hardware según el
  montaje de la imagen del estudio) + `EEGGround: Fpz`. JSON/electrodes de M31 T2 apartados
  a `data/bids/{raw,electrodes}/_descartados/` (nada se borra).
- `data/bids/subgroups.csv` (EDSS, DFIS T1 = col. D transversal, DFIS T1+T2 = col. E
  longitudinal, respuesta) exportado de la `Hoja1`. No mezclar D y E.
- Archivos de Javier reubicados en `data/clinical/`; PDF del mapa de dependencias y
  montaje de 32 canales en `docs/`.
- **Pendiente (requiere código del pipeline):** el cargador aún NO aplica `tmin_s/tmax_s`;
  hasta entonces no se debe lanzar el batch sobre M4/M5/M6 T1 (procesaría el raw
  completo como EC y como EO). Ver pendientes en el resumen de la sesión.

### 2026-10-07 — graphify versionado, parche Julia y CI

- `graphify-out/` versionado (grafo 1471 nodos / 3015 aristas / 64 comunidades,
  etiquetas en español) para compartirlo entre los dos equipos; se ignoran solo
  los archivos dependientes de máquina/ejecución y las copias fechadas.
  `.gitattributes` registra el merge driver de `graph.json`; hooks git
  post-commit/post-checkout reconstruyen el AST en segundo plano.
- graphify 0.9.79 ignoraba las funciones Julia con tipo de retorno
  (`function f(x)::T`, 253/762 funciones, p. ej. todo `Epochs.jl` y `FDR.jl`):
  parche local en `graphify-out/patches/` (aplicar en cada equipo).
- `CLAUDE.md` con sección graphify + hook PreToolUse (`.claude/settings.json`,
  ruta `$HOME/.local/bin/graphify`, portable entre equipos);
  `.cursor/rules/graphify.mdc` para Cursor.
- CI con GitHub Actions (PR #6): chequeo de datos versionados, syntax check y
  `test/runtests.jl` con Julia 1.12. Suite local: 449/449 pass.
- Plan de continuación (Report_Pre + revisión de resultados): `docs/PLAN_continuacion_2026-10.md`.

### 2026-07-28 (cont.) — Visor transversal con estadística única de producción

- Contrato transversal v2 obligatorio en `statistics_contract.toml`, con
  procedencia, ámbitos y tamaños FDR, bootstrap, cuantiles, Mann–Whitney,
  `effect_rrb` y `effect_d_pooled`; CSV/JSON guardan precisión completa y el
  redondeo queda exclusivamente en la presentación.
- `plot_transversal.jl` valida contrato, columnas, familias y artefactos antes
  de cargar. Los resultados antiguos o incompletos se muestran como
  incompatibles y deben regenerarse; no existen fallbacks estadísticos ni de
  matrices en el navegador.
- El productor persiste medias/IC, medianas/cuartiles, diferencia/IC, r_rb,
  probabilidad de superioridad y d pooled/IC. El visor separa aristas
  evaluadas, nominales, Top-N y FDR, y consume sus flags/rangos persistidos.
- Tras auditar consumidores longitudinales y transversales, se retiraron de
  `viewer_common.js` los estimadores estadísticos compartidos que ya no tenían
  llamadas. Potencia conserva BH entre canales dentro de cada banda y muestra
  cobertura por canal.

### 2026-07-28 (cont.) — Visor longitudinal con estadística única de producción

- Contrato longitudinal v2 obligatorio en `statistics_contract.toml`: registra
  fuente estadística, ámbitos/tamaños FDR, bootstrap, cuantiles, `effect_rrb` y
  `effect_dz`. Los artefactos se escriben con precisión completa y nombres
  canónicos (`effect_dz`, `effect_dz_ci_*`, `median_diff`, `q1_diff`, `q3_diff`);
  el redondeo queda exclusivamente en la presentación.
- `plot_longitudinal.jl` valida contrato, columnas y artefactos antes de cargar.
  Resultados antiguos/incompletos se muestran como incompatibles y deben
  regenerarse: no se reconstruyen matrices ni estadísticos como fallback.
- `Longitudinal.jl::paired_change_summary` centraliza cuantiles, `r_rb` pareado,
  Cohen dz e IC95% bootstrap. `network_global_statistics.csv` materializa todos
  esos valores, su expresión mean wPLI normalizada y la familia BH de 7 bandas;
  el visor longitudinal ya no recalcula estadísticos en JavaScript.
- Mean wPLI se presenta como normalización exacta de mean strength
  (`mean_strength=(n_channels−1)×mean_wPLI`), no como estimando independiente.
  ALPHA-EC queda fijado por regresión: `r_rb≈−0.1618`, mediana
  `≈−0.0051785`, 351 aristas, 20 nominales, Top-20 con solapamiento 17 y 0 FDR.
- La UI separa recuentos total/nominal/Top-N/FDR, explicita las tres familias
  de corrección, usa −log10(p) con color FDR en el volcano y representa en el
  topograma de potencia la cobertura por canal (BH entre 31 canales/banda).

### 2026-07-28 (cont.) — Suite de tests consolidada

- Eliminada la duplicidad `test/` + `tests/`: la suite completa vive ahora en
  la ubicación canónica Julia `test/runtests.jl`; ya no existe un wrapper que
  delegue a una segunda carpeta.
- `Project.toml` declara `Test` mediante `[extras]` y `[targets]`, por lo que
  tanto `Pkg.test()` como la ejecución directa usan el mismo archivo sin
  modificar manualmente `LOAD_PATH`. `Printf`, usado por los módulos de cohorte,
  queda además declarado explícitamente en `[deps]` (antes lo ocultaba ese
  `LOAD_PATH` global).
- Encabezado de la suite rehecho con el formato documental del proyecto,
  comandos de ejecución e inventario de los 36 `@testset` con su finalidad.

### 2026-07-28 — Auditoría y endurecimiento longitudinal

- Selección T1–T2 independiente por condición: un sujeto ya no queda fuera de
  EC por carecer de un par EO. La cohorte real pasa a 16 pares EC y conserva
  18 EO; M07 queda correctamente recuperado en EC.
- Wilcoxon signed-rank exacto condicional mediante programación dinámica para
  N≤30, con fallback asintótico documentado y método registrado en CSV, JSON y
  snapshot. La arista marginal que habría superado FDR al aplicar el exacto
  sobre el antiguo N=15 deja de hacerlo al recuperar M07: 0 aristas y 0 celdas
  de potencia superan FDR en EC/EO.
- `subject_band_means.csv` y el gráfico pareado usan ahora el mismo montaje
  común que el estimando C (`mean_wpli_equiv`), con 27 canales EC / 28 EO.
  `n_bands_ok` conserva el recuento real también para excluidos, y las tablas
  añaden los aliases explícitos `n_edges` y `mean_abs_dz`.
- Nueva figura `power_effect_heatmap_channel_band.png`, canales en orden
  anatómico y marcas naranjas para N<70 % del máximo de la condición; los
  topomapas incorporan la misma advertencia. Fp2 queda marcado con N=8/16 EC
  y N=4/18 EO.
- Forest, cambios individuales, matrices ALPHA y gráfico pareado ajustados para
  evitar recortes y hacer explícitos N, montaje común, prueba exacta y leyenda
  T1/T2. Suite completa y ejecución end-to-end verificadas.

### 2026-07-27 (cont.) — Figuras de síntesis nuevas en Transversal.jl

- `fig_raincloud_band`: generaliza el raincloud de wPLI-por-sujeto (antes solo
  ALPHA-EC de manuscrito, vía `fig_raincloud_alpha` que se conserva intacta)
  a las 7 bandas, en ambas condiciones — `figures/raincloud_{BAND}.png`.
- `fig_power_effect_heatmap`: heatmap canal×banda de Cohen d para potencia
  (EM−Control), con overlay de celdas FDR-significativas — resume los 14
  topomapas Δ-potencia individuales en una sola figura por condición —
  `figures/power_effect_heatmap_channel_band.png`.
- `fig_interaction_ec_eo`: interacción grupo×condición (EC−EO). Lee
  `subject_band_means.csv` de `eyesclosed/` y `eyesopen/` (ya escritos por el
  bucle principal — mismo principio CSV→figura del resto del módulo), los
  empareja por sujeto+banda y compara Δ(EC−EO) entre EM y Control
  (Mann–Whitney + FDR-BH entre las 7 bandas + bootstrap CI). Como no
  pertenece a ninguna condición sola, sale a una tercera carpeta
  `results/transversal/combined/{tables,figures}/` (forest del diff-of-diff +
  slope plots por sujeto), generada una vez tras procesar ambas condiciones.
- Manifiesto de trazabilidad (`figures_manifest.tsv`) ahora se escribe para
  **ambas** condiciones (antes solo `eyesclosed`, restringido a las 4
  figuras de manuscrito) — la variable `manifest` se inicializa al principio
  de cada iteración del bucle `for cond in [...]` en vez de dentro del
  bloque `if cond == "eyesclosed"`.
- Verificado contra datos reales: paridad numérica exacta con el estado
  previo (mismos p/q/d por banda en EC/EO); interacción ALPHA
  (ΔEC−EO=0.106 EM vs 0.063 Control) coincide con el hallazgo descriptivo
  que motivó la petición; ninguna banda sobrevive FDR en la interacción
  (q mínimo 0.786) — resultado honesto, no forzado.

### 2026-07-27 — Simplificación radical de la capa de cohorte (auditoría + rediseño)

Reestructuración en dos pasadas el mismo día: una relocalización inicial
preservando comportamiento (`src/viz/` → `src/visualization/`, flags
`--figures-only`/`--manuscript`), seguida de una auditoría completa
(evidencia por grep + 2 subagentes de exploración) y un rediseño desde
principios tras 3 ideas del usuario. Solo el estado final se documenta aquí
en detalle — ver historial de git para los pasos intermedios.

- **Estadística + figuras unificadas en un único módulo por dominio**,
  generadas siempre juntas en cada lanzamiento (sin modo `--figures-only`
  separado): `src/transversal/Transversal.jl` y `src/longitudinal/Longitudinal.jl`
  sustituyen `run_{transversal,longitudinal}_analysis.jl` (lógica) +
  `{Transversal,Longitudinal}Figures.jl` + `{Transversal,Longitudinal}Manuscript.jl`
  + `src/visualization/{GroupVizCommon,PublicationCommon,PublicationTheme,SummaryFigures}.jl`
  (8 archivos, ~1865 líneas con ~40 % de duplicación real, sustituidos por 2
  módulos de dominio de ~1250 líneas cada uno). `scripts/run_{transversal,longitudinal}_analysis.jl`
  quedan como lanzadores finos (`include` + `Transversal.run(config)`).
- **Duplicación aceptada a propósito entre los dos módulos de dominio**
  (primitivas de dibujo + constantes `CH_POS`/`BAND_ORDER`, ~150 líneas): es
  el precio de tener exactamente 2 archivos en vez de 2 + un módulo
  compartido — decisión explícita del usuario, documentada en el propio
  código como "nota de diseño".
- **`results/{transversal,longitudinal}/EC|EO/` → `{eyesclosed,eyesopen}/{tables,figures}/`**,
  alineado con la convención que ya usa `results/subjects/` (vía `norm_cond()`).
  Todos los CSV bajo `tables/` (antes sueltos en la raíz + `tables/spectral/`
  + `tables/network/` — ahora plano); todas las PNG bajo `figures/`
  (exploratorias + antiguas "manuscript", ya sin distinción de carpeta).
  Tocó también `src/webapp/App.jl` (paneles 13/14, más superficie de la
  esperada) y los dos visores interactivos.
- **`results/publication/` eliminada por completo** (Idea 3): sin PDF (solo
  PNG); las figuras de manuscrito viven en `results/{dominio}/eyesclosed/figures/`
  junto a las exploratorias; la única figura cross-dominio
  (`summary_transversal_longitudinal_effects`) vive en `results/summary/`
  (carpeta nueva y mínima), generada por el segundo de los dos análisis que
  se lance, comprobando en disco si el otro ya tiene resultados.
- **`reconstruct_longitudinal_C!` eliminada de raíz** (no solo relocalizada):
  `mean_strength` por sujeto×banda×T1/T2 se captura una única vez, dentro del
  propio bucle que ya calcula `network_global_statistics.csv`
  (`tables/mean_strength_scores.csv`) — ya no hay una segunda pasada que
  recalcule desde wPLI crudo ni asserts para verificar que dos caminos
  independientes coinciden, porque ahora solo hay un camino.
- **Visores interactivos reubicados a `src/interactive/`** (Idea 1): no son
  lanzadores, son una tercera categoría de artefacto (HTTP+JS en vivo) —
  `scripts/aux_viewers/` → `src/interactive/aux/`, `plot_{transversal,longitudinal}.jl`
  + `group_viewer_common.js` (→ `viewer_common.js`) desde `src/{transversal,longitudinal}/`
  y `src/visualization/`. Nuevo `src/interactive/viewer_support.jl` (ex
  `GroupVizCommon.jl`, sección IO/JSON) — detectado durante la migración que
  los visores dependían de ese archivo antes de que se borrara; corregido
  antes de completar el borrado.
- **Código muerto podado** (single-subject, confirmado por grep exhaustivo,
  0 llamadas): `src/visualization/ClinicalPlots.jl` (archivo entero),
  `Heatmaps.jl::plot_group_comparison`, `Spectra.jl::plot_spectrum`,
  `GraphPlots.jl::plot_graph_metrics`, `GraphPlots.jl::plot_surrogate_distribution`
  (duplicada por `SingleSubjectPipeline.jl::_plot_surrogate_null_band`, que sí
  se usa — se conserva esa). `src/report/HTMLReport.jl` (335 líneas, huérfano,
  0 llamadas, sin relación con `Report_Pre/`) y `src/dashboard/Dashboard.jl`
  (stub de 6 líneas nunca incluido) eliminados enteros.
- **Duplicación estadística de `_assign_ranks` (4 copias: `GroupStats.jl`
  testeada-pero-no-usada, `run_transversal_analysis.jl`, `run_longitudinal_analysis.jl`,
  `PublicationCommon.jl`) reducida a 2** (una por módulo de dominio) — no
  eliminada del todo a propósito: consolidar con `GroupStats.jl`/`FDR.jl`
  queda fuera de esta ronda (toca lógica estadística, `AGENTS.md` §12 pide no
  tocarla sin verificar contra `EEG_Julia/`).
- Limpieza de 7 PNG sin uso en `Report_Pre/figures/plots/` (6 exploratorios +
  1 QC legacy, confirmados por grep de `\includegraphics` en cada capítulo).
  Los 8 PNG de manuscrito que sí se usan (verificado 8/8 en LaTeX) quedan sin
  tocar pero **desactualizados** — la sincronización con `Report_Pre/` ya no
  es automática (`copy_to_report_plots` eliminado, Idea 3): es un paso manual
  (`cp`/`rsync` desde `results/{transversal,longitudinal}/eyesclosed/figures/`
  y `results/summary/` a `Report_Pre/figures/plots/`).
- Verificado extremo a extremo contra los datos reales del dataset (no solo
  sintéticos): ambos lanzadores, paridad numérica exacta con el estado previo
  a la sesión (ALPHA transversal ctrl=0.196/MS=0.260, 12 sig; longitudinal EC
  15 pares, asserts de reconstrucción implícitos por construcción), ambos
  visores cargando desde las nuevas rutas.

### 2026-07-25 (cont.) — Manuscript publication figures (EC)

- `scripts/export_publication_figures.jl` + `src/viz/publication/`: eight
  static manuscript figures (PNG+PDF) from existing transversal/longitudinal
  CSVs — forest plots, raincloud ALPHA, connectivity matrices with FDR,
  network+power, longitudinal individual Δ (estimand C), and a
  transversal–longitudinal summary. Bootstrap CI (B=5000, seed fixed);
  estimand A transversal / C longitudinal; reconstructs per-subject
  `mean_strength` with assert `ΔwPLI = Δstrength/(n_ch−1)`. Integrated in
  `Report_Pre/chapters/{12,13,10}_*.tex`.

### 2026-07-25 (cont.) — Audit, Report_Pre group chapters, aux viewers versioned

- Full audit of the Phase C batch and both group analyses:
  `docs/audits/batch_and_group_analysis_audit_2026-07-25.md` — channel hard-intersect
  quantified (24/31 transversal EC, 27–28/31 longitudinal), a longitudinal
  candidate-filter bug (`audit_full_dataset.jl:411` requires EC *and* EO complete before
  per-condition QC), `GroupStats.jl`/`FDR.jl` (the only tested statistics module) not
  actually used by `run_transversal_analysis.jl`/`run_longitudinal_analysis.jl`, an
  uncurated duplicate recording (M16 T1 EC), and cohort-count mismatches vs this file.
- Applied the audit's low-risk fixes: this file's §7 clarified (dynamic QC exclusion is
  independent of `exclude_fp2`), A&S citation corrected (7.1.26, not 26.2.17), `load_wpli`
  wrapped in try/catch in both group scripts, EO cohort figures regenerated
  (`regenerate_group_figures.jl EO` — previously only EC had been run).
- `figures/aux/` (12 single-subject viewers, previously inside `results/`, gitignored and
  at risk of being lost on regeneration) versioned to
  `scripts/aux_viewers/sub-M05_ses-T2_eyesclosed/`, with each script's path constants
  rewritten to resolve `results/` from the new location.
- `Report_Pre/chapters/12_transversal.tex` and `13_longitudinal.tex` added — the report's
  first group-level chapters (previous chapters covered only sub-M05).
  `10_conclusiones.tex` updated: dropped now-obsolete limitations ("batch not launched",
  "no group contrast yet citable") and added ones the audit surfaced (channel-intersection
  variability, the longitudinal filter bug, missing clinical metadata).

### 2026-07-25 — Phase C + cohort visualization overhaul
- Batch: 201 OK / 5 SKIP / 0 ERR (~51 min, surrogates OFF); transversal + longitudinal regenerated.
- `src/viz/GroupVizCommon.jl` + `group_viewer_common.js`: shared PNG/JS for cohort plots.
- Viewers (`plot_longitudinal.jl` :8780, `plot_transversal.jl` :8781): full `edge_stats` tooltips, shared color scales, FDR banner, volcano, SEM on strip/spaghetti, `Reverse(:RdBu)` aligned with UI.
- Analysis scripts write `heatmap_triplet_*`, topo `sig_network_*`, and `explore_network_topN_*` when FDR is empty; honest `best_band` (empty if `n_total_sig=0`).
- `scripts/regenerate_group_figures.jl` regenerates PNGs from CSV without re-running stats.

### 2026-05-26 — wPLI multi-method + script fixes
- `src/connectivity/wPLI.jl`: multi-method architecture (Hilbert/FourierCSD/Multitaper) using abstract types + multiple dispatch. Surrogates use the same estimator.
- `config/batch_pipeline.toml` and `config/single_subject.toml`: new `[connectivity]` block with `wpli_method`, `[connectivity.fourier_csd]`, `[connectivity.multitaper]`.
- `scripts/run_transversal_analysis.jl` and `run_longitudinal_analysis.jl`: fix `erf` → A&S polynomial approximation 7.1.26.
- `scripts/run_longitudinal_analysis.jl`: same `erf` fix.
- `scripts/run_longitudinal_analysis.jl`: fix `longitudinal_pairs.csv` reader (new boolean flag format).
- Segmentation: `profile = "default"`, `segment_length_seconds = 2.0` in batch config.

### 2026-05-25 — dwPLI, valid surrogates, GraphMetrics, QC v2
- Debiased dwPLI (Vinck 2011) added to `wPLI.jl`.
- `_circular_shift_surrogate` — methodologically correct surrogate (independent shift per channel).
- Monte Carlo correction (+1): p-values never 0.
- `qc_decision_table.csv` — automatic QC decision table.
- Fp2 exclusion across all subjects.

### 2026-05-24 — Full dataset, BrainVision, Batch pipeline
- `src/io/BrainVisionLoader.jl`: native reader IEEE_FLOAT_32 / INT_16.
- `scripts/audit_full_dataset.jl`: 212 vhdr audited, 0 errors.
- `scripts/build_bids_full.jl`: 205 metadata JSON created.
- `scripts/run_batch_pipeline.jl`: batch pipeline with CLI filters and CSV log.

### 2026-05-23 — Full dashboard (panels 0–14), PSD fix
- Panels 11–14 implemented (Results, Export, Transversal, Longitudinal).
- PSD normalization fixed: `Pseg ./= (n_samp * fs)`.
- Panel 15: MNE-Python validation — 3 stages, band power comparison, 6-metric per-channel table.
