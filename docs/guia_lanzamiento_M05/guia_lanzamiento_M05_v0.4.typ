// ═══════════════════════════════════════════════════════════════
//  Cuaderno de ejecución, interpretación y validación · sub-M05 / ses-T2 / eyesclosed
//  Versión v0.4 (10-11-2026): v0.3 + guardado de tablas y figuras en el Modo B (REPL).
//  10-11-2026: + Anexo C (cómo ver cada figura de M05 que cita el libro).
//  Compilar:  typst compile docs/guia_lanzamiento_M05/guia_lanzamiento_M05_v0.4.typ
//  Las figuras de figuras/ proceden de la corrida de prueba del 10-11-2026 (ver § Procedencia)
//  y NO se versionan (.gitignore local): son derivados de un paciente.
// ═══════════════════════════════════════════════════════════════

#let azul   = rgb("#1f3a5f")
#let azul2  = rgb("#3b6ea5")
#let azulf  = rgb("#eaf1fa")
#let gris   = rgb("#f3f4f6")
#let borde  = rgb("#d1d5db")
#let ambar  = rgb("#b45309")
#let ambarf = rgb("#fef3c7")
#let verde  = rgb("#166534")
#let verdef = rgb("#dcfce7")
#let rojo   = rgb("#991b1b")
#let rojof  = rgb("#fee2e2")
#let violeta  = rgb("#5b21b6")
#let violetaf = rgb("#ede9fe")

#set document(title: "Cuaderno de ejecución e interpretación · sub-M05", author: "NeuroMIND")
#set page(
  paper: "a4",
  margin: (x: 2cm, top: 2.2cm, bottom: 2cm),
  header: context {
    if counter(page).get().first() > 1 [
      #set text(8pt, fill: luma(110))
      NeuroMIND · Cuaderno de ejecución e interpretación · sub-M05 / ses-T2 / eyesclosed
      #h(1fr) v0.4 · 10-11-2026
      #v(-0.6em)
      #line(length: 100%, stroke: 0.4pt + borde)
    ]
  },
  footer: context {
    set text(8pt, fill: luma(110))
    h(1fr)
    counter(page).display("1 / 1", both: true)
  },
)
#set text(font: ("Helvetica Neue", "Helvetica"), size: 9.6pt, lang: "es")
#set par(justify: false, leading: 0.62em, spacing: 0.95em)
#show heading.where(level: 1): it => {
  set text(17pt, fill: azul, weight: "bold")
  block(above: 2.2em, below: 0.8em, sticky: true)[
    #it.body
    #v(-0.55em)
    #line(length: 100%, stroke: 1.2pt + azul2)
  ]
}
#show heading.where(level: 2): it => {
  set text(12pt, fill: azul, weight: "bold")
  block(above: 1.25em, below: 0.65em, it.body)
}
#show heading.where(level: 3): it => {
  set text(10pt, fill: azul2, weight: "bold")
  block(above: 1.0em, below: 0.5em, it.body)
}
#show raw: set text(font: ("Menlo", "DejaVu Sans Mono"), size: 7.9pt)
#show raw.where(block: false): it => box(
  fill: gris, inset: (x: 2.5pt), outset: (y: 2.5pt), radius: 2pt, it,
)
#show raw.where(block: true): it => block(
  width: 100%, fill: gris, inset: 7pt, radius: 3pt, stroke: 0.5pt + borde, it,
)
#show link: set text(fill: azul2)
#show math.equation.where(block: true): set block(above: 0.8em, below: 0.8em)
#set figure(supplement: [Figura], gap: 0.6em)
#show figure.caption: it => block(width: 100%, inset: (x: 4pt))[
  #set text(8.4pt, fill: luma(60))
  #set par(justify: false)
  *#it.supplement #context it.counter.display(it.numbering).* #it.body
]

#let caja(color, fondo, rotulo, cuerpo) = block(
  width: 100%, fill: fondo, inset: 8pt, radius: 3pt, breakable: true,
  stroke: (left: 3pt + color),
)[#text(fill: color, weight: "bold", rotulo) #h(0.3em) #cuerpo]
#let nota(c)    = caja(azul2, azulf, "Nota.", c)
#let aviso(c)   = caja(ambar, ambarf, "Atención.", c)
#let ok(c)      = caja(verde, verdef, "Comprobado.", c)
#let fallo(c)   = caja(rojo, rojof, "Problema conocido.", c)
#let lectura(c) = caja(violeta, violetaf, "Cómo leer la figura.", c)

// Estados de validación (leyenda en «Cómo usar este cuaderno»)
#let estados = (
  ok:   ("Ejecutado y comprobado",     verde, verdef),
  adv:  ("Ejecutado con advertencias", ambar, ambarf),
  pend: ("Pendiente de validación",    azul2, azulf),
  no:   ("No ejecutado",               rojo,  rojof),
)
#let estado(k) = {
  let (t, c, f) = estados.at(k)
  box(fill: f, stroke: 0.6pt + c, inset: (x: 4pt, y: 2.5pt), radius: 2pt,
      text(fill: c, weight: "bold", size: 7.8pt, t))
}

#let tabla(cols, cab, ..celdas) = {
  set text(8.4pt)
  set par(justify: false)
  table(
    columns: cols,
    stroke: 0.4pt + borde,
    inset: (x: 5pt, y: 4pt),
    fill: (_, y) => if y == 0 { azul } else if calc.even(y) { gris } else { white },
    ..cab.map(c => text(fill: white, weight: "bold", c)),
    ..celdas,
  )
}

// Cabecera de cada etapa: estado, paso, capítulo del libro y rutinas
#let cabecera(est, paso, libro, rutinas, inspeccion) = {
  set text(8.6pt)
  set par(justify: false)
  table(
    columns: (2.7cm, 1fr),
    stroke: 0.4pt + borde,
    inset: (x: 5pt, y: 4pt),
    fill: (x, _) => if x == 0 { azulf } else { white },
    [*Estado*], estado(est),
    [*Paso del pipeline*], paso,
    [*Libro*], libro,
    [*Rutinas*], rutinas,
    [*Inspección*], inspeccion,
  )
}

#let figura(ruta, ancho, pie) = figure(image(ruta, width: ancho), caption: pie)

// ─── Portada ─────────────────────────────────────────────────
#page(header: none, footer: none)[
  #v(2.6cm)
  #text(11pt, fill: azul2, weight: "bold", tracking: 1.5pt)[NEUROMIND · CONECTIVIDAD EEG EN ESCLEROSIS MÚLTIPLE]
  #v(0.5cm)
  #text(27pt, fill: azul, weight: "bold")[Cuaderno de ejecución,\ interpretación y validación]
  #v(0.3cm)
  #text(13.5pt, fill: luma(60))[Sujeto de referencia sub-M05: cómo lanzar cada etapa, qué calcula,\ qué significa su salida y hasta dónde está validada]
  #v(1.1cm)
  #line(length: 100%, stroke: 1.2pt + azul2)
  #v(0.4cm)
  #grid(
    columns: (4.4cm, 1fr), row-gutter: 0.75em,
    text(fill: luma(90))[Sujeto de referencia], [*sub-M05 · ses-T2 · eyesclosed* (EC)],
    text(fill: luma(90))[Código], [NeuroMIND, rama `feat/aclaraciones-clinicas` (commit `6a2089e`, maestro clínico vigente)],
    text(fill: luma(90))[Configuración], [`config/pipeline.toml` (perfiles `eeg_julia`, wPLI Hilbert clásico)],
    text(fill: luma(90))[Destino de resultados], [`~/Proyectos/NeuroMIND-results/{results, mne_results}`],
    text(fill: luma(90))[Estructura], [Capítulos del libro Quarto (`NeuroMIND-quarto/temario.yml`)],
    text(fill: luma(90))[Versión], [v0.4 · 10 de noviembre de 2026 (v0.3 + guardado de salidas en el Modo B)],
  )
  #v(1fr)
  #block(width: 100%, fill: gris, inset: 10pt, radius: 3pt)[
    #set text(8.8pt)
    *Qué es y qué no es.* Este cuaderno acompaña al libro: no repite su teoría, la enlaza. Para
    cada etapa del pipeline de un sujeto explica qué pregunta resuelve, qué objeto entra y sale,
    cómo se lanza cada subrutina, la ecuación que implementa el código, qué da en M05 y qué está
    comprobado y qué no. Los valores y figuras de M05 proceden de una corrida de prueba del
    10-11-2026 con el código y la configuración actuales; la corrida oficial en
    `NeuroMIND-results` *todavía no se ha lanzado*.
  ]
]

#outline(title: [Índice], depth: 2, indent: 1.2em)
#pagebreak()

// ═════════════════════════════════════════════════════════════
= Cómo usar este cuaderno

== Plantilla de cada etapa

Cada etapa del pipeline (capítulos 5 a 14 del libro) sigue la misma secuencia:

#tabla(
  (4.2cm, 1fr),
  ([Apartado], [Qué responde]),
  [1 · Objetivo científico], [Qué pregunta resuelve la etapa y por qué hace falta.],
  [2 · Entrada y salida], [Tipo de objeto, dimensiones, unidades y ficheros que deja.],
  [3 · Rutina], [Llamada mínima en el REPL y parámetros que importan; comando de inspección.],
  [4 · Fundamento], [La ecuación que implementa *el código* (no la de un libro de texto) y
    dónde se desarrolla en el libro.],
  [5 · Resultado de M05], [Tabla breve y, cuando existe, una figura real del sujeto.],
  [6 · Interpretación y validación], [Qué significa, qué se ha comprobado, qué no y por qué.],
)

== Leyenda de estados

#tabla(
  (5cm, 1fr),
  ([Estado], [Significado en este cuaderno]),
  estado("ok"), [Se ha ejecutado con el código actual y el resultado coincide con una
    referencia independiente (la corrida archivada de julio o una comprobación numérica).],
  estado("adv"), [Se ha ejecutado y produce salida, pero hay un problema conocido que limita
    su interpretación (se explica en la etapa).],
  estado("pend"), [Se ha ejecutado, pero falta la prueba que lo validaría (estabilidad,
    sensibilidad, N suficiente…). El protocolo está descrito.],
  estado("no"), [No se ha ejecutado en esta revisión, o no puede ejecutarse tal como está.],
)

#nota[«Reproducible» no es lo mismo que «validado». Que el pipeline dé el mismo número dos
veces demuestra reproducibilidad; que ese número mida lo que se pretende exige una validación
independiente (otra implementación, simulaciones, pruebas de sensibilidad). La columna de
estados separa las dos cosas.]

== Estado de cada etapa para M05

#tabla(
  (3.6cm, 4.4cm, 1fr),
  ([Etapa], [Estado], [Motivo principal]),
  [Carga · `[1/8]`], estado("ok"), [Señal idéntica a la de julio (|Δ| ≤ 10⁻⁴ µV, redondeo).],
  [QC · `[2/8]`], estado("ok"), [z-RMS idéntico a julio; Fp2 único canal marcado.],
  [Filtrado · `[3/8]`], estado("adv"), [Igual que julio; el Atlas registra 1 incumplimiento
    y 2 comprobaciones pendientes de `filter_recording`.],
  [ICA · `[4/8]`], estado("adv"), [No converge; clasificación
    relativa y cerca del umbral. Estabilidad: #estado("pend")],
  [Segmentación, línea base y rechazo · `[5/8]`], estado("ok"), [Misma época rechazada que
    en julio (la 12).],
  [Espectral · `[6/8]`], estado("adv"), [El «valor de banda» es una densidad media
    (µV²/Hz) etiquetada como µV²; bandas con solape y hueco.],
  [wPLI · `[7/8]`], estado("adv"), [DELTA por debajo del mínimo de ciclos; filtro por época
    de 1 s. Sensibilidad: #estado("pend")],
  [Surrogates · `[7/8]`], estado("pend"), [Solo N = 20 (prueba de llamada). N = 500 sin
    lanzar en esta revisión.],
  [Réplica MNE (fases 3–7)], estado("no"), [Sin entorno `.venv`; no se ha ejecutado.],
  [Validación Julia ↔ MNE], estado("no"), [`run_phase8_validation.py` lee rutas antiguas:
    fallará.],
  [Análisis de grupo], estado("no"), [Requiere la cohorte completa.],
)

== Procedencia de los valores y figuras

#tabla(
  (4cm, 1fr),
  ([Origen], [Detalle]),
  [Corrida de prueba], [10-11-2026, 17:32. Código `9e391e4`; `config/pipeline.toml` con dos
    cambios: rutas absolutas hacia una carpeta temporal y surrogates activos con N = 20. Su
    `pipeline_log.txt` y `config_snapshot.toml` están en `figuras/` junto a las figuras.],
  [REPL paso a paso], [Los fragmentos de este cuaderno se ejecutaron el mismo día y dieron
    los mismos números que la corrida.],
  [Topomapas IC7 e IC13], [Redibujados el 10-11-2026 con `plot_topomap` ya corregido, a
    partir de `ica_mixing_matrix.csv` de la corrida de prueba (no de la corrida oficial).],
  [Archivo de julio (`$ARCHIVO`)], [Corrida del 25-07-2026 (BIDS ligero, surrogates N = 500,
    ICA leída de una caché anterior). Solo lectura; se cita como referencia.],
)

#aviso[Ninguna figura de este cuaderno pertenece todavía a `NeuroMIND-results`. Cuando se
lance la corrida oficial, hay que regenerarlas desde allí y comprobar que coinciden.]

// ═════════════════════════════════════════════════════════════
= Preparación y lanzamiento

== Dónde se escriben los resultados

`load_ss_config` toma como raíz del proyecto la carpeta *dos niveles por encima* del TOML
(`config/pipeline.toml` → `NeuroMIND/`) y resuelve `[paths]` contra ella. Los enlaces
simbólicos del repositorio apuntan a la carpeta local `NeuroMIND-results`:

#tabla(
  (3.4cm, 1fr, 3.6cm),
  ([Enlace en el repo], [Destino real], [Contenido]),
  [`results`], [`~/Proyectos/NeuroMIND-results/results`], [Pipeline Julia (vacío)],
  [`mne_brain/results`], [`~/Proyectos/NeuroMIND-results/mne_results`], [Réplica MNE (vacío)],
  [`data`], [iCloud · `BRAIN/04_NeuroMIND/data`], [BIDS y `sourcedata` (solo lectura)],
  [—], [iCloud · `BRAIN/04_NeuroMIND/results`], [Archivo de julio (solo lectura)],
)

```bash
cd ~/Proyectos/NeuroMIND
readlink results               # → …/Proyectos/NeuroMIND-results/results
readlink mne_brain/results     # → …/Proyectos/NeuroMIND-results/mne_results
ls data/bids/sub-M05/ses-T2/eeg/
export M05=results/subjects/sub-M05/ses-T2/eyesclosed
export ARCHIVO="$HOME/Library/Mobile Documents/com~apple~CloudDocs/Universidad/\
Investigación/BRAIN/04_NeuroMIND/results/subjects/sub-M05/ses-T2/eyesclosed"
```

== Entornos (una vez)

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'include("src/NeuroMIND.jl"); println("OK")'
# réplica MNE (opcional):
cd mne_brain && uv sync && cd ..
```

== Configuración variante para surrogates

`config/pipeline.toml` ya apunta a M05/T2/eyesclosed con surrogates desactivados. La variante
con surrogates se genera fuera del repo, con rutas absolutas (necesarias porque la raíz se
deduce de la ubicación del TOML). Ya existe en `NeuroMIND-results/config/M05_surrogates.toml`.

```bash
mkdir -p ~/Proyectos/NeuroMIND-results/config
sed -e "s|^bids_root = .*|bids_root = \"$HOME/Proyectos/NeuroMIND/data/bids\"|" \
    -e "s|^results   = .*|results   = \"$HOME/Proyectos/NeuroMIND-results/results\"|" \
    -e '/^\[surrogates\]/,/^\[/ s|^enabled *= *false|enabled      = true|' \
    config/pipeline.toml > ~/Proyectos/NeuroMIND-results/config/M05_surrogates.toml
```

== Modo A · lanzamiento completo

```bash
# corrida base (sin surrogates), ~50 s + compilación
julia --project=. scripts/run_single_subject.jl
# corrida con surrogates N = 500, ~15 min más
julia --project=. scripts/run_single_subject.jl \
      --config ~/Proyectos/NeuroMIND-results/config/M05_surrogates.toml
```

#aviso[`--force` aparece en la documentación pero el script no lo lee. La caché de ICA se
invalida sola si cambian `[ica]` o `[filtering]`; para forzarla, borra `$M05/cache/`. Tampoco
hay limpieza previa de la carpeta: para una corrida limpia, `rm -rf $M05` antes de lanzar.]

== Modo B · paso a paso en el REPL

El orquestador es monolítico; para estudiar cada etapa se llaman sus rutinas en el mismo
orden, en una sola sesión de Julia. Las rutinas de cálculo solo devuelven objetos en memoria:
*no escriben tablas ni figuras*. Para grabarlas hay que llamar a las funciones de guardado del
orquestador (§ siguiente).

```julia
# julia --project=.   (desde ~/Proyectos/NeuroMIND)
using NeuroMIND, Statistics, DataFrames
cfg = load_ss_config("config/pipeline.toml")
```

Los fragmentos de cada etapa reutilizan las variables de la anterior.

== Modo B · guardar tablas y figuras desde el REPL

Quien graba las salidas en el Modo A es `run_single_subject_pipeline`, mediante funciones
internas `_save_*` (no exportadas: se llaman con el prefijo `NeuroMIND.`). En el REPL se
invocan a mano, al terminar cada etapa, sobre la misma carpeta que usaría el Modo A. Primero
se prepara el destino y el registro:

```julia
const NM = NeuroMIND
export_dir = joinpath(NM.results_dir(cfg), "subjects", "sub-M05", "ses-T2", "eyesclosed")
mkpath(joinpath(export_dir, "cache"))
log_io = open(joinpath(export_dir, "pipeline_log.txt"), "w")
dpi    = Int(get(cfg.export_cfg, "figure_dpi", 150))
val    = validate_channels(rec, NM.bids_electrodes_path(NM.bids_root_dir(cfg), "M05", "T2"))
```

Después, cada etapa graba lo suyo con las variables de sus fragmentos:

#tabla(
  (3.3cm, 1fr),
  ([Tras la etapa], [Llamada (deja en \$M05)]),
  [Raw (cap. 5)], [`NM._save_raw_signal(rec, export_dir, log_io)` → `tables/raw_signal.csv`],
  [Filtrado (cap. 6)], [`for (k, r) in pasos; NM._save_filtered_signal(r, export_dir, k, log_io); end`
    → `tables/filtered_signal_{notch,bandreject,highpass,lowpass}.csv`],
  [ICA (cap. 8)], [`NM._save_ica_results(ica, rec_filt, rec_ica, export_dir, t_ica, log_io;`
    `eval_df = ev, rejection_source = "auto",`
    `artifact_thresh = Float64(cfg.ica["artifact_threshold"]))` → `tables/ica/`, topomapas en
    `figures/ica/`, `json/ica_summary.json`. `t_ica` es la duración en segundos:
    `t_ica = @elapsed ica = run_ica(rec_filt, cfg)`],
  [Espectro (cap. 12)], [`NM._save_spectral_phase_results(sp, cfg, export_dir, dpi, log_io)` →
    PSD, potencia por banda y sus figuras],
  [wPLI (cap. 13)], [`NM._save_connectivity_phase_results(conn, cfg, export_dir, dpi, log_io)` →
    `tables/connectivity/`, `figures/connectivity/`],
  [Final], [`NM._save_all_results(rec, rec_filt, qc, bad, sp, conn, val, export_dir, "EC", cfg, dpi, log_io)`
    y `close(log_io)` → `channel_statistics`, `qc_summary`, `overview`, resto de figuras],
)

#ok[Probado el 10-11-2026 en M05/T2/eyesclosed sobre una carpeta temporal: con la cadena
raw → filtrado → ICA → épocas → PSD → wPLI y estas llamadas se escriben 23 tablas, 40 figuras y
3 JSON, y el registro.]

#aviso[El Modo B *no* graba: las tablas de segmentación y rechazo (`segments_table`,
`rejected_segments`…), porque dependen de objetos intermedios que solo construye el
orquestador; los surrogates; `config_snapshot.toml`; ni los índices globales
(`subjects_index.csv`, `qc/qc_decision_table.csv`). Para el conjunto completo y auditable,
usar el Modo A: el Modo B sirve para inspeccionar variables, no sustituye a la corrida oficial.
Si se combinan, el Modo A sobrescribe lo grabado en el B.]

#aviso[Los topomapas ICA corregidos (§ ICA) requieren el arreglo de `plot_topomap` de la PR #12
(`fix/topomap-transpuesto`), que todavía no está en `main` ni en `feat/aclaraciones-clinicas`.
Sin él, `figures/ica/ica_topomap_*.png` salen transpuestos.]

== Árbol de salida (corrida con surrogates)

```text
results/
├── subjects_index.csv · qc/qc_decision_table.csv       ← índices globales
└── subjects/sub-M05/ses-T2/eyesclosed/                  ← $M05
    ├── pipeline_log.txt · config_snapshot.toml · ica_labels_auto.csv
    ├── tables/   raw_signal · channel_statistics · qc_summary · overview
    │             filtered_signal_{notch,bandreject,highpass,lowpass}
    │             segments_table · channel_coverage · rejected_segments
    │             channel_artifact_summary · psd_by_channel · band_power_summary
    │             regional_psd · spectral_indices
    │   ├── ica/           ica_components · ica_component_features · ica_activations
    │   │                  ica_{mixing,unmixing}_matrix · ica_signal_{before,after}
    │   ├── connectivity/  wpli_{BANDA} · connectivity_edges · network_metrics
    │   └── surrogate/     wpli_{observed,pvalues,qvalues,significant}_{BANDA}
    │                      surrogate_null_stats_{BANDA} · significant_connections
    │                      surrogate_quality
    ├── json/     {ica,segmentation,artifact_rejection,spectral,connectivity,
    │              surrogate}_summary.json
    ├── figures/  filtered_signal_preview · psd_all_channels · band_power_summary
    │             ica/ica_topomap_NNN · connectivity/wpli_{BANDA}
    │             surrogate/surrogate_null_{BANDA}
    └── cache/    ica_result.jls · ica_config.hash · surrogate/null_distribution_{BANDA}.jls
```

// ═════════════════════════════════════════════════════════════
= Capítulos 1–4 · Contexto y datos

Los capítulos 1 (esclerosis múltiple) y 2 (estado del arte) no tienen rutinas. El capítulo 3
(fundamentos) solo pide conocer el montaje: 31 canales 10–20, referencia FCz en adquisición,
Fp2 en la última posición.

== Capítulo 4 · Organización del dataset (BIDS)

#cabecera("ok",
  [Fase 0, previa al pipeline (fases A y B).],
  [Cap. 4 · `sec-bids-organizacion`; Anexo A · flujo de ejecución.],
  [`audit_full_dataset.jl`, `build_bids_full.jl`, `build_clinical_master.jl`.],
  [`bids-validator-deno data/bids` (0 errores el 15-10-2026).],
)

El BIDS de M05 ya existe y no hace falta regenerarlo. Estas fases escriben en `data/` (iCloud),
no en `NeuroMIND-results`, y la fase A recorre siempre toda la cohorte.

```bash
julia --project=. scripts/audit_full_dataset.jl
julia --project=. scripts/build_bids_full.jl --subjects M05
bids-validator-deno data/bids
# comprobar M05 sin imprimir la columna 1 (nombre crudo con iniciales y fecha):
awk -F, 'NR>1 && $3=="M05" {print $7, $9, $12}' data/tables/inventory.csv
```

Resultado: cuatro grabaciones válidas (T1 EC, T1 EO, T2 EC, T2 EO). Las de T1 proceden de un
único registro recortado en la fase B, antes de cualquier filtro.

// ═════════════════════════════════════════════════════════════
= Capítulo 5 · Señal raw y control de calidad

#cabecera("ok",
  [`[1/8]` Carga y `[2/8]` QC de canales.],
  [Cap. 5 · `sec-raw-estadisticas`, `sec-rms-zscore`, `sec-raw-hfnoise`, `sec-raw-snr`.],
  [`load_single_subject`, `validate_channels`, `compute_channel_stats`,
    `compute_channel_spectral_qc`, `compute_correlation_summary`, `flag_bad_channels`.],
  [Visores `plot_raw.jl` (:8765), `plot_raw_butterfly.jl` (:8766),
    `plot_raw_histogram.jl` (:8768).],
)

=== 1 · Objetivo científico
Leer la señal tal como se registró y detectar canales cuyo comportamiento global es anómalo
*antes* de procesar nada, para no propagar un mal contacto a ICA, PSD y conectividad.

=== 2 · Entrada y salida
*Entrada:* tríptico BrainVision del árbol BIDS. *Salida:* `EEGRecording` de 31 × 50 180
muestras en µV (500 Hz, 100,36 s) y una tabla por canal. Ficheros: `tables/raw_signal.csv`,
`tables/channel_statistics.csv` (= `qc_summary.csv`), `tables/overview.csv`.

=== 3 · Rutina
```julia
rec = load_single_subject(cfg, "M05", "T2", "eyesclosed", 1)
qc  = compute_channel_stats(rec)                       # μ, RMS, σ, rango, z-RMS…
bad = flag_bad_channels(rec; z_threshold = cfg.qc["bad_channel_zscore_threshold"])
σm  = mean(std(rec.data[c, :]) for c in 1:rec.meta.n_channels)    # σ media cruda
```

=== 4 · Fundamento
Para cada canal $c$ con muestras $x_c [n]$, el criterio de canal anómalo es un z-score del
RMS calculado *entre canales del mismo registro*:

$ "RMS"_c = sqrt(1/N sum_n x_c [n]^2), quad
  z_c = ("RMS"_c - overline("RMS")) / s_"RMS", quad
  "malo" <=> abs(z_c) > 3 $

El aviso de amplitud compara la desviación típica media entre canales,
$overline(sigma) = 1/C sum_c sigma_c$, con 20 µV (posible filtro online apagado). Las métricas
extendidas están en el diccionario (Anexo A).

=== 5 · Resultado de M05
#tabla(
  (5.2cm, 1fr),
  ([Comprobación], [Resultado]),
  [Señal], [31 canales · 50 180 muestras · 500 Hz · 100,36 s],
  [Primer valor de Fz], [−7,19 µV (igual que en `$ARCHIVO`)],
  [$overline(sigma)$ cruda], [14,4 µV → `amplitude_warning = false`],
  [Canales con |z| > 3], [*Fp2* (z = 4,63), único],
)

=== 6 · Interpretación y validación
El QC *marca* Fp2 pero no lo elimina aquí; la exclusión se decide al construir el montaje de
análisis (capítulos 9–11). El z-score es relativo al propio registro: con 31 canales, un canal
muy ruidoso infla $s_"RMS"$ y puede enmascarar a un segundo canal moderadamente malo
(*masking*). Coincide con la corrida de julio canal a canal.

// ═════════════════════════════════════════════════════════════
= Capítulo 6 · Filtrado

#cabecera("adv",
  [`[3/8]` Filtrado de la señal *continua*.],
  [Cap. 6 · `sec-butterworth`, `sec-filtfilt`, `sec-cascada`; Atlas B · Filtrado.],
  [`filter_recording` (con `on_step`), `describe_filter_chain`.],
  [Visor `plot_filtered_vs_raw.jl` (:8769); figura `filtered_signal_preview.png`.],
)

=== 1 · Objetivo científico
Quitar la deriva lenta (< 0,5 Hz), la interferencia de red (50 Hz) y su armónico (100 Hz) y el
contenido por encima de 150 Hz, sin desplazar en el tiempo la actividad que interesa.

=== 2 · Entrada y salida
`EEGRecording` cruda → `EEGRecording` filtrada (mismas dimensiones y unidades). Ficheros:

#tabla(
  (7.2cm, 1fr),
  ([Fichero (bajo \$M05/)], [Contenido]),
  [`tables/filtered_signal_notch.csv`], [Señal tras el notch de 50 Hz.],
  [`tables/filtered_signal_bandreject.csv`], [Tras el bandreject de 100 Hz.],
  [`tables/filtered_signal_highpass.csv`], [Tras el pasa-alta de 0,5 Hz.],
  [`tables/filtered_signal_lowpass.csv`], [Señal final, tras el pasa-baja de 150 Hz.],
  [`figures/filtered_signal_preview.png`], [Cruda frente a filtrada, 5 canales, 10 s.],
)

=== 3 · Rutina
```julia
describe_filter_chain(cfg)               # orden, banda, orden del filtro y método
pasos = Dict{String,Any}()
rec_filt = filter_recording(rec, cfg; on_step = (k, r) -> pasos[k] = r)
```

=== 4 · Fundamento
Cascada Butterworth de orden 4 (perfil `eeg_julia`): notch 49,5–50,5 Hz y bandreject
99,5–100,5 Hz aplicados con `filt` (*causales*: introducen desfase, pero solo cerca de 50 y
100 Hz), y pasa-alta 0,5 Hz y pasa-baja 150 Hz con `filtfilt` (*fase cero*: la señal pasa
hacia delante y hacia atrás, el módulo efectivo es $abs(H(f))^2$ y el desfase se anula). La
teoría y las simulaciones de respuesta, bordes y escalón están en el capítulo 6
(Simulaciones A y C).

=== 5 · Resultado de M05
#figura("figuras/filtered_signal_preview_recorte.png", 46%)[Señal cruda (gris) y filtrada (azul) en
los cinco primeros canales, 0–10 s. Corrida de prueba del 10-11-2026.]

#lectura[El efecto dominante es el pasa-alta: en FT9 la señal cruda arranca en unos 120 µV y
deriva hacia cero, y en Fz hay un descenso lento de −20 µV; ambas componentes desaparecen sin
retrasar la actividad rápida (los picos de la curva azul caen donde los grises). Las
deflexiones de FT9 entre 4 y 7 s *sobreviven* al filtrado porque están en la banda útil: no
son un fallo del filtro, son un transitorio que tendrá que tratar la ICA (componente IC13 del
capítulo 8). La figura no muestra la respuesta en frecuencia ni los bordes: diez segundos y
cinco canales son una comprobación visual, no una validación.]

=== 6 · Interpretación y validación
La salida coincide con la de julio salvo el redondeo del CSV (|Δ| ≤ 10⁻⁴ µV). El estado es
«con advertencias» porque la ficha de `filter_recording` en el Atlas registra 1 comprobación
superada, 1 con incumplimiento y 2 pendientes (ver `atlas/filtrado.qmd`), y esa validación no
se ha cerrado. El `filter_order = 8` de `[connectivity]` es el del filtro de banda de wPLI
(capítulo 13), no el de esta cascada.

// ═════════════════════════════════════════════════════════════
= Capítulo 7 · Origen de los artefactos

Capítulo sin paso propio. En M05 se ven tres tipos, cada uno detectado en una etapa distinta:

#tabla(
  (4.2cm, 1fr, 3.8cm),
  ([Artefacto], [Evidencia en M05], [Dónde se trata]),
  [Canal ruidoso (Fp2)], [z-RMS = 4,63; además domina la columna de mezcla de IC7 (peso
    −11,8 frente a < 2,8 del resto).], [QC (marca) → ICA → montaje (lo excluye)],
  [Transitorio focal (FT9)], [Deriva de ~120 µV y deflexiones a 4–7 s; domina IC13 (peso
    −5,5 frente a 0,53 del siguiente).], [ICA],
  [Alta frecuencia lateral (T8, F4/FC6)], [IC15 e IC31 con potencia relativa alta por encima
    de 30 Hz.], [ICA (etiquetados `line_noise`)],
  [Amplitud puntual en C4], [Época 12: máx 72,2 µV > 70 µV.], [Rechazo por amplitud],
)

// ═════════════════════════════════════════════════════════════
= Capítulo 8 · ICA

#cabecera("adv",
  [`[4/8]` ICA sobre la señal continua filtrada (siempre antes de segmentar).],
  [Cap. 8 (previsto en el temario). Atlas: sin página todavía.],
  [`run_ica`, `compute_ica_features`, `evaluate_ica_components`, `load_ica_labels`,
    `write_ica_labels_auto`, `apply_ica_rejection`.],
  [Visores `plot_ica_components.jl` (:8770), `plot_ica_before_after.jl` (:8771).],
)

=== 1 · Objetivo científico
Separar la señal en fuentes estadísticamente independientes para eliminar las que son
artefacto (ojos, músculo, electrodo) conservando la actividad cerebral en todos los canales.

=== 2 · Entrada y salida
`rec_filt` (31 × 50 180) → `ICAResult` con la matriz de mezcla $A$ (31 × 31), la de desmezcla
$W$ y las activaciones $S$ (31 × 50 180), y la señal reconstruida sin los componentes
rechazados. Ficheros en `tables/ica/`, `json/ica_summary.json`, `figures/ica/` y
`ica_labels_auto.csv`.

=== 3 · Rutina
```julia
t_ica = @elapsed ica = run_ica(rec_filt, cfg)   # eeg_julia: 31 comp., 512 iter., tol 1e-7
ica.diagnostics                      # converged, n_iter, final_error…
feat = compute_ica_features(Float64.(ica.mixing_matrix), Float64.(ica.activations),
                            rec_filt.meta.fs, rec_filt.meta.channel_names)
ev   = evaluate_ica_components(feat; artifact_thresh = cfg.ica["artifact_threshold"])
rej  = Int.(ev.component[ev.artifact_type .!= "brain"])
rec_ica = apply_ica_rejection(rec_filt, ica, rej)
```

Precedencia del rechazo: un `$M05/ica_labels.csv` manual (columnas `component`, `label`) gana
siempre; si no existe y `auto_reject = true`, se rechazan los componentes con
`artifact_type ≠ brain`.

=== 4 · Fundamento
Modelo lineal $bold(x) = A bold(s)$. Tras centrar y blanquear con PCA
($bold(z) = Lambda^(-1\/2) E^T bold(x)$), FastICA simétrico con $g = tanh$ itera

$ W^+ = EE[g(W bold(z)) bold(z)^T] - "diag"(EE[g'(W bold(z))]) W, quad
  W <- (W W^T)^(-1\/2) W $

y se detiene cuando $max_k abs(1 - abs(bold(w)_k^+ dot bold(w)_k)) < 10^(-7)$. La señal limpia
es $hat(bold(x)) = A_"keep" bold(s)_"keep"$.

La clasificación calcula 7 rasgos por componente, los convierte en *z-scores entre los 31
componentes* y combina cuatro puntuaciones:
ocular $= 0,4 z_"front" + 0,4 z_"blink" + 0,2 z_"kurt"$;
muscular $= 0,5 z_"emg" + 0,3 z_"temp" + 0,2 z_"kurt"$;
línea $= z_"line"$; salto $= 0,5 z_"kurt" + 0,5 z_"extreme"$.
Si la mayor supera 1,5, el componente recibe la etiqueta de esa puntuación.

=== 5 · Resultado de M05
#tabla(
  (0.9cm, 1.6cm, 2.7cm, 1fr, 2.2cm),
  ([IC], [Var.], [Canal dominante (peso)], [Rasgos decisivos], [Etiqueta]),
  [7], [11,45 %], [Fp2 (−11,8)], [blink_ratio 13,6 · curtosis 21,6 · ocular 2,04 · salto *2,87*],
    [`jump`],
  [13], [1,41 %], [FT9 (−5,5)], [frontal_ratio 7,9 · blink 8,5 · curtosis 28,3 · ocular 3,68 ·
    salto *4,26*], [`jump`],
  [15], [0,84 %], [T8 (−2,8)], [line_ratio 0,0027 · muscular 1,40 · línea *2,02*],
    [`line_noise`],
  [31], [0,46 %], [F4 (2,3), FC6 (−1,6)], [emg_ratio 1,24 · muscular 1,77 · línea *3,01*],
    [`line_noise`],
)
Se elimina el 14,2 % de la varianza (retenida: 85,8 %). La FastICA *no converge*: agota las
512 iteraciones con error 3,2·10⁻⁷ frente a la tolerancia 10⁻⁷.

#grid(columns: (1fr, 1fr), gutter: 10pt,
  figura("figuras/topomap_IC7_corregido.png", 92%)[IC7 con `plot_topomap` corregido: el
    máximo cae sobre *Fp2*, el electrodo que domina su columna de mezcla (−11,8).],
  figura("figuras/topomap_IC13_corregido.png", 92%)[IC13 corregido: foco único sobre *FT9*
    (−5,5), frontotemporal izquierdo; el resto del cuero cabelludo, casi a cero.],
)

#lectura[Ambos mapas muestran un foco en *un solo electrodo*, sin el gradiente amplio que
produciría una fuente cortical ni el patrón bifrontal simétrico de un parpadeo: es la firma de un
artefacto de contacto o de canal, coherente con el z-RMS de Fp2 y con la deriva de FT9 del
capítulo 6. El signo es arbitrario en ICA (A y S pueden cambiar de signo a la vez). Los mapas se
han redibujado con la matriz de mezcla de la corrida de prueba; IC15 e IC31 no tienen figura en
esta versión.]

#caja(verde, verdef, "Corregido el 10-11-2026.")[Hasta esta fecha los topomapas salían
*transpuestos*: `plot_topomap` (`src/visualization/Topomaps.jl`) llenaba la rejilla como
`Z[j, i]` (fila = y) y `heatmap!` la lee como `Z[x, y]`, así que el campo quedaba reflejado
respecto a la diagonal mientras las etiquetas se dibujaban bien (IC7 aparecía junto a FT10/T8;
IC13, abajo a la derecha). La interpolación está ahora en `_topomap_idw_grid`, con la convención
`Z[i, j]` ↔ `(xg[i], yg[j])` y el test `Topomap_IDW_Orientation`. *Todos los
`ica_topomap_NNN.png` generados antes del 10-11-2026 están transpuestos* (también los de
`$ARCHIVO`); se corrigen al relanzar el pipeline. Las tablas y la clasificación no cambian: no
usan las figuras.]

=== 6 · Interpretación y validación
*Por qué se rechaza cada componente.* IC7 e IC13 son artefactos *de un solo electrodo*: casi
todo su peso cae en Fp2 (el canal ya marcado por el QC) y en FT9 (el que deriva en la figura
del capítulo 6), con potencia concentrada por debajo de 4 Hz y curtosis muy alta. La etiqueta
`jump` y no `eye/blink` sale de que la curtosis entra en las dos puntuaciones y la de salto
gana por poco (2,87 frente a 2,04 en IC7); para el rechazo da igual, para describir el
artefacto no. IC15 e IC31 se etiquetan `line_noise` porque su `line_ratio` es el mayor *de
este registro*, pero en valor absoluto la banda 48–52 Hz es solo el 0,3 % de su potencia: tras
el notch, lo que los distingue es la alta frecuencia lateral (puntuación muscular 1,4–1,8), más
propia de actividad muscular que de red.

*Fragilidad de la decisión.* Las puntuaciones son relativas (z-scores entre componentes): un
registro limpio también tendrá «el componente más raro». Hay componentes justo por debajo del
umbral que se conservan: IC9 (5,8 % de varianza, línea 1,44), IC3 (línea 1,41), IC28 (muscular
1,27) e IC14 (4,4 %, blink_ratio 12,3 pero ocular 0,97). Con un umbral de 1,4 en vez de 1,5
se eliminaría un 6,6 % más de varianza.

*Convergencia.* Que el error final (3,2·10⁻⁷) esté cerca de la tolerancia no demuestra que la
solución sea estable. La corrida de julio leyó una ICA antigua de caché (convergida en 290
iteraciones, 5 rechazados, 82,6 % retenido): las dos no son comparables componente a
componente, porque el orden de las fuentes de FastICA es arbitrario.

#caja(azul2, azulf, "Protocolo de estabilidad (pendiente).")[Repetir la ICA con el perfil
`default` (31 componentes, tol 10⁻⁷, `max_attempts = 1`) con 512 y 5000 iteraciones y semillas
1234, 1, 7, 42 y 2026. Para cada repetición: (a) convergencia e iteraciones; (b) emparejar cada
componente rechazado de referencia con el de mayor |correlación| de activaciones y anotar su
etiqueta; (c) diferencia relativa de la señal reconstruida, $norm(hat(x) - hat(x)_"ref") \/ norm(hat(x)_"ref")$;
(d) correlación y diferencia media absoluta de las aristas wPLI por banda. Añadir como cota la
corrida *sin* ICA. Criterio propuesto: los cuatro rechazos se reproducen (|r| > 0,9) y la wPLI
ALPHA cambia menos de 0,02 de media. El script está preparado y no se ha lanzado.]

// ═════════════════════════════════════════════════════════════
= Capítulos 9–11 · Segmentación, línea base y rechazo

#cabecera("ok",
  [`[5/8]`: montaje → `segment_recording` → `apply_baseline` → `reject_artifacts` →
    `apply_baseline`.],
  [Caps. 9, 10 y 11 (previstos en el temario).],
  [`segment_recording`, `apply_baseline`, `reject_artifacts`, `compute_epoch_quality_report`,
    `compute_channel_coverage`.],
  [Visores `plot_epochs.jl` (:8772) y `plot_baseline.jl` (:8773); necesitan la caché de ICA.],
)

=== 1 · Objetivo científico
Cortar la señal limpia en épocas cortas y estacionarias sobre las que promediar el espectro y
la fase, y descartar las que contienen amplitudes no fisiológicas.

=== 2 · Entrada y salida
`rec_ica` (31 canales) → montaje de análisis (30 canales) → `EpochSet` de 30 × 500 × 100
→ 99 épocas válidas. Ficheros: `json/segmentation_summary.json`,
`json/artifact_rejection_summary.json`, `tables/segments_table.csv`,
`tables/channel_coverage.csv`, `tables/rejected_segments.csv`,
`tables/channel_artifact_summary.csv`.

=== 3 · Rutina
```julia
# montaje: exclude_channels − Fp2 (porque exclude_fp2 = false) ∪ canales malos del QC
excl = String.(cfg.montage["exclude_channels"])
cfg.montage["exclude_fp2"] || (excl = filter(!=("Fp2"), excl))
excl = unique(vcat(excl, bad))                          # → ["Fp2"]
keep = findall(c -> c ∉ excl, rec_ica.meta.channel_names)
m = rec_ica.meta
rec_seg = EEGRecording(
    RecordingMeta(m.subject_id, m.session_id, m.condition, m.run, m.fs, length(keep),
                  m.channel_names[keep], m.channel_positions, m.bids_path),
    rec_ica.data[keep, :], rec_ica.times)
ep0 = segment_recording(rec_seg, cfg)    # cap. 9: 1,0 s sin solape (forzado por eeg_julia)
ep1 = apply_baseline(ep0, cfg)           # cap. 10: pasada 1
ep2 = reject_artifacts(ep1, cfg)         # cap. 11: ±70 µV
ep  = apply_baseline(ep2, cfg)           # cap. 10: pasada 2, solo épocas válidas
```

=== 4 · Fundamento
Época $e$ de longitud $L = 500$ muestras: $x_(c,e)[n] = x_c [(e-1)L + n]$. Línea base
`first_window_mean`: se resta la media de las 50 primeras muestras (0–100 ms),
$x_(c,e)[n] <- x_(c,e)[n] - 1/50 sum_(m=1)^(50) x_(c,e)[m]$. Rechazo: la época se descarta si
algún canal sale de $[-70, 70]$ µV en cualquier muestra. La calidad que lista
`segments_table.csv` es $q_e = max(0, 1 - 0,5 dot max abs(x) \/ 70)$ y es solo descriptiva.

=== 5 · Resultado de M05
#tabla(
  (5.2cm, 1fr),
  ([Comprobación], [Resultado]),
  [Montaje de análisis], [30 canales: Fp2 excluido por el QC dinámico (no por configuración)],
  [Épocas], [100 totales · *99 válidas* · 1 rechazada],
  [Época rechazada], [12 (11–12 s), peor canal C4, máx 72,18 µV, calidad 0,484],
  [Coincide con `$ARCHIVO`], [Sí: misma época y canal (julio: 72,56 µV, por la ICA distinta)],
)

=== 6 · Interpretación y validación
*Montaje.* `exclude_fp2 = false` conserva Fp2 por configuración, pero los canales con |z| > 3
se excluyen siempre: el análisis se hace con 30 canales y 435 aristas por banda, no 31 y 465.

*Línea base en reposo.* En EEG de reposo no hay estímulo: restar la media de los primeros 100
ms desplaza cada época por un valor arbitrario. No afecta a la PSD (que resta la media de cada
época antes de la FFT) ni, en principio, a wPLI (el filtro de banda elimina la continua), pero
*sí* a la decisión de rechazo: una época puede cruzar ±70 µV o no según cómo empiece. Es una
herencia de `EEG_Julia` documentada, no una elección óptima para reposo.

*Rechazo.* El umbral se aplica a los `n_channels_used` primeros canales *por posición*; con 30
canales se aplican todos. Falta una figura de amplitud por época (la genera el visor :8772).

// ═════════════════════════════════════════════════════════════
= Capítulo 12 · Análisis espectral

#cabecera("adv",
  [`[6/8]` PSD por canal promediada sobre las épocas válidas y valor por banda.],
  [Cap. 12 (previsto). Cap. 5 · `sec-raw-psd` para la PSD sin filtrar.],
  [`compute_psd`, `band_power`.],
  [Visores `plot_spectral.jl` y `plot_spectral_PSD.jl` (ambos :8767, no a la vez).],
)

=== 1 · Objetivo científico
Describir cuánta varianza tiene la señal en cada frecuencia y resumirla en bandas fisiológicas,
para caracterizar el ritmo dominante y comparar sujetos.

=== 2 · Entrada y salida
`EpochSet` 30 × 500 × 99 → `SpectralResult`: PSD de 30 × 513 (0–250 Hz, Δf = 0,488 Hz) en
*µV²/Hz* y un valor por canal y banda. Ficheros: `tables/psd_by_channel.csv`,
`tables/band_power_summary.csv`, `tables/regional_psd.csv`, `tables/spectral_indices.csv`,
`json/spectral_summary.json`, `figures/psd_all_channels.png`,
`figures/band_power_summary.png`.

=== 3 · Rutina
```julia
sp = compute_psd(ep, cfg)                       # nfft = 1024, ventana Hamming en un 10 %
sp.freqs[2] - sp.freqs[1], length(sp.freqs)     # (0.488, 513)
Dict(b => mean(v) for (b, v) in sp.band_power)  # valor medio de banda en los 30 canales
```

=== 4 · Fundamento
Para cada canal y época se resta la media, se aplica una ventana $w$ que es rectangular salvo
un 5 % en cada extremo (coseno de Hamming), se rellena con ceros hasta $N_"FFT" = 1024$ y

$ P_(c,e)(f_k) = 2 / (f_s N overline(w^2)) abs(sum_(n=0)^(N-1) w_n x_(c,e)[n]
  e^(-i 2 pi k n \/ N_"FFT"))^2 , quad N = 500, f_k = k f_s \/ N_"FFT" $

(sin el factor 2 en $k = 0$ y en Nyquist), que se promedia sobre las 99 épocas. El valor de
banda es la *media* de los bins de la banda, no su integral:

$ B_c ("banda") = 1 / abs(K) sum_(k in K) P_c (f_k), quad
  K = {k : f_"inf" <= f_k < f_"sup"} quad ["µV²/Hz"] $

La potencia integrada sería $approx B_c dot abs(K) dot Delta f$ [µV²].

#aviso[*Unidades.* La PSD y el valor de banda son densidades (µV²/Hz). La columna
`power_uv2` de `psd_by_channel.csv`, la tabla `band_power_summary.csv`, la figura
`band_power_summary.png` (eje «µV²») y la v0.1 de esta guía los etiquetan como µV²: es un error de
etiqueta, no de cálculo. La normalización $2 \/ (f_s N overline(w^2))$ sí da µV²/Hz y no
depende del relleno con ceros (corrección del 23-05-2026).]

=== 5 · Resultado de M05
#tabla(
  (2.3cm, 1.9cm, 1fr, 1fr, 1.2cm, 1.9cm),
  ([Banda], [Límites (Hz)], [Valor medio (µV²/Hz)], [Bins · ancho efectivo], [Ciclos en 1 s], [Integral aprox. (µV²)]),
  [DELTA], [0,5–4], [3,787], [7 · 3,42 Hz], [0,5–4], [12,9],
  [THETA], [4–8], [1,961], [8 · 3,91 Hz], [4–8], [7,7],
  [ALPHA], [7,8–11,7], [3,593], [8 · 3,91 Hz], [7,8–11,7], [14,0],
  [BETA_LOW], [12–15], [0,909], [6 · 2,93 Hz], [12–15], [2,7],
  [BETA_MID], [15–18], [0,884], [6 · 2,93 Hz], [15–18], [2,6],
  [BETA_HIGH], [18–30], [1,061], [25 · 12,2 Hz], [18–30], [13,0],
  [GAMMA], [30–50], [0,188], [41 · 20,0 Hz], [30–50], [3,8],
)
#text(8pt, fill: luma(90))[Media de 30 canales y 99 épocas. «Integral aprox.» = valor medio ×
nº de bins × Δf: es una estimación para leer la tabla, *no* la calcula el pipeline.]

#figura("figuras/psd_all_channels.png", 100%)[PSD de cada canal, 0–50 Hz, escala lineal y
eje vertical propio en cada panel (µV²/Hz, aunque el pipeline no lo rotule). Corrida de
prueba del 10-11-2026.]

#lectura[Tres rasgos. (1) Un pico α hacia 10 Hz, máximo en la región parieto-occipital
(O1, Oz, O2, P7, P8, TP9): es el ritmo esperado con ojos cerrados. (2) Un segundo pico hacia
20 Hz en los centrales (C3, CP1, CP2, Cz): actividad β sensoriomotora. (3) Un pico a ~1 Hz en
*todos* los canales y enorme en C4 (≈ 60 µV²/Hz, una escala diez veces mayor que la de sus
vecinos): C4 es también el canal de la época rechazada. Ese pico es actividad lenta que
sobrevive al pasa-alta de 0,5 Hz y explica que DELTA tenga el mayor valor medio. Como los ejes
son independientes, la figura sirve para ver la *forma* de cada espectro, no para comparar
amplitudes entre canales.]

=== 6 · Interpretación y validación
*Por qué FFT de 1024 con épocas de 500 muestras.* Rellenar con ceros reduce el espaciado de
la rejilla (Δf = 0,488 Hz) pero no añade información: la resolución real la fija la
duración de la época, $f_s \/ N = 1$ Hz, ensanchada por la ventana. Se mantiene por
coherencia con BrainVision Analyzer y `EEG_Julia`. Con un 90 % de ventana rectangular la fuga
espectral es alta: el pico de 1 Hz «derrama» potencia a los bins vecinos.

*Bandas.* Los intervalos son semiabiertos y no son contiguos: el bin de 7,81 Hz cuenta en
THETA *y* en ALPHA (solape 7,8–8 Hz) y el de 11,72 Hz no cuenta en ninguna (hueco 11,7–12 Hz).

*Comparar bandas.* Los valores medios no se pueden comparar entre bandas de anchura distinta
sin más: BETA_HIGH tiene un valor medio bajo (1,06) pero, integrado sobre sus 12 Hz, acumula
tanta potencia como DELTA. Ninguno de estos números tiene, por sí solo, significado clínico: un
sujeto no tiene grupo de comparación.

// ═════════════════════════════════════════════════════════════
= Capítulo 13 · Conectividad funcional (wPLI)

#cabecera("adv",
  [`[7/8]` wPLI en espacio de sensores (`use_csd = false`).],
  [Cap. 13 (previsto); cap. 2 · definición breve de PLI y wPLI.],
  [`compute_wpli` (estimador `hilbert`; alternativas `fourier_csd` y `multitaper`).],
  [Visor `plot_connectivity.jl` (:8774).],
)

=== 1 · Objetivo científico
Medir si dos canales mantienen una relación de fase *consistente y con desfase distinto de
cero*. El desfase nulo es lo que produce la conducción de volumen (una misma fuente vista por
dos electrodos); al descartarlo, wPLI mide acoplamiento menos contaminado por ella.

=== 2 · Entrada y salida
`EpochSet` 30 × 500 × 99 → `ConnectivityMatrix`: 7 matrices simétricas 30 × 30 en [0, 1],
diagonal 0, 435 aristas por banda. Ficheros: `tables/connectivity/wpli_{BANDA}.csv`,
`tables/connectivity/connectivity_edges.csv` (3045 filas), `network_metrics.csv`,
`json/connectivity_summary.json`, `figures/connectivity/wpli_{BANDA}.png`.

=== 3 · Rutina
```julia
conn = compute_wpli(ep, cfg)          # avisa: DELTA 0,5 ciclos/época < 4
W = conn.matrices["ALPHA"]
# comparar estimadores sin tocar el TOML:
c2 = deepcopy(cfg); c2.connectivity["wpli_method"] = "multitaper"
compute_wpli(ep, c2).matrices["ALPHA"]
```

=== 4 · Fundamento
Estimador `hilbert`: cada época de cada canal se filtra en la banda con un Butterworth de orden
8 aplicado con `filtfilt`, y se obtiene la señal analítica
$z_c (t) = x_c^"banda" (t) + i cal(H){x_c^"banda"}(t)$. Para cada par $(i, j)$:

$ I_(i j)(e, t) = "Im"(z_i (e,t) z_j^ast (e,t)) = A_i A_j sin(phi_i - phi_j), quad
  "wPLI"_(i j) = abs(sum_(e,t) I_(i j)(e,t)) / (sum_(e,t) abs(I_(i j)(e,t))) $

sumando sobre todas las muestras de todas las épocas. Vale 1 si el desfase tiene siempre el
mismo signo y 0 si los signos se compensan; las muestras con desfase cercano a 0 o π pesan
poco porque su seno es pequeño. Es la definición de Vinck et al. (2011) aplicada muestra a
muestra en lugar de sobre el espectro cruzado por época. `use_dwpli = true` daría el
estimador sin sesgo (rango −1…1).

=== 5 · Resultado de M05
#tabla(
  (2.4cm, 1fr, 1fr, 1fr, 1fr, 1fr, 1fr, 1fr),
  ([wPLI], [DELTA], [THETA], [ALPHA], [β low], [β mid], [β high], [GAMMA]),
  [media], [0,109], [0,129], [*0,338*], [0,143], [0,116], [0,087], [0,085],
  [máximo], [0,563], [0,383], [*0,864*], [0,541], [0,347], [0,366], [0,232],
)

#figura("figuras/wpli_ALPHA.png", 78%)[Matriz wPLI ALPHA (7,8–11,7 Hz) de M05, 30 × 30.
Cada celda $(i, j)$ es la consistencia del signo del desfase entre los dos canales. Corrida de
prueba del 10-11-2026.]

#lectura[El eje vertical está invertido respecto a la convención habitual (Fz abajo) y la
escala de color llega a 0,86, no a 1. El bloque claro une O1/Oz con O2, P4, P8 y CP6: un
acoplamiento α parieto-occipital *interhemisférico* (máximo 0,864 en O1–P4 y Oz–P4), coherente
con el pico α de la PSD. Fíjate en la celda O1–Oz, casi 0: son electrodos vecinos que ven la
misma fuente con desfase nulo, y wPLI la descarta. Es el comportamiento buscado, no ausencia
de actividad. La figura no dice qué aristas son significativas (capítulo 14) ni si el patrón
es patológico: no hay grupo de comparación.]

=== 6 · Interpretación y validación
*Por qué Hilbert.* No hay en el repositorio una comparación que lo justifique frente a las
alternativas: es el «método original» de NeuroMIND, conservado por continuidad con `EEG_Julia`
y con los resultados de julio. El propio código lo describe como «válido para épocas
suficientemente largas, sensible a artefactos de borde» y señala el multitaper como el mejor
compromiso entre varianza y sesgo. Es una decisión de reproducibilidad, no de optimalidad.

*Épocas de 1 s.* El filtro de banda se aplica *dentro* de cada época. Un ciclo de 0,5 Hz dura
2 s: con 1 s, DELTA tiene entre 0,5 y 4 ciclos y THETA entre 4 y 8 (el código avisa por debajo
de 4). La respuesta al impulso de un Butterworth de orden 8 con borde a 0,5 Hz dura varios
segundos, más que la época: en DELTA buena parte de la señal analítica es transitorio de borde.
En la matriz DELTA (`wpli_DELTA.png`) la arista máxima (C4–FC1, 0,563) aparece aislada y
toca a C4, el canal del pico lento anómalo: sospechosa de artefacto. *Recomendación:* no
interpretar DELTA, o activar `exclude_unreliable_bands = true`.

#caja(azul2, azulf, "Protocolo de sensibilidad (pendiente).")[(a) Duración de época: perfil de
segmentación `default` con $L$ = 1, 2 y 4 s, sin solape, misma ICA, y comparar media y máximo
de wPLI por banda y número de épocas válidas; si DELTA/THETA cambian mucho con $L$ y ALPHA no,
las bandas lentas dependen del borde. (b) Estimador: `hilbert`, `fourier_csd` y `multitaper`
sobre las mismas épocas, con la correlación de aristas frente a Hilbert por banda. El script
está preparado y no se ha lanzado.]

// ═════════════════════════════════════════════════════════════
= Capítulo 14 · Datos subrogados y estadística

#cabecera("pend",
  [`[7/8]`, bloque opcional (solo si `[surrogates] enabled = true`).],
  [Cap. 14 (previsto).],
  [`surrogate_test` (desplazamiento circular), FDR de Benjamini–Hochberg (`fdr_correction`).],
  [Visor `plot_surrogate.jl` (:8775).],
)

=== 1 · Objetivo científico
Decidir, arista por arista, si la wPLI observada es mayor de lo que daría el azar con las
mismas señales pero sin relación de fase entre canales, controlando los falsos positivos de
435 contrastes por banda.

=== 2 · Entrada y salida
Épocas + matrices observadas → por banda: p, q, máscara de significación, media y σ nulas y
z. Ficheros en `tables/surrogate/`, `json/surrogate_summary.json`,
`figures/surrogate/surrogate_null_{BANDA}.png` y `cache/surrogate/` (lo lee el visor).

=== 3 · Rutina
```julia
sr = surrogate_test(ep, conn, "ALPHA", cfg)     # N, α y semilla de cfg.surrogates
sr.p_values, sr.sig_mask, sr.fdr_threshold
```
```bash
julia --project=. scripts/run_single_subject.jl \
      --config ~/Proyectos/NeuroMIND-results/config/M05_surrogates.toml
```

=== 4 · Fundamento
*Nula.* En cada repetición, cada canal de cada época se desplaza circularmente un número
aleatorio de muestras entre $L\/5$ y $4L\/5$: se conservan el espectro y la autocorrelación de
cada canal, pero se destruye la relación de fase entre canales. Se recalcula wPLI con *el mismo
estimador*. *p de Monte Carlo:*
$p_(i j) = (\#{W_(i j)^"nula" >= W_(i j)^"obs"} + 1) \/ (N + 1)$, que nunca es 0 y tiene un
mínimo de $1\/(N+1)$.

*FDR (Benjamini–Hochberg, por banda).* Con los $m = 435$ p ordenados, se busca el mayor $k$ con
$p_((k)) <= k alpha \/ m$ y se declaran significativas las $k$ aristas de menor p. La familia es
cada banda por separado; una corrección conjunta de las 7 bandas sería más estricta.

=== 5 · Resultado de M05
#tabla(
  (3.4cm, 1fr),
  ([Corrida], [Resultado]),
  [Prueba, N = 20 (10-11-2026)], [0 aristas significativas en las 7 bandas; ≈ 5 s por banda.],
  [Archivo de julio, N = 500 (solo referencia)], [ALPHA 246/435 (56,6 %), GAMMA 174 (40,0 %),
    THETA 39, BETA_HIGH 29, resto 0. Umbral FDR en ALPHA: p ≤ 0,028.],
  [Estimación con N = 500 hoy], [≈ 2 min por banda, ≈ 15 min en total.],
)

=== 6 · Interpretación y validación
*Por qué N = 20 no da nada.* El p mínimo es $1\/21 approx 0,048$. Para que BH declare
significativa alguna arista con ese p hace falta $k >= 0,048 dot 435 \/ 0,05 approx 414$: al
menos 415 de las 435 aristas tendrían que tener la nula siempre por debajo de lo observado. No
es imposible, pero exige un efecto casi total; en la práctica, con N = 20 y 435 contrastes por
banda el procedimiento no tiene potencia. Con N = 500 el p mínimo es ≈ 0,002 y bastan unas 18
aristas con ese p para superar el umbral.

*Leer los resultados de julio con cautela.* (1) Proceden de otra ICA. (2) Un 40 % de aristas
GAMMA significativas no indica conectividad gamma cerebral: la actividad muscular residual
también produce relaciones de fase consistentes (IC15 e IC31 sugieren músculo lateral). (3) La
significación es intra-sujeto: dice que una arista no es azar en *este* registro, no que difiera
de un control.

*Estado.* Pendiente hasta lanzar la corrida con N = 500 sobre el código actual.

// ═════════════════════════════════════════════════════════════
= Capítulos 15–17 · Grupo y conclusiones

#cabecera("no",
  [Post-hoc, sobre la cohorte procesada por lotes.],
  [Caps. 15 (transversal) y 16 (longitudinal), previstos.],
  [`Transversal.jl`, `Longitudinal.jl`.],
  [Visores `plot_transversal.jl` (:8781) y `plot_longitudinal.jl` (:8780).],
)

No se pueden ejecutar con un sujeto. Con M05 solo se puede mirar su par T1→T2 procesando sus
cuatro grabaciones:

```bash
julia --project=. scripts/run_batch_pipeline.jl --subjects M05 --dry-run
julia --project=. scripts/run_batch_pipeline.jl --subjects M05
```

Orden para la cohorte (referencia):

```bash
julia --project=. scripts/run_batch_pipeline.jl --skip-done
julia --project=. scripts/run_transversal_analysis.jl
julia --project=. scripts/run_longitudinal_analysis.jl
```

El capítulo 17 (conclusiones) no tiene rutinas.

// ═════════════════════════════════════════════════════════════
= Cierre · comprobaciones y validación cruzada

== Lista de comprobación tras el modo A

```bash
ls $M05/tables/overview.csv $M05/config_snapshot.toml $M05/pipeline_log.txt
grep -E 'exclude_fp2|use_dwpli|n_channels_used|wpli_method' $M05/config_snapshot.toml
# 31 columnas = 1 de nombres + 30 canales
head -1 $M05/tables/connectivity/wpli_ALPHA.csv | awk -F, '{print NF}'
# 435 aristas por banda
grep -c ',ALPHA,' $M05/tables/connectivity/connectivity_edges.csv
grep '^M05,' results/qc/qc_decision_table.csv
grep 'Convergencia' $M05/pipeline_log.txt
```

== Regresión frente al archivo de julio

```bash
diff <(cut -d, -f1,6 "$ARCHIVO/tables/channel_statistics.csv") \
     <(cut -d, -f1,6 $M05/tables/channel_statistics.csv) && echo "z-RMS idéntico"
diff "$ARCHIVO/tables/rejected_segments.csv" $M05/tables/rejected_segments.csv
```

#tabla(
  (4.8cm, 1fr),
  ([Comparación (prueba del 10-11-2026)], [Resultado frente a \$ARCHIVO]),
  [`raw_signal.csv`, `filtered_signal_*.csv`], [Iguales salvo redondeo (|Δ| ≤ 10⁻⁴ µV).],
  [`channel_statistics.csv`], [z-RMS idéntico; Fp2 único canal malo.],
  [`qc_decision_table.csv`], [Misma fila: `include`, 30 canales, 99/100 épocas.],
  [ICA], [Distinta: julio leyó una ICA antigua de caché.],
  [`rejected_segments.csv`], [Misma época y canal; amplitud 72,18 frente a 72,56 µV.],
  [PSD y wPLI], [Diferencias pequeñas heredadas de la ICA.],
)

== Validación cruzada Julia ↔ MNE

#cabecera("no",
  [Réplica independiente de las fases 3–7 en MNE-Python y comparación (fase 8).],
  [Dashboard, panel 15; `validation/` del libro.],
  [`run_full_pipeline.py`, `run_phase{3..7}_m05.py`, `run_phase8_validation.py`,
    `compare_neuromind_mne_rejections.py`.],
  [Informe en `mne_brain/validation/reports/`.],
)

```bash
cd mne_brain
uv run python scripts/run_full_pipeline.py --subject M05 --session T2 --condition EC
uv run python scripts/run_phase8_validation.py
cd ..
```

#fallo[`run_phase8_validation.py` busca las tablas de Julia en rutas planas
(`…/eyesclosed/wpli_{BANDA}.csv`, `psd_by_channel.csv`, `band_power_summary.csv`), pero el
pipeline las escribe en `tables/` y `tables/connectivity/`: fallará con `FileNotFoundError`.
Además escribe su informe dentro del repo (`mne_brain/validation/reports/`, ignorado), no en
`NeuroMIND-results`. Mientras no se corrija, *no hay validación independiente* de PSD ni de
wPLI para el código actual: el pipeline es reproducible, pero no está validado frente a MNE.]

== Dashboard y libro

`scripts/launch_dashboard.jl` (http://localhost:8080) solo lee `results/`; los paneles de señal
cruda y PSD siguen buscando `*_eeg_data.tsv` y pueden salir vacíos. El libro registra la
ejecución `sub-M05_ses-T2_eyesclosed_2026-07-25` con el sha256 de cada figura y tabla citada
(`metadata/figuras-resultados.json`): tras la corrida oficial hay que registrar la nueva y
revisar el texto que las describe. Cómo ver cada figura de M05 que cita el libro: Anexo C.

== Antes de dar por definitivo este cuaderno

#tabla(
  (0.6cm, 1fr, 3.2cm),
  ([], [Cuestión], [Estado]),
  [1], [Estabilidad de la ICA (protocolo del capítulo 8). Topomapas: corregidos el 10-11-2026.],
    estado("pend"),
  [2], [Validación Julia ↔ MNE: actualizar rutas de la fase 8 y ejecutarla.], estado("no"),
  [3], [Justificación de los parámetros de conectividad: sensibilidad a la duración de época y
    al estimador (capítulo 13); decisión sobre DELTA.], estado("pend"),
  [4], [Corrida oficial en `NeuroMIND-results` con N = 500 y regeneración de las figuras.],
    estado("no"),
  [5], [Etiquetas de unidades (µV²/Hz) en tablas y figuras espectrales.], estado("adv"),
)

// ═════════════════════════════════════════════════════════════
= Anexo A · Diccionario de variables

#{
  set text(7.9pt)
  tabla(
    (3.4cm, 1fr, 1.3cm, 3.8cm),
    ([Variable (fichero)], [Definición y fórmula], [Unidades], [Significado y precaución]),
    [`rms_uv` (channel_statistics)], [$sqrt(1/N sum x^2)$ sobre la señal cruda.], [µV],
      [Amplitud global del canal.],
    [`rms_zscore`], [$(R M S_c - overline(R M S)) \/ s_(R M S)$ entre canales.], [—],
      [> 3 en valor absoluto ⇒ `is_bad`. Relativo al registro.],
    [`hfnoise`], [$P_(25-45) \/ P_(1-45)$ (Welch sobre la cruda).], [—],
      [Fracción de alta frecuencia; alto sugiere músculo o mal contacto.],
    [`snr_db`], [$10 log_10 (P_"ALPHA" \/ P_"GAMMA")$.], [dB],
      [«Señal» = α; en ojos abiertos pierde sentido.],
    [`mean_abs_corr`], [Media de |r| con el resto de canales.], [—],
      [Muy baja: canal desconectado; muy alta: puente.],
    [$overline(sigma)$ (`raw_sigma_mean_uv`)], [Media entre canales de la σ cruda.], [µV],
      [> 20 µV ⇒ `amplitude_warning`.],
    [`variance_pct` (ica_components)], [Varianza explicada por el componente.], [%],
      [Coste de rechazarlo.],
    [`frontal_ratio`, `temporal_ratio`], [Media de |A| en canales F\*/T\* entre la del resto.], [—],
      [FT9/FT10 cuentan como frontales (regex `^F`).],
    [`blink_ratio`], [$P_(0,5-4) \/ P_(4-40)$ de la activación.], [—], [Potencia lenta.],
    [`emg_ratio`], [$P_(30-80) \/ P_(1-30)$.], [—], [Potencia rápida (músculo).],
    [`line_ratio`], [$P_(48-52) \/ P_(0,5-100)$.], [—], [Tras el notch suele ser < 0,5 %.],
    [`kurtosis`, `extreme_frac`], [Curtosis (no en exceso; normal = 3) y fracción de muestras
      a más de 5σ.], [—], [Transitorios.],
    [`*_score`, `artifact_score`], [Combinaciones de z-scores entre componentes (cap. 8).],
      [—], [> 1,5 ⇒ artefacto. Relativas.],
    [`quality` (segments_table)], [$max(0, 1 - 0,5 dot max abs(x) \/ 70)$.], [—],
      [Descriptiva; no decide el rechazo.],
    [`max_amp_uv`, `p2p_uv`], [Máximo y pico a pico de la época antes del rechazo.], [µV],
      [p2p solo es diagnóstico.],
    [`coverage_pct`], [% de épocas válidas por canal.], [%], [],
    [`power_uv2` (psd_by_channel)], [$P_c (f_k)$, media sobre épocas (cap. 12).], [*µV²/Hz*],
      [Mal rotulada como µV².],
    [`band_power_` `summary`], [Media de los bins de la banda.], [*µV²/Hz*],
      [No es potencia integrada (× nº bins × Δf).],
    [`spectral_indices`], [Cocientes α/θ, β/α, θ/β, γ/α de valores de banda; pico α (Hz).], [—],
      [Adimensionales.],
    [`wpli` (`wpli_*`, edges)], [$abs(sum I) \/ sum abs(I)$, $I = "Im"(z_i z_j^ast)$.], [—, 0–1],
      [Consistencia del signo del desfase.],
    [`rank` (connectivity_edges)], [Orden por wPLI descendente dentro de la banda.], [—], [],
    [`strength` (network_metrics)], [Suma de wPLI del canal, promediada sobre las 7 bandas.],
      [—], [Mezcla bandas.],
    [`degree`], [Nº de aristas con wPLI > 0,1, promediado sobre bandas.], [—],
      [Usa `[graph] density = 0,1` como umbral *absoluto*, no como densidad.],
    [`p`, `q` (surrogate)], [p de Monte Carlo (+1); q de BH por banda.], [—],
      [p mínimo = 1/(N+1).],
    [`null_mean`, `null_std`, `z_score`], [Media y σ de la nula; $(W - mu_0) \/ sigma_0$.],
      [—], [z orientativo; la decisión usa p.],
  )
}

// ═════════════════════════════════════════════════════════════
= Anexo B · Resumen capítulo → paso → salida

#tabla(
  (1.1cm, 1.3cm, 1fr, 1fr, 3.5cm),
  ([Cap.], [Paso], [Rutina principal], [Salida clave (bajo \$M05/)], [Estado]),
  [4], [Fase 0], [`build_bids_full.jl`], [`data/bids/sub-M05/`], estado("ok"),
  [5], [1–2/8], [`load_single_subject`, `flag_bad_channels`], [`tables/raw_signal.csv`,
    `tables/channel_statistics.csv`], estado("ok"),
  [6], [3/8], [`filter_recording`], [`tables/filtered_signal_*.csv`], estado("adv"),
  [8], [4/8], [`run_ica`, `apply_ica_rejection`], [`tables/ica/`,
    `json/ica_summary.json`], estado("adv"),
  [9–11], [5/8], [`segment_recording`, `apply_baseline`, `reject_artifacts`],
    [`tables/segments_table.csv`, `tables/rejected_segments.csv`], estado("ok"),
  [12], [6/8], [`compute_psd`], [`tables/psd_by_channel.csv`], estado("adv"),
  [13], [7/8], [`compute_wpli`], [`tables/connectivity/wpli_*.csv`], estado("adv"),
  [14], [7/8], [`surrogate_test`], [`tables/surrogate/`], estado("pend"),
  [15–16], [grupo], [`Transversal.jl`, `Longitudinal.jl`], [`results/{transversal,
    longitudinal}/`], estado("no"),
)

// ═════════════════════════════════════════════════════════════
= Anexo C · Figuras del libro: cómo verlas

Para cada figura de M05 que cita el libro (`NeuroMIND-quarto`), cómo verla desde este
repositorio. La lista sigue la «Procedencia de las figuras» del libro, capítulos 5–13. Hay tres
vías:

#tabla(
  (3.2cm, 1fr),
  ([Vía], [Qué es]),
  [Modo A], [La figura es una salida del pipeline: tras la corrida, `open $M05/figures/…`. Es
    la única vía cuyas figuras puede enlazar el libro (con su `pipeline_log.txt` y
    `config_snapshot.toml`).],
  [Visor], [`julia --project=. src/interactive/aux/sub-M05_ses-T2_eyesclosed/<visor>.jl` y
    abrir `http://127.0.0.1:<puerto>/`. Lee las tablas de `$M05`, así que necesita antes el
    Modo A. «Guardar PNG» escribe en la carpeta del visor (ignorada por git): sirve para mirar,
    pero el libro *no* enlaza esas figuras.],
  [REPL], [Una rutina del paquete sobre las variables del Modo B. Tampoco es enlazable.],
)

#aviso[Si el libro está servido en local en el puerto 8765, `plot_raw.jl` no arranca (usa el
mismo puerto): detener ese servidor antes. Los visores marcados como aproximación muestran los
mismos datos, pero no reproducen la figura del informe.]

#nota[«No reproducida» quiere decir que el libro describe la figura del informe en un recuadro,
sin mostrarla, porque no es una salida trazable del pipeline. «Informe» quiere decir que el
capítulo todavía no se ha migrado al libro. El libro solo enlazará las figuras de la vía Modo A.]

#let figuras-cap(titulo, ..filas) = {
  block(sticky: true, above: 1.1em, below: 0.45em,
        text(9.4pt, fill: azul, weight: "bold", titulo))
  set text(7.6pt)
  tabla((3.4cm, 2.7cm, 1fr), ([Figura del libro (M05)], [En el libro], [Cómo verla]), ..filas)
}

#figuras-cap([Capítulo 5 · Señal raw y control de calidad],
    [Fp2 frente a Cz y Cz frente a F3, registro completo], [No reproducida.], [Script del informe: `walkthrough_M05.jl`, no versionado. Visor
      `plot_raw.jl` (:8765): canales Fp2 y Cz
      (o Cz y F3), ventana completa. Aproximación.],
    [Trazas apiladas, 31 canales, 0–10 s], [No reproducida.], [Visor `plot_raw_butterfly.jl`
      (:8766): todos los canales, 0–10 s.],
    [PSD media de la señal sin filtrar], [No reproducida.], [Sin vía. El libro la atribuye a
      `plot_raw_psd_average` (`src/visualization/qc_figures.jl`), que no aparece en ningún
      commit de NeuroMIND; `plot_spectral*.jl` muestran la PSD del paso `[6/8]`, no la cruda.],
    [Histogramas de amplitud por canal (láminas 1–16 y 17–31)], [No reproducida.],
      [Visor `plot_raw_histogram.jl` (:8768): canales 1–16 y después 17–31.],
    [Topografía de la desviación típica por canal], [No reproducida.], [REPL (después del
      bloque de la § Capítulo 5):
      ```julia
      using CairoMakie
      fig = plot_topomap(Float64.(qc.std_uv),
                         rec.meta.channel_names,
                         rec.meta.channel_positions;
                         style = :power, colorbar_label = "µV")
      f = joinpath(tempdir(), "sigma_M05.png")
      save(f, fig); run(`open $f`)
      ```
      En Cursor, `display(fig)` la muestra en el panel de gráficos.
      Sin el arreglo de la rama `fix/topomap-transpuesto` sale transpuesta (el máximo de Fp2
      aparece junto a FT10/T8). Con la versión de esa rama, el máximo cae en Fp2 (probado el
      10-11-2026).],
    [Panel de métricas QC por canal (5 paneles)], [No reproducida.], [Sin vía gráfica. Los
      valores son los de `$M05/tables/channel_statistics.csv`, que el libro muestra como
      tabla.]
)

#figuras-cap([Capítulo 6 · Filtrado],
    [Vista previa cruda frente a filtrada, 5 canales, 10 s], [*Enlazada* (ejecución del 25-07-2026).],
      [Modo A: `open $M05/figures/filtered_signal_preview.png`.],
    [Cz y F3, cruda frente a filtrada, registro completo], [No reproducida.],
      [Script del informe: `walkthrough_M05.jl`. Visor `plot_filtered_vs_raw.jl` (:8769): canal Cz (o F3), etapas crudo y `lowpass`.
      Aproximación.],
    [PSD cruda frente a filtrada], [No reproducida.], [Sin vía (script del informe: `walkthrough_M05.jl`).]
)

#figuras-cap([Capítulo 8 · ICA (previsto)],
    [Topografías de IC5, IC7, IC9 e IC13], [Informe.], [Modo A:
      `open $M05/figures/ica/ica_topomap_{005,007,009,013}.png`. Transpuestas sin la rama
      `fix/topomap-transpuesto`. El número de componente no se conserva entre corridas.],
    [Puntuaciones de artefacto por componente], [Informe.], [Sin figura en el código. Visor
      `plot_ica_components.jl` (:8770) muestra los rasgos; valores en
      `$M05/tables/ica/ica_component_features.csv`.],
    [Varianza por canal antes y después de la ICA], [Informe.], [Sin vía.],
    [IC13: curso temporal y espectro], [Informe.], [Visor `plot_ica_components.jl` (:8770),
      componente 13 activo: activación y PSD de los primeros 10 s exportados, no del registro
      completo. Aproximación.]
)

#figuras-cap([Capítulos 9–11 · Segmentación, línea base y rechazo (previstos)],
    [Histograma de calidad, épocas superpuestas en Oz y época 12 rechazada], [Informe.],
      [Visor `plot_epochs.jl` (:8772): regenera las tres. Necesita `cache/ica_result.jls`,
      `ica_labels_auto.csv`, `segments_table.csv` y `rejected_segments.csv`.],
    [Antes y después de la línea base, C4], [Informe (época 89).], [Visor `plot_baseline.jl`
      (:8773): por defecto C4 y la época con mayor desplazamiento medio en la ventana de
      línea base.],
    [Histograma de pico a pico], [Informe.], [Sin figura en el código. Los 20 intervalos están
      en `$M05/json/artifact_rejection_summary.json` (`p2p_histogram`) y cada valor en la
      columna `p2p_uv` de `$M05/tables/segments_table.csv`.]
)

#figuras-cap([Capítulo 12 · Análisis espectral (previsto)],
    [PSD de cada canal], [Informe.], [Modo A: `open $M05/figures/psd_all_channels.png`. Visor
      `plot_spectral_PSD.jl` (:8767).],
    [Potencia por banda y canal], [Informe.], [Modo A:
      `open $M05/figures/band_power_summary.png` (eje «µV²»: son µV²/Hz, § Capítulo 12). Visor
      `plot_spectral.jl` (:8767, no a la vez que el anterior).],
    [Potencia por banda (TikZ del informe)], [Informe.], [Visor `plot_spectral.jl`, modo
      `bands_abs`, con los datos de `$M05/tables/band_power_summary.csv`. Aproximación.],
    [Topografía de potencia por banda], [Informe.], [Sin vía: `band_topomap_grid` no aparece
      en el código.]
)

#figuras-cap([Capítulo 13 · Conectividad (previsto)],
    [Matrices wPLI ALPHA y DELTA], [Informe.], [Modo A:
      `open $M05/figures/connectivity/wpli_{ALPHA,DELTA}.png`. Visor `plot_connectivity.jl`
      (:8774).]
)

