# Technical Design Document — Sistema de ejecución de rutinas

**Proyecto:** NeuroMIND · **Autor:** Rafael Castro Triguero · **Fecha:** 2026-07-22
**Estado:** propuesta de diseño (sin implementar) · **Destino:** convertir en tareas para Claude Code

---

## 0. Resumen ejecutivo

NeuroMIND ejecuta hoy sus rutinas (`run_single_subject`, `run_batch_pipeline`,
`run_transversal_analysis`, `run_longitudinal_analysis`) con `println` sueltos,
sin limpieza determinista de resultados, sin manifiesto de ejecución y sin un
sistema de logs persistente y estructurado. Funciona, pero no es trazable ni
robusto frente a cambios de código.

Este documento propone una capa de **runtime** (ejecución) transversal a todas
las rutinas, con seis responsabilidades separadas: gestión de resultados,
manifiesto de ejecución, logging, consola, progreso y resumen final.

**Decisión de dependencias — clave:** las cuatro librerías maduras necesarias
(`ProgressMeter`, `Crayons`, `PrettyTables`, `LoggingExtras`) **ya están en el
`Manifest.toml`** del proyecto. No se añade ninguna dependencia nueva. Esto
elimina el principal riesgo de adopción.

**Recomendaciones de una línea:**

| Tema | Decisión |
|------|----------|
| Limpieza de resultados | **Clean-slate por unidad hoja** (borrar la carpeta del sujeto/análisis antes de regenerarla) + **manifiesto** de ficheros producidos |
| Metadatos de ejecución | **`run_manifest.json`** (máquina) + conservar `pipeline_log.txt` (humano) y `config_snapshot.toml` (parámetros) |
| Logging | **`LoggingExtras`** con `TeeLogger` (consola + fichero), niveles estándar |
| Progreso | **`ProgressMeter`** |
| Colores | **`Crayons`** con 4 niveles semánticos (INFO/WARN/ERROR/SUCCESS) |
| Tablas en terminal | **`PrettyTables`** |

---

## 1. Gestión de la carpeta de resultados

### El problema central

Si en una versión nueva del código cambia el nombre de una figura o tabla
(p. ej. `wpli_matrix_ALPHA_EC.csv` → `wpli_ALPHA.csv`), la estrategia actual
—escribir/sobrescribir sin limpiar— deja el fichero antiguo en la carpeta.
Resultado: la carpeta contiene una mezcla de dos versiones, un análisis
posterior o un lector de dashboard puede tomar el fichero viejo, y nadie se
entera. Es exactamente lo que ya ocurrió en este proyecto con el árbol heredado
y los sufijos `_EC`. **Cualquier diseño que no resuelva esto es inaceptable.**

### Alternativas

| Estrategia | Ventajas | Inconvenientes |
|-----------|----------|----------------|
| **A. Borrar todo `results/` antes de cada lanzamiento** | Nunca hay ficheros obsoletos | Catastrófico en batch: procesar un sujeto borraría los otros 204. Rompe análisis de grupo que dependen de resultados previos. Inaceptable |
| **B. Sobrescribir solo los ficheros existentes** | Simple, rápido | **Es el bug actual**: un fichero renombrado deja huérfano al antiguo. No detecta obsoletos |
| **C. Mantener resultados antiguos** (nunca borrar) | Cero pérdida | La carpeta se llena de versiones mezcladas; imposible saber qué es actual. Empeora el problema |
| **D. Carpeta nueva por ejecución** (`results/run_2026-07-22_1030/`) | Aislamiento total, histórico completo | Rompe rutas estables (el análisis de grupo espera `results/subjects/sub-M05/…` fijo). Explosión de disco (GB por corrida). Complica encontrar "el resultado actual" |
| **E. Clean-slate por unidad hoja + manifiesto** ⭐ | Cada unidad (sujeto/sesión/tarea, o análisis de grupo) borra **solo su propia carpeta** antes de regenerarla → imposible dejar huérfanos. No toca a las demás. Rutas estables. Manifiesto lista lo producido | Requiere disciplina: definir bien la "unidad" y que el borrado sea atómico |

### Recomendación: E — clean-slate por unidad hoja

Es la práctica madura en pipelines científicos reproducibles (equivale al patrón
de Snakemake/Nextflow: cada *rule/process* posee su directorio de salida y lo
regenera por completo). Concreta para NeuroMIND:

- **Unidad del pipeline individual:** `results/subjects/sub-{ID}/ses-{SES}/{task}/`.
  Antes de escribir, se **vacía esa carpeta** (no `results/subjects/` entero).
  Reprocesar M05/EC no puede dejar ni un fichero de la corrida anterior, y no
  toca a M07 ni al resto.
- **Unidad de grupo:** `results/transversal/{EC|EO}/` y `results/longitudinal/{EC|EO}/`.
  Mismo criterio: vaciar la carpeta de la condición antes de regenerarla.
- **Índices globales** (`subjects_index.csv`, `qc/qc_decision_table.csv`): NO se
  borran; se actualizan por clave (upsert de la fila del sujeto). Son acumulativos
  por diseño.

**Procedimiento de escritura seguro (por unidad):**

1. Escribir en un **directorio de staging** hermano y oculto: `…/{task}.tmp/`.
2. Al completar sin error, **borrar** el directorio final y **renombrar**
   (`mv`) el staging al nombre final. `mv` dentro del mismo sistema de ficheros
   es atómico: o está la versión nueva completa, o la vieja intacta. Nunca un
   estado a medias si el proceso muere a mitad.
3. Escribir el **manifiesto** (`run_manifest.json`, §2) como último paso, con la
   lista exacta de ficheros producidos.

Esta combinación (staging atómico + clean-slate + manifiesto) da las tres
garantías: sin huérfanos, sin corrupción por interrupción, y trazabilidad de qué
se generó.

**Salvaguarda para la unidad "borrar":** el borrado solo puede apuntar a rutas
bajo `results/subjects/`, `results/transversal/` o `results/longitudinal/`
verificadas con `startswith(realpath(...))`. Nunca acepta una ruta arbitraria.
Una función `_safe_clean(dir)` que aborta si `dir` no está dentro de `results/`
previene un `rm -rf` accidental por un `subject_id` mal formado.

---

## 2. Información de la ejecución (`run_manifest.json`)

### Formato: JSON, no TXT

| Criterio | TXT | JSON ⭐ | TOML |
|----------|-----|--------|------|
| Legible por humano | ✅ | ✅ (indentado) | ✅ |
| Parseable por máquina | ❌ | ✅ | ✅ |
| Consultable (`jq`, filtros, agregación batch) | ❌ | ✅ | parcial |
| Anidamiento (parámetros, listas de ficheros) | ❌ | ✅ | limitado |
| Ecosistema Julia (`JSON3` stdlib-adyacente) | — | ✅ | ✅ |

**Decisión:** `run_manifest.json` por unidad de ejecución. Razón: se quiere poder
**agregar** los manifiestos de las 205 grabaciones (p. ej. "¿cuántas generaron
menos figuras de lo esperado?", "¿cuál tardó más?") — eso exige máquina, y `jq`
sobre JSON es el estándar. TXT no escala a análisis batch.

**No se reinventa lo que ya existe:** `config_snapshot.toml` (parámetros exactos)
y `pipeline_log.txt` (traza cronológica humana) se conservan. El manifiesto los
**complementa**, no los sustituye — es el índice de la corrida.

### Contenido de `run_manifest.json`

```
{
  "schema_version": 1,
  "project":     { "name": "NeuroMIND", "version": "0.2.0" },
  "run": {
    "unit":        "sub-M05/ses-T2/eyesclosed",
    "started_at":  "2026-07-22T10:30:00",
    "finished_at": "2026-07-22T10:38:42",
    "duration_s":  521.7,
    "status":      "success"          // success | warning | error
  },
  "environment": {
    "julia_version": "1.11.2",
    "os":            "macOS 14.5 (arm64)",
    "hostname":      "…", "user": "rafa",
    "git_commit":    "04dc8a9", "git_branch": "main", "git_dirty": false
  },
  "config": {
    "path":   "config/pipeline.toml",
    "sha256": "…",                     // hash del config efectivo
    "seed":   1234                     // seed de ICA / surrogates
  },
  "inputs":  { "n_subjects": 1, "vhdr": "…/sub-M05_…​.vhdr", "n_channels": 31, "fs_hz": 500.0 },
  "outputs": {
    "n_files": 98, "n_tables": 69, "n_figures": 21, "n_json": 6,
    "files": [ "overview.csv", "wpli_ALPHA.csv", "figures/psd_all_channels.png", … ]
  },
  "warnings": [ "DELTA: 0.5 ciclos/época < mínimo 4.0" ],
  "errors":   []
}
```

`git_dirty` (árbol de trabajo con cambios sin commitear) y `config.sha256` son
los dos campos que más protegen la reproducibilidad: detectan "este resultado se
generó con código/config no versionado o modificado".

---

## 3. Organización del árbol de resultados

La estructura vigente (subjects/transversal/longitudinal/qc/logs) es correcta.
El diseño la **completa** por unidad hoja, sin cambiar la raíz:

```
results/
├── README.md
├── subjects_index.csv              ← índice global (upsert)
├── qc/
│   └── qc_decision_table.csv       ← índice global de decisiones QC
├── logs/
│   ├── batch_run_{timestamp}.csv   ← resumen de cada corrida en lote
│   └── run_{timestamp}.log         ← log persistente (§5)
├── subjects/
│   └── sub-{ID}/ses-{SES}/{task}/  ← UNIDAD hoja (clean-slate, ver §3.1)
├── transversal/{EC|EO}/            ← UNIDAD hoja (clean-slate)
│   └── run_manifest.json + tablas + figures/
└── longitudinal/{EC|EO}/           ← UNIDAD hoja (clean-slate)
    └── run_manifest.json + tablas + figures/
```

**Creación automática:** solo la unidad hoja de la corrida actual y `logs/`,
`qc/`. No se pre-crean carpetas vacías de sujetos que no se procesan (evita el
ruido de directorios vacíos que generaba `ensure_dirs`). `cache/` se excluye
del clean-slate opcionalmente (para no invalidar la caché de ICA en
re-corridas de pasos posteriores) — decisión configurable.

**Nota sobre `temp/` y `exports/`:** `temp/` se materializa como el staging
atómico `{unit}.tmp/` y se autoelimina; no queda como carpeta permanente.
`exports/` (informes agregados, PDF) tiene sentido a nivel raíz solo si se genera
un entregable global; se propone `results/exports/` creado bajo demanda por el
script de informe, no por el pipeline.

### 3.1 Estructura interna de la unidad hoja (sujeto)

Hoy `{task}/` mezcla ~45 CSV, varios JSON, un TOML y un log sueltos junto a
`figures/` — difícil de navegar. Se organiza por **tipo de salida**, con
subcarpeta por fase **solo donde el volumen lo justifica**:

```
sub-{ID}/ses-{SES}/{task}/
├── run_manifest.json              ← §2
├── config_snapshot.toml           ← parámetros exactos
├── pipeline_log.txt               ← traza humana
├── tables/
│   ├── overview.csv
│   ├── qc_summary.csv
│   ├── channel_statistics.csv
│   ├── band_power_summary.csv
│   ├── raw_signal.csv
│   ├── ica_signal_before.csv
│   ├── ica_signal_after.csv
│   └── connectivity/              ← única subcarpeta por fase (~35 ficheros)
│       ├── connectivity_edges.csv
│       ├── significant_connections.csv
│       ├── wpli_{band}.csv                  (×7)
│       ├── wpli_pvalues_{band}.csv          (×7)
│       ├── wpli_qvalues_{band}.csv          (×7)
│       ├── wpli_significant_{band}.csv      (×7)
│       └── surrogate_null_stats_{band}.csv  (×7)
├── figures/
│   └── connectivity/              ← mismo criterio si las figuras de wPLI proliferan por banda
├── json/
│   ├── ica_summary.json
│   └── surrogate_summary.json
└── cache/
    └── ica_result.jls             ← excluido del manifiesto y del clean-slate
```

**Regla de decisión (repetible, no ad-hoc):** tres categorías fijas en la raíz
de la unidad — `tables/`, `figures/`, `json/` — más los metadatos de la
corrida sueltos en la raíz (manifiesto, config, log), porque describen la
corrida entera, no un tipo de salida. Dentro de `tables/`/`figures/`, solo se
crea subcarpeta por fase cuando el volumen de ficheros lo justifica: hoy eso
es únicamente `connectivity/` (wPLI × 4 variantes × 7 bandas + surrogates).
QC, espectral e ICA generan 1–3 ficheros cada uno — meterlos en subcarpetas
de fase añadiría profundidad sin reducir el ruido real.

**Alternativa descartada:** una subcarpeta por cada una de las 8 fases del
pipeline, incluso con un solo fichero dentro (`qc/qc_summary.csv`,
`spectral/band_power_summary.csv`). Se descarta por generar carpetas de un
solo fichero que no ayudan a navegar.

**Alcance:** este desglose aplica a la unidad hoja de sujeto
(`subjects/sub-{ID}/ses-{SES}/{task}/`). Las unidades de grupo
(`transversal/`, `longitudinal/`) generan muchos menos ficheros y se
mantienen planas por ahora; si su volumen crece, se les aplica el mismo
criterio.

**Dependencias a actualizar (fuera de esta fase):** el escritor
(`_save_all_results` en `SingleSubjectPipeline.jl`) y las rutas de lectura del
dashboard (`src/webapp/App.jl`, `web/views/dashboard.html`) — se abordan en
una fase posterior, ya acordada como no bloqueante para adoptar esta
estructura.

---

## 4. Lanzamiento de rutinas desde terminal

### 4.1 Encabezado

Un banner al inicio, idéntico para las cuatro rutinas, generado por el módulo de
consola. Muestra: proyecto + versión, fecha/hora, versión de Julia, SO,
directorio de trabajo, config en uso, unidad(es) a procesar. Es el equivalente
en pantalla de la cabecera del `run_manifest.json`.

### 4.2 Comprobaciones previas (preflight)

Antes de tocar datos, una batería de checks con salida `OK`/`ERROR` y **fallo
temprano** (abortar si algo crítico falla, en vez de morir a mitad del paso 6):

| Check | Criticidad |
|-------|-----------|
| `config/pipeline.toml` existe y parsea | crítico |
| `data/bids/` y el `.vhdr`/metadata del sujeto existen | crítico |
| Claves obligatorias del config presentes (validación de esquema) | crítico |
| Permisos de escritura en `results/` | crítico |
| Espacio libre en disco > umbral (p. ej. 2 GB) | aviso |
| `electrodes.tsv` presente (topomaps) | aviso |
| Git limpio (sin cambios sin commitear) | aviso |

El preflight imprime una tabla resumen y, si hay un `ERROR` crítico, se detiene
con código de salida ≠ 0 (importante para CI y para encadenar en shell).

### 4.3 Durante la ejecución

Mostrar por paso: etapa actual (`[4/8] ICA`), sujeto en curso, tiempo del paso.
Para el **lote** (205 grabaciones) una **barra de progreso** con ETA es muy
valiosa; para un sujeto suelto (8 pasos) basta la línea por paso que ya existe.

**Librería:** `ProgressMeter` (ya instalada). Es el estándar de facto en Julia,
soporta ETA, throughput y descripción dinámica. Se usa la barra en batch
(progreso sobre N sujetos) y opcionalmente por-paso.

### 4.4 Colores

Sí, con moderación. Cuatro niveles semánticos con `Crayons` (ya instalada):

- `INFO` (azul/cian) — progreso normal
- `SUCCESS` (verde) — paso u operación completada
- `WARNING` (amarillo) — algo revisable, no bloqueante
- `ERROR` (rojo) — fallo

Regla de robustez: **detectar si la salida es un terminal** (`isa(stdout, Base.TTY)`)
y desactivar colores si se redirige a fichero o corre en CI, para no ensuciar los
logs con códigos ANSI. `Crayons` + esa comprobación cubre esto.

### 4.5 Resumen final

Bloque final claro, alimentado por el mismo estado que el manifiesto:

```
─────────── Resumen ───────────
  Unidad:      sub-M05/ses-T2/eyesclosed
  Estado:      ✔ SUCCESS
  Tiempo:      8 min 42 s
  Tablas:      69   Figuras: 21   JSON: 6
  Warnings:    1  (DELTA poco fiable)
  Errores:     0
  Salida:      results/subjects/sub-M05/ses-T2/eyesclosed/
────────────────────────────────
```

En batch, un resumen agregado adicional: procesados / saltados / con error, y
ruta del log y del CSV de la corrida.

---

## 5. Registro de logs

**Problema con el estado actual:** los `println` van solo a pantalla y se pierden;
`pipeline_log.txt` es manual por sujeto. Falta un log persistente de la corrida
completa con niveles.

**Diseño:** `LoggingExtras` (ya instalada) + `Logging` de la stdlib.

- Un `TeeLogger` envía cada mensaje a **dos destinos**: consola (formateado, con
  color, nivel ≥ `Info`) y `results/logs/run_{timestamp}.log` (texto plano, nivel
  ≥ `Debug`, sin color). El código emite `@info`/`@warn`/`@error`/`@debug`
  estándar; no hay `println` de estado.
- **Niveles:** `DEBUG` (detalle interno, solo a fichero) · `INFO` (progreso) ·
  `WARNING` (revisable) · `ERROR` (fallo). El nivel de consola es configurable
  por flag `--verbose`/`--quiet`.
- **Rotación:** por timestamp de corrida (un fichero por ejecución) es suficiente
  y sencilla; evita la complejidad de rotación por tamaño. Una tarea de limpieza
  opcional (`--keep-logs N`) conserva los N más recientes.
- `pipeline_log.txt` por sujeto se conserva (traza local autocontenida junto a
  sus resultados); el `run_{timestamp}.log` es la traza global de la corrida.

---

## 6. Tablas en terminal (`PrettyTables`)

`PrettyTables` (ya instalada) para todas. Las realmente útiles en un pipeline EEG:

1. **Preflight** — check / estado (OK·ERROR) / detalle.
2. **Datos cargados** — sujeto, canales, fs, duración, nº épocas, condición.
3. **Estado de procesamiento (batch)** — la más valiosa:
   ```
   Subject   ICA    Segment   wPLI    Figures   Status
   M05        ✔      99/100    ✔        21       include
   M07        ✔      88/95     ✔        21       warning
   ```
4. **Tiempo por etapa** — Load/QC/Filter/ICA/Segment/Spectral/wPLI/Save + total.
   Identifica cuellos de botella (ICA suele dominar).
5. **Ficheros generados** — recuento por tipo (CSV/PNG/JSON) desde el manifiesto.
6. **Ranking de métricas** (opcional, científico) — top sujetos por potencia de
   banda, canales con mayor conectividad media. Útil pero secundario; no bloquea.

Las tablas 1–5 son de infraestructura (aportan control de calidad del *run*); la
6 es analítica. Priorizar 1–5; la 6 como extra.

---

## 7. Calidad visual de la terminal

Estilo **moderno pero sobrio**, reutilizando la base ya presente (`═`, `▶`, `✓`,
`⚠`, `·`):

- **Separadores** y **cajas** con caracteres de dibujo Unicode (`─═│├└`) para el
  banner y el resumen.
- **Iconos** semánticos discretos: `✔` éxito, `⚠` aviso, `✖` error, `▶` inicio de
  etapa. Unicode, no emojis a color (se ven mal en muchos terminales y en logs).
- **Sin emojis decorativos** en la salida de infraestructura. Un proyecto
  científico gana en credibilidad con sobriedad; los iconos monocromo bastan.
- **Alineación** de columnas y cifras (la da `PrettyTables`).
- **Colores** solo para los 4 niveles semánticos (§4.4), no decorativos.
- **Degradación elegante:** si no hay TTY o el terminal no soporta Unicode
  (`ENV["TERM"]`), caer a ASCII (`+--+`, `[OK]`, `[!]`) automáticamente.

---

## 8. Arquitectura propuesta

Una capa de runtime en `src/runtime/`, transversal a las cuatro rutinas, con
responsabilidad única por módulo:

```
src/runtime/
├── RunContext.jl     — objeto que atraviesa toda la corrida
├── ResultsManager.jl — clean-slate por unidad, staging atómico, _safe_clean
├── RunManifest.jl    — construye y escribe run_manifest.json
├── RunLogging.jl     — TeeLogger (consola + fichero) sobre LoggingExtras
├── Console.jl        — banner, cajas, colores (Crayons), iconos, degradación TTY
├── Progress.jl       — envoltura fina sobre ProgressMeter
├── Preflight.jl      — checks previos + tabla OK/ERROR + fallo temprano
└── Summary.jl        — resumen final (sujeto) y agregado (batch)
```

| Módulo | Responsabilidad | Depende de |
|--------|-----------------|-----------|
| **RunContext** | Estado de la corrida: unidad, tiempos, contadores (figuras/tablas/warnings), config, entorno. Se pasa a todos los pasos; es la fuente única del manifiesto y del resumen | — |
| **ResultsManager** | Resolver la ruta de la unidad, `_safe_clean`, crear staging `{unit}.tmp/`, swap atómico, excluir `cache/` | RunContext |
| **RunManifest** | Recolectar del RunContext + listar ficheros reales del directorio → `run_manifest.json` (vía JSON3) | RunContext, ResultsManager |
| **RunLogging** | Configurar el TeeLogger al inicio, cerrar al final; API `@info/@warn/...` estándar | LoggingExtras |
| **Console** | Todo lo visual: banner, separadores, colores semánticos, iconos, detección de TTY/Unicode | Crayons |
| **Progress** | Barra de progreso batch y por-paso; silenciable | ProgressMeter |
| **Preflight** | Batería de checks, tabla resultado, abortar si crítico | Console, PrettyTables |
| **Summary** | Render del bloque final desde RunContext; versión batch | Console, PrettyTables |

**Punto de integración:** cada script (`run_single_subject.jl`, etc.) crea un
`RunContext`, llama a `Preflight`, abre `RunLogging`, ejecuta la rutina
(que va rellenando el context y usando `Progress`/`Console`), y cierra con
`RunManifest` + `Summary`. El núcleo científico (ICA, wPLI…) **no cambia**: solo
recibe el context y emite `@info` en vez de `println`.

---

## 9. Roadmap de implementación (para Claude Code)

Ordenado por dependencia y riesgo. Cada fase es un bloque de trabajo verificable
de forma independiente.

### Fase 1 — Fundamentos de reproducibilidad (mayor valor, menor riesgo)
1. `ResultsManager.jl`: `_safe_clean` (con guardas de ruta), staging atómico,
   integración en `_save_all_results` y en los dos scripts de grupo.
2. `RunManifest.jl` + `RunContext.jl`: escribir `run_manifest.json` por unidad.
3. Ajustar la estructura de carpetas a §3/§3.1 (`tables/figures/json` + subcarpeta
   `connectivity/`; sin `ensure_dirs` heredado).
4. **Verificación:** reprocesar M05 dos veces; comprobar que no quedan ficheros
   huérfanos y que el manifiesto lista exactamente lo que hay en disco.

### Fase 2 — Observabilidad (logging + preflight + progreso)
5. `RunLogging.jl`: TeeLogger consola/fichero; sustituir los `println` de estado
   por `@info`/`@warn` en `SingleSubjectPipeline.jl`.
6. `Preflight.jl`: checks previos con tabla OK/ERROR y fallo temprano.
7. `Progress.jl`: barra en `run_batch_pipeline.jl`.
8. **Verificación:** una corrida deja `results/logs/run_*.log` completo; un
   preflight con datos ausentes aborta con mensaje claro y código ≠ 0.

### Fase 3 — Experiencia de terminal (consola + tablas + resumen)
9. `Console.jl`: banner, colores semánticos (Crayons), iconos, degradación TTY/ASCII.
10. `Summary.jl` + tablas `PrettyTables` (preflight, datos, estado batch, tiempos).
11. **Verificación:** salida coherente con y sin TTY (redirigida a fichero sin
    códigos ANSI); resumen final correcto en sujeto y en batch.

### Fase 4 — Consolidación y escala
12. Refactor: las cuatro rutinas comparten el mismo esqueleto `RunContext`.
13. Agregación batch: consolidar los `run_manifest.json` en un informe de corrida
    (`results/logs/batch_{timestamp}.json`) con `jq`-abilidad.
14. Documentar la capa runtime en README y AGENTS; tests de `_safe_clean` y del
    swap atómico (los dos puntos con riesgo de pérdida de datos).

### Riesgos y mitigaciones

| Riesgo | Mitigación |
|--------|-----------|
| Un `clean` mal dirigido borra datos válidos | `_safe_clean` con `startswith(realpath, results/)`; tests unitarios dedicados; staging + swap en vez de borrado in-place |
| Swap atómico entre sistemas de ficheros distintos no es atómico | Forzar que staging y destino compartan FS (staging hermano del destino) |
| Colores/Unicode rompen logs o CI | Detección de TTY y `TERM`; degradación a ASCII; fichero de log sin color |
| Regresión al tocar el pipeline científico | La capa runtime no toca la lógica de ICA/wPLI; solo I/O y presentación. Validar contra el caso M05 tras cada fase |
| Caché ICA invalidada por clean-slate | Excluir `cache/` del borrado (configurable) |

### Buenas prácticas transversales
- Ninguna dependencia nueva (todo ya en `Manifest.toml`).
- El núcleo científico permanece intacto; la capa es I/O + presentación.
- Todo lo que se muestra en pantalla existe también en el manifiesto/log (una
  sola fuente de verdad: el `RunContext`).
- Validar cada fase contra M05 antes de pasar a la siguiente.
- `results/` y `deprecated/results/` siguen fuera de git.

---

*Fin del documento. No incluye código: es la base para generar tareas de
implementación con Claude Code, fase por fase.*
