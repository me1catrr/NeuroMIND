# NeuroMIND report

Informe científico-técnico en LaTeX del proyecto **NeuroMIND**. Documenta el fundamento EEG, la metodología de conectividad funcional, el pipeline Julia y una lectura prudente de los resultados generados por las rutinas del repositorio.

## Estructura del proyecto

- `main_es.tex`: versión principal en español del informe oficial.
- `main.tex`: versión inglesa heredada, conservada como fuente secundaria.
- `sections_es/`: contenido modular de la versión española.
  - `summary.tex`: resumen ejecutivo, mapa visual y propósito del cuaderno.
  - `intro.tex`: fundamentos EEG, conectividad funcional y base matemática.
  - `data.tex`: datos, diseño experimental y estructura de la ejecución.
  - `methods.tex`: teoría operativa del pipeline, criterios técnicos y contratos por fase.
  - `results.tex`: resultados por fase del pipeline NeuroMIND e interpretación técnica.
  - `notes.tex`: estado de outputs, decisiones y tareas pendientes.
  - `appendix.tex`: apéndices BIDS, Julia, dashboard, LaTeX, Git/GitHub e IA.
- `sections/`: contenido modular heredado de la versión inglesa.
  - `intro.tex`
  - `data.tex`
  - `methods.tex`
  - `notes.tex`
  - `appendix.tex`
- `report.cls`: clase principal del informe.
- `refs.bib`: base de datos bibliográfica.
- `assets/`: recursos del proyecto.
  - `logo_EPSC.pdf`, `logo_EPSC.png`, `ARR.pdf`
  - `fonts/`
  - `figures/` y `tables/` (recursos enlazados desde resultados de Julia)
  - `figures/results/`: copias estables de figuras y resúmenes usados por la sección 5.
- `archive/`: archivos legacy no usados en compilación.
- `Makefile`, `.latexmkrc`: flujo de compilación.

## Requisitos

- `pdflatex`
- `latexmk`
- `biber`

## Versión principal

La versión principal del proyecto es `main_es.tex`. El PDF final esperado es `build/pdf/main_es.pdf`.

## Compilación

Desde `NeuroMIND/report/`:

```bash
make build-es   # Compila main_es.tex
make watch-es   # Recompilación continua de la versión castellana
make clean      # Limpia auxiliares de ambas versiones
```

## Nota

La compilación oficial genera `build/pdf/main_es.pdf`. La carpeta `build/pdf` debe mantenerse limpia, con artefactos `main_es.*` como salida de referencia.

La sección de resultados se basa actualmente en la ejecución `sub-M05/ses-T2/eyesclosed`. La lectura es técnica y prudente: documenta QC, ICA, segmentación, espectro, wPLI, métricas de red y surrogates/FDR, pero no formula conclusiones clínicas hasta disponer de una ejecución de referencia congelada y resultados agregados.
