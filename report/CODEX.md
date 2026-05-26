# CODEX.md

## Estado actual

Proyecto LaTeX del informe científico-técnico oficial de NeuroMIND sobre conectividad funcional EEG en esclerosis múltiple. La versión activa es `main_es.tex` y el PDF principal esperado es `build/pdf/main_es.pdf`. Integrado en el repositorio Julia, este directorio debe vivir como `NeuroMIND/report/`.

El informe integra ahora tres niveles:

- Base teórica: EEG, bandas, conectividad funcional, números complejos, fase, wPLI, CSD, surrogates y FDR.
- Metodología/pipeline: BIDS, filtrado, ICA, segmentación, baseline, rechazo de artefactos, FFT/CSD/wPLI e inferencia.
- Resultados técnicos: primera lectura de una ejecución NeuroMIND (`sub-M05`, `ses-T2`, `eyesclosed`) con QC, PSD, potencia por banda, matrices wPLI descriptivas y surrogates/FDR.

## Cambios recientes

- Añadida la sección `Resultados del análisis de conectividad funcional EEG` en `sections_es/results.tex`.
- Insertada la sección en `main_es.tex` antes del seguimiento operativo.
- Adaptada portada, resumen e introducción para presentar el documento como informe oficial de NeuroMIND.
- Incorporadas figuras de resultados exportadas por NeuroMIND y copiadas localmente a `assets/figures/results/`.
- Cambiada la versión española a una columna real; las notas laterales se convierten en notas al pie para mejorar legibilidad y evitar solapes.
- Reconstruidas en TikZ varias figuras problemáticas detectadas en revisión visual: perspectivas de conectividad, flujo QC raw y flujo FFT.
- Reubicadas las definiciones matemáticas de QC de Fase 0 en una caja dentro del cuerpo del texto, vinculada al gráfico ilustrativo de estadísticos básicos.
- Reubicada la lógica de scores ICA en una caja de control QC con gráfico de umbral y decisión por máximo.
- Archivada la figura vectorial sobre números complejos, fase y wPLI por redundancia con la figura TikZ `\NSWPLIIntuition`.
- Corregidas las anclas internas de líneas de algoritmos para evitar identificadores PDF duplicados.
- Reorganizada `assets/figures/` en carpetas semánticas.
- Archivadas figuras heredadas no usadas en `assets/figures/archive/`.
- Renombrado el archivo de estado anterior a `CODEX.md`.
- Preparada la integración como `NeuroMIND/report/`, manteniendo rutas internas relativas para que el informe compile desde la nueva ubicación.
- Refinado el informe por fases: coherencia CSD/resultados, justificación wPLI, conceptos de conectividad, protocolo RS-MIND, motivación clínica EM, separación de apéndices operativos y cierre de próximos pasos.
- Añadida la sección `Conclusiones operativas y próximos pasos`.

## Organización de figuras

Estructura activa:

- `assets/figures/theory/`: fundamentos teóricos, cerebro, oscilaciones, números complejos/fase.
- `assets/figures/methodology/`: figuras de diseño, Fase 0, QC raw, segmentación, baseline y rechazo de artefactos.
- `assets/figures/connectivity/`: figuras conceptuales de conectividad.
- `assets/figures/results/`: figuras y copia de metadatos seleccionados desde NeuroMIND.
- `assets/figures/generated/`: figuras TikZ generadas desde cero para el informe.
- `assets/figures/archive/`: material heredado no usado directamente, preservado para no perder recursos.

Figuras incorporadas al informe desde NeuroMIND:

- `results/neuromind_m05_t2_ec_senal_raw_filtrada.png`
- `results/neuromind_m05_t2_ec_senal_raw_filtrada_crop.png`
- `results/neuromind_m05_t2_ec_psd_canales.png`
- `results/neuromind_m05_t2_ec_potencia_bandas.png`
- `results/neuromind_m05_t2_ec_wpli_alpha.png`
- `results/neuromind_m05_t2_ec_wpli_theta.png`

Figuras metodológicas/conceptuales reconstruidas en TikZ:

- `\NSConnectivityPerspectives`
- `\NSRawQCFlow`
- `\NSFFTFlow`
- `\NSWPLIIntuition`

Figuras descartadas para esta versión:

- Heatmaps wPLI de delta, beta baja/media/alta y gamma: calidad suficiente, pero redundantes para una primera lectura técnica.
- Topomaps ICA individuales: útiles para depuración, pero demasiado granulares para el flujo narrativo actual.
- Familias vectoriales completas de wPLI y complejos: archivadas como material reutilizable; se priorizaron figuras TikZ o piezas puntuales coherentes con el informe.
- Figuras sustituidas por TikZ: `conectividad_tres_perspectivas.pdf`, `qc_raw_umbrales_flujo_control.png`, `fase0_timeline_procesamiento.png` y `complejos_fase_wpli_idea_clave.pdf`; conservadas en `assets/figures/archive/replaced_in_tikz/`.

## Decisiones técnicas

- Las salidas de NeuroMIND se interpretan como resultados descriptivos de una ejecución individual, no como conclusiones clínicas.
- La ejecución revisada usa `use_csd=false`; por tanto, las matrices wPLI actuales son sensor--sensor preliminares.
- CSD permanece como etapa metodológica prevista/opcional, no como propiedad de la ejecución preliminar ya incorporada.
- La sección incorpora surrogates de fase aleatorizada y FDR: en la ejecución `M05/T2/EC` no sobreviven aristas significativas (`n_sig_total = 0`).
- Las rutas del informe apuntan a copias locales dentro de `assets/figures/results/`, evitando depender de rutas externas de NeuroMIND para compilar.
- La carpeta `report/` contiene el informe reproducible; `src/report/` queda reservado para código Julia de generación/exportación de reportes.
- Los recursos no usados se archivaron en lugar de borrarse definitivamente.

## Pendientes

- Relanzar NeuroMIND con configuración final congelada: CSD activo si procede, número definitivo de surrogates, semilla y snapshot de configuración.
- Exportar resultados comparables para más sujetos, condiciones EO/EC y grupos.
- Generar tablas agregadas por banda/condición/grupo antes de formular conclusiones científicas.
- Revisar visualmente el PDF completo tras cada reorganización grande de figuras.
- Resolver warnings menores de maquetación (`underfull`) cuando se estabilice el contenido.

## Compilación

Comando principal:

```bash
make build-es
```

Salida esperada:

```text
build/pdf/main_es.pdf
```

Comprobaciones útiles:

```bash
rg -n "LaTeX Error|File .* not found|Undefined control sequence|Reference .* undefined|Citation .* undefined" build/pdf/main_es.log
rg -n "NeuroMIND/|Cerebro/|Julia/" main_es.tex sections_es
```

`build/pdf` debe conservar únicamente los artefactos de la compilación española (`main_es.*`) tras la limpieza final.
