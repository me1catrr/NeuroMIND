# Changelog

## 2026-05-23 - Ajuste de portada

- Compactada la portada para evitar que autores y fecha saltaran fuera de la primera página.
- Añadida geometría local de portada y reducidos logo, visual TikZ y tamaño del título.
- PDF español recompilado y portada verificada visualmente.

## 2026-05-23 - Refinamiento científico por fases

- Ajustado el título y el resumen para evitar contradicción entre CSD metodológico y resultados actuales con `use_csd=false`.
- Añadida una justificación de wPLI basada en Vinck et al. 2011, comparándolo con coherencia, PLV, ImC y PLI.
- Añadida una tabla de conceptos de conectividad para fijar que NeuroMIND trabaja con conectividad funcional, no dirigida, estática y a nivel de sensores.
- Reforzada la integración del protocolo RS-MIND en ICA, segmentación, baseline, rechazo de artefactos y FFT.
- Ampliada la motivación clínica de EM con hipótesis prudentes y sin conclusiones grupales no sustentadas.
- Añadida la sección `Conclusiones operativas y próximos pasos`.
- Añadidas referencias reales de Vinck 2011, Friston 2011 y Cao 2022; retiradas entradas bibliográficas de plantilla o incompletas.

## 2026-05-23 - Integración de apuntes internos de ICA

- Revisados los apuntes ICA escaneados aportados como material interno de apoyo.
- Ampliada la subsección de Fase 2 con una lectura más completa de centrado, whitening/sphering, FastICA, demixing/unmixing y back-projection.
- Añadida una tabla operativa de ICA y una caja de cautelas sobre supuestos, no gaussianidad e indeterminación de componentes.
- Registrada la revisión en `../docs/ICA_NOTES_INTEGRATION.md`.
- No se incorporaron capturas manuscritas al informe; se usaron para reforzar el texto y conservar una estética vectorial/LaTeX.

## 2026-05-23 - Figuras reales del pipeline incorporadas

- Incorporado un panel de topomaps ICA seleccionados en la sección de resultados.
- Incorporado un atlas complementario de matrices wPLI por bandas a partir de imágenes exportadas por NeuroMIND.
- Copiadas las figuras seleccionadas a `assets/figures/results/ica/` y `assets/figures/results/wpli/`.
- Actualizada la tabla de selección de artefactos para diferenciar figuras incluidas, resumidas y reservadas.
- Mantenida la cautela inferencial: las figuras amplían la trazabilidad visual, pero no modifican el resultado FDR de la ejecución individual.

## 2026-05-23 - Secciones alineadas con pipeline/dashboard

- Renombrada la sección 4 a `Teoría operativa del pipeline por fases`.
- Renombrada la sección 5 a `Resultados por fase del pipeline NeuroMIND`.
- Reetiquetadas las subsecciones de metodología y resultados para que sigan las fases NeuroMIND y sus paneles del dashboard.
- Actualizadas las tablas de trazabilidad y contratos operativos para usar nombres F0--F8 coherentes.
- Conservada la teoría ya redactada, pero agrupada bajo títulos que diferencian teoría de resultados por fase.

## 2026-05-23 - Consolidación documental

- Actualizado `README.md` del informe para describir la estructura española actual por secciones.
- Documentada la sección 5 como resultados fase a fase basados en la ejecución `sub-M05/ses-T2/eyesclosed`.
- Explicitado que la lectura de resultados es técnica y prudente hasta disponer de ejecución congelada y agregados por cohorte.
- Alineado el README del informe con el apéndice dashboard y los contratos operativos de metodología.

## 2026-05-23 - Rediseño de resultados fase a fase

- Reorganizada la sección `Resultados del análisis de conectividad funcional EEG` para leer la ejecución NeuroMIND por fases.
- Añadida una tabla de lectura fase-dashboard con artefactos revisados, resultado observado, interpretación técnica y panel asociado.
- Añadida una subsección de ICA/segmentación/rechazo de artefactos con indicadores exportados por NeuroMIND.
- Añadida una subsección de lectura de red derivada con nodos de mayor fuerza normalizada.
- Añadida una tabla final de selección de artefactos incluidos, resumidos o reservados fuera del PDF.
- Conservado el tono prudente: no se formulan conclusiones clínicas a partir de una ejecución individual.

## 2026-05-23 - Contratos operativos por fase

- Actualizada la tabla de trazabilidad metodológica para conectar fases del informe con módulos Julia actuales y outputs persistidos.
- Añadida la subsección `Contratos operativos por fase` en `sections_es/methods.tex`.
- Incorporado un mapa fase a fase de entradas, transformaciones, salidas QC y paneles del dashboard.
- Añadida una caja de control QC que formaliza cuándo una fase puede alimentar conclusiones de la fase siguiente.
- Pulida la maquetación de tablas para mantener legibilidad en una sola columna.

## 2026-05-23 - Apéndice técnico del dashboard NeuroMIND

- Añadido al apéndice español un bloque específico sobre el dashboard NeuroMIND.
- Documentada la arquitectura `scripts/launch_dashboard.jl` -> Genie API (`src/webapp/App.jl`) -> interfaz HTML (`web/views/dashboard.html`) -> outputs en `results/`.
- Añadido un diagrama TikZ de arquitectura dashboard/pipeline/informe.
- Añadida una tabla de componentes y una tabla de paneles 0--14 con endpoints principales y uso técnico.
- Explicitada la diferencia entre numeración metodológica del informe y numeración granular del dashboard.
- Añadidas cautelas para evitar interpretar capturas o vistas interactivas como evidencia estadística independiente.

## 2026-05-23 - Integración como informe oficial de NeuroMIND

- Adaptada la portada de `main_es.tex` para presentar el documento como informe científico-técnico oficial del proyecto NeuroMIND.
- Actualizados resumen ejecutivo e introducción para conectar explícitamente el informe con el pipeline Julia, sus salidas y su función como reporte vivo del proyecto.
- Actualizado `README.md` del informe para explicar que la versión principal es `main_es.tex` y que la compilación oficial se realiza desde `NeuroMIND/report/`.
- Actualizado `CODEX.md` para documentar la ubicación esperada `NeuroMIND/report/`, la relación con `src/report/` y las decisiones de integración.
- La integración conserva rutas internas relativas y copias locales de figuras/resultados para que el PDF compile sin depender de rutas externas.

## 2026-05-23 - Resultados wPLI con surrogates/FDR

- Revisados los nuevos outputs de NeuroMIND para `sub-M05`, `ses-T2`, condición `eyesclosed`, incluyendo matrices wPLI observadas, matrices de p-valores, q-valores, máscaras significativas, estadísticas nulas surrogate y resúmenes JSON/CSV.
- Ampliada `sections_es/results.tex` con una subsección específica de inferencia por surrogates de fase aleatorizada y corrección FDR.
- Añadido un gráfico TikZ/PGFPlots de resumen por banda que compara wPLI observado medio, media nula surrogate y p-valor medio.
- Actualizada la tabla de aristas wPLI destacadas para mostrar p-valores empíricos y q-valores FDR asociados.
- Documentado el resultado inferencial de la ejecución individual: 0/465 aristas significativas en todas las bandas y `n_sig_total = 0`, sin convertir esta lectura en conclusión clínica.
- Copiados a `assets/figures/results/` los resúmenes `neuromind_m05_t2_ec_surrogate_summary.json`, `neuromind_m05_t2_ec_surrogate_quality.csv` y `neuromind_m05_t2_ec_significant_connections.csv` para trazabilidad local.
- Actualizadas las notas metodológicas, el resumen ejecutivo, el apéndice Julia y `CODEX.md` para no tratar la fase surrogate/FDR como inexistente; queda marcada como activa exploratoria hasta congelar configuración final con CSD y análisis comparativo.

## 2026-05-22 - Ajuste a una columna y limpieza visual por capturas

- Cambiada la maquetación española a una sola columna real, sin margen lateral de notas; las antiguas notas laterales se renderizan ahora como notas al pie para evitar solapes con figuras, fórmulas y captions.
- Reemplazada la figura raster/PDF de tres perspectivas de conectividad por una reconstrucción TikZ (`\NSConnectivityPerspectives`) integrada en `assets/figures/generated/neurosmart_generated_figures.tex`.
- Reemplazado el flujo de control QC raw por una figura TikZ nueva (`\NSRawQCFlow`) y archivada la imagen previa de menor calidad.
- Reemplazado el diagrama FFT blanco/negro por una versión TikZ coherente con la paleta del informe (`\NSFFTFlow`).
- Retirada del flujo narrativo la figura duplicada `theory/complejos_fase_wpli_idea_clave.pdf`; se conserva archivada porque la explicación equivalente queda mejor integrada con `\NSWPLIIntuition`.
- Añadidas anclas únicas para líneas de pseudocódigo, eliminando los avisos de identificadores duplicados de `algorithmicx`/`hyperref`.
- Archivados en `assets/figures/archive/replaced_in_tikz/` los recursos sustituidos:
  - `connectivity/conectividad_tres_perspectivas.pdf`
  - `methodology/qc_raw_umbrales_flujo_control.png`
  - `methodology/fase0_timeline_procesamiento.png`
  - `theory/complejos_fase_wpli_idea_clave.pdf`
- Eliminada la repetición del mapa del pipeline dentro de metodología; la sección remite al mapa visual global para evitar redundancia y problemas de composición.
- Reubicadas las ecuaciones de QC de Fase 0 desde notas al pie al cuerpo principal: ahora aparecen como caja matemática junto al gráfico ilustrativo de media, desviación estándar, extremos, RMS, asimetría y curtosis.
- Simplificado el mapeo Pluto -> src para evitar notas al pie largas y cortes visuales entre páginas.
- Reubicada la explicación de `ocular_score`, `muscle_score`, `line_score` y `jump_score` desde notas al pie a una caja `Control QC` con gráfico de barras del máximo y umbral ICA.

## 2026-05-22 - Resultados NeuroMIND y reorganización de assets

### Nueva sección de resultados

- Añadida `sections_es/results.tex` con la sección `Resultados del análisis de conectividad funcional EEG`.
- Integradas figuras exportadas por NeuroMIND y copiadas a `assets/figures/results/`:
  - `neuromind_m05_t2_ec_senal_raw_filtrada.png`
  - `neuromind_m05_t2_ec_senal_raw_filtrada_crop.png`
  - `neuromind_m05_t2_ec_psd_canales.png`
  - `neuromind_m05_t2_ec_potencia_bandas.png`
  - `neuromind_m05_t2_ec_wpli_alpha.png`
  - `neuromind_m05_t2_ec_wpli_theta.png`
- La interpretación se limita a resultados descriptivos de una ejecución individual (`sub-M05`, `ses-T2`, `eyesclosed`) y diferencia resultado observado, interpretación técnica y limitaciones.
- Se documenta explícitamente que la ejecución revisada usa `use_csd=false`; la inferencia surrogate/FDR queda actualizada en la entrada del 2026-05-23.

### Figuras y organización

- Reorganizada `assets/figures/` en:
  - `theory/`
  - `methodology/`
  - `connectivity/`
  - `results/`
  - `generated/`
  - `archive/`
- Revisada la figura vectorial `theory/complejos_fase_wpli_idea_clave.pdf`; finalmente se archivó por redundancia con la figura TikZ de wPLI integrada en metodología.
- Movidas a `archive/` las familias heredadas no usadas directamente: NeuroMIND, Cerebro y Julia.
- Actualizadas las rutas en `sections_es/intro.tex`, `sections_es/data.tex` y `sections_es/methods.tex`.

### Documentación y limpieza

- Renombrado el archivo de estado anterior a `CODEX.md` y reescrito como guía operativa actual.
- Limpieza de `build/pdf`: eliminados los artefactos de la compilación inglesa `main.*`, conservando solo `main_es.*`.
- Eliminados `.DS_Store` encontrados en `assets/figures/`.

## 2026-05-22 - Revisión de figuras wPLI vectoriales

- Revisada la carpeta `assets/figures/NeuroMIND/wpli_vector_figures/`.
- Se verificó que contiene pares SVG/PDF vectoriales de buena base, aunque algunas versiones PDF muestran textos cortados o solapados al renderizarse en el documento.
- Añadida una figura TikZ nueva, `\NSWPLIIntuition`, inspirada en esa colección pero reconstruida desde cero con la paleta y estilos del informe.
- Incorporada la figura en la sección 7.3 para explicar desfase cero, desfase no nulo, componente imaginaria y matriz wPLI sin introducir dependencias SVG ni `shell-escape`.
- Retirada de esa zona una figura lateral previa de conectividad funcional porque competía visualmente con el nuevo bloque a ancho completo.

## 2026-05-22 - Mejora visual y estructural del cuaderno español

### Cambios de estructura

- Añadida una sección inicial `Resumen ejecutivo y propósito del cuaderno`.
- Reordenada la narrativa de la versión española para abrir con propósito, mapa visual, diseño experimental y pipeline global antes de la introducción teórica.
- Reescrito el apéndice español para agrupar mejor BIDS, Julia, LaTeX, Git/GitHub e IA, con tablas más limpias y texto unificado en castellano.

### Figuras añadidas

Se añadieron figuras vectoriales generadas desde cero en `assets/figures/generated/neurosmart_generated_figures.tex`:

- Portada visual con cerebro/red EEG y líneas de señal.
- Mapa visual del estudio.
- Diseño experimental transversal y longitudinal.
- Pipeline global Fase 0--8.
- Flujo BrainVision -> BIDS.
- Montaje EEG 10--20 con regiones, REF FCz y GND Fpz.
- Cascada de filtrado.
- Flujo ICA.
- Segmentación, baseline y rechazo de artefactos.
- CSD + wPLI con fórmula.
- Surrogates + FDR.
- Dashboard de outputs por ejecución.

### Paquetes y macros

- Añadido `macros/neurosmart_style.tex` con paleta, estilos de títulos/captions y cajas:
  - `ideaclave`
  - `decisionmetodologica`
  - `controlqc`
  - `parametros`
  - `pendiente`
- Añadida la librería TikZ `decorations.pathreplacing` para llaves de agrupación en diagramas.

### Problemas detectados

- La fase de surrogates/FDR estaba pendiente en esta revisión, pero queda actualizada en la entrada del 2026-05-23 tras revisar los nuevos outputs de NeuroMIND.
- El documento conserva algunas figuras rasterizadas previas de NeuroMIND como apoyo conceptual; las nuevas figuras principales del rediseño son vectoriales/TikZ.

### Tareas pendientes

- Revisar visualmente el PDF completo página por página para ajustar floats o notas laterales que puedan quedar desplazadas.
- Congelar la fase surrogate/FDR en una ejecución de referencia con CSD, número final de surrogates, semilla y snapshot de configuración.
- Congelar una ejecución de referencia con `run_<fecha>_<branch>_<cfg>`, commit hash y snapshot de configuración.

## 2026-05-22 - Auditoría de figuras

### Cambios

- Revisada la carpeta `assets/figures/` por formato, resolución y uso real en LaTeX.
- Añadida una figura de alta resolución sobre flujo de control de calidad raw:
  - `NeuroMIND/02_EEG_extracts_curated_highres/02_calidad_raw_EEG/08_umbrales_y_flujo_de_control_de_calidad.png`
- Sustituidas por PDF vectorial varias figuras que antes se incluían como PNG:
  - `04_misma_red_tres_perspectivas.pdf`
  - `02_concepto_epoch_ventana_2s.pdf`
  - `06_proceso_correccion_baseline.pdf`
  - `09_revision_rechazo_artefactos.pdf`
  - `07_eeg_conectividad_funcional.pdf`

### Limpieza

- Eliminados los PNG duplicados no usados tras migrar esas inclusiones a PDF vectorial.
- No se eliminaron las familias highres restantes: no se detectaron problemas de resolución en los PNG reales y algunas sirven como reserva visual para futuras secciones.
