# Revisión visual del PDF — Fase 6

Fecha: 2026-05-23  
PDF revisado: `report/build/pdf/main_es.pdf`  
Páginas: 47  
Render temporal usado para QA: `/private/tmp/neuromind_pdf_review/`

## Alcance

La revisión se centró en las páginas con mayor riesgo tras las fases 2--5:

- portada e índice inicial;
- metodología y contratos operativos;
- sección 5 de resultados fase a fase;
- tablas de ICA, segmentación, red, surrogates/FDR y selección de artefactos;
- apéndice del dashboard y tabla de endpoints 0--14;
- estado de `build/pdf`.

## Páginas inspeccionadas visualmente

- `1--3`: portada e índice inicial.
- `14--15`: resumen del pipeline, trazabilidad y contratos operativos.
- `30--37`: entrada de resultados, tabla fase-dashboard, QC, ICA, espectral, wPLI, red e inferencia.
- `39--45`: seguimiento operativo, apéndices BIDS/Julia/dashboard y tabla de endpoints.

## Resultado

- No se observaron solapes de texto, figuras o tablas.
- Las tablas nuevas son densas pero legibles en A4.
- La sección 5 mantiene una lectura secuencial y separa resultado observado, interpretación técnica y limitaciones.
- El apéndice del dashboard queda coherente con la numeración 0--14 y con los endpoints documentados.
- `build/pdf` conserva únicamente artefactos `main_es.*`.

## Validación LaTeX

- `make build-es` queda actualizado.
- No hay errores LaTeX.
- No hay referencias indefinidas.
- No hay `Overfull \hbox`.
- Permanecen algunos `Underfull \hbox` menores, sin impacto visual relevante.

## Decisión

La versión actual del PDF queda aceptada como cierre de revisión visual técnica. La próxima revisión visual debería repetirse después de incorporar nuevos sujetos, CSD activo o resultados grupales.
