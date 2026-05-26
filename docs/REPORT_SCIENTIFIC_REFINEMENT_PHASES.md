# Refinamiento científico del informe por fases

Fecha: 2026-05-23

## Fuentes revisadas

| Fuente | Extracción | Uso |
|---|---|---|
| `2011_wPLI.pdf` | Texto completo extraíble | Base principal para justificar wPLI frente a coherencia, PLV, ImC y PLI. |
| `Functional and Effective Connectivity A Review.pdf` | Texto extraíble | Distinción conceptual entre conectividad funcional y efectiva. |
| `Brain Functional and Effective Connectivity based on EEG.pdf` | Texto extraíble | Contexto EEG: medidas funcionales/efectivas, sensor/fuente y visualización. |
| `Protocolo_RS_MIND.pdf` | Escaneo sin OCR útil | Se usó como protocolo operativo interno según las decisiones ya especificadas por el proyecto. |
| `Apuntes_Conectividad.pdf`, `Apuntes_FFT_*`, `Mapa_EEG.pdf`, `Teoría_BrainVision_Synchro.pdf` | Escaneos sin OCR útil | Mantienen valor como apoyo interno; no se incorporan como citas formales. |
| `Functional Connectivity in MS Recents Findings and Future.pdf` | Escaneo sin OCR útil | Se usa solo como orientación temática prudente sobre EM, sin afirmaciones clínicas fuertes. |

## Cambios por fase

### Fase 1 - Coherencia interna urgente

- Ajustado el título/subtítulo para no presentar CSD como resultado ya activo.
- Aclarado que CSD es etapa metodológica prevista y que la ejecución preliminar actual se interpreta con `use_csd=false` en espacio de sensores.
- Corregida la errata `Hipotesis` -> `Hipótesis`.
- Eliminadas referencias bibliográficas de plantilla o incompletas (`Smith...`, `sun2012`, `bakhshayesh2019`) y añadidas referencias reales para wPLI y conectividad.

### Fase 2 - Justificación de wPLI

- Añadida una subsección breve sobre wPLI frente a coherencia, PLV, ImC y PLI.
- Incorporada cita de Vinck et al. 2011.
- Añadida cautela sobre tamaño muestral, longitud de segmentos, número de épocas y calidad de señal.

### Fase 3 - Conceptos de conectividad

- Añadida tabla de ejes conceptuales: estructural, funcional, efectiva, dirigida/no dirigida, sensor/fuente y estática/dinámica.
- Explicitado que NeuroMIND trabaja actualmente con conectividad funcional, no dirigida, estática y a nivel de sensores.

### Fase 4 - Protocolo RS-MIND

- Reforzada la sección de ICA con número de componentes igual al número de canales seleccionados cuando procede y revisión semiautomática/manual.
- Reforzadas segmentación de 1 s sin solapamiento, baseline 0--100 ms, rechazo ±70 µV, FFT por canal/segmento y promediado posterior.
- Cuando el pipeline actual difiere del ideal CSD+wPLI, se documenta como adaptación exploratoria y no como error.

### Fase 5 - Motivación clínica EM

- Ampliada ligeramente la motivación clínica sobre alteraciones funcionales de red en EM.
- Incluidas hipótesis prudentes sobre conectividad interhemisférica, bandas alfa/beta, compensación temprana, control top--down y relación futura con EDSS/MSFC.

### Fase 6 - Separación ciencia/documentación interna

- Añadidas notas en apéndices para mantener dashboard, Git/GitHub e IA como trazabilidad, no como evidencia científica principal.
- Se remite la operación diaria a `CODEX.md` y `CHANGELOG.md`.

### Fase 7 - Cierre operativo

- Añadida la sección `Conclusiones operativas y próximos pasos`.
- Se documenta que el pipeline es funcional y reproducible, pero que los resultados actuales son preliminares y no clínicamente concluyentes.

## Validación

- `report/main_es.tex` se compiló con `make build-es`.
- PDF generado: `report/build/pdf/main_es.pdf`.
- Revisión de log sin errores LaTeX, figuras ausentes, referencias/citas indefinidas ni overfull graves.
