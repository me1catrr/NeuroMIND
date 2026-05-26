# Integración de imágenes de resultados en el informe

Fecha: 2026-05-23  
Ejecución fuente: `results/subjects/sub-M05/ses-T2/eyesclosed/`  
Informe destino: `report/main_es.tex`

## Principio de selección

Las imágenes incorporadas al informe deben cumplir al menos una de estas funciones:

- sostener una decisión QC del pipeline;
- mostrar un resultado interpretable por fase;
- conectar una figura con una tabla o resumen exportado;
- evitar depender de capturas del dashboard;
- mejorar la lectura técnica sin introducir conclusiones clínicas nuevas.

No se incorporan imágenes solo por disponibilidad. Los originales permanecen en `results/`; el informe usa copias estables en `report/assets/figures/results/`.

## Fase 1 — Inventario

Se revisaron:

- señales y PSD: `filtered_signal_preview.png`, `psd_all_channels.png`, `band_power_summary.png`;
- matrices wPLI: `wpli_DELTA.png`, `wpli_THETA.png`, `wpli_ALPHA.png`, `wpli_BETA_LOW.png`, `wpli_BETA_MID.png`, `wpli_BETA_HIGH.png`, `wpli_GAMMA.png`;
- topomaps ICA: `figures/ica_topomap_001.png` a `figures/ica_topomap_031.png`;
- resúmenes CSV/JSON: `ica_component_features.csv`, `connectivity_summary.json`, `surrogate_summary.json`, `surrogate_quality.csv`.

## Fase 2 — Copia estable al informe

Se copiaron al informe:

```text
report/assets/figures/results/wpli/
  neuromind_m05_t2_ec_wpli_delta.png
  neuromind_m05_t2_ec_wpli_beta_low.png
  neuromind_m05_t2_ec_wpli_beta_mid.png
  neuromind_m05_t2_ec_wpli_beta_high.png
  neuromind_m05_t2_ec_wpli_gamma.png

report/assets/figures/results/ica/
  neuromind_m05_t2_ec_ica_ic05_jump.png
  neuromind_m05_t2_ec_ica_ic13_jump.png
  neuromind_m05_t2_ec_ica_ic14_line_noise.png
  neuromind_m05_t2_ec_ica_ic21_line_noise.png
```

Las figuras alfa/theta, señal filtrada, PSD y potencia por banda ya estaban copiadas previamente en `report/assets/figures/results/`.

## Fase 3 — Incorporación narrativa

Se añadieron dos bloques nuevos en `report/sections_es/results.tex`:

- panel ICA de componentes con mayor evidencia automática de artefacto;
- atlas complementario de matrices wPLI por bandas no mostradas en la figura principal.

## Fase 4 — Figuras no incorporadas

Quedan fuera del PDF principal, pero conservadas en `results/`:

- matrices completas de p-valores, q-valores y máscaras significativas;
- series temporales ICA completas;
- todos los topomaps ICA no seleccionados;
- duplicados bajo la ruta legacy `results/M05/T2/`.

Motivo: aportan trazabilidad o material de dashboard, pero saturan el informe si no hay una pregunta específica que los requiera.

## Cautela interpretativa

La inclusión de más figuras no cambia la conclusión técnica actual: la ejecución individual `M05/T2/EC` sirve para validar el pipeline y revisar patrones observados, pero no produce aristas wPLI significativas tras FDR. Las imágenes deben leerse como soporte visual de QC y resultados observados, no como evidencia clínica concluyente.
