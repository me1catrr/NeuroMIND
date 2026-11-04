# NeuroSmart-EEG dentro de NeuroMIND

Esta carpeta es el repositorio **me1catrr/NeuroSmart-EEG** importado con todo su historial (`git subtree`) el 2026-11-04. El repositorio original se **archivará** (solo lectura) en GitHub.

- **Qué contiene:** los 8 cuadernos Pluto didácticos (`Pluto/`), su exportación HTML (`docs/`, `exports/`), el módulo antiguo `EEG_Julia` (`src/`) y las herramientas de publicación (`tools/pages/`).
- **Estado:** histórico. El pipeline vigente es el de NeuroMIND (`src/` en la raíz del repo). Los cuadernos usan su propio entorno (`Project.toml`/`Manifest.toml` de esta carpeta) y el módulo de `src/` de esta misma carpeta.
- **Web:** no se publica. Las exportaciones HTML (`docs/`, `exports/`) se guardan solo como copia y se pueden abrir en local con el navegador. El GitHub Pages de NeuroMIND queda reservado para el libro Quarto del informe.
- **Datos:** como en el repo original, `data/`, `results/` y `Javier_results/` no se versionan.

Para abrir un cuaderno:

```bash
cd legacy/neurosmart
julia --project=. -e 'using Pkg; Pkg.instantiate(); using Pluto; Pluto.run()'
```
