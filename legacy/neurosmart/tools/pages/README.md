# Publicacion GitHub Pages (`tools/pages`)

Guia principal y unica del flujo de publicacion web para este repositorio Julia + Pluto.

## Que se publica

- GitHub Pages publica desde `main:/docs`.
- Portada: `docs/index.html`.
- Modulos HTML: `docs/Pluto/<Modulo>/<Notebook>.html`.

## Que NO se publica

- Carpetas de resultados internos o no publicables.
- `Javier_results` no forma parte de la web publicada.
- Notebooks fuente `.jl` no se publican directamente como pagina web.

## Ubicacion de cada tipo de archivo

- Fuente cientifica: `Pluto/<Modulo>/<Modulo>.jl`
- Staging de HTML exportado manualmente: `exports/Pluto/<Modulo>/` (carpeta completa)
- Destino publicado en Pages: `docs/Pluto/<Modulo>/<Notebook>.html`

## Scripts

- `tools/pages/build.jl`
  - Busca exportaciones reales en `exports/Pluto/...` (o en una ruta custom opcional).
  - Si existe exportacion real del modulo, copia la carpeta completa a `docs/Pluto/...`.
  - Valida referencias locales `./...` del `index.html` (JS/CSS/assets).
  - Si faltan assets requeridos, genera una pagina de diagnostico de exportacion incompleta.
  - Si no existe, genera un placeholder claro con enlace de vuelta a portada.
  - Regenera `docs/pages_manifest.json`.

- `tools/pages/publish.jl`
  - Ejecuta `build.jl`.
  - Muestra `git status --short`.
  - Si no hay cambios: termina correctamente.
  - Si hay cambios: `git add .` -> `git commit` -> `git push`.
  - Fallback unico para upstream ausente: `git push --set-upstream origin <rama_actual>`.

- `tools/pages/export_pluto.jl`
  - Exporta notebooks Pluto a HTML estatico sin usar la web.
  - Salida: `exports/Pluto/<Modulo>/<Notebook>.html` (+ assets).
  - Si no indicas modulos, exporta todos.

## Comandos

```bash
julia tools/pages/build.jl
```

```bash
julia tools/pages/build.jl "/ruta/custom/exports/Pluto"
```

```bash
julia tools/pages/publish.jl
```

```bash
julia tools/pages/publish.jl "mensaje"
```

```bash
julia --project=. tools/pages/export_pluto.jl
```

```bash
julia --project=. tools/pages/export_pluto.jl ICA Spectral Connectivity
```

## Flujo recomendado

1. Editar notebook fuente en `Pluto/<Modulo>/<Modulo>.jl`.
2. Exportar HTML (CLI recomendado):
   - `julia --project=. tools/pages/export_pluto.jl`
3. Guardar la exportacion completa en `exports/Pluto/<Modulo>/` (incluyendo `index.html` y assets).
4. Ejecutar `julia tools/pages/build.jl`.
5. Verificar `docs/index.html` y `docs/Pluto/<Modulo>/index.html`.
6. Ejecutar `julia tools/pages/publish.jl "mensaje"`.

## Nota sobre exportacion de Pluto

Este repositorio no fuerza una exportacion automatica de Pluto desde `build.jl`.
La exportacion HTML se considera un paso manual y verificable; `build.jl` solo sincroniza a `docs/` o genera placeholders cuando falta exportacion real.  
Si el `index.html` referencia archivos locales que no estan en `exports/Pluto/<Modulo>/`, el build publica un diagnostico de exportacion incompleta para evitar una pagina rota en GitHub Pages.
