@AGENTS.md

## Privacidad (regla absoluta)

- **Nunca escribir nombres de archivo crudos de la adquisición** (llevan iniciales del paciente y fecha de registro) en código, comentarios, documentación, tests, mensajes de commit ni en `graphify-out/`. Se escriben anonimizados: `M12_T2_2_*_Ojoscerrados`. Si un script necesita los nombres reales, los lee en ejecución de `data/`.
- Al auditar, contar coincidencias y dar rutas de archivo; **nunca imprimir los valores**.
- `data`, `results` y `mne_brain/results` son enlaces simbólicos (privados, no versionar nada de ahí): `data` → `BRAIN/04_NeuroMIND/data` (iCloud); `results` y `mne_brain/results` → `~/Proyectos/NeuroMIND-results/{results,mne_results}` (local). Los resultados antiguos quedan en `BRAIN/04_NeuroMIND/` como archivo (solo lectura). Ver AGENTS.md §8.
- Hook de privacidad: `git config core.hooksPath scripts/hooks`; auditoría manual: `scripts/hooks/pre-commit --all` (debe dar 0). No usar `--no-verify`.

## graphify

This project has a knowledge graph at graphify-out/ with god nodes, community structure, and cross-file relationships.

Rules:
- For codebase questions, first run `graphify query "<question>"` when graphify-out/graph.json exists. Use `graphify path "<A>" "<B>"` for relationships and `graphify explain "<concept>"` for focused concepts. These return a scoped subgraph, usually much smaller than GRAPH_REPORT.md or raw grep output.
- If graphify-out/wiki/index.md exists, use it for broad navigation instead of raw source browsing.
- Read graphify-out/GRAPH_REPORT.md only for broad architecture review or when query/path/explain do not surface enough context.
- After modifying code, run `graphify update .` to keep the graph current (AST-only, no API cost).
