# GitHub workflow for NeuroMIND

NeuroMIND should be versioned as its own repository, separate from the original
EEG_Julia/NeuroSmart-EEG project.

## What belongs in Git

- Julia source code in `src/`
- command-line entry points in `scripts/`
- example configuration in `config/`
- tests in `tests/`
- dashboard templates and static assets in `web/`
- documentation in `README.md` and `docs/`
- `Project.toml` and `Manifest.toml` for reproducible Julia environments

## What must stay local

- real EEG recordings in `data/`
- generated subject results in `results/`
- clinical derivatives, reports, logs, and exports
- machine-local editor or assistant settings
- credentials, tokens, API keys, and `.env` files

## Suggested branch model

- `main`: stable, runnable code
- `dev`: integration branch for ongoing work
- `codex/<task>` or `claude/<task>`: local AI-assisted changes
- short feature branches for scientific changes, for example:
  - `codex/align-original-pipeline`
  - `codex/csd-wpli-mode`
  - `claude/dashboard-ica-panel`

Before merging scientific changes, compare outputs against a known local dataset
without committing the data or generated result files.
