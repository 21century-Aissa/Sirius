# SIRIUS_MOFA

Assistant "Sirius" avec UI MOFA (CustomTkinter) + pipeline MOFA (R) + chat Ollama.

## Lancer (Windows)

1. Installer les dependances Python:

```bat
pip install -r requirements.txt
```

2. Lancer:

```bat
run_sirius_mofa_windows.bat
```

## Notes

- Le pipeline R necessite `Rscript` sur le PATH et les packages R requis (MOFA2, etc.).
- Les documents du `work_dir` (PDF/CSV/TSV/TXT) sont automatiquement fournis au chat (PDF via extraction de texte `pypdf`).
