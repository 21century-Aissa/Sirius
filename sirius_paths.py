"""Chemins racine du paquet SIRIUS_MOFA (autonome, tout le necessaire sous ce dossier)."""
from __future__ import annotations

from pathlib import Path

SIRIUS_ROOT: Path = Path(__file__).resolve().parent

SIRIUS_EMBEDDED_ROOT: Path = SIRIUS_ROOT / "SIRIUS_embedded"
SIRIUS_WORKSPACE_DIR: Path = SIRIUS_ROOT / "SIRIUS_workspace"
SIRIUS_DIALOGS_DIR: Path = SIRIUS_ROOT / "SIRIUS_dialogs"
SIRIUS_CONVERSATIONS_DIR: Path = SIRIUS_ROOT / "SIRIUS_conversations"

MOFA_PIPELINE_DIR: Path = SIRIUS_EMBEDDED_ROOT / "Mofa_pipeline"
MOFA_PIPELINE_R: Path = MOFA_PIPELINE_DIR / "mofa_pipeline.R"
MOFA_PLOTS_R: Path = MOFA_PIPELINE_DIR / "mofa_plots.R"
MOFA_SURVIVAL_R: Path = MOFA_PIPELINE_DIR / "mofa_survival.R"
MOFA_CLUSTERING_R: Path = MOFA_PIPELINE_DIR / "mofa_clustering.R"
MOFA_GSEA_R: Path = MOFA_PIPELINE_DIR / "mofa_gsea.R"

DEFAULT_WORKSPACE: Path = SIRIUS_WORKSPACE_DIR / "default"
