"""Execution analytique SIRIUS_MOFA: pipeline R MOFA (scripts embarques)."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple

from sirius_paths import (
    MOFA_CLUSTERING_R,
    MOFA_GSEA_R,
    MOFA_PIPELINE_R,
    MOFA_PLOTS_R,
    MOFA_SURVIVAL_R,
    SIRIUS_CONVERSATIONS_DIR,
)


def find_rscript() -> str | None:
    return shutil.which("Rscript")


def ensure_conversation_dir() -> Path:
    SIRIUS_CONVERSATIONS_DIR.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    d = SIRIUS_CONVERSATIONS_DIR / ts
    d.mkdir(parents=True, exist_ok=True)
    return d


def run_mofa_pipeline(
    work_dir: str,
    selected_omics: Dict[str, str],
    num_factors: int,
    max_r2: float | None,
    has_metadata: bool,
    metadata_file: str,
    rscript_exe: str | None = None,
    timeout_s: int = 6 * 3600,
) -> Tuple[int, str]:
    root = Path(work_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)

    rscript = rscript_exe or find_rscript()
    if not rscript:
        return -1, "Rscript introuvable dans le PATH. Installez R."
    if not MOFA_PIPELINE_R.is_file():
        return -1, f"Script R introuvable: {MOFA_PIPELINE_R}"

    cmd = [
        rscript,
        str(MOFA_PIPELINE_R),
        f"work_dir={str(root)}",
        f"num_factors={int(num_factors)}",
        *([] if max_r2 is None else [f"max_r2={float(max_r2)}"]),
        f"has_metadata={'TRUE' if has_metadata else 'FALSE'}",
        f"metadata_file={metadata_file}",
        f"python_path={sys.executable}",
    ]
    for name, path in selected_omics.items():
        cmd.append(f"omic_{name}={path}")

    try:
        proc = subprocess.run(
            cmd,
            cwd=str(root),
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            env={**os.environ},
            timeout=timeout_s,
        )
    except subprocess.TimeoutExpired:
        return -2, f"Pipeline R interrompu : timeout {timeout_s} s depasse."
    except Exception as e:
        return -3, f"Erreur lancement Rscript : {e}"

    log = (proc.stdout or "") + "\n" + (proc.stderr or "")
    return proc.returncode, log


def run_mofa_plots(
    work_dir: str,
    model_path: str,
    plot_type: str,
    out_name: str,
    extra_params: Dict[str, str],
    rscript_exe: str | None = None,
    timeout_s: int = 2 * 3600,
) -> Tuple[int, str]:
    root = Path(work_dir).resolve()
    rscript = rscript_exe or find_rscript()
    if not rscript:
        return -1, "Rscript introuvable dans le PATH."
    if not MOFA_PLOTS_R.is_file():
        return -1, f"Script R introuvable: {MOFA_PLOTS_R}"

    cmd = [
        rscript,
        str(MOFA_PLOTS_R),
        f"model_path={model_path}",
        f"work_dir={str(root)}",
        f"plot_type={plot_type}",
        f"out_name={out_name}",
    ]
    for k, v in extra_params.items():
        cmd.append(f"{k}={v}")

    try:
        proc = subprocess.run(
            cmd,
            cwd=str(root),
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            env={**os.environ},
            timeout=timeout_s,
        )
    except subprocess.TimeoutExpired:
        return -2, f"Plots R interrompus : timeout {timeout_s} s depasse."
    except Exception as e:
        return -3, f"Erreur lancement Rscript : {e}"

    log = (proc.stdout or "") + "\n" + (proc.stderr or "")
    return proc.returncode, log


def run_mofa_survival(
    work_dir: str,
    model_path: str,
    analysis: str,
    out_name: str,
    extra_params: Dict[str, str],
    rscript_exe: str | None = None,
    timeout_s: int = 2 * 3600,
) -> Tuple[int, str]:
    root = Path(work_dir).resolve()
    rscript = rscript_exe or find_rscript()
    if not rscript:
        return -1, "Rscript introuvable dans le PATH."
    if not MOFA_SURVIVAL_R.is_file():
        return -1, f"Script R introuvable: {MOFA_SURVIVAL_R}"

    cmd = [
        rscript,
        str(MOFA_SURVIVAL_R),
        f"model_path={model_path}",
        f"work_dir={str(root)}",
        f"analysis={analysis}",
        f"out_name={out_name}",
    ]
    for k, v in extra_params.items():
        cmd.append(f"{k}={v}")

    try:
        proc = subprocess.run(
            cmd,
            cwd=str(root),
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            env={**os.environ},
            timeout=timeout_s,
        )
    except subprocess.TimeoutExpired:
        return -2, f"Survival R interrompu : timeout {timeout_s} s depasse."
    except Exception as e:
        return -3, f"Erreur lancement Rscript : {e}"

    log = (proc.stdout or "") + "\n" + (proc.stderr or "")
    return proc.returncode, log


def run_mofa_clustering(
    model_path: str,
    output_dir: str,
    analysis: str,
    method: str = "kmeans",
    k: int = 3,
    var_threshold: float = 0,
    factor_x: int = 1,
    factor_y: int = 2,
    covariates: str | None = None,
    rscript_exe: str | None = None,
    timeout_s: int = 2 * 3600,
) -> Tuple[int, str]:
    root = Path(output_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)

    rscript = rscript_exe or find_rscript()
    if not rscript:
        return -1, "Rscript introuvable dans le PATH."
    if not MOFA_CLUSTERING_R.is_file():
        return -1, f"Script R introuvable: {MOFA_CLUSTERING_R}"

    cmd = [
        rscript,
        str(MOFA_CLUSTERING_R),
        "--model", model_path,
        "--output_dir", str(root),
        "--analysis", analysis,
        "--method", str(method),
        "--k", str(int(k)),
        "--var_threshold", str(var_threshold),
        "--factor_x", str(int(factor_x)),
        "--factor_y", str(int(factor_y)),
    ]
    if covariates is not None and str(covariates).strip() != "":
        if isinstance(covariates, (list, tuple)):
            covariates = ",".join(str(c) for c in covariates)
        cmd.extend(["--covariates", str(covariates)])

    try:
        proc = subprocess.run(
            cmd,
            cwd=str(root),
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            env={**os.environ},
            timeout=timeout_s,
        )
    except subprocess.TimeoutExpired:
        return -2, f"Clustering R interrompu : timeout {timeout_s} s depasse."
    except Exception as e:
        return -3, f"Erreur lancement Rscript : {e}"

    log = (proc.stdout or "") + "\n" + (proc.stderr or "")
    return proc.returncode, log


def run_mofa_gsea(
    work_dir: str,
    model_path: str,
    analysis: str,
    out_name: str,
    extra_params: Dict[str, str],
    rscript_exe: str | None = None,
    timeout_s: int = 3 * 3600,
) -> Tuple[int, str]:
    root = Path(work_dir).resolve()
    rscript = rscript_exe or find_rscript()
    if not rscript:
        return -1, "Rscript introuvable dans le PATH."
    if not MOFA_GSEA_R.is_file():
        return -1, f"Script R introuvable: {MOFA_GSEA_R}"

    cmd = [
        rscript,
        str(MOFA_GSEA_R),
        f"model_path={model_path}",
        f"work_dir={str(root)}",
        f"analysis={analysis}",
        f"out_name={out_name}",
    ]
    for k, v in extra_params.items():
        cmd.append(f"{k}={v}")

    try:
        proc = subprocess.run(
            cmd,
            cwd=str(root),
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            env={**os.environ},
            timeout=timeout_s,
        )
    except subprocess.TimeoutExpired:
        return -2, f"GSEA R interrompu : timeout {timeout_s} s depasse."
    except Exception as e:
        return -3, f"Erreur lancement Rscript : {e}"

    log = (proc.stdout or "") + "\n" + (proc.stderr or "")
    return proc.returncode, log


def list_result_files(work_dir: str, max_files: int = 300) -> List[Path]:
    root = Path(work_dir).resolve()
    if not root.is_dir():
        return []
    exts = {".pdf", ".png", ".jpg", ".jpeg", ".webp", ".txt", ".tsv", ".csv", ".log", ".hdf5"}
    found: List[Path] = []
    for p in sorted(root.rglob("*")):
        if p.is_file() and p.suffix.lower() in exts:
            found.append(p)
        if len(found) >= max_files:
            break
    return found
