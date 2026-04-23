from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from pypdf import PdfReader


@dataclass(frozen=True)
class DocumentBlob:
    path: Path
    text: str


def _safe_read_text(path: Path, max_chars: int) -> str:
    try:
        raw = path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        try:
            raw = path.read_text(encoding="latin-1", errors="replace")
        except Exception as e:
            return f"[ERROR reading {path.name}] {e!s}"
    return raw[:max_chars]


def extract_text_from_pdf(path: Path, max_chars: int = 300_000) -> str:
    try:
        reader = PdfReader(str(path))
        parts: list[str] = []
        for page in reader.pages:
            try:
                t = page.extract_text() or ""
            except Exception:
                t = ""
            if t:
                parts.append(t)
        text = "\n\n".join(parts)
        text = re.sub(r"\n{3,}", "\n\n", text)
        return text[:max_chars]
    except Exception as e:
        return f"[ERROR extracting PDF {path.name}] {e!s}"


def iter_workspace_documents(
    work_dir: Path,
    max_files: int = 80,
    max_chars_per_file: int = 50_000,
) -> Iterable[DocumentBlob]:
    """Collecte les documents du work_dir pour les injecter au LLM.

    - PDF : extraction via pypdf
    - CSV/TSV/TXT/LOG : lecture texte
    """
    work_dir = work_dir.resolve()
    if not work_dir.is_dir():
        return

    exts_text = {".txt", ".tsv", ".csv", ".log"}
    exts_pdf = {".pdf"}

    files: list[Path] = []
    for p in sorted(work_dir.rglob("*")):
        if not p.is_file():
            continue
        suf = p.suffix.lower()
        if suf in exts_text or suf in exts_pdf:
            files.append(p)
        if len(files) >= max_files:
            break

    for p in files:
        suf = p.suffix.lower()
        if suf in exts_pdf:
            yield DocumentBlob(path=p, text=extract_text_from_pdf(p, max_chars=max_chars_per_file))
        else:
            yield DocumentBlob(path=p, text=_safe_read_text(p, max_chars=max_chars_per_file))
