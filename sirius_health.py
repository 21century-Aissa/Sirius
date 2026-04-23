"""Verifications prealable (Ollama, Rscript) pour AI Sirius MOFA."""
from __future__ import annotations

import http.client
import json
import os
import shutil
import socket
from urllib.parse import urlparse


def _parse(base_url: str) -> tuple[str, str, int]:
    u = urlparse(base_url)
    host = u.hostname or "127.0.0.1"
    port = u.port or (443 if u.scheme == "https" else 11434)
    scheme = u.scheme or "http"
    return scheme, host, port


def check_ollama(
    base: str | None = None,
    timeout: float = 5.0,
    model: str = "llama3:latest",
    gen_timeout: float = 120.0,
) -> tuple[bool, str]:
    """Verifie que le daemon Ollama repond et qu'une generation courte fonctionne.

    Utilise http.client directement pour eviter tout proxy systeme (registre Windows).
    """
    raw = (base or os.environ.get("OLLAMA_HOST", "http://127.0.0.1:11434")).rstrip("/")
    scheme, host, port = _parse(raw)
    Conn = http.client.HTTPSConnection if scheme == "https" else http.client.HTTPConnection

    # 1) Ping /api/tags (rapide)
    conn = Conn(host, port, timeout=timeout)
    try:
        conn.request("GET", "/api/tags")
        resp = conn.getresponse()
        resp.read()
        if resp.status != 200:
            return False, f"Ollama HTTP {resp.status}"
    except (TimeoutError, socket.timeout):
        return False, "Ollama : delai depasse (tags)."
    except (ConnectionError, OSError) as e:
        return False, f"Ollama injoignable ({e!s}). Lancez `ollama serve`."
    finally:
        try:
            conn.close()
        except Exception:
            pass

    # 2) Generation minimale (force le chargement du modele)
    payload = json.dumps(
        {
            "model": model,
            "prompt": "health-check",
            "stream": False,
            "options": {"num_predict": 1},
        }
    ).encode("utf-8")
    conn = Conn(host, port, timeout=gen_timeout)
    try:
        conn.request(
            "POST",
            "/api/generate",
            body=payload,
            headers={"Content-Type": "application/json"},
        )
        resp = conn.getresponse()
        raw_body = resp.read().decode("utf-8", errors="replace")
        if resp.status != 200:
            details = ""
            try:
                parsed = json.loads(raw_body) if raw_body else {}
                details = parsed.get("error") or raw_body
            except Exception:
                details = raw_body
            msg = details.strip() if details else f"HTTP {resp.status}"
            return False, f"Ollama generation indisponible: {msg}"
        return True, ""
    except (TimeoutError, socket.timeout):
        return False, "Ollama : delai depasse (chargement du modele)."
    except (ConnectionError, OSError) as e:
        return False, f"Ollama : erreur de connexion ({e!s})."
    except Exception as e:
        return False, f"Ollama : {e!s}"
    finally:
        try:
            conn.close()
        except Exception:
            pass


def check_rscript() -> tuple[bool, str]:
    if shutil.which("Rscript"):
        return True, ""
    return False, "Rscript introuvable dans le PATH (installez R ou activez Conda avec R)."
