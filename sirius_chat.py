from __future__ import annotations

import json
import os
import time
import socket
import http.client
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from langchain_core.messages import AIMessage, HumanMessage, SystemMessage

from sirius_documents import iter_workspace_documents
from sirius_prompt_builder import build_system_prompt
from sirius_pubmed import PubMedClient, build_pubmed_context, looks_like_pubmed_request


_DEFAULT_MODELS = {
    "ollama": "llama3:latest",
    "claude": "claude-3-5-sonnet-latest",
    "deepseek": "deepseek-chat",
}


@dataclass
class ChatConfig:
    provider: str = "ollama"           # "ollama" | "claude" | "deepseek"
    model: str = "llama3:latest"
    temperature: float = 0.2
    host: str | None = None            # used only for ollama
    request_timeout_s: float = 240.0
    num_predict: int = 1024            # max output tokens
    api_key: str | None = None         # overrides env var when provided


class SiriusChat:
    def __init__(self, config: ChatConfig | None = None):
        self.config = config or ChatConfig()
        self.base_url = (self.config.host or os.environ.get("OLLAMA_HOST", "http://127.0.0.1:11434")).rstrip("/")
        self._opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
        self.history: list[HumanMessage | AIMessage | SystemMessage] = []

    def set_provider(self, provider: str, model: str | None = None) -> None:
        """Switch LLM provider. Resets history to avoid cross-provider contamination."""
        provider = (provider or "ollama").lower().strip()
        if provider not in _DEFAULT_MODELS:
            raise ValueError(f"Unknown provider: {provider}")
        self.config.provider = provider
        self.config.model = model or _DEFAULT_MODELS[provider]
        self.history = []

    def _get_api_key(self) -> str:
        if self.config.api_key:
            return self.config.api_key
        env_var = {
            "claude": "ANTHROPIC_API_KEY",
            "deepseek": "DEEPSEEK_API_KEY",
        }.get(self.config.provider, "")
        if not env_var:
            return ""
        key = os.environ.get(env_var, "").strip()
        if not key:
            raise RuntimeError(
                f"Cle API introuvable. Definissez {env_var} avant de lancer l'application."
            )
        return key

    def _split_system_and_messages(self, messages: list) -> tuple[str, list[dict]]:
        """Separate system messages (merged) from user/assistant messages."""
        system_parts: list[str] = []
        chat_msgs: list[dict] = []
        for m in messages:
            if isinstance(m, SystemMessage):
                system_parts.append(str(getattr(m, "content", "")))
            else:
                role = "assistant" if isinstance(m, AIMessage) else "user"
                chat_msgs.append({"role": role, "content": str(getattr(m, "content", ""))})
        return "\n\n".join(system_parts), chat_msgs

    def _to_ollama_messages(self, messages: list) -> list[dict[str, str]]:
        out: list[dict[str, str]] = []
        for m in messages:
            role = "user"
            if isinstance(m, SystemMessage):
                role = "system"
            elif isinstance(m, AIMessage):
                role = "assistant"
            out.append({"role": role, "content": str(getattr(m, "content", ""))})
        return out

    def _parse_base_url(self) -> tuple[str, int, str]:
        from urllib.parse import urlparse
        u = urlparse(self.base_url)
        host = u.hostname or "127.0.0.1"
        port = u.port or (443 if u.scheme == "https" else 11434)
        scheme = u.scheme or "http"
        return scheme, host, port

    def _invoke_with_timeout(self, messages: list) -> str:
        provider = (self.config.provider or "ollama").lower()
        if provider == "ollama":
            return self._invoke_ollama(messages)
        if provider == "claude":
            return self._invoke_claude(messages)
        if provider == "deepseek":
            return self._invoke_deepseek(messages)
        raise RuntimeError(f"Unknown provider: {provider}")

    def _invoke_ollama(self, messages: list) -> str:
        start = time.monotonic()
        payload = json.dumps(
            {
                "model": self.config.model,
                "messages": self._to_ollama_messages(messages),
                "stream": False,
                "options": {
                    "temperature": self.config.temperature,
                    "num_predict": self.config.num_predict,
                },
            }
        ).encode("utf-8")
        scheme, host, port = self._parse_base_url()
        Conn = http.client.HTTPSConnection if scheme == "https" else http.client.HTTPConnection
        conn = Conn(host, port, timeout=self.config.request_timeout_s)
        try:
            conn.request(
                "POST",
                "/api/chat",
                body=payload,
                headers={"Content-Type": "application/json"},
            )
            resp = conn.getresponse()
            raw = resp.read().decode("utf-8", errors="replace")
            if resp.status != 200:
                details = ""
                try:
                    parsed = json.loads(raw) if raw else {}
                    details = parsed.get("error") or raw
                except Exception:
                    details = raw
                raise RuntimeError(
                    f"Ollama HTTP error: {details.strip() if details else f'HTTP {resp.status}'}"
                )
            data = json.loads(raw) if raw else {}
            return data.get("message", {}).get("content", "")
        except (TimeoutError, socket.timeout) as e:
            elapsed_s = time.monotonic() - start
            raise TimeoutError(
                "Ollama n'a pas repondu dans le delai imparti "
                f"(~{elapsed_s:.1f}s ecoulees, limite {self.config.request_timeout_s:.0f}s)."
            ) from e
        except (http.client.RemoteDisconnected, ConnectionError, OSError) as e:
            raise RuntimeError(f"Ollama connection error: {e!s}") from e
        finally:
            try:
                conn.close()
            except Exception:
                pass

    def _invoke_claude(self, messages: list) -> str:
        start = time.monotonic()
        api_key = self._get_api_key()
        system_text, chat_msgs = self._split_system_and_messages(messages)
        if not chat_msgs:
            chat_msgs = [{"role": "user", "content": "(empty)"}]
        payload = json.dumps({
            "model": self.config.model,
            "max_tokens": int(self.config.num_predict),
            "temperature": float(self.config.temperature),
            "system": system_text,
            "messages": chat_msgs,
        }).encode("utf-8")
        conn = http.client.HTTPSConnection(
            "api.anthropic.com", 443, timeout=self.config.request_timeout_s
        )
        try:
            conn.request("POST", "/v1/messages", body=payload, headers={
                "Content-Type": "application/json",
                "x-api-key": api_key,
                "anthropic-version": "2023-06-01",
            })
            resp = conn.getresponse()
            raw = resp.read().decode("utf-8", errors="replace")
            if resp.status != 200:
                details = raw
                try:
                    parsed = json.loads(raw) if raw else {}
                    err = parsed.get("error", {})
                    details = err.get("message") or raw
                except Exception:
                    pass
                raise RuntimeError(f"Claude HTTP {resp.status}: {details[:500]}")
            data = json.loads(raw) if raw else {}
            blocks = data.get("content", [])
            return "".join(b.get("text", "") for b in blocks if b.get("type") == "text")
        except (TimeoutError, socket.timeout) as e:
            elapsed_s = time.monotonic() - start
            raise TimeoutError(
                f"Claude timeout (~{elapsed_s:.1f}s, limite {self.config.request_timeout_s:.0f}s)."
            ) from e
        except (http.client.RemoteDisconnected, ConnectionError, OSError) as e:
            raise RuntimeError(f"Claude connection error: {e!s}") from e
        finally:
            try:
                conn.close()
            except Exception:
                pass

    def _invoke_deepseek(self, messages: list) -> str:
        start = time.monotonic()
        api_key = self._get_api_key()
        system_text, chat_msgs = self._split_system_and_messages(messages)
        if system_text:
            chat_msgs = [{"role": "system", "content": system_text}] + chat_msgs
        if not any(m["role"] == "user" for m in chat_msgs):
            chat_msgs.append({"role": "user", "content": "(empty)"})
        payload = json.dumps({
            "model": self.config.model,
            "messages": chat_msgs,
            "temperature": float(self.config.temperature),
            "max_tokens": int(self.config.num_predict),
            "stream": False,
        }).encode("utf-8")
        conn = http.client.HTTPSConnection(
            "api.deepseek.com", 443, timeout=self.config.request_timeout_s
        )
        try:
            conn.request("POST", "/v1/chat/completions", body=payload, headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {api_key}",
            })
            resp = conn.getresponse()
            raw = resp.read().decode("utf-8", errors="replace")
            if resp.status != 200:
                details = raw
                try:
                    parsed = json.loads(raw) if raw else {}
                    err = parsed.get("error", {})
                    details = err.get("message") or raw
                except Exception:
                    pass
                raise RuntimeError(f"DeepSeek HTTP {resp.status}: {details[:500]}")
            data = json.loads(raw) if raw else {}
            choices = data.get("choices", [])
            if not choices:
                return ""
            return choices[0].get("message", {}).get("content", "")
        except (TimeoutError, socket.timeout) as e:
            elapsed_s = time.monotonic() - start
            raise TimeoutError(
                f"DeepSeek timeout (~{elapsed_s:.1f}s, limite {self.config.request_timeout_s:.0f}s)."
            ) from e
        except (http.client.RemoteDisconnected, ConnectionError, OSError) as e:
            raise RuntimeError(f"DeepSeek connection error: {e!s}") from e
        finally:
            try:
                conn.close()
            except Exception:
                pass

    def reset(self) -> None:
        self.history = []

    def _route_intent(self, user_text: str) -> str:
        t = user_text.lower()
        if looks_like_pubmed_request(user_text):
            return "PUBMED_EVIDENCE"
        if any(k in t for k in ("error", "traceback", "failed", "failure", "exception", "rscript")):
            return "PIPELINE_ERROR_HELP"
        return "MOFA_ASSISTANCE"

    @staticmethod
    def _build_file_index(work_dir: Path, max_entries: int = 500) -> str:
        """Liste recursive de TOUS les fichiers du work_dir (nom, taille, mtime).

        Sert a donner au LLM une vision complete du contenu du dossier,
        y compris les fichiers non textuels (PNG, HDF5, RDS, ...).
        """
        import datetime as _dt
        rows: list[str] = []
        total = 0
        try:
            for p in sorted(work_dir.rglob("*")):
                if not p.is_file():
                    continue
                total += 1
                if len(rows) >= max_entries:
                    continue
                try:
                    st = p.stat()
                    size_kb = st.st_size / 1024.0
                    mtime = _dt.datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M")
                except Exception:
                    size_kb, mtime = 0.0, "?"
                rel = p.relative_to(work_dir).as_posix()
                rows.append(f"  - {rel}  ({size_kb:.1f} KB, {mtime})")
        except Exception as e:
            return f"[ERROR listing work_dir] {e!s}"
        header = f"Total files: {total}"
        if total > len(rows):
            header += f" (showing first {len(rows)})"
        return header + "\n" + "\n".join(rows)

    def ask(
        self,
        user_text: str,
        work_dir: str | None = None,
        max_docs: int = 8,
        pubmed_enabled: bool = False,
        pubmed_email: str = "",
        pubmed_api_key: str = "",
        pubmed_retmax: int = 8,
    ) -> str:
        intent = self._route_intent(user_text)
        blobs = []
        file_index = ""
        pubmed_text = ""

        if work_dir:
            wd = Path(work_dir)
            if wd.is_dir():
                file_index = self._build_file_index(wd)
                blobs = list(iter_workspace_documents(wd, max_files=max_docs))

        if pubmed_enabled and intent == "PUBMED_EVIDENCE":
            try:
                client = PubMedClient(email=pubmed_email, api_key=pubmed_api_key)
                pmids = client.search_pmids(user_text, retmax=pubmed_retmax)
                articles = client.fetch_articles(pmids)
                pubmed_text = build_pubmed_context(articles)
            except Exception:
                pubmed_text = ""

        msgs: list = [
            SystemMessage(
                content=build_system_prompt(
                    intent=intent,
                    has_docs=bool(blobs),
                    has_pubmed=bool(pubmed_text),
                )
            ),
            SystemMessage(
                content=(
                    "Intent router decision: "
                    f"{intent}. Keep answers structured, clinically cautious, and grounded in available evidence."
                )
            ),
        ]

        if work_dir and file_index:
            wd = Path(work_dir)
            msgs.append(SystemMessage(content=(
                f"Index complet du work_dir={wd.resolve().as_posix()} "
                "(tous fichiers, y compris PNG/HDF5/RDS non extraits) :\n\n"
                f"{file_index}"
            )))

        if work_dir and blobs:
            wd = Path(work_dir)
            docs_text = "\n\n".join(
                [f"[FILE] {b.path.name}\n{b.text}" for b in blobs if b.text]
            )
            msgs.append(SystemMessage(content=(
                f"Contenu textuel extrait (PDF/CSV/TSV/TXT/LOG, max {max_docs} fichiers) "
                f"du work_dir={wd.resolve().as_posix()} :\n\n{docs_text}"
            )))

        if pubmed_enabled and intent == "PUBMED_EVIDENCE":
            if pubmed_text:
                msgs.append(
                    SystemMessage(
                        content=(
                            "PubMed grounding context (NCBI E-utilities):\n\n"
                            f"{pubmed_text}\n\n"
                            "When using this context, cite PMID in your final answer (e.g. [PMID:12345678])."
                        )
                    )
                )
            else:
                msgs.append(
                    SystemMessage(
                        content=(
                            "PubMed search was enabled but no matching article was retrieved for this query. "
                            "State this clearly if literature evidence is requested."
                        )
                    )
                )

        # Historique court pour eviter de surcharger
        msgs.extend(self.history[-12:])
        msgs.append(HumanMessage(content=user_text))

        content = self._invoke_with_timeout(msgs)

        if intent == "PUBMED_EVIDENCE" and pubmed_enabled and pubmed_text and "[PMID:" not in content:
            retry_msgs = list(msgs)
            retry_msgs.append(
                SystemMessage(
                    content=(
                        "Quality check failed: include PMID citations for literature-backed claims. "
                        "Regenerate with explicit citations in the Evidence section."
                    )
                )
            )
            retry_content = self._invoke_with_timeout(retry_msgs)
            if "[PMID:" in retry_content:
                content = retry_content

        # Màj historique
        self.history.append(HumanMessage(content=user_text))
        self.history.append(AIMessage(content=content))
        return content
