"""Sirius — Intelligent Agent for Automated Multi-Omics Analysis and Clinical Data integration for precision medicine.

- UI "Mofa-like": boutons pour lancer analyses
- LLM: Ollama (langchain_ollama)
- Documents: collecte depuis work_dir, y compris extraction PDF via pypdf
"""

from __future__ import annotations

import os
import subprocess
import threading
import json as _json
import time as _time
import urllib.request as _urlreq
import urllib.error as _urlerr

for _var in ("NO_PROXY", "no_proxy"):
    _existing = os.environ.get(_var, "")
    _entries = [e.strip() for e in _existing.split(",") if e.strip()]
    for _host in ("localhost", "127.0.0.1"):
        if _host not in _entries:
            _entries.insert(0, _host)
    os.environ[_var] = ",".join(_entries)

from pathlib import Path

import customtkinter as ctk
import tkinter as tk
from tkinter import filedialog, messagebox

from sirius_chat import SiriusChat
from sirius_executor import list_result_files, run_mofa_pipeline
from sirius_health import check_ollama, check_rscript
from sirius_paths import MOFA_PIPELINE_R
from plots_window import PlotsWindow


LOGO_PATH = Path(__file__).resolve().parent / "SIRIUS_dialogs" / "sirius_logo.png"
ICON_PATH = Path(__file__).resolve().parent / "SIRIUS_dialogs" / "sirius.ico"


# ── Windows : forcer un AppUserModelID propre pour que la barre des tâches
#    affiche l'icône de l'app au lieu de celle de python.exe.
#    Doit être appelé AVANT la création de la fenêtre Tk.
def _set_windows_app_user_model_id(app_id: str = "Sirius.MOFA.App.1") -> None:
    if os.name != "nt":
        return
    try:
        import ctypes
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(app_id)
    except Exception:
        pass


_set_windows_app_user_model_id()


ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")


OMIC_TYPES = ["Mutations", "Methylation", "mRNA", "miRNA", "Protein", "Metabolom"]
OMIC_LABELS = {
    "Mutations": "Mutations (DNA)",
    "Methylation": "Methylation (DNA)",
    "mRNA": "mRNA",
    "miRNA": "miRNA",
    "Protein": "Proteomique",
    "Metabolom": "Metabolom",
}

ACCENT = "#00C9A7"
ACCENT2 = "#845EC2"
BG_DARK = "#0F1117"
BG_MID = "#1A1D27"
BG_CARD = "#21253A"
TEXT_DIM = "#8892A4"


class OmicFileCard(ctk.CTkFrame):
    def __init__(self, master, omic_key: str, label: str, **kwargs):
        super().__init__(master, fg_color=BG_CARD, corner_radius=10, **kwargs)
        self.omic_key = omic_key
        self.filepath = tk.StringVar(value="")
        self.enabled = tk.BooleanVar(value=False)

        self.chk = ctk.CTkCheckBox(
            self,
            text=label,
            variable=self.enabled,
            font=ctk.CTkFont(family="Consolas", size=13, weight="bold"),
            text_color=ACCENT,
            fg_color=ACCENT,
            hover_color=ACCENT2,
            command=self._toggle,
        )
        self.chk.grid(row=0, column=0, padx=12, pady=(10, 4), sticky="w")

        self.entry = ctk.CTkEntry(
            self,
            textvariable=self.filepath,
            placeholder_text="No file selected…",
            font=ctk.CTkFont(family="Consolas", size=11),
            fg_color="#12151F",
            text_color="white",
            border_color="#2D3348",
            border_width=1,
            state="disabled",
            width=340,
        )
        self.entry.grid(row=1, column=0, padx=12, pady=(0, 4), sticky="ew")

        self.btn = ctk.CTkButton(
            self,
            text="Browse",
            width=90,
            font=ctk.CTkFont(size=12),
            fg_color="#2D3348",
            hover_color=ACCENT2,
            state="disabled",
            command=self._browse,
        )
        self.btn.grid(row=1, column=1, padx=(4, 12), pady=(0, 4))

        self.lbl_fmt = ctk.CTkLabel(
            self,
            text="Format: .csv .tsv .txt .rds",
            font=ctk.CTkFont(size=10),
            text_color=TEXT_DIM,
        )
        self.lbl_fmt.grid(row=2, column=0, columnspan=2, padx=12, pady=(0, 8), sticky="w")

        self.columnconfigure(0, weight=1)

    def _toggle(self):
        state = "normal" if self.enabled.get() else "disabled"
        self.entry.configure(state=state)
        self.btn.configure(state=state)
        if not self.enabled.get():
            self.filepath.set("")

    def _browse(self):
        path = filedialog.askopenfilename(
            title=f"Select file for {self.omic_key}",
            filetypes=[("Omics data", "*.csv *.tsv *.txt *.rds"), ("All", "*.*")],
        )
        if path:
            self.filepath.set(path)

    def get_value(self):
        if self.enabled.get() and self.filepath.get():
            return (self.omic_key, self.filepath.get())
        return None


class AISiriusMofaApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("AI Sirius MOFA — MOFA2 Desktop")
        self.geometry("1280x860")
        self.configure(fg_color=BG_DARK)
        self.minsize(1100, 760)

        self._logo_image = None
        self._icon_image = None
        self._logo_label = None
        self._try_set_window_icon()

        # Charger les cles API persistees (depuis %APPDATA%/sirius_mofa/api_keys.json)
        self._apply_saved_api_keys_to_env()

        self.chat = SiriusChat()

        self.work_dir = tk.StringVar()
        self.has_metadata = tk.BooleanVar(value=False)
        self.metadata_path = tk.StringVar()
        self.num_factors = tk.IntVar(value=15)
        self.max_r2 = tk.DoubleVar(value=15.0)
        self.max_r2_use_default = tk.BooleanVar(value=True)
        self.pubmed_enabled = tk.BooleanVar(value=False)
        self.pubmed_email = tk.StringVar()
        self.pubmed_api_key = tk.StringVar()
        self.pubmed_retmax = tk.IntVar(value=8)

        self.omic_cards: dict[str, OmicFileCard] = {}

        self._build_ui()
        self.after(500, self._run_health_checks)

        self.work_dir.trace_add("write", lambda *_: self._on_workdir_changed())

    def _build_ui(self):
        hdr = ctk.CTkFrame(self, fg_color=BG_MID, corner_radius=0, height=58)
        hdr.pack(fill="x")

        ctk.CTkLabel(
            hdr,
            text="  AI Sirius MOFA ",
            font=ctk.CTkFont(family="Consolas", size=18, weight="bold"),
            text_color=ACCENT,
        ).pack(side="left", padx=16, pady=14)

        self.btn_viz = ctk.CTkButton(
            hdr,
            text="Visualization",
            font=ctk.CTkFont(size=12, weight="bold"),
            fg_color=ACCENT2,
            hover_color="#6A4BAE",
            width=140,
            height=34,
            command=self._open_plots_window,
        )
        self.btn_viz.pack(side="right", padx=(10, 10), pady=12)

        self.btn_ollama_start = ctk.CTkButton(
            hdr,
            text="Start Ollama",
            font=ctk.CTkFont(size=12, weight="bold"),
            fg_color="#2D3348",
            hover_color=ACCENT2,
            width=130,
            height=34,
            command=self._start_ollama_manually,
        )
        self.btn_ollama_start.pack(side="right", padx=(0, 10), pady=12)

        self.btn_ollama_test = ctk.CTkButton(
            hdr,
            text="Test Ollama",
            font=ctk.CTkFont(size=12, weight="bold"),
            fg_color="#2D3348",
            hover_color=ACCENT2,
            width=120,
            height=34,
            command=self._test_ollama_direct,
        )
        self.btn_ollama_test.pack(side="right", padx=(0, 10), pady=12)

        self._lbl_health = ctk.CTkLabel(
            hdr,
            text="Health: ...",
            font=ctk.CTkFont(size=12),
            text_color=TEXT_DIM,
        )
        self._lbl_health.pack(side="right", padx=16)

        self._try_add_header_logo(hdr)

        body = ctk.CTkFrame(self, fg_color=BG_DARK)
        body.pack(fill="both", expand=True, padx=14, pady=12)
        body.columnconfigure(0, weight=3)
        body.columnconfigure(1, weight=2)
        body.rowconfigure(0, weight=1)

        self.left = ctk.CTkScrollableFrame(body, fg_color=BG_DARK)
        self.left.grid(row=0, column=0, sticky="nsew", padx=(0, 10))

        self.right = ctk.CTkFrame(body, fg_color=BG_CARD, corner_radius=12)
        self.right.grid(row=0, column=1, sticky="nsew")
        self.right.rowconfigure(2, weight=1)
        self.right.columnconfigure(0, weight=1)

        self._build_pipeline_panel()
        self._build_chat_panel()

    def _try_set_window_icon(self):
        # Icône native Windows (.ico) pour la barre de titre / barre des tâches
        if ICON_PATH.is_file():
            try:
                self.iconbitmap(default=str(ICON_PATH))
            except Exception:
                pass
        # Fallback / complément PNG via iconphoto (multi-plateforme)
        if not LOGO_PATH.is_file():
            return
        try:
            try:
                from PIL import Image, ImageTk  # type: ignore

                img = Image.open(LOGO_PATH)
                self._icon_image = ImageTk.PhotoImage(img)
            except Exception:
                self._icon_image = tk.PhotoImage(file=str(LOGO_PATH))
            self.iconphoto(True, self._icon_image)
        except Exception:
            return

    def _try_add_header_logo(self, header: ctk.CTkFrame):
        if not LOGO_PATH.is_file():
            return
        try:
            from PIL import Image  # type: ignore

            img = Image.open(LOGO_PATH)
            self._logo_image = ctk.CTkImage(light_image=img, dark_image=img, size=(34, 34))
            self._logo_label = ctk.CTkLabel(header, text="", image=self._logo_image)
            self._logo_label.pack(side="right", padx=(0, 14), pady=12)
        except Exception:
            try:
                self._logo_image = tk.PhotoImage(file=str(LOGO_PATH))
                self._logo_label = tk.Label(
                    header, image=self._logo_image, bg=BG_MID, bd=0, highlightthickness=0
                )
                self._logo_label.pack(side="right", padx=(0, 14), pady=12)
            except Exception:
                self._logo_image = None

    def _label_section(self, parent, text: str):
        ctk.CTkLabel(
            parent,
            text=text,
            font=ctk.CTkFont(family="Consolas", size=14, weight="bold"),
            text_color="#C8D0E0",
            anchor="w",
        ).pack(anchor="w", pady=(14, 4))

    def _build_pipeline_panel(self):
        self._label_section(self.left, "Working directory")
        frm = ctk.CTkFrame(self.left, fg_color=BG_CARD, corner_radius=10)
        frm.pack(fill="x", pady=(0, 14))

        ctk.CTkEntry(
            frm,
            textvariable=self.work_dir,
            placeholder_text="Select working directory…",
            font=ctk.CTkFont(family="Consolas", size=12),
            fg_color="#12151F",
            border_color="#2D3348",
            border_width=1,
        ).pack(side="left", fill="x", expand=True, padx=(12, 6), pady=12)

        ctk.CTkButton(
            frm,
            text="Browse",
            width=100,
            fg_color=ACCENT2,
            hover_color="#6A4BAE",
            command=self._browse_workdir,
        ).pack(side="right", padx=(0, 12), pady=12)

        self._label_section(self.left, "Metadata (optional)")
        meta = ctk.CTkFrame(self.left, fg_color=BG_CARD, corner_radius=10)
        meta.pack(fill="x", pady=(0, 14))

        chk = ctk.CTkCheckBox(
            meta,
            text="Provide a metadata file",
            variable=self.has_metadata,
            font=ctk.CTkFont(size=13, weight="bold"),
            text_color=ACCENT,
            fg_color=ACCENT,
            hover_color=ACCENT2,
            command=self._toggle_metadata,
        )
        chk.grid(row=0, column=0, columnspan=2, padx=12, pady=(10, 4), sticky="w")

        self.meta_entry = ctk.CTkEntry(
            meta,
            textvariable=self.metadata_path,
            placeholder_text="Metadata file (.csv / .tsv / .rds)…",
            font=ctk.CTkFont(family="Consolas", size=11),
            fg_color="#12151F",
            border_color="#2D3348",
            border_width=1,
            state="disabled",
        )
        self.meta_entry.grid(row=1, column=0, padx=(12, 6), pady=(0, 10), sticky="ew")

        self.meta_btn = ctk.CTkButton(
            meta,
            text="Browse",
            width=100,
            fg_color="#2D3348",
            hover_color=ACCENT2,
            state="disabled",
            command=self._browse_metadata,
        )
        self.meta_btn.grid(row=1, column=1, padx=(0, 12), pady=(0, 10))
        meta.columnconfigure(0, weight=1)

        self._label_section(self.left, "Omics data")
        ctk.CTkLabel(
            self.left,
            text=(
                "Check and select the input files for each data type.\n"
                "Features in rows, samples in columns."
            ),
            font=ctk.CTkFont(size=11),
            text_color=TEXT_DIM,
            justify="left",
        ).pack(anchor="w", pady=(0, 8))

        grid = ctk.CTkFrame(self.left, fg_color="transparent")
        grid.pack(fill="x")

        for i, key in enumerate(OMIC_TYPES):
            card = OmicFileCard(grid, omic_key=key, label=OMIC_LABELS[key])
            r, c = divmod(i, 2)
            card.grid(row=r, column=c, padx=6, pady=5, sticky="nsew")
            self.omic_cards[key] = card

        grid.columnconfigure(0, weight=1)
        grid.columnconfigure(1, weight=1)

        self._label_section(self.left, "Model parameters")
        prm = ctk.CTkFrame(self.left, fg_color=BG_CARD, corner_radius=10)
        prm.pack(fill="x", pady=(0, 14))

        row = ctk.CTkFrame(prm, fg_color="transparent")
        row.pack(fill="x", padx=12, pady=(12, 8))
        ctk.CTkLabel(
            row,
            text="Number of latent factors:",
            font=ctk.CTkFont(size=13, weight="bold"),
            text_color="white",
        ).pack(side="left")

        ctk.CTkButton(
            row,
            text="-",
            width=34,
            fg_color="#2D3348",
            hover_color=ACCENT2,
            command=lambda: self._step_factor(-1),
        ).pack(side="right")
        self.factor_entry = ctk.CTkEntry(
            row,
            textvariable=self.num_factors,
            width=70,
            justify="center",
            font=ctk.CTkFont(family="Consolas", size=16, weight="bold"),
            fg_color="#12151F",
            border_color=ACCENT,
            border_width=1,
            text_color=ACCENT,
        )
        self.factor_entry.pack(side="right", padx=6)
        ctk.CTkButton(
            row,
            text="+",
            width=34,
            fg_color="#2D3348",
            hover_color=ACCENT2,
            command=lambda: self._step_factor(1),
        ).pack(side="right")

        row2 = ctk.CTkFrame(prm, fg_color="transparent")
        row2.pack(fill="x", padx=12, pady=(0, 12))
        ctk.CTkLabel(
            row2,
            text="Max R² for plot_variance_explained:",
            font=ctk.CTkFont(size=13, weight="bold"),
            text_color="white",
        ).pack(side="left")

        self._r2_val_label = ctk.CTkLabel(row2, text="NULL", text_color=TEXT_DIM, width=38)
        self._r2_val_label.pack(side="right")

        self.r2_slider = ctk.CTkSlider(
            row2,
            from_=5,
            to=50,
            variable=self.max_r2,
            fg_color="#2D3348",
            progress_color=ACCENT2,
            button_color=ACCENT,
            state="disabled",
            command=lambda v: self._r2_val_label.configure(text=str(int(v))),
        )
        self.r2_slider.pack(side="right", fill="x", expand=True, padx=(8, 6))

        def _toggle_r2():
            if self.max_r2_use_default.get():
                self.r2_slider.configure(state="disabled")
                self._r2_val_label.configure(text="NULL")
            else:
                self.r2_slider.configure(state="normal")
                self._r2_val_label.configure(text=str(int(self.max_r2.get())))

        ctk.CTkCheckBox(
            row2,
            text="NULL (default)",
            variable=self.max_r2_use_default,
            font=ctk.CTkFont(size=12),
            text_color=TEXT_DIM,
            fg_color=ACCENT2,
            hover_color=ACCENT2,
            border_color="#2D3348",
            command=_toggle_r2,
            width=120,
        ).pack(side="right", padx=(8, 6))

        self._label_section(self.left, "Run")
        run = ctk.CTkFrame(self.left, fg_color="transparent")
        run.pack(fill="x", pady=(6, 10))

        self.btn_run = ctk.CTkButton(
            run,
            text="RUN MOFA PIPELINE",
            font=ctk.CTkFont(size=15, weight="bold"),
            fg_color=ACCENT,
            hover_color="#00A884",
            height=46,
            corner_radius=12,
            command=self._validate_and_run,
        )
        self.btn_run.pack(fill="x")

        self.progress = ctk.CTkProgressBar(self.left, mode="indeterminate", fg_color="#2D3348", progress_color=ACCENT)
        self.progress.pack(fill="x", pady=(0, 8))
        self.progress.set(0)

        self._label_section(self.left, "Console")
        log_frm = ctk.CTkFrame(self.left, fg_color="#0A0D14", corner_radius=10)
        log_frm.pack(fill="both", expand=True, pady=(0, 20))

        self.log_text = tk.Text(
            log_frm,
            height=16,
            bg="#0A0D14",
            fg="#C8FFD4",
            font=("Consolas", 11),
            bd=0,
            relief="flat",
            insertbackground=ACCENT,
            wrap="word",
        )
        self.log_text.pack(fill="both", expand=True, padx=10, pady=10)

        sb = ctk.CTkScrollbar(log_frm, command=self.log_text.yview)
        sb.pack(side="right", fill="y")
        self.log_text.configure(yscrollcommand=sb.set)

        self._log("AI Sirius MOFA is ready. Configure your data and run the pipeline.\n")

    def _build_chat_panel(self):
        title = ctk.CTkLabel(
            self.right,
            text="Sirius chat",
            font=ctk.CTkFont(family="Consolas", size=14, weight="bold"),
            text_color=ACCENT,
        )
        title.grid(row=0, column=0, sticky="w", padx=14, pady=(14, 6))

        self._lbl_docs = ctk.CTkLabel(
            self.right,
            text="Documents: (none)",
            font=ctk.CTkFont(size=11),
            text_color=TEXT_DIM,
            justify="left",
            wraplength=420,
        )
        self._lbl_docs.grid(row=1, column=0, sticky="we", padx=14, pady=(0, 10))

        self.chat_out = tk.Text(
            self.right,
            height=18,
            bg="#0A0D14",
            fg="#d8e9ff",
            font=("Consolas", 11),
            bd=0,
            relief="flat",
            wrap="word",
        )
        self.chat_out.grid(row=2, column=0, sticky="nsew", padx=14)

        # --- LLM provider selection (Ollama / Claude / DeepSeek) ---
        self.llm_provider = tk.StringVar(value="Ollama (llama3)")
        self._llm_provider_map = {
            "Ollama (llama3)":    ("ollama",   "llama3:latest"),
            "Claude (Anthropic)": ("claude",   "claude-3-5-sonnet-latest"),
            "DeepSeek":           ("deepseek", "deepseek-chat"),
        }
        llm_row = ctk.CTkFrame(self.right, fg_color="transparent")
        llm_row.grid(row=3, column=0, sticky="we", padx=14, pady=(8, 0))
        ctk.CTkLabel(
            llm_row, text="LLM:",
            text_color=TEXT_DIM, font=ctk.CTkFont(size=11, weight="bold"),
        ).pack(side="left", padx=(0, 8))
        ctk.CTkComboBox(
            llm_row,
            variable=self.llm_provider,
            values=["Ollama (llama3)", "Claude (Anthropic)", "DeepSeek"],
            width=220,
            command=self._on_llm_provider_change,
            fg_color="#12151F",
            border_color="#2D3348",
            button_color=ACCENT2,
            dropdown_fg_color="#1A1D27",
        ).pack(side="left")
        ctk.CTkButton(
            llm_row, text="API key", width=80,
            fg_color=ACCENT, hover_color="#00A88C",
            command=self._open_api_key_dialog,
        ).pack(side="left", padx=(8, 0))
        self._llm_status = ctk.CTkLabel(
            llm_row, text="(local Ollama)",
            text_color=TEXT_DIM, font=ctk.CTkFont(size=10),
        )
        self._llm_status.pack(side="left", padx=(10, 0))

        pub = ctk.CTkFrame(self.right, fg_color="transparent")
        pub.grid(row=4, column=0, sticky="we", padx=14, pady=(8, 4))
        pub.columnconfigure(1, weight=1)

        ctk.CTkCheckBox(
            pub,
            text="PubMed grounding",
            variable=self.pubmed_enabled,
            font=ctk.CTkFont(size=11, weight="bold"),
            text_color=ACCENT,
            fg_color=ACCENT,
            hover_color=ACCENT2,
        ).grid(row=0, column=0, sticky="w", padx=(0, 8))

        ctk.CTkEntry(
            pub,
            textvariable=self.pubmed_email,
            placeholder_text="NCBI email (recommended)",
            font=ctk.CTkFont(size=11),
            fg_color="#12151F",
            border_color="#2D3348",
            border_width=1,
        ).grid(row=0, column=1, sticky="we", padx=(0, 8))

        ctk.CTkEntry(
            pub,
            textvariable=self.pubmed_retmax,
            placeholder_text="retmax",
            width=70,
            font=ctk.CTkFont(size=11),
            fg_color="#12151F",
            border_color="#2D3348",
            border_width=1,
        ).grid(row=0, column=2, sticky="e", padx=(0, 8))

        ctk.CTkEntry(
            pub,
            textvariable=self.pubmed_api_key,
            placeholder_text="NCBI API key (optional)",
            font=ctk.CTkFont(size=11),
            fg_color="#12151F",
            border_color="#2D3348",
            border_width=1,
        ).grid(row=1, column=1, columnspan=2, sticky="we", pady=(6, 0))

        row = ctk.CTkFrame(self.right, fg_color="transparent")
        row.grid(row=5, column=0, sticky="we", padx=14, pady=14)
        row.columnconfigure(0, weight=1)

        self.chat_in = ctk.CTkEntry(
            row,
            placeholder_text="Ask a question about the results (PDF/CSV/TSV)...",
            font=ctk.CTkFont(size=12),
            fg_color="#12151F",
            border_color="#2D3348",
            border_width=1,
        )
        self.chat_in.grid(row=0, column=0, sticky="we", padx=(0, 8))
        self.chat_in.bind("<Return>", lambda _e: self._send_chat())

        self.btn_send = ctk.CTkButton(
            row,
            text="Send",
            width=90,
            fg_color=ACCENT2,
            hover_color="#6A4BAE",
            command=self._send_chat,
        )
        self.btn_send.grid(row=0, column=1)

    def _run_health_checks(self):
        def run():
            ok_o, msg_o = check_ollama()
            ok_r, msg_r = check_rscript()
            parts = [
                "Ollama OK" if ok_o else f"Ollama NOK: {msg_o}",
                "Rscript OK" if ok_r else f"Rscript NOK: {msg_r}",
            ]
            self.after(0, lambda: self._lbl_health.configure(text=" | ".join(parts)))
        threading.Thread(target=run, daemon=True).start()

    def _start_ollama_manually(self):
        self.btn_ollama_start.configure(state="disabled", text="Starting...")
        self._lbl_health.configure(text="Health: starting Ollama...")

        def run():
            try:
                creationflags = getattr(subprocess, "CREATE_NEW_CONSOLE", 0)
                subprocess.Popen(
                    ["ollama", "serve"],
                    stdin=subprocess.DEVNULL,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    creationflags=creationflags,
                )
            except FileNotFoundError:
                self.after(
                    0,
                    lambda: self._lbl_health.configure(
                        text="Ollama NOK: commande `ollama` introuvable dans le PATH."
                    ),
                )
            except Exception as e:
                self.after(0, lambda: self._lbl_health.configure(text=f"Ollama NOK: {e!s}"))
            finally:
                def finish_ui():
                    self.btn_ollama_start.configure(state="normal", text="Start Ollama")
                    self._run_health_checks()
                    self.after(2500, self._run_health_checks)
                    self.after(6000, self._run_health_checks)

                self.after(0, finish_ui)

        threading.Thread(target=run, daemon=True).start()

    def _test_ollama_direct(self):
        self.btn_ollama_test.configure(state="disabled", text="Testing...")
        base = self.chat.base_url
        model = self.chat.config.model
        self._chat_log(f"\n[TEST] POST {base}/api/chat  model={model}\n")

        def run():
            import http.client as _http
            from urllib.parse import urlparse as _parse
            start = _time.monotonic()
            status = "?"
            detail = ""
            conn = None
            try:
                u = _parse(base)
                host = u.hostname or "127.0.0.1"
                port = u.port or (443 if u.scheme == "https" else 11434)
                Conn = _http.HTTPSConnection if u.scheme == "https" else _http.HTTPConnection
                conn = Conn(host, port, timeout=120)
                payload = _json.dumps(
                    {
                        "model": model,
                        "messages": [{"role": "user", "content": "ping"}],
                        "stream": False,
                        "options": {"num_predict": 16},
                    }
                ).encode("utf-8")
                conn.request(
                    "POST", "/api/chat", body=payload,
                    headers={"Content-Type": "application/json"},
                )
                resp = conn.getresponse()
                raw = resp.read().decode("utf-8", errors="replace")
                if resp.status != 200:
                    status = f"HTTP {resp.status}"
                    detail = raw[:300]
                else:
                    data = _json.loads(raw) if raw else {}
                    content = data.get("message", {}).get("content", "")
                    status = "OK"
                    detail = content.strip()[:200] or "(empty content)"
            except Exception as e:
                status = type(e).__name__
                detail = str(e)[:300]
            finally:
                try:
                    if conn is not None:
                        conn.close()
                except Exception:
                    pass
            elapsed = _time.monotonic() - start

            def done():
                self._chat_log(f"[TEST] {status} in {elapsed:.1f}s :: {detail}\n")
                self.btn_ollama_test.configure(state="normal", text="Test Ollama")

            self.after(0, done)

        threading.Thread(target=run, daemon=True).start()

    def _browse_workdir(self):
        p = filedialog.askdirectory()
        if p:
            self.work_dir.set(p)

    def _on_workdir_changed(self):
        self._refresh_docs_label()
        self._auto_detect_omics()

    def _auto_detect_omics(self):
        wd = self.work_dir.get().strip()
        root = Path(wd) if wd else None
        if not root or not root.is_dir():
            return
        exts = {".csv", ".tsv", ".txt", ".rds"}
        candidates: dict[str, Path] = {}
        for p in root.iterdir():
            if p.is_file() and p.suffix.lower() in exts:
                candidates[p.stem.lower()] = p
        for key, card in self.omic_cards.items():
            match = candidates.get(key.lower())
            if match:
                card.enabled.set(True)
                card.filepath.set(str(match))
                card.entry.configure(state="normal")
                card.btn.configure(state="normal")
                self._log(f"Auto-detected: {key} → {match.name}\n")

    def _toggle_metadata(self):
        state = "normal" if self.has_metadata.get() else "disabled"
        self.meta_entry.configure(state=state)
        self.meta_btn.configure(state=state)

    def _browse_metadata(self):
        p = filedialog.askopenfilename(filetypes=[("Metadata", "*.csv *.tsv *.txt *.rds"), ("All", "*.*")])
        if p:
            self.metadata_path.set(p)

    def _step_factor(self, delta: int):
        val = int(self.num_factors.get()) + delta
        if val >= 1:
            self.num_factors.set(val)

    def _log(self, msg: str):
        self.log_text.insert("end", msg)
        self.log_text.see("end")

    def _chat_log(self, msg: str):
        self.chat_out.insert("end", msg)
        self.chat_out.see("end")

    def _api_keys_config_path(self) -> Path:
        """Local file to persist API keys (user profile, outside the repo)."""
        base = Path(os.environ.get("APPDATA", str(Path.home()))) / "sirius_mofa"
        base.mkdir(parents=True, exist_ok=True)
        return base / "api_keys.json"

    def _load_saved_api_keys(self) -> dict:
        p = self._api_keys_config_path()
        if not p.is_file():
            return {}
        try:
            import json as _json
            return _json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            return {}

    def _save_api_keys_to_disk(self, data: dict) -> None:
        p = self._api_keys_config_path()
        import json as _json
        p.write_text(_json.dumps(data, indent=2), encoding="utf-8")

    def _apply_saved_api_keys_to_env(self):
        """Load persisted keys into os.environ at startup (if any)."""
        saved = self._load_saved_api_keys()
        for k, v in saved.items():
            if v and not os.environ.get(k, "").strip():
                os.environ[k] = str(v)

    def _open_api_key_dialog(self):
        import tkinter as _tk
        dlg = ctk.CTkToplevel(self)
        dlg.title("API keys — Claude / DeepSeek")
        dlg.geometry("520x280")
        dlg.configure(fg_color=BG_DARK)
        dlg.transient(self)
        dlg.grab_set()

        ctk.CTkLabel(
            dlg, text="Cles API (laisser vide pour conserver la valeur existante)",
            text_color=TEXT_DIM, font=ctk.CTkFont(size=11),
        ).pack(pady=(14, 10), padx=16, anchor="w")

        saved = self._load_saved_api_keys()
        cur_claude = os.environ.get("ANTHROPIC_API_KEY", "") or saved.get("ANTHROPIC_API_KEY", "")
        cur_deep   = os.environ.get("DEEPSEEK_API_KEY", "")  or saved.get("DEEPSEEK_API_KEY", "")

        var_claude = _tk.StringVar(value=cur_claude)
        var_deep   = _tk.StringVar(value=cur_deep)
        var_save   = _tk.BooleanVar(value=bool(saved))

        def mk_row(label, var):
            row = ctk.CTkFrame(dlg, fg_color="transparent")
            row.pack(fill="x", padx=16, pady=4)
            ctk.CTkLabel(row, text=label, width=150, anchor="w",
                         text_color="#d8e9ff", font=ctk.CTkFont(size=11)).pack(side="left")
            entry = ctk.CTkEntry(
                row, textvariable=var, show="*",
                fg_color="#12151F", border_color="#2D3348", border_width=1,
            )
            entry.pack(side="left", fill="x", expand=True)
            return entry

        mk_row("ANTHROPIC_API_KEY", var_claude)
        mk_row("DEEPSEEK_API_KEY",  var_deep)

        ctk.CTkCheckBox(
            dlg, text="Enregistrer sur ce poste (fichier local, non chiffre)",
            variable=var_save,
            font=ctk.CTkFont(size=10), text_color=TEXT_DIM,
            fg_color=ACCENT, hover_color=ACCENT2,
        ).pack(pady=(10, 4), padx=16, anchor="w")

        info = ctk.CTkLabel(
            dlg,
            text=f"Stockage local: {self._api_keys_config_path()}",
            text_color=TEXT_DIM, font=ctk.CTkFont(size=9),
        )
        info.pack(padx=16, anchor="w")

        btns = ctk.CTkFrame(dlg, fg_color="transparent")
        btns.pack(fill="x", padx=16, pady=14)

        def on_save():
            k_claude = var_claude.get().strip()
            k_deep   = var_deep.get().strip()
            # Applique immediatement pour la session
            if k_claude:
                os.environ["ANTHROPIC_API_KEY"] = k_claude
            if k_deep:
                os.environ["DEEPSEEK_API_KEY"] = k_deep
            # Persistance optionnelle
            if var_save.get():
                data = self._load_saved_api_keys()
                if k_claude: data["ANTHROPIC_API_KEY"] = k_claude
                if k_deep:   data["DEEPSEEK_API_KEY"]  = k_deep
                try:
                    self._save_api_keys_to_disk(data)
                    self._chat_log(f"\n[INFO] Cles API sauvegardees : {self._api_keys_config_path()}\n")
                except Exception as e:
                    self._chat_log(f"\n[ERROR] Sauvegarde cles : {e}\n")
            else:
                self._chat_log("\n[INFO] Cles API definies pour cette session uniquement.\n")
            # Re-apply current provider to pick up new key
            self._on_llm_provider_change(self.llm_provider.get())
            dlg.destroy()

        ctk.CTkButton(btns, text="Save", width=100,
                      fg_color=ACCENT, hover_color="#00A88C",
                      command=on_save).pack(side="right", padx=(8, 0))
        ctk.CTkButton(btns, text="Cancel", width=100,
                      fg_color="#3A3F52", hover_color="#4A4F62",
                      command=dlg.destroy).pack(side="right")

    def _on_llm_provider_change(self, value: str):
        """Switch LLM backend when user changes the dropdown."""
        provider, model = self._llm_provider_map.get(value, ("ollama", "llama3:latest"))
        try:
            self.chat.set_provider(provider, model=model)
            status_map = {
                "ollama":   "(local Ollama)",
                "claude":   "(API · ANTHROPIC_API_KEY)",
                "deepseek": "(API · DEEPSEEK_API_KEY)",
            }
            self._llm_status.configure(text=status_map.get(provider, ""))
            self._chat_log(f"\n[INFO] Provider={value}  model={model}\n")
            if provider in ("claude", "deepseek"):
                env_var = "ANTHROPIC_API_KEY" if provider == "claude" else "DEEPSEEK_API_KEY"
                if not os.environ.get(env_var, "").strip():
                    self._chat_log(
                        f"[WARN] Variable {env_var} non definie. "
                        f"Definissez-la avant d'envoyer un message.\n"
                    )
        except Exception as e:
            self._chat_log(f"\n[ERROR] Provider switch: {e}\n")

    def _refresh_docs_label(self):
        wd = self.work_dir.get().strip()
        if not wd or not os.path.isdir(wd):
            self._lbl_docs.configure(text="Documents: (none)")
            return
        files = list_result_files(wd, max_files=50)
        if not files:
            self._lbl_docs.configure(text="Documents: (none)")
            return
        names = "\n".join([f"- {p.name}" for p in files[:12]])
        more = "" if len(files) <= 12 else f"\n(+{len(files)-12} more)"
        self._lbl_docs.configure(text=f"Detected documents:\n{names}{more}")

    def _open_plots_window(self):
        wd = self.work_dir.get().strip()
        if not wd or not os.path.isdir(wd):
            messagebox.showerror("Error", "Select a valid working directory.")
            return
        rds = Path(wd) / "model_with_metadata.rds"
        hdf = Path(wd) / "model.hdf5"
        mp = str(rds) if rds.is_file() else (str(hdf) if hdf.is_file() else "")
        PlotsWindow(self, work_dir=wd, model_path=mp)

    def _validate_and_run(self):
        errors: list[str] = []
        work_dir = self.work_dir.get().strip()
        if not work_dir or not os.path.isdir(work_dir):
            errors.append("- Select a valid working directory.")

        selected_omics: dict[str, str] = {}
        for key, card in self.omic_cards.items():
            val = card.get_value()
            if val:
                selected_omics[val[0]] = val[1]
        if len(selected_omics) < 2:
            errors.append("- Select at least 2 omics views.")

        for name, path in selected_omics.items():
            if not os.path.isfile(path):
                errors.append(f"- File not found for {name}: {path}")

        try:
            nf = int(self.num_factors.get())
            if nf < 1:
                raise ValueError
        except Exception:
            errors.append("- The number of factors must be a positive integer.")
            nf = 15

        has_meta = bool(self.has_metadata.get())
        meta_path = self.metadata_path.get().strip() if has_meta else ""
        if has_meta and (not meta_path or not os.path.isfile(meta_path)):
            errors.append("- Metadata is enabled but the file was not found.")

        if not MOFA_PIPELINE_R.is_file():
            errors.append(f"- R script not found: {MOFA_PIPELINE_R}")

        if errors:
            messagebox.showerror("Errors", "\n".join(errors))
            return

        if not messagebox.askyesno(
            "Confirm",
            f"Working directory: {work_dir}\nViews: {', '.join(selected_omics.keys())}\nFactors: {nf}\nMax R2: {'NULL (default)' if self.max_r2_use_default.get() else f'{self.max_r2.get():.0f}'}\n\nRun the pipeline?",
        ):
            return

        self.btn_run.configure(state="disabled", text="Pipeline running...")
        self.progress.start()
        self.log_text.delete("1.0", "end")

        t = threading.Thread(
            target=self._run_pipeline_thread,
            args=(work_dir, selected_omics, nf, None if self.max_r2_use_default.get() else float(self.max_r2.get()), has_meta, meta_path),
            daemon=True,
        )
        t.start()

    def _run_pipeline_thread(self, work_dir: str, selected_omics: dict[str, str], nf: int, max_r2: float | None, has_meta: bool, meta_path: str):
        code, log = run_mofa_pipeline(
            work_dir=work_dir,
            selected_omics=selected_omics,
            num_factors=nf,
            max_r2=max_r2,
            has_metadata=has_meta,
            metadata_file=meta_path,
        )

        def done():
            self.progress.stop()
            self.btn_run.configure(state="normal", text="RUN MOFA PIPELINE")
            if code == 0:
                self.progress.set(1)
                self._log("\nOK: pipeline finished.\n")
                self._refresh_docs_label()
                messagebox.showinfo(
                    "Pipeline finished",
                    "The MOFA pipeline has finished. You can now use the chat to interpret the results.",
                )
            else:
                self.progress.set(0)
                self._log(f"\nERROR: pipeline failed (code {code}).\n")
                messagebox.showerror("Pipeline error", f"The pipeline failed (code {code}). See the console.")

        self.after(0, self._log, log)
        self.after(0, done)

    def _send_chat(self):
        text = self.chat_in.get().strip()
        if not text:
            return
        self.chat_in.delete(0, "end")
        wd = self.work_dir.get().strip()
        self._chat_log(f"\nYOU> {text}\n")
        self.btn_send.configure(state="disabled")

        def run_llm():
            try:
                try:
                    retmax = int(self.pubmed_retmax.get())
                except Exception:
                    retmax = 8
                answer = self.chat.ask(
                    text,
                    work_dir=wd if os.path.isdir(wd) else None,
                    pubmed_enabled=bool(self.pubmed_enabled.get()),
                    pubmed_email=self.pubmed_email.get().strip(),
                    pubmed_api_key=self.pubmed_api_key.get().strip(),
                    pubmed_retmax=retmax,
                )
            except TimeoutError:
                provider = (self.chat.config.provider or "ollama").lower()
                try:
                    # Retry en conservant le work_dir, les documents, l'historique
                    # et PubMed : les documents sont legers (figures exclues, PDF/TXT
                    # plafonnes), donc ce n'est pas le contexte qu'il faut sacrifier.
                    answer = self.chat.ask(
                        text,
                        work_dir=wd if os.path.isdir(wd) else None,
                        pubmed_enabled=bool(self.pubmed_enabled.get()),
                        pubmed_email=self.pubmed_email.get().strip(),
                        pubmed_api_key=self.pubmed_api_key.get().strip(),
                        pubmed_retmax=retmax,
                    )
                    answer = (
                        "[RECOVERY] Premiere requete en timeout, reponse obtenue au second essai "
                        "(contexte workspace preserve).\n\n"
                        f"{answer}"
                    )
                except Exception:
                    if provider == "ollama":
                        answer = (
                            "[ERROR] Ollama n'a pas repondu a temps, meme apres retry. "
                            "Cliquez sur 'Start Ollama' puis reessayez."
                        )
                    else:
                        answer = (
                            f"[ERROR] Le fournisseur '{provider}' n'a pas repondu a temps, "
                            "meme apres retry. Reessayez dans un instant ou augmentez le "
                            "timeout (ChatConfig.request_timeout_s)."
                        )
            except Exception as e:
                msg = str(e)
                lower = msg.lower()
                if "not found" in lower and "model" in lower:
                    answer = (
                        "[ERROR] Le modele Ollama configure est introuvable. "
                        "Verifiez `llama3:latest` via `ollama list` et `ollama pull llama3:latest`."
                    )
                elif "connection" in lower or "connect" in lower or "refused" in lower:
                    answer = (
                        "[ERROR] Connexion Ollama impossible. "
                        "Verifiez `OLLAMA_HOST` (ex: http://127.0.0.1:11434) et que `ollama serve` est actif."
                    )
                else:
                    answer = f"[ERROR] {msg}"

            def done():
                self._chat_log(f"AI> {answer}\n")
                self.btn_send.configure(state="normal")

            self.after(0, done)

        threading.Thread(target=run_llm, daemon=True).start()


if __name__ == "__main__":
    app = AISiriusMofaApp()
    app.mainloop()
