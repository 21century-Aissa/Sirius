"""Visualization window (plots / survival / GSEA) for AI Sirius MOFA."""

from __future__ import annotations

import os
import threading
from pathlib import Path

import customtkinter as ctk
import tkinter as tk
from tkinter import messagebox

from sirius_executor import run_mofa_gsea, run_mofa_plots, run_mofa_survival


ACCENT = "#00C9A7"
ACCENT2 = "#845EC2"
BG_DARK = "#0F1117"
BG_MID = "#1A1D27"
BG_CARD = "#21253A"
BG_INPUT = "#12151F"
TEXT_DIM = "#8892A4"


class PlotsWindow(ctk.CTkToplevel):

    _PLOT_FIELDS: dict[str, list] = {
        "factor_values": [
            ("factors",    "Factors (e.g. 1,2)",    "entry", "1,2"),
        ],
        "top_weights": [
            ("view",      "View (name, CSV, or ALL)", "entry", "ALL"),
            ("factor",    "Factor",       "entry", "1"),
            ("nfeatures", "N features",   "entry", "25"),
            ("scale",     "scale",        "combo", "TRUE", ["TRUE", "FALSE"]),
        ],
        "factor_violin": [
            ("factors",  "Factors (e.g. 1,2)", "entry", "1,2"),
            ("color_by", "color_by",           "entry", ""),
        ],
        "feature_weights": [
            ("view",      "View (omics)", "entry", "mRNA"),
            ("factor",    "Factor",       "entry", "1"),
            ("nfeatures", "N features",   "entry", "25"),
        ],
        "data_scatter": [
            ("view",     "View (omics)",          "entry", "mRNA"),
            ("factor",   "Factor",                "entry", "1"),
            ("features", "Features (comma-sep)",  "entry", ""),
            ("sign",     "sign",                  "combo", "positive", ["positive", "negative"]),
            ("color_by", "color_by",              "entry", ""),
            ("y_label",  "y_label",               "entry", ""),
        ],
        "data_heatmap": [
            ("view",          "View (omics)",  "entry", "mRNA"),
            ("factor",        "Factor",         "entry", "1"),
            ("features",      "N features",     "entry", "20"),
            ("cluster_rows",  "cluster_rows",   "combo", "TRUE",   ["TRUE", "FALSE"]),
            ("cluster_cols",  "cluster_cols",   "combo", "TRUE",   ["TRUE", "FALSE"]),
            ("show_rownames", "show_rownames",  "combo", "TRUE",   ["TRUE", "FALSE"]),
            ("show_colnames", "show_colnames",  "combo", "FALSE",  ["TRUE", "FALSE"]),
            ("scale_mode",    "scale",           "combo", "row",    ["row", "column", "none"]),
        ],
        "factors_scatter": [
            ("factor_x", "Factor X",  "entry", "1"),
            ("factor_y", "Factor Y",  "entry", "2"),
            ("color_by", "color_by",  "entry", ""),
            ("shape_by", "shape_by",  "entry", ""),
        ],
        "variance_by_factor": [
            ("min_r2",  "min_r2 (optional)",  "entry", ""),
            ("max_r2",  "max_r2 (optional)",  "entry", ""),
            ("x",       "X axis",  "combo", "view",   ["view", "factor", "group"]),
            ("y",       "Y axis",  "combo", "factor", ["factor", "view", "group"]),
            ("palette", "Color palette", "combo", "", ["", "Blues", "Reds", "YlOrRd", "Greens", "Purples", "RdYlBu"]),
        ],
        "variance_total": [
            ("factors",         "Factors (ex: 1,2,3 or all)", "entry", "all"),
            ("label_threshold", "Label threshold (%)",        "entry", "1"),
            ("palette", "Color palette", "combo", "", [
                "",
                "Oranges", "Reds", "Blues", "Greens", "Purples", "Greys",
                "YlOrRd", "YlOrBr", "OrRd",
                "YlGnBu", "BuGn", "PuBu", "GnBu",
                "BuPu", "PuRd", "RdPu",
                "RdYlBu", "RdYlGn", "Spectral",
            ]),
        ],
        "variance_detailed": [
            ("min_r2",  "min_r2 (optional)",  "entry", ""),
            ("max_r2",  "max_r2 (optional)",  "entry", ""),
            ("x",       "X axis",  "combo", "view",   ["view", "factor", "group"]),
            ("y",       "Y axis",  "combo", "factor", ["factor", "view", "group"]),
            ("palette", "Color palette", "combo", "", ["", "Blues", "Reds", "YlOrRd", "Greens", "Purples", "RdYlBu"]),
        ],
        "factor_correlation": [],
        "data_overview":      [],
        "enrichment_heatmap": [
            ("view",   "View (omics)", "entry", "mRNA"),
            ("factor", "Factor",       "entry", "1"),
        ],
        "weights_overview": [
            ("view",   "View (omics)", "entry", "mRNA"),
            ("factor", "Factor",       "entry", "1"),
        ],
        "dimred": [
            ("method",      "Method",      "combo", "UMAP", ["UMAP", "TSNE"]),
            ("color_by",    "color_by",    "entry", ""),
            ("n_neighbors", "n_neighbors", "entry", "15"),
            ("min_dist",    "min_dist",    "entry", "0.1"),
        ],
        "association_analysis": [
            ("covariates", "Covariates (comma-sep, e.g. Gender,age,died)", "entry", ""),
            ("plot",       "Plot type", "combo", "log_pval", ["log_pval", "r"]),
        ],
    }

    def __init__(self, master, work_dir: str, model_path: str):
        super().__init__(master)
        self.title("MOFA — Visualization / Survival / GSEA")
        self.geometry("1400x720")
        self.configure(fg_color=BG_DARK)

        self.work_dir_var = tk.StringVar(value=work_dir)
        self.model_path_var = tk.StringVar(value=model_path)
        self._preview_image: object = None
        self._surv_preview_image: object = None
        self._gsea_preview_image: object = None
        self._last_plot_params: dict | None = None
        self._last_gsea_plot_params: dict | None = None
        self._last_surv_params: dict | None = None

        self._build_ui()

    def _build_ui(self):
        hdr = ctk.CTkFrame(self, fg_color=BG_MID, corner_radius=0, height=56)
        hdr.pack(fill="x")

        ctk.CTkLabel(
            hdr,
            text="  MOFA Visualization",
            font=ctk.CTkFont(family="Consolas", size=16, weight="bold"),
            text_color=ACCENT,
        ).pack(side="left", padx=16, pady=14)

        body = ctk.CTkFrame(self, fg_color=BG_DARK)
        body.pack(fill="both", expand=True, padx=14, pady=14)
        body.rowconfigure(1, weight=1)
        body.columnconfigure(0, weight=1)

        top = ctk.CTkFrame(body, fg_color=BG_CARD, corner_radius=10)
        top.grid(row=0, column=0, sticky="ew", pady=(0, 12))
        top.columnconfigure(1, weight=1)

        ctk.CTkLabel(top, text="Work dir:", text_color=TEXT_DIM).grid(row=0, column=0, padx=12, pady=(12, 6), sticky="w")
        ctk.CTkEntry(
            top,
            textvariable=self.work_dir_var,
            font=ctk.CTkFont(family="Consolas", size=11),
            fg_color=BG_INPUT,
            border_color="#2D3348",
            border_width=1,
        ).grid(row=0, column=1, padx=12, pady=(12, 6), sticky="ew")

        ctk.CTkLabel(top, text="MOFA Model:", text_color=TEXT_DIM).grid(row=1, column=0, padx=12, pady=(0, 12), sticky="w")
        ctk.CTkEntry(
            top,
            textvariable=self.model_path_var,
            font=ctk.CTkFont(family="Consolas", size=11),
            fg_color=BG_INPUT,
            border_color="#2D3348",
            border_width=1,
        ).grid(row=1, column=1, padx=12, pady=(0, 12), sticky="ew")

        tabs = ctk.CTkTabview(body, fg_color=BG_DARK)
        tabs.grid(row=1, column=0, sticky="nsew")

        tab_plots = tabs.add("Plots")
        tab_surv = tabs.add("Survival")
        tab_gsea = tabs.add("GSEA")

        self._build_tab_plots(tab_plots)
        self._build_tab_survival(tab_surv)
        self._build_tab_gsea(tab_gsea)

        self.log = tk.Text(
            body,
            height=10,
            bg="#0A0D14",
            fg="#C8FFD4",
            font=("Consolas", 10),
            bd=0,
            relief="flat",
            wrap="word",
        )
        self.log.grid(row=2, column=0, sticky="nsew", pady=(12, 0))

    def _log(self, msg: str):
        self.log.insert("end", msg)
        self.log.see("end")

    def _validate_paths(self) -> tuple[bool, str]:
        wd = self.work_dir_var.get().strip()
        mp = self.model_path_var.get().strip()
        if not wd or not os.path.isdir(wd):
            return False, "Invalid working directory."
        if not mp or not os.path.isfile(mp):
            return False, "Model .hdf5 not found."
        return True, ""

    def _run_in_thread(self, fn, *args, **kwargs):
        ok, err = self._validate_paths()
        if not ok:
            messagebox.showerror("Error", err)
            return

        def run():
            code, out = fn(*args, **kwargs)

            def done():
                self._log(out + "\n")
                if code == 0:
                    self._log("\nOK\n")
                else:
                    self._log(f"\nERROR code={code}\n")

            self.after(0, done)

        threading.Thread(target=run, daemon=True).start()

    def _build_tab_plots(self, tab):
        tab.configure(fg_color=BG_DARK)

        self.plot_type = tk.StringVar(value="factor_values")
        self.plot_out = tk.StringVar(value="plot")
        self._plot_vars: dict[str, tk.StringVar] = {}

        # 2-column layout: controls left, preview right
        main = ctk.CTkFrame(tab, fg_color=BG_DARK)
        main.pack(fill="both", expand=True)
        main.columnconfigure(0, weight=2)
        main.columnconfigure(1, weight=3)
        main.rowconfigure(0, weight=1)

        left = ctk.CTkFrame(main, fg_color=BG_DARK)
        left.grid(row=0, column=0, sticky="nsew", padx=(0, 6))

        right = ctk.CTkFrame(main, fg_color=BG_CARD, corner_radius=10)
        right.grid(row=0, column=1, sticky="nsew")

        # --- Controls (left) ---
        top_card = ctk.CTkFrame(left, fg_color=BG_CARD, corner_radius=10)
        top_card.pack(fill="x", padx=8, pady=(8, 4))
        top_card.columnconfigure(1, weight=1)

        ctk.CTkLabel(top_card, text="plot_type", text_color=TEXT_DIM).grid(row=0, column=0, padx=12, pady=(12, 6), sticky="w")
        ctk.CTkComboBox(
            top_card,
            variable=self.plot_type,
            values=list(self._PLOT_FIELDS.keys()),
            fg_color=BG_INPUT,
            border_color="#2D3348",
            button_color=ACCENT2,
            dropdown_fg_color=BG_MID,
            command=lambda _: self._refresh_plot_fields(),
        ).grid(row=0, column=1, padx=12, pady=(12, 6), sticky="ew")

        ctk.CTkLabel(top_card, text="out_name", text_color=TEXT_DIM).grid(row=1, column=0, padx=12, pady=(0, 12), sticky="w")
        ctk.CTkEntry(
            top_card,
            textvariable=self.plot_out,
            font=ctk.CTkFont(family="Consolas", size=11),
            fg_color=BG_INPUT,
            border_color="#2D3348",
            border_width=1,
        ).grid(row=1, column=1, padx=12, pady=(0, 12), sticky="ew")

        self._dyn_card = ctk.CTkFrame(left, fg_color=BG_CARD, corner_radius=10)
        self._dyn_card.pack(fill="x", padx=8, pady=(0, 4))
        self._dyn_card.columnconfigure(1, weight=1)

        btn_row = ctk.CTkFrame(left, fg_color=BG_DARK)
        btn_row.pack(fill="x", padx=8, pady=(6, 8))
        btn_row.columnconfigure(0, weight=1)
        btn_row.columnconfigure(1, weight=1)

        self._btn_preview = ctk.CTkButton(
            btn_row,
            text="Preview",
            fg_color=ACCENT2,
            hover_color="#6B4DB0",
            height=40,
            command=self._run_preview,
        )
        self._btn_preview.grid(row=0, column=0, padx=(0, 4), sticky="ew")

        self._btn_save_pdf = ctk.CTkButton(
            btn_row,
            text="Save PDF",
            fg_color=ACCENT,
            hover_color="#00A884",
            height=40,
            state="disabled",
            command=self._save_pdf,
        )
        self._btn_save_pdf.grid(row=0, column=1, padx=(4, 0), sticky="ew")

        # --- Preview panel (right) ---
        _hdr = ctk.CTkFrame(right, fg_color="transparent")
        _hdr.pack(fill="x", padx=8, pady=(10, 2))
        ctk.CTkLabel(
            _hdr,
            text="Preview",
            font=ctk.CTkFont(size=12, weight="bold"),
            text_color=ACCENT,
        ).pack(side="left")
        ctk.CTkButton(
            _hdr, text="✕ Clear", width=70, height=24,
            fg_color=BG_MID, hover_color="#2D3348", text_color=TEXT_DIM,
            font=ctk.CTkFont(size=11),
            command=self._clear_preview,
        ).pack(side="right")

        self._preview_label = ctk.CTkLabel(
            right,
            text="No preview yet.\nGenerate a plot to see it here.",
            text_color=TEXT_DIM,
            font=ctk.CTkFont(size=11),
        )
        self._preview_label.pack(fill="both", expand=True, padx=8, pady=(4, 12))

        self._refresh_plot_fields()

    def _refresh_plot_fields(self):
        for widget in self._dyn_card.winfo_children():
            widget.destroy()
        self._plot_vars.clear()

        fields = self._PLOT_FIELDS.get(self.plot_type.get(), [])

        if not fields:
            ctk.CTkLabel(
                self._dyn_card,
                text="No additional parameters required for this plot type.",
                text_color=TEXT_DIM,
                font=ctk.CTkFont(size=11),
            ).grid(row=0, column=0, columnspan=2, padx=12, pady=12, sticky="w")
            return

        for i, field in enumerate(fields):
            key, label, ftype, default = field[0], field[1], field[2], field[3]
            options = field[4] if len(field) > 4 else None
            pady_top = 12 if i == 0 else 0
            pady_bot = 12 if i == len(fields) - 1 else 6

            ctk.CTkLabel(self._dyn_card, text=label, text_color=TEXT_DIM).grid(
                row=i, column=0, padx=12, pady=(pady_top, pady_bot), sticky="w"
            )
            var = tk.StringVar(value=default)
            self._plot_vars[key] = var

            if ftype == "combo" and options:
                ctk.CTkComboBox(
                    self._dyn_card,
                    variable=var,
                    values=options,
                    fg_color=BG_INPUT,
                    border_color="#2D3348",
                    button_color=ACCENT2,
                    dropdown_fg_color=BG_MID,
                ).grid(row=i, column=1, padx=12, pady=(pady_top, pady_bot), sticky="ew")
            else:
                ctk.CTkEntry(
                    self._dyn_card,
                    textvariable=var,
                    font=ctk.CTkFont(family="Consolas", size=11),
                    fg_color=BG_INPUT,
                    border_color="#2D3348",
                    border_width=1,
                ).grid(row=i, column=1, padx=12, pady=(pady_top, pady_bot), sticky="ew")

    def _collect_params(self) -> dict:
        return {
            "wd":        self.work_dir_var.get().strip(),
            "mp":        self.model_path_var.get().strip(),
            "plot_type": self.plot_type.get().strip(),
            "out_name":  self.plot_out.get().strip() or "plot",
            "extra":     {k: v.get().strip() for k, v in self._plot_vars.items()},
        }

    def _run_preview(self):
        ok, err = self._validate_paths()
        if not ok:
            messagebox.showerror("Error", err)
            return

        p = self._collect_params()
        self._last_plot_params = p
        self._btn_preview.configure(state="disabled")
        self._btn_save_pdf.configure(state="disabled")
        self._preview_label.configure(text="Generating preview...", image=None)
        self._preview_image = None

        def run():
            extra = {**p["extra"], "preview_only": "TRUE"}
            code, out = run_mofa_plots(
                work_dir=p["wd"], model_path=p["mp"],
                plot_type=p["plot_type"], out_name=p["out_name"],
                extra_params=extra,
            )

            def done():
                self._log(out + "\n")
                if code == 0:
                    self._log("\nPreview OK — click Save PDF to export.\n")
                    self._load_preview(p["wd"], p["out_name"])
                    self._btn_save_pdf.configure(state="normal")
                else:
                    self._log(f"\nERROR code={code}\n")
                    self._preview_label.configure(text="Plot failed. See log.", image=None)
                self._btn_preview.configure(state="normal")

            self.after(0, done)

        threading.Thread(target=run, daemon=True).start()

    def _save_pdf(self):
        if self._last_plot_params is None:
            return
        p = dict(self._last_plot_params)
        p["out_name"] = self.plot_out.get().strip() or "plot"
        self._btn_save_pdf.configure(state="disabled")
        self._log("\nSaving PDF...\n")

        def run():
            code, out = run_mofa_plots(
                work_dir=p["wd"], model_path=p["mp"],
                plot_type=p["plot_type"], out_name=p["out_name"],
                extra_params=p["extra"],
            )

            def done():
                self._log(out + "\n")
                if code == 0:
                    self._log(f"\n✓ PDF saved: {p['wd']}/{p['out_name']}.pdf\n")
                else:
                    self._log(f"\nERROR code={code}\n")
                self._btn_save_pdf.configure(state="normal")

            self.after(0, done)

        threading.Thread(target=run, daemon=True).start()

    def _clear_preview(self):
        self._preview_image = None
        self._safe_configure_label(self._preview_label, text="No preview yet.\nGenerate a plot to see it here.", image=None)
        self._btn_save_pdf.configure(state="disabled")

    def _clear_surv_preview(self):
        self._surv_preview_image = None
        self._safe_configure_label(self._surv_preview_label, text="No preview yet.\nRun survival analysis to see it here.", image=None)
        if hasattr(self, "_btn_surv_save_pdf"):
            self._btn_surv_save_pdf.configure(state="disabled")

    def _clear_gsea_preview(self):
        self._gsea_preview_image = None
        self._safe_configure_label(self._gsea_preview_label, text="No preview yet.\nRun GSEA to see it here.", image=None)
        if hasattr(self, "_btn_gsea_save_pdf"):
            self._btn_gsea_save_pdf.configure(state="disabled")

    def _clear_label(self, label, default_text: str):
        self._safe_configure_label(label, text=default_text, image=None)

    def _safe_configure_label(self, label, *, text=None, image=None):
        """Configure a CTkLabel robustly, bypassing stale pyimage references."""
        # Access underlying tk.Label directly to bypass CTkLabel's state machine
        inner = getattr(label, "_label", None)
        try:
            if image is not None:
                # Retrieve tk PhotoImage from CTkImage
                if hasattr(image, "_get_current_image"):
                    tk_img = image._get_current_image(ctk.get_appearance_mode())
                else:
                    tk_img = image
                if inner is not None:
                    inner.configure(image=tk_img, text=text if text is not None else "")
                label._image = image
                label._text = text if text is not None else ""
            else:
                # Clear image
                if inner is not None:
                    inner.configure(image="", text=text if text is not None else "")
                label._image = None
                label._text = text if text is not None else ""
        except Exception as e:
            # Final fallback: just try inner with text only
            try:
                if inner is not None:
                    inner.configure(text=text if text is not None else "")
            except Exception:
                pass

    def _load_preview(self, work_dir: str, out_name: str, label=None):
        if label is None:
            label = self._preview_label
        png_path = str(Path(work_dir) / f"{out_name}_preview.png")
        if not os.path.isfile(png_path):
            self._safe_configure_label(label, text="Preview not available.", image=None)
            return
        try:
            from PIL import Image, ImageTk
            with Image.open(png_path) as raw:
                img = raw.copy()
            self.update_idletasks()
            w = label.winfo_width()
            h = label.winfo_height()
            w = w if w > 50 else 580
            h = h if h > 50 else 400
            img.thumbnail((w - 16, h - 16), Image.LANCZOS)
            tk_img = ImageTk.PhotoImage(img)
            # Keep strong reference
            if label is self._preview_label:
                self._preview_image = tk_img
            elif label is self._surv_preview_label:
                self._surv_preview_image = tk_img
            else:
                self._gsea_preview_image = tk_img
            label._cascade_img_ref = tk_img
            # Directly configure inner tk.Label
            inner = getattr(label, "_label", None)
            if inner is not None:
                inner.configure(image=tk_img, text="")
                label._text = ""
        except Exception as e:
            self._safe_configure_label(label, text=f"Preview error: {e}", image=None)

    def _build_tab_survival(self, tab):
        tab.configure(fg_color=BG_DARK)

        self.surv_analysis = tk.StringVar(value="Plot Hazard ratios")
        self._surv_analysis_map = {
            "Plot Hazard ratios": "cox_hr",
            "Cox summary (bubble)": "cox_summary",
            "Kaplan Meier plot": "kaplan_meier",
            "Kaplan Meier (all factors)": "kaplan_meier_all",
        }
        self.surv_out = tk.StringVar(value="survival")
        self.surv_time = tk.StringVar(value="time")
        self.surv_event = tk.StringVar(value="event")

        pane = ctk.CTkFrame(tab, fg_color=BG_DARK)
        pane.pack(fill="both", expand=True)
        pane.columnconfigure(0, weight=1)
        pane.columnconfigure(1, weight=2)
        pane.rowconfigure(0, weight=1)

        left = ctk.CTkFrame(pane, fg_color=BG_DARK)
        left.grid(row=0, column=0, sticky="nsew", padx=(0, 6))
        right = ctk.CTkFrame(pane, fg_color=BG_CARD, corner_radius=10)
        right.grid(row=0, column=1, sticky="nsew")

        card = ctk.CTkFrame(left, fg_color=BG_CARD, corner_radius=10)
        card.pack(fill="x", padx=8, pady=8)
        card.columnconfigure(1, weight=1)

        ctk.CTkLabel(card, text="analysis", text_color=TEXT_DIM).grid(row=0, column=0, padx=12, pady=(12, 6), sticky="w")
        ctk.CTkComboBox(
            card,
            variable=self.surv_analysis,
            values=["Plot Hazard ratios", "Cox summary (bubble)", "Kaplan Meier plot", "Kaplan Meier (all factors)"],
            fg_color=BG_INPUT,
            border_color="#2D3348",
            button_color=ACCENT2,
            dropdown_fg_color=BG_MID,
        ).grid(row=0, column=1, padx=12, pady=(12, 6), sticky="ew")

        ctk.CTkLabel(card, text="out_name", text_color=TEXT_DIM).grid(row=1, column=0, padx=12, pady=(0, 6), sticky="w")
        ctk.CTkEntry(card, textvariable=self.surv_out, fg_color=BG_INPUT, border_color="#2D3348").grid(
            row=1, column=1, padx=12, pady=(0, 6), sticky="ew"
        )

        ctk.CTkLabel(card, text="col_time", text_color=TEXT_DIM).grid(row=2, column=0, padx=12, pady=(0, 6), sticky="w")
        ctk.CTkEntry(card, textvariable=self.surv_time, fg_color=BG_INPUT, border_color="#2D3348").grid(
            row=2, column=1, padx=12, pady=(0, 6), sticky="ew"
        )

        ctk.CTkLabel(card, text="col_event", text_color=TEXT_DIM).grid(row=3, column=0, padx=12, pady=(0, 12), sticky="w")
        ctk.CTkEntry(card, textvariable=self.surv_event, fg_color=BG_INPUT, border_color="#2D3348").grid(
            row=3, column=1, padx=12, pady=(0, 12), sticky="ew"
        )

        btn_row = ctk.CTkFrame(left, fg_color=BG_DARK)
        btn_row.pack(fill="x", padx=8, pady=(6, 8))
        btn_row.columnconfigure(0, weight=1)
        btn_row.columnconfigure(1, weight=1)

        self._btn_surv_preview = ctk.CTkButton(
            btn_row,
            text="Preview",
            fg_color=ACCENT2,
            hover_color="#6B4DB0",
            height=40,
            command=self._run_survival_preview,
        )
        self._btn_surv_preview.grid(row=0, column=0, padx=(0, 4), sticky="ew")

        self._btn_surv_save_pdf = ctk.CTkButton(
            btn_row,
            text="Save PDF",
            fg_color=ACCENT,
            hover_color="#00A884",
            height=40,
            state="disabled",
            command=self._save_survival_pdf,
        )
        self._btn_surv_save_pdf.grid(row=0, column=1, padx=(4, 0), sticky="ew")

        _shdr = ctk.CTkFrame(right, fg_color="transparent")
        _shdr.pack(fill="x", padx=8, pady=(10, 2))
        ctk.CTkLabel(
            _shdr, text="Preview",
            font=ctk.CTkFont(size=12, weight="bold"),
            text_color=ACCENT,
        ).pack(side="left")
        ctk.CTkButton(
            _shdr, text="✕ Clear", width=70, height=24,
            fg_color=BG_MID, hover_color="#2D3348", text_color=TEXT_DIM,
            font=ctk.CTkFont(size=11),
            command=self._clear_surv_preview,
        ).pack(side="right")
        self._surv_preview_label = ctk.CTkLabel(
            right,
            text="No preview yet.\nRun survival analysis to see it here.",
            text_color=TEXT_DIM,
            font=ctk.CTkFont(size=11),
        )
        self._surv_preview_label.pack(fill="both", expand=True, padx=8, pady=(4, 12))

    def _run_survival_preview(self):
        wd = self.work_dir_var.get().strip()
        mp = self.model_path_var.get().strip()
        analysis_display = self.surv_analysis.get().strip()
        analysis = self._surv_analysis_map.get(analysis_display, analysis_display)
        out_name = self.surv_out.get().strip() or "survival"
        extra = {
            "col_time": self.surv_time.get().strip(), 
            "col_event": self.surv_event.get().strip(),
            "preview_only": "TRUE",
        }

        ok, err = self._validate_paths()
        if not ok:
            messagebox.showerror("Error", err)
            return
        
        p = {
            "wd": wd,
            "mp": mp,
            "analysis": analysis,
            "out_name": out_name,
            "extra": extra,
        }
        self._last_surv_params = p
        self._btn_surv_preview.configure(state="disabled")
        self._btn_surv_save_pdf.configure(state="disabled")
        self._surv_preview_label.configure(text="Generating preview...", image=None)
        self._surv_preview_image = None

        def run():
            code, out = run_mofa_survival(work_dir=wd, model_path=mp, analysis=analysis, out_name=out_name, extra_params=extra)

            def done():
                self._log(out + "\n")
                if code == 0:
                    self._log("\nPreview OK — click Save PDF to export.\n")
                    self._load_preview(wd, out_name, label=self._surv_preview_label)
                    self._btn_surv_save_pdf.configure(state="normal")
                else:
                    self._log(f"\nERROR code={code}\n")
                    self._surv_preview_label.configure(text="Analysis failed. See log.", image=None)
                self._btn_surv_preview.configure(state="normal")

            self.after(0, done)

        threading.Thread(target=run, daemon=True).start()

    def _save_survival_pdf(self):
        if self._last_surv_params is None:
            return
        p = dict(self._last_surv_params)
        p["out_name"] = self.surv_out.get().strip() or "survival"
        p["extra"] = {k: v for k, v in p["extra"].items() if k != "preview_only"}
        self._btn_surv_save_pdf.configure(state="disabled")
        self._log("\nSaving PDF...\n")

        def run():
            code, out = run_mofa_survival(
                work_dir=p["wd"], 
                model_path=p["mp"], 
                analysis=p["analysis"], 
                out_name=p["out_name"], 
                extra_params=p["extra"],
            )

            def done():
                self._log(out + "\n")
                if code == 0:
                    self._log(f"\n✓ PDF saved: {p['wd']}/{p['out_name']}.pdf\n")
                else:
                    self._log(f"\nERROR code={code}\n")
                self._btn_surv_save_pdf.configure(state="normal")

            self.after(0, done)

        threading.Thread(target=run, daemon=True).start()

    def _build_tab_gsea(self, tab):
        tab.configure(fg_color=BG_DARK)

        self.gsea_analysis = tk.StringVar(value="gsea_heatmap")
        self.gsea_out = tk.StringVar(value="gsea")
        self.gsea_view = tk.StringVar(value="mRNA")
        self.gsea_factors = tk.StringVar(value="1:3")
        self.gsea_sign = tk.StringVar(value="positive")
        self.gsea_test = tk.StringVar(value="parametric")
        self.gsea_source = tk.StringVar(value="reactome")
        self.gsea_result_type = tk.StringVar(value="positive")
        self.gsea_factor_plot = tk.StringVar(value="5")
        self.gsea_max_pathways = tk.StringVar(value="15")

        pane = ctk.CTkFrame(tab, fg_color=BG_DARK)
        pane.pack(fill="both", expand=True)
        pane.columnconfigure(0, weight=1)
        pane.columnconfigure(1, weight=2)
        pane.rowconfigure(0, weight=1)

        left = ctk.CTkFrame(pane, fg_color=BG_DARK)
        left.grid(row=0, column=0, sticky="nsew", padx=(0, 6))
        right = ctk.CTkFrame(pane, fg_color=BG_CARD, corner_radius=10)
        right.grid(row=0, column=1, sticky="nsew")

        # Single GSEA Section
        card = ctk.CTkFrame(left, fg_color=BG_CARD, corner_radius=10)
        card.pack(fill="x", padx=8, pady=8)
        card.columnconfigure(1, weight=1)

        def row(label, var, values=None, r=0):
            ctk.CTkLabel(card, text=label, text_color=TEXT_DIM).grid(row=r, column=0, padx=12, pady=(12 if r == 0 else 0, 6), sticky="w")
            if values:
                ctk.CTkComboBox(
                    card,
                    variable=var,
                    values=values,
                    fg_color=BG_INPUT,
                    border_color="#2D3348",
                    button_color=ACCENT2,
                    dropdown_fg_color=BG_MID,
                ).grid(row=r, column=1, padx=12, pady=(12 if r == 0 else 0, 6), sticky="ew")
            else:
                ctk.CTkEntry(card, textvariable=var, fg_color=BG_INPUT, border_color="#2D3348").grid(
                    row=r, column=1, padx=12, pady=(12 if r == 0 else 0, 6), sticky="ew"
                )

        row("plot_type", self.gsea_analysis, ["gsea_heatmap", "gsea_top_pathways", "gsea_detailed"], 0)
        row("out_name", self.gsea_out, None, 1)
        row("view_gsea", self.gsea_view, None, 2)
        row("factors_gsea", self.gsea_factors, None, 3)
        row("sign_gsea", self.gsea_sign, ["positive", "negative"], 4)
        row("stat_test", self.gsea_test, ["parametric", "cor.adj.parametric", "permutation"], 5)
        row("geneset_source", self.gsea_source, ["reactome", "msigdb_c2_human", "msigdb_c5_human", "msigdb_c2_mouse", "msigdb_c5_mouse", "custom"], 6)
        row("factor_plot", self.gsea_factor_plot, None, 7)
        row("max_pathways", self.gsea_max_pathways, None, 8)

        btn_row = ctk.CTkFrame(left, fg_color=BG_DARK)
        btn_row.pack(fill="x", padx=8, pady=(6, 8))
        btn_row.columnconfigure(0, weight=1)
        btn_row.columnconfigure(1, weight=1)

        self._btn_gsea_preview = ctk.CTkButton(
            btn_row,
            text="Preview",
            fg_color=ACCENT2,
            hover_color="#6B4DB0",
            height=40,
            command=self._run_gsea_preview,
        )
        self._btn_gsea_preview.grid(row=0, column=0, padx=(0, 4), sticky="ew")

        self._btn_gsea_save_pdf = ctk.CTkButton(
            btn_row,
            text="Save PDF",
            fg_color=ACCENT,
            hover_color="#00A884",
            height=40,
            state="disabled",
            command=self._save_gsea_pdf,
        )
        self._btn_gsea_save_pdf.grid(row=0, column=1, padx=(4, 0), sticky="ew")

        _ghdr = ctk.CTkFrame(right, fg_color="transparent")
        _ghdr.pack(fill="x", padx=8, pady=(10, 2))
        ctk.CTkLabel(
            _ghdr, text="Preview",
            font=ctk.CTkFont(size=12, weight="bold"),
            text_color=ACCENT,
        ).pack(side="left")
        ctk.CTkButton(
            _ghdr, text="✕ Clear", width=70, height=24,
            fg_color=BG_MID, hover_color="#2D3348", text_color=TEXT_DIM,
            font=ctk.CTkFont(size=11),
            command=self._clear_gsea_preview,
        ).pack(side="right")
        self._gsea_preview_label = ctk.CTkLabel(
            right,
            text="No preview yet.\nRun GSEA to see it here.",
            text_color=TEXT_DIM,
            font=ctk.CTkFont(size=11),
        )
        self._gsea_preview_label.pack(fill="both", expand=True, padx=8, pady=(4, 12))

    def _run_gsea_preview(self):
        wd = self.work_dir_var.get().strip()
        mp = self.model_path_var.get().strip()
        analysis = self.gsea_analysis.get().strip()
        out_name = self.gsea_out.get().strip() or "gsea"
        sign = self.gsea_sign.get().strip()

        # Preview toujours: run enrichment avec le signe choisi + génère le plot demandé
        extra = {
            "geneset_source": self.gsea_source.get().strip(),
            "custom_gs_path": "",
            "view_gsea": self.gsea_view.get().strip(),
            "factors_gsea": self.gsea_factors.get().strip(),
            "sign_gsea": sign,
            "stat_test": self.gsea_test.get().strip(),
            "stat_test2": "",
            "gene_case": "as_is",
            "load_result": "FALSE",
            "factor_plot": self.gsea_factor_plot.get().strip(),
            "max_pathways": self.gsea_max_pathways.get().strip(),
            "preview_only": "TRUE",
        }

        ok, err = self._validate_paths()
        if not ok:
            messagebox.showerror("Error", err)
            return
        
        p = {
            "wd": wd,
            "mp": mp,
            "analysis": analysis,
            "out_name": out_name,
            "extra": extra,
        }
        self._last_gsea_plot_params = p
        self._btn_gsea_preview.configure(state="disabled")
        self._btn_gsea_save_pdf.configure(state="disabled")
        self._gsea_preview_label.configure(text="Generating preview...", image=None)
        self._gsea_preview_image = None

        def run():
            code, out = run_mofa_gsea(work_dir=wd, model_path=mp, analysis=analysis, out_name=out_name, extra_params=extra)

            def done():
                self._log(out + "\n")
                if code == 0:
                    self._log("\nPreview OK — click Save PDF to export.\n")
                    self._load_preview(wd, out_name, label=self._gsea_preview_label)
                    self._btn_gsea_save_pdf.configure(state="normal")
                else:
                    self._log(f"\nERROR code={code}\n")
                    self._gsea_preview_label.configure(text="Analysis failed. See log.", image=None)
                self._btn_gsea_preview.configure(state="normal")

            self.after(0, done)

        threading.Thread(target=run, daemon=True).start()

    def _save_gsea_pdf(self):
        if self._last_gsea_plot_params is None:
            return
        p = dict(self._last_gsea_plot_params)
        p["out_name"] = self.gsea_out.get().strip() or "gsea_plot"
        p["extra"] = {k: v for k, v in p["extra"].items() if k != "preview_only"}
        self._btn_gsea_save_pdf.configure(state="disabled")
        self._log("\nSaving PDF...\n")

        def run():
            code, out = run_mofa_gsea(
                work_dir=p["wd"], 
                model_path=p["mp"], 
                analysis=p["analysis"], 
                out_name=p["out_name"], 
                extra_params=p["extra"],
            )

            def done():
                self._log(out + "\n")
                if code == 0:
                    self._log(f"\n✓ PDF saved: {p['wd']}/{p['out_name']}.pdf\n")
                else:
                    self._log(f"\nERROR code={code}\n")
                self._btn_gsea_save_pdf.configure(state="normal")

            self.after(0, done)

        threading.Thread(target=run, daemon=True).start()  # Added parentheses here
