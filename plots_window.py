"""Visualization window (plots / survival / GSEA) for AI Sirius MOFA."""

from __future__ import annotations

import os
import threading
from pathlib import Path

import customtkinter as ctk
import tkinter as tk
from tkinter import filedialog, messagebox

from sirius_executor import (
    run_mofa_clustering,
    run_mofa_gsea,
    run_mofa_plots,
    run_mofa_survival,
)


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
        self._clust_preview_image: object = None
        self._last_plot_params: dict | None = None
        self._last_gsea_plot_params: dict | None = None
        # Last produced PDF path per tab — used by the "Open PDF" buttons.
        self._last_clust_pdf: str | None = None
        self._last_surv_pdf: str | None = None
        # Preview labels are created lazily by _build_tab_*; default to None
        # so legacy helpers (e.g. _load_preview) don't raise AttributeError.
        self._clust_preview_label = None
        self._surv_preview_label = None

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
        tab_clust = tabs.add("Clustering")
        tab_surv = tabs.add("Survival")
        tab_gsea = tabs.add("GSEA")

        self._build_tab_plots(tab_plots)
        self._build_tab_clustering(tab_clust)
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

    def _clear_clust_preview(self):
        self._clust_preview_image = None
        if self._clust_preview_label is not None:
            self._safe_configure_label(
                self._clust_preview_label,
                text="No preview yet.\nRun an analysis to see it here.",
                image=None,
            )

    def _open_last_clust_pdf(self):
        if not self._last_clust_pdf or not os.path.isfile(self._last_clust_pdf):
            self._set_cluster_error("No PDF available — run an analysis first.")
            return
        try:
            os.startfile(self._last_clust_pdf)  # type: ignore[attr-defined]
        except Exception as e:
            self._set_cluster_error(f"Cannot open PDF: {e}")

    def _open_last_surv_pdf(self):
        if not self._last_surv_pdf or not os.path.isfile(self._last_surv_pdf):
            self._set_surv_status("No PDF available — run an analysis first.")
            return
        try:
            os.startfile(self._last_surv_pdf)  # type: ignore[attr-defined]
        except Exception as e:
            self._set_surv_status(f"Cannot open PDF: {e}")

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
            elif label is self._clust_preview_label:
                self._clust_preview_image = tk_img
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

    def _build_tab_clustering(self, tab):
        tab.configure(fg_color=BG_DARK)

        self.method_var       = tk.StringVar(value="kmeans")
        self.k_var            = tk.StringVar(value="3")
        self.var_thresh_var   = tk.StringVar(value="0")
        self.factor_x_var     = tk.StringVar(value="1")
        self.factor_y_var     = tk.StringVar(value="2")
        self.covariates_var   = tk.StringVar(value="")

        pane = ctk.CTkFrame(tab, fg_color=BG_DARK)
        pane.pack(fill="both", expand=True)
        pane.columnconfigure(0, weight=2)
        pane.columnconfigure(1, weight=3)
        pane.rowconfigure(0, weight=1)

        left = ctk.CTkFrame(pane, fg_color=BG_DARK)
        left.grid(row=0, column=0, sticky="nsew", padx=(0, 6))
        right = ctk.CTkFrame(pane, fg_color=BG_CARD, corner_radius=10)
        right.grid(row=0, column=1, sticky="nsew")

        # --- Right preview pane ---
        _hdr_cl = ctk.CTkFrame(right, fg_color="transparent")
        _hdr_cl.pack(fill="x", padx=8, pady=(10, 2))
        ctk.CTkLabel(
            _hdr_cl, text="Preview",
            font=ctk.CTkFont(size=12, weight="bold"),
            text_color=ACCENT,
        ).pack(side="left")
        ctk.CTkButton(
            _hdr_cl, text="Open PDF", width=90, height=24,
            fg_color=BG_MID, hover_color="#2D3348", text_color=TEXT_DIM,
            font=ctk.CTkFont(size=11),
            command=self._open_last_clust_pdf,
        ).pack(side="right", padx=(0, 6))
        ctk.CTkButton(
            _hdr_cl, text="✕ Clear", width=70, height=24,
            fg_color=BG_MID, hover_color="#2D3348", text_color=TEXT_DIM,
            font=ctk.CTkFont(size=11),
            command=self._clear_clust_preview,
        ).pack(side="right")
        self._clust_preview_label = ctk.CTkLabel(
            right,
            text="No preview yet.\nRun an analysis to see it here.",
            text_color=TEXT_DIM,
            font=ctk.CTkFont(size=11),
        )
        self._clust_preview_label.pack(fill="both", expand=True, padx=8, pady=(4, 12))

        # --- Section "Clustering settings" ---
        card_cl = ctk.CTkFrame(left, fg_color=BG_CARD, corner_radius=10)
        card_cl.pack(fill="x", padx=8, pady=(8, 4))
        card_cl.columnconfigure(1, weight=1)

        ctk.CTkLabel(
            card_cl, text="Clustering settings",
            font=ctk.CTkFont(size=12, weight="bold"),
            text_color=ACCENT,
        ).grid(row=0, column=0, columnspan=2, padx=12, pady=(12, 8), sticky="w")

        def _row_entry(parent, label, var, r):
            ctk.CTkLabel(parent, text=label, text_color=TEXT_DIM).grid(
                row=r, column=0, padx=12, pady=(0, 6), sticky="w"
            )
            ctk.CTkEntry(
                parent,
                textvariable=var,
                font=ctk.CTkFont(family="Consolas", size=11),
                fg_color=BG_INPUT,
                border_color="#2D3348",
                border_width=1,
            ).grid(row=r, column=1, padx=12, pady=(0, 6), sticky="ew")

        ctk.CTkLabel(card_cl, text="Method", text_color=TEXT_DIM).grid(
            row=1, column=0, padx=12, pady=(0, 6), sticky="w"
        )
        ctk.CTkOptionMenu(
            card_cl,
            variable=self.method_var,
            values=["kmeans", "pam", "hclust_ward"],
            fg_color=BG_INPUT,
            button_color=ACCENT2,
            dropdown_fg_color=BG_MID,
        ).grid(row=1, column=1, padx=12, pady=(0, 6), sticky="ew")

        _row_entry(card_cl, "Number of clusters (k)",        self.k_var,         2)
        _row_entry(card_cl, "Min variance explained (%)",    self.var_thresh_var, 3)
        _row_entry(card_cl, "Factor X axis",                  self.factor_x_var,  4)
        _row_entry(card_cl, "Factor Y axis",                  self.factor_y_var,  5)

        btn_cl = ctk.CTkFrame(left, fg_color=BG_DARK)
        btn_cl.pack(fill="x", padx=8, pady=(6, 8))
        btn_cl.columnconfigure(0, weight=1)
        btn_cl.columnconfigure(1, weight=1)

        self._btn_suggest_k = ctk.CTkButton(
            btn_cl,
            text="Suggest k",
            fg_color=ACCENT2,
            hover_color="#6B4DB0",
            height=40,
            command=self._run_suggest_k,
        )
        self._btn_suggest_k.grid(row=0, column=0, padx=(0, 4), sticky="ew")

        self._btn_run_clustering = ctk.CTkButton(
            btn_cl,
            text="Run Clustering",
            fg_color=ACCENT,
            hover_color="#00A884",
            height=40,
            command=self._run_clustering_action,
        )
        self._btn_run_clustering.grid(row=0, column=1, padx=(4, 0), sticky="ew")

        # --- Section "Survival on clusters" ---
        card_surv = ctk.CTkFrame(left, fg_color=BG_CARD, corner_radius=10)
        card_surv.pack(fill="x", padx=8, pady=(4, 4))
        card_surv.columnconfigure(1, weight=1)

        ctk.CTkLabel(
            card_surv, text="Survival on clusters",
            font=ctk.CTkFont(size=12, weight="bold"),
            text_color=ACCENT,
        ).grid(row=0, column=0, columnspan=2, padx=12, pady=(12, 8), sticky="w")

        _row_entry(
            card_surv,
            "Covariates (comma-separated, e.g. age,stage)",
            self.covariates_var,
            1,
        )

        # Time / event columns reuse the survival tab vars if defined; else default.
        ctk.CTkLabel(card_surv, text="col_time", text_color=TEXT_DIM).grid(
            row=2, column=0, padx=12, pady=(0, 6), sticky="w"
        )
        self.cluster_time_var = tk.StringVar(value="time")
        ctk.CTkEntry(
            card_surv, textvariable=self.cluster_time_var,
            fg_color=BG_INPUT, border_color="#2D3348",
        ).grid(row=2, column=1, padx=12, pady=(0, 6), sticky="ew")

        ctk.CTkLabel(card_surv, text="col_event", text_color=TEXT_DIM).grid(
            row=3, column=0, padx=12, pady=(0, 12), sticky="w"
        )
        self.cluster_event_var = tk.StringVar(value="event")
        ctk.CTkEntry(
            card_surv, textvariable=self.cluster_event_var,
            fg_color=BG_INPUT, border_color="#2D3348",
        ).grid(row=3, column=1, padx=12, pady=(0, 12), sticky="ew")

        btn_sv = ctk.CTkFrame(left, fg_color=BG_DARK)
        btn_sv.pack(fill="x", padx=8, pady=(6, 8))
        btn_sv.columnconfigure(0, weight=1)
        btn_sv.columnconfigure(1, weight=1)
        btn_sv.columnconfigure(2, weight=1)

        self._btn_km_clusters = ctk.CTkButton(
            btn_sv,
            text="Kaplan-Meier by cluster",
            fg_color=ACCENT2,
            hover_color="#6B4DB0",
            height=40,
            command=self._run_km_clusters,
        )
        self._btn_km_clusters.grid(row=0, column=0, padx=(0, 4), sticky="ew")

        self._btn_cox_clusters = ctk.CTkButton(
            btn_sv,
            text="Cox by cluster",
            fg_color=ACCENT,
            hover_color="#00A884",
            height=40,
            command=self._run_cox_clusters,
        )
        self._btn_cox_clusters.grid(row=0, column=1, padx=4, sticky="ew")

        self._btn_cox_diag = ctk.CTkButton(
            btn_sv,
            text="Cox diagnostics",
            fg_color=ACCENT2,
            hover_color="#6B4DB0",
            height=40,
            command=self._run_cox_diagnostics,
        )
        self._btn_cox_diag.grid(row=0, column=2, padx=(4, 0), sticky="ew")

        # Error / status label
        self._cluster_error_label = ctk.CTkLabel(
            left,
            text="",
            text_color="#FF6B6B",
            font=ctk.CTkFont(size=11),
            anchor="w",
            justify="left",
        )
        self._cluster_error_label.pack(fill="x", padx=12, pady=(4, 8))

        self._cluster_buttons = [
            self._btn_suggest_k,
            self._btn_run_clustering,
            self._btn_km_clusters,
            self._btn_cox_clusters,
            self._btn_cox_diag,
        ]

    # -------- Clustering actions --------
    def _set_cluster_buttons_state(self, state: str):
        for b in self._cluster_buttons:
            try:
                b.configure(state=state)
            except Exception:
                pass

    def _set_cluster_error(self, text: str, *, error: bool = True):
        try:
            self._cluster_error_label.configure(
                text=text,
                text_color="#FF6B6B" if error else TEXT_DIM,
            )
        except Exception:
            pass

    def _open_pdf(self, pdf_path: str):
        try:
            if os.path.isfile(pdf_path):
                os.startfile(pdf_path)  # type: ignore[attr-defined]
            else:
                self._set_cluster_error(f"File not found: {pdf_path}")
        except Exception as e:
            self._set_cluster_error(f"Cannot open PDF: {e}")

    def _parse_int(self, s: str, default: int) -> int:
        try:
            return int(str(s).strip())
        except Exception:
            return default

    def _parse_float(self, s: str, default: float) -> float:
        try:
            return float(str(s).strip())
        except Exception:
            return default

    def _run_clustering_in_thread(self, *, analysis: str, busy_text: str,
                                  output_pdf: str, use_survival: bool = False,
                                  surv_analysis: str | None = None):
        ok, err = self._validate_paths()
        if not ok:
            messagebox.showerror("Error", err)
            return

        wd = self.work_dir_var.get().strip()
        mp = self.model_path_var.get().strip()
        method = self.method_var.get().strip() or "kmeans"
        k = self._parse_int(self.k_var.get(), 3)
        var_t = self._parse_float(self.var_thresh_var.get(), 0.0)
        fx = self._parse_int(self.factor_x_var.get(), 1)
        fy = self._parse_int(self.factor_y_var.get(), 2)
        covariates = self.covariates_var.get().strip() or None

        self._set_cluster_buttons_state("disabled")
        self._set_cluster_error(busy_text, error=False)

        def run():
            if use_survival:
                extra = {
                    "col_time": self.cluster_time_var.get().strip() or "time",
                    "col_event": self.cluster_event_var.get().strip() or "event",
                }
                if covariates:
                    extra["covariates"] = covariates
                code, out = run_mofa_survival(
                    work_dir=wd,
                    model_path=mp,
                    analysis=surv_analysis or analysis,
                    out_name="cluster_surv",
                    extra_params=extra,
                )
            else:
                code, out = run_mofa_clustering(
                    model_path=mp,
                    output_dir=wd,
                    analysis=analysis,
                    method=method,
                    k=k,
                    var_threshold=var_t,
                    factor_x=fx,
                    factor_y=fy,
                    covariates=covariates,
                )

            def done():
                self._log(out + "\n")
                if code == 0:
                    self._log("\nOK\n")
                    self._set_cluster_error("", error=False)
                    pdf_path = str(Path(wd) / output_pdf)
                    self._last_clust_pdf = pdf_path
                    # Load inline preview (PNG produced alongside the PDF)
                    basename = output_pdf[:-4] if output_pdf.lower().endswith(".pdf") else output_pdf
                    if self._clust_preview_label is not None:
                        self._load_preview(wd, basename, label=self._clust_preview_label)
                    # Refresh cluster-dependent survival buttons (clusters.rds may now exist)
                    if hasattr(self, "_refresh_surv_cluster_buttons"):
                        try:
                            self._refresh_surv_cluster_buttons()
                        except Exception:
                            pass
                else:
                    self._log(f"\nERROR code={code}\n")
                    self._set_cluster_error(
                        f"Failed (code={code}). See log for details."
                    )
                self._set_cluster_buttons_state("normal")

            self.after(0, done)

        threading.Thread(target=run, daemon=True).start()

    def _run_suggest_k(self):
        self._run_clustering_in_thread(
            analysis="choose_k",
            busy_text="Computing WSS / Silhouette / Gap...",
            output_pdf="choose_k.pdf",
        )

    def _run_clustering_action(self):
        self._run_clustering_in_thread(
            analysis="run_clustering",
            busy_text="Running clustering...",
            output_pdf="cluster_scatter.pdf",
        )

    def _run_km_clusters(self):
        self._run_clustering_in_thread(
            analysis="kaplan_meier_clusters",
            busy_text="Kaplan-Meier on clusters...",
            output_pdf="km_clusters.pdf",
            use_survival=True,
            surv_analysis="kaplan_meier_clusters",
        )

    def _run_cox_clusters(self):
        self._run_clustering_in_thread(
            analysis="cox_clusters",
            busy_text="Cox regression on clusters...",
            output_pdf="cox_clusters.pdf",
            use_survival=True,
            surv_analysis="cox_clusters",
        )

    def _run_cox_diagnostics(self):
        self._run_clustering_in_thread(
            analysis="cox_diagnostics",
            busy_text="Cox PH diagnostics...",
            output_pdf="cox_zph.pdf",
            use_survival=True,
            surv_analysis="cox_diagnostics",
        )

    @staticmethod
    def check_clusters_exist(work_dir: str) -> bool:
        """Return True if clusters.rds is present in `work_dir`."""
        try:
            return bool(work_dir) and os.path.exists(os.path.join(work_dir, "clusters.rds"))
        except Exception:
            return False

    def _build_tab_survival(self, tab):
        tab.configure(fg_color=BG_DARK)

        # ── Variable definitions (all StringVar / BooleanVar local to this tab)
        # NB: model_path_var and work_dir_var are shared with the global header.
        self.surv_out_name_var     = tk.StringVar(value="survival_output")
        self.surv_col_time_var     = tk.StringVar(value="time")
        self.surv_col_event_var    = tk.StringVar(value="event")
        self.surv_factor_index_var = tk.StringVar(value="1")
        self.surv_km_label_low_var = tk.StringVar(value="")
        self.surv_km_label_high_var = tk.StringVar(value="")
        self.surv_title_var        = tk.StringVar(value="")
        self.surv_xlab_var         = tk.StringVar(value="Time")
        self.surv_conf_int_var     = tk.BooleanVar(value=True)
        self.surv_pval_show_var    = tk.BooleanVar(value=True)
        self.surv_covariates_var   = tk.StringVar(value="")

        self._surv_buttons: list[ctk.CTkButton] = []
        self._surv_cluster_buttons: list[ctk.CTkButton] = []

        pane = ctk.CTkFrame(tab, fg_color=BG_DARK)
        pane.pack(fill="both", expand=True)
        pane.columnconfigure(0, weight=2)
        pane.columnconfigure(1, weight=3)
        pane.rowconfigure(0, weight=1)

        left = ctk.CTkScrollableFrame(pane, fg_color=BG_DARK)
        left.grid(row=0, column=0, sticky="nsew", padx=(0, 6))
        right = ctk.CTkFrame(pane, fg_color=BG_CARD, corner_radius=10)
        right.grid(row=0, column=1, sticky="nsew")

        # --- Right preview pane ---
        _hdr_sv = ctk.CTkFrame(right, fg_color="transparent")
        _hdr_sv.pack(fill="x", padx=8, pady=(10, 2))
        ctk.CTkLabel(
            _hdr_sv, text="Preview",
            font=ctk.CTkFont(size=12, weight="bold"),
            text_color=ACCENT,
        ).pack(side="left")
        ctk.CTkButton(
            _hdr_sv, text="Open PDF", width=90, height=24,
            fg_color=BG_MID, hover_color="#2D3348", text_color=TEXT_DIM,
            font=ctk.CTkFont(size=11),
            command=self._open_last_surv_pdf,
        ).pack(side="right", padx=(0, 6))
        ctk.CTkButton(
            _hdr_sv, text="✕ Clear", width=70, height=24,
            fg_color=BG_MID, hover_color="#2D3348", text_color=TEXT_DIM,
            font=ctk.CTkFont(size=11),
            command=self._clear_surv_preview,
        ).pack(side="right")
        self._surv_preview_label = ctk.CTkLabel(
            right,
            text="No preview yet.\nRun an analysis to see it here.",
            text_color=TEXT_DIM,
            font=ctk.CTkFont(size=11),
        )
        self._surv_preview_label.pack(fill="both", expand=True, padx=8, pady=(4, 12))

        # =========== Card 1 — Required parameters ===========
        card_req = ctk.CTkFrame(left, fg_color=BG_CARD, corner_radius=10)
        card_req.pack(fill="x", padx=8, pady=(8, 4))
        card_req.columnconfigure(1, weight=1)

        ctk.CTkLabel(
            card_req, text="Required parameters",
            font=ctk.CTkFont(size=12, weight="bold"),
            text_color=ACCENT,
        ).grid(row=0, column=0, columnspan=3, padx=12, pady=(12, 8), sticky="w")

        def _entry(parent, var, r, col=1, columnspan=1):
            e = ctk.CTkEntry(
                parent,
                textvariable=var,
                font=ctk.CTkFont(family="Consolas", size=11),
                fg_color=BG_INPUT,
                border_color="#2D3348",
                border_width=1,
            )
            e.grid(row=r, column=col, columnspan=columnspan,
                   padx=12, pady=(0, 6), sticky="ew")
            return e

        def _label(parent, text, r):
            ctk.CTkLabel(parent, text=text, text_color=TEXT_DIM).grid(
                row=r, column=0, padx=12, pady=(0, 6), sticky="w"
            )

        # Model path with Browse
        _label(card_req, "Model path (.hdf5 / .rds)", 1)
        ctk.CTkEntry(
            card_req,
            textvariable=self.model_path_var,
            font=ctk.CTkFont(family="Consolas", size=11),
            fg_color=BG_INPUT, border_color="#2D3348", border_width=1,
        ).grid(row=1, column=1, padx=(12, 4), pady=(0, 6), sticky="ew")
        ctk.CTkButton(
            card_req, text="Browse", width=80, height=28,
            fg_color=BG_MID, hover_color="#2D3348",
            command=self._browse_model_path,
        ).grid(row=1, column=2, padx=(0, 12), pady=(0, 6), sticky="e")

        # Work dir with Browse
        _label(card_req, "Output directory", 2)
        ctk.CTkEntry(
            card_req,
            textvariable=self.work_dir_var,
            font=ctk.CTkFont(family="Consolas", size=11),
            fg_color=BG_INPUT, border_color="#2D3348", border_width=1,
        ).grid(row=2, column=1, padx=(12, 4), pady=(0, 6), sticky="ew")
        ctk.CTkButton(
            card_req, text="Browse", width=80, height=28,
            fg_color=BG_MID, hover_color="#2D3348",
            command=self._browse_work_dir,
        ).grid(row=2, column=2, padx=(0, 12), pady=(0, 6), sticky="e")

        _label(card_req, "Output name", 3)
        _entry(card_req, self.surv_out_name_var, 3, columnspan=2)

        _label(card_req, "Time column (in metadata)", 4)
        _entry(card_req, self.surv_col_time_var, 4, columnspan=2)

        _label(card_req, "Event column (in metadata, 0/1)", 5)
        _entry(card_req, self.surv_col_event_var, 5, columnspan=2)

        # Spacer
        ctk.CTkLabel(card_req, text="", height=4).grid(row=6, column=0)

        # =========== Card 2 — Analysis buttons ===========
        card_btn = ctk.CTkFrame(left, fg_color=BG_CARD, corner_radius=10)
        card_btn.pack(fill="x", padx=8, pady=(4, 4))
        card_btn.columnconfigure(0, weight=1)

        ctk.CTkLabel(
            card_btn, text="Analyses",
            font=ctk.CTkFont(size=12, weight="bold"),
            text_color=ACCENT,
        ).grid(row=0, column=0, padx=12, pady=(12, 8), sticky="w")

        def _make_btn(parent, text, analysis, color, r, *, cluster=False):
            b = ctk.CTkButton(
                parent,
                text=text,
                fg_color=color,
                hover_color="#6B4DB0" if color == ACCENT2 else "#00A884",
                height=38,
                command=lambda a=analysis: self._run_survival_action(a),
            )
            b.grid(row=r, column=0, padx=12, pady=4, sticky="ew")
            self._surv_buttons.append(b)
            if cluster:
                b.configure(state="disabled")
                self._surv_cluster_buttons.append(b)
            return b

        self._btn_surv_cox_hr        = _make_btn(card_btn, "Cox HR — Forest plot",         "cox_hr",                ACCENT2, 1)
        self._btn_surv_cox_summary   = _make_btn(card_btn, "Cox Summary — Bubble plot",    "cox_summary",           ACCENT2, 2)
        self._btn_surv_km            = _make_btn(card_btn, "Kaplan-Meier (1 factor)",      "kaplan_meier",          ACCENT,  3)
        self._btn_surv_km_all        = _make_btn(card_btn, "Kaplan-Meier (all factors)",   "kaplan_meier_all",      ACCENT,  4)
        self._btn_surv_km_clust      = _make_btn(card_btn, "Kaplan-Meier by cluster",      "kaplan_meier_clusters", ACCENT2, 5, cluster=True)
        self._btn_surv_cox_clust     = _make_btn(card_btn, "Cox by cluster",               "cox_clusters",          ACCENT,  6, cluster=True)
        self._btn_surv_cox_diag      = _make_btn(card_btn, "Cox diagnostics (Schoenfeld)", "cox_diagnostics",       ACCENT2, 7, cluster=True)

        # Bottom spacer
        ctk.CTkLabel(card_btn, text="", height=4).grid(row=8, column=0)

        # =========== Card 3 — Kaplan-Meier options ===========
        card_km = ctk.CTkFrame(left, fg_color=BG_CARD, corner_radius=10)
        card_km.pack(fill="x", padx=8, pady=(4, 4))
        card_km.columnconfigure(1, weight=1)

        ctk.CTkLabel(
            card_km, text="Kaplan-Meier options",
            font=ctk.CTkFont(size=12, weight="bold"),
            text_color=ACCENT,
        ).grid(row=0, column=0, columnspan=2, padx=12, pady=(12, 8), sticky="w")

        _label(card_km, "Factor index",        1); _entry(card_km, self.surv_factor_index_var, 1)
        _label(card_km, "Label LOW group",     2); _entry(card_km, self.surv_km_label_low_var,  2)
        _label(card_km, "Label HIGH group",    3); _entry(card_km, self.surv_km_label_high_var, 3)
        _label(card_km, "Plot title",          4); _entry(card_km, self.surv_title_var,         4)
        _label(card_km, "X axis label",        5); _entry(card_km, self.surv_xlab_var,          5)

        ctk.CTkCheckBox(
            card_km, text="Show confidence interval",
            variable=self.surv_conf_int_var,
            fg_color=ACCENT, hover_color="#00A884",
            text_color=TEXT_DIM,
        ).grid(row=6, column=0, columnspan=2, padx=12, pady=(6, 4), sticky="w")

        ctk.CTkCheckBox(
            card_km, text="Show p-value",
            variable=self.surv_pval_show_var,
            fg_color=ACCENT, hover_color="#00A884",
            text_color=TEXT_DIM,
        ).grid(row=7, column=0, columnspan=2, padx=12, pady=(0, 12), sticky="w")

        # =========== Card 4 — Cluster survival options ===========
        card_cl = ctk.CTkFrame(left, fg_color=BG_CARD, corner_radius=10)
        card_cl.pack(fill="x", padx=8, pady=(4, 4))
        card_cl.columnconfigure(1, weight=1)

        ctk.CTkLabel(
            card_cl, text="Cluster survival options",
            font=ctk.CTkFont(size=12, weight="bold"),
            text_color=ACCENT,
        ).grid(row=0, column=0, columnspan=2, padx=12, pady=(12, 8), sticky="w")

        _label(card_cl, "Covariates (comma-separated, e.g. age,stage)", 1)
        _entry(card_cl, self.surv_covariates_var, 1)

        ctk.CTkLabel(
            card_cl,
            text="⚠ Cluster analyses require clustering to be run first",
            text_color="#E0B43E",
            font=ctk.CTkFont(size=11),
        ).grid(row=2, column=0, columnspan=2, padx=12, pady=(2, 12), sticky="w")

        # =========== Status / error label ===========
        self._surv_status_label = ctk.CTkLabel(
            left, text="",
            text_color="#FF6B6B",
            font=ctk.CTkFont(size=11),
            anchor="w", justify="left",
        )
        self._surv_status_label.pack(fill="x", padx=12, pady=(4, 8))

        # ── Trace work_dir for cluster button availability
        try:
            self.work_dir_var.trace_add("write", lambda *a: self._refresh_surv_cluster_buttons())
        except Exception:
            # Fallback for older tkinter
            self.work_dir_var.trace("w", lambda *a: self._refresh_surv_cluster_buttons())
        self._refresh_surv_cluster_buttons()

    # -------- Survival helpers --------
    def _browse_model_path(self):
        path = filedialog.askopenfilename(
            title="Select MOFA model",
            filetypes=[
                ("MOFA models", "*.hdf5 *.rds"),
                ("HDF5", "*.hdf5"),
                ("RDS", "*.rds"),
                ("All files", "*.*"),
            ],
        )
        if path:
            self.model_path_var.set(path)

    def _browse_work_dir(self):
        path = filedialog.askdirectory(title="Select output directory")
        if path:
            self.work_dir_var.set(path)

    def _refresh_surv_cluster_buttons(self):
        wd = self.work_dir_var.get().strip()
        enabled = self.check_clusters_exist(wd)
        state = "normal" if enabled else "disabled"
        for b in getattr(self, "_surv_cluster_buttons", []):
            try:
                b.configure(state=state)
            except Exception:
                pass

    def _set_surv_buttons_state(self, state: str):
        for b in getattr(self, "_surv_buttons", []):
            try:
                b.configure(state=state)
            except Exception:
                pass
        # After a global re-enable, re-apply cluster availability gating
        if state == "normal":
            self._refresh_surv_cluster_buttons()

    def _set_surv_status(self, text: str, *, error: bool = True):
        try:
            self._surv_status_label.configure(
                text=text,
                text_color="#FF6B6B" if error else TEXT_DIM,
            )
        except Exception:
            pass

    def _build_survival_extra(self, analysis: str) -> dict[str, str]:
        """Build the extra_params dict mirroring the spec's --key value list."""
        extra: dict[str, str] = {}

        if analysis == "kaplan_meier":
            extra["factor_index"] = self.surv_factor_index_var.get().strip() or "1"
            extra["conf_int"]  = "TRUE" if bool(self.surv_conf_int_var.get())  else "FALSE"
            extra["pval_show"] = "TRUE" if bool(self.surv_pval_show_var.get()) else "FALSE"
            low  = self.surv_km_label_low_var.get().strip()
            high = self.surv_km_label_high_var.get().strip()
            ttl  = self.surv_title_var.get().strip()
            xlab = self.surv_xlab_var.get().strip()
            if low:  extra["km_label_low"]  = low
            if high: extra["km_label_high"] = high
            if ttl:  extra["title_surv"]    = ttl
            if xlab: extra["xlab_surv"]     = xlab

        elif analysis == "kaplan_meier_all":
            extra["conf_int"]  = "TRUE" if bool(self.surv_conf_int_var.get())  else "FALSE"
            extra["pval_show"] = "TRUE" if bool(self.surv_pval_show_var.get()) else "FALSE"

        elif analysis in ("cox_clusters", "cox_diagnostics"):
            cov = self.surv_covariates_var.get().strip()
            if cov:
                extra["covariates"] = cov

        return extra

    # Mapping analysis → output PDF basename (relative to work_dir)
    _SURV_OUTPUT_PDF: dict[str, str] = {
        "cox_hr":                None,            # uses out_name + ".pdf"
        "cox_summary":           None,
        "kaplan_meier":          None,
        "kaplan_meier_all":      None,
        "kaplan_meier_clusters": "km_clusters.pdf",
        "cox_clusters":          "cox_clusters.pdf",
        "cox_diagnostics":       "cox_zph.pdf",
    }

    def _surv_output_pdf_path(self, analysis: str, work_dir: str, out_name: str) -> str:
        fixed = self._SURV_OUTPUT_PDF.get(analysis)
        return str(Path(work_dir) / (fixed if fixed else f"{out_name}.pdf"))

    def _run_survival_action(self, analysis: str):
        # Re-check cluster availability for cluster-based analyses
        if analysis in ("kaplan_meier_clusters", "cox_clusters", "cox_diagnostics"):
            self._refresh_surv_cluster_buttons()
            wd_check = self.work_dir_var.get().strip()
            if not self.check_clusters_exist(wd_check):
                messagebox.showerror(
                    "Error",
                    "clusters.rds not found in the output directory.\n"
                    "Run clustering first (Clustering tab → Run Clustering).",
                )
                return

        ok, err = self._validate_paths()
        if not ok:
            messagebox.showerror("Error", err)
            return

        wd       = self.work_dir_var.get().strip()
        mp       = self.model_path_var.get().strip()
        out_name = self.surv_out_name_var.get().strip() or "survival_output"
        col_time  = self.surv_col_time_var.get().strip()  or "time"
        col_event = self.surv_col_event_var.get().strip() or "event"

        if not col_time or not col_event:
            messagebox.showerror("Error", "Time and Event column names are required.")
            return

        extra = {
            "col_time":  col_time,
            "col_event": col_event,
            **self._build_survival_extra(analysis),
        }

        self._set_surv_buttons_state("disabled")
        self._set_surv_status(f"Running '{analysis}'...", error=False)
        self._log(f"\n[survival] Running analysis: {analysis}\n")

        def run():
            code, out = run_mofa_survival(
                work_dir=wd,
                model_path=mp,
                analysis=analysis,
                out_name=out_name,
                extra_params=extra,
            )

            def done():
                self._log(out + "\n")
                if code == 0:
                    self._log("\nOK\n")
                    self._set_surv_status("", error=False)
                    pdf_path = self._surv_output_pdf_path(analysis, wd, out_name)
                    self._last_surv_pdf = pdf_path
                    # Load inline preview (PNG produced alongside the PDF)
                    pdf_name = os.path.basename(pdf_path)
                    basename = pdf_name[:-4] if pdf_name.lower().endswith(".pdf") else pdf_name
                    if self._surv_preview_label is not None:
                        self._load_preview(wd, basename, label=self._surv_preview_label)
                    if not os.path.isfile(pdf_path):
                        self._set_surv_status(f"PDF not produced: {pdf_path}")
                else:
                    self._log(f"\nERROR code={code}\n")
                    self._set_surv_status(
                        f"'{analysis}' failed (code={code}). See log."
                    )
                self._set_surv_buttons_state("normal")

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
