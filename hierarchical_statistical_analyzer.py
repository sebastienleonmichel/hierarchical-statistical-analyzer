
"""
Hierarchical Statistical Analyzer
=================================

GUI-based framework for:
- stratified nonparametric testing
- linear mixed-effects modeling
- replicate-aware statistical analysis
- skewed biological datasets

Supports:
- blocked experimental designs
- hierarchical inference
- multiple testing correction
- contextual data mapping

Author: Sébastien Terreau
Year: 2026
Version: 2.1.1
"""


import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from pathlib import Path

import numpy as np
import pandas as pd

from scipy.stats import rankdata

import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests


APP_VERSION = "v2.1.1"


# ============================================================
# RESHAPE
# ============================================================

def reshape_long(df_wide, col_to_group, col_to_repl):

    rows = []

    for col in df_wide.columns:

        grp = col_to_group[col]
        rep = col_to_repl[col]

        vals = pd.to_numeric(
            df_wide[col],
            errors="coerce"
        ).dropna()

        vals = vals[vals > 0]

        for v in vals:

            rows.append({

                "value": float(v),
                "group": grp,
                "replicate_id": rep,
                "column": col

            })

    return pd.DataFrame(rows)


# ============================================================
# WILCOXON UTILITIES
# ============================================================

def mannwhitney_u(x, gA_mask):

    r = rankdata(x, method="average")

    n1 = int(gA_mask.sum())

    U = r[gA_mask].sum() - n1 * (n1 + 1) / 2

    return float(U)


def van_elteren_test(
    df_long,
    groupA,
    groupB,
    n_perm=10000,
    seed=0
):

    rng = np.random.default_rng(seed)

    observed = 0.0

    for _, sub in df_long.groupby("replicate_id"):

        sub = sub[
            sub["group"].isin([groupA, groupB])
        ]

        if len(sub) == 0:
            continue

        x = sub["value"].to_numpy()

        gA = (
            sub["group"].to_numpy() == groupA
        )

        if gA.sum() == 0:
            continue

        observed += mannwhitney_u(x, gA)

    perms = []

    for _ in range(n_perm):

        total = 0.0

        for _, sub in df_long.groupby("replicate_id"):

            sub = sub[
                sub["group"].isin([groupA, groupB])
            ]

            if len(sub) == 0:
                continue

            x = sub["value"].to_numpy()

            nA = np.sum(
                sub["group"].to_numpy() == groupA
            )

            idx = np.arange(len(x))

            rng.shuffle(idx)

            mask = np.zeros(
                len(x),
                dtype=bool
            )

            mask[idx[:nA]] = True

            total += mannwhitney_u(x, mask)

        perms.append(total)

    perms = np.array(perms)

    center = perms.mean()

    p = (
        np.sum(
            np.abs(perms - center)
            >= np.abs(observed - center)
        ) + 1
    ) / (n_perm + 1)

    return observed, p


# ============================================================
# LMM
# ============================================================

def run_lmm(
    df_long,
    contrasts,
    log_transform=False
):

    df = df_long.copy()

    if log_transform:

        df["response"] = np.log(df["value"])

    else:

        df["response"] = df["value"]

    rows = []

    raw_p = []

    for A, B in contrasts:

        sub = df[
            df["group"].isin([A, B])
        ].copy()

        # --------------------------------------------
        # REQUIRE MULTIPLE REPLICATES
        # --------------------------------------------

        if sub["replicate_id"].nunique() < 2:

            rows.append({

                "Contrast":
                    f"{A} vs {B}",

                "Method":
                    "LMM failed: <2 replicate levels"

            })

            raw_p.append(np.nan)

            continue

        sub["group"] = (
            sub["group"]
            .astype("category")
        )

        try:

            model = smf.mixedlm(
                "response ~ group",
                data=sub,
                groups=sub["replicate_id"]
            )

            fit = model.fit(reml=False)

            # ----------------------------------------
            # ROBUST EXTRACTION
            # ----------------------------------------

            beta = float(
                fit.params.iloc[1]
            )

            pval = float(
                fit.pvalues.iloc[1]
            )

            raw_p.append(pval)

            rows.append({

                "Contrast":
                    f"{A} vs {B}",

                "Method":
                    "Linear Mixed Model",


                "p_raw":
                    pval

            })

        except Exception as e:

            rows.append({

                "Contrast":
                    f"{A} vs {B}",

                "Method":
                    f"LMM failed: {e}"

            })

            raw_p.append(np.nan)

    # ========================================================
    # MULTIPLE TESTING CORRECTIONS
    # ========================================================

    raw_p_array = np.array(raw_p)

    valid = np.isfinite(raw_p_array)

    holm_adj = np.full(len(raw_p), np.nan)
    bonf_adj = np.full(len(raw_p), np.nan)
    fdr_adj = np.full(len(raw_p), np.nan)

    if np.sum(valid) > 0:

        _, p_holm, _, _ = multipletests(
            raw_p_array[valid],
            method="holm"
        )

        _, p_bonf, _, _ = multipletests(
            raw_p_array[valid],
            method="bonferroni"
        )

        _, p_fdr, _, _ = multipletests(
            raw_p_array[valid],
            method="fdr_bh"
        )

        holm_adj[valid] = p_holm
        bonf_adj[valid] = p_bonf
        fdr_adj[valid] = p_fdr

    for i in range(len(rows)):

        rows[i]["p_bonferroni"] = bonf_adj[i]
        rows[i]["p_holm"] = holm_adj[i]
        rows[i]["p_fdr_bh"] = fdr_adj[i]

    return pd.DataFrame(rows)


# ============================================================
# SCROLLABLE GUI FRAME
# ============================================================

class ScrollableFrame(ttk.Frame):
    """
    A reusable scrollable frame for GUI tabs containing many rows.

    This is used for the Mapping and Contrasts tabs so that datasets
    with many conditions or many pairwise contrasts remain accessible.
    """

    def __init__(self, parent):

        super().__init__(parent)

        self.canvas = tk.Canvas(
            self,
            borderwidth=0,
            highlightthickness=0
        )

        self.v_scrollbar = ttk.Scrollbar(
            self,
            orient="vertical",
            command=self.canvas.yview
        )

        self.scrollable_frame = ttk.Frame(
            self.canvas
        )

        self.window_id = self.canvas.create_window(
            (0, 0),
            window=self.scrollable_frame,
            anchor="nw"
        )

        self.canvas.configure(
            yscrollcommand=self.v_scrollbar.set
        )

        self.canvas.grid(
            row=0,
            column=0,
            sticky="nsew"
        )

        self.v_scrollbar.grid(
            row=0,
            column=1,
            sticky="ns"
        )

        self.rowconfigure(
            0,
            weight=1
        )

        self.columnconfigure(
            0,
            weight=1
        )

        self.scrollable_frame.bind(
            "<Configure>",
            self._update_scroll_region
        )

        self.canvas.bind(
            "<Configure>",
            self._resize_inner_frame
        )

        self.canvas.bind(
            "<Enter>",
            self._bind_mousewheel
        )

        self.canvas.bind(
            "<Leave>",
            self._unbind_mousewheel
        )

    def _update_scroll_region(self, event=None):

        self.canvas.configure(
            scrollregion=self.canvas.bbox("all")
        )

    def _resize_inner_frame(self, event):

        self.canvas.itemconfigure(
            self.window_id,
            width=event.width
        )

    def _bind_mousewheel(self, event=None):

        self.canvas.bind_all(
            "<MouseWheel>",
            self._on_mousewheel
        )

        self.canvas.bind_all(
            "<Button-4>",
            self._on_mousewheel
        )

        self.canvas.bind_all(
            "<Button-5>",
            self._on_mousewheel
        )

    def _unbind_mousewheel(self, event=None):

        self.canvas.unbind_all(
            "<MouseWheel>"
        )

        self.canvas.unbind_all(
            "<Button-4>"
        )

        self.canvas.unbind_all(
            "<Button-5>"
        )

    def _on_mousewheel(self, event):

        if getattr(event, "num", None) == 4:

            delta = -1

        elif getattr(event, "num", None) == 5:

            delta = 1

        else:

            delta = -int(event.delta / 120)

            if delta == 0:

                delta = -1 if event.delta > 0 else 1

        self.canvas.yview_scroll(
            delta,
            "units"
        )


# ============================================================
# GUI
# ============================================================

class HierarchicalStatisticalAnalyzerGUI(tk.Tk):

    def __init__(self):

        super().__init__()

        self.title(
            f"Hierarchical Statistical Analyzer {APP_VERSION}"
        )

        self.geometry("1300x900")

        self.df = None
        self.input_path = None

        self.group_vars = {}
        self.repl_vars = {}
        self.contrast_vars = {}

        self.engine_var = tk.StringVar(
            value="Both"
        )

        self.log_var = tk.BooleanVar(
            value=False
        )

        self._build_gui()

    def _build_gui(self):

        nb = ttk.Notebook(self)

        nb.pack(fill=tk.BOTH, expand=True)

        self.tab_input = ttk.Frame(nb)
        self.tab_mapping = ttk.Frame(nb)
        self.tab_contrast = ttk.Frame(nb)
        self.tab_options = ttk.Frame(nb)
        self.tab_run = ttk.Frame(nb)

        nb.add(self.tab_input, text="1) Input")
        nb.add(self.tab_mapping, text="2) Mapping")
        nb.add(self.tab_contrast, text="3) Contrasts")
        nb.add(self.tab_options, text="4) Options")
        nb.add(self.tab_run, text="5) Run")

        self.build_input()
        self.build_mapping()
        self.build_contrasts()
        self.build_options()
        self.build_run()

    def build_input(self):

        ttk.Button(
            self.tab_input,
            text="Open data file (.csv or .xlsx)",
            command=self.open_data_file
        ).pack(padx=20, pady=(20, 5))

        ttk.Label(
            self.tab_input,
            text=f"Supported input formats: .csv and .xlsx | {APP_VERSION}"
        ).pack(padx=20, pady=(0, 20))

        self.columns_box = tk.Listbox(
            self.tab_input,
            width=100,
            height=25
        )

        self.columns_box.pack(
            padx=20,
            pady=20
        )

    def open_data_file(self):

        path = filedialog.askopenfilename(
            title="Open CSV or XLSX file",
            filetypes=[
                ("CSV or Excel files", "*.csv *.xlsx"),
                ("CSV files", "*.csv"),
                ("Excel files", "*.xlsx"),
                ("All files", "*.*")
            ]
        )

        if not path:
            return

        file_path = Path(path)
        suffix = file_path.suffix.lower()

        try:

            if suffix == ".csv":

                self.df = pd.read_csv(file_path)

            elif suffix == ".xlsx":

                self.df = pd.read_excel(
                    file_path,
                    engine="openpyxl"
                )

            else:

                messagebox.showerror(
                    "Unsupported file",
                    "Please select a .csv or .xlsx file."
                )

                return

        except ImportError:

            messagebox.showerror(
                "Missing dependency",
                "Reading .xlsx files requires openpyxl. Install it with: pip install openpyxl"
            )

            return

        except Exception as e:

            messagebox.showerror(
                "File loading error",
                f"Could not load the selected file:\n{e}"
            )

            return

        self.input_path = file_path

        self.columns_box.delete(0, tk.END)

        for c in self.df.columns:

            self.columns_box.insert(tk.END, c)

        self.refresh_mapping()
        self.refresh_contrasts()

    def open_csv(self):
        # Backwards-compatible alias for older callback names.
        self.open_data_file()

    def build_mapping(self):

        self.mapping_scroll = ScrollableFrame(
            self.tab_mapping
        )

        self.mapping_scroll.pack(
            fill=tk.BOTH,
            expand=True,
            padx=20,
            pady=20
        )

        self.mapping_frame = self.mapping_scroll.scrollable_frame

    def refresh_mapping(self):

        for w in self.mapping_frame.winfo_children():
            w.destroy()

        self.group_vars = {}
        self.repl_vars = {}

        if self.df is None:
            return

        ttk.Label(
            self.mapping_frame,
            text="Input column",
            width=40
        ).grid(row=0, column=0, sticky="w", padx=5, pady=(0, 8))

        ttk.Label(
            self.mapping_frame,
            text="Group / condition",
            width=20
        ).grid(row=0, column=1, sticky="w", padx=5, pady=(0, 8))

        ttk.Label(
            self.mapping_frame,
            text="Replicate ID",
            width=20
        ).grid(row=0, column=2, sticky="w", padx=5, pady=(0, 8))

        for i, col in enumerate(self.df.columns, start=1):

            ttk.Label(
                self.mapping_frame,
                text=col,
                width=40
            ).grid(row=i, column=0, sticky="w", padx=5, pady=2)

            g = tk.StringVar(
                value=col.split("_")[0]
            )

            r = tk.StringVar(
                value=col.split("_")[-1]
            )

            self.group_vars[col] = g
            self.repl_vars[col] = r

            ttk.Entry(
                self.mapping_frame,
                textvariable=g,
                width=20
            ).grid(row=i, column=1, sticky="w", padx=5, pady=2)

            ttk.Entry(
                self.mapping_frame,
                textvariable=r,
                width=20
            ).grid(row=i, column=2, sticky="w", padx=5, pady=2)

    def build_contrasts(self):

        ttk.Button(
            self.tab_contrast,
            text="Update contrasts from mapping",
            command=self.refresh_contrasts
        ).pack(anchor="w", padx=20, pady=(20, 5))

        self.contrast_scroll = ScrollableFrame(
            self.tab_contrast
        )

        self.contrast_scroll.pack(
            fill=tk.BOTH,
            expand=True,
            padx=20,
            pady=20
        )

        self.contrast_frame = self.contrast_scroll.scrollable_frame

    def refresh_contrasts(self):

        for w in self.contrast_frame.winfo_children():
            w.destroy()

        self.contrast_vars = {}

        if self.df is None:
            return

        groups = []

        for col in self.df.columns:

            group_name = self.group_vars[col].get().strip()

            if group_name and group_name not in groups:

                groups.append(group_name)

        groups = sorted(groups)

        row = 0

        for i in range(len(groups)):

            for j in range(i + 1, len(groups)):

                A = groups[i]
                B = groups[j]

                var = tk.BooleanVar(
                    value=True
                )

                self.contrast_vars[(A, B)] = var

                ttk.Checkbutton(
                    self.contrast_frame,
                    text=f"{A} vs {B}",
                    variable=var
                ).grid(row=row, column=0, sticky="w")

                row += 1

    def build_options(self):

        frm = self.tab_options

        ttk.Label(
            frm,
            text="Statistical engine"
        ).pack(anchor="w", padx=20, pady=10)

        ttk.Combobox(
            frm,
            textvariable=self.engine_var,
            values=[
                "Stratified Wilcoxon",
                "Linear Mixed Model",
                "Both"
            ],
            state="readonly"
        ).pack(anchor="w", padx=20)

        ttk.Checkbutton(
            frm,
            text="Log-transform values",
            variable=self.log_var
        ).pack(anchor="w", padx=20, pady=20)

    def build_run(self):

        ttk.Button(
            self.tab_run,
            text="Run Analysis",
            command=self.run_analysis
        ).pack(anchor="w", padx=20, pady=20)

        output_frame = ttk.Frame(
            self.tab_run
        )

        output_frame.pack(
            fill=tk.BOTH,
            expand=True,
            padx=20,
            pady=20
        )

        self.output = tk.Text(
            output_frame,
            width=180,
            height=40,
            wrap="none"
        )

        output_y_scroll = ttk.Scrollbar(
            output_frame,
            orient="vertical",
            command=self.output.yview
        )

        output_x_scroll = ttk.Scrollbar(
            output_frame,
            orient="horizontal",
            command=self.output.xview
        )

        self.output.configure(
            yscrollcommand=output_y_scroll.set,
            xscrollcommand=output_x_scroll.set
        )

        self.output.grid(
            row=0,
            column=0,
            sticky="nsew"
        )

        output_y_scroll.grid(
            row=0,
            column=1,
            sticky="ns"
        )

        output_x_scroll.grid(
            row=1,
            column=0,
            sticky="ew"
        )

        output_frame.rowconfigure(
            0,
            weight=1
        )

        output_frame.columnconfigure(
            0,
            weight=1
        )

    def run_analysis(self):

        if self.df is None:

            messagebox.showerror(
                "No input file",
                "Please open a CSV or Excel file first."
            )

            return

        col_to_group = {
            c: self.group_vars[c].get().strip()
            for c in self.df.columns
        }

        col_to_repl = {
            c: self.repl_vars[c].get().strip()
            for c in self.df.columns
        }

        if any(not group for group in col_to_group.values()):

            messagebox.showerror(
                "Missing group name",
                "All columns must have a group name in the Mapping tab."
            )

            return

        if any(not repl for repl in col_to_repl.values()):

            messagebox.showerror(
                "Missing replicate ID",
                "All columns must have a replicate ID in the Mapping tab."
            )

            return

        df_long = reshape_long(
            self.df,
            col_to_group,
            col_to_repl
        )

        contrasts = []

        current_groups = set(col_to_group.values())

        for (A, B), var in self.contrast_vars.items():

            if var.get():

                if A not in current_groups or B not in current_groups:

                    messagebox.showerror(
                        "Outdated contrasts",
                        "The contrast list does not match the Mapping tab. Click 'Update contrasts from mapping' and run the analysis again."
                    )

                    return

                contrasts.append((A, B))

        if not contrasts:

            messagebox.showerror(
                "No contrast selected",
                "Please select at least one contrast before running the analysis."
            )

            return

        all_results = []

        engine = self.engine_var.get()

        # ====================================================
        # STRATIFIED WILCOXON
        # ====================================================

        if engine in [
            "Stratified Wilcoxon",
            "Both"
        ]:

            rows = []

            raw_p = []

            for A, B in contrasts:

                U, p = van_elteren_test(
                    df_long,
                    A,
                    B
                )

                raw_p.append(p)

                rows.append({

                    "Contrast":
                        f"{A} vs {B}",

                    "Method":
                        "Stratified Wilcoxon",

                    "p_raw":
                        p

                })

            raw_p_array = np.array(raw_p)

            valid = np.isfinite(raw_p_array)

            holm_adj = np.full(len(raw_p), np.nan)
            bonf_adj = np.full(len(raw_p), np.nan)
            fdr_adj = np.full(len(raw_p), np.nan)

            if np.sum(valid) > 0:

                _, p_holm, _, _ = multipletests(
                    raw_p_array[valid],
                    method="holm"
                )

                _, p_bonf, _, _ = multipletests(
                    raw_p_array[valid],
                    method="bonferroni"
                )

                _, p_fdr, _, _ = multipletests(
                    raw_p_array[valid],
                    method="fdr_bh"
                )

                holm_adj[valid] = p_holm
                bonf_adj[valid] = p_bonf
                fdr_adj[valid] = p_fdr

            for i in range(len(rows)):

                rows[i]["p_bonferroni"] = bonf_adj[i]
                rows[i]["p_holm"] = holm_adj[i]
                rows[i]["p_fdr_bh"] = fdr_adj[i]

            wilcox_df = pd.DataFrame(rows)

            all_results.append(wilcox_df)

        # ====================================================
        # LMM
        # ====================================================

        if engine in [
            "Linear Mixed Model",
            "Both"
        ]:

            lmm_df = run_lmm(
                df_long,
                contrasts,
                log_transform=self.log_var.get()
            )

            all_results.append(lmm_df)

        results = pd.concat(
            all_results,
            ignore_index=True
        )

        if self.input_path is not None:

            output_path = (
                self.input_path.parent
                / f"{self.input_path.stem}_statistical-analysis.csv"
            )

        else:

            output_path = Path("statistical-analysis.csv")

        results.to_csv(
            output_path,
            index=False
        )

        self.output.delete("1.0", tk.END)

        self.output.insert(
            tk.END,
            results.to_string(index=False)
        )

        messagebox.showinfo(
            "Done",
            f"Analysis completed. Results saved as:\n{output_path}"
        )


if __name__ == "__main__":

    app = HierarchicalStatisticalAnalyzerGUI()

    app.mainloop()
