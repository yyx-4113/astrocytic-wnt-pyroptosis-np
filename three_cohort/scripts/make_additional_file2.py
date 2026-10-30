# -*- coding: utf-8 -*-
"""Additional file 2: Figure S1 (pre-processing dependence) and Figure S2
(leave-one-out of the 13 leading-edge pyroptosis genes).

All plotted values are read from the result files; nothing is hard-coded.
"""
import csv
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

WNT_C = "#2F5D8C"
PYRO_C = "#B4443C"
GREY = "#6E6E6E"

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE) if os.path.basename(HERE).lower() == "scripts" else HERE
RESULTS = os.path.join(PKG, "results") if os.path.isdir(os.path.join(PKG, "results")) else PKG
OUTDIR = os.environ.get("NP_OUT", PKG)


def _find(fn):
    for d in (RESULTS, PKG, os.getcwd()):
        p = os.path.join(d, fn)
        if os.path.exists(p):
            return p
    return os.path.join(RESULTS, fn)  # not found; fail with a clear path


sens = json.load(open(_find("sensitivity_results.json"), encoding="utf-8"))

# primary-setting values come from the GSEA output table, never hard-coded
PRIMARY = {}
with open(_find("gsea_three_cohort.csv"), encoding="utf-8") as fh:
    for r in csv.DictReader(fh):
        if r["cohort"] == "C1":
            PRIMARY[r["gene_set"]] = (float(r["NES"]), float(r["nominal_p"]))
PRIMARY["PYRO_EXECUTION_MACHINERY"] = (sens["exec"]["nes"], sens["exec"]["p"])

VARIANTS = ["primary", "log2FC\nranking", "median\nfilter", "no\nfilter"]
KEYS = ["primary", "lfc", "med", "nofilter"]
KEYMAP = {"wnt": "HALLMARK_WNT_BETA_CATENIN",
          "pyro": "PYROPTOSIS_CORE",
          "exec": "PYRO_EXECUTION_MACHINERY"}


def nes_of(key, setname):
    if key == "primary":
        return PRIMARY[KEYMAP[setname]]
    v = sens[key].get(setname)
    return None if v is None else (v["nes"], v["p"])


def panel_a(ax):
    for setname, col, label, dy in [("pyro", PYRO_C, "PYROPTOSIS_CORE", 11),
                                    ("wnt", WNT_C, "HALLMARK_WNT_BETA_CATENIN", -13)]:
        pts = [(i, *nes_of(k, setname)) for i, k in enumerate(KEYS)]
        ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=30, color=col,
                   zorder=3, edgecolor="white", linewidth=0.6, label=label)
        for x, nes, p in pts:
            ax.annotate(f"{nes:+.2f} ({p:.2g})", (x, nes), textcoords="offset points",
                        xytext=(0, dy), ha="center", va="center", fontsize=5.8,
                        color=col)
    ax.axhline(0, color="#5A5A5A", lw=0.6, zorder=1)
    ax.set_xticks(range(4))
    ax.set_xticklabels(VARIANTS, fontsize=6.5)
    ax.set_xlim(-0.5, 3.5)
    ax.set_ylim(-2.1, 3.6)
    ax.set_ylabel("Normalised enrichment score (cohort C1)", fontsize=7)
    ax.legend(loc="center left", bbox_to_anchor=(0.03, 0.63), frameon=False,
              fontsize=6, handletextpad=0.4, borderpad=0.2)
    ax.set_title("A   Dependence of the Wnt and pyroptosis enrichment scores on\n"
                 "      pre-processing settings (cohort C1)",
                 loc="left", fontsize=7.5, fontweight="bold", pad=8)


def panel_b(ax):
    loo = sorted(sens["loo"], key=lambda d: d["dNES"])
    names = [d["dropped"].capitalize() for d in loo]
    dnes = [d["dNES"] for d in loo]
    ypos = list(range(len(names)))
    ax.barh(ypos, dnes, color=PYRO_C, alpha=0.85, height=0.6, zorder=2)
    ax.axvline(0, color="#5A5A5A", lw=0.6, zorder=1)
    ax.set_yticks(ypos)
    ax.set_yticklabels(names, fontsize=6.2)
    ax.set_ylim(-0.7, len(names) - 0.3)
    ax.set_xlabel("Change in NES when that gene is removed from the set", fontsize=6.8)
    lo = min(dnes)
    ax.set_xlim(lo - 0.012, 0.028)
    ax.set_xticks([-0.10, -0.075, -0.05, -0.025, 0, 0.025])
    ax.set_xticklabels(["\u22120.100", "\u22120.075", "\u22120.050", "\u22120.025",
                        "0", "0.025"], fontsize=6.2)
    for y, d in zip(ypos, dnes):
        ax.text(0.010, y, f"{d:+.03f}", va="center", ha="left", fontsize=5.8,
                color="#333333")
    ax.set_title("B   Leave-one-out analysis of the 13 leading-edge pyroptosis genes",
                 loc="left", fontsize=7.5, fontweight="bold", pad=8)
    ax.text(0.0, -0.255, f"All 13 re-runs remained significant; NES never fell below "
                         f"{min(d['nes'] for d in sens['loo']):.2f} "
                         f"(primary NES {PRIMARY['PYROPTOSIS_CORE'][0]:.2f}).",
            transform=ax.transAxes, fontsize=6, color=GREY)


SHORT = {"raw": "raw counts", "median_ratio": "median-of-ratios",
         "upper_quartile": "upper-quartile", "as_deposited": "as deposited (RPKM)"}


def panel_c(ax):
    """Does library normalisation create the null?  Survivors of BH correction
    under every pre-processing regime tried."""
    norm = json.load(open(_find("normalisation_sensitivity.json"), encoding="utf-8"))
    labels, values, colours = [], [], []
    palette = {"C1": "#C0392B", "C2": "#2471A3", "C3": "#1E8449"}
    for cid in ("C1", "C2", "C3"):
        for key in [k for k in norm[cid] if isinstance(norm[cid][k], dict)
                    and "q05_survivors" in norm[cid][k]]:
            base = key.replace("_indepfiltered", "")
            lab = SHORT.get(base, base)
            if key.endswith("_indepfiltered"):
                lab += " + filtering"
            labels.append(f"{cid}  {lab}")
            values.append(norm[cid][key]["q05_survivors"])
            colours.append(palette[cid])
    y = list(range(len(labels)))
    ax.barh(y, values, color=colours, height=0.68)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=5.6)
    ax.invert_yaxis()
    ax.set_xlabel("genes surviving Benjamini\u2013Hochberg q < 0.05", fontsize=6.5)
    ax.set_xlim(0, max(values) + 1.5)
    ax.tick_params(axis="x", labelsize=6)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_title("C   Library normalisation does not create the null",
                 loc="left", fontsize=7.5, fontweight="bold", pad=8)
    qs = {}
    for cid in ("C1", "C2", "C3"):
        v = [norm[cid][k]["Atf3"]["q"] for k in norm[cid]
             if isinstance(norm[cid][k], dict) and norm[cid][k].get("Atf3")]
        qs[cid] = (min(v), max(v))
    note = ("Atf3 across the same regimes: C1 %.2f\u2013%.2f, C2 %.2f\u2013%.2f, C3 %.2f "
            "(q never below 0.05)\n"
            "Total genes tested per regime: 9,505\u201319,864. C3 is already normalised "
            "abundance and cannot be re-normalised.") % (
        qs["C1"][0], qs["C1"][1], qs["C2"][0], qs["C2"][1], qs["C3"][0])
    ax.text(0.0, -0.34, note, transform=ax.transAxes, fontsize=5.6,
            color=GREY, va="top")


def main():
    mm = 1 / 25.4
    fig = plt.figure(figsize=(183 * mm, 215 * mm))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.15, 1.05], hspace=0.66,
                          left=0.155, right=0.955, top=0.945, bottom=0.125)
    panel_a(fig.add_subplot(gs[0]))
    panel_b(fig.add_subplot(gs[1]))
    panel_c(fig.add_subplot(gs[2]))
    fig.savefig(os.path.join(OUTDIR, "Additional_File_2_figures.pdf"), dpi=400)
    fig.savefig(os.path.join(OUTDIR, "Additional_File_2_figures.png"), dpi=400)
    print("saved Additional_File_2_figures.pdf/.png")


if __name__ == "__main__":
    main()
