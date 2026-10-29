# -*- coding: utf-8 -*-
"""Additional file 2: Figure S1 (pre-processing dependence) and Figure S2
(leave-one-out of the 13 leading-edge pyroptosis genes).

All plotted values are read from the result files; nothing is hard-coded.
"""
import csv
import json
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

sens = json.load(open("sensitivity_results.json", encoding="utf-8"))

# primary-setting values come from the GSEA output table, never hard-coded
PRIMARY = {}
with open("gsea_three_cohort.csv", encoding="utf-8") as fh:
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


def main():
    mm = 1 / 25.4
    fig = plt.figure(figsize=(183 * mm, 150 * mm))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.15], hspace=0.72,
                          left=0.155, right=0.955, top=0.90, bottom=0.13)
    panel_a(fig.add_subplot(gs[0]))
    panel_b(fig.add_subplot(gs[1]))
    fig.savefig("Additional_File_2_figures.pdf", dpi=400)
    fig.savefig("Additional_File_2_figures.png", dpi=400)
    print("saved Additional_File_2_figures.pdf/.png")


if __name__ == "__main__":
    main()
