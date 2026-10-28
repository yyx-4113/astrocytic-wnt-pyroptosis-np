# -*- coding: utf-8 -*-
"""
Generate Figure 10 (Structural Analysis) for the honest manuscript.
Panels:
  A. Four PDBFixer-remediated structures (PDB ID, target, MW) -- prepared & inspected.
  B. Published experimental binding affinities (Kd, nM) from literature.
  C. Molecular-dynamics status: 0/4 completed; no new binding energies computed.
No 3D surface rendering and no novel structural results are claimed (transparency).
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

OUT = r"D:\省自然\生信分析\results\Figure10_Structural_Analysis.png"

structs = [
    ("7PZC", "NLRP3 NACHT domain", "118 kDa"),
    ("6N9O", "CTNNB1 / beta-catenin", "86 kDa"),
    ("1JDH", "TNKS2 tankyrase-2", "71 kDa"),
    ("3UA9", "GSDMD pore", "53 kDa"),
]

fig, axes = plt.subplots(1, 3, figsize=(13, 4.4))

# ---- Panel A: structure cards ----
ax = axes[0]
ax.axis("off")
ax.set_title("A. PDB structures prepared\n(PDBFixer-remediated, inspected)", fontsize=10, weight="bold")
for i, (pdb, tgt, mw) in enumerate(structs):
    y = 0.86 - i * 0.23
    ax.add_patch(FancyBboxPatch((0.04, y - 0.14), 0.92, 0.19,
                 boxstyle="round,pad=0.01", fc="#eef3fb", ec="#3b6ea5", lw=1.2,
                 transform=ax.transAxes))
    ax.text(0.08, y - 0.03, pdb, fontsize=10.5, weight="bold", transform=ax.transAxes, color="#1f3b66")
    ax.text(0.34, y - 0.01, tgt, fontsize=8.5, transform=ax.transAxes)
    ax.text(0.34, y - 0.10, "MW " + mw, fontsize=8.5, color="#555", transform=ax.transAxes)
ax.text(0.04, 0.02, "All four remediated with PDBFixer (Eastman et al., 2017);\n"
                     "3D ligand-binding-site surface rendering not performed.",
        fontsize=7.5, color="#444", transform=ax.transAxes)

# ---- Panel B: literature binding affinities ----
ax = axes[1]
ax.set_title("B. Published binding affinities", fontsize=10, weight="bold")
labels = ["MCC950-NLRP3", "XAV939-TNKS2"]
kd = [5.0, 46.0]
bars = ax.barh(labels, kd, color=["#c0392b", "#2980b9"])
ax.set_xlabel("Kd (nM, literature)")
ax.invert_yaxis()
for b, v in zip(bars, kd):
    ax.text(b.get_width() + 1.0, b.get_y() + b.get_height() / 2, f"{v:g} nM",
            va="center", fontsize=9)
ax.set_xlim(0, 62)
ax.text(0.0, -0.32, "Sources: Coll et al., 2019; Huang et al., 2009.",
        fontsize=7, color="#555", transform=ax.transAxes)

# ---- Panel C: MD status ----
ax = axes[2]
ax.axis("off")
ax.set_title("C. Molecular-dynamics status", fontsize=10, weight="bold")
ax.text(0.05, 0.72, "MD completed: 0 / 4", fontsize=12, weight="bold",
        color="#c0392b", transform=ax.transAxes)
ax.text(0.05, 0.50,
        "No new binding-energy calculations performed.\n"
        "Reported energies are from published\n"
        "experimental literature\n"
        "(Coll et al., 2019; Huang et al., 2009).",
        fontsize=9, transform=ax.transAxes)
ax.text(0.05, 0.10,
        "Transparency instrument only:\nno novel structural results claimed.",
        fontsize=7.5, color="#444", style="italic", transform=ax.transAxes)

fig.suptitle("Figure 10. Structural Analysis of Drug-Target Complexes", fontsize=13, weight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(OUT, dpi=160, bbox_inches="tight")
print("Saved", OUT)
