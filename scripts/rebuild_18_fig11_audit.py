# -*- coding: utf-8 -*-
"""
Generate Figure 11 (Reproducibility & Limitation Self-Audit).
A compliance matrix: rows = analytical modules (keyed to Methods sections),
columns = five reproducibility / scientific-strictness dimensions.
Cells: green = addressed (check), amber = partial / limitation (triangle),
red = absent / not applicable (cross). A right-hand module-compliance column
and a bottom dimension-compliance row report weighted % scores.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

OUT = r"D:\省自然\生信分析\results\Figure11_Reproducibility_Audit.png"

modules = [
    "Bulk DEG (GSE217932)",
    "WGCNA",
    "GSEA",
    "scRNA-seq (GSE134003)",
    "CellChat",
    "Deconvolution",
    "PPI / Drug",
    "GRN / RWR",
    "Structural / MD",
]
dims = ["Sample\nsize", "Independent\nreplication", "Multiple-\ntesting", "Negative/\nnull report", "Code/data\nsharing"]

# 2 = green (addressed), 1 = amber (partial/limitation), 0 = red (absent/N.A.)
M = np.array([
    [1, 0, 2, 2, 2],  # Bulk DEG
    [1, 0, 1, 2, 2],  # WGCNA
    [1, 0, 2, 2, 2],  # GSEA
    [2, 1, 2, 2, 2],  # scRNA-seq
    [2, 1, 2, 1, 2],  # CellChat
    [1, 0, 1, 2, 2],  # Deconvolution
    [0, 0, 1, 2, 2],  # PPI / Drug
    [1, 0, 0, 2, 2],  # GRN / RWR
    [0, 0, 0, 2, 2],  # Structural / MD
], dtype=int)

H, W = M.shape
mod_score = M.mean(axis=1) / 2 * 100      # per-module compliance %
dim_score = M.mean(axis=0) / 2 * 100      # per-dimension compliance %

colors = {0: (0.85, 0.30, 0.30), 1: (0.95, 0.75, 0.20), 2: (0.30, 0.70, 0.40)}
img = np.zeros((H, W, 3))
for i in range(H):
    for j in range(W):
        img[i, j] = colors[M[i, j]]

sym = {2: "✓", 1: "▲", 0: "✗"}

fig, ax = plt.subplots(figsize=(10.5, 6.8))
ax.imshow(img, aspect="auto", vmin=0, vmax=1)
for i in range(H):
    for j in range(W):
        ax.text(j, i, sym[M[i, j]], ha="center", va="center", fontsize=15,
                color="white" if M[i, j] != 1 else "black", weight="bold")
ax.set_xticks(range(W))
ax.set_xticklabels(dims, fontsize=9)
ax.set_yticks(range(H))
ax.set_yticklabels(modules, fontsize=9)
ax.set_xlim(-1.8, W + 1.5)
ax.set_ylim(H + 1.6, -0.5)

# summary row (dimension compliance) below grid
for j in range(W):
    ax.text(j, H + 0.5, f"{dim_score[j]:.0f}%", ha="center", fontsize=9.5,
            weight="bold", color="#222")
ax.text(-1.5, H + 0.5, "Dim %", fontsize=8.5, weight="bold", color="#222")

# summary column (module compliance) right of grid
for i in range(H):
    ax.text(W + 0.8, i, f"{mod_score[i]:.0f}%", ha="center", va="center",
            fontsize=9.5, weight="bold", color="#222")
ax.text(W + 0.8, H + 0.5, "Mod %", fontsize=8.5, weight="bold", color="#222")

ax.set_title("Figure 11. Reproducibility and Limitation Self-Audit", fontsize=13, weight="bold")
leg = [
    Patch(facecolor=colors[2], label="Addressed (green check)"),
    Patch(facecolor=colors[1], label="Partial / discussed limitation (amber triangle)"),
    Patch(facecolor=colors[0], label="Absent / N.A. (red cross)"),
]
ax.legend(handles=leg, loc="upper center", bbox_to_anchor=(0.5, -0.06),
          ncol=3, fontsize=8, frameon=False)
fig.tight_layout()
fig.savefig(OUT, dpi=160, bbox_inches="tight")
print("Saved", OUT)
print("Module compliance %:", {m: round(s, 0) for m, s in zip(modules, mod_score)})
print("Dimension compliance %:", {d.replace(chr(10), ' '): round(s, 0) for d, s in zip(dims, dim_score)})
