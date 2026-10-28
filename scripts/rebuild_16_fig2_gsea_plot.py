# -*- coding: utf-8 -*-
"""Regenerate Figure 2 GSEA plot from the GENUINE GSE217932 ranking.
The previous Figure2_Wnt_pathway.png showed a fabricated positive enrichment
(NES=1.68, FDR=0.012). The honest recomputation gives NES ~= -1.6 (depleted
in the CCI/injury direction, nominal p~0.02, not significant after correction).
This script plots the running enrichment score and ranked log2FC with the
TRUE values so the figure matches the corrected text.
"""
import csv, math, random
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEG = r"D:\省自然\生信分析\results\GSE217932_DEGs.csv"
OUT = r"D:\省自然\生信分析\results\Figure2_Wnt_pathway.png"

HALLMARK = {
    "MYC","CCND1","CTNNB1","LEF1","TCF7","FOSL1","JAG2","NOTUM","PPP2R5B","WIF1",
    "AXIN2","NKD1","TBL1XR1","CD44","KLF4","PAX2","CXXC4","SOX17","TCF7L1","TLE1",
    "USP34","GLIS2","SFRP1","DKK1","DKK4","FZD7","CER1","FST","ASCL2","LGR4","LRP6",
    "OLFM1","PITX2","SPRY2","T","WNT11","WNT3","WNT5A"
}
HUP = {h.upper() for h in HALLMARK}

sym_lfc = {}
with open(DEG, newline="") as f:
    for row in csv.DictReader(f):
        try:
            s = row["symbol"].strip(); v = float(row["log2FC"])
        except: 
            continue
        if s and s != "NA":
            sym_lfc.setdefault(s, []).append(v)
mean_lfc = {s: sum(v)/len(v) for s, v in sym_lfc.items()}
rank = sorted(mean_lfc, key=lambda g: -mean_lfc[g])
S = [g for g in rank if g.upper() in HUP]
Sset = {g.upper() for g in S}

# running enrichment score
Nh = len(S); N = len(rank)
hit_w = [abs(mean_lfc[g]) if g.upper() in Sset else 0.0 for g in rank]
sum_hit = sum(hit_w) or 1.0
running = np.zeros(N)
cum = 0.0
for i, g in enumerate(rank):
    if g.upper() in Sset:
        cum += hit_w[i] / sum_hit
    else:
        cum -= 1.0 / (N - Nh)
    running[i] = cum
ES_obs = running[np.argmax(np.abs(running))]
# NES via permutation
random.seed(7); nperm = 500
es_perm = np.empty(nperm)
for p in range(nperm):
    perm = rank[:]; random.shuffle(perm)
    hw = [abs(mean_lfc[g]) if g.upper() in Sset else 0.0 for g in perm]
    sh = sum(hw) or 1.0; c = 0.0; r = np.zeros(N)
    for i, g in enumerate(perm):
        if g.upper() in Sset: c += hw[i]/sh
        else: c -= 1.0/(N-Nh)
        r[i]=c
    es_perm[p] = r[np.argmax(np.abs(r))]
nes = ES_obs / np.mean(np.abs(es_perm))
pval = (np.sum(es_perm <= ES_obs)+1)/(nperm+1) if ES_obs < 0 else (np.sum(es_perm >= ES_obs)+1)/(nperm+1)

print(f"ES={ES_obs:.3f} NES={nes:.3f} nominal p={pval:.3f} (n_hallmark={Nh})")

# ---- plot ----
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True,
                                gridspec_kw={"height_ratios": [1, 1.3]})
# top: ranked log2FC
lfc_arr = np.array([mean_lfc[g] for g in rank])
is_hit = np.array([g.upper() in Sset for g in rank])
x = np.arange(N)
ax1.scatter(x[~is_hit], lfc_arr[~is_hit], s=0.5, c="lightgray", label="Background")
ax1.scatter(x[is_hit], lfc_arr[is_hit], s=4, c="crimson", label="Wnt/β-cat hallmark")
ax1.axhline(0, color="k", lw=0.5)
ax1.set_ylabel("log2FC (CCI vs Sham)")
ax1.set_title("Figure 2. GSEA — Hallmark Wnt/β-catenin (GSE217932, genuine)")
ax1.legend(fontsize=7, loc="upper right")

# bottom: running ES
ax2.plot(x, running, color="navy", lw=1.0)
ax2.axhline(0, color="k", lw=0.5)
ymin, ymax = running.min(), running.max()
ax2.vlines([0, N], ymin, ymax, color="gray", lw=0.5)
ax2.set_ylim(ymin*1.1, ymax*1.1)
ax2.set_xlabel("Ranked genes")
ax2.set_ylabel("Running ES")
ax2.set_title(f"Running enrichment score  (NES ≈ {nes:.2f}, nominal p ≈ {pval:.2f}; "
              f"depleted in CCI direction, not significant after correction)")
# mark hallmark hit positions lightly
ax2.vlines(np.where(is_hit)[0], ymin, ymax, color="red", alpha=0.05, lw=0.3)

fig.tight_layout()
fig.savefig(OUT, dpi=160, bbox_inches="tight")
print("Saved", OUT)
