# -*- coding: utf-8 -*-
"""Task 14 — Honest rebuild of the WGCNA deliverable.

The original MC5_WGCNA_Preservation.csv was FORGED: it compared GSE217932
(reference) against GSE186237 (test), but GSE186237 is a byte-identical
duplicate of GSE217932. That made mean_corr_ref == mean_corr_test,
conn_corr = 1.0, mm_corr = 1.0, and Z_summary = median(0,1,1) = 1.0 for
every module — exactly the artefact Reviewer 1 flagged.

This script rebuilds the deliverable HONESTLY:
  * The genuine per-gene module assignment (gene -> module) computed on the
    single real GSE217932 cohort is retained (MC5_GeneModuleMembership.csv).
  * No cross-dataset Z_summary is reported, because no independent bulk cohort
    exists (the manuscript explicitly states preservation was not assessed).
  * Output: MC5_WGCNA_Modules.csv (Wnt/pyroptosis gene module assignments +
    module sizes) and an honest MC5_WGCNA_Summary.json.
  * The forged MC5_WGCNA_Preservation.csv is removed.
"""
import csv, json, os
from collections import Counter, defaultdict

RES = r"D:\省自然\生信分析\results"
MODFILE = os.path.join(RES, "MC5_GeneModuleMembership.csv")   # genuine, 5000 genes
OLD = os.path.join(RES, "MC5_WGCNA_Preservation.csv")          # forged
NEW = os.path.join(RES, "MC5_WGCNA_Modules.csv")               # new honest file
SUM = os.path.join(RES, "MC5_WGCNA_Summary.json")              # new honest json

# ---- Wnt / pyroptosis gene sets (rat symbols) ----
WNT = ["Wnt1","Wnt2","Wnt2b","Wnt3","Wnt3a","Wnt4","Wnt5a","Wnt5b","Wnt6","Wnt7a",
       "Wnt7b","Wnt8a","Wnt8b","Wnt9a","Wnt9b","Wnt10a","Wnt10b","Wnt11","Wnt16",
       "Fzd1","Fzd2","Fzd3","Fzd4","Fzd5","Fzd6","Fzd7","Fzd8","Fzd9","Fzd10",
       "Lrp5","Lrp6","Ctnnb1","Tcf4","Tcf3","Lef1","Axin1","Axin2","Gsk3b",
       "Dvl1","Dvl2","Dvl3","Sfrp1","Sfrp2","Dkk1","Dkk2","Dkk3","Ror1","Ror2",
       "Ryk","Ccnd1","Myc"]
PYRO = ["Nlrp3","Pycard","Casp1","Gsdmd","Il1b","Il18","Nlrp1","Nlrc4","Aim2",
        "Casp4","Casp5","Gsdme","Il1a","Hmgb1","Tnf","Nfkb1"]
wnt_set, pyro_set = set(WNT), set(PYRO)

# ---- Load genuine module assignment ----
gene_module = {}
with open(MODFILE, newline="") as f:
    for row in csv.DictReader(f):
        gene_module[row["gene"]] = row["module"]
print(f"Loaded {len(gene_module)} gene->module assignments")

# ---- Module sizes (genuine) ----
mod_size = Counter(gene_module.values())
print(f"Modules: {len(mod_size)}; sizes: {dict(sorted(mod_size.items()))}")

# ---- Build honest Wnt/pyro module-assignment table ----
rows = []
for g in sorted(wnt_set | pyro_set):
    if g in gene_module:
        rows.append({
            "gene": g,
            "pathway": "Wnt" if g in wnt_set else "Pyroptosis",
            "module": gene_module[g],
            "module_size": mod_size[gene_module[g]],
        })
print(f"Wnt/pyroptosis genes present in WGCNA set: {len(rows)}")

with open(NEW, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["gene","pathway","module","module_size"])
    w.writeheader()
    w.writerows(rows)
print(f"Wrote {NEW}")

# ---- Honest summary (NO forged Z_summary) ----
n_wnt = sum(1 for g in wnt_set if g in gene_module)
n_pyro = sum(1 for g in pyro_set if g in gene_module)
# distribution of Wnt/pyro genes across modules
mod_counts = Counter(r["module"] for r in rows)
summary = {
    "dataset": "GSE217932",
    "note": ("WGCNA module ASSIGNMENT on the single genuine rat bulk cohort "
             "(n=12, n=4 per primary contrast). Module PRESERVATION across an "
             "independent dataset was NOT assessed because no independent bulk "
             "NP cohort exists; the Langfelder et al. (2011) Z-summary therefore "
             "was not computed. The earlier MC5_WGCNA_Preservation.csv (Z_summary "
             "= 1.0 for all modules) was forged by comparing GSE217932 against its "
             "own duplicate GSE186237 and has been deleted."),
    "n_genes_analyzed": len(gene_module),
    "n_modules": len(mod_size),
    "module_sizes": dict(sorted(mod_size.items())),
    "n_wnt_genes_in_set": n_wnt,
    "n_pyroptosis_genes_in_set": n_pyro,
    "wnt_pyroptosis_module_distribution": dict(sorted(mod_counts.items())),
    "method": "Signed WGCNA, Spearman correlation, soft-threshold beta=6, top 5000 variable genes, min module size=30 (original pipeline_mc5 run on GSE217932 only)",
    "interpretation": "Exploratory; modules are descriptive and not replicated.",
    "caution": "n=12 total (n=4 per contrast) limits reliability; interpret qualitatively."
}
with open(SUM, "w") as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)
print(f"Wrote {SUM}")

# ---- Remove the forged preservation file ----
if os.path.exists(OLD):
    os.remove(OLD)
    print(f"Deleted forged file: {OLD}")
else:
    print(f"(forged file not present: {OLD})")

print("\n=== Task 14 WGCNA rebuild complete ===")
