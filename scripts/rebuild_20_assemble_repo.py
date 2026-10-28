# -*- coding: utf-8 -*-
"""
Assemble a CLEAN public reproducibility repository at
  D:/省自然/生信分析/astrocytic-wnt-pyroptosis-np/
from the honest manuscript + genuine artifacts only.

Included:
  - manuscript_final.md / manuscript_final_with_figures.docx
  - figures/  : the 12 genuine PNGs
  - data/     : genuine processed CSVs + JSON summaries (NO MR / cross-dataset / GSE186237)
  - scripts/  : rebuild_*.py (the honest regeneration pipeline)
  - README.md, LICENSE (MIT), CITATION.cff, MANIFEST.md
Excluded: any MR, cross-dataset, or GSE186237-derived file (retired/withdrawn).
"""
import os, shutil, glob

BASE = r"D:/省自然/生信分析"
RESULTS = os.path.join(BASE, "results")
PKG = os.path.join(BASE, "submission", "REPRODUCIBILITY_PACKAGE")
REPO = os.path.join(BASE, "astrocytic-wnt-pyroptosis-np")

FIGS = [
    "Figure0_Graphical_Abstract.png", "Figure1_DEG_overview.png",
    "Figure2_Wnt_pathway.png", "Figure3_Pyroptosis_Wnt_crosstalk.png",
    "Figure4_scRNAseq_cell_atlas.png", "Figure5_CellChat_Pseudotime.png",
    "Figure6_Deconvolution.png", "Figure7_PPI_Hub_TCF4.png",
    "Figure8_Drug_ThreeTier_Classification.png", "Figure9_GRN_VirtualKO.png",
    "Figure10_Structural_Analysis.png", "Figure11_Reproducibility_Audit.png",
]
DATA_CSV = [
    "GSE217932_DEGs.csv", "MC5_GeneModuleMembership.csv", "MC5_WGCNA_Modules.csv",
    "MC2_Drug_ThreeTier.csv", "Figure5_communication_SNI.csv",
    "Figure6_deconv_GSE217932.csv", "Figure7_HubGenes.csv",
    "Figure9_RWR_Scores.csv", "Figure9_VirtualKO.csv", "Figure9_Perturbation_OE_KD.csv",
    "Figure8_DrugConnectivity.csv", "Figure8_TargetEnrichment.csv", "MC4_Pathway_Comparison.csv",
]
DATA_JSON = [
    "Figure7_summary.json", "Figure8_summary.json", "Figure9_summary.json",
    "Figure10_summary.json", "MC5_WGCNA_Summary.json",
]

def mkdir(p):
    os.makedirs(p, exist_ok=True)

mkdir(REPO)
for sub in ["figures", "data", "scripts"]:
    mkdir(os.path.join(REPO, sub))

# --- figures ---
print("== figures ==")
for f in FIGS:
    src = os.path.join(RESULTS, f)
    if os.path.exists(src):
        shutil.copy2(src, os.path.join(REPO, "figures", f))
        print("  copied", f)
    else:
        print("  MISSING", f)

# --- data ---
print("== data ==")
for f in DATA_CSV + DATA_JSON:
    src = os.path.join(RESULTS, f)
    if os.path.exists(src):
        shutil.copy2(src, os.path.join(REPO, "data", f))
        print("  copied", f)
    else:
        print("  MISSING", f)

# --- scripts (rebuild_*.py only) ---
print("== scripts ==")
n = 0
for src in glob.glob(os.path.join(BASE, "rebuild_*.py")):
    shutil.copy2(src, os.path.join(REPO, "scripts", os.path.basename(src)))
    n += 1
print("  copied", n, "rebuild_*.py scripts")

# --- manuscript + manifest ---
shutil.copy2(os.path.join(BASE, "submission", "manuscript_final.md"), os.path.join(REPO, "manuscript_final.md"))
shutil.copy2(os.path.join(BASE, "submission", "manuscript_final_with_figures.docx"), os.path.join(REPO, "manuscript_final_with_figures.docx"))
shutil.copy2(os.path.join(PKG, "MANIFEST.md"), os.path.join(REPO, "MANIFEST.md"))
print("copied manuscript + MANIFEST")

# --- LICENSE (MIT) ---
LICENSE = """MIT License

Copyright (c) 2026 Yang Yong-Xin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
with open(os.path.join(REPO, "LICENSE"), "w", encoding="utf-8") as fh:
    fh.write(LICENSE)

# --- CITATION.cff ---
CITATION = """cff-version: 1.2.0
message: "If you use this work, please cite the manuscript and this repository."
title: "Computational re-analysis of convergent Wnt/beta-catenin and pyroptosis signaling in neuropathic pain, with single-cell resolution of cell-type specificity"
authors:
  - family-names: Yang
    given-names: Yong-Xin
    orcid: "https://orcid.org/0009-0004-9698-6552"
repository-code: "https://github.com/yyx-4113/astrocytic-wnt-pyroptosis-np"
license: MIT
keywords:
  - neuropathic pain
  - Wnt/beta-catenin
  - pyroptosis
  - bioinformatics
  - reproducible research
"""
with open(os.path.join(REPO, "CITATION.cff"), "w", encoding="utf-8") as fh:
    fh.write(CITATION)

# --- README.md ---
README = """# astrocytic-wnt-pyroptosis-np

**Honest computational re-analysis of convergent Wnt/β-catenin and pyroptosis
signaling in neuropathic pain, with single-cell resolution of cell-type
specificity — a hypothesis-generating study.**

Single-author manuscript (Yang Yong-Xin). All analyses use publicly available,
de-identified datasets; no new experimental data were generated.

## Important honesty disclosures (read before reuse)
- **Single bulk cohort.** Only one bulk RNA-seq cohort exists (GSE217932, rat
  SNI/CCI spinal cord, n=4 per condition). There is **no independent bulk
  cohort**; all cross-dataset replication claims from earlier drafts were
  withdrawn. `GSE186237` was found to be a byte-identical duplicate of
  `GSE217932` (identical MD5) and is excluded.
- **Mendelian-randomization analysis withdrawn.** The MR component of earlier
  drafts was removed; no MR results, instruments, or causal estimates are
  reported.
- **WGCNA module preservation not assessed.** The earlier `Z_summary = 1.0`
  preservation file was forged (self-comparison) and has been deleted. Module
  *assignment* is reported descriptively as exploratory.
- **Negative/null findings are reported transparently** (see the manuscript
  and Figure 11 self-audit).

## Repository structure
```
manuscript_final.md                  # source manuscript (honest)
manuscript_final_with_figures.docx   # submission-ready Word file (12 figures embedded)
figures/                             # 12 genuine PNG figures (Figure0–Figure11)
data/                                # genuine processed CSVs + JSON summaries
scripts/                             # rebuild_*.py — honest regeneration pipeline
MANIFEST.md                          # file inventory + checksums + provenance
LICENSE                              # MIT
CITATION.cff                         # citation metadata
```

## Data sources
- Bulk RNA-seq: GSE217932 (GEO; https://www.ncbi.nlm.nih.gov/geo/)
- Single-cell RNA-seq: GSE134003 (GEO)
- Protein structures: PDB 7PZC, 6N9O, 1JDH, 3UA9 (RCSB PDB)
- PPI: STRING; Drug targets: DrugBank

## How to reproduce
Run the scripts in `scripts/` (Python 3, numpy, matplotlib, python-docx):
```
python scripts/rebuild_19_build_docx.py   # rebuilds the .docx from the .md + figures
python scripts/rebuild_16_fig2_gsea_plot.py
python scripts/rebuild_17_fig10_structural.py
python scripts/rebuild_18_fig11_audit.py
...
```
Each script regenerates its figure/data from the genuine inputs.

## License
MIT — see LICENSE.
"""
with open(os.path.join(REPO, "README.md"), "w", encoding="utf-8") as fh:
    fh.write(README)

print("\nRepo assembled at:", REPO)
print("Tree:")
for root, dirs, files in os.walk(REPO):
    rel = os.path.relpath(root, REPO)
    if rel == ".":
        rel = ""
    for f in sorted(files):
        print(f"  {rel}/{f}")
