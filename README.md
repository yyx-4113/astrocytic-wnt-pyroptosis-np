# astrocytic-wnt-pyroptosis-np

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
