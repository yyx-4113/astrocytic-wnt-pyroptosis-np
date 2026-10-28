# -*- coding: utf-8 -*-
"""A-task: rebuild Figure 0 (Graphical Abstract) honestly.
Removes: (1) fabricated 'n=2 datasets', MR, and MD-simulation claims;
         (2) depiction of UNMEASURED events (caspase-1 cleavage, GSDMD-N pore)
             as if observed. Reflects the corrected single-cohort + full-transcriptome
         + descriptive-drug-catalog + hypothesis-generating framing."""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patches as mpatches

OUT = r'D:\省自然\生信分析\results\Figure0_Graphical_Abstract.png'
C_WNT='#1f77b4'; C_PYRO='#d62728'; C_ASTRO='#ff7f0e'; C_DRUG='#9467bd'; C_GRAY='#666666'; C_BG='#fafafa'

fig, ax = plt.subplots(1, 1, figsize=(24, 13))
ax.set_xlim(0, 28); ax.set_ylim(0, 14); ax.axis('off'); ax.set_facecolor(C_BG)

ax.text(14, 13.2, "Graphical Abstract: Wnt/β-Catenin and Pyroptosis Signaling in Neuropathic Pain",
        ha='center', fontsize=15, fontweight='bold', color='#222')
ax.text(14, 12.5, "A computational, hypothesis-generating re-analysis  (single bulk RNA-seq cohort + full-transcriptome scRNA-seq + descriptive drug catalog)",
        ha='center', fontsize=8.5, fontstyle='italic', color=C_GRAY)

def box(x, y, w, h, title, lines, ec, tc=C_GRAY):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.15",
                 facecolor='white', edgecolor=ec, linewidth=2))
    ax.text(x+w/2, y+h-0.35, title, ha='center', fontsize=10.5, fontweight='bold', color=ec)
    for i, ln in enumerate(lines):
        ax.text(x+w/2, y+h-0.95-i*0.62, ln, ha='center', fontsize=7.6, color=tc)

# Box 1: NP model -> Wnt (computational inference)
box(0.8, 7.0, 5.2, 3.4, "Neuropathic pain model (CCI/SNI)",
    ["Rodent spinal cord", "Bulk RNA-seq: GSE217932", "Sham vs CCI (n=4/group)",
     "Wnt/β-catenin genes:", "  Wnt2/Wnt7a, Ctnnb1, Tcf4", "  ↑ vs sham (nominal)"],
    C_WNT)
# Box 2: Wnt -> pyroptosis crosstalk (hypothesis)
box(7.0, 7.0, 5.6, 3.4, "Wnt -> pyroptosis crosstalk (hypothesis)",
    ["Inferred from correlation", "network + motif analysis", "NOT demonstrated:",
     "  pathway activation,", "  cell death, or pain mechanism", "(hypothesis-generating only)"],
    C_PYRO)
# Box 3: pyroptosis localization (scRNA, corrected)
box(13.8, 7.0, 5.6, 3.4, "Pyroptosis localization (scRNA, full transcriptome)",
    ["GSE134003 SNI spinal cord", "Pyroptosis genes ↑ in", "  MICROGLIA / MACROPHAGES",
     "Astrocyte pyroptosis:", "  LOW (25% detection),", "  unchanged SNI vs SI"],
    C_ASTRO)
# Box 4: drug catalog
box(20.6, 7.0, 5.6, 3.4, "Descriptive drug candidate catalog",
    ["Tier 1 (candidate targets):", "  MCC950 (NLRP3),", "  disulfiram (GSDMD),",
     "  XAV939 (TNKS2) ...", "Descriptive only —", "no statistical ranking"],
    C_DRUG)

# arrows
ax.add_patch(FancyArrowPatch((6.0,8.7),(7.0,8.7), arrowstyle='->', color=C_WNT, lw=2.5, mutation_scale=18))
ax.add_patch(FancyArrowPatch((12.6,8.7),(13.8,8.7), arrowstyle='->', color=C_PYRO, lw=2.5, mutation_scale=18))
ax.add_patch(FancyArrowPatch((19.4,8.7),(20.6,8.7), arrowstyle='->', color=C_ASTRO, lw=2.5, mutation_scale=18))

# bottom caveat bar
ax.add_patch(FancyBboxPatch((0.8, 2.2), 25.4, 2.0, boxstyle="round,pad=0.15",
             facecolor='#fff3cd', edgecolor=C_GRAY, linewidth=1.2))
ax.text(14, 3.7, "Limitations (transparent reporting)", ha='center', fontsize=10, fontweight='bold', color='#7a5c00')
ax.text(14, 3.1, "• Single bulk cohort, NO independent replication   • scRNA-seq from a different species/model than bulk",
        ha='center', fontsize=7.6, color='#7a5c00')
ax.text(14, 2.7, "• Weak bulk signal (2 genes significant after FDR)   • No MD simulations completed   • All findings are computational predictions requiring validation",
        ha='center', fontsize=7.6, color='#7a5c00')

ax.text(14, 1.0, "Computational predictions only — not measured protein events (no caspase-1 cleavage / GSDMD-N pore assays performed).",
        ha='center', fontsize=7.2, style='italic', color=C_GRAY)

fig.tight_layout()
fig.savefig(OUT, dpi=160, bbox_inches='tight')
print('saved', OUT, flush=True)
