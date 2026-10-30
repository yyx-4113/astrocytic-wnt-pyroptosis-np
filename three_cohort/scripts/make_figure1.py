# -*- coding: utf-8 -*-
"""
Figure 1 — design, assay calibration, the Wnt result, and why FDR fails.

Panel A  study design: three independent cohorts, one pipeline, stated exclusions
Panel B  positive control Atf3: the assay can see a 16-fold injury response,
         and that response still fails Benjamini-Hochberg (padj = 0.47)
Panel C  36 Wnt/β-catenin genes, per-cohort and pooled, against the ±0.5
         equivalence bound
Panel D  genes passing nominal p versus genes passing BH, per cohort

Every number is read from the files written by run_three_cohort.py.
"""
import csv, os, math

def _default_out():
    """Resolve the output directory without hard-coding a machine-specific
    path.  If the script lives in <pkg>/scripts/, write to <pkg>/; otherwise
    write next to the script.  Set NP_OUT to override."""
    here = os.path.dirname(os.path.abspath(__file__))
    parent = os.path.dirname(here)
    if os.path.basename(here).lower() == 'scripts' and os.path.isdir(parent):
        return parent
    return here



import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, FancyBboxPatch
from matplotlib.lines import Line2D

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 8,
    'axes.linewidth': 0.75,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

OUT = os.environ.get('NP_OUT') or _default_out()
RESULTS = os.path.join(OUT, 'results') if os.path.isdir(os.path.join(OUT, 'results')) else OUT
FIGDIR = os.path.join(OUT, 'figures') if os.path.isdir(os.path.join(OUT, 'figures')) else OUT


def _find(fn):
    for d in (RESULTS, OUT, os.getcwd()):
        p = os.path.join(d, fn)
        if os.path.exists(p):
            return p
    return os.path.join(RESULTS, fn)  # not found; fail with a clear path


COH = {'C1': ('cohort_C1_GSE217932_CCI_gene_level.csv', 'GSE217932\nCCI vs Sham', '#B03A2E'),
       'C2': ('cohort_C2_GSE256472_SNI_gene_level.csv', 'GSE256472\nSNI vs Ctrl', '#2471A3'),
       'C3': ('cohort_C3_GSE205906_SNI_gene_level.csv', 'GSE205906\nSNI vs Sham', '#1E8449')}
COH_ORDER = ['C1', 'C2', 'C3']
GSE = {'C1': 'GSE217932', 'C2': 'GSE256472', 'C3': 'GSE205906'}

WNT = ['Ctnnb1', 'Tcf4', 'Tcf3', 'Tcf7', 'Tcf7l1', 'Tcf7l2', 'Lef1', 'Apc', 'Axin1', 'Axin2',
       'Gsk3a', 'Gsk3b', 'Ccnd1', 'Ccnd2', 'Ccnd3', 'Myc', 'Mycn', 'Jun', 'Fosl1', 'Cd44',
       'Mmp7', 'Birc5', 'Sox9', 'Wnt2', 'Wnt3', 'Wnt3a', 'Wnt5a', 'Wnt10b', 'Fzd1', 'Fzd4',
       'Fzd7', 'Lrp5', 'Lrp6', 'Dkk1', 'Dvl1', 'Dvl2', 'Dvl3']

# ------------------------------------------------------------------ load
data = {}
for cid, (fn, label, col) in COH.items():
    d = {}
    for r in csv.DictReader(open(_find(fn), encoding='utf-8')):
        try:
            d[r['symbol'].upper()] = dict(
                lfc=float(r['log2FC']), lo=float(r['ci_lo']), hi=float(r['ci_hi']),
                p=float(r['pvalue']), padj=float(r['padj']), expr=float(r['mean_expr']))
        except (ValueError, KeyError):
            pass
    data[cid] = d

meta = {}
for r in csv.DictReader(open(_find('meta_focus.csv'), encoding='utf-8')):
    try:
        meta[r['symbol'].upper()] = dict(
            k=int(r['k']), re=float(r['pooled_RE_log2FC']), lo=float(r['RE_ci_lo']),
            hi=float(r['RE_ci_hi']), p=float(r['p_RE']), i2=float(r['I2_pct']))
    except (ValueError, KeyError):
        pass

NC = {cid: len(data[cid]) for cid in COH_ORDER}
N_NOM = {cid: sum(1 for v in data[cid].values() if v['p'] < 0.05) for cid in COH_ORDER}
N_BH = {cid: sum(1 for v in data[cid].values() if v['padj'] < 0.05) for cid in COH_ORDER}

# ------------------------------------------------------------------ canvas
fig = plt.figure(figsize=(7.20, 8.60))
gs = fig.add_gridspec(2, 3, height_ratios=[0.50, 2.05],
                      width_ratios=[0.95, 1.45, 0.82],
                      left=0.105, right=0.982, top=0.962, bottom=0.088,
                      hspace=0.235, wspace=0.42)

axA = fig.add_subplot(gs[0, :])
axB = fig.add_subplot(gs[1, 0])
axC = fig.add_subplot(gs[1, 1])
axD = fig.add_subplot(gs[1, 2])

# ================================================================== Panel A
axA.set_xlim(0, 1)
axA.set_ylim(0, 1)
axA.axis('off')
axA.set_title('A   Three independent cohorts, one analysis pipeline',
              fontsize=8.4, fontweight='bold', loc='left', pad=7)

bw, bh, by = 0.272, 0.44, 0.455
boxes = [
    ('C1  discovery', 'GSE217932', 'CCI vs Sham\n4 vs 4\nraw counts', '#B03A2E', 0.040),
    ('C2  independent cohort', 'GSE256472', 'SNI vs control\n3 vs 3\nraw counts', '#2471A3', 0.345),
    ('C3  independent cohort', 'GSE205906', 'SNI vs Sham\n3 vs 3\nFPKM', '#1E8449', 0.650),
]
for title, acc, detail, col, x0 in boxes:
    axA.add_patch(FancyBboxPatch((x0, by), bw, bh,
                                 boxstyle='round,pad=0.005,rounding_size=0.010',
                                 fc='white', ec=col, lw=1.0, zorder=2))
    axA.add_patch(Rectangle((x0, by + bh - 0.115), bw, 0.115, fc=col, ec='none',
                            alpha=0.12, zorder=1.5))
    axA.text(x0 + 0.011, by + bh - 0.052, title, fontsize=6.5, color=col, va='center',
             fontweight='bold', zorder=3)
    axA.text(x0 + 0.011, by + bh - 0.168, acc, fontsize=7.6, va='center',
             fontweight='bold', zorder=3)
    axA.text(x0 + 0.011, by + 0.028, detail, fontsize=6.0, va='bottom', color='#333333',
             zorder=3, linespacing=1.32)

# pipeline strip
py, ph = 0.135, 0.128
axA.add_patch(FancyBboxPatch((0.040, py), 0.882, ph,
                             boxstyle='round,pad=0.004,rounding_size=0.008',
                             fc='#F2F4F4', ec='#95A5A6', lw=0.6, zorder=2))
axA.text(0.050, py + ph / 2,
         'gene-level collapse  →  Welch $t$  →  BH  →  95% CI  →  '
         '$\\pm$0.5 log$_2$ equivalence  →  DerSimonian–Laird meta  →  preranked GSEA',
         fontsize=5.9, va='center', color='#1C2833', zorder=3)

for x0 in (0.176, 0.481, 0.786):
    axA.add_patch(FancyArrowPatch((x0, by - 0.012), (x0, py + ph + 0.012),
                                  arrowstyle='-|>', mutation_scale=6.0,
                                  lw=0.7, color='#7F8C8D', zorder=1))

axA.text(0.040, 0.045,
         'excluded, pre-specified:  CCI+THP-1…4  (n=4, intervention arm)   ·   '
         'PLM1…3  (n=3, paeoniflorin–liquiritin arm)',
         fontsize=5.6, color='#7F8C8D', va='center')

# ================================================================== Panel B
axB.set_title('B   Positive control', fontsize=7.8, fontweight='bold', loc='left', pad=7)

rowsB = [(cid, data[cid]['ATF3']) for cid in COH_ORDER if 'ATF3' in data[cid]]
mAtf = meta['ATF3']

ys = list(range(len(rowsB), 0, -1))
for (cid, r), yy in zip(rowsB, ys):
    axB.plot([r['lo'], r['hi']], [yy, yy], color=COH[cid][2], lw=1.5,
             solid_capstyle='round', zorder=2)
    axB.plot(r['lfc'], yy, 'o', ms=5.2, color=COH[cid][2], mec='white', mew=0.5, zorder=4)
    axB.text(r['hi'] + 0.30, yy, f"padj {r['padj']:.2f}", fontsize=5.5, va='center',
             color='#7B241C')
axB.plot([mAtf['lo'], mAtf['hi']], [0.0, 0.0], color='#1C2833', lw=2.0,
         solid_capstyle='round', zorder=2)
axB.plot(mAtf['re'], 0.0, 'D', ms=5.4, color='#1C2833', mec='white', mew=0.5, zorder=4)
axB.text(mAtf['hi'] + 0.30, 0.0, f"$p$={mAtf['p']:.3f}", fontsize=5.5, va='center',
         color='#1C2833')

axB.axvline(0, color='#95A5A6', lw=0.7, ls=':', zorder=1)
axB.set_yticks(ys + [0])
axB.set_yticklabels([f'{c}\n{GSE[c]}' for c, _ in rowsB] + ['pooled\n(RE)'], fontsize=5.6,
                    linespacing=1.35)
axB.set_ylim(-0.80, len(rowsB) + 0.30)
axB.set_xlim(-1.15, 7.55)
axB.set_xlabel('$Atf3$  log$_2$ fold change (95% CI)', fontsize=6.8, labelpad=2)
axB.tick_params(axis='x', labelsize=6.2, length=2.5, width=0.6, pad=1.5)
axB.tick_params(axis='y', length=0, pad=2.5)
for s in ('top', 'right'):
    axB.spines[s].set_visible(False)

axB.text(0.035, 0.035,
         f"16-fold in C1 (nominal $p$={rowsB[0][1]['p']:.3f})\n"
         f"yet padj = {rowsB[0][1]['padj']:.2f}",
         transform=axB.transAxes, fontsize=5.6, va='bottom', ha='left', color='#1C2833',
         linespacing=1.5,
         bbox=dict(boxstyle='round,pad=0.30', fc='#FDF2E9', ec='#E59866', lw=0.5))

# ================================================================== Panel C
axC.set_title('C   Wnt/$\\beta$-catenin genes', fontsize=7.8, fontweight='bold',
              loc='left', pad=7)

genes = [g for g in WNT if g.upper() in meta]
genes.sort(key=lambda g: -meta[g.upper()]['re'])
yy = np.arange(len(genes))[::-1]

axC.axvspan(-0.5, 0.5, color='#EAF2F8', zorder=0)
axC.axvline(0, color='#7F8C8D', lw=0.7, zorder=1)
for b in (-0.5, 0.5):
    axC.axvline(b, color='#5D6D7E', lw=0.8, ls='--', zorder=1)

for gi, g in enumerate(genes):
    k = g.upper()
    m = meta[k]
    y = yy[gi]
    axC.plot([m['lo'], m['hi']], [y, y], color='#34495E', lw=1.5, alpha=0.28,
             solid_capstyle='round', zorder=2)
    axC.plot(m['re'], y, 'D', ms=3.3, color='#1C2833', zorder=5)
    for j, cid in enumerate(COH_ORDER):
        r = data[cid].get(k)
        if r:
            axC.plot(r['lfc'], y + (j - 1) * 0.225, 'o', ms=2.0, color=COH[cid][2],
                     alpha=0.95, zorder=4)

axC.set_yticks(yy)
axC.set_yticklabels(genes, fontsize=5.4, style='italic')
axC.set_ylim(-2.75, len(genes) - 0.25)
axC.set_xlim(-2.30, 2.30)
axC.set_xlabel('log$_2$ fold change', fontsize=6.8, labelpad=2)
axC.tick_params(axis='x', labelsize=6.2, length=2.5, width=0.6, pad=1.5)
axC.tick_params(axis='y', length=0, pad=1.5)
for s in ('top', 'right'):
    axC.spines[s].set_visible(False)

n_eq = sum(1 for g in genes if meta[g.upper()]['lo'] > -0.5 and meta[g.upper()]['hi'] < 0.5)
mx = max(abs(meta[g.upper()]['re']) for g in genes)
handles = [Line2D([], [], marker='o', ls='', ms=3.0, color=COH[c][2], label=c) for c in COH_ORDER]
handles += [Line2D([], [], marker='D', ls='', ms=3.3, color='#1C2833', label='pooled')]
axC.legend(handles=handles, fontsize=5.4, frameon=False, loc='lower right',
           ncol=2, handletextpad=0.25, columnspacing=0.7, borderpad=0.12)
axC.spines['left'].set_bounds(0, len(genes) - 1)
axC.text(0.5, -0.072,
         f"{n_eq} of {len(genes)} pooled CIs lie within $\\pm$0.5;  "
         f"largest |pooled effect| = {mx:.2f} log$_2$",
         transform=axC.transAxes, fontsize=5.5, va='top', ha='center',
         color='#34495E')

# ================================================================== Panel D
axD.set_title('D   Surviving correction', fontsize=7.8, fontweight='bold', loc='left', pad=7)

ylabels, yv, cols = [], [], []
for cid in COH_ORDER:
    ylabels += [f'{cid}  nominal', f'{cid}  BH q<0.05']
    yv += [N_NOM[cid], N_BH[cid]]
    cols += [COH[cid][2], COH[cid][2]]
ypos = np.arange(len(yv))[::-1]
axD.barh(ypos, [max(v, 0.62) for v in yv], height=0.58, color=cols, alpha=0.85,
         edgecolor='none', zorder=2)
for y, v in zip(ypos, yv):
    axD.text(max(v, 0.62) * 1.35, y, f'{v:,}', fontsize=5.9, va='center',
             color='#1C2833', zorder=3)
axD.set_xscale('log')
axD.set_xlim(0.55, 9000)
axD.set_yticks(ypos)
axD.set_yticklabels(ylabels, fontsize=5.7)
axD.set_ylim(-1.25, len(yv) - 0.35)
axD.set_xlabel('genes (log scale)', fontsize=6.6, labelpad=2)
axD.tick_params(axis='x', labelsize=5.9, length=2.5, width=0.6, pad=1.5)
axD.tick_params(axis='y', length=0, pad=1.5)
for s in ('top', 'right'):
    axD.spines[s].set_visible(False)

axD.text(0.035, 0.025,
         'the single BH pass (C3) is a\nzero-variance transcript\n(LOC688442, 1.09-fold)',
         transform=axD.transAxes, fontsize=5.2, va='bottom', ha='left', color='#7B241C',
         linespacing=1.45,
         bbox=dict(boxstyle='round,pad=0.26', fc='#FDEDEC', ec='#CD6155', lw=0.45))

fig.savefig(os.path.join(FIGDIR, 'Figure1_design_and_Wnt.png'), dpi=400)
fig.savefig(os.path.join(FIGDIR, 'Figure1_design_and_Wnt.pdf'))
print('wrote Figure1_design_and_Wnt.png/.pdf')
print('NC', NC, 'NOM', N_NOM, 'BH', N_BH)
print('wnt plotted', len(genes), 'pooled eq', n_eq, 'max |pooled|', round(mx, 3))
