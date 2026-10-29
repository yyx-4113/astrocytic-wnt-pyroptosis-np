# -*- coding: utf-8 -*-
"""Figure: three-cohort forest plot + detectability calibration."""
import csv, math, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
import matplotlib.font_manager as fm

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 8,
    'axes.linewidth': 0.8,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

OUT = r'D:/省自然/生信分析/submission/RESCUE_2026-10-29'
COH = {'C1': ('cohort_C1_GSE217932_CCI_gene_level.csv', 'GSE217932  CCI vs Sham', '#C0392B'),
       'C2': ('cohort_C2_GSE256472_SNI_gene_level.csv', 'GSE256472  SNI vs Ctrl', '#2471A3'),
       'C3': ('cohort_C3_GSE205906_SNI_gene_level.csv', 'GSE205906  SNI vs Sham', '#1E8449')}

data = {}
for cid, (fn, label, col) in COH.items():
    d = {}
    for r in csv.DictReader(open(os.path.join(OUT, fn), encoding='utf-8')):
        try:
            d[r['symbol'].upper()] = dict(lfc=float(r['log2FC']), lo=float(r['ci_lo']),
                                          hi=float(r['ci_hi']), p=float(r['pvalue']),
                                          padj=float(r['padj']), expr=float(r['mean_expr']))
        except (ValueError, KeyError):
            pass
    data[cid] = d

meta = {}
for r in csv.DictReader(open(os.path.join(OUT, 'meta_focus.csv'), encoding='utf-8')):
    try:
        meta[r['symbol'].upper()] = dict(k=int(r['k']), re=float(r['pooled_RE_log2FC']),
                                        lo=float(r['RE_ci_lo']), hi=float(r['RE_ci_hi']),
                                        p=float(r['p_RE']), i2=float(r['I2_pct']))
    except (ValueError, KeyError):
        pass


def forest(ax, genes, title, note=None, xlim=(-1.6, 1.6), hide_ylabel=False):
    y = 0
    ticks, labels = [], []
    shades = ['#FDFEFE', '#F4F6F7']
    for gi, g in enumerate(genes):
        k = g.upper()
        base = y
        for j, cid in enumerate(['C1', 'C2', 'C3']):
            r = data[cid].get(k)
            yy = base - j * 0.22
            if r:
                ax.plot([r['lo'], r['hi']], [yy, yy], color=COH[cid][2], lw=1.0, solid_capstyle='butt')
                ax.plot(r['lfc'], yy, 'o', ms=3.0, color=COH[cid][2], zorder=3)
            else:
                ax.text(0.0, yy, 'n/a', fontsize=5, ha='center', va='center', color='#999999')
        m = meta.get(k)
        ym = base - 3 * 0.22 - 0.10
        if m:
            ax.plot([m['lo'], m['hi']], [ym, ym], color='black', lw=1.1)
            ax.plot(m['re'], ym, 'D', ms=3.6, color='black', zorder=4)
            star = '***' if m['p'] < 0.001 else '**' if m['p'] < 0.01 else '*' if m['p'] < 0.05 else 'ns'
            ax.text(xlim[1] - 0.02, ym, star, fontsize=6, ha='right', va='center',
                    color='black' if star != 'ns' else '#888888')
        ticks.append(base - 3 * 0.22 - 0.10)
        labels.append(g)
        y = base - (3 * 0.22 + 0.10) - 0.30
    ax.axvline(0, color='black', lw=0.7, zorder=1)
    ax.axvspan(-0.5, 0.5, color='#D5D8DC', alpha=0.35, zorder=0)
    if hide_ylabel:
        ax.set_yticks([])
    else:
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels, fontsize=7.5, style='italic')
    ax.set_xlim(*xlim)
    ax.set_xlabel('log$_2$ fold change (case vs control)', fontsize=7.5)
    ax.set_title(title, fontsize=8.5, pad=6, loc='left')
    ax.tick_params(axis='x', labelsize=6.5, length=2)
    ax.tick_params(axis='y', length=0)
    ax.spines[['top', 'right']].set_visible(False)
    if note:
        ax.text(0.99, 0.015, note, transform=ax.transAxes, fontsize=6,
                ha='right', va='bottom', color='#555555')


fig = plt.figure(figsize=(7.2, 8.4))
gs = fig.add_gridspec(2, 2, height_ratios=[0.85, 1.5], hspace=0.34, wspace=0.30)

axA = fig.add_subplot(gs[0, 0])
forest(axA, ['Atf3'], 'A   Atf3 — positive control (injury-response gene)',
       xlim=(-1.5, 6.6), hide_ylabel=True)

axD = fig.add_subplot(gs[0, 1])
# calibration: CI width vs mean expression, C1
c1 = data['C1']
xs, ys, eqs = [], [], []
for k, r in c1.items():
    w = r['hi'] - r['lo']
    if not math.isfinite(w) or w > 8:
        continue
    xs.append(r['expr']); ys.append(w)
    eqs.append(1 if (r['lo'] > -0.5 and r['hi'] < 0.5) else 0)
xs = np.array(xs); ys = np.array(ys); eqs = np.array(eqs)
axD.scatter(xs[eqs == 0], ys[eqs == 0], s=2.2, c='#C8CDD2', alpha=0.45, lw=0, rasterized=True)
axD.scatter(xs[eqs == 1], ys[eqs == 1], s=2.2, c='#2C3E50', alpha=0.55, lw=0, rasterized=True)
axD.axhline(1.0, color='#C0392B', lw=0.9, ls='--', zorder=5)
axD.set_xlabel('mean log$_2$ expression', fontsize=7.5)
axD.set_ylabel('95% CI width of log$_2$FC', fontsize=7.5)
axD.set_title('B   Detectability is set by expression level', fontsize=8.5, pad=6, loc='left')
axD.set_ylim(0, 6)
axD.tick_params(labelsize=6.5, length=2)
axD.spines[['top', 'right']].set_visible(False)
_eq_all = sum(1 for _r in c1.values() if _r['lo'] > -0.5 and _r['hi'] < 0.5)
_n_all = len(c1)
_skip = _n_all - len(eqs)
axD.text(0.025, 0.965,
         f'GSE217932 (n=4 vs 4)\n{_eq_all:,} of {_n_all:,} genes ({100*_eq_all/_n_all:.0f}%) have a\n'
         f'95% CI entirely within $\\pm$0.5\n({_skip} genes with CI width >8 not drawn)',
         transform=axD.transAxes, fontsize=5.8, ha='left', va='top', color='#2C3E50',
         linespacing=1.45,
         bbox=dict(boxstyle='round,pad=0.35', fc='white', ec='none', alpha=0.92))


axB = fig.add_subplot(gs[1, 0])
forest(axB, ['Ctnnb1', 'Axin2', 'Ccnd1', 'Ccnd2', 'Tcf4', 'Lef1', 'Gsk3b', 'Apc', 'Jun'],
       'C   Wnt/$\\beta$-catenin pathway (37 genes tested)', xlim=(-2.5, 1.2))

axC = fig.add_subplot(gs[1, 1])
forest(axC, ['Nlrp3', 'Pycard', 'Casp1', 'Gsdmd', 'Il1b', 'Il6', 'Il1a', 'Tnf'],
       'D   Pyroptosis machinery (22 genes tested)', xlim=(-1.9, 2.1))

handles = [plt.Line2D([], [], color=COH[c][2], lw=1.6, label=COH[c][1]) for c in ['C1', 'C2', 'C3']]
handles += [plt.Line2D([], [], color='black', marker='D', lw=1.2, ms=4, label='Random-effects pooled (k=2-3)')]
fig.legend(handles=handles, loc='lower center', ncol=2, fontsize=6.8, frameon=False,
           bbox_to_anchor=(0.5, 0.002))
fig.text(0.5, 0.985, 'Wnt and pyroptosis transcriptional programmes are not induced in neuropathic-pain spinal cord',
         ha='center', va='top', fontsize=9.2, fontweight='bold')
fig.text(0.5, 0.958, 'Three independent rat bulk RNA-seq cohorts, identical pipeline; grey band = $\\pm$0.5 log$_2$ (1.41-fold) equivalence margin',
         ha='center', va='top', fontsize=6.6, color='#555555')

for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(OUT, f'Figure_three_cohort_validation.{ext}'), dpi=400,
                bbox_inches='tight', facecolor='white')
print('saved Figure_three_cohort_validation.png/.pdf')
print(f'panel B: equivalence {int((eqs==1).sum())}/{len(eqs)} in C1')
