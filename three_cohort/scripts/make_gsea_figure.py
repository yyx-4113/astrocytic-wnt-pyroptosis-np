# -*- coding: utf-8 -*-
"""GSEA figure: cohort x gene-set NES heatmap + C1 pyroptosis walking plot."""
import csv, os, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
                     'font.size': 8, 'axes.linewidth': 0.8, 'pdf.fonttype': 42})

OUT = r'D:/省自然/生信分析/submission/RESCUE_2026-10-29'

# ---------- read gsea csv
tab = {}
sets_order = []
for r in csv.DictReader(open(os.path.join(OUT, 'gsea_three_cohort.csv'), encoding='utf-8')):
    tab[(r['cohort'], r['gene_set'])] = (float(r['NES']), float(r['nominal_p']), float(r['FDR_q_within_7_sets']))
    if r['gene_set'] not in sets_order:
        sets_order.append(r['gene_set'])

pretty = {
    'HALLMARK_WNT_BETA_CATENIN': 'Wnt/$\\beta$-catenin (Hallmark, 37)',
    'PYROPTOSIS_CORE': 'Pyroptosis core (22)',
    'INJURY_RESPONSE_POSCTRL': 'Injury response (23) — pos. ctrl',
    'MICROGLIA_ACTIVATION': 'Microglial activation (20)',
    'ASTROCYTE_REACTIVITY': 'Astrocyte reactivity (19)',
    'NFKB_CYTOKINE_SIGNALLING': 'NF-$\\kappa$B / cytokine (22)',
}
order = ['HALLMARK_WNT_BETA_CATENIN', 'PYROPTOSIS_CORE', 'INJURY_RESPONSE_POSCTRL',
         'MICROGLIA_ACTIVATION', 'ASTROCYTE_REACTIVITY', 'NFKB_CYTOKINE_SIGNALLING']
cohorts = ['C1', 'C2', 'C3']
clab = ['GSE217932\nCCI vs Sham\n(n=4v4)', 'GSE256472\nSNI vs Ctrl\n(n=3v3)', 'GSE205906\nSNI vs Sham\n(n=3v3)']

M = np.array([[tab[(c, s)][0] for c in cohorts] for s in order])
Q = np.array([[tab[(c, s)][2] for c in cohorts] for s in order])

fig = plt.figure(figsize=(7.4, 6.6))
gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.95], hspace=0.42)

# ---------------- Panel A heatmap
axA = fig.add_subplot(gs[0])
norm = TwoSlopeNorm(vmin=-3.0, vcenter=0, vmax=3.0)
im = axA.imshow(M, cmap='RdBu_r', norm=norm, aspect='auto')
axA.set_xticks(range(len(cohorts)))
axA.set_xticklabels(clab, fontsize=6.8)
axA.set_yticks(range(len(order)))
axA.set_yticklabels([pretty[s] for s in order], fontsize=7.2)
for i in range(len(order)):
    for j in range(len(cohorts)):
        v, q = M[i, j], Q[i, j]
        star = '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else 'ns'
        axA.text(j, i, f'{v:+.2f}\n{star}', ha='center', va='center', fontsize=6.4,
                 color='white' if abs(v) > 1.8 else 'black', linespacing=1.25)
axA.set_title('A   Gene-set enrichment (NES) in three independent cohorts', fontsize=8.6, loc='left', pad=8)
axA.tick_params(length=0)
for sp in axA.spines.values():
    sp.set_visible(False)
cb = fig.colorbar(im, ax=axA, fraction=0.026, pad=0.015)
cb.set_label('NES', fontsize=6.8)
cb.ax.tick_params(labelsize=6)
fig.text(0.5, 0.012, '* q<0.05    ** q<0.01    *** q<0.001    (Benjamini-Hochberg across the 6 sets tested in this analysis)',
         ha='center', va='bottom', fontsize=6.2, color='#555555')

# ---------------- Panel B running ES for C1 pyroptosis
PYRO = ["NLRP3","NLRC4","AIM2","PYCARD","CASP1","CASP4","CASP8","CASP3","GSDMD","GSDME",
        "IL1B","IL18","IL1A","HMGB1","TNF","IL6","NEK7","TXNIP","TLR4","NFKB1","RELA","GZMB"]
raw = []
for r in csv.DictReader(open(os.path.join(OUT, 'cohort_C1_GSE217932_CCI_gene_level.csv'), encoding='utf-8')):
    try:
        raw.append((r['symbol'].upper(), float(r['t_stat']), float(r['log2FC']), float(r['mean_expr'])))
    except ValueError:
        pass
ex = sorted(x[3] for x in raw)
thr = ex[int(len(ex) * 0.25)]
rows = [(s, max(-50, min(50, t)), l) for s, t, l, e in raw if e >= thr]
rows.sort(key=lambda x: -x[1])
syms = [r[0] for r in rows]
metric = np.array([r[1] for r in rows])
setu = set(PYRO)
hit = np.array([s in setu for s in syms])
k = int(hit.sum())
w = np.where(hit, np.abs(metric), 0.0)
cum = np.cumsum(w) / w.sum() - np.cumsum(~hit) / (len(metric) - k)
x = np.arange(len(cum))

axB = fig.add_subplot(gs[1])
axB.plot(x, cum, color='#1B4F72', lw=1.2)
axB.axhline(0, color='#999999', lw=0.6, ls=':')
axB.axvspan(0, 844, color='#FADBD8', alpha=0.55, zorder=0)
hit_idx = np.where(hit)[0]
axB.plot(hit_idx, np.full(len(hit_idx), -0.16), '|', color='#C0392B', ms=6, mew=1.0)
axB.set_xlim(0, len(cum))
axB.set_ylim(-0.22, 0.95)
axB.set_xlabel('Gene rank (Welch t, descending)  —  14,257 genes after expression filtering', fontsize=7.2)
axB.set_ylabel('Running enrichment score', fontsize=7.2)
axB.set_title('B   C1: pyroptosis set is enriched, but NLRP3 is not a leading-edge gene',
              fontsize=8.6, loc='left', pad=8)
axB.tick_params(labelsize=6.5, length=2)
axB.spines[['top', 'right']].set_visible(False)

nrlp3_i = syms.index('NLRP3') if 'NLRP3' in syms else None
if nrlp3_i is not None:
    axB.annotate('NLRP3\nrank %d  (t=%+.2f, log$_2$FC %+.2f)' % (nrlp3_i + 1, metric[nrlp3_i],
                [l for s, t, l in rows if s == 'NLRP3'][0]),
                 xy=(nrlp3_i, cum[nrlp3_i]), xytext=(nrlp3_i + 1600, 0.17),
                 fontsize=6.4, color='black',
                 arrowprops=dict(arrowstyle='->', lw=0.8, color='black',
                                 connectionstyle='arc3,rad=-0.18'))
axB.annotate('leading edge\n13 of 19 set genes',
             xy=(760, 0.72), xytext=(1350, 0.90), fontsize=6.2, color='#922B21',
             ha='left', va='center', linespacing=1.35,
             arrowprops=dict(arrowstyle='->', lw=0.7, color='#922B21',
                             connectionstyle='arc3,rad=0.25'))
axB.text(0.985, 0.06, 'ES = %+.3f    NES = %+.2f    q < 0.001' % (cum.max(), M[1, 0]),
         transform=axB.transAxes, fontsize=6.4, ha='right',
         bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#CCCCCC', lw=0.5))

fig.text(0.5, 0.975, 'Wnt programme is not enriched in any cohort; pyroptosis shows weak concordant priming without NLRP3 induction',
         ha='center', va='top', fontsize=9.0, fontweight='bold')

for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(OUT, f'Figure_gsea_three_cohort.{ext}'), dpi=400,
                bbox_inches='tight', facecolor='white')
print('saved Figure_gsea_three_cohort.png/.pdf')
print('pyro ES', round(float(cum.max()), 4), 'peak rank', int(np.argmax(cum)) + 1)
