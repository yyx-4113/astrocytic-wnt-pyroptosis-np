# -*- coding: utf-8 -*-
"""A-task: rebuild Figure 1 as a SINGLE-cohort volcano (GSE217932 only).
Removes the fabricated GSE186237 volcano panel."""
import csv, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

SRC = r'D:\省自然\生信分析\results\GSE217932_DEGs.csv'
OUT = r'D:\省自然\生信分析\results\Figure1_DEG_overview.png'

rows = list(csv.DictReader(open(SRC)))
lg2 = np.array([float(r['log2FC']) for r in rows])
p   = np.array([float(r['pvalue']) for r in rows])
sym = [r['symbol'] for r in rows]
logp = -np.log10(np.clip(p, 1e-300, 1.0))

sig = (p < 0.05) & (np.abs(lg2) > 1)
up  = sig & (lg2 > 0)
dn  = sig & (lg2 < 0)

fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(lg2[~sig], logp[~sig], s=6, c='#bdc3c7', alpha=0.4, label='ns')
ax.scatter(lg2[up],  logp[up],  s=10, c='#e74c3c', alpha=0.7, label=f'up (n={up.sum()})')
ax.scatter(lg2[dn],  logp[dn],  s=10, c='#3498db', alpha=0.7, label=f'down (n={dn.sum()})')
ax.axhline(-np.log10(0.05), ls='--', c='k', lw=0.8)
ax.axvline(-1, ls=':', c='gray'); ax.axvline(1, ls=':', c='gray')

# label top genes by |log2FC|
order = np.argsort(-np.abs(lg2))[:8]
for i in order:
    if sym[i]:
        ax.annotate(sym[i], (lg2[i], logp[i]), fontsize=6, xytext=(2,2),
                    textcoords='offset points')

ax.set_xlabel('log2 fold change (Sham vs CCI)')
ax.set_ylabel('-log10 nominal p')
ax.set_title(f'Figure 1. GSE217932 bulk RNA-seq DEG volcano (single cohort)\n'
             f'n_genes={len(rows)}; nominal p<0.05 & |log2FC|>1: up={up.sum()}, down={dn.sum()}\n'
             f'(No second/independent bulk cohort; cross-dataset replication NOT claimed)')
ax.legend(fontsize=8, loc='upper right')
fig.tight_layout()
fig.savefig(OUT, dpi=160, bbox_inches='tight')
print('saved', OUT, '| sig-up', int(up.sum()), 'sig-down', int(dn.sum()), flush=True)
