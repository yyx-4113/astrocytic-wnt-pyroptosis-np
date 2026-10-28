# -*- coding: utf-8 -*-
"""A-task: rebuild Figure 6 as SINGLE-cohort CIBERSORTx deconvolution (GSE217932 only).
Removes the fabricated GSE186237 panel and the false cross-dataset r=0.89 claim.
Uses the REAL deconvolution values from Figure6_deconv_GSE217932.csv
(astrocyte fraction ~0.10, NOT the 38.2-62.1% figure quoted in older text)."""
import csv, numpy as np
from scipy import stats
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

SRC = r'D:\省自然\生信分析\results\Figure6_deconv_GSE217932.csv'
OUT = r'D:\省自然\生信分析\results\Figure6_Deconvolution.png'

with open(SRC) as f:
    r = csv.reader(f)
    header = next(r)
    data = [row for row in r if row]
cts = header[1:]
sham = np.array([[float(x) for x in row[1:]] for row in data if row[0].startswith('Sham')])
cci  = np.array([[float(x) for x in row[1:]] for row in data if row[0].startswith('CCI')])
sham_m = sham.mean(0); cci_m = cci.mean(0)

# astrocyte sham vs CCI test
astro_i = cts.index('Astrocytes')
_, p_astro = stats.mannwhitneyu(sham[:, astro_i], cci[:, astro_i])

fig, ax = plt.subplots(figsize=(13, 5))
x = np.arange(len(cts)); w = 0.38
ax.bar(x - w/2, sham_m*100, w, label=f'Sham (n={len(sham)})', color='#3498db')
ax.bar(x + w/2, cci_m*100, w, label=f'CCI (n={len(cci)})', color='#e74c3c')
ax.set_xticks(x); ax.set_xticklabels(cts, rotation=45, ha='right', fontsize=7)
ax.set_ylabel('Estimated cell proportion (%)')
ax.set_title(f'Figure 6. CIBERSORTx deconvolution of GSE217932 (single cohort)\n'
             f'Astrocyte: Sham {sham_m[astro_i]*100:.1f}% vs CCI {cci_m[astro_i]*100:.1f}% '
             f'(Mann-Whitney p={p_astro:.3f}, not significant). No cross-dataset validation.')
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUT, dpi=160, bbox_inches='tight')
print('saved', OUT, '| astro Sham%%=%.2f CCI%%=%.2f p=%.3f' %
      (sham_m[astro_i]*100, cci_m[astro_i]*100, p_astro), flush=True)
