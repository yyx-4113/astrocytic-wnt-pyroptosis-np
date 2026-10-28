# -*- coding: utf-8 -*-
"""
A-task (Figure 4 rebuild): regenerate Figure 4 from the FULL-transcriptome
GSE134003 data, reflecting the corrected scRNA conclusion (B-task):
  - pyroptosis signal is in microglia/macrophages, NOT astrocytes
  - astrocyte pyroptosis is low-detection and unchanged SNI vs SI
No pandas; numpy + matplotlib only.
"""
import gzip, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

META = r'D:\省自然\生信分析\data\scRNAseq\GSE134003_ClusterMetadata.txt.gz'
EXPR = r'D:\省自然\生信分析\data\scRNAseq\GSE134003_SI_SNI_normalizedcounts.txt.gz'
OUT  = r'D:\省自然\生信分析\results\Figure4_scRNAseq_cell_atlas.png'

def map_ct(s):
    s = str(s)
    if 'Oligo' in s and 'myl' not in s: return 'Oligodendrocytes'
    elif 'mylOligo' in s: return 'Myelinating_Oligo'
    elif 'Micro' in s: return 'Microglia'
    elif 'Astrocyte' in s: return 'Astrocytes'
    elif 'Endothelial' in s: return 'Endothelial'
    elif 'Pericyte' in s: return 'Pericytes'
    elif 'smVasc' in s: return 'Vascular_SMC'
    elif 'Leptomeningeal' in s: return 'Leptomeningeal'
    elif 'OPC' in s: return 'OPCs'
    elif 'Excit' in s: return 'Excitatory_Neurons'
    elif 'Inhib' in s: return 'Inhibitory_Neurons'
    elif 'Schwann' in s and 'myl' in s: return 'Myelinating_Schwann'
    elif 'Schwann' in s and 'nm' in s: return 'NM_Schwann'
    elif 'Macro' in s: return 'Macrophages'
    else: return 'Other'

print('load meta', flush=True)
bc_info = {}
samples = {}
with gzip.open(META, 'rt') as f:
    f.readline()
    for line in f:
        p = line.rstrip('\n').split('\t')
        bc, sample, cond, subtype = p[0], p[3], p[4], p[6]
        bc_info[bc] = (sample, cond, map_ct(subtype))
        samples[sample] = cond

pyro = ['Nlrp3','Nlrc4','Aim2','Pycard','Casp1','Casp4','Gsdmd','Il1b','Il18',
        'Hmgb1','Tnf','Il6','Nek7','Txnip','Tlr4','Nfkb1','Rela','Casp8','Casp3','Gzmb']
wnt = ['Wnt1','Wnt3','Wnt3a','Wnt4','Wnt5a','Wnt7a','Wnt10b','Fzd1','Fzd2','Fzd7',
       'Lrp5','Lrp6','Ctnnb1','Tcf4','Tcf7','Tcf7l2','Lef1','Apc','Axin1','Gsk3b','Myc']
gs = set(pyro + wnt)
expr = {}
print('extract genes', flush=True)
with gzip.open(EXPR, 'rt') as f:
    header = f.readline().rstrip('\n').split('\t')
    cells = header[1:]
    for line in f:
        parts = line.rstrip('\n').split('\t')
        if parts[0] in gs:
            expr[parts[0]] = np.array(parts[1:], dtype=np.float32)
            gs.discard(parts[0])
            if not gs: break

ct_of_cell = np.array([bc_info.get(c, ('?','?','Other'))[2] for c in cells])
cond_of_cell = np.array([bc_info.get(c, ('?','?','Other'))[1] for c in cells])
animal_of_cell = np.array([bc_info.get(c, ('?','?','Other'))[0] for c in cells])

ct_order = ['Astrocytes','Microglia','Macrophages','Endothelial','OPCs',
            'Oligodendrocytes','Excitatory_Neurons','Inhibitory_Neurons','Other']
ct_order = [c for c in ct_order if c in set(ct_of_cell)]

def mod_score(genes, mask):
    rows = [expr[g][mask] for g in genes if g in expr]
    if not rows: return np.array([])
    return np.vstack(rows).mean(axis=0)

# ---- 4A composition ----
comp = {}
for cond in ['SI','SNI']:
    m = cond_of_cell == cond
    cnt = np.array([(ct_of_cell[m]==c).sum() for c in ct_order], dtype=float)
    comp[cond] = cnt / cnt.sum() * 100

# ---- 4B/4C dot matrices (mean in expressed, pct>0) ----
def dot_matrix(genelist):
    mean_m = np.zeros((len(genelist), len(ct_order)))
    pct_m = np.zeros_like(mean_m)
    for j, ct in enumerate(ct_order):
        m = ct_of_cell == ct
        for i, g in enumerate(genelist):
            if g in expr:
                v = expr[g][m]
                pos = v[v > 0]
                mean_m[i, j] = pos.mean() if len(pos) else 0
                pct_m[i, j] = 100 * (v > 0).mean()
    return mean_m, pct_m

wnt_mean, _ = dot_matrix(wnt)
pyro_mean, _ = dot_matrix(pyro)

# ---- 4D pyro module score per cell type, animal-level mean +/- SEM ----
pyro_means_per_animal = {ct: [] for ct in ct_order}
for ct in ct_order:
    m = ct_of_cell == ct
    sc = mod_score(pyro, m)
    if len(sc)==0: continue
    ani = animal_of_cell[m]
    for a in np.unique(ani):
        pyro_means_per_animal[ct].append(sc[ani==a].mean())
pyro_bar = np.array([np.mean(pyro_means_per_animal[c]) if pyro_means_per_animal[c] else 0 for c in ct_order])
pyro_sem = np.array([np.std(pyro_means_per_animal[c], ddof=1)/np.sqrt(len(pyro_means_per_animal[c])) if len(pyro_means_per_animal[c])>1 else 0 for c in ct_order])

print('plotting', flush=True)
fig, axes = plt.subplots(2, 2, figsize=(15, 11))
# 4A
ax = axes[0,0]
x = np.arange(len(ct_order)); w=0.38
ax.bar(x-w/2, comp['SI'], w, label='SI', color='#3498db')
ax.bar(x+w/2, comp['SNI'], w, label='SNI', color='#e74c3c')
ax.set_xticks(x); ax.set_xticklabels(ct_order, rotation=45, ha='right', fontsize=8)
ax.set_ylabel('Cell proportion (%)'); ax.set_title('A. Cell-type composition (SI vs SNI)'); ax.legend(fontsize=8)

# 4B Wnt heatmap
ax = axes[0,1]
im = ax.imshow(wnt_mean, aspect='auto', cmap='YlOrRd')
ax.set_yticks(range(len(wnt))); ax.set_yticklabels(wnt, fontsize=6)
ax.set_xticks(range(len(ct_order))); ax.set_xticklabels(ct_order, rotation=45, ha='right', fontsize=7)
ax.set_title('B. Wnt signaling expression (mean in expressed)')
fig.colorbar(im, ax=ax, shrink=0.7)

# 4C Pyro heatmap
ax = axes[1,0]
im2 = ax.imshow(pyro_mean, aspect='auto', cmap='YlOrRd')
ax.set_yticks(range(len(pyro))); ax.set_yticklabels(pyro, fontsize=6)
ax.set_xticks(range(len(ct_order))); ax.set_xticklabels(ct_order, rotation=45, ha='right', fontsize=7)
ax.set_title('C. Pyroptosis expression (enriched in Microglia/Macrophages)')
fig.colorbar(im2, ax=ax, shrink=0.7)

# 4D pyro module score bar
ax = axes[1,1]
colors = ['#e74c3c' if c=='Astrocytes' else '#3498db' for c in ct_order]
ax.bar(range(len(ct_order)), pyro_bar, yerr=pyro_sem, color=colors, alpha=0.85, capsize=3)
ax.set_xticks(range(len(ct_order))); ax.set_xticklabels(ct_order, rotation=45, ha='right', fontsize=7)
ax.axhline(0, color='k', lw=0.5)
ax.set_ylabel('Pyroptosis module score (animal-level mean ± SEM)')
ax.set_title('D. Cell-type pyroptosis (Microglia/Macrophage >> Astrocyte)')

fig.suptitle('Figure 4: GSE134003 SNI spinal cord, full-transcriptome (14,700 genes; n=19 animals)',
             fontsize=12, fontweight='bold')
fig.tight_layout(rect=[0,0,1,0.97])
fig.savefig(OUT, dpi=160, bbox_inches='tight')
print('saved', OUT, flush=True)
