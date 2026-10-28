# -*- coding: utf-8 -*-
"""
B-task: scRNA reconciliation for GSE134003 (neuropathic-pain spinal cord).
Goal (Reviewer 1 #5/#6): re-evaluate whether astrocyte-enriched pyroptosis
expression is supported when using the FULL transcriptome (NOT the 141-gene
pre-filtered matrix) and accounting for BETWEEN-ANIMAL variation.

Inputs (unrestricted):
  data/scRNAseq/GSE134003_SI_SNI_normalizedcounts.txt.gz  (gene x cell, log-norm)
  data/scRNAseq/GSE134003_ClusterMetadata.txt.gz          (cell -> Subtype,Condition,Sample)
No pandas; pure numpy + gzip streaming.
"""
import gzip, numpy as np
from scipy import stats

META = r'D:\省自然\生信分析\data\scRNAseq\GSE134003_ClusterMetadata.txt.gz'
EXPR = r'D:\省自然\生信分析\data\scRNAseq\GSE134003_SI_SNI_normalizedcounts.txt.gz'

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

print('=== STEP 1: read metadata ===', flush=True)
bc_info = {}   # barcode -> (sample, condition, celltype)
samples = {}
with gzip.open(META, 'rt') as f:
    header = f.readline().rstrip('\n').split('\t')
    # header: ['', Num.Genes, Num.Transcripts, Sample, Condition, InitialCluster, Subtype]
    for line in f:
        p = line.rstrip('\n').split('\t')
        bc, sample, cond, subtype = p[0], p[3], p[4], p[6]
        ct = map_ct(subtype)
        bc_info[bc] = (sample, cond, ct)
        samples.setdefault(sample, cond)
print('total cells in meta:', len(bc_info), flush=True)
print('unique samples:', len(samples), flush=True)
from collections import Counter
print('condition per sample:', dict(Counter(samples.values())), flush=True)

# CellType x Condition cell counts
ct_cond = Counter()
for bc, (s, c, ct) in bc_info.items():
    ct_cond[(ct, c)] += 1
print('\n--- CellType x Condition cell counts ---', flush=True)
for (ct, c) in sorted(ct_cond):
    print(f'  {ct:22s} {c:4s}: {ct_cond[(ct,c)]:5d}', flush=True)

# animals per (CellType, Condition)
animal_ct = {}
for bc, (s, c, ct) in bc_info.items():
    animal_ct.setdefault((ct, c), set()).add(s)
print('\n--- Animals per (CellType, Condition) ---', flush=True)
for (ct, c) in sorted(animal_ct):
    print(f'  {ct:22s} {c:4s}: {len(animal_ct[(ct,c)])} animals', flush=True)

print('\n=== STEP 2: extract gene sets from FULL transcriptome ===', flush=True)
pyro = ['Nlrp3','Nlrc4','Aim2','Pycard','Casp1','Casp4','Gsdmd','Il1b','Il18',
        'Hmgb1','Tnf','Il6','Nek7','Txnip','Tlr4','Nfkb1','Rela','Casp8','Casp3','Gzmb']
wnt_sig = ['Ctnnb1','Tcf4','Tcf3','Tcf7','Tcf7l1','Tcf7l2','Lef1','Apc','Axin1',
           'Axin2','Gsk3a','Gsk3b','Ccnd1','Ccnd2','Ccnd3','Myc','Mycn','Jun','Fosl1','Cd44','Mmp7','Birc5','Sox9']
astro_mk = ['Gfap','Aldh1l1','S100b','Aqp4','Slc1a3','Slc1a2']
micro_mk = ['Aif1','Cx3cr1','Tmem119','Ctss','Itgam','Fcgr1']
neu_mk = ['Snap25','Syp','Map2']

all_genes = set(pyro + wnt_sig + astro_mk + micro_mk + neu_mk)
gene_expr = {}
with gzip.open(EXPR, 'rt') as f:
    header = f.readline().rstrip('\n').split('\t')
    cells = header[1:]
    # align cells to bc_info
    cell_idx = np.arange(len(cells))
    for line in f:
        parts = line.rstrip('\n').split('\t')
        g = parts[0]
        if g in all_genes:
            gene_expr[g] = np.array(parts[1:], dtype=np.float32)
            all_genes.discard(g)
            if not all_genes:
                break
print('genes found:', sorted(set(pyro+wnt_sig+astro_mk+micro_mk+neu_mk) - all_genes), flush=True)
print('genes MISSING:', sorted(all_genes), flush=True)

# build per-cell arrays grouped by celltype+condition+animal
# index cells
ct_of_cell = np.array([bc_info.get(c, ('?','?','Other'))[2] for c in cells])
cond_of_cell = np.array([bc_info.get(c, ('?','?','Other'))[1] for c in cells])
animal_of_cell = np.array([bc_info.get(c, ('?','?','Other'))[0] for c in cells])

def module_score(genes, cell_mask):
    """mean expression across genes for cells in mask (nan if none expressed)."""
    sub = None
    for g in genes:
        if g in gene_expr:
            arr = gene_expr[g][cell_mask]
            sub = arr if sub is None else np.vstack([sub, arr])
    if sub is None: return np.array([])
    return sub.mean(axis=0)  # per-cell mean across genes

print('\n=== STEP 3: cell-type level pyroptosis expression (FULL transcriptome) ===', flush=True)
cts_of_interest = ['Astrocytes','Microglia','Oligodendrocytes','OPCs',
                   'Excitatory_Neurons','Inhibitory_Neurons','Macrophages','Endothelial']
print(f'{"CellType":22s} {"n":>6s} {"pyro_mean_all":>13s} {"pyro_mean_pos":>14s} {"pyro_pct>0":>10s}', flush=True)
ct_pyro_mean = {}
for ct in cts_of_interest:
    mask = ct_of_cell == ct
    n = mask.sum()
    if n == 0: 
        print(f'  {ct:22s} {0:6d} (none)', flush=True); continue
    sc = module_score(pyro, mask)
    mean_all = np.mean(sc)
    pos = sc[sc > 0]
    mean_pos = np.mean(pos) if len(pos) else 0.0
    pct = 100.0 * (sc > 0).mean()
    ct_pyro_mean[ct] = mean_all
    print(f'  {ct:22s} {n:6d} {mean_all:13.4f} {mean_pos:14.4f} {pct:10.1f}', flush=True)

print('\n--- Pyroptosis module score RANKING (higher = more pyroptosis expr) ---', flush=True)
for ct, v in sorted(ct_pyro_mean.items(), key=lambda x: -x[1]):
    print(f'  {ct:22s} {v:8.4f}', flush=True)
print('  => Microglia vs Astrocytes ratio:',
      f"{ct_pyro_mean.get('Microglia',0)/ct_pyro_mean.get('Astrocytes',1e-9):.2f}x", flush=True)

print('\n=== STEP 4: ASTROCYTE pyroptosis genes SNI vs SI (gene-level) ===', flush=True)
a_mask = ct_of_cell == 'Astrocytes'
for g in pyro:
    if g not in gene_expr: continue
    arr = gene_expr[g][a_mask]
    si = arr[cond_of_cell[a_mask] == 'SI']
    sni = arr[cond_of_cell[a_mask] == 'SNI']
    if len(si) >= 3 and len(sni) >= 3:
        fc = np.mean(sni) - np.mean(si)
        _, p = stats.mannwhitneyu(sni, si)
        sig = '***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else 'ns'
        print(f'  {g:10s} astro SI={np.mean(si):7.3f} SNI={np.mean(sni):7.3f} delta={fc:+.3f} {sig} (nSI={len(si)},nSNI={len(sni)})', flush=True)

print('\n=== STEP 5: MICROGLIA pyroptosis genes SNI vs SI (contrast) ===', flush=True)
m_mask = ct_of_cell == 'Microglia'
for g in pyro:
    if g not in gene_expr: continue
    arr = gene_expr[g][m_mask]
    si = arr[cond_of_cell[m_mask] == 'SI']
    sni = arr[cond_of_cell[m_mask] == 'SNI']
    if len(si) >= 3 and len(sni) >= 3:
        fc = np.mean(sni) - np.mean(si)
        _, p = stats.mannwhitneyu(sni, si)
        sig = '***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else 'ns'
        print(f'  {g:10s} micro SI={np.mean(si):7.3f} SNI={np.mean(sni):7.3f} delta={fc:+.3f} {sig}', flush=True)

print('\n=== STEP 6: ANIMAL-LEVEL inference (between-animal variation) ===', flush=True)
def animal_level(genes, ct):
    """For each animal, mean pyroptosis module score across its cells of this ct;
       then compare SI animals vs SNI animals (animal is the replicate)."""
    mask = ct_of_cell == ct
    sc = module_score(genes, mask)  # per-cell pyro module
    if len(sc) == 0: return None
    ani = animal_of_cell[mask]
    out = {}
    for a in np.unique(ani):
        out[a] = np.mean(sc[ani == a])
    si_vals = [v for a,v in out.items() if samples[a]=='SI']
    sni_vals = [v for a,v in out.items() if samples[a]=='SNI']
    return si_vals, sni_vals, out

for ct in ['Astrocytes','Microglia']:
    r = animal_level(pyro, ct)
    if not r: 
        print(f'  {ct}: no cells', flush=True); continue
    si_v, sni_v, out = r
    print(f'\n  [{ct}] animal-level pyroptosis module (mean per animal):', flush=True)
    for a in sorted(out):
        print(f'    {a} ({samples[a]}): {out[a]:.4f}', flush=True)
    if len(si_v) >= 2 and len(sni_v) >= 2:
        _, p = stats.mannwhitneyu(sni_v, si_v)
        diff = np.mean(sni_v) - np.mean(si_v)
        sig = '***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else 'ns'
        print(f'    ANIMAL-LEVEL SI={np.mean(si_v):.4f} SNI={np.mean(sni_v):.4f} delta={diff:+.4f} {sig} '
              f'(nSI={len(si_v)},nSNI={len(sni_v)})', flush=True)
    else:
        print(f'    insufficient animals (nSI={len(si_v)},nSNI={len(sni_v)}) for test', flush=True)

print('\n=== DONE ===', flush=True)
